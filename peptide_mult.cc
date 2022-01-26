
#include <unistd.h>
#include <assert.h>
#include <vector>
#include <map>
#include <algorithm>
#include <sstream>
#include "keyword_tree.h"
#include "keyword_tree.t"
#include "shift_and.h"
#include "fasta_io.h"
#include "fasta_io.t"
#include "char_io.t"
#include "primer_alignment.h"
#include "util.h"
#include "types.h"
#include "select.t"
#include "math.h"

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

char release_tag[] = "$Name:  $";

#if defined(__alpha) 
// Hack to get around problems with mmap
const unsigned long __sbrk_override = 1;
#endif

void newfail_handler() {
  cerr << "Fatal Error: Attempt to allocate memory with new() failed.\n" 
       << endl;
  exit(1);
  return;
}

struct Options {
  std::string indfile;
  std::string database;
  std::string massfile;
  double tol;
  bool reltol;
  int miscl;
  ostream *out;
  bool delout;
  char eos_char;
  bool verbose;
  bool memmap;
  int dbind;
  bool xwild;
  
  Options(int argc, char *argv[]);
  ~Options();
  void usage(char *msg=NULL);
};

Options::Options(int argc, char *argv[]) : out(&cout) {
  signed char c;
  optarg = NULL;
  dbind = 0;
  eos_char = '\n';
  memmap = true;
  delout = false;
  verbose = false;
  tol = 2;
  reltol = false;
  miscl = 1;
  xwild = false;
    
  while ((c = getopt(argc, argv, "i:o:hm:I:BD:e:rC:Xv")) != -1)
    switch (c) {
    case 'i':
      database = optarg;
      break;
    case 'o':
      if (delout) {
	delete out;
      }
      out = new ofstream(optarg);
      delout = true;
      break;
    case 'm':
      massfile = optarg;
      break;
    case 'I':
      indfile = optarg;
      break;
    case 'D':
      dbind = atoi(optarg);
      break;
    case 'C':
      miscl = atoi(optarg);
      break;
    case 'e':
      tol = atof(optarg);
      break;
    case 'E': {
      int eos_char0;
      if (!sscanf(optarg,"%i",&eos_char0)) {
	usage("Invalid end-of-sequence specification.\n");
      }
      eos_char = (char)eos_char0;
    }
    break;
    case 'r':
       reltol= true;
      break; 
    case 'v':
      verbose = true;
      break; 
    case 'B':
      memmap = false;
      break; 
    case 'X':
      xwild = true;
      break; 
    case 'h':
    default :
      usage();
    }
  if ((indfile == "" || database == "" || massfile == "")&&!verbose) usage("One of protein indices, sequence database, or mass file is missing.");
  if (dbind < 0 || dbind> 4) usage("Invalid integer for fasta database indexing (-D).");
}


Options::~Options() {
  if (delout) {
    delete out;
  }
}

void Options::usage(char *message) {
  if (message != NULL && strlen(message) > 0) {
    cerr << message << endl;
    cerr << endl;
  }
  cerr << "Usage: peptide_mult [options] \n\n";
  cerr << "Options: \n";
  cerr << "  -I <protein-indices>   Indices of proteins, and the queries' \n";
  cerr << "                         mw for multiplicities. Required. \n";
  cerr << "  -i <sequence-database> Input sequence database. Required.\n";
  cerr << "  -m <mass-file>         File of masses. Required\n";
  cerr << "  -o <output-file>       Output file. Defaults to standard out.\n";
  cerr << "  -e <tolerance>         Query MW error tolerance. Default 2.\n";
  cerr << "  -r                     Error tolerance is relative, not absolute.\n";
  cerr << "  -C                     Number of missed cleavages permitted.\n";
  cerr << "  -X                     Treat X as wildcard for other aas.\n";
  cerr << "  -E <int>               End-of-sequence character. Default is \'\\n\'\n";
  cerr << "  -B                     Don\'t use memmap for I/O, use buffered I/O instead.\n";
  cerr << "  -D (0|1|2|3|4)         Fasta database indexing and preprocessing.\n";
  cerr << "                         0: Auto, 1: None, 2: Indexed, 3: Normalized,\n";
  cerr << "                         4: Compressed. Default: Auto.\n";
  cerr << "  -v                     Verbose (version & diagnostic) output.\n";
  cerr << "  -h                     Command line option help.\n";
  cerr << "\n";
  exit(1);
}

struct query_elt {
  long unsigned int query;
  long unsigned int rank;
};

struct pep_elt {
  FILE_POSITION_TYPE spos;
  FILE_POSITION_TYPE epos;
  long unsigned int miscl;
  char sub[2];
};

#define max_obs(m,t,r) ((r)?((m)*(1+(t))):((m)+(t)))
#define min_obs(m,t,r) ((r)?((m)*(1-(t))):((m)-(t)))

int main(int argc,char *argv[]) {

  set_new_handler(newfail_handler);
  
  Options opt(argc,argv);

  // checkpoint;

  FastaFile<Lazy_Header_SI>* ff;
  fasta_file_seq_params ffp;
  ffp.upper_case = true;
  ffp.eos_char = opt.eos_char;
  ffp.eos_start = true;
  ff = pick_fasta_file<Lazy_Header_SI>(opt.database,opt.dbind,opt.memmap,
				       /*alignments=*/true,ffp,opt.verbose);
  
  // checkpoint;

  std::vector<double> masses(ff->size(),0.0);
  double nterm=0.0,cterm=0.0;
  ifstream mifs(opt.massfile.c_str());
  char symbol;
  double mass;
  while (mifs) {
    mifs >> symbol >> mass;
    if (symbol == 'n') {
      nterm = mass;
    } else if (symbol == 'c') {
      cterm = mass;
    } else if (symbol >= 'A' && symbol <= 'Z') {
      if (ff->nch(symbol) >= 0) {
	masses[ff->nch(symbol)] = mass;
      }
    }
  }
  mifs.close();

  if (opt.xwild) {
    masses[ff->nch('X')] = 0.0;
    masses[ff->nch('B')] = 0.0;
    masses[ff->nch('Z')] = 0.0;
  }

  // checkpoint;

  unsigned char neos = ff->nch(opt.eos_char);
  unsigned char nx = ff->nch('X');
  unsigned char nb = ff->nch('B');
  unsigned char nz = ff->nch('Z');
  unsigned char nm = ff->nch('M');

  long unsigned int const sz=ff->size();

  std::vector<bool> trpair(ff->size()*ff->size(),false);
#define trpair_index(i,j) ((i)*sz+(j))
  for (long unsigned int i=0;i<ff->size();i++) {
    trpair[trpair_index(neos,i)] = true;
    trpair[trpair_index(i,neos)] = true;
    trpair[trpair_index(ff->nch('K'),i)] = true;
    trpair[trpair_index(ff->nch('R'),i)] = true;
  }
  trpair[trpair_index(ff->nch('K'),ff->nch('P'))] = false;
  trpair[trpair_index(ff->nch('R'),ff->nch('P'))] = false;

#define tryptic(ch1,ch2) (trpair[trpair_index((ch1),(ch2))])

  // checkpoint;

//   for (int i=0;i<masses.size();i++) {
//     cerr << i << " " << ff->ch(i) << " " << masses[i] << endl;
//   }
  
  //   cerr << (int)neos << endl;
  
  unsigned char c,lc;
  bool delete_stream = false;
  istream * ifs(&cin);
  if (opt.indfile != "-") {
    ifs = new ifstream(opt.indfile.c_str());
    delete_stream = true;
  }
  unsigned long index;
  const int BUFSIZE=8*1024;
  char *buffer = new char[BUFSIZE];
  long unsigned int buffer1size=BUFSIZE;
  unsigned char *buffer1=new unsigned char[buffer1size];
  long unsigned int input_chars;
  do {
    
    // Read in query...

    // checkpoint;

    sortedvector<double,query_elt> mws;
    sortedvector<double,query_elt>::iterator mwit,mwit1;

    ifs->getline(buffer,BUFSIZE-1);
    if (ifs->fail()) break;
    assert(ifs->gcount()<BUFSIZE-1);
    buffer[ifs->gcount()]='\0';
    istrstream ss(buffer);
    ss >> index;
    if (index == 0) break;
    while (ss) {
      query_elt qe;
      double molwt=0;
      ss >> qe.query >> qe.rank >> molwt;
      if (molwt <= 0) break;
      mws.push_back(molwt,qe);
    }
    mws.normalize();
    
    // checkpoint;
    //     cerr << index << " " 
    // 	 << mws.rbegin()->key() << " " 
    // 	 << mws.begin()->key() << " " 
    // 	 << endl;

    // Read in sequence...

    FILE_POSITION_TYPE len;
    FILE_POSITION_TYPE pos=0;
    ff->fasta_pos(index-1,0);
    const Lazy_Header_SI & h = ff->get_header_data(ff->pos());
    buffer1[pos++] = neos;
    while (1) {
      while ((c=ff->getnch()) != neos && pos < buffer1size-3) {
	buffer1[pos++] = c;
      }
      buffer1[pos++] = c;
      // for (long unsigned int i=1;i<pos;i++) {
      // cerr << ff->ch(buffer1[i]);
      // }
      // cerr << endl;
      if ( pos < buffer1size-3 ) {
	break;
      } 
      unsigned char *buffer2=new unsigned char[buffer1size*2];
      memcpy(buffer2,buffer1,buffer1size);
      buffer1size*=2;
      if (buffer1) delete [] buffer1;
      buffer1 = buffer2;
    }
    // checkpoint;
    buffer1[pos++] = '\0';

//     checkpoint;
//     cerr << ">" << h.header() << endl;
//     for (long unsigned int i=1;i<pos;i++) {
//       cerr << ff->ch(buffer1[i]);
//     }
//     cerr << endl;
    
    typedef sortedvector<double,pep_elt> peplist;
    peplist pepmw;
    peplist::iterator pit,pit1;
    long unsigned int spos=1,epos=1;
    do {
      epos = spos-1;
      // checkpoint;
      // cerr << "spos: " << spos << " epos: " << epos << endl;
      long int mcl=-1;
      if (spos == 1 /* && buffer1[spos] == nm */) {
	mcl--;
      }
      mass = nterm;
      int numx = 0,numb = 0,numz = 0;
      long signed int wcpos=-1;
      // mass += masses[buffer1[spos]];
      while (1) {
	// checkpoint;
	epos++;
	if (buffer1[epos] == neos) break;
	mass += masses[buffer1[epos]];
	if (buffer1[epos] == nx) {
	  numx++;
	  wcpos = epos;
	}
	if (buffer1[epos] == nb) {
	  numb++;
	  wcpos = epos;
	}
	if (buffer1[epos] == nz) {
	  numz++;
	  wcpos = epos;
	}
// 	checkpoint;
// 	cerr << "Adding " << masses[buffer1[epos]] 
// 	     << " (" << ff->ch(buffer1[epos]) << ")" 
// 	     << " = " << mass
// 	     << endl;
	// cerr << "spos: " << spos << " epos: " << epos << endl;
	while (epos > 1 && !tryptic(buffer1[epos],buffer1[epos+1]) && 
	       buffer1[epos+1]!= neos) {
	  epos++;
// 	  checkpoint;
// 	  cerr << "Adding " << masses[buffer1[epos]] 
// 	       << " (" << ff->ch(buffer1[epos]) << ")" 
// 	       << " = " << mass
// 	       << endl;
	  mass += masses[buffer1[epos]];
	  if (buffer1[epos] == nx) {
	    numx++;
	    wcpos = epos;
	  }
	  if (buffer1[epos] == nb) {
	    numb++;
	    wcpos = epos;
	  }
	  if (buffer1[epos] == nz) {
	    numz++;
	    wcpos = epos;
	  }
	}
// 	if (opt.verbose) {
// 	  // checkpoint;
// 	  cerr << mass+cterm << '\t';
// 	  cerr << spos << '\t';
// 	  cerr << epos << '\t';
// 	  cerr << mcl+1 << '\t';
// 	  for (long unsigned int i=spos;i<=epos;i++) {
// 	    cerr << ff->ch(buffer1[i]);
// 	  }
// 	  cerr << endl;
// 	}
	mcl++;
	if (mcl > opt.miscl) {
	  break;
	}
	if (!opt.xwild || numx + numb + numz == 0) {
	  if (mass + cterm > max_obs(mws.rbegin()->key(),opt.tol,opt.reltol)) {
	    break;
	  }
	  if (mass + cterm < min_obs(mws.begin()->key(),opt.tol,opt.reltol)) {
	    continue;
	  }
	  pep_elt pe;
	  pe.spos = spos;
	  pe.epos = epos;
	  pe.miscl = mcl;
	  pe.sub[0] = 0;
	  pe.sub[1] = 0;
	  pepmw.push_back(mass+cterm,pe);
	} else if (numx + numb + numz == 1) {
	  static int const naa = 20;
	  static char aas_for_x[naa+1]= "ACDEFGHIKLMNPQRSTVWY";
	  static char aas_for_b[3]= "ND";
	  static char aas_for_z[3]= "EQ";
	  int limit = naa;
	  if (numx == 0) limit = 2;
	  for (int i=0;i<limit;i++) {
	    double m  = mass + cterm;
	    if (numx > 0) m += masses[ff->nch(aas_for_x[i])];
	    if (numb > 0) m += masses[ff->nch(aas_for_b[i])];
	    if (numz > 0) m += masses[ff->nch(aas_for_z[i])];
	    if (m > max_obs(mws.rbegin()->key(),opt.tol,opt.reltol)) {
	      continue;
	    } 
	    if (m < min_obs(mws.begin()->key(),opt.tol,opt.reltol)) {
	      continue;
	    } 
	    if (numx > 0 && 
		(aas_for_x[i] == 'K' || aas_for_x[i] == 'R') && 
		(buffer1[wcpos+1] != ff->nch('P')) &&
		mcl+1 > opt.miscl) {
	      continue;
	    }
	    pep_elt pe;
	    if (numx > 0) {
	      pe.sub[0] = 'X';
	      pe.sub[1] = aas_for_x[i];
	    }
	    if (numb > 0) {
	      pe.sub[0] = 'B';
	      pe.sub[1] = aas_for_b[i];
	    }
	    if (numz > 0) {
	      pe.sub[0] = 'Z';
	      pe.sub[1] = aas_for_z[i];
	    }
	    pe.spos = spos;
	    pe.epos = epos;
	    pe.miscl = mcl;
	    pepmw.push_back(m,pe);
	  }
	}
      }
      // checkpoint;
      // cerr << "spos: " << spos << " epos: " << epos << endl;
      spos++;
      while (((spos == 2 && buffer1[spos-1]!=nm) || spos > 2) && 
	     !tryptic(buffer1[spos-1],buffer1[spos]) && 
	     buffer1[spos]!=neos)
	spos++;
      // checkpoint;
      // cerr << "spos: " << spos << " epos: " << epos << endl;
    } while (buffer1[spos] != neos);
    // checkpoint;
    // cerr << "spos: " << spos << " epos: " << epos << endl;
    pepmw.normalize();

    if (opt.verbose) {
      pit = pepmw.begin();
      while (pit != pepmw.end()) {
	cerr << pit->key() << '\t';
	cerr << pit->value().spos << '\t';
	cerr << pit->value().epos << '\t';
	cerr << pit->value().miscl << '\t';
	for (long unsigned int i=pit->value().spos;i<=pit->value().epos;i++) {
	  cerr << ff->ch(buffer1[i]);
	}
	if (pit->value().sub[0] != 0) {
	  cerr << '\t';
	  cerr << pit->value().sub[0] << "=" << pit->value().sub[1];
	}
	cerr << endl;
	++pit;
      }
      mwit = mws.begin();
      while (mwit != mws.end()) {
	cerr << mwit->key() << " \\in "
	     << "[" << min_obs(mwit->key(),opt.tol,opt.reltol)
	     << "," << max_obs(mwit->key(),opt.tol,opt.reltol)
	     << "]\n";
	++mwit;
      }
    }

    // checkpoint;

    *opt.out << index;

    mwit = mws.begin();
    pit = pepmw.end();
    while (mwit != mws.end()) {
      // checkpoint;
      if (pit == pepmw.end()) {
	// checkpoint;
	// cerr << mwit->key() << endl;
	try {
	  pit1 = pit = pepmw.locate_first_at_least(min_obs(mwit->key(),opt.tol,opt.reltol));
	} 
	catch (peplist::KeyOutOfRange const &) {
	  // checkpoint;
	  pit1 = pit = pepmw.end();
	}
      } else {
	// checkpoint;
	// cerr << mwit->key() << endl;
	// cerr << min_obs(mwit->key(),opt.tol,opt.reltol) << endl;
	try {
	  pit1 = pit = pepmw.finger_locate_first_at_least(pit,min_obs(mwit->key(),opt.tol,opt.reltol));
	} 
	catch (peplist::KeyOutOfRange const &) {
	  // checkpoint;
	  pit1 = pit = pepmw.end();
	}
	catch (peplist::InvalidFinger const & e) {
	  // checkpoint
	  pit1 = pit = pepmw.end();
	}
      }
      // checkpoint;
      long unsigned int mult = 0;
      while (pit1!=pepmw.end() && 
	     pit1->key() <= max_obs(mwit->key(),opt.tol,opt.reltol)) {
	// cerr << "found: " << pit1->key() << endl;
	++mult;
	++pit1;
      }
      // checkpoint;
      *opt.out << '\t' << mwit->value().query 
	       << '\t' << mwit->value().rank
	       // << '\t' << mwit->key()
	       << '\t' << mult;
      if (mult == 0 && 1/*opt.verbose*/) {
	cerr << "Query " << mwit->value().query << ", Rank "
	     << mwit->value().rank << ", has multiplicity " << mult
	     << " for protein entry " << h.index() << ":\n>" << h.header() << "\n";
      }
      // checkpoint;
      ++mwit;
    }
    *opt.out << endl;
  } while (!ifs->eof());
  if (delete_stream) {
    delete ifs;
  }

  return 0;
}

