
#include <unistd.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include "keyword_tree.h"
#include "keyword_tree.t"
#include "shift_and_inexact.h"
#include "shift_and.h"
#include "fasta_io.h"
#include "fasta_io.t"
#include "char_io.t"
#include "primer_alignment.h"
#include "util.h"
#include "types.h"
#include "select.h"
#include "select.t"

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
  bool rev_comp;
  bool pattern_file;
  bool fasta_pattern_file;
  std::string patterns;
  std::string database;
  int nmismatch;
  char eos_char;
  std::string alignformat;
  bool alignments;
  ostream *out;
  bool delout;
  bool verbose;
  bool memmap;
  int dbind;
  int node;
  bool wc;
  bool tn;
  int minmotifcount;
  int minmotiflen;
  Options(int argc, char *argv[]);
  ~Options();
  void usage(char *msg=NULL);
};

Options::Options(int argc, char *argv[]) : out(&cout) {
  signed char c;
  optarg = NULL;
  rev_comp = false;
  pattern_file = false;
  fasta_pattern_file = false;
  eos_char = '\n';
  nmismatch = 0;
  delout = false;
  alignformat = ">%h\n %T %s %e\n %A\n %Q %i%R\n";
  alignments = true;
  verbose = false;
  memmap = true;
  dbind = 0;
  wc = false;
  tn = false;
  minmotifcount = -1;
  minmotiflen = -1;
  node =0;
    
  while ((c = getopt(argc, argv, "p:i:o:E:hrvA:BD:wWN:c:l:")) != -1)
    switch (c) {
    case 'p':
      patterns = optarg;
      break;
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
    case 'D':
      dbind = atoi(optarg);
      break;
    case 'E': {
      int eos_char;
      if (!sscanf(optarg,"%i",&eos_char)) {
	usage("Invalid end-of-sequence specification.\n");
      }
      eos_char = (char)eos_char;
    }
    break;
    case 'v':
      verbose = true;
      break; 
    case 'B':
      memmap = false;
      break; 
    case 'h':
    default :
      usage();
    }
  if (patterns == "" || database == "") usage();
  if (dbind < 0 || dbind > 4) usage("Invalid integer for fasta database indexing (-D).");
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
  cerr << "Usage: tandem_match [options] \n\n";
  cerr << "Options: \n";
  cerr << "  -p <sequences>         Tandem repeat motifs, separated by whitespace.\n";
  cerr << "  -i <sequence-database> Input sequence database. Required.\n";
  cerr << "  -o <output-file>       Output file. Defaults to standard out.\n";
  cerr << "  -E <int>               End-of-sequence character. Default is \'\\n\'\n";
  cerr << "  -v                     Verbose (version & diagnostic) output.\n";
  cerr << "  -B                     Don\'t use memmap for I/O, use buffered I/O instead.\n";
  cerr << "  -D (0|1|2|3|4)         Fasta database indexing and preprocessing.\n";
  cerr << "                         0: Auto, 1: None, 2: Indexed, 3: Normalized,\n";
  cerr << "                         4: Compressed. Default: 0.\n";
  cerr << "  -h                     Command line option help.\n";
  cerr << "\n";
  exit(1);
}

void alignformat(ostream & os,
		 std::string const & alignformat, 
		 FILE_POSITION_TYPE const & s,
		 FILE_POSITION_TYPE const & e,
		 FILE_POSITION_TYPE const & five,
		 FILE_POSITION_TYPE const & three,
		 FILE_POSITION_TYPE const & S,
		 FILE_POSITION_TYPE const & E,
		 long unsigned int const & i,
		 unsigned int const & d,
		 std::string const & p,
		 std::string const & P,
		 std::string const & q,
		 std::string const & Q,
		 std::string const & r,
		 std::string const & R,
		 std::string const & t,
		 std::string const & T,
		 std::string const & A,
		 std::string const & h,
		 std::string const & H,
		 long unsigned int const & f) {

  //
  // The following codes get transformed into positional markers for a
  // printf format.
  // 
  // %s    start of alignment from the start of the fasta entry
  // %e    end of alignment from the start of the fasta entry
  // %S    start of alignment in compressed sequence file
  // %E    end of alignment in compressed sequence file
  // %i    index of pattern
  // %d    edit distance of alignment
  // %p    tandem repeat motif
  // %Q    aligned tandem repeat with insertion characters
  // %t    aligned text
  // %T    aligned text with alignment characters
  // %A    alignment string
  // %h    fasta header of hit sequence
  // %H    first word of fasta header of hit sequence
  // %%    percent character
  //

  int pos=0;
  while (pos<alignformat.length()) {
    if (alignformat[pos] == '%') {
      pos++;
      if (pos < alignformat.length()) {
	switch(alignformat[pos]) {
	case 's': 
	  os << s;
	  break;
	case 'e':
	  os << e;
	  break;
	case 'l':
	  os << e-s;
	  break;
	case 'S': 
	  os << S;
	  break;
	case 'E':
	  os << E;
	  break;
	case 'i':
	  os << i;
	  break;
	case 'd':
	  os << d;
	  break;
	case 'p':
	  os << p;
	  break;
	case 'P':
	  os << P;
	  break;
	case 'q':
	  os << q;
	  break;
	case 'Q':
	  os << Q;
	  break;
	case 'r':
	  os << r;
	  break;
	case 'R':
	  os << R;
	  break;
	case 't':
	  os << t;
	  break;
	case 'T':
	  os << T;
	  break;
	case 'A':
	  os << A;
	  break;
	case 'h':
	  os << h;
	  break;
	case 'H':
	  os << H;
	  break;
	case 'f':
	  os << f;
	  break;
	case '%':
	  os << "%";
	  break;
	default:
	  os << alignformat[pos];
	}
      } else {
	os << "%";
      }
    } else if (alignformat[pos] == '\\') {
      pos++;
      if (pos < alignformat.length()) {
	switch(alignformat[pos]) {
	case 'n': 
	  os << '\n';
	  break;
	case 't': 
	  os << '\t';
	  break;
	case '\\': 
	  os << '\\';
	  break;
	default:
	  os << alignformat[pos];
	}
      } else {
	os << '\\';
      }
    } else {
      os << alignformat[pos];
    }
    pos++;
  }
}

int main(int argc,char *argv[]) {
  
  set_new_handler(newfail_handler);
  assert(sizeof(FILE_POSITION_TYPE)>=8);

  Options opt(argc,argv);

  if (opt.verbose) {
    ostrstream ss;
    ss << "Release Tag: " << release_tag << ends;
    std::string v(ss.str());
    timestamp(v.c_str());
  }

  std::list<std::string> patterns;
  std::list<std::string> patdeflines;
  std::list<std::string>::iterator patit;
  std::list<std::string>::iterator patdefit;

  istrstream sis(opt.patterns.c_str());
  std::string pattern;
  while (sis >> pattern) {
    patterns.push_back(pattern);
  }
  
  if (patterns.empty()) {
    exit(0);
  }

  for (patit=patterns.begin();patit!=patterns.end();++patit) {
    uppercase(*patit);
  }

  FastaFile<Header_SI>* ff;
  ff = pick_fasta_file<Header_SI>(opt.database,opt.dbind,opt.memmap,
				  opt.alignments,/*ucdict=*/true,
				  opt.eos_char,opt.verbose);
  
  //    PatternMatch* pm = new shift_and_inexact(opt.nmismatch,
  //  					   opt.eos_char,
  //  					   opt.wc,opt.tn);
  
  // PatternMatch* pm = new shift_and(opt.wc,opt.tn);
  // 

  PatternMatch* pm;
  pm = pick_pattern_index(ff,opt.node,opt.nmismatch,0,0,-1,opt.wc,opt.tn,true,opt.eos_char,opt.verbose);
  
  std::string pattern;
  std::list<std::string>::size_type i=1,n=patterns.size();
  unsigned long int N1=((!opt.rev_comp)?1:2)*n;
  std::vector<std::string> patarray(N1+1); /* indices 1..n, index 0 wasted */
  std::vector<std::string> patdefarray;
  if (opt.fasta_pattern_file) {
    patdefarray.resize(n+1);
    patdefit=patdeflines.begin();
  }

  for (patit=patterns.begin();patit!=patterns.end();++patit) {
    if (opt.verbose) {
      ostrstream ss;
      ss << "Tandem motif: " << setw(3) << i << " > " 
	 << *patit << ends;
      std::string v(ss.str());
      timestamp(v.c_str());
    }
    patarray[i]=*patit;
    if (opt.fasta_pattern_file) {
      patdefarray[i]=*patdefit;
      patdefit++;
    }
    
    pm->add_pattern(patarray[i],i);

    if (opt.rev_comp) {
      patarray[i+n]=reverse_comp(*patit);

      if (opt.verbose) {
	ostrstream ss;
	ss << "Tandem motif: " << setw(3) << i << " < " 
	   << patarray[i+n] << ends;
	std::string v(ss.str());
	timestamp(v.c_str());
      }
      pm->add_pattern(patarray[i+n],i+n);
    }
    i++;
  }
  
  if (opt.verbose) {
    pm->progress_interval(*ff);
  }
  pm->init(*ff);

  pattern_hit_vector l(1000);
  pattern_hit_vector::iterator it;
  bool more;
  // checkpoint;
  while ((more=pm->find_patterns(*ff,l,1000))||!l.empty()) {
    // checkpoint;
    FILE_POSITION_TYPE oldcharspos;
    oldcharspos = ff->pos();
    it = l.begin();
    // checkpoint;
    while (it != l.end()) {
      pattern_hit* ph = it->value();
      // checkpoint;
      if (ph && ph->pattern_id()) {
	// checkpoint;
	// cerr << "ph: " << ph->pattern_id() << " " << ph->pos() << endl;
	std::list<pattern_hit*> repeats;
	std::list<std::list<pattern_hit*>::iterator> repeatsits;
	repeats.push_back(ph);
	repeatsits.push_back(it);
	std::list<pattern_hit*>::iterator it1(it);
	++it1;
	bool det=false;
	unsigned int motiflen;
	FILE_POSITION_TYPE next_pos;
	motiflen = patarray[ph->pattern_id()].length();
	next_pos = ph->pos()+motiflen;
	pattern_hit* ph1;
	// checkpoint;
	while (it1 != l.end()) {
	  // checkpoint;
	  if ((ph1=*it1)!=NULL) {
	    // checkpoint;
	    // cerr << next_pos << endl;
	    if (ph1->pos() > next_pos) {
	      // checkpoint;
	      det=true;
	      break;
	    }
	    if (ph1->pattern_id() == ph->pattern_id()) {
	      if (ph1->pos() == next_pos) {
		repeats.push_back(ph1);
		repeatsits.push_back(it1);
		next_pos = ph1->pos()+motiflen;
	      } else if (ph1->pos() < next_pos) {
		repeats.push_back(ph1);
		repeatsits.push_back(it1);
	      }
	    }
	  }
	  it1++;
	}
	// checkpoint;
	if (!det&&more) {
	  // checkpoint;
	  ++it;
	  continue;
	}
	// checkpoint;
	// At this point, repeats contains all the motif hits, one
	// after the other in order along the sequence...
	pattern_alignment *pa = 
	  new editdist_alignment(*repeats.back(),0,opt.eos_char,
				 opt.wc,opt.tn,0,0);
	// checkpoint;
	// We need to figure out what to align against. 
	std::string motif = patarray[ph->pattern_id()];
	int copyn = (repeats.back()->pos()
		     -repeats.front()->pos()
		     +motif.length())
	            /motif.length();
	// checkpoint;
	// cerr << copyn << endl;
	if (copyn >= opt.minmotifcount && 
	    repeats.back()->pos()-repeats.front()->pos() + motif.length() 
	    >= opt.minmotiflen) {
	  // checkpoint;
	  std::string alignstr=motif;
	  int i=1;
	  while (i < copyn) {
	    alignstr+=motif;
	    i++;
	  }
	  // checkpoint;
	  pa->align(*ff,alignstr);	
	  // checkpoint;
	  // cerr << pa->editdist() << endl;
	  if (pa->editdist() == 0) {
	    // checkpoint;
	    const FILE_POSITION_TYPE spe = ff->get_seq_pos(pa->end());
	    const FILE_POSITION_TYPE sps = spe-pa->length()+1;
	    const FILE_POSITION_TYPE pe = pa->end();
	    const FILE_POSITION_TYPE ps = pe-pa->length()+1;
	    const bool rc = (ph->pattern_id()>n);
	    const long unsigned int ind = repeats.front()->pattern_id()-(rc?n:0);
	    const unsigned int ed = pa->editdist();
	    const Header_SI & h = ff->get_header_data(pa->end());
	    std::string patdef("");
	    if (opt.fasta_pattern_file) {
	      patdef = patdefarray[ind];
	    }
	    alignformat(*opt.out,
			opt.alignformat,sps,spe,(rc?spe:sps),
			(rc?sps:spe),ps,pe,ind,ed,patarray[ind],
			patdef,alignstr,pa->alignment_pattern(alignstr),
			(rc?"R":"F"),(rc?" REVCOMP":""),
			pa->matching_text(),pa->alignment_text(),
			pa->alignment_string(),
			h.header(),h.short_header(),h.index());
	    // checkpoint;
	  }
	}
	delete pa;
	// checkpoint;
	// Delete repeats from list...
	std::list<std::list<pattern_hit*>::iterator>::iterator repitit;
	repitit=repeatsits.begin();
	// First one is special...
	++repitit;
	while (repitit!=repeatsits.end()) {
	  l.erase(*repitit);
	  ++repitit;
	}
	// checkpoint;
	std::list<pattern_hit*>::iterator repit;
	repit = repeats.begin();
	while (repit != repeats.end()) {
	  delete *repit;
	  ++repit;
	}
      }
      // checkpoint;
      // This takes care of the first one...
      it = l.erase(it);
      // checkpoint;
    }
    // checkpoint;
    ff->pos(oldcharspos);
    // checkpoint;
  }
  return 0;
}

