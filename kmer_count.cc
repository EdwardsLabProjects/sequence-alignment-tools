
#include <unistd.h>
#include <assert.h>
#include <vector>
#include <map>
#include <algorithm>
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
  std::string database;
  ostream *out;
  bool delout;
  bool memmap;
  int node;
  int dbind;
  bool wc;
  bool tn;
  char eos_char;
  bool verbose;
  int k;
  double t;
  int T;
  int m;
  bool aggregate;
  bool nmer;
  bool binary;
  bool nonacgtmer;
  bool addrc;
  Options(int argc, char *argv[]);
  ~Options();
  void usage(const char *msg=NULL);
};

Options::Options(int argc, char *argv[]) : out(&cout) {
  signed char c;
  optarg = NULL;
  eos_char = '\n';
  memmap = true;
  dbind = 0;
  verbose = false;
  k = 1;
  t = 0;
  T = 0;
  m = MAXINT;
  aggregate = false;
  nmer = false;
  nonacgtmer = false;
  binary = false;
  addrc = false;
    
  while ((c = getopt(argc, argv, "i:o:E:hBD:wWvk:t:T:M:anNbr")) != -1)
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
    case 'w':
      wc = true;
      tn = false;
      break;
    case 'W':
      wc = true;
      tn = true;
      break;
    case 'k':
      k = atoi(optarg);
      break;
      /*    case 't':
      t = atof(optarg);
      break;
    case 'T':
      T = atoi(optarg);
      break; */
    case 'M':
      m = atoi(optarg);
      break;
    case 'D':
      dbind = atoi(optarg);
      break;
    case 'E': {
      int eos_char0;
      if (!sscanf(optarg,"%i",&eos_char0)) {
	usage("Invalid end-of-sequence specification.\n");
      }
      eos_char = (char)eos_char0;
    }
    break;
    case 'v':
      verbose = true;
      break; 
    case 'r':
      addrc = true;
      break; 
    case 'n':
      nmer = true;
      break; 
    case 'N':
      nonacgtmer = true;
      nmer = true;
      break; 
    case 'a':
      aggregate = true;
      break; 
    case 'b':
      binary = true;
      break; 
    case 'B':
      memmap = false;
      break; 
    case 'h':
    default :
      usage();
    }
  if (database == "") usage();
  if (dbind < 0 || dbind> 4) usage("Invalid integer for fasta database indexing (-D).");
  if (binary && ! aggregate) usage("Cannot output binary data in non-aggregate mode.");
}


Options::~Options() {
  if (delout) {
    delete out;
  }
}

void Options::usage(const char *message) {
  if (message != NULL && strlen(message) > 0) {
    cerr << message << endl;
    cerr << endl;
  }
  cerr << "Usage: kmer_count [options] \n\n";
  cerr << "Options: \n";
  cerr << "  -i <sequence-database> Input sequence database. Required.\n";
  cerr << "  -o <output-file>       Output file. Defaults to standard out.\n";
  cerr << "  -k <mer-size>          k-mer size.\n";
  cerr << "  -n                     Include Ns in k-mers.\n";
  cerr << "  -N                     Consider any non ACGT as N.\n";
  cerr << "  -a                     Aggregate counts.\n";
  cerr << "  -b                     Binary output for aggregate counts.\n";
  cerr << "  -r                     Aggregate forward and reverse complement counts.\n";
  cerr << "  -M <max-output>        Max number of mers to output.\n";
  cerr << "  -E <int>               End-of-sequence character. Default is \'\\n\'\n";
  cerr << "  -B                     Don\'t use memmap for I/O, use buffered I/O instead.\n";
  cerr << "  -D (0|1|2|3|4)         Fasta database indexing and preprocessing.\n";
  cerr << "                         0: Auto, 1: None, 2: Indexed, 3: Normalized,\n";
  cerr << "                         4: Compressed. Default: 0.\n";
  cerr << "  -v                     Verbose (version & diagnostic) output.\n";
  cerr << "  -h                     Command line option help.\n";
  cerr << "\n";
  exit(1);
}

long unsigned int hash(const char *mer, const long unsigned int len, 
		       const int base, const int charmap[]) {
  long unsigned int h=0;
  const char *mer1 = mer;
  for (FILE_POSITION_TYPE i=0;i<len;i++) {
    int m = charmap[*mer1];
    assert(m >= 0 && m < base);
    h = h*base+m;
    mer1++;
  }
  return h;
}

long unsigned int hashrc(const char *mer, const long unsigned int len, 
			 const int base, const int rccharmap[]) {
  long unsigned int h=0;
  const char *mer1 = mer;
  for (FILE_POSITION_TYPE i=0;i<len;i++) {
    int m = rccharmap[*mer1];
    assert(m >= 0 && m < base);
    h = h*base+m;
    mer1--;
  }
  return h;
}

char *unhash(long unsigned int h, int len, int base, 
	     char invcharmap[], char *mer) {
  for (int i=len-1;i>=0;i--) {
    mer[i] = invcharmap[h%base];
    h/=base;
  }
  return mer;
}

class mypair : public std::pair<long unsigned int, long unsigned int> {
public:
  mypair() {};
  mypair(long unsigned int f, long unsigned int s) 
    : std::pair<long unsigned int, long unsigned int>(f,s) {};
  mypair(std::pair<long unsigned int, long unsigned int> const & a) 
    : std::pair<long unsigned int, long unsigned int>(a) {};
};

bool operator<(const mypair & a, const mypair & b) {
  if (a.second < b.second) return true;
  if (a.second > b.second) return false;
  if (a.first < b.first) return true;
  if (a.first > b.first) return false;
  return false;
};

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

  FastaFile<Header_SI> *ff;
  fasta_file_seq_params ffp;
  ffp.upper_case = true;
  ffp.eos_char = opt.eos_char;
  ff = pick_fasta_file<Header_SI>(opt.database,opt.dbind,opt.memmap,
				  /*alignments=*/true,ffp,opt.verbose);

  // checkpoint;
  
  long unsigned int fe=0;
  char ch;
  char *mer = new char[(opt.k+1)*opt.k];
  for (int i=0;i<opt.k;i++) {
    mer[i*(opt.k+1)+opt.k] = '\0';
  }
  int *pos= new int[opt.k];
  for (int i=0;i<opt.k;i++) {
    pos[i]=-i;
  }
  long unsigned int h;
  std::map<long unsigned int,unsigned int> counts;
  std::map<long unsigned int,unsigned int>::iterator ctit;
  std::multimap<long unsigned int,long unsigned int> orderedcounts;
  std::multimap<long unsigned int,long unsigned int>::iterator orctit;
  std::multimap<long unsigned int,long unsigned int>::reverse_iterator orctrit;
  
  int charmap[256];
  int rccharmap[256];
  for (int i=0;i<256;i++) {
    charmap[i] = -1;
    rccharmap[i] = -1;
  }
  int nchar;
  char *invcharmap;
  if (opt.nmer) {
    charmap['A'] = 0;  rccharmap['T'] = 0;
    charmap['C'] = 1;  rccharmap['G'] = 1;
    charmap['G'] = 2;  rccharmap['C'] = 2;
    charmap['T'] = 3;  rccharmap['A'] = 3;
    charmap['N'] = 4;  rccharmap['N'] = 4;
    nchar= 5; 
    invcharmap = new char[nchar];
    invcharmap[0] = 'A';
    invcharmap[1] = 'C';
    invcharmap[2] = 'G';
    invcharmap[3] = 'T';
    invcharmap[4] = 'N';
  } else {
    charmap['A'] = 0;  rccharmap['T'] = 0;
    charmap['C'] = 1;  rccharmap['G'] = 1;
    charmap['G'] = 2;  rccharmap['C'] = 2;
    charmap['T'] = 3;  rccharmap['A'] = 3;
    nchar= 4; 
    invcharmap = new char[nchar];
    invcharmap[0] = 'A';
    invcharmap[1] = 'C';
    invcharmap[2] = 'G';
    invcharmap[3] = 'T';
  }

  double totkmers=pow((double)nchar,opt.k);
  double q=((totkmers-1)/totkmers);
  double expectation;
  double stddev;
  FILE_POSITION_TYPE seqlen=0;
  // checkpoint;
  while (!ff->eof()) {
    // checkpoint;
    ch = ff->getch();
    if (opt.nonacgtmer && ch != 'A' && ch != 'C' && 
	ch != 'G' && ch != 'T' && ch != opt.eos_char) {
      ch = 'N';
    } 
    // cerr << ch << endl;
    if (ch == opt.eos_char) {
      /* print out stats */
      if (!opt.aggregate) {
	bool firstresult=true;
	expectation=((seqlen/opt.k)/totkmers);
	stddev=((seqlen/opt.k)*q/totkmers);
	// cout << "Expected count: " << expectation << endl;
	// cout << "Stddev count: " << stddev << endl;
	for (ctit=counts.begin();ctit!=counts.end();++ctit) {
	  if (ctit->second // && 
	      // ctit->second > expectation+stddev*opt.t &&
	      /*ctit->second > opt.T*/) {
	    orderedcounts.insert(std::pair<const unsigned long,unsigned long>((const unsigned long)ctit->second,ctit->first));
	  }
	}
	int numberoutput=0;
	for (orctrit=orderedcounts.rbegin();orctrit!=orderedcounts.rend();++orctrit) {
	  if (firstresult) {
	    Header_SI const & h=ff->get_header_data(ff->pos()-1);
	    *opt.out << '>' << h.header() << endl;
	    firstresult=false;
	  }
	  *opt.out << unhash(orctrit->second,opt.k,nchar,invcharmap,mer) << " " 
	       << orctrit->first << endl;
	  numberoutput++;
	  if (numberoutput >= opt.m) break;
	}
	counts.clear();
	orderedcounts.clear();
      }
      for (int i=0;i<opt.k;i++) {
	pos[i]=-i;
      }
      seqlen=0;
      continue;
    }
    // checkpoint;
    if (charmap[ch] < 0) {
      for (int i=0;i<opt.k;i++) {
	pos[i]=-i;
      }
      continue;
    }
    // checkpoint;
    for (int i=0;i<opt.k;i++) {
      if (pos[i]>=0) {
	mer[i*(opt.k+1)+pos[i]] = ch;
      }
    }
    seqlen++;
    // checkpoint;
    for (int i=0;i<opt.k;i++) {
      if (pos[i]==(opt.k-1)) {
	h = hash(&mer[i*(opt.k+1)],opt.k,nchar,charmap);
	if ((ctit=counts.find(h)) != counts.end()) {
	  ctit->second++;
	} else {
	  counts[h]=1;
	}
	if (opt.addrc) {
	  h = hashrc(&mer[i*(opt.k+1)],opt.k,nchar,rccharmap);
	  if ((ctit=counts.find(h)) != counts.end()) {
	    ctit->second++;
	  } else {
	    counts[h]=1;
	  }	  
	}
	pos[i] = 0;
      } else {
	pos[i]++;
      }
    }
  }
  
  if (opt.aggregate) {
    ctit = counts.begin();
    while (ctit!=counts.end()) {
      if (!opt.binary) {
	if (ctit->second > 0) {
	  *opt.out << unhash(ctit->first,opt.k,nchar,invcharmap,mer) << " " 
		   << ctit->second << endl;
	}
      } else {
	if (ctit->second > 0) {
	  long unsigned int b1;
	  unsigned int b2;
	  (*opt.out).write((const char *)(&ctit->first),sizeof(long unsigned int));
	  (*opt.out).write((const char *)(&ctit->second),sizeof(unsigned int));
	}
      }
      ctit++;
    }
  }
  delete ff;
  
  return 0;
}

