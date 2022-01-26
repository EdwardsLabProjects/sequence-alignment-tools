
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
  char eos_char;
  double min;
  double max;
  double randprob;
  
  Options(int argc, char *argv[]);
  ~Options();
  void usage(char *msg=NULL);
};

Options::Options(int argc, char *argv[]) : out(&cout) {
  signed char c;
  optarg = NULL;
  min = 100;
  max = 3000;
  randprob = 1;
    
  while ((c = getopt(argc, argv, "i:o:hm:M:r:")) != -1)
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
      min = atof(optarg);
      break;
    case 'M':
      max = atof(optarg);
      break;
    case 'r':
      randprob = atof(optarg);
      break;
    case 'h':
    default :
      usage();
    }
  if (database == "") usage();
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
  cerr << "Usage: aacomp [options] \n\n";
  cerr << "Options: \n";
  cerr << "  -i <sequence-database> Input sequence database. Required.\n";
  cerr << "  -o <output-file>       Output file. Defaults to standard out.\n";
  cerr << "  -m <mass>              Min peptide mass to index.\n";
  cerr << "  -M <mass>              Max peptide mass to index.\n";
  cerr << "  -r <prob>              Random sampling probability.\n";
  cerr << "  -h                     Command line option help.\n";
  cerr << "\n";
  exit(1);
}

int main(int argc,char *argv[]) {

  set_new_handler(newfail_handler);
  
  Options opt(argc,argv);
  
  long unsigned int i=0;
  FILE_POSITION_TYPE pos=0;
  typedef sortedvector<float,std::pair<long unsigned int,FILE_POSITION_TYPE> > massindex;
  massindex mass2pos(10000);
  massindex::iterator m2pit;

  srand48(time(NULL));
  bool delete_stream = false;
  ifstream ifs(opt.database.c_str());
  fasta_entry f;
  while (ifs >> f) {
    if (f.sequence() == "") break;
    if (drand48() > opt.randprob) continue;
    timestampi("Fasta entry: ",i);
    char const * const seq=f.sequence().c_str();
    char const * seq1=seq;
    long unsigned int seqlen=f.sequence().length();
    FILE_POSITION_TYPE s,l;
    float m0,m1 = 0;
    for (s=0;s<seqlen;s++) {
      m0 = 0;
      seq1 = seq+s;
      for (l=0;((l<opt.max/50)&&(s+l<seqlen));l++) {
	// cerr << seq1[l] << " " << monomolwt(seq1[l]) << " ";
	if ((m1=monomolwt(seq1[l])) >= 0) {
	  m0 += m1;
	  // cerr << m0 << " " << i << endl;
	  if (m0 >= opt.min && m0 <= opt.max) {
	    mass2pos.push_back(m0,std::make_pair(i,pos+s+l/2));
	  }
	} else {
	  break;
	}
      }
    }
    pos += seqlen;
    i = i+1;
  }
  mass2pos.normalize_strict();
  mass2pos.bwrite(*opt.out);

  return 0;
}

