

#include <unistd.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include "keyword_tree.h"
#include "keyword_tree.t"
#include "shift_and_inexact.h"
#include "rand_hash_table.h"
#include "shift_and.h"
#include "fasta_io.h"
#include "fasta_io.t"
#include "char_io.t"
#include "primer_alignment.h"
#include "util.h"
#include "types.h"
#include "select.h"
#include "select.t"
#include "gs_hash_table.h"

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

char release_tag[] = "$Name:  $";

#if defined(__alpha) 
// Hack to get around problems with mmap
const unsigned long __sbrk_override = 1;
#endif

void newfail_handler() {
  timestamp("Fatal Error: Attempt to allocate memory with new() failed.");
  exit(1);
  return;
}

struct Options {
  ostream *out;
  bool delout;
  bool verbose;
  FILE_POSITION_TYPE length;
  Options(int argc, char *argv[]);
  ~Options();
  void usage(char *msg=NULL);
};

Options::Options(int argc, char *argv[]) : out(&cout) {
  signed char c;
  optarg = NULL;
  length = 0;
  delout = false;
  verbose = false;
  while ((c = getopt(argc, argv, "o:E:l:hv")) != -1)
    switch (c) {
    case 'l':
      length = atol(optarg);
      break;
    case 'o':
      if (delout) {
	delete out;
      }
      out = new ofstream(optarg);
      delout = true;
      break;
    case 'v':
      verbose = true;
      break; 
    case 'h':
    default :
      usage();
    }
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
  cerr << "Usage: genome_simulation [options] \n\n";
  cerr << "Options: \n";
  cerr << "  -l <length>            Number of characters.\n";
  cerr << "  -o <output-file>       Output file. Defaults to standard out.\n";
  cerr << "  -E <int>               End-of-sequence character. Default is \'\\n\'\n";
  cerr << "  -v                     Verbose (version & diagnostic) output.\n";
  cerr << "  -h                     Command line option help.\n";
  cerr << "\n";
  exit(1);
}

int main(int argc,char *argv[]) {
  
  set_new_handler(newfail_handler);
#if defined(PROFILE)
  set_profile_signal_handler();
#endif

  Options opt(argc,argv);

  srand48(time(NULL));
  FILE_POSITION_TYPE i=0;
  *opt.out << ((unsigned char)4);
  while (i++ < opt.length) {
    unsigned char s0 = (unsigned char)floor(DRAND48*4);
    *opt.out << s0;
  }
  *opt.out << ((unsigned char)4);
}

