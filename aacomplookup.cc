
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
  std::string index;
  ostream *out;
  bool delout;
  bool count;
  bool aggregate;
  double tol;
  Options(int argc, char *argv[]);
  ~Options();
  void usage(char *msg=NULL);
};

Options::Options(int argc, char *argv[]) : out(&cout) {
  signed char c;
  index = "";
  optarg = NULL;
  tol = 0.1;
  count = false;
  aggregate = false;
    
  while ((c = getopt(argc, argv, "i:o:hct:a")) != -1)
    switch (c) {
    case 'i':
      index = optarg;
      break;
    case 'o':
      if (delout) {
	delete out;
      }
      out = new ofstream(optarg);
      delout = true;
      break;
    case 't':
      tol = atof(optarg);
      break;
    case 'c':
      count = true;
      break;
    case 'a':
      aggregate = true;
      break;
    case 'h':
    default :
      usage();
    }
  if (index == "") usage();
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
  cerr << "Usage: aacomplookup [options] \n\n";
  cerr << "Options: \n";
  cerr << "  -i <aa-comp-index>     AA composition index. Required.\n";
  cerr << "  -o <output-file>       Output file. Defaults to standard out.\n";
  cerr << "  -t <mass-tolerance>    Mass tolerance.\n";
  cerr << "  -c                     Output counts.\n";
  cerr << "  -a                     Output aggregate counts.\n";
  cerr << "  -h                     Command line option help.\n";
  cerr << "\n";
  exit(1);
}

int main(int argc,char *argv[]) {

  set_new_handler(newfail_handler);
  
  Options opt(argc,argv);
  
  typedef sortedvector<float,std::pair<long unsigned int,FILE_POSITION_TYPE> > massindex;
  massindex mass2pos;
  massindex::iterator m2pit;
  
  ifstream ifs(opt.index.c_str());
  mass2pos.bread(ifs);
  ifs.close();

  if (mass2pos.empty()) {
    exit(0);
  }

  sortedvector<long unsigned int,long unsigned int> protstats;
  sortedvector<long unsigned int,long unsigned int>::iterator pstit;

  sortedvector<long unsigned int,char> protset;
  sortedvector<long unsigned int,char>::iterator psit;

  int i=0;
  double m;
  m2pit = mass2pos.end();
  while (cin >> m) {
    if (m <= 0) break;
    protset.clear();
    try {
      if (m2pit == mass2pos.end()) {
	m2pit = mass2pos.locate_first_at_least(m-opt.tol);
      } else {
	m2pit = mass2pos.finger_locate_first_at_least(m2pit,m-opt.tol);
      }
    } 
    catch (massindex::KeyOutOfRange const & e) {
      m2pit = mass2pos.end();
    }
    catch (massindex::InvalidFinger const & e) {
      m2pit = mass2pos.end();
    }
    while (m2pit != mass2pos.end() && m2pit->key() <= m+opt.tol) {
      *opt.out << i << " " << m2pit->value().first << " " << m2pit->value().second << endl;
      ++m2pit;
    }
    i++;
  }

  return 0;
}

