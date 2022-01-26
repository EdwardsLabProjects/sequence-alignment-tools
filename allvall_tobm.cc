

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

#include "hash.h"
#include "perfposht.h"
#include "bitmap.h"

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
  std::string database;
  std::string datfile;
  int mersize;
  int chunksize;
  uint32 offset;
  bool exclude;
  Options(int argc, char *argv[]);
  ~Options();
  void usage(char *msg=NULL);
};
Options::Options(int argc, char *argv[]) {
  signed char c;
  optarg = NULL;
  mersize = 0; 
  exclude = false;
  offset = 0;
  chunksize = 0;

  while ((c = getopt(argc, argv, "i:C:m:d:no:")) != -1)
    switch (c) {
    case 'm':
      mersize = atoi(optarg);
      break;
    case 'i':
      database = optarg;
      break;
    case 'C':
      chunksize = atoi(optarg);
      break;
    case 'd':
      datfile = optarg;
      break;
    case 'n':
      exclude = true;
      break;
    case 'o':
      offset = atoi(optarg);
      break;
    case 'h':
    default :
      usage();
    }
  if (database == "" || datfile == "" || mersize == 0) usage();
}


Options::~Options() {
  ;
}

void Options::usage(char *message) {
  if (message != NULL && strlen(message) > 0) {
    cerr << message << endl;
    cerr << endl;
  }
  cerr << "Usage: allvall_tobm [options] \n\n";
  cerr << "Options: \n";
  cerr << "  -m <int>               Mersize of mers.\n";
  cerr << "  -i <sequence-database> Input sequence database. Required.\n";
  cerr << "  -d <datfile>           Input bitmap file. Required.\n";
  cerr << "  -C <chunksize>         Chunksize.\n";
  cerr << "  -n                     Mark mer positions not included in ranges.\n";
  cerr << "  -o                     Position offset of bitmap vs ranges.\n";
  cerr << "  -h                     Command line option help.\n";
  cerr << "\n";
  exit(1);
}

int main(int argc,char *argv[]) {
  
  set_new_handler(newfail_handler);
#if defined(PROFILE)
  set_profile_signal_handler();
#endif
  assert(sizeof(uint64)==8);
  assert(sizeof(uint32)==4);
  assert(sizeof(int64)==8);
  assert(sizeof(int32)==4);

  Options opt(argc,argv);

  int block=0;
  int p = opt.database.rfind('.');
  if (p != std::string::npos) {
    block = atoi(opt.database.substr(p+1).c_str());
  }
  uint32 offset = block*opt.chunksize;
  if (opt.offset) {
    offset = opt.offset;
  }

  Normalized<BufferedFileChars> cp(opt.database+".seq",'$',0,10000000);

  bitmap match(cp.length()+1);

  ifstream match_in(opt.datfile.c_str());
  uint32 posin = 0;
  uint32 spanin = 0;
  if (match_in) {
    string l;
    assert(getline(match_in,l) && l == "BEGIN");
    match_in >> spanin >> posin;
    assert(getline(match_in,l));
    match.read(match_in);
    assert(getline(match_in,l) && l == "END");
  }
  match_in.close();

  uint32 nextmark=0;
  if (opt.exclude) {
    for (int i=0;i<opt.mersize-1;i++) {
      // cerr << "Mark position " << i << endl;
      match.set(i);
    }
    nextmark=opt.mersize-1;
  }

  int32 first;
  uint32 count;
  while (true) {
    cin >> first >> count;
    // cerr << first << " " << count << endl; 
    if (cin.eof()) {
      break;
    }
    first -= offset;
    if (opt.exclude) {
      if (first < opt.mersize-1 || first + count >= (cp.length()+1)) {
	continue;
      }
      for (int i=0;i<first-nextmark;i++) {
	if (nextmark+i < opt.mersize-1 || nextmark+i >= (cp.length()+1)) {
	  continue;
	}
	// cerr << "Mark position " << nextmark+i << endl;
	match.set(nextmark+i);
      }
      nextmark=first+count;
    } else {
      if (first + count < opt.mersize-1 || first >= (cp.length()+1)) {
	continue;
      }
      for (int i=0;i<count;i++) {
	// cerr << "Mark position " << first+i << endl;
	match.set(first+i);
      }
    }
  }
  if (opt.exclude) {
    for (int i=0;i<cp.length()+1-nextmark;i++) {
      if (nextmark+i < opt.mersize-1 || nextmark+i >= (cp.length()+1)) {
	continue;
      }
      // cerr << "Mark position " << nextmark+i << endl;
      match.set(nextmark+i);
    }
  }

  ostrstream ss;
  ss << "BEGIN" << endl;
  ss << 0 << " " << 0 << endl;
  match.write(ss);
  ss << "END" << endl;
  
  string s(ss.str());
  
  // Atomic write of match state
  ofstream ofs(opt.datfile.c_str());
  ofs.write(s.c_str(),s.length());
  ofs.close();

  return 0;
}

