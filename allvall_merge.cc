

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
  std::string outfile;
  std::list<std::string> bmin;
  bool verbose;
  bool ascii;
  Options(int argc, char *argv[]);
  ~Options();
  void usage(char *msg=NULL);
};
Options::Options(int argc, char *argv[]) {
  signed char c;
  optarg = NULL;
  verbose = false;
  outfile = "";
  ascii = false;

  while ((c = getopt(argc, argv, "vAo:")) != -1)
    switch (c) {
    case 'o':
      outfile = optarg;
      break;
    case 'v':
      verbose = true;
      break;
    case 'A':
      ascii = true;
      break;
    case 'h':
    default :
      usage();
    }
  if (outfile == "") usage();

  for (int i=optind;i<argc;i++) {
    bmin.push_back(argv[i]);
  }
  
  if (bmin.size() < 1) usage();

}


Options::~Options() {
  ;
}

void Options::usage(char *message) {
  if (message != NULL && strlen(message) > 0) {
    cerr << message << endl;
    cerr << endl;
  }
  cerr << "Usage: allvall_merge [options] bitmap files \n";
  cerr << "Options: \n";
  cerr << "  -o <output-bitmap>     Output bitmap, - implies stdout. Required.\n";
  cerr << "  -A                     Ascii bitmap out. Default: False.\n";
  cerr << "  -v                     Verbose. Default: False.\n";
  cerr << "  -h                     Command line option help.\n";
  cerr << "\n";
  exit(1);
}

bool readbm(bitmap & match, std::string filename) {
  istream* match_in;
  if (filename != "-") {
    match_in = new ifstream(filename.c_str());
  } else {
    match_in = &cin;
  }
  match_in->peek();
  match_in->peek();
  FILE_POSITION_TYPE posin = 0;
  uint32 spanin = 0;
  bool good = false;
  if (*match_in) {
    string l;
    assert(getline(*match_in,l) && l == "BEGIN");
    (*match_in) >> spanin >> posin;
    assert(getline(*match_in,l));
    match.read(*match_in);
    assert(getline(*match_in,l) && l == "END");
    good = true;
  }
  if (filename != "-") {
    ((ifstream*)match_in)->close();
  }
  return good;
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

  std::list<std::string>::const_iterator it=opt.bmin.begin();

  bitmap match;
  match.set_nocheck();
  if (opt.verbose) {
    timestamps("Reading in bitmap: ",it->c_str());
  }
  readbm(match,*it);

  if (opt.verbose) {
    timestampli("Bitmap bits unset: ",match.nunset());
  }

  if (*it != "-") {
    ++it;
  }
  while (it != opt.bmin.end()) {
    bitmap match1;
    match1.set_nocheck();

    if (opt.verbose) {
      timestamps("Reading in bitmap: ",it->c_str())
    }
    
    bool good = readbm(match1,*it);
    if (!good) {
      break;
    }
    if (opt.verbose) {
      timestampli("Bitmap bits unset: ",match1.nunset());
    }

    match |= match1;

    if (opt.verbose) {
      timestampli("Total unset:       ",match.nunset());
    }

    if (*it != "-") {
      ++it;
    }
  }
  timestampli("Final unset:       ",match.nunset());

  ostrstream ss;
  ss << "BEGIN" << endl;
  ss << 0 << " " << 0 << endl;
  match.write(ss,opt.ascii);
  ss << "END" << endl;
  
  string s(ss.str());
  
  // Atomic write of match state
  if (opt.outfile != "-") {
    ofstream ofs(opt.outfile.c_str());
    ofs.write(s.c_str(),s.length());
    ofs.close();
  } else {
    cout.write(s.c_str(),s.length());
  }

  return 0;
}

