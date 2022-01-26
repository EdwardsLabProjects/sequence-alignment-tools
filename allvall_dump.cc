

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
  int chunksize;
  int mersize;
  bool verbose;
  Options(int argc, char *argv[]);
  ~Options();
  void usage(char *msg=NULL);
};
Options::Options(int argc, char *argv[]) {
  signed char c;
  optarg = NULL;
  mersize = 0;
  chunksize = 0;
  verbose = false;

  while ((c = getopt(argc, argv, "i:m:d:C:vh")) != -1)
    switch (c) {
    case 'C':
      chunksize = atoi(optarg);
      break;
    case 'm':
      mersize = atoi(optarg);
      break;
    case 'i':
      database = optarg;
      break;
    case 'd':
      datfile = optarg;
      break;
    case 'v':
      verbose++;
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
  cerr << "Usage: allvall_dump [options] \n\n";
  cerr << "Options: \n";
  cerr << "  -m <int>               Mersize of mers.\n";
  cerr << "  -i <sequence-database> Input sequence database. Required.\n";
  cerr << "  -d <datfile>           Input bitmap file. Required.\n";
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

  int dbblock=0;
  int p = opt.database.rfind('.');
  if (p != std::string::npos) {
    dbblock = atoi(opt.database.substr(p+1).c_str());
  }

  FILE_POSITION_TYPE dboffset = ((FILE_POSITION_TYPE)opt.chunksize)*dbblock;
  if (opt.verbose)
    timestampli("database offset: ",dboffset);

  FastaFile<Lazy_Header_SI>* ff, *ff1;
  fasta_file_seq_params ffp;
  ffp.upper_case = true;
  ffp.eos_char = '$';
  ffp.check_params = false;
  ffp.translate = false;
  ffp.eos_start = true;
  ffp.offset = dboffset;
  ff  = pick_fasta_file<Lazy_Header_SI>(opt.database,0,false,true,ffp,opt.verbose);

  FastaFile<Lazy_Header_SI> & cp  = *ff;

  bitmap match(cp.length());

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

  bitmap::runs_t runs;
  match.runs(runs,false);

  for (int i=0;i<runs.size();i++) {
    FILE_POSITION_TYPE pos = runs[i].first;
    uint32 len = runs[i].second;
    // cerr << pos << " " << len << endl;
    if (len == 0) {
      continue;
    }
    cp.pos(pos+dboffset-opt.mersize);
    string seq = cp.getstr(len+opt.mersize-1);
    const Lazy_Header_SI & h  = cp.get_header_data(pos+dboffset);
    FILE_POSITION_TYPE st = cp.get_seq_pos(pos+dboffset)-opt.mersize;
    FILE_POSITION_TYPE ed = st+len+opt.mersize-1;
    cout << ">" << h.header() << " /run=" << dbblock << "." << i << " /pos=" << pos+dboffset << " /index=" << h.index() << " /start=" << st << " /end=" << ed << " /len=" << len+opt.mersize-1 << endl;
    cout << seq << endl;
  }
  /*
  for (int i=0;i<cp.length();i++) {
    if (!match[i]) {
      cp.pos(i-opt.mersize);
      cout << i << " " << cp.getstr((unsigned int)(opt.mersize)) << endl;
    }
  }
  */
  return 0;
}

