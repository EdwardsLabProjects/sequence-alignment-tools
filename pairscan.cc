

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
  std::string output;
  int verbose;
  int mersize;
  string temp;
  bool ignore;
  int distmin;
  int distmax;
  int exitthresh;
  Options(int argc, char *argv[]);
  ~Options();
  void usage(char *msg=NULL);
};
Options::Options(int argc, char *argv[]) {
  signed char c;
  optarg = NULL;
  verbose = 0;
  mersize = 0;
  distmin = 0;
  distmax = 1000;
  ignore = false;
  exitthresh = -1;
  while ((c = getopt(argc, argv, "i:o:hvm:l:d:D:e:")) != -1)
    switch (c) {
    case 'm':
      mersize = atoi(optarg);
      break;
    case 'd':
      distmin = atoi(optarg);
      break;
    case 'D':
      distmax = atoi(optarg);
      break;
    case 'I':
      ignore = true;
      break;
    case 'i':
      database = optarg;
      break;
    case 'o':
      output = optarg;
      break;
    case 'e':
      exitthresh = atoi(optarg);
      break;
    case 'l':
      stderr = fopen(optarg,"a");
      setlinebuf(stderr);
      break;
    case 'v':
      verbose++;
      break; 
    case 'h':
    default :
      usage();
    }
  if (database == "" || output == "" || mersize == 0) usage();
}


Options::~Options() {
  ;
}

void Options::usage(char *message) {
  if (message != NULL && strlen(message) > 0) {
    cerr << message << endl;
    cerr << endl;
  }
  cerr << "Usage: pairscan [options] \n\n";
  cerr << "Options: \n";
  cerr << "  -m <int>               Mersize of mers.\n";
  cerr << "  -d <int>               Min. distance between 3' ends of mers\n";
  cerr << "  -D <int>               Max. distance between 3' ends of mers\n";
  cerr << "  -i <sequence-database> Input sequence database. Required.\n";
  cerr << "  -o <output-file>       Output file. Required.\n";
  cerr << "  -l <log-file>          Redirect stderr.\n";
  cerr << "  -v                     Verbose (version & diagnostic) output.\n";
  cerr << "  -h                     Command line option help.\n";
  cerr << "\n";
  exit(1);
}

void
progress(CharacterProducer const & cp1, bitmap const & first, bitmap const & again, int interval, bool force=false) {
  static int nextprogress=0;
  if (nextprogress == 0 || force) {
    nextprogress = time(NULL);
  }
  if (time(NULL) >= nextprogress) {
    ostrstream ss;
    ss.setf(ios::fixed);
    ss << "Progress:";
    ss << setprecision(3) << setw(7) << cp1.progress()*100;
    ss << "%";
    ss << " Seen:"; 
    ss << setprecision(3) << setw(7) << 100.0*first.nset()/first.size();
    ss << "%";
    ss << " Dupe:"; 
    ss << setprecision(3) << setw(7) << 100.0*again.nset()/again.size();
    ss << "% ";
    ss << setw(10) << again.nunset();
    ss << " left.";
    ss << ends;
    std::string v(ss.str());
    timestamp(v.c_str());
    nextprogress += interval;
  }
}

void 
write_chkpnt(int span, uint32 pos, bitmap const & first, bitmap const & again, std::string & filename, int interval, bool force=false) {
  static int nextprogress=0;
  if (nextprogress == 0) {
    nextprogress = time(NULL)+interval;
  }
  if (time(NULL) >= nextprogress || force) {
    ostrstream ss;
    ss << "BEGIN" << endl;
    ss << span << " " << pos << endl;
    first.write(ss);
    again.write(ss);
    ss << "END" << endl;
    
    string s(ss.str());

    // Atomic write of match state
    ofstream ofs(filename.c_str());
    ofs.write(s.c_str(),s.length());
    ofs.close();
    
    nextprogress += interval;
  }
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

  Normalized<BufferedFileChars> cp(opt.database+".seq",'$',0,10000000);
  
  long unsigned int hashsize = (((hash_t)1)<<(opt.mersize*4));
  int interval = 5;
  int chkpnt_write_interval = 60;
  
  bitmap first(hashsize);
  bitmap again(hashsize);

  ifstream first_in(opt.output.c_str());
  uint32 posin = 0;
  uint32 spanin = opt.distmin;
  if (first_in) {
    string l;
    assert(getline(first_in,l) && l == "BEGIN");
    first_in >> spanin >> posin;
    assert(getline(first_in,l));
    first.read(first_in);
    again.read(first_in);
    assert(getline(first_in,l) && l == "END");
  }
  first_in.close();

  if (opt.ignore) {
    posin = 0;
    spanin = opt.distmin;
  }

  for (int i=spanin;i<=opt.distmax;i++) {

    char *temp=new char[i+2*opt.mersize+1];
    int p=0;
    for (int j=0;j<opt.mersize;j++) {
      temp[p++] = '1';
    }
    for (int j=0;j<i;j++) {
      temp[p++] = '0';
    }
    for (int j=0;j<opt.mersize;j++) {
      temp[p++] = '1';
    }
    temp[p] = 0;
    
    hash * hptr = hashselect(cp,4,temp);
    hash & h = *hptr;

    if (i == spanin && posin > 0) {
      h.reset(posin+1);
    } else {
      h.reset();
    }

    if (opt.verbose) {
      progress(cp,first,again,interval,true);
    }

    while (h.next()) {
      hash_t v = h.value();
      hash_t v1 = h.rcvalue();

      if (!first.get(v)) {
	first.set(v);
	if (opt.verbose >= 2) {
	  string vs = h.str(v);
	  fprintf(stderr," . %s %lu %d %s - %s F\n",
		  opt.database.c_str(),h.pos(),i,vs.substr(0,opt.mersize).c_str(),
		  vs.substr(opt.mersize,opt.mersize).c_str());
	}
      } else {
	if (!again.get(v)) {
	  again.set(v);
	  if (opt.verbose >= 2) {
	    string vs = h.str(v);
	    fprintf(stderr," + %s %lu %d %s - %s F\n",
		    opt.database.c_str(),h.pos(),i,vs.substr(0,opt.mersize).c_str(),
		    vs.substr(opt.mersize,opt.mersize).c_str());
	  }
	}
      }

      v = h.rcvalue();
      if (!first.get(v)) {
	first.set(v);
	if (opt.verbose >= 2) {
	  string vs = h.str(v);
	  fprintf(stderr," . %s %lu %d %s - %s R\n",
		  opt.database.c_str(),h.pos(),i,vs.substr(0,opt.mersize).c_str(),
		  vs.substr(opt.mersize,opt.mersize).c_str());
	}
      } else {
	if (!again.get(v)) {
	  again.set(v);
	  if (opt.verbose >= 2) {
	    string vs = h.str(v);
	    fprintf(stderr," + %s %lu %d %s - %s R\n",
		    opt.database.c_str(),h.pos(),i,vs.substr(0,opt.mersize).c_str(),
		    vs.substr(opt.mersize,opt.mersize).c_str());
	  }
	}
      }

      if (opt.verbose) {
	progress(cp,first,again,interval);
      } 
      write_chkpnt(i,cp.pos(),first,again,opt.output,chkpnt_write_interval);
    }

    if (opt.verbose) {
      progress(cp,first,again,interval,true);
    }

    delete hptr;

    if (opt.exitthresh > 0 && again.nunset() < opt.exitthresh) {
      write_chkpnt(0,0,first,again,opt.output,chkpnt_write_interval,true);
      exit(2);
    }

    write_chkpnt(i+1,0,first,again,opt.output,chkpnt_write_interval,true);

  }

  if (opt.verbose) {
    progress(cp,first,again,interval,true);
  }

  write_chkpnt(0,0,first,again,opt.output,chkpnt_write_interval,true);

  if (opt.exitthresh > 0 && again.nunset() < opt.exitthresh) {
    exit(2);
  }

  return 0;
}

