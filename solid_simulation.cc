

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
  std::string database;
  int reads;
  char eos_char;
  ostream *out;
  bool delout;
  bool verbose;
  bool memmap;
  int mersize;
  int dbind;
  bool noindex;
  bool rc;
  int samples;
  int rounds;
  int period;
  std::string errprob;
  Options(int argc, char *argv[]);
  ~Options();
  void usage(char *msg=NULL);
};

Options::Options(int argc, char *argv[]) : out(&cout) {
  signed char c;
  optarg = NULL;
  eos_char = '\n';
  delout = false;
  verbose = false;
  rc = false;
  memmap = true;
  noindex = false;
  mersize = 0;
  samples = 0;
  dbind = 0;
  errprob = "0 0 0 0 0.1";
  rounds = 5;
  period = 5;
  while ((c = getopt(argc, argv, "i:o:r:p:E:hvBIm:RS:e:")) != -1)
    switch (c) {
    case 'm':
      mersize = atoi(optarg);
      break;
    case 'r':
      rounds = atoi(optarg);
      break;
    case 'p':
      period = atoi(optarg);
      break;
    case 'i':
      database = optarg;
      break;
    case 'e':
      errprob = optarg;
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
    case 'S':
      samples = atoi(optarg);
      break;
    case 'E': {
      int eos;
      if (!sscanf(optarg,"%i",&eos)) {
	usage("Invalid end-of-sequence specification.\n");
      }
      eos_char = (char)eos;
    }
    break;
    case 'v':
      verbose = true;
      break; 
    case 'R':
      rc = true;
      break; 
    case 'B':
      memmap = false;
      break; 
    case 'I':
      noindex = true;
      break; 
    case 'h':
    default :
      usage();
    }
  if (database == "" || mersize == 0) usage();
  if (dbind < 0 || dbind > 4) usage("Invalid integer for fasta database indexing (-D).");
  if (rounds * period != mersize) {
    usage("Inconsistent parameters");
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
  cerr << "Usage: xmers [options] \n\n";
  cerr << "Options: \n";
  cerr << "  -m <int>               Mersize of mers.\n";
  cerr << "  -R                     Reverse complement too.\n";
  cerr << "  -i <sequence-database> Input sequence database. Required.\n";
  cerr << "  -o <output-file>       Output file. Defaults to standard out.\n";
  cerr << "  -E <int>               End-of-sequence character. Default is \'\\n\'\n";
  cerr << "  -I                     Do not load fasta database index.\n";
  cerr << "  -v                     Verbose (version & diagnostic) output.\n";
  cerr << "  -B                     Don\'t use memmap for I/O, use buffered I/O instead.\n";
  cerr << "  -D (0|1|2|3|4)         Fasta database indexing and preprocessing.\n";
  cerr << "                         0: Auto, 1: None, 2: Indexed, 3: Normalized,\n";
  cerr << "                         4: Compressed. Default: 0.\n";
  cerr << "  -h                     Command line option help.\n";
  cerr << "\n";
  exit(1);
}

char getcschar(char f, char s) {
  // checkpoint;
  // cerr << f << " " << s << endl;
  switch (f) {
  case 'A':
    switch (s) {
    case 'A':
      return '0';
    case 'C':
      return '1';
    case 'G':
      return '2';
    case 'T':
      return '3';
    }
    break;
  case 'C':
    switch (s) {
    case 'A':
      return '1';
    case 'C':
      return '0';
    case 'G':
      return '3';
    case 'T':
      return '2';
    }
    break;
  case 'G':
    switch (s) {
    case 'A':
      return '2';
    case 'C':
      return '3';
    case 'G':
      return '0';
    case 'T':
      return '1';
    }
    break;
  case 'T':
    switch (s) {
    case 'A':
      return '3';
    case 'C':
      return '2';
    case 'G':
      return '1';
    case 'T':
      return '0';
    }
    break;
  }
  assert(0);
}

void tocs(char* a, int n, char *b) {
  b[0] = 'G';
  b[1] = getcschar('G',a[0]);
  for (int i=1;i<n;i++) {
    b[i+1] = getcschar(a[i-1],a[i]);
  }
  b[n+1] = 0;
}

char rndcs(char c) {
  char c1=c;
  while (c1 == c) {
    c1 = '0'+(char)(floor(DRAND48*4));
  }
  return c1;
}

void applyerror(char* a, float *errprob, int rounds, int period) {
  for (int r=0;r<rounds;r++) {
    float ep = errprob[r];
    if (ep > 0) {
      for (int p=0;p<period;p++) {
	if (DRAND48 < ep) {
	  char c = a[1+r*period+p];
	  a[1+r*period+p] = rndcs(c);
	} 
      }
    }
  }
}

int main(int argc,char *argv[]) {
  
  set_new_handler(newfail_handler);
#if defined(PROFILE)
  set_profile_signal_handler();
#endif
  // assert(sizeof(FILE_POSITION_TYPE)>=8);
  // assert(sizeof(bigword)>=8);
  assert(sizeof(FILE_POSITION_TYPE)>=8);

  Options opt(argc,argv);

  if (opt.verbose) {
    ostrstream ss;
    ss << "Release Tag: " << release_tag << ends;
    std::string v(ss.str());
    timestamp(v.c_str());
  }

  float eprob[opt.rounds];
  istrstream epstr(opt.errprob);

  float p;
  int i=0;
  while (i < opt.rounds && epstr >> p) {
    eprob[i] = p;
    i++;
  }

  // checkpoint;
  //   for (int i=0;i<opt.rounds;i++) {
  //     cerr << eprob[i] << " ";
  //   }
  //   cerr << endl;

  FastaFile<Lazy_Header_SI>* ff;
  fasta_file_seq_params ffp;
  ffp.upper_case = true;
  ffp.eos_char = opt.eos_char;
  ffp.check_params = true;
  ff = pick_fasta_file<Lazy_Header_SI>(opt.database,opt.dbind,opt.memmap,
				       !opt.noindex,ffp,opt.verbose);
  
  char neos=ff->nch(opt.eos_char);

  double s0;
  FILE_POSITION_TYPE s;
  FILE_POSITION_TYPE size = ff->length();
  char *buffer = new char[opt.mersize+1];
  char *buffer1 = new char[opt.mersize+2];
  char *buffer2 = new char[opt.mersize+2];

  srand48(time(NULL));
  i=0;
  while (i < opt.samples) {
    s0 = DRAND48;
    s = ((FILE_POSITION_TYPE)(s0*size));
    // cerr << s << endl;
    if (opt.rc && DRAND48 > 0.5) {
      if (s < opt.mersize) {
	continue;
      }
      bool good = true;
      ff->pos(s-opt.mersize);
      for (int j=opt.mersize-1;j>=0;j--) {
	if (ff->eof()) {
	  good = false;
	  break;
	}
	int nch = ff->getnch();
	if (nch > 3) {
	  good = false;
	  break;
	}
	buffer[j] = iupac_revcomp(ff->ch(nch));
      }
      if (good) {
	buffer[opt.mersize] = '\0';
	tocs(buffer,opt.mersize,buffer1);
	// checkpoint;
	// cerr << buffer1 << endl;
	applyerror(buffer1,eprob,opt.rounds,opt.period);
	// checkpoint;
	// cerr << buffer1 << endl;
	*opt.out << ">" << -i << " " << buffer << "\n" << buffer1 << endl;
	i++;
      }
    } else {
      bool good = true;
      ff->pos(s);
      for (int j=0;j<opt.mersize;j++) {
	if (ff->eof()) {
	  good = false;
	  break;
	}
	int nch = ff->getnch();
	if (nch > 3) {
	  good = false;
	  break;
	}
	buffer[j] = ff->ch(nch);
      }
      if (good) {
	buffer[opt.mersize] = '\0';
	tocs(buffer,opt.mersize,buffer1);
	// checkpoint;
	// cerr << buffer1 << endl;
	applyerror(buffer1,eprob,opt.rounds,opt.period);
	// checkpoint;
	// cerr << buffer1 << endl;
	*opt.out << ">" << i << " " << buffer << "\n" << buffer1 << endl;
	i++;
      }
    }
  }

  return 0;
}

