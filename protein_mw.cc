
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
  std::string indfile;
  std::string database;
  std::string massfile;
  ostream *out;
  bool delout;
  char eos_char;
  bool verbose;
  bool memmap;
  int dbind;
  
  Options(int argc, char *argv[]);
  ~Options();
  void usage(char *msg=NULL);
};

Options::Options(int argc, char *argv[]) : out(&cout) {
  signed char c;
  optarg = NULL;
  dbind = 0;
  eos_char = '\n';
  memmap = true;
  delout = false;
  verbose = false;
    
  while ((c = getopt(argc, argv, "i:o:hm:I:BD:v")) != -1)
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
      massfile = optarg;
      break;
    case 'I':
      indfile = optarg;
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
  if ((indfile == "" || database == "" || massfile == "")&&!verbose) usage("One of protein indices, sequence database, or mass file is missing.");
  if (dbind < 0 || dbind> 4) usage("Invalid integer for fasta database indexing (-D).");
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
  cerr << "Usage: protein_mw [options] \n\n";
  cerr << "Options: \n";
  cerr << "  -I <protein-indices>   Indices of proteins to compute MW for. Required. \n";
  cerr << "  -i <sequence-database> Input sequence database. Required.\n";
  cerr << "  -m <mass-file>         File of masses. Required\n";
  cerr << "  -o <output-file>       Output file. Defaults to standard out.\n";
  cerr << "  -E <int>               End-of-sequence character. Default is \'\\n\'\n";
  cerr << "  -B                     Don\'t use memmap for I/O, use buffered I/O instead.\n";
  cerr << "  -D (0|1|2|3|4)         Fasta database indexing and preprocessing.\n";
  cerr << "                         0: Auto, 1: None, 2: Indexed, 3: Normalized,\n";
  cerr << "                         4: Compressed. Default: Auto.\n";
  cerr << "  -v                     Verbose (version & diagnostic) output.\n";
  cerr << "  -h                     Command line option help.\n";
  cerr << "\n";
  exit(1);
}

int main(int argc,char *argv[]) {

  set_new_handler(newfail_handler);
  
  Options opt(argc,argv);

  // checkpoint;

  FastaFile<Lazy_Header_SI>* ff;
  fasta_file_seq_params ffp;
  ffp.upper_case = true;
  ffp.eos_char = opt.eos_char;
  ffp.eos_start = true;
  ff = pick_fasta_file<Lazy_Header_SI>(opt.database,opt.dbind,opt.memmap,
				       /*alignments=*/true,ffp,opt.verbose);
  
  // checkpoint;

  std::vector<double> masses(ff->size());

  ifstream mifs(opt.massfile.c_str());
  unsigned char symbol;
  double mass;
  while (mifs) {
    mifs >> symbol >> mass;
    if (symbol > 31 && symbol <= 127) {
      if (ff->nch(symbol) >= 0) {
	masses[ff->nch(symbol)] = mass;
      }
    }
  }
  mifs.close();

  // checkpoint;

  // for (int i=0;i<masses.size();i++) {
  //   cerr << i << " " << ff->ch(i) << " " << masses[i] << endl;
  //

  unsigned char neos = ff->nch(opt.eos_char);
  // cerr << (int)neos << endl;
  unsigned char c;
  bool delete_stream = false;
  istream * ifs(&cin);
  if (opt.indfile != "-") {
    ifs = new ifstream(opt.indfile.c_str());
    delete_stream = true;
  }
  unsigned long index;
  while ((*ifs) >> index) {
    // cerr << index << endl;
    ff->fasta_pos(index-1,0);
    mass = 0;
    while ((c=ff->getnch()) != neos) {
      // cerr << ff->ch(c) << " " << (int)c << " " << masses[c] << " " << mass << endl;
      mass += masses[c];
    }
    *opt.out << index << " " << setprecision(10) << mass << endl;
  }
  if (delete_stream) {
    delete ifs;
  }

  return 0;
}

