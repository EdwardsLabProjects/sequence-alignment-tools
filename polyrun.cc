
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
  int dbind;
  char eos_char;
  std::string alignformat;
  bool outputn;
  bool outputanynonacgt;
  bool verbose;
  unsigned long int l;
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
  l = 20;
  outputn = false;
  outputanynonacgt = false;
  alignformat = ">%h\n %s %e %t x %l\n";
    
  while ((c = getopt(argc, argv, "i:o:E:hBD:vl:nNA:")) != -1)
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
    case 'l':
      l = atoi(optarg);
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
    case 'n':
      outputn = true;
      break; 
    case 'N':
      outputanynonacgt = true;
      outputn = true;
      break; 
    case 'A':
      alignformat = optarg;
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
  cerr << "Usage: polyrun [options] \n\n";
  cerr << "Options: \n";
  cerr << "  -i <sequence-database> Input sequence database. Required.\n";
  cerr << "  -o <output-file>       Output file. Defaults to standard out.\n";
  cerr << "  -l <length>            Minimum length poly-nucleotide run.\n";
  cerr << "  -n                     Output \'N\' runs instead.\n";
  cerr << "  -N                     Output non-ACGT runs instead.\n";
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

void alignformat(ostream & os,
		 std::string const & alignformat, 
		 FILE_POSITION_TYPE const & s,
		 FILE_POSITION_TYPE const & e,
		 FILE_POSITION_TYPE const & S,
		 FILE_POSITION_TYPE const & E,
		 char const & t,
		 std::string const & h,
		 std::string const & H,
		 long unsigned int const & f) {

  //
  // The following codes get transformed into positional markers for a
  // printf format.
  // 
  // %s    start of alignment from the start of the fasta entry
  // %e    end of alignment from the start of the fasta entry
  // %S    start of alignment in compressed sequence file
  // %E    end of alignment in compressed sequence file
  // %t    aligned text with alignment characters
  // %h    fasta header of hit sequence
  // %H    first word of fasta header of hit sequence
  // %f    index of header
  // %%    percent character
  //

  int pos=0;
  while (pos<alignformat.length()) {
    if (alignformat[pos] == '%') {
      pos++;
      if (pos < alignformat.length()) {
	switch(alignformat[pos]) {
	case 's': 
	  os << s;
	  break;
	case 'e':
	  os << e;
	  break;
	case 'l':
	  os << e-s;
	  break;
	case 'S': 
	  os << S;
	  break;
	case 'E':
	  os << E;
	  break;
	case 't':
	  os << t;
	  break;
	case 'h':
	  os << h;
	  break;
	case 'H':
	  os << H;
	  break;
	case 'f':
	  os << f;
	  break;
	case '%':
	  os << "%";
	  break;
	default:
	  os << alignformat[pos];
	}
      } else {
	os << "%";
      }
    } else if (alignformat[pos] == '\\') {
      pos++;
      if (pos < alignformat.length()) {
	switch(alignformat[pos]) {
	case 'n': 
	  os << '\n';
	  break;
	case 't': 
	  os << '\t';
	  break;
	case '\\': 
	  os << '\\';
	  break;
	default:
	  os << alignformat[pos];
	}
      } else {
	os << '\\';
      }
    } else {
      os << alignformat[pos];
    }
    pos++;
  }
}

int main(int argc,char *argv[]) {

  set_new_handler(newfail_handler);
  // assert(sizeof(FILE_POSITION_TYPE)>=8);
  
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

  char ch;
  char eos=ff->nch(opt.eos_char);
  signed char lastchar=-1;
  char charn=ff->nch('N');
  char chara=ff->nch('A');
  char charc=ff->nch('C');
  char charg=ff->nch('T');
  char chart=ff->nch('G');
  FILE_POSITION_TYPE polylen=0;
  // checkpoint;
  // cerr << polylen << " " << opt.l << endl;
  while (!ff->eof()) {
    ch = ff->getnch();
    if (opt.outputanynonacgt && ch != chara && ch != charc && 
	                         ch != charg && ch != chart && ch != eos) {
      ch = charn;
    } 
    // cerr << "ch: " << (int)ch << " lastchar " << (int) lastchar << endl;
    if (ch == lastchar) {
      polylen++;
    } else {
      if (polylen >= opt.l && lastchar != eos && 
	  ((opt.outputn && lastchar == charn) || 
	   (!opt.outputn && lastchar != charn))) {
	// checkpoint;
	// cerr << ff->pos() << " " << polylen << endl;
	const FILE_POSITION_TYPE pe = ff->pos()-1;
	const FILE_POSITION_TYPE ps = pe-polylen;
	const FILE_POSITION_TYPE spe = ff->get_seq_pos(pe);
	const FILE_POSITION_TYPE sps = spe-polylen;
	const Header_SI & h = ff->get_header_data(pe);
	alignformat(*opt.out,opt.alignformat,sps,spe,ps,pe,ff->ch(lastchar),
		    h.header(),h.short_header(),h.index());
      }
      lastchar = ch;
      polylen=1;
    }
  }
  
  delete ff;
  
  return 0;
}

