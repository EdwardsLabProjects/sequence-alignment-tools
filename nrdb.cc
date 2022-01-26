
#include <unistd.h>
#include <assert.h>
#include <vector>
#include "keyword_tree.h"
#include "keyword_tree.t"
#include "hash_table.h"
#include "shift_and.h"
#include "shift_and_inexact.h"
#include "fasta_io.h"
#include "fasta_io.t"
#include "char_io.t"
#include "primer_alignment.h"
#include "util.h"
#include "types.h"
#include "select.h"
#include "select.t"
#include "sts_io.h"

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
  bool ucdict;
  int seedlen;
  char eos_char;
  ostream *out;
  bool delout;
  unsigned long int report_interval;
  bool verbose;
  bool memmap;
  int dbind;
  int contained;
  bool bareout;
  bool noshort;
  Options(int argc, char *argv[]);
  ~Options();
  void usage(char msg[]=NULL);
};

Options::Options(int argc, char *argv[]) : out(&cout) {
  signed char c;
  optarg = NULL;
  ucdict = false;
  seedlen = 6;
  eos_char = '\n';
  report_interval = 1000;
  delout = false;
  verbose = false;
  memmap = true;
  contained = false;
  bareout = false;
  noshort = false;
  dbind = 0;
    
  while ((c = getopt(argc, argv, "i:o:x:SCbE:huvR:BD:")) != -1) {
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
    case 'x':
      seedlen = atoi(optarg);
      break;
    case 'u':
      ucdict = true;
      break; 
    case 'b':
      bareout = true;
      break; 
    case 'S':
      noshort = true;
      break; 
    case 'R':
      report_interval = atoi(optarg);
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
    case 'C':
      contained = true;
      break; 
    case 'B':
      memmap = false;
      break; 
    case 'h':
    default :
      usage();
    }
  }
  
  if (database == ""&&!verbose) usage();
  if (seedlen == 0) usage();
  if (dbind < 0 || dbind> 4) usage("Invalid integer for fasta database indexing (-D).");
}

Options::~Options() {
  if (delout) {
    delete out;
  }
}

void Options::usage(char message[]) {
  if (message != NULL && strlen(message) > 0) {
    cerr << message << endl;
    cerr << endl;
  }
  cerr << "Usage: nrdb [options] \n\n";
  cerr << "Options: \n";
  cerr << "  -i <sequence-database> Input sequence database. Required.\n";
  cerr << "  -o <output-file>       Output file. Defaults to standard out.\n";
  cerr << "  -b                     Bare, sequence only, output format. Default: false.\n";
  cerr << "  -x <#-chars>           Number of characters for e(x)act seed (word size).\n";
  cerr << "                         Default is 6.\n";
  cerr << "  -S                     Do not suppress short ( < seedlen ) sequences. Default: false\n";
  cerr << "  -C                     Suppress contained sequences too.\n";
  cerr << "  -E <int>               End-of-sequence character. Default is \'\\n\'\n";
  cerr << "  -R <int>               Alignment report interval. Default is 1000.\n";
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
#if defined(PROFILE)
  set_profile_signal_handler();
#endif
  // assert(sizeof(FILE_POSITION_TYPE)>=8);
  // assert(sizeof(bigword)>=8);
  
  Options opt(argc,argv);
  
  if (opt.verbose) {
    ostrstream ss;
    ss << "Release Tag: " << release_tag << ends;
    std::string v(ss.str());
    timestamp(v.c_str());
  }

  if (opt.database == "") {
    exit(1);
  }
  
  // checkpoint;
  FastaFile<Lazy_Header_SI>* ffd;
  FastaFile<Lazy_Header_SI>* ffq;
  fasta_file_seq_params ffp;
  ffp.upper_case = opt.ucdict;
  ffp.eos_char = opt.eos_char;
  ffd = pick_fasta_file<Lazy_Header_SI>(opt.database,opt.dbind,opt.memmap,
					true,ffp,opt.verbose);
  ffq = pick_fasta_file<Lazy_Header_SI>(opt.database,opt.dbind,opt.memmap,
					true,ffp,opt.verbose);
  PatternMatch* pm;
  pm = new keyword_tree<ktnode_list>;

  int i=0;
  char* buffer = new char[opt.seedlen+1];
  buffer[opt.seedlen] = 0;
  while (ffq->fasta_pos(i,0)) {
    int j=0;
    for (j=0;j<opt.seedlen;j++) {
      buffer[j] = ffq->getch();
      if (buffer[j] == opt.eos_char) {
	break;
      }
    }
    buffer[j] = 0;
    if (j < opt.seedlen) {
      if (opt.verbose) {
        cerr << "Short pattern " << i+1 << " : " << buffer << endl;
      }
      if (!opt.noshort) {
        pm->add_pattern(std::string(buffer),i+1);
      }
    } else {
      pm->add_pattern(std::string(buffer),i+1);
    }
    i++;
  }
  delete [] buffer;

  // checkpoint;

  std::vector<std::list<long unsigned int> > contains;
  contains.resize(i);
  std::vector<bool> contained;
  contained.resize(i);
  fill(contained.begin(),contained.end(),false);

  pm->init(*ffd);

  pattern_hit_vector l(opt.report_interval*2);
  pattern_hit_vector::iterator it;
  // checkpoint;
  while (pm->find_patterns(*ffd,l,opt.report_interval)) {
    // checkpoint;
    FILE_POSITION_TYPE oldcharspos;
    oldcharspos = ffd->pos();
    it = l.begin();
    char qch,dch;
    char neos=ffq->nch(opt.eos_char);
    while (it != l.end()) {
      long unsigned int feq = it->value()->id();      
      feq--;
      FILE_POSITION_TYPE pos = it->key();
      long unsigned int fed = ffd->get_header_data(pos).index();
      fed--;
      if (feq != fed) {
	assert(ffq->fasta_pos(feq,it->value()->pattern().length()));
	ffd->pos(pos);
	qch = ffq->getnch();
	dch = ffd->getnch();
	while (qch == dch && 
	       qch != neos && 
	       dch != neos) {
	  qch = ffq->getnch();
	  dch = ffd->getnch();
	}
	if (qch == neos) {
	  if (dch == neos && ffd->get_seq_pos(pos) == it->value()->pattern().length()) {
	    // cerr << fed << " equals " << feq << endl;
	    if (feq > fed) {
	      contains[fed].push_back(feq);
	      contained[feq] = true;
	    }
	  } else if (opt.contained) {
	    // cerr << fed << " contains " << feq << endl;
	    contains[fed].push_back(feq);
	    contained[feq] = true;
	  }
	} 
      }
      ++it;
    }
    l.clear();
    ffd->pos(oldcharspos);
  }
  // checkpoint;
  std::vector<bool> output;
  output.resize(contains.size());
  
  fill(output.begin(),output.end(),true);
  std::list<long unsigned int> stck;
  // checkpoint;
  for (long unsigned int i=0; i < contains.size(); i++) {
    if (output[i] && !contained[i]) {
      output[i] = false;
      if (!opt.bareout) {
	ffq->fasta_pos(i,0);
	FILE_POSITION_TYPE pos = ffq->pos();
	const Lazy_Header_SI & h = ffq->get_header_data(pos);
	std::string header="";
	// cout << "Output: " << i << endl;
	stck.clear();
	stck.push_back(i);
	while (!stck.empty()) {
	  long unsigned int j = stck.front();
	  stck.pop_front();
	  std::list<long unsigned int>::iterator it;
	  it = contains[j].begin();
	  while (it != contains[j].end()) {
	    if (output[*it]) {
	      stck.push_back(*it);
	    }
	    ++it;
	  }
	  output[j] = false;
	  ffq->fasta_pos(j,0);
	  FILE_POSITION_TYPE pos = ffq->pos();
	  const Lazy_Header_SI & h = ffq->get_header_data(pos);
	  if (header != "") {
	    header+=';';
	  }
	  header+=h.header();
	}
	*opt.out << ">" << header;
      }
      ffq->fasta_pos(i,0);
      char ch;
      FILE_POSITION_TYPE p = 0;
      while ((ch=ffq->getch())!=opt.eos_char) {
	if (!opt.bareout && p%60==0) *opt.out << endl;
	*opt.out << ch;	
	p++;
      }
      *opt.out << endl;
    }
  }
  return 0;
}

