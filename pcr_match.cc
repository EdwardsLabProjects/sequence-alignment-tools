
#include <unistd.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include <map>
#include "keyword_tree.h"
#include "keyword_tree.t"
#include "shift_and_inexact.h"
#include "shift_and.h"
#include "sts_io.h"
#include "fasta_io.t"
#include "char_io.t"
#include "primer_alignment.h"
#include "util.h"
#include "types.h"
#include "select.h"
#include "select.t"
#include "pattern_match.h"
#include "pattern_alignment.h"
#include "sortedvector.t"

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
  abort();
  return;
}

struct Options {
  bool pattern_file;
  bool sts_pattern_file;
  bool fasta_pattern_file;
  std::string patterns;
  std::string database;
  int nmismatch;
  char eos_char;
  std::string default_alignformat;
  std::string alignformat;
  bool alignments;
  ostream *out;
  bool delout;
  bool verbose;
  bool veryverbose;
  bool memmap;
  int dbind;
  int node;
  bool wc;
  bool tn;
  bool ucdict;
  bool allorient;
  bool rev_comp;
  int maxdist;
  int mindist;
  int deviation;
  int stlen;
  int edlen;
  int fplen;
  int tplen;
  int seedlen;
  bool indels;
  bool betweenlen;
  unsigned long int report_interval;
  Options(int argc, char *argv[]);
  ~Options();
  void usage(const char *msg=NULL);
};

Options::Options(int argc, char *argv[]) : out(&cout) {
  signed char c;
  optarg = NULL;
  allorient = false;
  rev_comp = false;
  pattern_file = false;
  sts_pattern_file = false;
  fasta_pattern_file = false;
  eos_char = '\n';
  nmismatch = 0;
  delout = false;
  default_alignformat = ">%h\\n %>T %>s ... %l ... %<e %<T\\n %>A  %!>s    %!l    %!<e  %<A\\n %>Q %>r%!>s    %!l    %!<e%<r %<Q %a%R\\n";
  alignformat = default_alignformat;
  alignments = true;
  verbose = false;
  veryverbose = false;
  memmap = true;
  dbind = 0;
  wc = false;
  tn = false;
  node =0;
  report_interval=1000;
  mindist=0;
  maxdist=2000;
  deviation=-1;
  stlen =  0;
  edlen = 0;
  tplen = 0;
  fplen = 0;
  seedlen=0;
  indels = true;
  betweenlen = false;
    
  while ((c = getopt(argc, argv, "p:i:o:P:S:F:E:R:k:K:s:e:5:3:x:hrvVubaA:BD:wWN:M:m:d:")) != -1)
    switch (c) {
    case 'p':
      patterns = optarg;
      pattern_file = false;
      break;
    case 'P':
      patterns = optarg;
      pattern_file = true;
      break;
    case 'S':
      patterns = optarg;
      sts_pattern_file = true;
      break;
     case 'F':
      patterns = optarg;
      fasta_pattern_file = true;
      break;
   case 'i':
      database = optarg;
      break;
    case 'o':
      if (delout) {
	delete out;
      }
      out = new ofstream(optarg,(ios::out|ios::app|ios::ate));
      delout = true;
      break;
    case 'k':
      nmismatch = atoi(optarg);
      indels=true;
      break;
    case 'K':
      nmismatch = atoi(optarg);
      indels=false;
      break;
    case '3':
      tplen = atoi(optarg);
      if (optarg[0] == '~') {
	tplen = -atoi(optarg+1);
      } else {
	tplen = atoi(optarg);
      }
      break;
    case '5':
      fplen = atoi(optarg);
      if (optarg[0] == '~') {
	fplen = -atoi(optarg+1);
      } else {
	fplen = atoi(optarg);
      }
      break;
    case 's':
      if (optarg[0] == '~') {
	stlen = -atoi(optarg+1);
      } else {
	stlen = atoi(optarg);
      }
      break;
    case 'e':
      if (optarg[0] == '~') {
	edlen = -atoi(optarg+1);
      } else {
	edlen = atoi(optarg);
      }
      break;
    case 'x':
      seedlen = atoi(optarg);
      break;
    case 'R':
      report_interval = atoi(optarg);
      break;
    case 'A':
      alignformat = optarg;
      alignments = true;
      break;
    case 'w':
      wc = true;
      tn = false;
      break;
    case 'W':
      wc = true;
      tn = true;
      break;
    case 'u':
      ucdict = true;
      break;
    case 'D':
      dbind = atoi(optarg);
      break;
    case 'N':
      node = atoi(optarg);
      break;
    case 'M':
      maxdist = atoi(optarg);
      break;
    case 'd':
      deviation = atoi(optarg);
      break;
    case 'm':
      mindist = atoi(optarg);
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
    case 'b':
      betweenlen = true;
      break; 
    case 'V':
      veryverbose = true;
      verbose = true;
      break; 
    case 'r':
      rev_comp = true;
      break; 
    case 'a':
      allorient = true;
      break; 
    case 'B':
      memmap = false;
      break; 
    case 'h':
    default :
      usage();
    }
  if ((patterns == "" || database == "")&&!verbose) usage();
  if (nmismatch < 0) usage("Number of mismatches (-k) must be at least 0");
  if (dbind < 0 || dbind > 4) usage("Invalid integer for fasta database indexing (-D).");
}

Options::~Options() {
  if (delout) {
    delete out;
  }
}

std::string spaces(FILE_POSITION_TYPE fp) {
  std::string ret=" ";
  while (fp/=10) ret+=' ';
  return ret;
}

std::string spaces(std::string const & instr) {
  std::string ret=instr;
  for (int i=0;i<instr.length();i++) {
    ret[i] = ' ';
  }
  return ret;
}

void Options::usage(const char *message) {
  if (message != NULL && strlen(message) > 0) {
    cerr << message << endl;
    cerr << endl;
  }
  cerr << "Usage: pcr_match [options] \n\n";
  cerr << "Options: \n";
  cerr << "  -i <sequence-database> Input sequence database. Required.\n";
  cerr << "  -p <sequences>         Primer pairs, separated by whitespace.\n";
  cerr << "  -P <sequence-file>     Primer pairs, separated by whitespace.\n";
  cerr << "                         \"-\" indicates standard input.\n";
  cerr << "  -S <sequence-file>     Primer pairs in UniSTS format.\n";
  cerr << "                         \"-\" indicates standard input.\n";
  cerr << "  -F <sequence-file>     Primer pairs in FASTA format.\n";
  cerr << "                         \"-\" indicates standard input.\n";
  cerr << "                         One of -p, -P or -S is required.\n";
  cerr << "  -o <output-file>       Output file. Defaults to standard out.\n";
  cerr << "                         Will append if output-file exists.\n";
  cerr << "  -k <#-edits>           Number of insertions, deletions \n";
  cerr << "                         and substitutions permitted. Defaults to 0.\n";
  cerr << "  -K <#-mismatches>      Number of substitutions (mismatches) permitted.\n";
  cerr << "                         No insertions or deletions allowed. Defaults to 0.\n";
  cerr << "                         At most one of -k and -K may be specified.\n";
  cerr << "  -r                     Reverse reverse complement primer.\n";
  cerr << "                         Note: -r is automatically set for \n";
  cerr << "                         UniSTS format primer-pair input.\n";
  cerr << "  -a                     Output all primer-pair orientations.\n";
  cerr << "  -x <#-chars>           Number of characters for e(x)act primer seed (word size).\n";
  cerr << "  -s <#-chars>           Number of characters from start of pattern\n";
  cerr << "                         (or reverse complement) that must match exactly\n";
  cerr << "                         (inexactly if preceeded by ~).\n";
  cerr << "  -e <#-chars>           Number of characters from end of pattern\n";
  cerr << "                         (or reverse complement) that must match exactly\n";
  cerr << "                         (inexactly if preceeded by ~).\n";
  cerr << "  -5 <#-chars>           Number of characters from 5' end of pattern\n";
  cerr << "                         (or reverse complement) that must match exactly\n";
  cerr << "                         (inexactly if preceeded by ~).\n";
  cerr << "  -3 <#-chars>           Number of characters from 3' end of pattern\n";
  cerr << "                         (or reverse complement) that must match exactly\n";
  cerr << "                         (inexactly if preceeded by ~).\n";
  cerr << "  -u                     Upper-case all primers.\n";
  cerr << "  -w                     Honor IUPAC wildcards in pattern & text,\n";
  cerr << "                         no text N wildcard.\n";
  cerr << "  -W                     Honor IUPAC wildcards in pattern & text\n";
  cerr << "                         with text N wildcard.\n";
  cerr << "  -E <int>               End-of-sequence character. Default is \'\\n\'\n";
  cerr << "  -m <int>               Minimum amplicon length. Default: 0.\n";
  cerr << "  -M <int>               Maximum amplicon length. Default: 2000.\n";
  cerr << "  -d <int>               Deviation from reported amplicon length. \n";
  cerr << "                         STS format primers required. \n";
  cerr << "                         Default: No constraint.\n";
  cerr << "  -b                     Ignore primers in amplicon length computation.\n";                       
  cerr << "                         Default: False.\n";
  cerr << "  -A <format>            Alignment output format.\n";
  cerr << "                         Default: " << default_alignformat << endl;
  cerr << "  -R <int>               Alignment report interval. Default is 2000.\n";
  cerr << "  -N <int>               Data structure for primer index.\n";
  cerr << "                         Default: Auto.\n";
  cerr << "  -B                     Don\'t use memmap for I/O, use buffered I/O instead.\n";
  cerr << "  -D (0|1|2|3|4)         Fasta database indexing and preprocessing.\n";
  cerr << "                         0: Auto, 1: None, 2: Indexed, 3: Normalized,\n";
  cerr << "                         4: Compressed. Default: Auto.\n";
  cerr << "  -v                     Verbose (version & diagnostic) output.\n";
  cerr << "  -h                     Command line option help.\n";
  cerr << "\n";
  exit(1);
}

void alignformat(ostream & os,
		 std::string const & alignformat, 
		 FILE_POSITION_TYPE const & s,
		 FILE_POSITION_TYPE const & s1,
		 FILE_POSITION_TYPE const & e,
		 FILE_POSITION_TYPE const & e1,
		 FILE_POSITION_TYPE const & five,
		 FILE_POSITION_TYPE const & five1,
		 FILE_POSITION_TYPE const & three,
		 FILE_POSITION_TYPE const & three1,
		 FILE_POSITION_TYPE const & S,
		 FILE_POSITION_TYPE const & S1,
		 FILE_POSITION_TYPE const & E,
		 FILE_POSITION_TYPE const & E1,
		 long unsigned int const & i,
		 unsigned int const & d,
		 unsigned int const & d1,
		 std::string const & p,
		 std::string const & p1,
		 sts_entry const & sts,
		 std::string const & patdeff,
		 std::string const & patdefr,
		 std::string const & q,
		 std::string const & q1,
		 std::string const & Q,
		 std::string const & Q1,
		 std::string const & r,
		 std::string const & r1,
		 std::string const & R,
		 std::string const & R1,
		 bool const & ppo,
		 std::string const & t,
		 std::string const & t1,
		 std::string const & T,
		 std::string const & T1,
		 std::string const & A,
		 std::string const & A1,
		 std::string const & h,
		 std::string const & H,
		 long unsigned int const & f,
		 char * const & a,
		 long unsigned int const & ncount) {

  //
  // The following codes get transformed into positional markers for a
  // printf format.
  // 
  // %s    start of alignment from the start of the fasta entry
  // %e    end of alignment from the start of the fasta entry
  // %S    start of alignment in compressed sequence file
  // %E    end of alignment in compressed sequence file
  // %i    index of pattern
  // %I    sts id number
  // %d    edit distance of alignment
  // %p    tandem repeat motif
  // %Q    aligned tandem repeat with insertion characters
  // %t    aligned text
  // %T    aligned text with alignment characters
  // %A    alignment string
  // %h    fasta header of hit sequence
  // %H    first word of fasta header of hit sequence
  // %%    percent character
  //

  int pos=0;
  while (pos<alignformat.length()) {
    if (alignformat[pos] == '%') {
      pos++;
      if (pos < alignformat.length()) {
	bool widthonly=false;
	if (alignformat[pos] == '!') {
	  widthonly = true;
	  pos++;
	}
	int dirn=0;
	if (alignformat[pos] == '>') dirn = 1;	
	if (alignformat[pos] == '<') dirn = -1;
	if (dirn!=0) pos++;
	switch(alignformat[pos]) {
	case 's': 
	  if (dirn > 0) {
	    if (!widthonly) {
	      os << s;
	    } else {
	      os << spaces(s);
	    }
	  } else if (dirn < 0) {
	    if (!widthonly) {
	      os << s1;
	    } else {
	      os << spaces(s1);
	    }
	  }
	  break;
	case 'e':
	  if (dirn > 0) {
	    if (!widthonly) {
	      os << e;
	    } else {
	      os << spaces(e);
	    }
	  } else if (dirn < 0) {
	    if (!widthonly) {
	      os << e1;
	    } else {
	      os << spaces(e1);
	    }
	  } 
	  break;
	case 'l':
	  if (dirn > 0) {
	    os << e-s;
	  } else if (dirn < 0) {
	    os << e1-s1;
	  } else {
	    if (!widthonly) {
	      os << e1-s;
	    } else {
	      os << spaces(e1-s);
	    }
	  }
	  break;
	case 'S': 
	  if (dirn > 0) {
	    os << S;
	  } else if (dirn < 0) {
	    os << S1;
	  }
	  break;
	case 'E':
	  if (dirn > 0) {
	    os << E;
	  } else if (dirn < 0) {
	    os << E1;
	  }
	  break;
	case 'i':
	  os << i;
	  break;
	case 'd':
	  if (dirn > 0) {
	    os << d;
	  } else if (dirn < 0) {
	    os << d1;
	  }
	  break;
	case 'p':
	  if (dirn > 0) {
	    os << p;
	  } else if (dirn < 0) {
	    os << p1;
	  }
	  break;
	case 'P':
	  if (dirn > 0) {
	    os << patdeff;
	  } else if (dirn < 0) {
	    os << patdefr;
	  }
	  break;
	case 'I':
	  os << sts.id();
	  break;
	case 'L':
	  if (sts.sizeub() != sts.sizelb()) {
	    if (dirn > 0) {
	      os << sts.sizelb();
	    } else if (dirn < 0) {
	      os << sts.sizeub();
	    } else {
	      os << sts.sizelb() << "-" << sts.sizeub();
	    }
	  } else {
	    os << sts.sizelb();
	  }
	  break;
	case 'D': 
	  {
	    int amplen = e1-s;
	    int deviance = 0;
	    if (amplen > sts.sizeub()) {
	      deviance = amplen - sts.sizeub();
	    } else if (amplen < sts.sizelb()) {
	      deviance = sts.sizelb() - amplen;
	    }
	    os << deviance;
	  } 
	  break;
	case 'a':
	  os << sts.accession();
	  break;	  
	case 'O':
	  os << sts.species();
	  break;	  
	case '&':
	  os << sts.altacc();
	  break;	  
	case 'X':
	  os << sts.chrom();
	  break;	  
	case 'q':
	  if (dirn > 0) {
	    os << q;
	  } else if (dirn < 0) {
	    os << q1;
	  }
	  break;
	case 'Q':
	  if (dirn > 0) {
	    if (!widthonly) {
	      os << Q;
	    } else {
	      os << spaces(Q);
	    }
	  } else if (dirn < 0) {
	    if (!widthonly) {
	      os << Q1;
	    } else {
	      os << spaces(Q1);
	    }

	  }
	  break;
	case 'r':
	  if (dirn > 0) {
	    os << r;
	  } else if (dirn < 0) {
	    os << r1;
	  } else {
	    if (ppo) {
	      os << "F";
	    } else {
	      os << "R";
	    }

	  }
	  break;
	case 'R':
	  if (dirn > 0) {
	    os << R;
	  } else if (dirn < 0) {
	    os << R1;
	  } else {
	    if (ppo) {
	      os << "";
	    } else {
	      os << " REVERSE-STRAND";
	    }
	  }
	  break;
	case 't':
	  if (dirn > 0) {
	    os << t;
	  } else if (dirn < 0) {
	    os << t1;
	  }
	  break;
	case 'T':
	  if (dirn > 0) {
	    os << T;
	  } else if (dirn < 0) {
	    os << T1;
	  }
	  break;
	case 'A':
	  if (dirn > 0) {
	    if (!widthonly) {
	      os << A;
	    } else {
	      os << spaces(A);
	    }
	  } else if (dirn < 0) {
	    if (!widthonly) {
	      os << A1;
	    } else {
	      os << spaces(A1);
	    }
	  }
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
	case '@':
	  os << a;
	  break;
	case '*':
	  os << ((ppo)?a:reverse_comp(a));
	  break;
	case 'N':
	  os << ncount;
	  break;
	case '%':
	  os << "%";
	  break;
	case '0':
	  os << H << " " 
	     << s+1 << ".." << e1 << '\t' 
	     << sts.id() << '\t';
	  if (sts.accession() != "") {
	    os << '\t' << sts.accession();
	    if (sts.chrom() != "") {
	      os << '\t' << sts.chrom();
	      if (sts.altacc() != "") {
		os << '\t' << sts.altacc();
		if (sts.species() != "") {
		  os << '\t' << sts.species();
		}
	      }
	    }
	  }
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
	  os << endl;
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
#if defined(PROFILE)
  set_profile_signal_handler();
#endif
  // assert(sizeof(FILE_POSITION_TYPE)>=8);

  Options opt(argc,argv);

  if (opt.verbose) {
    ostrstream ss;
    ss << "Release Tag: " << release_tag << ends;
    std::string v(ss.str());
    timestamp(v.c_str());
  }
  std::list<std::string> patterns;
  std::list<std::string>::iterator patit;
  std::list<sts_entry> sts;
  sts_entry null_sts_entry;
  std::list<sts_entry>::iterator stsit;
  std::list<std::string> patdeflines;
  std::list<std::string>::iterator patdefit;

  if (opt.pattern_file) {
    bool delete_stream = false;
    istream * ifs(&cin);
    if (opt.patterns != "-") {
      ifs = new ifstream(opt.patterns.c_str());
      delete_stream = true;
    }
    std::string pattern;
    while ((*ifs) >> pattern) {
      patterns.push_back(pattern);
    }
    if (delete_stream) {
      delete ifs;
    }
  } else if (opt.sts_pattern_file) {
    bool delete_stream = false;
    istream * ifs(&cin);
    if (opt.patterns != "-") {
      ifs = new ifstream(opt.patterns.c_str());
      delete_stream = true;
    }
    sts_entry s;
    while ((*ifs) >> s) {
      if (s.forward_primer() == "") break;
      sts.push_back(s);
      // cerr << s << endl;
      patterns.push_back(s.forward_primer());
      patterns.push_back(s.reverse_primer());
    }
    if (delete_stream) {
      delete ifs;
    }
  } else if (opt.fasta_pattern_file) {
    bool delete_stream = false;
    istream * ifs(&cin);
    if (opt.patterns != "-") {
      ifs = new ifstream(opt.patterns.c_str());
      delete_stream = true;
    }
    fasta_entry f;
    while ((*ifs) >> f) {
      if (f.sequence() == "") break;
      patdeflines.push_back(f.defline());
      // cerr << s << endl;
      patterns.push_back(f.sequence());
    }
    if (delete_stream) {
      delete ifs;
    }
  } else {
    istrstream sis(opt.patterns.c_str());
    std::string pattern;
    while (sis >> pattern) {
      patterns.push_back(pattern);
    }
  }
  if (patterns.empty()) {
    exit(0);
  } else if (patterns.size() % 2 != 0) {
    opt.usage("Odd number of primers!");
  }

  if (opt.ucdict) {
    for (patit=patterns.begin();patit!=patterns.end();++patit) {
      uppercase(*patit);
    }
  }
  
  if (opt.rev_comp || opt.sts_pattern_file) {
    opt.rev_comp = true;
    int i=0;
    for (patit=patterns.begin();patit!=patterns.end();++patit) {
      if (i%2==1) {
	*patit = reverse_comp(*patit);
      }
      i++;
    }
  }

  std::string pattern;
  std::list<std::string>::size_type i=1,j=1,k=1,n=patterns.size();
  unsigned long int N1=2*n;
  std::vector<std::string> patarray(N1+1); /* indices 1..2*n, index 0 wasted */
  std::vector<std::pair<int,int> > patconst(N1+1);
  std::vector<int> patlen(N1+1);
  std::vector<sts_entry> stsarray(sts.size()+1);
  std::vector<std::string> patdefarray(patdeflines.size()+1);
  int min_exact_const=MAXINT;

  for (patit=patterns.begin(),stsit=sts.begin(),patdefit=patdeflines.begin();
       patit!=patterns.end();++patit) {

    int fplen=opt.fplen;
    int tplen=opt.tplen;
    if (i%2==0) {
      fplen = opt.tplen;
      tplen = opt.fplen;
    }
    if (opt.verbose && patterns.size() < 100 || opt.veryverbose) {
      ostrstream ss;
      if (i%2==1) {
	ss << "[" << setw(4) << i << "] Forward primer: " 
	   << setw(3) << (i-1)/2+1 << " > " 
	   << *patit << ends;
      } else {
	ss << "[" << setw(4) << i << "] Reverse primer: " 
	   << setw(3) << (i-1)/2+1 << " > " 
	   << *patit << ends;
      }
      std::string v(ss.str());
      timestamp(v.c_str());
    }
    patarray[i]=*patit;
    patlen[i]=patit->length();
    if (opt.sts_pattern_file && i%2 == 1) {
      stsarray[j] = *stsit;
      j++;
      stsit++;
    }
    if (opt.fasta_pattern_file) {
      patdefarray[k] = *patdefit;
      k++;
      patdefit++;
    }
    if (opt.stlen > 0) {
      patconst[i].first = opt.stlen;
    } else {
      patconst[i].first = 0;
    }
    if (fplen > patconst[i].first) {
      patconst[i].first = fplen;
    }
    if (opt.edlen < 0 && patit->length() + opt.edlen > patconst[i].first) {
      patconst[i].first = patit->length()+ opt.edlen;
    }
    if (tplen < 0 && patit->length() + tplen > patconst[i].first) {
      patconst[i].first = patit->length()+ tplen;
    }    
    if (opt.edlen > 0) {
      patconst[i].second = opt.edlen;
    } else {
      patconst[i].second = 0;
    }
    if (tplen > patconst[i].second) {
      patconst[i].second = tplen;
    }
    if (opt.stlen < 0 && patit->length() + opt.stlen > patconst[i].second) {
      patconst[i].second = patit->length()+ opt.stlen;
    }
    if (fplen < 0 && patit->length() + fplen > patconst[i].second) {
      patconst[i].second = patit->length()+ fplen;
    }    
    patarray[i+n]=reverse_comp(*patit);
    patlen[i+n]=patarray[i+n].length();
    if (opt.verbose && patterns.size() < 100 || opt.veryverbose) {
      ostrstream ss;
      if (i%2==1) {
	ss << "[" << setw(4) << i+n << "] Forward primer: " 
	   << setw(3) << (i-1)/2+1 << " < " 
	   << patarray[i+n] << ends;
      } else {
	ss << "[" << setw(4) << i+n << "] Reverse primer: " 
	   << setw(3) << (i-1)/2+1 << " < " 
	   << patarray[i+n] << ends;
      }
      std::string v(ss.str());
      timestamp(v.c_str());
    }
    if (opt.stlen > 0) {
      patconst[i+n].first = opt.stlen;
    } else {
      patconst[i+n].first = 0;
    }
    if (tplen > patconst[i+n].first) {
      patconst[i+n].first = tplen;
    }
    if (opt.edlen < 0 && patit->length() + opt.edlen > patconst[i+n].first) {
      patconst[i+n].first = patit->length()+ opt.edlen;
    }
    if (fplen < 0 && patit->length() + fplen > patconst[i+n].first) {
      patconst[i+n].first = patit->length()+ fplen;
    }    
    if (opt.edlen > 0) {
      patconst[i+n].second = opt.edlen;
    } else {
      patconst[i+n].second = 0;
    }
    if (fplen > patconst[i+n].second) {
      patconst[i+n].second = fplen;
    }
    if (opt.stlen < 0 && patit->length() + opt.stlen > patconst[i+n].second) {
      patconst[i+n].second = patit->length()+ opt.stlen;
    }
    if (tplen < 0 && patit->length() + tplen > patconst[i+n].second) {
      patconst[i+n].second = patit->length()+ tplen;
    }    
    i++;
  }
  patterns.clear();
  sts.clear();
  patdeflines.clear();

  FastaFile<Lazy_Header_SI>* ff;
  fasta_file_seq_params ffp;
  ffp.upper_case = false;
  ffp.eos_char = opt.eos_char;
  ff = pick_fasta_file<Lazy_Header_SI>(opt.database,opt.dbind,opt.memmap,
				       /*alignments=*/true,ffp,opt.verbose);
  PatternMatch* pm;
  pm = pick_pattern_index(ff,opt.node,opt.nmismatch,&patconst,&patlen,
			  opt.seedlen,opt.wc,opt.tn,opt.indels,false,
			  opt.eos_char,opt.verbose);

  for (int i=1;i<=N1;i++) {
    pm->add_pattern(patarray[i],i,patconst[i].first,patconst[i].second);
  }
  if (opt.verbose) {
    pm->progress_interval(*ff);
  }
  pm->init(*ff);

  if (opt.verbose) {
    timestamp("Scanning sequence database...");
  }

  typedef sortedvector<FILE_POSITION_TYPE,pattern_hit_vector::iterator> pathitpq;
  typedef std::map<unsigned long int,pathitpq> pathitmap;
  typedef std::pair<unsigned long int,pathitpq> pathitmapelt;  

  pattern_hit_vector l(opt.report_interval*2);
  pattern_hit_vector::iterator it,it2;
  std::list<pattern_hit_vector::iterator> l1;
  std::list<pattern_hit_vector::iterator>::iterator it1;

  // checkpoint;
  bool more;
  while ((more=pm->find_patterns(*ff,l,opt.report_interval))||!l.empty()) {
    // checkpoint;
    // cerr << l.size() << endl;
    FILE_POSITION_TYPE oldcharspos;
    oldcharspos = ff->pos();
    
    pathitmap m;
    pathitmap::iterator mit,mit1,mit2;
    pathitpq::iterator pqit,pqit1;

    l.normalize();

    // checkpoint;
    it = l.begin();
    while (it != l.end()) {
      // checkpoint;
      if ((mit=m.find(it->value().first->id()))==m.end()) {
	// checkpoint;
	std::pair<pathitmap::iterator,bool> ret;
	ret = m.insert(pathitmapelt(it->value().first->id(),pathitpq()));
	mit=ret.first;
      } 
      mit->second.push_back(it->key(),it);
      ++it;
    }
    // checkpoint;
    it = l.begin();
    while (it != l.end()) {
      // checkpoint;
      long unsigned int pid=it->value().first->id();
      long unsigned int pos=it->key();
      // cerr << pid << " " << pos << endl;
      long unsigned int pid1=0,pid2=0;
      if (pid <= n && pid%2 == 1) {
	pid1 = pid+1;
	// cerr << "pid1 = " << pid1 << endl;
      } else if (pid > n && (pid-n)%2==0) {
	pid1 = pid-1;
	// cerr << "pid1 = " << pid1 << endl;
      } 
      if (opt.allorient) {
	if (pid <= n) {
	  if (pid%2 == 1) {
	    pid2 = pid+n+1;
	    // cerr << "pid2 = " << pid2 << endl;
	  } else {
	    pid1 = pid-1;
	    pid2 = pid+n-1;
	    // cerr << "pid1 = " << pid1 << endl;
	    // cerr << "pid2 = " << pid2 << endl;
	  }
	} else {
	  if (pid%2 == 0) {
	    pid2 = pid-n-1;
	    // cerr << "pid2 = " << pid2 << endl;
	  } else {
	    pid1 = pid+1;
	    pid2 = pid-n+1;
	    // cerr << "pid1 = " << pid1 << endl;
	    // cerr << "pid2 = " << pid2 << endl;
	  }
	}
      }

      // checkpoint;

      // checkpoint;
      // cerr << pid << " " << pid1 << " " << pid2 << endl;

      l1.erase(l1.begin(),l1.end());
      bool det;
      mit=m.find(pid);
      assert(mit != m.end());
      pqit = mit->second.locate_first_at_least(pos);
      assert(pqit->key() == pos);
      if (pid1>0) {
	mit1=m.find(pid1);
      } else {
	mit1=m.end();
      }
      if (pid2>0) {
	mit2=m.find(pid2);
      } else {
	mit2=m.end();
      }
      // checkpoint;
      long unsigned int pair = (pid-((pid>n)?n:0)+1)/2;
      FILE_POSITION_TYPE stretch_max=opt.maxdist;
      FILE_POSITION_TYPE stretch_min=opt.mindist;
      if (opt.betweenlen) {
	FILE_POSITION_TYPE plen=0;
	if (pid1 != 0) {
	  plen = patlen[pid1];
	} 
	if (pid2 != 0 && patlen[pid2]>plen) {
	  plen = patlen[pid2];
	}
	stretch_max += plen + patlen[pid];
      }
      if (opt.sts_pattern_file && opt.deviation >= 0) {
	if (stretch_max > stsarray[pair].sizeub() + opt.deviation) {
	  stretch_max = stsarray[pair].sizeub() + opt.deviation;
	} 
	if (stretch_min < stsarray[pair].sizelb() - opt.deviation) {
	  stretch_min = stsarray[pair].sizelb() - opt.deviation;
	}
      }
      stretch_max += pos - patlen[pid] + (opt.indels?opt.nmismatch:1);
      stretch_min += pos - patlen[pid] - (opt.indels?opt.nmismatch:1);
      if (oldcharspos < stretch_max && more) {
	det=false;
      } else {
	det=true;
	// checkpoint;
	if (mit1 != m.end()) {
	  try {
	    pqit1 = mit1->second.locate_first_at_least(stretch_min);
	  } 
	  catch (pathitpq::KeyOutOfRange const & e) {
	    // checkpoint;
	    pqit1 = mit1->second.end();
	  }
	  if (pqit1!=mit1->second.end()) {
	    while (pqit1!=mit1->second.end() &&
		   pqit1->value()->key() <= stretch_max) {
	      if (pqit1->value()->key()) {
		l1.push_back(pqit1->value());
	      }
	      ++pqit1;
	    }
	  }
	}
	// checkpoint;
	if (mit2 != m.end()) {
	  try {
	    pqit1 = mit2->second.locate_first_at_least(stretch_min);
	  } 
	  catch (pathitpq::KeyOutOfRange const & e) {
	    pqit1 = mit2->second.end();
	  }
	  if (pqit1!=mit2->second.end()) {
	    while (pqit1!=mit2->second.end() &&
		   pqit1->value()->key() <= stretch_max) {
	      if (pqit1->value()->key()) {
		l1.push_back(pqit1->value());
	      }
	      ++pqit1;
	    }
	  }
	}
      }

      // checkpoint;
      if (!det&&more) {
	++it;
	continue;
      }
      // checkpoint;
      for (it1=l1.begin();it1!=l1.end();++it1) {
	// checkpoint;
	pattern_alignment *pa = 
	  new editdist_alignment(pos,pos,opt.nmismatch,opt.eos_char,
				 opt.wc,opt.tn,opt.indels,false,
				 patconst[pid].first,
				 patconst[pid].second,false);
	long unsigned int pid1 = (*it1)->value().first->id();
	FILE_POSITION_TYPE pos1 = (*it1)->key();
	// cerr << pid1 << " " << pos1 << endl;
	pattern_alignment *pa1 = 
	  new editdist_alignment(pos1,pos1,opt.nmismatch,opt.eos_char,
				 opt.wc,opt.tn,opt.indels,false,
				 patconst[pid1].first,
				 patconst[pid1].second,false);
	// checkpoint;
	// cerr << pid << " " << pos << " " << patarray[pid] << endl;
	// checkpoint;
	// cerr << pid1 << " " << pos1 << " " << patarray[pid1] << endl;
	pa->align(*ff,patarray[pid]);	
	// checkpoint;
	pa1->align(*ff,patarray[pid1]);	
	// checkpoint;
	if (pa->editdist() <= opt.nmismatch && 
	    pa1->editdist() <= opt.nmismatch ) {
	  FILE_POSITION_TYPE spe,spe1;
	  FILE_POSITION_TYPE sps,sps1;
	  FILE_POSITION_TYPE pe,pe1;
	  FILE_POSITION_TYPE ps,ps1;
	  bool rc,rc1;
	  long unsigned int ind,ind1,pind;
	  unsigned int ed,ed1;
	  spe = ff->get_seq_pos(pa->end());
	  spe1 = ff->get_seq_pos(pa1->end());
	  sps = spe-pa->length()+1;
	  sps1 = spe1-pa1->length()+1;
	  pe = pa->end();
	  pe1 = pa1->end();
	  ps = pe-pa->length()+1;
	  ps1 = pe1-pa1->length()+1;
	  rc = (pid>n);
	  rc1 = (pid1>n);
	  ind = pid-(rc?n:0);
	  ind1 = pid1-(rc1?n:0);
	  ed = pa->editdist();
	  ed1 = pa1->editdist();
	  sts_entry & stsref = null_sts_entry;
	  if (ind<ind1) {
	    pind = ind/2+1;
	  } else {
	    pind = ind1/2+1;
	  }
	  if (opt.sts_pattern_file) {
	    stsref = stsarray[pind];
	  }
	  std::string patdef("");
	  if (opt.fasta_pattern_file) {
	    patdef = patdefarray[ind];
	  }
	  std::string patdef1("");
	  if (opt.fasta_pattern_file) {
	    patdef1 = patdefarray[ind1];
	  }
	  if (opt.rev_comp) {
	    if (ind%2 == 0) {
	      rc = !rc;
	    } else if (ind1%2 == 0) {
	      rc1 = !rc1;
	    }
	  }
	  long signed int amplicon_len;
	  if (!opt.betweenlen) {
	    amplicon_len = pe1;
	    amplicon_len -= ps;
	  } else {
	    amplicon_len = ps1;
	    amplicon_len -= pe;	    
	  }
	  if (ff->is_subseq(ps,pe1) && 
	      amplicon_len <= opt.maxdist && amplicon_len >= opt.mindist && 
	      (!opt.sts_pattern_file || opt.deviation < 0 || 
	       ((amplicon_len + opt.deviation >= stsref.sizelb()) &&
		(amplicon_len <= ((int)stsref.sizeub()) + opt.deviation)))) {
	    const Lazy_Header_SI & h = ff->get_header_data(pa->end());
	    char *buffer = new char[amplicon_len+1];
	    long unsigned int ncount=0;
	    int i=0;
	    ff->pos(ps);
	    while (i<amplicon_len) {
	      buffer[i] = ff->getch();
	      if (buffer[i] == 'N' || buffer[i] == 'n') {
		ncount++;
	      }
	      i++;
	    }
	    buffer[i] = '\0';
	    alignformat(*opt.out,
			opt.alignformat,sps,sps1,spe,spe1,
			(rc?spe:sps),(rc1?spe1:sps1),
			(rc?sps:spe),(rc1?sps1:spe1),
			ps,ps1,pe,pe1,pind,ed,ed1,patarray[ind],
			patarray[ind1],
			stsref,patdef,patdef1,patarray[pid],
			patarray[pid1],
			pa->alignment_pattern(patarray[pid]),
			pa1->alignment_pattern(patarray[pid1]),
			(rc?"R":"F"),(rc1?"R":"F"),
			(rc?" REVCOMP":""),(rc1?" REVCOMP":""),
			(ind<ind1),
			pa->matching_text(),pa1->matching_text(),
			pa->alignment_text(),pa1->alignment_text(),
			pa->alignment_string(),pa1->alignment_string(),
			h.header(),h.short_header(),h.index(),buffer,ncount);
	    delete [] buffer;
	    // checkpoint;
	  } 
	  // checkpoint;
	}
	delete pa;
	delete pa1;
      }
      it->key() = 0;
      ++it;
    }
    m.erase(m.begin(),m.end());
    // checkpoint;
    // fix up l...
    // Move everything with valid key to beginning...
    // int i,i1;
    it=l.begin();
    it2=l.begin();
    // i=0;
    // i1=0;
    long unsigned int newsize=0;
    while (it!=l.end()) {
      // cerr << i << " " << i1 << endl;
      // cerr << l.key(it) << " " << l.value(it)->id() << endl;
      if (it->key() != 0) {
	if (it != it2) {
	  it2->key() = it->key();
	  it2->value() = it->value();
	}
	++newsize;
	++it2;
	// i1++;
      }
      ++it;
      // i++;
    }
    // cerr << newsize << endl;
    l.resize(newsize);
    ff->pos(oldcharspos);
    // checkpoint;
  }

  if (opt.verbose) {
    timestamp("Done.");
  }
  return 0;
}
