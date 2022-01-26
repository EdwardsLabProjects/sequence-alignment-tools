
#include <unistd.h>
#include <assert.h>
#include <vector>
#include "keyword_tree.h"
#include "keyword_tree.t"
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
#include "math.h"

#if defined(PRIMER3TM)
extern "C" {
#include "oligotm.h"
}
#endif

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
  bool rev_comp;
  bool pattern_file;
  bool fasta_pattern_file;
  bool sts_pattern_file;
  std::string patterns;
  std::string database;
  bool ucdict;
  int tplen;
  int fplen;
  int stlen;
  int edlen;
  int seedlen;
  char eos_char;
  int nmismatch;
  int maxcount;
  std::string default_alignformat;
  std::string alignformat;
  bool alignments;
  std::string default_countformat;
  std::string countformat;
  bool counts;
  ostream *out;
  bool delout;
  unsigned long int report_interval;
  bool verbose;
  bool veryverbose;
  bool memmap;
  int node;
  int dbind;
  bool dbindex;
  bool wc;
  bool tn;
  bool indels;
  bool aggregate;
  bool dna_mutations;
  bool translate;
  Options(int argc, char *argv[]);
  ~Options();
  void usage(const char *msg=NULL);
};

Options::Options(int argc, char *argv[]) : out(&cout) {
  signed char c;
  optarg = NULL;
  pattern_file = false;
  fasta_pattern_file = false;
  sts_pattern_file = false;
  rev_comp = false;
  ucdict = false;
  tplen = 0;
  fplen = 0;
  stlen = 0;
  edlen = 0;
  seedlen = 0;
  eos_char = '\n';
  nmismatch = 0;
  maxcount = 0;
  report_interval = 1000;
  delout = false;
  default_alignformat = ">%h\\n %T %s %e %d\\n %A\\n %Q %i%R\\n";
  alignformat = default_alignformat;
  alignments = true;
  default_countformat = "%i %r %q %c%+ ( %C )\\n";
  countformat = default_countformat;
  counts = false;
  verbose = false;
  veryverbose = false;
  memmap = true;
  node = 0;
  dbind = 0;
  dbindex = true;
  wc = false;
  tn = false;
  indels = true;
  aggregate = false;
  dna_mutations = false;
  translate = false;
    
  while ((c = getopt(argc, argv, "p:i:o:P:F:S:M:k:K:s:e:3:5:x:E:hrucavA:C:R:BN:D:IwWT")) != -1)
    switch (c) {
    case 'p':
      patterns = optarg;
      pattern_file = false;
      break;
    case 'P':
      patterns = optarg;
      pattern_file = true;
      break;
    case 'F':
      patterns = optarg;
      fasta_pattern_file = true;
      break;
    case 'S':
      patterns = optarg;
      sts_pattern_file = true;
      rev_comp = true;
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
    case 'k':
      if (optarg[0] == '.') {
	nmismatch = atoi(optarg+1);
	dna_mutations = true;
      } else {
	nmismatch = atoi(optarg);
      }
      indels = true;
      break;
    case 'K':
      if (optarg[0] == '.') {
	nmismatch = atoi(optarg+1);
	dna_mutations = true;
      } else {
	nmismatch = atoi(optarg);
      }
      indels = false;
      break;
    case 'r':
      rev_comp = true;
      break;
    case 'c':
      counts = true;
      alignments = false;
      break;
    case 'M':
      maxcount = atoi(optarg);
      break;
    case 'x':
      seedlen = atoi(optarg);
      break;
    case 'A':
      if (strlen(optarg)>0)
	alignformat = optarg;
      alignments = true;
      break;
    case 'C':
      if (strlen(optarg)>0)
	countformat = optarg;
      counts = true;
      break;
    case 'u':
      ucdict = true;
      break; 
    case 'a':
      aggregate = true;
      break; 
    case 'T':
      translate = true;
      break; 
    case 'w':
      wc = true;
      tn = false;
      break;
    case 'W':
      wc = true;
      tn = true;
      break;
    case 'R':
      report_interval = atoi(optarg);
      break;
    case 'N':
      node = atoi(optarg);
      break;
    case 'D':
      dbind = atoi(optarg);
      break;
    case 'E': {
      int ec;
      if (!sscanf(optarg,"%i",&ec)) {
	usage("Invalid end-of-sequence specification.\n");
      }
      eos_char = (char)ec;
    }
    break;
    case 'v':
      verbose = true;
      break; 
    case 'I':
      dbindex = false;
      break; 
    case 'V':
      verbose = true;
      veryverbose = true;
      break; 
    case 'B':
      memmap = false;
      break; 
    case 'h':
    default :
      usage();
    }
  if ((patterns == "" || database == "")&&!verbose) usage("No primers and/or no sequence database supplied.");
  if (nmismatch < 0) usage("Number of mismatches (-k) must be >= 0.");
  if (maxcount > 0 && !counts) usage("Can''t use maxcount (-M) without counts (-c or -C).");
  if (aggregate && !counts) usage("Can''t use aggregate (-a) without counts (-c or -C).");
  if (dbind < 0 || dbind> 4) usage("Invalid integer for fasta database indexing (-D).");
  if (dna_mutations && nmismatch > 3) usage("Can''t model > 3 DNA mutations");
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
  cerr << "Usage: primer_match [options] \n\n";
  cerr << "Options: \n";
  cerr << "  -i <sequence-database> Input sequence database. Required.\n";
  cerr << "  -p <sequences>         Primer sequences, separated by whitespace.\n";
  cerr << "  -P <sequence-file>     Primer sequences, separated by whitespace.\n";
  cerr << "                         \"-\" indicates standard input.\n";
  cerr << "  -F <sequence-file>     Primer sequences in FASTA format.\n";
  cerr << "                         \"-\" indicates standard input.\n";
  cerr << "  -S <sequence-file>     Primer sequences in UniSTS format.\n";
  cerr << "                         \"-\" indicates standard input.\n";
  cerr << "                         One of -p, -P, -F, or -S is required.\n";
  cerr << "                         -S automatically sets the -r option.\n";
  cerr << "  -o <output-file>       Output file. Defaults to standard out.\n";
  cerr << "                         Will append if output-file exists.\n";
  cerr << "  -k <#-mismatches>      Number of insertions, deletions \n";
  cerr << "                         and substitutions permitted. Defaults to 0.\n";
  cerr << "  -K <#-mismatches>      Number of substitutions (mismatches) permitted.\n";
  cerr << "                         No insertions or deletions allowed. Defaults to 0.\n";
  cerr << "                         At most one of -k and -K may be specified.\n";
  cerr << "  -r                     Match reverse complement of primers too.\n";
  cerr << "  -x <#-chars>           Number of characters for e(x)act seed (word size).\n";
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
  cerr << "  -u                     Uppercase pattern sequences.\n";
  cerr << "  -w                     Honor IUPAC wildcards in pattern & text,\n";
  cerr << "                         no text N wildcard.\n";
  cerr << "  -W                     Honor IUPAC wildcards in pattern & text\n";
  cerr << "                         with text N wildcard.\n";
  cerr << "  -E <int>               End-of-sequence character. Default is \'\\n\'\n";
  cerr << "  -c                     Output counts (only).\n";
  cerr << "  -a                     Output the aggregate of forward & \n";
  cerr << "                         reverse complement counts.\n";
  cerr << "  -M <max-count>         Maximum number of occurances to count.\n";
  cerr << "  -A <format>            Alignment output format.\n";
  cerr << "                         Default: \"" << default_alignformat << "\".\n";
  cerr << "  -C <format>            Counts output format.\n";
  cerr << "                         Default: \"" << default_countformat << "\".\n";
  cerr << "  -R <int>               Alignment report interval. Default is 1000.\n";
  cerr << "  -N <int>               Data structure for primer index.\n";
  cerr << "                         Default: Auto.\n";
  cerr << "  -B                     Don\'t use memmap for I/O, use buffered I/O instead.\n";
  cerr << "  -D (0|1|2|3|4)         Fasta database indexing and preprocessing.\n";
  cerr << "                         0: Auto, 1: None, 2: Indexed, 3: Normalized,\n";
  cerr << "                         4: Compressed. Default: Auto.\n";
  cerr << "  -I                     Do not load fasta database index.\n";
  cerr << "  -v                     Verbose (version & diagnostic) output.\n";
  cerr << "  -h                     Command line option help.\n";
  cerr << "\n";
  exit(1);
}

void alignformat(ostream & os,
		 std::string const & alignformat, 
		 FILE_POSITION_TYPE const & s,
		 FILE_POSITION_TYPE const & e,
		 FILE_POSITION_TYPE const & five,
		 FILE_POSITION_TYPE const & three,
		 FILE_POSITION_TYPE const & S,
		 FILE_POSITION_TYPE const & E,
		 long unsigned int const & i,
		 unsigned int const & d,
		 std::string const & p,
		 std::string const & P,
		 std::string const & q,
		 std::string const & Q,
		 std::string const & r,
		 std::string const & R,
		 std::string const & t,
		 std::string const & T,
		 std::string const & A,
		 std::string const & h,
		 std::string const & H,
		 long unsigned int const & f,
		 sts_entry const & sts,
		 int frame=-1,std::string buffer="") {
  
  //
  // The following codes get transformed into positional markers for a
  // printf format.
  // 
  // %s    start of alignment from the start of the fasta entry
  // %e    end of alignment from the start of the fasta entry
  // %5    position of 5' end of the alignment from the start of fasta entry
  // %3    position of 3' end of the alignment from the start of fasta entry
  // %S    start of alignment in compressed sequence file
  // %E    end of alignment in compressed sequence file
  // %i    index of pattern
  // %d    edit distance of alignment
  // %p    pattern
  // %q    pattern or rev comp of pattern depending on which hit
  // %Q    aligned pattern or rev comp with insertion characters
  // %r    "F" or "R" depending on whether pattern or rev comp was hit
  // %R    "" or "REVCOMP" depending on whether pattern or rev comp was hit
  // %t    aligned text
  // %T    aligned text with alignment characters
  // %A    alignment string
  // %h    fasta header of hit sequence
  // %H    first work of fasta header of hit sequence
  // %%    percent character
  // 
  // %=    The usual default format, wrapped to approx 50 characters.

  int pos=0;
  bool sc=false;
  unsigned int ins=0;
  unsigned int del=0;
  unsigned int sub=0;
  unsigned int wcm=0;
  unsigned int mat=0;
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
	case '5': 
	  os << five;
	  break;
	case '3':
	  os << three;
	  break;
	case 'S': 
	  os << S;
	  break;
	case 'E':
	  os << E;
	  break;
	case 'i':
	  os << i;
	  break;
	case 'd':
	  os << d;
	  break;
	case 'D':
	  os << p.length()-(s-e);
	  break;
	case 'M': 
	  {
	    double monomw1 = 0.0;
	    double monomw2 = 0.0;
	    for (unsigned int i0=0;i0<p.length();i0++) {
	      monomw1 += monomolwt(p[i0]);
	    }
	    for (unsigned int i0=0;i0<q.length();i0++) {
	      monomw2 += monomolwt(t[i0]);
	    }
	    os << floor((monomw1-monomw2)*100)/100;
	  }
	  break;
	case 'p':
	  os << p;
	  break;
	case 'P':
	  os << P;
	  break;
	case 'q':
	  os << q;
	  break;
	case 'Q':
	  os << Q;
	  break;
	case 'r':
	  os << r;
	  break;
	case 'R':
	  os << R;
	  break;
	case 't':
	  os << t;
	  break;
	case 'T':
	  os << T;
	  break;
	case 'U':
	  os << ((r=="R")?reverse_comp(t):t);
	  break;
	case 'A':
	  os << A;
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
	case 'I':
	  os << sts.id();
	  break;
	case 'L':
	  if (sts.sizeub() != sts.sizelb()) {
	    os << sts.sizelb() << "-" << sts.sizeub();
	  } else {
	    os << sts.sizelb();
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
	case 'F':
	  os << frame;
	  break;
	case 'n':
	  os << buffer;
	  break;
#ifdef PRIMER3TM
	case 'm':
	case 'G':
	  char *dna;
	  char v;
	  v = alignformat[pos];
	  pos++;
	  if (alignformat[pos] == 'p') {
	    dna = strdup(p.c_str());
	  } else if (alignformat[pos] == 'q') {
	    dna = strdup(q.c_str());
	  } else if (alignformat[pos] == 't') {
	    dna = strdup(t.c_str());
	  } else if (alignformat[pos] == 'U') {
	    dna = strdup(((r=="R")?reverse_comp(t):t).c_str());
	  } else {
	    dna = strdup(t.c_str());
	    pos--;
	  }
	  int i,j;
	  i=j=0;
	  while (dna[j] != '\0') {
	    if (index("ACGT",dna[j])) {
	      if (i < j) {
		dna[i] = dna[j];
	      }
	      i += 1;
            }
	    j += 1;
	  }
	  dna[i] = '\0';
	  char buf[100];
	  if (v == 'm') { 
	    sprintf(buf,"%.2f",oligotm(dna,50.0,50.0,0.0,0.0,TM_METHOD_SANTALUCIA,SALT_CORRECTION_SANTALUCIA));
	  } else if (v == 'G') {
	    sprintf(buf,"%.2f",oligodg(dna,TM_METHOD_SANTALUCIA));
	  }
	  os << buf;
	  free(dna);
	  break;
#endif
	case '%':
	  os << "%";
	  break;
	case '|':
	  {
	    if (!sc) {
	      for (unsigned int i0=0;i0<A.length();i0++) {
		switch (A[i0]) {
		case '|': mat++; break;
		case '^': del++; break;
		case 'v': ins++; break;
		case '*': sub++; break;
		case '+': wcm++; break;
		}
	      }
	      sc = true;
	    }
	    os << mat;
	  }
	  break;
	case '^':
	  {
	    if (!sc) {
	      for (unsigned int i0=0;i0<A.length();i0++) {
		switch (A[i0]) {
		case '|': mat++; break;
		case '^': del++; break;
		case 'v': ins++; break;
		case '*': sub++; break;
		case '+': wcm++; break;
		}
	      }
	      sc = true;
	    }
	    os << del;
	  }
	  break;
	case 'v':
	  {
	    if (!sc) {
	      for (unsigned int i0=0;i0<A.length();i0++) {
		switch (A[i0]) {
		case '|': mat++; break;
		case '^': del++; break;
		case 'v': ins++; break;
		case '*': sub++; break;
		case '+': wcm++; break;
		}
	      }
	      sc = true;
	    }
	    os << ins;
	  }
	  break;
	case '*':
	  {
	    if (!sc) {
	      for (unsigned int i0=0;i0<A.length();i0++) {
		switch (A[i0]) {
		case '|': mat++; break;
		case '^': del++; break;
		case 'v': ins++; break;
		case '*': sub++; break;
		case '+': wcm++; break;
		}
	      }
	      sc = true;
	    }
	    os << sub;
	  }
	  break;
	case '+':
	  {
	    if (!sc) {
	      for (unsigned int i0=0;i0<A.length();i0++) {
		switch (A[i0]) {
		case '|': mat++; break;
		case '^': del++; break;
		case 'v': ins++; break;
		case '*': sub++; break;
		case '+': wcm++; break;
		}
	      }
	      sc = true;
	    }
	    os << wcm;
	  }
	  break;
	case '=': 
	  {
	    unsigned int len0=T.length();
	    const unsigned int width0=50;
	    unsigned int textchars_start=0;
	    for (unsigned int i0=0;i0<len0;i0+=width0) {
	      unsigned int nchars=width0;
	      if (i0+nchars > len0) {
		nchars = len0-i0;
	      }
	      unsigned int textchars_end=textchars_start+nchars;
	      unsigned int editcount0=nchars;
	      for (unsigned int j0=0;j0<nchars;j0++) {
		if (A[i0+j0] == '|' || A[i0+j0] == '+') editcount0--;
		if (A[i0+j0] == 'v') textchars_end--;
	      }
	      os << " " << T.substr(i0,width0) 
		 << " " << textchars_start
		 << " " << textchars_end
		 << " " << editcount0
		 << "\n" 
		 << " " << A.substr(i0,width0)
		 << "\n" 
		 << " " << Q.substr(i0,width0)
		 << " " << i << R 
		 << "\n";
	      if (len0-i0 > width0) {
		os << endl;
	      }
	      textchars_start = textchars_end;
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

void countformat(ostream & os,
		 std::string const & format, 
		 long unsigned int const & i,
		 std::string const & p,
		 std::string const & P,
		 std::string const & q,
		 std::string const & r,
		 std::string const & R,
		 long unsigned int const & c,
		 long unsigned int * const & C,
		 unsigned int const & k,
		 bool const & gtmax,
		 sts_entry const & sts) {
  //
  // The following codes get transformed into positional markers for a
  // printf format.
  // 
  // %i    index of pattern
  // %p    pattern
  // %q    pattern or rev comp of pattern depending on which hit
  // %r    "F" or "R" depending on whether pattern or rev comp was hit
  // %R    "" or "REVCOMP" depending on whether pattern or rev comp was hit
  // %c    count of pattern or rev comp
  // %C    whitespace separated vector of counts by editdist
  // %+    greater than max count marker "+"
  // %%    percent character
  //

  int pos=0;
  while (pos<format.length()) {
    if (format[pos] == '%') {
      pos++;
      if (pos < format.length()) {
	switch(format[pos]) {
	case 'i':
	  os << i;
	  break;
	case 'p':
	  os << p;
	  break;
	case 'P':
	  os << P;
	  break;
	case 'q':
	  os << q;
	  break;
	case 'r':
	  os << r;
	  break;
	case 'R':
	  os << R;
	  break;
	case 'c':
	  os << c;
	  break;
	case 'C':
	  for (unsigned int i=0;i<k;i++) {
	    os << C[i] << " ";
	  }
	  os << C[k];
	  break;
	case '+':
	  if (gtmax) {
	    os << "+";
	  }
	  break;
	case '%':
	  os << "%";
	  break;
	case 'I':
	  os << sts.id();
	  break;
	case 'L':
	  if (sts.sizeub() != sts.sizelb()) {
	    os << sts.sizelb() << "-" << sts.sizeub();
	  } else {
	    os << sts.sizelb();
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
	default:
	  os << format[pos];
	}
      } else {
	os << "%";
      }
    } else if (format[pos] == '\\') {
      pos++;
      if (pos < format.length()) {
	switch(format[pos]) {
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
	  os << format[pos];
	}
      } else {
	os << '\\';
      }
    } else {
      os << format[pos];
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
  // assert(sizeof(bigword)>=8);

  Options opt(argc,argv);

  if (opt.verbose) {
    ostrstream ss;
    ss << "Release Tag: " << release_tag << ends;
    std::string v(ss.str());
    timestamp(v.c_str());
  }

  std::list<std::string> patterns;
  std::list<std::string> patdeflines;
  std::list<sts_entry> sts;
  std::list<std::string>::iterator patit;
  std::list<std::string>::iterator patdefit;
  sts_entry null_sts_entry;
  std::list<sts_entry>::iterator stsit;

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
      patterns.push_back(f.sequence());
      patdeflines.push_back(f.defline());
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
      patterns.push_back(s.forward_primer());
      patterns.push_back(s.reverse_primer());
      sts.push_back(s);
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
  }

  if (opt.verbose) {
    timestamp("Read primers");
  }

  if (opt.ucdict) {
    for (patit=patterns.begin();patit!=patterns.end();++patit) {
      uppercase(*patit);
    }
    if (opt.verbose) {
      timestamp("Uppercase primers");
    }
  }

  std::string pattern;
  std::list<std::string>::size_type i=1,n=patterns.size();
  unsigned long int N1=((!opt.rev_comp&&!opt.translate)?1:2)*n;
  std::vector<std::string> patarray(N1+1); /* indices 1..N1, index 0 wasted */
  std::vector<std::pair<int,int> > patconst(N1+1);
  std::vector<int> patlen(N1+1);
  std::vector<std::string> patdefarray;
  std::vector<sts_entry> stsarray;
  if (opt.fasta_pattern_file) {
    patdefarray.resize(n+1);
    patdefit=patdeflines.begin();
  }
  if (opt.sts_pattern_file) {
    stsarray.resize(n/2+1);
    stsit=sts.begin();
  }

  /* This is only used if we are counting hits rather than outputing them */
  std::vector<long unsigned int> patcount;
  std::vector<bool> maxpatcount;
  if (opt.counts) {
    patcount.resize(N1*(opt.nmismatch+1));
    if (opt.maxcount>0) {
      maxpatcount.resize(N1+1);
    }
  }

#define index_patcount(i,k) (((i)-1)*(opt.nmismatch+1)+(k))

  for (patit=patterns.begin();patit!=patterns.end();++patit) {
    if (opt.verbose && patterns.size() < 100 || opt.veryverbose) {
      ostrstream ss;
      ss << "Pattern " << setw(3) << i << " > " 
	 << *patit << ends;
      std::string v(ss.str());
      timestamp(v.c_str());
    }
    patarray[i]=*patit;
    patlen[i]=patit->length();
    if (opt.fasta_pattern_file) {
      patdefarray[i]=*patdefit;
      patdefit++;
    }
    if (opt.sts_pattern_file && i%2 == 1) {
      stsarray[(i+1)/2]=*stsit;
      stsit++;
    }
    if (opt.stlen > 0) {
      patconst[i].first = opt.stlen;
    } else {
      patconst[i].first = 0;
    }
    if (opt.fplen > patconst[i].first) {
      patconst[i].first = opt.fplen;
    }
    if (opt.edlen < 0 && patit->length() + opt.edlen > patconst[i].first) {
      patconst[i].first = patit->length()+ opt.edlen;
    }
    if (opt.tplen < 0 && patit->length() + opt.tplen > patconst[i].first) {
      patconst[i].first = patit->length()+ opt.tplen;
    }    
    if (opt.edlen > 0) {
      patconst[i].second = opt.edlen;
    } else {
      patconst[i].second = 0;
    }
    if (opt.tplen > patconst[i].second) {
      patconst[i].second = opt.tplen;
    }
    if (opt.stlen < 0 && patit->length() + opt.stlen > patconst[i].second) {
      patconst[i].second = patit->length()+ opt.stlen;
    }
    if (opt.fplen < 0 && patit->length() + opt.fplen > patconst[i].second) {
      patconst[i].second = patit->length()+ opt.fplen;
    }    
    if (opt.counts) {
      for (int k=0;k<=opt.nmismatch;k++) {
	patcount[index_patcount(i,k)]=0;
      }
      if (opt.maxcount>0) {
	maxpatcount[i] = false;
      }
    }
    if (opt.rev_comp || opt.translate) {
      if (opt.translate) {
	patarray[i+n]=reverse(*patit);	
      } else {
	patarray[i+n]=reverse_comp(*patit);
      }
      patlen[i+n]=patarray[i+n].length();
      if (opt.stlen > 0) {
	patconst[i+n].first = opt.stlen;
      } else {
	patconst[i+n].first = 0;
      }
      if (opt.tplen > patconst[i+n].first) {
	patconst[i+n].first = opt.tplen;
      }
      if (opt.edlen < 0 && patit->length() + opt.edlen > patconst[i+n].first) {
	patconst[i+n].first = patit->length()+ opt.edlen;
      }
      if (opt.fplen < 0 && patit->length() + opt.fplen > patconst[i+n].first) {
	patconst[i+n].first = patit->length()+ opt.fplen;
      }    
      if (opt.edlen > 0) {
	patconst[i+n].second = opt.edlen;
      } else {
	patconst[i+n].second = 0;
      }
      if (opt.fplen > patconst[i+n].second) {
	patconst[i+n].second = opt.fplen;
      }
      if (opt.stlen < 0 && patit->length() + opt.stlen > patconst[i+n].second) {
	patconst[i+n].second = patit->length()+ opt.stlen;
      }
      if (opt.tplen < 0 && patit->length() + opt.tplen > patconst[i+n].second) {
	patconst[i+n].second = patit->length()+ opt.tplen;
      }    
      
      if (opt.verbose && patterns.size() < 100 || opt.veryverbose) {
	ostrstream ss;
	ss << "Pattern " << setw(3) << i << " < " 
	   << patarray[i+n] << ends;
	std::string v(ss.str());
	timestamp(v.c_str());
      }
      
      if (opt.counts) {
	for (int k=0;k<=opt.nmismatch;k++) {
	  patcount[index_patcount(i+n,k)]=0;
	}
	if (opt.maxcount>0) {
	  maxpatcount[i+n] = false;
	}
      }
    }
    i++;
  }

  if (opt.verbose) {
    timestamp("Put primers in an array");
  }

  patterns.clear();
  patdeflines.clear();
  sts.clear();

  FastaFile<Lazy_Header_SI>* ff;
  fasta_file_seq_params ffp;
  ffp.upper_case = opt.ucdict;
  ffp.eos_char = opt.eos_char;
  ffp.check_params = opt.dbindex;
  ffp.translate = opt.translate;
  ff = pick_fasta_file<Lazy_Header_SI>(opt.database,opt.dbind,opt.memmap,
				       opt.alignments&&opt.dbindex,
				       ffp,opt.verbose);
  
  PatternMatch* kt;
  kt = pick_pattern_index(ff,opt.node,opt.nmismatch,&patconst,&patlen,
                          opt.seedlen,opt.wc,opt.tn,opt.indels,opt.dna_mutations,
			  opt.eos_char,opt.verbose);

  for (int i=1;i<=N1;i++) {
    kt->add_pattern(patarray[i],i,patconst[i].first,patconst[i].second);
  }

  if (opt.verbose) {
    kt->progress_interval(*ff,1.0);
  }
  kt->init(*ff);

  pattern_hit_vector l(opt.report_interval*2);
  pattern_hit_vector::iterator it;
  bool more;
  // checkpoint;
  while ((more=kt->find_patterns(*ff,l,opt.report_interval))||!l.empty()) {
    // checkpoint;
    FILE_POSITION_TYPE oldcharspos;
    oldcharspos = ff->pos();
    it = l.begin();
    while (it != l.end()) {
      // checkpoint;
      long unsigned int pid = it->value().first->id();
      // checkpoint;
      // checkpoint;
      // cerr << pid << endl;
      if (pid &&
	  (opt.maxcount <= 0 || 
	   !maxpatcount[pid])) {
	pattern_alignment *pa(0);
	// checkpoint;
	// cerr << pid << endl;
	if (opt.nmismatch == 0) { 
	  if (!opt.wc) {
	    pa = new exact_alignment(it->key());
	  } else {
	    pa = new exact_wc_alignment(it->key(),opt.tn);	      
	  }
	} else {
	  // checkpoint;
	  pa = new editdist_alignment(it->key(),it->key(),
				      opt.nmismatch,opt.eos_char,
				      opt.wc,opt.tn,opt.indels,
				      opt.dna_mutations,
				      patconst[pid].first,
				      patconst[pid].second,false);
	}
	// checkpoint;
	pa->align(*ff,patarray[pid]);
	// checkpoint;
	if (pa->editdist() <= (unsigned int) opt.nmismatch) {
	  // checkpoint;
	  // cerr << pid << endl;
	  if (opt.alignments) {
	    // checkpoint;
	    FILE_POSITION_TYPE p = pa->end();
	    // cout << p << endl;
	    int frame = 0;
	    if (opt.translate) {
	      FILE_POSITION_TYPE p1=0;
	      ff->getbasepos(p,p1,frame);
	      p = p1;
	    }
	    // cout << p << endl;
	    FILE_POSITION_TYPE spe = ff->get_seq_pos(p);
	    // cout << spe << endl;
	    if (opt.translate) {
	      frame = spe % 3 + 1;
	      spe = spe / 3;	      
	    }
	    // cout << spe << endl;
	    FILE_POSITION_TYPE sps = spe-pa->length()+1;
	    FILE_POSITION_TYPE pe = pa->end();
	    FILE_POSITION_TYPE ps = pe-pa->length()+1;
	    const bool rc = (pid>n);
	    const long unsigned int ind = pid-(rc?n:0);
	    const unsigned int ed = pa->editdist();
	    const std::string pat = patarray[pid];
	    const Lazy_Header_SI & h = ff->get_header_data(p);
	    // checkpoint;
	    std::string patdef("");
	    if (opt.fasta_pattern_file) {
	      patdef = patdefarray[ind];
	    }
	    sts_entry & stsref = null_sts_entry;
	    if (opt.sts_pattern_file) {
	      stsref = stsarray[(ind+1)/2];
	    }
	    // checkpoint;
	    if (!opt.translate) {
	      alignformat(*opt.out,
			  opt.alignformat,sps,spe,(rc?spe:sps),
			  (rc?sps:spe),ps,pe,ind,ed,patarray[ind],
			  patdef,pat,pa->alignment_pattern(pat),
			  (rc?"R":"F"),(rc?" REVCOMP":""),
			  pa->matching_text(),pa->alignment_text(),
			  pa->alignment_string(),
			  h.header(),h.short_header(),h.index(),stsref);
	    } else {
	      std::string buffer((pa->length()-1)*3,' ');
	      int i=0;
	      ff->pos(ps);
	      while (i<(pa->length()-1)*3) {
		buffer[i] = ff->getbasech();
		i++;
	      }
	      if (!rc) {
		alignformat(*opt.out,
			    opt.alignformat,sps,spe,(rc?spe:sps),
			    (rc?sps:spe),ps,pe,ind,ed,patarray[ind],
			    patdef,pat,pa->alignment_pattern(pat),
			    (rc?"R":"F"),(rc?" REVCOMP":""),
			    pa->matching_text(),pa->alignment_text(),
			    pa->alignment_string(),
			    h.header(),h.short_header(),h.index(),stsref,
			    frame,buffer);
	      } else {
		frame = -frame;
		alignformat(*opt.out,
			    opt.alignformat,sps,spe,(rc?spe:sps),
			    (rc?sps:spe),ps,pe,ind,ed,patarray[ind],
			    patdef,reverse(pat),
			    reverse(pa->alignment_pattern(pat)),
			    "R"," REVSTRAND",
			    reverse(pa->matching_text()),
			    reverse(pa->alignment_text()),
			    reverse(pa->alignment_string()),
			    h.header(),h.short_header(),h.index(),stsref,
			    frame,reverse_comp(buffer));
	      }
	    }
	    // checkpoint;
	  }
	  if (opt.counts) {
	    patcount[index_patcount(pid,pa->editdist())]++;
	    if (opt.maxcount > 0) {
	      long unsigned int count=0;
	      for (int k=0;k<=opt.nmismatch;k++) {
		count += patcount[index_patcount(pid,k)];
	      }
	      if (count >= (unsigned int) opt.maxcount) {
		maxpatcount[pid] = true;
	      }
	    }
	  }
	} else {
	  timestamp("Bogus hit returned to primer_match main()");
	  if (opt.alignments) {
	    const Lazy_Header_SI & h = ff->get_header_data(it->key());
	    cerr << "Problem sequence is near:\n";
	    cerr << ">" << h.header() << endl;
	  } else {
	    cerr << "Approximate absolute sequence position:\n";
	    cerr << " " << it->key() << endl;
	  }
	  cerr << "Problem primer:\n";
	  cerr << " " << patarray[pid] << endl;
	  exit(1);
	}
	delete pa;
      }
      ++it;
    }
    l.clear();
    ff->pos(oldcharspos);
  }
  long unsigned int *counts = new long unsigned int [opt.nmismatch+1];
  if (opt.counts) {
    for (unsigned int i=1;i<=n;i++) {
      unsigned long int total=0;
      for (int k=0;k<=opt.nmismatch;k++) {
	counts[k]=patcount[index_patcount(i,k)];
	total += patcount[index_patcount(i,k)];
      }
      bool gtmax=false;
      if (opt.maxcount>0) {
	gtmax = maxpatcount[i];
      }
      std::string patdef("");
      if (opt.fasta_pattern_file) {
	patdef = patdefarray[i];
      }
      sts_entry & stsref = null_sts_entry;
      if (opt.sts_pattern_file) {
	stsref = stsarray[(i+1)/2];
      }
      if (!opt.aggregate) {
	countformat(*opt.out,opt.countformat,
		    i,patarray[i],patdef,
		    patarray[i],"F","",
		    total,counts,opt.nmismatch,gtmax,
		    stsref);
      }
      if (opt.rev_comp || opt.translate) {
	if (!opt.aggregate) {
	  total = 0;
	  for (int k=0;k<=opt.nmismatch;k++) {
	    counts[k] = 0;
	  }
	  gtmax=false;
	}
	for (int k=0;k<=opt.nmismatch;k++) {
	  counts[k]+=patcount[index_patcount(i+n,k)];
	  total += patcount[index_patcount(i+n,k)];
	}
	if (opt.maxcount>0) {
	  gtmax = gtmax || maxpatcount[i+n];
	}
	if (!opt.aggregate) {
	  countformat(*opt.out,opt.countformat,
		      i,patarray[i],patdef,
		      patarray[i+n],"R"," REVCOMP",
		      total,counts,opt.nmismatch,gtmax,
		      stsref);
	}
      }
      if (opt.aggregate) {
	countformat(*opt.out,opt.countformat,
		    i,patarray[i],patdef,
		    "","","",
		    total,counts,
		    opt.nmismatch,gtmax,
		    stsref);
      }
    }
  }
  // delete [] counts;
  if (opt.verbose) {
    timestamp("Done.");
  }
  return 0;
}

