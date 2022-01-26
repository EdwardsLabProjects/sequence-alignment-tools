
#include <unistd.h>
#include <assert.h>
#include <vector>
#include <algorithm>
#include <iomanip>
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
  bool pattern_file;
  std::string patterns;
  std::string database;
  char eos_char;
  ostream *out;
  bool delout;
  unsigned long int report_interval;
  bool verbose;
  bool veryverbose;
  bool tryptic;
  bool memmap;
  bool translate;
  int frame;
  int mapindex;
  int node;
  int dbind;
  int hashsize;
  int nmismatches;
  int contextlen;
  Options(int argc, char *argv[]);
  ~Options();
  void usage(const char *msg=NULL);
};

Options::Options(int argc, char *argv[]) : out(&cout) {
  signed char c;
  optarg = NULL;
  pattern_file = false;
  eos_char = '\n';
  report_interval = 1000;
  delout = false;
  verbose = false;
  veryverbose = false;
  memmap = true;
  node = 0;
  dbind = 0;
  tryptic = false;
  translate = false;
  frame = 0;
  hashsize = 4;
  nmismatches = 0;
  contextlen = 1;
  mapindex = 0;
    
  while ((c = getopt(argc, argv, "p:i:o:P:E:hvR:BN:D:tT:x:K:C:M:")) != -1)
    switch (c) {
    case 'p':
      patterns = optarg;
      pattern_file = false;
      break;
    case 'P':
      patterns = optarg;
      pattern_file = true;
      break;
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
    case 'R':
      report_interval = atoi(optarg);
      break;
    case 'N':
      node = atoi(optarg);
      break;
    case 'D':
      dbind = atoi(optarg);
      break;
    case 'K':
      nmismatches = atoi(optarg);
      break;
    case 'x':
      hashsize = atoi(optarg);
      break;
    case 'M':
      mapindex = atoi(optarg);
      break;
    case 'C':
      contextlen = atoi(optarg);
      break;
    case 'T':
      translate = true;
      switch (*optarg) {
      case 'A':	frame = 0; break;
      case 'F':	frame = 4; break;
	/* case 'R':	frame = -4; break;
	   case '1':	frame = 1; break;
	   case '2':	frame = 2; break;
	   case '3':	frame = 3; break;
	   case '4':	frame = -1; break;
	   case '5':	frame = -2; break;
	   case '6':	frame = -3; break; */
      }
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
    case 't':
      tryptic = true;
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
  if ((patterns == "" || database == "")&&!verbose) usage("No peptides and/or no sequence database supplied.");
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
  cerr << "Usage: peptide_scan [options] \n\n";
  cerr << "Options: \n";
  cerr << "  -i <sequence-database> Input sequence database. Required.\n";
  cerr << "  -p <sequences>         Peptide sequences, separated by whitespace.\n";
  cerr << "  -P <sequence-file>     Peptide sequences, separated by whitespace.\n";
  cerr << "                         \"-\" indicates standard input.\n";
  cerr << "                         One of -p or -P is required.\n";
  cerr << "  -T (A|F)               Translate DNA sequence.\n";
  cerr << "  -M <int>               Amino-acid symbol map index. 2: I/L; 3: I/L,K/Q.\n";
  cerr << "  -K <int>               Number of permitted DNA substitutions. Default: 0.\n";
  cerr << "  -x <int>               Hash size (in amino-acids). Default: 4.\n";
  cerr << "  -C <int>               Lengh of amino-acid context. Default: 1.\n";
  cerr << "  -o <output-file>       Output file. Defaults to standard out.\n";
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

  Options opt(argc,argv);

  if (opt.verbose) {
    ostrstream ss;
    ss << "Release Tag: " << release_tag << ends;
    std::string v(ss.str());
    timestamp(v.c_str());
  }

  std::list<std::string> patterns;
  std::list<std::string>::iterator patit;

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
  } else {
    istrstream sis(opt.patterns.c_str());
    std::string pattern;
    while (sis >> pattern) {
      patterns.push_back(pattern);
    }
  }

  if (opt.verbose) {
    timestamp("Finished reading peptides");
  }
  
  if (patterns.empty()) {
    exit(0);
  }

  for (patit=patterns.begin();patit!=patterns.end();++patit) {
    uppercase(*patit);
  }

  if (opt.verbose) {
    timestamp("Finished uppercasing peptides");
  }

  FastaFile<Lazy_Header_SI>* ff;
  fasta_file_seq_params ffp;
  ffp.upper_case = false;
  ffp.eos_char = opt.eos_char;
  ffp.eos_start = true;
  ffp.translate = opt.translate;
  ffp.frame = opt.frame;
  ffp.mapindex = opt.mapindex;
  ff = pick_fasta_file<Lazy_Header_SI>(opt.database,opt.dbind,opt.memmap,
				       /*alignments=*/true,ffp,opt.verbose);

  PatternMatch* kt;
  kt = pick_pattern_index(ff,opt.node,opt.nmismatches,0,0,
                          /*seedlen=*/opt.hashsize,/*wc=*/false,/*tn=*/false,
			  /*indels=*/false,/*dna_mut*/true,opt.eos_char,opt.verbose);

  int npatterns = patterns.size();
  for (patit=patterns.begin();patit!=patterns.end();++patit) {
    kt->add_pattern(*patit);
  }
  if (opt.translate && opt.frame <= 0) {
    for (patit=patterns.begin();patit!=patterns.end();++patit) {
      kt->add_pattern(reverse(*patit));
    }
  }
  
  if (opt.verbose) {
    kt->progress_interval(*ff);
  }
  kt->init(*ff);

  if (opt.verbose) {
    timestamp("Finished building peptide index");
  }

  pattern_hit_vector l(opt.report_interval*2);
  pattern_hit_vector::iterator it;
  bool more;
  // checkpoint;
  if (opt.verbose) {
    timestamp("Begin.");
  }
  while ((more=kt->find_patterns(*ff,l,opt.report_interval))||!l.empty()) {
    // checkpoint;
    FILE_POSITION_TYPE oldcharspos;
    oldcharspos = ff->pos();
    it = l.begin();
    while (it != l.end()) {
      // checkpoint;
      long unsigned int pid = it->value().first->id();
      bool rc=false;
      if (pid > npatterns) {
 	pid -= npatterns;
	rc = true;
      }
      
      pattern_alignment *pa=0;

      if (opt.nmismatches == 0) {
	pa = new exact_peptide_alignment(it->key());
	((exact_peptide_alignment*)pa)->set_context(opt.contextlen);
      } else {
	pa = new editdist_peptide_alignment(it->key(),it->key(),
					    opt.nmismatches,opt.eos_char,
					    false,false,false,
					    opt.translate,0,0,false,opt.translate);
	((editdist_peptide_alignment*)pa)->set_context(opt.contextlen);
      }
      std::string pepstr(it->value().first->pattern());
      pa->align(*ff,pepstr); 
      if (rc) {
        pepstr = reverse(pepstr);
      }
      FILE_POSITION_TYPE p = pa->end();
      // cout << p << endl;
      int frame = 0;
      if (opt.translate) {
	FILE_POSITION_TYPE p1=0;
	ff->getbasepos(p,p1,frame);
	p = p1;
      }
      if (pa->value() > opt.nmismatches || (rc && frame < 3) || (!rc && frame > 2)) {
	delete pa;
	++it;
	continue;
      }
      // cout << p << endl;
      FILE_POSITION_TYPE spe = ff->get_seq_pos(p);
      // cout << spe << endl;
      if (opt.translate) {
	frame = spe % 3 + 1;
	spe = spe / 3;	      
      }
      const FILE_POSITION_TYPE sps = spe-pa->length()+1;
      const FILE_POSITION_TYPE pe = pa->end();
      const FILE_POSITION_TYPE ps = pe-pa->length()+1;
      const Lazy_Header_SI & h = ff->get_header_data(p);
      std::string lcontext;
      std::string rcontext;
      if (opt.nmismatches == 0) {
	lcontext = ((exact_peptide_alignment*)pa)->lcontext();
	rcontext = ((exact_peptide_alignment*)pa)->rcontext();
      } else {
	lcontext = ((editdist_peptide_alignment*)pa)->lcontext();
	rcontext = ((editdist_peptide_alignment*)pa)->rcontext();
      }
      if (rc) {
	swap(lcontext,rcontext);
      }
      std::string::size_type pos = lcontext.rfind(opt.eos_char);
      if (pos != string::npos) {
	lcontext = std::string("-") + lcontext.substr(pos+1);
      }
      pos = rcontext.find(opt.eos_char);
      if (pos != string::npos) {
	rcontext = rcontext.substr(0,pos)+'-';
      }
      if (sps >= 0) {
	std::string buffer="-";
	if (opt.translate) {
	  buffer = std::string((pa->length()-1)*3,' ');
	  int i=0;
	  ff->pos(ps);
	  while (i<(pa->length()-1)*3) {
	    buffer[i] = ff->getbasech();
	    i++;
	  }
	  if (rc) {
	    buffer = reverse_comp(buffer);
	  }
	} else {
          buffer=pa->matching_text();
        }
	*opt.out << pid << " " 
		 << sps << " " 
		 << spe << " " 
		 << lcontext << " "
		 << pepstr << " " 
		 << rcontext << " "
	         << ff->get_seq_pos(p)-(pa->length()-1)*(opt.translate?3:1) << " "
	         << ff->get_seq_pos(p) << " "
	         << frame << " "
	         << (rc?"R":"F") << " "
	         << buffer << " "
		 << h.index() << " >" << h.header();
	if (pa->value() > 0) {
	  float delta = 0.0;
	  string at = pa->alignment_string();
	  string mt = pa->matching_text();
	  size_t p= at.find_first_not_of("|");
	  int j = 1;
	  while (p != string::npos) {
	    char to = pepstr[p];
	    char from = mt[p];
	    *opt.out << " /sub" << j << "=" << from << p+1 << "->" << to << "(" << aasubdist(from,to) << ")";
	    delta += (monomolwt(to)-monomolwt(from));
	    p = at.find_first_not_of("|",p+1);
	    j += 1;
	  }
	  char buf1[100];
	  sprintf(buf1,"%.2f",delta);
	  *opt.out << " /delta=" << buf1;
	}
	*opt.out << endl;
      }
      delete pa;
      ++it;
    }
    l.clear();
    ff->pos(oldcharspos);
  }
  if (opt.verbose) {
    timestamp("Done.");
  }
  return 0;
}

