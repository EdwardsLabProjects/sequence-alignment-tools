
#include <unistd.h>
#include <assert.h>
#include <vector>
#include "keyword_tree.h"
#include "shift_and.h"
#include "shift_and_inexact.h"
#include "fasta_io.h"
#include "fasta_io.t"
#include "char_io.t"
#include "pattern_alignment.h"
#include "select.t"

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

void usage(const char *message=NULL) {
  if (message != NULL && strlen(message) > 0) {
    cerr << message << endl;
    cerr << endl;
  }
  cerr << "Usage: inexact_match [options] \n\n";
  cerr << "Options: \n";
  cerr << "  -p <sequences>         Pattern sequences, separated by whitespace.\n";
  cerr << "  -P <sequence-file>     Pattern sequences, in FASTA format.\n";
  cerr << "                         One of -p or -P is required.\n";
  cerr << "  -k <#-mismatches>      Number of permitted character mismatches. Required.\n";
  cerr << "  -i <sequence-database> Input sequence database. Required.\n";
  cerr << "  -E <int>               End-of-sequence character.\n";
  cerr << "  -r                     Match reverse complement too.\n";
  cerr << "  -w                     Honor IUPAC wildcards in pattern & text, no text N wildcard.\n";
  cerr << "  -W                     Honor IUPAC wildcards in pattern & text with text N wildcard.\n";
  cerr << "  -u                     Uppercase pattern sequences.\n";
  cerr << "  -q                     Suppress diagnostic messages. Optional.\n";
  cerr << "  -h                     Command line option help.\n";
  cerr << "\n";
  exit(1);
}

struct Options {
  bool rev_comp;
  bool quiet;
  bool pattern_file;
  std::string patterns;
  std::string database;
  bool ucdict;
  int nmismatch;
  char eos_char;
  bool wc;
  bool tn;
  bool memmap;
  int dbind;
};

void options(int argc, char *argv[], Options & opt) {
  signed char c;
  optarg = NULL;
  opt.quiet = false;
  opt.rev_comp = false;
  opt.ucdict = false;
  opt.nmismatch = -1;
  opt.eos_char = '\n';
  opt.memmap = true;
  opt.dbind = 0;
  opt.wc = false;
  opt.tn = false;
  while ((c = getopt(argc, argv, "p:i:P:hruE:k:wWqD:B")) != -1)
  switch (c) {
  case 'p':
    opt.patterns = optarg;
    opt.pattern_file = false;
    break;
  case 'P':
    opt.patterns = optarg;
    opt.pattern_file = true;
    break;
  case 'i':
    opt.database = optarg;
    break;
  case 'q':
    opt.quiet = true;
    break;
  case 'r':
    opt.rev_comp = true;
    break;
  case 'u':
    opt.ucdict = true;
    break; 
  case 'w':
    opt.wc = true;
    opt.tn = false;
    break;
  case 'W':
    opt.wc = true;
    opt.tn = true;
    break;
  case 'k':
    opt.nmismatch = atoi(optarg);
    break;
  case 'E': {
    int eos_char;
    if (!sscanf(optarg,"%i",&eos_char)) {
      usage("Invalid end-of-sequence specification.\n");
    }
    opt.eos_char = (char)eos_char;
    }
    break;
  case 'B':
    opt.memmap = false;
    break; 
  case 'D':
    opt.dbind = atoi(optarg);
    break;
  case 'h':
  default :
    usage();
  }
  if (opt.patterns == "" || opt.database == "" || opt.nmismatch < 0) usage();
  if (opt.dbind < 0 || opt.dbind> 4) usage("Invalid integer for fasta database indexing (-D).");
}

int main(int argc,char *argv[]) {
  
#ifdef PROFILE
  set_profile_signal_handler();
#endif

  Options opt;
  options(argc,argv,opt);

  std::list<std::string> patterns;
  std::list<std::string>::iterator patit;
  if (opt.pattern_file) {
    ifstream ifs(opt.patterns.c_str());
    std::string pattern;
    while (ifs >> pattern) {
      patterns.push_back(pattern);
    }      
    if (pattern != "") {
      patterns.push_back(pattern);
    }
    ifs.close();
  } else {
    istrstream sis(opt.patterns.c_str());
    std::string pattern;
    while (sis >> pattern) {
      patterns.push_back(pattern);
    }
  }

  if (opt.ucdict) {
    for (patit=patterns.begin();patit!=patterns.end();++patit) {
      uppercase(*patit);
    }
  }
  
  PatternMatch *pm=new shift_and_inexact(opt.nmismatch,
					 opt.eos_char,opt.wc,opt.tn);
  std::string pattern;
  std::list<std::string>::size_type i=1,n=patterns.size();
  std::vector<fasta_entry> pattern_map(n+1); /* indices 1..n, 0 index wasted */
  
  for (patit=patterns.begin();patit!=patterns.end();++patit) {
    pm->add_pattern(*patit,i,0,0);
    pattern_map[i] = *patit;
    cerr << "[" << i << "] " << "Add pattern > " << *patit << endl;
    if (opt.rev_comp) {
      pm->add_pattern(reverse_comp(*patit),n+i,0,0);
      cerr << "[" << n+i << "] " << "Add pattern < " << reverse_comp(*patit) << endl;
    }
    i++;
  }

  
  FastaFile<Header>* chars;
  fasta_file_seq_params ffp;
  ffp.upper_case = opt.ucdict;
  ffp.eos_char = opt.eos_char;
  chars = pick_fasta_file<Header>(opt.database,opt.dbind,opt.memmap,
				  !opt.quiet,ffp,!opt.quiet);
  pm->progress_interval(*chars,1.0);
  pm->init(*chars);

  pattern_hit_vector l(1000);
  pattern_hit_vector::iterator it;
  while (pm->find_patterns(*chars,l,1000)) {
    it = l.begin();
    while (it != l.end()) {
      std::string kw;
      if (it->value().first->id() > (unsigned int) n) {
	kw = reverse_comp(pattern_map[it->value().first->id()-n].sequence());
	} else {
	  kw = pattern_map[it->value().first->id()].sequence();
	}
      if (!opt.quiet) {
	editdist_alignment ea(it->value().first->id(),it->key(),
			      opt.nmismatch,opt.eos_char,opt.wc,opt.tn,true,false,0,0);
	ea.align(*chars,kw);
	cout << '>' << chars->get_header_data(ea.end()).header() << endl;
	std::string const & mt = ea.matching_text();
	int p=0;
	cout << " ";
	for (int i=0; i<ea.size(); i++) {
	  if (ea[i] != alignment_deletion) {
	    cout << mt[p];
	    p++;
	  } else {
	    cout << "-";
	  }
	}
	cout << " " << ea.start() << " " << it->key() << " " << ea.editdist();
	cout << endl;
	cout << " ";
	for (int i=0; i<ea.size(); i++) {
	  if (ea[i] == alignment_equal) {
	    cout << "|";
	  } else if (ea[i] == alignment_wildcard_equal) {
	    cout << "+";
	  } else if (ea[i] == alignment_substitution) {
	    cout << "*";
	  } else if (ea[i] == alignment_insertion) {
	    cout << "^";
	  } else if (ea[i] == alignment_deletion) {
	    cout << "v";
	  } else {
	    assert(0);
	  }
	}
	cout << endl;
	p = 0;
	// cout << kw << endl;
	cout << " ";
	for (int i=0; i<ea.size(); i++) {
	  if (ea[i] != alignment_insertion) {
	    cout << kw[p];
	    p++;
	  } else {
	    cout << "-";
	  }
	}
	cout << " " << it->value().first->id();
	cout << endl;
      } else {
	cout << it->value().first->id() << " " 
	     << kw << " "
	     << it->key()
	     << endl;
      }
      ++it;
    }
    l.clear();
  }
  
  delete pm;

}




