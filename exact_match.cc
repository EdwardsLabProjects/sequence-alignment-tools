
#include <unistd.h>
#include <assert.h>
#include <vector>
#include <list>
#include "keyword_tree.h"
#include "keyword_tree.t"
#include "shift_and.h"
#include "fasta_io.h"
#include "fasta_io.t"
#include "char_io.t"
#include "types.h"

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

void usage(char *message=NULL) {
  if (message != NULL && strlen(message) > 0) {
    cerr << message << endl;
    cerr << endl;
  }
  cerr << "Usage: exact_match [options] \n\n";
  cerr << "Options: \n";
  cerr << "  -p <sequences>         Pattern sequences, separated by whitespace.\n";
  cerr << "  -P <sequence-file>     Pattern sequences, in FASTA format.\n";
  cerr << "                         One of -p or -p is required.\n";
  cerr << "  -i <sequence-database> Input sequence database. Required.\n";
  cerr << "  -r                     Match reverse complement too.\n";
  cerr << "  -u                     Uppercase pattern sequences.\n";
  cerr << "  -k                     Use keyword tree.\n";
  cerr << "  -b                     Use bitvector.\n";
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
  bool keyword_tree;
  bool bitvec;
  bool ucdict;
};

void options(int argc, char *argv[], Options & opt) {
  signed char c;
  optarg = NULL;
  opt.quiet = false;
  opt.rev_comp = false;
  opt.keyword_tree = false;
  opt.bitvec = false;
  opt.ucdict = false;
  while ((c = getopt(argc, argv, "p:ri:P:hrukbq")) != -1)
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
  case 'k':
    opt.keyword_tree = true;
    break;
  case 'b':
    opt.bitvec = true;
    break;
  case 'h':
  default :
    usage();
  }
  if (opt.patterns == "" || opt.database == "") usage();
}

int main(int argc,char *argv[]) {
  
  Options opt;
  options(argc,argv,opt);

  std::list<std::string> patterns;
  std::list<std::string>::iterator patit;

  if (opt.pattern_file) {
    ifstream ifs(opt.patterns.c_str());
    std::string p;
    while (ifs >> p) {
      patterns.push_back(p);
    }
    ifs.close();
  } else {
    istrstream sis(opt.patterns.c_str());
    std::string p;
    while (sis >> p) {
      patterns.push_back(p);
    }
  }

  
  if (opt.ucdict) {
    for(patit=patterns.begin();patit!=patterns.end();++patit) {
      uppercase(*patit);
    }
  }

  PatternMatch * pm=NULL;
  if (opt.keyword_tree) {
    pm = new keyword_tree<ktnode_jtable>;
  } else if (opt.bitvec) {
    pm = new shift_and;
  } else {
    // Some sort of default behaviour or error here...
    pm = new keyword_tree<ktnode_jtable>;
  }
  
  std::string pattern;
  std::list<std::string>::size_type i=1,n=patterns.size();
  std::vector<std::string> pattern_map(n+1); /* indices 1..n, 0 index wasted */
  for(patit=patterns.begin();patit!=patterns.end();++patit) {
    if (!opt.quiet) {
      cerr << "Adding pattern:" << endl;
      cerr << *patit << endl;
    }
    
    pm->add_pattern(*patit,i);
    pattern_map[i] = *patit;
    if (opt.rev_comp) {
      if (!opt.quiet) {
	cerr << "Adding reverse complement:" << endl;
	cerr << reverse_comp(*patit) << endl;
      }
      pm->add_pattern(reverse_comp(*patit),n+i);
    }
    i++;
  }

  fasta_file_seq_params ffp;
  IndexedFastaFile<Normalized<MapFileChars> ,Header> chars(opt.database,!opt.quiet,ffp);
  pm->init(chars);

  pattern_hit_vector l;
  while (pm->find_patterns(chars,l,100)) {
    FILE_POSITION_TYPE oldcharspos;
    oldcharspos = chars.pos();
    l.normalize();
    pattern_hit_vector::iterator it;
    it = l.begin();
    while (it != l.end()) {
      std::string kw;
      if (it->value().first->id() > (unsigned int) n) {
	kw = pattern_map[it->value().first->id()-n];
      } else {
	kw = pattern_map[it->value().first->id()];
      }
      if (!opt.quiet) {
	exact_alignment ea(it->value().first->id());
	ea.align(chars,kw);
	Header const & h = chars.get_header_data(ea.end());
	FILE_POSITION_TYPE seq_end = chars.get_seq_pos(ea.end());
	cout << '>' << h.header() << endl;
	cout << it->value().first->id() << " " 
	     << kw << " "
	     << ea.matching_text() << " "
	     << seq_end - ea.matching_text().length() << " "
	     << seq_end 
	     << endl;
      } else {
	cout << it->key() << " "
	     << kw << " ";
	if (it->value().first->id() > (unsigned int) n) {
	  cout << "REV";
	}
	cout << endl;
      }
      ++it;
    }
    l.clear();
    chars.pos(oldcharspos);
  }
  delete pm;
}




