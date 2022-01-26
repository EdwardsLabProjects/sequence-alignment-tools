
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
  std::string mertable;
  std::string seqdb;
  ostream *out;
  bool delout;
  bool verbose;
  int k;
  bool nmer;
  bool binary;
  Options(int argc, char *argv[]);
  ~Options();
  void usage(char *msg=NULL);
};

Options::Options(int argc, char *argv[]) : out(&cout) {
  signed char c;
  optarg = NULL;
  verbose = false;
  k = 1;
  nmer = false;
  binary = false;
    
  while ((c = getopt(argc, argv, "m:s:o:hvk:nb")) != -1)
    switch (c) {
    case 'm':
      mertable = optarg;
      break;
    case 's':
      seqdb = optarg;
      break;
    case 'o':
      if (delout) {
	delete out;
      }
      out = new ofstream(optarg);
      delout = true;
      break;
    case 'k':
      k = atoi(optarg);
      break;
    case 'b':
      binary = true;
      break;
    case 'v':
      verbose = true;
      break; 
    case 'n':
      nmer = true;
      break; 
    case 'h':
    default :
      usage();
    }
  if (mertable == "" || seqdb == "") usage();
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
  cerr << "Usage: kmer_annotate [options] \n\n";
  cerr << "Options: \n";
  cerr << "  -m <mer-table>         Input mer-table. Required.\n";
  cerr << "  -s <seqdb>             Sequence database to annotate. Required.\n";
  cerr << "  -o <output-file>       Output file. Defaults to standard out.\n";
  cerr << "  -k <mer-size>          k-mer size.\n";
  cerr << "  -n                     Include Ns in k-mers.\n";
  cerr << "  -v                     Verbose (version & diagnostic) output.\n";
  cerr << "  -h                     Command line option help.\n";
  cerr << "\n";
  exit(1);
}

long unsigned int hash(const char *mer, 
		       FILE_POSITION_TYPE start, FILE_POSITION_TYPE end, 
		       int base, int charmap[]) {
  long unsigned int h=0;
  const char *mer1 = mer+start;
  const long unsigned int len=end-start;
  for (FILE_POSITION_TYPE i=0;i<len;i++) {
    int m = charmap[*mer1];
    assert(m >= 0 && m < base);
    h = h*base+m;
    mer1++;
  }
  return h;
}

long unsigned int hashrc(const char *mer, 
			 FILE_POSITION_TYPE start, FILE_POSITION_TYPE end, 
			 int base, int rccharmap[]) {
  long unsigned int h=0;
  const char *mer1 = mer+end-1;
  const long unsigned int len=end-start;
  for (FILE_POSITION_TYPE i=0;i<len;i++) {
    int m = rccharmap[*mer1];
    assert(m >= 0 && m < base);
    h = h*base+m;
    mer1--;
  }
  return h;
}

char *unhash(long unsigned int h, int len, int base, 
	     char invcharmap[], char *mer) {
  for (int i=len-1;i>=0;i--) {
    mer[i] = invcharmap[h%base];
    h/=base;
  }
  return mer;
}

class mypair : public std::pair<long unsigned int, long unsigned int> {
public:
  mypair() {};
  mypair(long unsigned int f, long unsigned int s) 
    : std::pair<long unsigned int, long unsigned int>(f,s) {};
  mypair(std::pair<long unsigned int, long unsigned int> const & a) 
    : std::pair<long unsigned int, long unsigned int>(a) {};
};

bool operator<(const mypair & a, const mypair & b) {
  if (a.second < b.second) return true;
  if (a.second > b.second) return false;
  if (a.first < b.first) return true;
  if (a.first > b.first) return false;
  return false;
};

int main(int argc,char *argv[]) {

  set_new_handler(newfail_handler);
  
  Options opt(argc,argv);
  
  if (opt.verbose) {
    ostrstream ss;
    ss << "Release Tag: " << release_tag << ends;
    std::string v(ss.str());
    timestamp(v.c_str());
  }

  long unsigned int h;
  std::map<long unsigned int,unsigned int> counts;
  std::map<long unsigned int,unsigned int>::iterator ctit;
  
  int charmap[256];
  int rccharmap[256];
  for (int i=0;i<256;i++) {
    charmap[i] = -1;
    rccharmap[i] = -1;
  }
  int nchar;
  char *invcharmap;
  if (opt.nmer) {
    charmap['A'] = 0;  rccharmap['T'] = 0;
    charmap['C'] = 1;  rccharmap['G'] = 1;
    charmap['G'] = 2;  rccharmap['C'] = 2;
    charmap['T'] = 3;  rccharmap['A'] = 3;
    charmap['N'] = 4;  rccharmap['N'] = 4;
    nchar= 5; 
    invcharmap = new char[nchar];
    invcharmap[0] = 'A';
    invcharmap[1] = 'C';
    invcharmap[2] = 'G';
    invcharmap[3] = 'T';
    invcharmap[4] = 'N';
  } else {
    charmap['A'] = 0;  rccharmap['T'] = 0;
    charmap['C'] = 1;  rccharmap['G'] = 1;
    charmap['G'] = 2;  rccharmap['C'] = 2;
    charmap['T'] = 3;  rccharmap['A'] = 3;
    nchar= 4; 
    invcharmap = new char[nchar];
    invcharmap[0] = 'A';
    invcharmap[1] = 'C';
    invcharmap[2] = 'G';
    invcharmap[3] = 'T';
  }

  if (!opt.binary) {
    
    ifstream ifs(opt.mertable.c_str());
    std::string mer;
    long unsigned int count;
    // cerr << opt.mertable.c_str() << endl;
    while (ifs) {
      ifs >> mer >> count;
      // cerr << mer << count << endl;
      if (mer.length() != 0) {
	if (count > 0) {
	  h = hash(mer.c_str(),0,opt.k,nchar,charmap);
	  counts[h] = count;
	  // cerr << h << " " << count << endl;
	}
      }
    }
    ifs.close();

  } else {
    
    long unsigned int b1;
    unsigned int b2;

    ifstream ifs(opt.mertable.c_str());
    ifs.read((char *)(&b1),sizeof(long unsigned int));
    ifs.read((char *)(&b2),sizeof(unsigned int));
    while (!ifs.fail()) {
      counts[b1] = b2;
      // cerr << buffer[0] << " " << buffer[1] << endl;
      ifs.read((char *)(&b1),sizeof(long unsigned int));
      ifs.read((char *)(&b2),sizeof(unsigned int));
    }

  }

  char *buffer1 = new char[opt.k+1];
  buffer1[opt.k] = '\0';

  ifstream ifs1(opt.seqdb.c_str());
  fasta_entry fe;
  while (ifs1) {
    ifs1 >> fe;
    if (fe.sequence() != "") {
      *opt.out << fe.defline() << '\t';
      const char * cstr = fe.sequence().c_str();
      for (FILE_POSITION_TYPE k=opt.k;k<=fe.sequence().length();k++) {
	long unsigned int count=0;
	// cerr << fe.sequence().substr(k-opt.k,opt.k) << " ";
	h = hash(cstr,k-opt.k,k,nchar,charmap);
	// cerr << h << " " << unhash(h,opt.k,nchar,invcharmap,buffer1) 
	// << " " << counts[h] << " ";
	count += counts[h];
	h = hashrc(cstr,k-opt.k,k,nchar,rccharmap);
	// cerr << h << " " << unhash(h,opt.k,nchar,invcharmap,buffer1) 
	// << " " << counts[h] << " ";
	count += counts[h];
	*opt.out << count << " ";
	// cerr << endl;
      }
      *opt.out << endl;
    }
  }
  ifs1.close();

  checkpoint;
  
  return 0;
}

