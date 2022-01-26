

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
  std::string patterns;
  std::string database;
  char eos_char;
  ostream *out;
  bool delout;
  bool verbose;
  bool memmap;
  int dbind;
  int mersize;
  int u;
  bool noindex;
  bool indels;
  bool rc;
  int node1;
  int node2;
  int nmismatch;
  int blocksize;
  bool delpat;
  bool uniq;
  bool exonly;
  int hashtablesize;
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
  dbind = 0;
  mersize = 0;
  u = 0;
  noindex = false;
  indels = false;
  nmismatch = 0;
  blocksize = 50000;
  delpat = false;
  node1 = 0;
  node2 = 0;
  uniq = false;
  exonly = false;
  hashtablesize = -1;
  while ((c = getopt(argc, argv, "r:i:o:E:hvBID:m:k:K:u:b:l:dRn:N:UXH:")) != -1)
    switch (c) {
    case 'm':
      mersize = atoi(optarg);
      break;
    case 'k':
      nmismatch = atoi(optarg);
      indels = true;
      break;
    case 'K':
      nmismatch = atoi(optarg);
      indels = false;
      break;
    case 'u':
      u = atoi(optarg);
      break;
    case 'n':
      node1 = atoi(optarg);
      break;
    case 'N':
      node2 = atoi(optarg);
      break;
    case 'r':
      patterns = optarg;
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
    case 'l':
      stderr = fopen(optarg,"a");
      setlinebuf(stderr);
      break;
    case 'D':
      dbind = atoi(optarg);
      break;
    case 'b':
      blocksize = atoi(optarg);
      break;
    case 'E': {
      int eos;
      if (!sscanf(optarg,"%i",&eos)) {
	usage("Invalid end-of-sequence specification.\n");
      }
      eos_char = (char)eos;
    }
    break;
    case 'H':
      hashtablesize = atoi(optarg);
      break;
    case 'v':
      verbose = true;
      break; 
    case 'd':
      delpat = true;
      break; 
    case 'U':
      uniq = true;
      exonly = false;
      break; 
    case 'X':
      exonly = true;
      uniq = false;
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
  if (patterns == "" || database == "" || mersize == 0) usage();
  if (dbind < 0 || dbind > 4) usage("Invalid integer for fasta database indexing (-D).");
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
  cerr << "  -r <sequences>         Regular expressions for mers, separated by whitespace.\n";
  cerr << "  -R                     Reverse complement too.\n";
  cerr << "  -k <int>               Edit distance.\n";
  cerr << "  -K <int>               Hamming distance.\n";
  cerr << "  -b <int>               Inexact search pattern batch size. Default: 50000\n";
  cerr << "  -d                     Delete patterns in phase 1 once xmers. Default: Do not delete.\n";
  cerr << "  -U                     Establish Unique xmers only. Default: Count all matches.\n";
  cerr << "  -X                     Establish xmers with eXact matches only. Default: Count all matches.\n";
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

// A merelt can be in one of 3 states...
// 1. Unique
// 2. non-unique, but not xmer
// 3. xmer
//
// If unique, need to keep l and r chars, and single position.
// If non-unique, but not xmer, need to keep l and r chars, and multiple positions.
// If xmer, need to keep xmer index, and all r chars
//

class merelt {
public:
  static long unsigned int xmerind;
  static long unsigned int nxmers;
  static long unsigned int totmers;
  char *str;
  unsigned char u_:1;
  unsigned char x_:1;
  unsigned char nrl_:2;
  unsigned char nrr_:2;
  unsigned char oep_:1;
  union {
    struct {
      char l;
      char r;
      FILE_POSITION_TYPE p0;
    } unx;
    struct {
      char l;
      char r;
      tinylist<FILE_POSITION_TYPE> *pl;
    } nunx;
    struct {
      char l;
      char r;
      long unsigned int xindex;
    } ux;
    struct {
      u_int32_t l:8;
      u_int32_t r:24; // Have alphabet sizes bigger than 16...
      long unsigned int xindex;
    } nux;
  };
  bool x() const {
    return x_;
  }
  void setx() {
    if (u()) {
      char r0 = unx.r;
      ux.r = r0;
      ux.xindex = xmerind+1;
    } else {
      char r0 = nunx.r;
      nux.r = r0;
      nux.xindex = xmerind+1;
    }
    x_ = true;
    u_ = true;
    xmerind++;
    nxmers++;
  }
  void unsetx() {
    x_ = false;
  }
  bool u() const {
    return u_;
  }
  void unsetu() {
    u_ = false;
  }
  void setu() {
    u_ = true;
  }
  void zero_left_char_count() {
    nrl_ = 0;
  }
  void inc_left_char_count() {
    if ( nrl_ < 2 ) nrl_++;
  }
  void zero_right_char_count() {
    nrr_ = 0;
  }
  void inc_right_char_count() {
    if ( nrr_ < 2 ) nrr_++;    
  }
  int get_left_char_count() {
    return nrl_;
  }
  int get_right_char_count() {
    return nrr_;
  }
  void unset_one_exact_pos() {
    oep_ = 0;
  }
  void set_one_exact_pos() {
    oep_ = 1;
  }
  bool one_exact_pos() {
    return (oep_ == 1);
  }
  merelt(char *s,char lin,char rin,char neos,FILE_POSITION_TYPE p) {
    setu();
    unsetx();
    zero_left_char_count();
    if (lin != neos) {
      inc_left_char_count();
    }
    zero_right_char_count();    
    if (rin != neos) {
      inc_right_char_count();
    }
    set_one_exact_pos();
    str = s;
    unx.l = lin; unx.r = rin;
    unx.p0 = p;
    totmers++;
  }
  int strcmp(const char *s,int k) {
    for (int i=0;i<k;i++) {
      if (s[i] != str[i]) {
	return (s[i]-str[i]);
      }
    }
    return 0;
  }
  bool rchar(int j) {
    assert(x());
    if (u()) {
      return j == ux.r;
    } else {
      return ((((long unsigned int)1)<<j)&nux.r);
    }
  }
  void dump(ostream *o,FastaFile<Lazy_Header_SI> & ff, int k) {
    char *buffer1=new char[k+1];
    for (int i=0;i<k;i++) {
      buffer1[i] = ff.ch(str[i]);
    }
    buffer1[k] = 0;
    if (x()) {
      if (u()) {
	(*o) << "  UX: "
	     << "seq:" << buffer1 << " "
	     << "ind:" << ux.xindex << " "
	     << "nlc:" << get_left_char_count() << " "
	     << "nrc:" << get_right_char_count() << " "
	     << "rch:" << ff.ch(ux.r) << endl;
      } else {
	(*o) << " NUX: "
	     << "seq:" << buffer1 << " "
	     << "ind:" << nux.xindex << " " 
	     << "nlc:" << get_left_char_count() << " "
	     << "nrc:" << get_right_char_count() << " "
	     << "rch:" ;
	for (int i=0;i<3*8;i++) {
	  if ((((long unsigned int)1)<<i)&nux.r) {
	    (*o) << ff.ch(i) << " ";
	  }
	}
	(*o) << endl;
      }
    } else {
      if (u()) {
	(*o) << " UNX: "
	     << "seq:" << buffer1 << " "
	     << "pos:" << unx.p0 << " "
	     << "nlc:" << get_left_char_count() << " "
	     << "nrc:" << get_right_char_count() << " "
	     << "lch:" << ff.ch(unx.l) << " " 
	     << "rch:" << ff.ch(unx.r) << endl;	
      } else {
	(*o) << "NUNX: "
	     << "seq:" << buffer1 << " pos:";
	tinylist<FILE_POSITION_TYPE>::iterator it;
	for (it=nunx.pl->begin();it!=nunx.pl->end();++it) {
	  (*o) << *it << " ";
	}
	(*o) << "nlc:" << get_left_char_count() << " "
	     << "nrc:" << get_right_char_count() << " "
	     << "lch:" << ff.ch(unx.l) << " " 
	     << "rch:" << ff.ch(unx.r) << endl;	
      }
    }
  }
  bool update(char lin,char rin,char neos,int phase,bool indel, 
	      FILE_POSITION_TYPE p, 
	      tinylist<FILE_POSITION_TYPE> * & lout,char & rout) {
    if (x()) {
      if (phase == 0) { // record additional rins
	if (u()) {
	  if (lin != neos) {
	    if (ux.l == neos) {
	      ux.l = lin;
	      inc_left_char_count();
	    } else if (lin != ux.l) {
	      inc_left_char_count();
	    }
	  }
	  if (rin != ux.r) { // transition from ux to nux
	    unsetu();
	    char r0 = ux.r;
	    char l0 = ux.l;
	    nux.r  = (((unsigned)1)<<r0);
	    nux.r |= (((unsigned)1)<<rin);
	    nux.l = l0;
	    if (rin != neos) {
	      inc_right_char_count();
	    }
	  }
	} else { // !u() just add addition(?) rin
	  if (lin != neos) {
	    if (nux.l == neos) {
	      nux.l = lin;
	      inc_left_char_count();
	    } else if (lin != nux.l) {
	      inc_left_char_count();
	    }
	  }
	  if (rin != neos) {
	    if (!((((unsigned)1)<<rin) & nux.r)) {
	      inc_right_char_count();
	      nux.r |= (((unsigned)1)<<rin);
	    }
	  }
	}
	lout = 0;
	return true;
      } else { // !phase == 0 ... do nothing...
	return false;
      }
    } else { // !x()
      if (u()) {
	if (phase == 0) {
	  if (lin != neos) {
	    if (unx.l == neos) {
	      unx.l = lin;
	      inc_left_char_count();
	    } else if (lin != unx.l) {
	      inc_left_char_count();
	    }
	  }
	  if (rin != neos) {
	    if (unx.r == neos) {
	      unx.r = rin;
	      inc_right_char_count();
	    } else if (rin != unx.r) {
	      inc_right_char_count();
	    }
	  }
	}
	if (lin != unx.l || rin != unx.r || (phase != 0 && indel)) {
	  // checkpoint;
	  // cerr << unx.r << endl;
	  // transition to xmer
	  lout = new tinylist<FILE_POSITION_TYPE>;
	  lout->push_front(unx.p0);
	  char r0 = unx.r;
	  char l0 = unx.l;
	  rout = r0;
 	  if (rin != r0 && phase == 0) {
	    // checkpoint;
	    // transition to nux
	    setx();
	    nux.r  = (((unsigned)1)<<r0);
	    nux.r |= (((unsigned)1)<<rin);
	    unsetu();
	    nux.l = l0;
	  } else { // rin == unx.r
	    // transition to ux
	    setx();
	    ux.r = r0;
	    ux.l = l0;
	  }
	  return true;
	} else if (phase == 0) { // now we are not unique, but not xmer
	  FILE_POSITION_TYPE p0 = unx.p0;
	  nunx.pl = new tinylist<FILE_POSITION_TYPE>();
	  nunx.pl->push_front(p0);
	  nunx.pl->push_front(p);
	  unsetu();
	  unset_one_exact_pos();
	  return false;
	}
      } else { // !u()
	if (phase == 0) {
	  if (lin != neos) {
	    if (nunx.l == neos) {
	      nunx.l = lin;
	      inc_left_char_count();
	    } else if (lin != nunx.l) {
	      inc_left_char_count();
	    }
	  }
	  if (rin != neos) {
	    if (nunx.r == neos) {
	      nunx.r = rin;
	      inc_right_char_count();
	    } else if (rin != nunx.r) {
	      inc_right_char_count();
	    }
	  }
	}
	if (lin != nunx.l || rin != nunx.r || (phase != 0 && indel)) {
	  // transition to xmer
	  lout = nunx.pl;
	  char r0 = nunx.r;
	  char l0 = nunx.l;
	  rout = r0;
 	  if (rin != r0 && phase == 0) {
	    // transition to nux
	    setx();
	    nux.r  = (((unsigned)1)<<r0);
	    nux.r |= (((unsigned)1)<<rin);
	    nux.l = l0;
	    unsetu();
	  } else { // rin == unx.r
	    // transition to ux
	    setx();
	    ux.r = r0;
	    ux.l = l0;
	    setu();
	  }
	  return true;
	} else if (phase == 0) { // we are (still) not unique, but not xmer
	  nunx.pl->push_front(p);
	  unset_one_exact_pos();
	  return false;
	}
      }
    }
    return false;
  }
  long unsigned int ind() {
    if (u()) {
      return ux.xindex;
    } else {
      return nux.xindex;
    }
    return 0;
  }
  
};

long unsigned int merelt::xmerind = 0;
long unsigned int merelt::nxmers = 0;
long unsigned int merelt::totmers = 0;

typedef tinylist<merelt> htelt;

long unsigned int hash(const char *buff, int l, unsigned long int base, long unsigned int p) {
  long unsigned h = 0;
  for (int i=0;i<l;i++) {
    h = (h*base+buff[i])%p;
  }
  return h;
}

/* 
#include "primegen.h"
void random_primes_lt(unsigned long int m, std::vector<long unsigned int> &p) {
  primegen pg;
  primegen_init(&pg);
  uint64 np = primegen_count(&pg,(uint64)m);
  uint64 pr;
  
  // checkpoint;
  // cerr << m << " " << np << endl;
  
  srand48(time(NULL));
  DRAND48;
  std::vector<long unsigned int> pi(p.size());
  for (int i=0;i<p.size();i++) {
    pi[i] = (long unsigned int) floor(np*DRAND48);
    // checkpoint;
    // cerr << pi[i] << endl;
  }
  sort(pi.begin(),pi.end());
  primegen_init(&pg);
  int j = 0,k = 0;
  while (j<=pi[p.size()-1]) {
    // checkpoint;
    pr = primegen_next(&pg);
    // cerr << j << " " << pr << endl;
    while (j == pi[k] && k < p.size()) {
      p[k] = pr;
      k++;
    }
    j++;
  }

}
*/

bool htlookup(std::vector<htelt> & htable, 
	      long unsigned int h, const char *s, int k,
	      tinylist<merelt>::iterator & mitret) {

  bool found = false;
  tinylist<merelt>::iterator mit=htable[h].begin();
  tinylist<merelt>::iterator mit0=htable[h].end();
  while (mit != htable[h].end()) {
    int cmp = mit->strcmp(s,k);
    if (cmp == 0) {
      found = true;
      break;
    } else if (cmp < 0) {
      break;
    }
    mit0 = mit;
    ++mit;
  }
  if (found) {
    mitret = mit;
  } else {
    mitret = mit0;
  }
  return found;
}

void dump_xmer_cnt(ostream *o,
		   long unsigned int u,
		   long unsigned int xind,
		   char *s,
		   int mersize,
		   char r,
		   long signed int cnt,
		   FastaFile<Lazy_Header_SI>* ff,
		   bool noindex) {
  if (noindex) {
    (*o) << u << ":"
	 << xind << " "
	 << ff->ch(r) << " - "
	 << cnt << endl;
  } else {
    (*o) << u << ":"
	 << xind << " "
	 << (int)r << " ";
    for (int i=0;i<mersize;i++) {
      (*o) << ff->ch(s[i]);
    }
    (*o) << " " << ff->ch(r); 
    (*o) << " " << cnt << endl;
  }
}

void dump_xmer_loc(ostream *o,
		   long unsigned int u,
		   long unsigned int xind,
		   char *s,
		   FILE_POSITION_TYPE e,
		   int mersize,
		   char r,
		   FastaFile<Lazy_Header_SI>* ff,
		   bool noindex,
		   char *actual_seq,
		   int phase=0) {
  if (noindex) {
      (*o) << u << ":"
	   << xind << " "
	   << ff->ch(r) << " ";
      if (e >= 0) {
	(*o) << e+1 << endl;
      } else {
	(*o) << e << endl;
    }
  } else {
    if (e >= 0) {
      const Lazy_Header_SI & h = ff->get_header_data(e);
      (*o) << u << ":"
	   << xind << " "
	   << (int)r << " "
	   << e-mersize << " "
	   << e+1 << " ";
      for (int i=0;i<mersize;i++) {
	(*o) << ff->ch(s[i]);
      }
      (*o) << " " << ff->ch(r); 
      (*o) << " " << h.index() << " "
	   << h.short_header() << endl;
    } else {
      const Lazy_Header_SI & h = ff->get_header_data(-e);
      (*o) << u << ":"
	   << xind << " "
	   << (int)r << " "
	   << -e << " "
	   << -e-(mersize+1) << " ";
      for (int i=0;i<mersize;i++) {
	(*o) << ff->ch(s[i]);
      }
      (*o) << " " << ff->ch(r); 
      (*o) << " " << h.index() << " "
	   << h.short_header() << endl;      
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

  std::list<std::string> patterns;
  std::list<std::string>::iterator patit;

  istrstream sis(opt.patterns.c_str());
  std::string pattern;
  while (sis >> pattern) {
    if (pattern.length() != opt.mersize) {
      cerr << "Bad pattern: " << pattern << " has length "
	   << pattern.length() << " != " << opt.mersize << endl;
      exit(1);
    }
    patterns.push_back(pattern);
  }
  
  if (patterns.empty()) {
    exit(0);
  }

  for (patit=patterns.begin();patit!=patterns.end();++patit) {
    uppercase(*patit);
  }
  
  FastaFile<Lazy_Header_SI>* ff;
  fasta_file_seq_params ffp;
  ffp.upper_case = true;
  ffp.eos_char = opt.eos_char;
  ffp.check_params = true;
  ff = pick_fasta_file<Lazy_Header_SI>(opt.database,opt.dbind,opt.memmap,
				       !opt.noindex,ffp,opt.verbose);

  
  char neos=ff->nch(opt.eos_char);

  // Set up hash table...
  long unsigned int maxp = (((unsigned long int)1)<<27);
  std::vector<long unsigned int> primes(1);
  rand_hash_table::random_primes_lt(maxp,primes);
  long unsigned int prime = primes[0];
  if (opt.hashtablesize != -1) {
    prime = opt.hashtablesize;
  }
  // long unsigned int prime = 37;
  // if (opt.verbose) {
  //  checkpoint;
  //  timestampli("Hash table size: ",prime);
  // }
  std::vector<htelt> htable(prime);
  std::vector<char*> htablekeys;
  std::vector<char*>::iterator htkit;
  char *buffer = new char[opt.mersize];
  char *buffer1 = new char[opt.mersize+2];
  char left,right;

  PatternMatch* pm;

  int rcstart = 0;

  int phase=0;
  int block=0;
  int nblocks=0;
  int blockstart=0;
  int blockstart0=0;
  while (phase < ((opt.nmismatch==0)?1:3)) {
    ff->reset();
    
    if (opt.verbose) {
      ostrstream ss;
      if (nblocks >= 0) {
	ss << "Starting batch " << block+1 << "/" << ((nblocks==0)?1:nblocks) << " of phase " << phase << ends;
      } else {
	ss << "Starting batch " << block+1 << "/?? of phase " << phase << ends;      
      }
      std::string v(ss.str());
      timestamp(v.c_str());
    }

    std::vector<long> counts(0);
    if (phase == 2) {
      counts.resize(opt.blocksize,-1);
      // fill(counts.begin(),counts.end(),-1);
    }
    int patcnt=0,patcntbl=0,skipped=0;
    switch (phase) {
    case 0:
      pm = new shift_and(false,false,true,opt.eos_char);
      for (patit=patterns.begin();patit!=patterns.end();++patit) {
	// checkpoint;
	// cerr << "Add pattern: " << *patit << endl;
	pm->add_pattern(*patit);
	if (opt.rc) {
	  // cerr << "Add pattern: " << reverse_comp(*patit) << endl;
	  pm->add_pattern(reverse_comp(*patit));
	}
	patcnt++;
      }
      patcntbl = patcnt;
      break;
    case 1: 
    case 2: 
      {
	std::vector<int> patlen(2,opt.mersize+((phase==2)?1:0));
	std::vector<std::pair<int,int> > patconst(2,make_pair(0,0));
	pm = pick_pattern_index(ff,(phase==1?opt.node1:opt.node2),
				opt.nmismatch,&patconst,&patlen,
				0,false,false,opt.indels,false,
				opt.eos_char,opt.verbose);
	tinylist<merelt>::iterator mit;
	long unsigned int h;
	for (htkit=htablekeys.begin();htkit!=htablekeys.end();++htkit) {
	  long unsigned int h = hash(*htkit,opt.mersize,ff->size(),prime);
	  if (htlookup(htable,h,*htkit,opt.mersize,mit)) {
	    if (phase == 1) {
	      if (!mit->x()) {
		if (patcnt < blockstart) {
		  patcnt++;
		  continue;
		}
		if (patcntbl >= opt.blocksize) {
		  break;
		}
		for (int i=0;i<opt.mersize;i++) {
		  buffer1[i] = ff->ch((*htkit)[i]);
		}
		buffer1[opt.mersize] = 0;
		// checkpoint;
		// cerr << "Add pattern: " << buffer1 << endl;
		pm->add_pattern(buffer1);
		if (opt.rc) {
		  // cerr << "Add pattern: " << reverse_comp(buffer1) << endl;
		  pm->add_pattern(reverse_comp(buffer1));
		}
		patcntbl++;
	      }
	      patcnt++;
	    } else if (phase == 2 && mit->x()) {
	      for (int i=0;i<opt.mersize;i++) {
		buffer1[i] = ff->ch((*htkit)[i]);
	      }
	      buffer1[opt.mersize+1] = 0;
	      for (int j=0;j<ff->size();j++) {
		if (mit->rchar(j) && j != neos) {
		  if (patcnt < blockstart) {
		    patcnt++;
		    continue;
		  }
		  if (patcntbl >= opt.blocksize) {
		    break;
		  }
		  buffer1[opt.mersize] = ff->ch(j);
		  if (!(opt.uniq || opt.exonly) ||
		      (mit->get_left_char_count() <= 1 &&
		       mit->get_right_char_count() <= 1 &&
		       (opt.exonly || (opt.uniq && mit->one_exact_pos())))) {
		    pm->add_pattern(buffer1,patcntbl*((opt.rc)?2:1)+1);
		    if (opt.rc) {
		      pm->add_pattern(reverse_comp(buffer1),patcntbl*((opt.rc)?2:1)+2);
		    }
		    counts[patcntbl]=0;
		  } else {
		    skipped++;
		  }
		  patcntbl++;
		  patcnt++;
		}
	      }
	    }
	  }
	}
	blockstart0 = blockstart;
	blockstart = patcnt;
      }
      break;
    }

    if (opt.verbose) {
      ostrstream ss;
      ss << (patcntbl-skipped) << " patterns" << ends;
      std::string v(ss.str());
      timestamp(v.c_str());
    }

    // checkpoint;

    std::vector<int> ignore(0);
    if (phase == 1) {
      ignore.resize(patcntbl,0);
      // fill(ignore.begin(),ignore.end(),0);
    }
    if (phase == 2) {
      // Must be smaller...
      counts.resize(patcntbl);
      // fill(counts.begin(),counts.end(),-1);
    }
    // checkpoint;

    if (opt.verbose) {
      pm->progress_interval(*ff);
    }

    // checkpoint;

    pm->init(*ff);
  
    // checkpoint;

    pattern_hit_vector l(1000);
    pattern_hit_vector::iterator it; 
    while ((patcntbl>0)&&(pm->find_patterns(*ff,l,1000)||!l.empty())) {
      // checkpoint;
      FILE_POSITION_TYPE oldcharspos;
      oldcharspos = ff->pos();
      it = l.begin();
      while (it != l.end()) {
	int hitid = it->value().first->id()-1;
	int hitrc = false;
	FILE_POSITION_TYPE hitpos = it->key();
	// checkpoint;
	// cerr << hitid << " "
	// << hitrc << " " 
	// << hitpos << endl;
	if (opt.rc) {
	  if (hitid%2==1) {
	    hitrc = true;
	    hitpos = -hitpos;
	  }
	  hitid = hitid/2;
	}
	// checkpoint;
	// cerr << hitid << " "
	//      << hitrc << " " 
	//      << hitpos << endl;
	if (phase == 2) {
	  if (!opt.exonly || it->value().second > 0) {
	    // checkpoint;
	    // cerr << (int)it->value().second << endl;
	    // cerr << hitid << endl;
	    counts[hitid]++;
	  }
	  ++it;
	  continue;
	} else if (phase == 1 && ignore[hitid]>0) {
	  if (ignore[hitid] == 2) {
	    checkpoint;
	    cerr << it->value().first->pattern() << endl;
	    exit(1);
	  }
	  ++it;
	  continue;
	}
	long unsigned int h;
	bool found;
	tinylist<merelt>::iterator mit;
	std::string pat = it->value().first->pattern();
	// checkpoint;
	// cerr << hitid << " "
	// << hitrc << " " 
	// << hitpos << " "
	// << pat << endl;
	if (hitrc) {
	  pat = reverse_comp(pat);
	}
	// checkpoint;
	// cerr << hitid << " "
	// << hitrc << " " 
	// << hitpos << " "
	// << pat << endl;
	if (phase == 1) {
	  for (int i=0;i<opt.mersize;i++) {
	    buffer[i] = ff->nch(pat[i]);
	  }
	  h = hash(buffer,opt.mersize,ff->size(),prime);
	  found = htlookup(htable,h,buffer,opt.mersize,mit);
	}
	ff->pos(it->key()-opt.mersize-1);
	if (!hitrc) {
	  left = ff->getnch();
	  for (int i=0;i<opt.mersize;i++) {
	    buffer[i] = ff->getnch();
	  }
	  right = ff->getnch();
	} else {
	  right = ff->nch(iupac_revcomp(ff->getch()));
	  for (int i=opt.mersize-1;i>=0;i--) {
	    buffer[i] = ff->nch(iupac_revcomp(ff->getch()));
	  }
	  left = ff->nch(iupac_revcomp(ff->getch()));
	}
	if (phase == 0) {
	  h = hash(buffer,opt.mersize,ff->size(),prime);
	  found = htlookup(htable,h,buffer,opt.mersize,mit);
	} 
	
	//checkpoint;
	//for (int i=0;i<opt.mersize;i++) {
	//buffer1[i] = ff->ch(buffer[i]);
	// }
	// buffer1[opt.mersize] = 0;
	// cerr << it->key() << " " << pat << " " << buffer1 << " " << ff->ch(left) << " " << ff->ch(right) << endl;
	
	tinylist<merelt>::iterator mit1;
	if (!found) {
	  assert(phase==0);
	  char *str = new char[opt.mersize];
	  memcpy(str,buffer,opt.mersize);
	  htablekeys.push_back(str);
	  // checkpoint;
	  // cerr << *(htablekeys.rbegin()) << endl;
	  if (mit == htable[h].end()) {
	    mit1 = htable[h].push_front(merelt(str,left,right,neos,hitpos));
	  } else {
	    mit1 = htable[h].insert_after(mit,merelt(str,left,right,neos,hitpos));
	  }
	  mit = mit1;
	  // checkpoint;
	  // mit->dump(opt.out,(*ff),opt.mersize);
	  if (left == neos || right == neos) { 
	    mit->setx();
	    dump_xmer_loc(opt.out,opt.u,mit->ind(),mit->str,hitpos,opt.mersize,right,ff,opt.noindex,buffer,phase);
	    // checkpoint;
	    // mit->dump(opt.out,(*ff),opt.mersize);
	    // (*opt.out) << endl;
	  }
	} else if (phase == 0 || mit->strcmp(buffer,opt.mersize)!=0){
	  tinylist<FILE_POSITION_TYPE> *pos=0;
	  char rout;
	  // checkpoint;
	  // mit->dump(opt.out,(*ff),opt.mersize);
	  if (mit->update(left,right,neos,phase,opt.indels,hitpos,pos,rout)) {
	    if (pos) {
	      tinylist<FILE_POSITION_TYPE>::iterator lit=pos->begin();
	      while (lit != pos->end()) {
		dump_xmer_loc(opt.out,opt.u,mit->ind(),mit->str,(*lit),opt.mersize,rout,ff,opt.noindex,buffer,phase);
		++lit;
	      }
	      delete pos;
	      if (phase == 1) {
		ignore[hitid] = 1;
		// checkpoint;
		// cerr << "setting ignore["<<it->value().first->pattern()<<"]"<<endl;
		// pm->del_pattern(*ff,it->value().first->pattern());
	      }
	    }
	    if (phase == 0) {
	      dump_xmer_loc(opt.out,opt.u,mit->ind(),mit->str,hitpos,opt.mersize,right,ff,opt.noindex,buffer,phase);
	    }
	  } 
	  // checkpoint;
	  // mit->dump(opt.out,(*ff),opt.mersize);
	  // (*opt.out) << endl;
	}
	++it;
      }
      if (opt.delpat && (phase == 1 || (phase == 2 && (opt.exonly || opt.uniq)))) {
	it = l.begin();
	while (it != l.end()) {
	  int hitid = it->value().first->id()-1;
	  if (opt.rc) {
	    hitid = hitid/2;
	  }
 	  if ((phase == 1 && ignore[hitid]==1)||
	      (phase == 2 && ((opt.uniq && counts[hitid]>1)||
			      (opt.exonly && counts[hitid]>0)))) {
	    if (phase == 1) {
	      ignore[hitid] = 2;
	    } else if (phase == 2) {
	      if (opt.uniq) {
		counts[hitid] = 2;
	      } else if (opt.exonly) {
		counts[hitid] = 1;
	      }
	    }
	    // checkpoint;
	    // cerr << "Calling del_pattern for " << it->value().first->pattern() << endl;
	    pm->del_pattern(*ff,it->value().first->pattern());
	    if (opt.rc) {
	      // cerr << "Calling del_pattern for " << reverse_comp(it->value().first->pattern()) << endl;
	      pm->del_pattern(*ff,reverse_comp(it->value().first->pattern()));
	    }
	  }
	  ++it;
	}
      }
      l.clear();
      ff->pos(oldcharspos);
    }
    if (phase == 0) {
      block=0;
      blockstart=0;
      nblocks = (int)ceil((((double)merelt::totmers)-merelt::nxmers)/opt.blocksize);
      phase++;
      if (opt.verbose) {
	ostrstream ss;
	ss << "At end of phase 0, total mers: "
	   << merelt::totmers 
	   << " xmers: "
	   << merelt::nxmers 
	   << " non-xmers: " 
	   << merelt::totmers-merelt::nxmers << ends;
	std::string v(ss.str());
	timestamp(v.c_str());
      }
    } else if (phase == 1) {
      block++;
      if (block>=nblocks) {
	block=0;
	blockstart=0;
	nblocks = -1;
	phase++;
	if (opt.verbose) {
	  ostrstream ss;
	  ss << "At end of phase 1, total mers: "
	     << merelt::totmers 
	     << " xmers: "
	     << merelt::nxmers 
	     << " non-xmers: " 
	     << merelt::totmers-merelt::nxmers << ends;
	  std::string v(ss.str());
	  timestamp(v.c_str());
	}
      }
    } else if (phase == 2) {
      block++;
      if (patcntbl < opt.blocksize) {
	phase++;
      }
      // We should dump out the counts...
      int patcnt=0,patcntbl=0;
      tinylist<merelt>::iterator mit;
      long unsigned int h;
      for (htkit=htablekeys.begin();htkit!=htablekeys.end();++htkit) {
	long unsigned int h = hash(*htkit,opt.mersize,ff->size(),prime);
	if (htlookup(htable,h,*htkit,opt.mersize,mit)) {
	  if (mit->x()) {
	    for (int j=0;j<ff->size();j++) {
	      if (mit->rchar(j) && j != neos) {
		if (patcnt < blockstart0) {
		  patcnt++;
		  continue;
		}
		if (patcntbl >= opt.blocksize) {
		  break;
		}
		long cnt = counts[patcntbl];
		if (cnt < 0) {
		  if (opt.uniq) {
		    cnt = 2;
		  } else if (opt.exonly) {
		    cnt = 1;
		  }
		}
		dump_xmer_cnt(opt.out,opt.u,mit->ind(),mit->str,opt.mersize,j,cnt,ff,opt.noindex);
		patcnt++;
		patcntbl++;
	      }
	      if (patcntbl >= opt.blocksize) {
		break;
	      }
	    }
	  }
	}
      }
    }
    delete pm;
  }
  return 0;
}

