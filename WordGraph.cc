/* posix compliant mmap */
#include "util.h"
# include <fcntl.h>
# include <errno.h>
# include <sys/types.h>
# include <unistd.h>
# include <stdlib.h>
# include <sys/stat.h>
# include <sys/shm.h>
# ifndef __USE_MISC
#   define __USE_MISC
# endif
# ifndef __MINGW32__
#   include <sys/mman.h>
#   if !defined(MAP_SHARED)
#     error MAP_SHARED
#   endif
#endif

#ifdef __alpha
#define MMAPFLAGS    (MAP_FILE | MAP_VARIABLE | MAP_SHARED)
#endif

#ifdef __linux
#define MMAPFLAGS    (MAP_FILE | MAP_SHARED)
#endif

#ifdef __FreeBSD__
#define MMAPFLAGS    (MAP_FILE | MAP_SHARED)
#endif

#ifdef _AIX
#define MMAPFLAGS    (MAP_FILE | MAP_SHARED)
#endif

#ifdef __CYGWIN__
#define MMAPFLAGS    (MAP_FILE | MAP_SHARED)
#endif

#ifdef __sun
#define MMAPFLAGS    (MAP_SHARED)
#endif

/* posix mmap */

#include <assert.h>
#include <sys/stat.h>

#include <cstdio>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
// #include <utils>
using namespace std;

static unsigned VERBOSE=0;
static char *oname = NULL;
static char *fname = NULL;
static char *xname = NULL;
static char *whoami;
static bool Amino=false;
static bool UC=false;
static unsigned MerSize=0;
static char *smap =0;
static long unsigned slen=0;
// static char MAP[256];
static bool collapse_nodes=true;
static bool distinguish_eos=false;
static bool keep_eos=false;
static bool aa_map=true;

// typedef map<pair<long unsigned int,long unsigned int>,int> ndcmpcache;
// static ndcmpcache merseqcmpcache;

void usage(FILE *f) {
  fprintf(f,
	  "%s: [-v] [-h] [-A] [-e] [-E] -x xspace_file -m mer_size -f fwd_file [ -o output ]\n",
	  whoami);
}

void unmap() {
  if (smap) {
#ifdef __CYGWIN__
    (void)munmap((char*)smap, slen);
#elif defined(__MINGW32__) 
    delete [] smap;
#else
    (void)munmap(smap, slen);
#endif
    smap = 0;
  }
}

static int get_args(int argc, const char **argv) {
  int errflg = 0;
  whoami = strdup(*argv);
  argv++;
#define OPT_ARG ( (argv[0][1])?(++(*argv)):(*(++argv)) )
  while(!errflg && *argv) {
    if (**argv == '-') {
      (*argv)++;
      while (!errflg && **argv) {
        switch (**argv) {
        case 'v':
          VERBOSE++;
          goto loopin;
        case 'E':
	  distinguish_eos = true;
          goto loopin;
        case 'e':
	  keep_eos = true;
          goto loopin;
        case 'A':
	  aa_map = true;
          goto loopin;
	  //         case 'A':
	  // 	  Amino = true;
	  //           goto loopin;
	  //         case 'U':
	  // 	  UC = true;
	  //           goto loopin;
        case 'm':
          MerSize = atoi(OPT_ARG);
          goto loopout;
        case 'o':
          oname = strdup(OPT_ARG);
          goto loopout;
        case 'f':
          fname = strdup(OPT_ARG);
          goto loopout;
        case 'x':
          xname = strdup(OPT_ARG);
          goto loopout;
        case 'h':
          usage(stdout);
          exit(0);
          goto loopin;
        default:
          cerr << whoami << ": unknown flag '-" << *argv <<"'"<< endl;
          errflg++;
        }
      loopin:
        (*argv)++;
      }
    loopout:
      argv++;
    }
    else {
      errflg++;
    }
  }
  if (!xname) ++errflg;
  if (!fname) ++errflg;
  if (MerSize==0) ++errflg;
  
  return errflg;
}

template < class basic_t, unsigned int n, unsigned int val_bits >
class largeword {
private:
# define basic_n (((n*val_bits)-1)/(sizeof(basic_t)*8)+1)
  basic_t data_[basic_n];
public:
  largeword() {
    assert(val_bits <= sizeof(basic_t)*8);
    for (unsigned int i=0;i<basic_n;i++) data_[i]=0;
  }
  largeword(const largeword & w) {
    for (unsigned int i=0;i<basic_n;i++) data_[i] = w.data_[i];
  }
  static unsigned int size() {
    return sizeof(basic_t)*basic_n*8/val_bits;
  }
  inline long unsigned int operator[](unsigned int i) const {
    return value(i);
  }
  long unsigned int value(unsigned int i) const {
    assert( i < size());
    // cerr << "I: " << i << endl;
    unsigned int start_bit = val_bits*i;
    unsigned int end_bit = start_bit+val_bits-1;
    unsigned int shift;
    unsigned int chunk;
    long unsigned int retval=0;
    shift = start_bit%(sizeof(basic_t)*8);
    chunk = start_bit/(sizeof(basic_t)*8);
    // cerr << "2: " << start_bit << " " << end_bit << endl;
    // cerr << "4: " << shift << " " << chunk << endl;
    if (shift!=0) {
      retval = data_[chunk]>>shift;
    } else {
      retval = data_[chunk];
    } 
    if (end_bit/(sizeof(basic_t)*8) > start_bit/(sizeof(basic_t)*8)) {
      chunk++;
      shift = (sizeof(basic_t)*8)-shift;
      retval |= data_[chunk]<<shift;
    }
    long unsigned int mask=(((long unsigned int)1)<<val_bits)-1;
    return (retval & mask);
  }
  void set(unsigned int i, long unsigned int v) {
    assert( i < size());
    // cerr << "I: " << i << " " << v << endl;
    unsigned int start_bit = val_bits*i;
    unsigned int end_bit = start_bit+val_bits-1;
    // cerr << "2: " << start_bit << " " << end_bit << endl;
    unsigned int shift;
    unsigned int chunk;
    long unsigned int mask=(((long unsigned int)1)<<val_bits)-1;
    long unsigned int val=v&mask;
    shift = start_bit%(sizeof(basic_t)*8);
    chunk = start_bit/(sizeof(basic_t)*8);
    // cerr << "3: " << mask << " " << val << endl;
    // cerr << "4: " << shift << " " << chunk << endl;
    // cerr << chunk << ": ";
    // cerr.width(3);
    // cerr << oct << (unsigned long int)data_[chunk] << endl;
    if (shift!=0) {
      data_[chunk] |= val<<shift;
    } else {
      data_[chunk] = val;
    } 
    // cerr << chunk << ": ";
    // cerr.width(3);
    // cerr << oct << (unsigned long int)data_[chunk] << endl;
    if (end_bit/(sizeof(basic_t)*8) > start_bit/(sizeof(basic_t)*8)) {
      chunk++;
      shift = (sizeof(basic_t)*8)-shift;
      // cerr << "5: " << shift << " " << chunk << endl;
      // cerr << chunk << ": ";
      // cerr.width(3);
      // cerr << oct << (unsigned long int)data_[chunk] << endl;
      data_[chunk] |= val>>shift;
      // cerr << chunk << ": ";
      // cerr.width(3);
      // cerr << oct << (unsigned long int)data_[chunk] << endl;
    } 
    // assert(value(i)==v);
  }
};

//typedef largeword<unsigned char, 33, 5> nodechca_t;
//typedef largeword<unsigned char, 20, 3> nodechca_t;
typedef largeword<unsigned char, 100, 3> nodechca_t;

class node {
private:
  static bool aamap_;
public:
  unsigned int label;
  nodechca_t chars;
  node(unsigned int l) : label(l) {
    assert(chars.size() >= MerSize+1);
  }
  void getseq(long unsigned int p) {
    if (chars[0] == 0) {
      if (p == 1) {
	p = slen;
      } 
      char *retval = smap+(p-MerSize);
      for (unsigned int i=0;i<MerSize;i++) {
	chars.set(i,map(retval[i]));
      }
    }
  };
  static char map(long unsigned int v) {
    if (aamap_) {
      assert(v >= 1 && v <= 27);
      if (v <= 26) {
	return ('A'-1+v);
      } else {
	return '$';
      }
    } else {
      switch (v) {
      case 1:
	return 'A';
	break;
      case 2:
	return 'C';
	break;
      case 3:
	return 'G';
	break;
      case 4:
	return 'T';
	break;
      case 5:
	return 'N';
	break;
      case 6:
	return '$';
	break;
      default:
	assert(0);
	break;
      }
      return '\0';
    }
  }
  static unsigned int map(char c) {
    if (aamap_) {
      assert (c >= 'A' && c <= 'Z' || c == '$');
      if (c >= 'A' && c <= 'Z') {
	return (c+1-'A');
      } else {
	return ('Z'+1-'A'+1);
      }
    } else {
      switch (c) {
      case 'A': 
	return ((unsigned int) 1); break;
      case 'C': 
	return ((unsigned int) 2); break;
      case 'G': 
      return ((unsigned int) 3); break;
      case 'T': 
      return ((unsigned int) 4); break;
      case 'N': 
	return ((unsigned int) 5); break;
      case '$': 
	return ((unsigned int) 6); break;
      default: assert(0);
      }
      return 0;
    }
  }
  static void set_aa() {
    aamap_ = true;
  }  
  static void set_nuc() {
    aamap_ = false;
  }
  bool contains_term() {
    if (node::map(chars[0]) == '$') return true;
    for (int i=MerSize-1; i>=1; i--) {
      if (node::map(chars[i]) == '$') return true;
    }
    return false;
  }
  bool operator<(const node & a) {
//     static unsigned long int cmpcount=0;
//     cmpcount++;
//     if (cmpcount % 100000 == 0) {
//       fprintf(stderr,".");
//       if (cmpcount % 1000000 == 0) fprintf(stderr," ");
//       if (cmpcount % 10000000 == 0) fprintf(stderr," %10.3e\n",(double)cmpcount);
//       fflush(stderr);
//     }
    
    for (unsigned int i=0; i<MerSize; i++) {
      if (chars[i] < a.chars[i]) return true;
      if (chars[i] > a.chars[i]) return false;
    }
    return false;
  }
  bool operator==(const node & a) {
    for (int i=MerSize-1; i>=0; i--) {
      if (chars[i] != a.chars[i]) return false;
    }
    return true;
  }
};

bool node::aamap_ = false;

vector<node> nodes;

class space {
public:
  unsigned long pos_;
  unsigned char ch_;
  unsigned int  node_;
  
  space(long unsigned s, unsigned int x, char c) {
    pos(s);
    node(x);
    ch(c);
  };
  inline void pos(long unsigned p) {
    pos_ = p;
  }
  inline long unsigned pos() const {
    return pos_;
  }
  inline void node(unsigned int x) {
    node_ = ((unsigned int)x);
  }
  inline unsigned int node() const {
    return node_;
  }
  inline char ch() const {
    return ch_;
  }
  inline void ch(char c) {
    ch_ = c;
  }
  inline bool operator<(const space &s) const {
//     static unsigned long int cmpcount=0;
//     cmpcount++;
//     if (cmpcount % 100000 == 0) {
//       fprintf(stderr,"o");
//       if (cmpcount % 1000000 == 0) fprintf(stderr," ");
//       if (cmpcount % 10000000 == 0) fprintf(stderr," %10.3e\n",(double)cmpcount);
//       fflush(stderr);
//     }
    if (pos() != s.pos()) return pos() < s.pos();
    if (nodes[node()].label != nodes[s.node()].label) 
      return nodes[node()].label < nodes[s.node()].label;
    return ch() < s.ch();
  }
  inline bool operator==(const space &s) const {
    if (pos() != s.pos()) return false;
    if (nodes[node()].label != nodes[s.node()].label) return false;
    return (ch() == s.ch());
  }
  inline char *getseq() const {
    long unsigned int p=pos();
    if (p == 1) {
      p = slen;
    } 
    return smap+(p-MerSize);
  }
};

vector<space> spaces;

class edge {
public:
  unsigned int i;
  unsigned int j;

  edge(unsigned int const & a, 
       unsigned int const & b): i(a),j(b) {};
  inline bool operator<(const edge &e) const {
//     static unsigned long int cmpcount=0;
//     cmpcount++;
//     if (cmpcount % 100000 == 0) {
//       fprintf(stderr,"x");
//       if (cmpcount % 1000000 == 0) fprintf(stderr," ");
//       if (cmpcount % 10000000 == 0) fprintf(stderr," %10.3e\n",(double)cmpcount);
//       fflush(stderr);
//     }
    if (nodes[spaces[i].node()].label != nodes[spaces[e.i].node()].label) return nodes[spaces[i].node()].label < nodes[spaces[e.i].node()].label;
    if (nodes[spaces[j].node()].label != nodes[spaces[e.j].node()].label) return nodes[spaces[j].node()].label < nodes[spaces[e.j].node()].label;
    return spaces[i].ch() < spaces[e.i].ch();
  }
  inline bool operator==(const edge &e) const {
    return 
      nodes[spaces[i].node()].label == nodes[spaces[e.i].node()].label && 
      nodes[spaces[j].node()].label == nodes[spaces[e.j].node()].label && 
      spaces[i].ch() == spaces[e.i].ch();
  }
};

vector<edge>  edges;

bool ndindlt(unsigned int const & a,
	     unsigned int const & b) {
  return (nodes[a] < nodes[b]);
}

bool contains_term(node & a) {
  return a.contains_term ();
}

void process_graph(FILE *f,FILE *o) {
  
checkpoint;

  static unsigned graphno=0;
  graphno++;
  long unsigned spc;
  char ch,term;
  int c=0,id=0,cnt=0;
  bool push_node;
  
  static unsigned long int cmpcount=0;

  for(int c=0; (c=fgetc(f)) != EOF; ++id) {
    //printf("read char %c(%d)\n",char(c),c);
    cnt = 0;
    while(c==' ' && c!=EOF) {
      fscanf(f,"%lu.%c",&spc,&ch);
      if (id==0)
	term = ch;
      else {
	spaces.push_back(space(spc,id-1,ch));
	if (cnt==0) {
	  nodes.push_back(node(id));
          cmpcount++;
          if (cmpcount % 100000 == 0) {
            fprintf(stderr,".");
            if (cmpcount % 1000000 == 0) fprintf(stderr," ");
            if (cmpcount % 10000000 == 0) fprintf(stderr," %10.3e\n",(double)cmpcount);
            fflush(stderr);
          }
        }
      }
      cnt++;
      c = fgetc(f);
    }
    //printf("end of line and c=%c\n",c);
    if (c != '\n') {
      fprintf(stderr,"format error at graph %u after node %d\n",graphno,id);
      fprintf(stderr,"c = %c(%d)\n",char(c),c);
      exit(33);
    }
    if (!cnt) break;
  }
  if (c == EOF || !id) {
    // nothing to process, else c would be '\n'
    return;
  }
  fprintf(stderr,"\n");
  checkpoint;

  sort(spaces.begin(),spaces.end());
  checkpoint;

  for (unsigned int i=0;i<spaces.size();++i) {
    nodes[spaces[i].node()].getseq(spaces[i].pos());
  }
  unmap();

  if (collapse_nodes) {
    
checkpoint;
  
    vector<unsigned int> nodeinds;
    nodeinds.resize(nodes.size());
    for (unsigned int i=0;i<nodeinds.size();i++) {
      nodeinds[i] = i;
    }
  
    stable_sort(nodeinds.begin(),nodeinds.end(),ndindlt);
    fprintf(stderr,"\n");
checkpoint;

    unsigned int next_node=1;
    nodes[nodeinds[0]].label = next_node;
    
    next_node++;
    for(unsigned i=1; i < nodeinds.size(); ++i) {
      if (nodes[nodeinds[i-1]] == nodes[nodeinds[i]]) {
	nodes[nodeinds[i]].label = nodes[nodeinds[i-1]].label;
      } else {
	nodes[nodeinds[i]].label = next_node;
	next_node++;
      }
    }
checkpoint;
    
//     for(unsigned i=0; i < nodeinds.size(); ++i) {
//       fprintf(stderr,"%10u: ",nodes[nodeinds[i]].label);
//       for (unsigned int j=0;j<MerSize;++j) {
// 	fprintf(stderr,"%c",node::map(nodes[nodeinds[i]].chars[j]));
//       }
//       fprintf(stderr," %s\n",contains_term(nodes[nodeinds[i]])?"true":"false");
//     }

  }

  unsigned j = 0;
  for(unsigned i=1; i < spaces.size(); ++i) {
    if (spaces[j].pos() == spaces[i].pos()) {
      /* if (!(spaces[j] == spaces[i])) {
	space & s0 = spaces[j];
	node & n0 = nodes[s0.node];
	fprintf(stderr,"%7d %7d %-7s\n%.*s\n",n0.label,s0.node,
		(contains_term(n0)?"true":"false"),
		MerSize,n0.getseq());
	fprintf(stderr,"space: %d %u %c\n",n0.label,s0.pos,s0.ch);
	space & s1 = spaces[i];
	node & n1 = nodes[s1.node];
	fprintf(stderr,"%7d %7d %-7s\n%.*s\n",n1.label,s1.node,
		(contains_term(n1)?"true":"false"),
		MerSize,n1.getseq());
	fprintf(stderr,"space: %d %u %c\n",n1.label,s1.pos,s1.ch);
	} */
    } else {
      node & ni = nodes[spaces[i].node()];
      node & nj = nodes[spaces[j].node()];
      if (!keep_eos) {
	if (!contains_term(nj) && 
	    !contains_term(ni)) {
	  edges.push_back(edge(j,i));
	}
      } else {
	if (!contains_term(nj) || 
	    !contains_term(ni)) {
	  edges.push_back(edge(j,i));
	}
      }
      j = i;
    }
  }

checkpoint;

  sort(edges.begin(),edges.end());
  fprintf(stderr,"\n");

checkpoint;

  unsigned int same=0;
  
  for(unsigned i=1; i < edges.size(); ++i) {
    if ( i!=0 && edges[i-1] == edges[i] && (!distinguish_eos || spaces[edges[i-1].i].ch() != term)) {
      // same;
      fprintf(o,"%c\t%u\t%u\t%lu\t%lu\n",
	      'e',
	      nodes[spaces[edges[i-1].i].node()].label,
	      nodes[spaces[edges[i-1].j].node()].label,
	      spaces[edges[i-1].i].pos(),
	      spaces[edges[i-1].j].pos());
      same++;
    } else {
      // write out edges[i-1];
      fprintf(o,"%c\t%u\t%u\t%lu\t%lu\t%u\n",
	      'E',
	      nodes[spaces[edges[i-1].i].node()].label,
	      nodes[spaces[edges[i-1].j].node()].label,
	      spaces[edges[i-1].i].pos(),
	      spaces[edges[i-1].j].pos(),
	      same+1);
      same = 0;
    }
  }
  fprintf(o,"%c\t%u\t%u\t%lu\t%lu\t%u\n",
	  'E',
	  nodes[spaces[edges[edges.size()-1].i].node()].label,
	  nodes[spaces[edges[edges.size()-1].j].node()].label,
	  spaces[edges[edges.size()-1].i].pos(),
	  spaces[edges[edges.size()-1].j].pos(),
	  same+1);
  fprintf(o,".\n");
checkpoint;

}

int main(int argc,const char **argv) {

  if (get_args(argc,argv)) {
    usage(stderr);
    exit(1);
  }
  FILE *stdverb = stdout;

checkpoint;

 if (aa_map) {
   node::set_aa();
 } 

 fprintf(stderr,"Size of space: %d\n",sizeof(space));
 fprintf(stderr,"Size of node:  %d\n",sizeof(node));
 fprintf(stderr,"Size of edge:  %d\n",sizeof(edge));

 checkpoint;

  FILE *o=stdout;
  if (oname) {
    if (!(o = fopen(oname,"w"))) {
      fprintf(stderr,"%s:main: unable to open %s\n",whoami,oname);
      exit(101);
    }
    stdverb = stderr;
  }
  FILE *xil;
  if (!(xil = fopen(xname,"r"))) {
    fprintf(stderr,"%s:main: unable to open %s\n",whoami,xname);
    exit(103);
  }
  errno = 0;
  
  int f = open(fname, O_RDONLY);
  if ((f < 0) || (errno)) {
    fprintf(stderr, "Couldn't open '%s'\n", fname);
    perror("open");
    exit(1);
  }
  
  struct stat  sb;
  fstat(f, &sb);
  if (errno) {
    fprintf(stderr, "Couldn't stat '%s'\n", fname);
    perror("fstat\n");
    exit(1);
  }
  
  slen = sb.st_size;
  
# ifndef __MINGW32__
  smap = (char *)mmap(0L, slen, PROT_READ, MMAPFLAGS, f, 0);
  // smap = (char *)shmat(f,0,SHM_MAP);
# else
  smap = new char[slen];
  read(f,smap,slen);
# endif
  
  if (errno) {
    fprintf(stderr, "Couldn't map '%s'\n", fname);
    // fprintf(stderr, "Couldn't shmat '%s'\n", fname);
    perror("mmap");
    // perror("shmat");
    exit(1);
  }
  
  close(f);
    
  // On Brian's advice, touch every page to ensure it is all
  // streamed into memory...
  
  // long unsigned int j=0;
  // for (long unsigned int i=0;
  // i<slen;i+=1024) {
  // if (smap[i] == '\0') j++;
  // }
checkpoint;
  // Hopefully, this won't get optimized out!
  
  //   if (Amino) {
  //     for(int i=0; i < 256; ++i) {
  //       MAP[i] = CharMap.aminoacid()[i];
  //     }
  //     // turn X's into mismatches
  //     for(int i=0; i < 256; ++i) if (MAP[i] == 'X') MAP[i] = CharMap.term1();
  //   } 
  //   else if (UC) {
  //     for(int i=0; i < 256; ++i) {
  //       MAP[i] = CharMap.uppercase()[i];
  //     }
  //     // turn non-letter's into mismatches
  //     for(int i=0; i < 256; ++i) if (MAP[i] == CharMap.term3()) MAP[i] = CharMap.term1();
  //   }
  //   else {
  //     for(int i=0; i < 256; ++i) {
  //       MAP[i] = CharMap.canonical()[i];
  //     }
  //     // turn n's into mismatches
  //     for(int i=0; i < 256; ++i) if (MAP[i] == 'N') MAP[i] = CharMap.term1();
  //   }

  int c;
  while(!feof(xil)) {
    process_graph(xil,o);
  }

  if (oname) fclose(o);
  if (xname) fclose(xil);
  unmap();
}
