#ifndef SUFFIX_TREE_H
#define SUFFIX_TREE_H

#include <iostream>
using namespace std;
#include <cstdlib>
#include <assert.h>
#include "mapFile.h"

//#define PRINT_CRAP
template <class Node,class Suffix> class suffix_tree;

class st_index {
  unsigned leaf:1;
  unsigned good:1;
  unsigned idx:30;

 public:
  st_index(): leaf(0),good(0),idx(0) {};
  st_index(const st_index & a): leaf(a.leaf), good(a.good), idx(a.idx) {
    // fprintf(stderr,"%s:%d st_index copy constructor\n",__FILE__,__LINE__);
  };

  inline bool is_leaf() const {return leaf;}
  inline bool is_noleaf() const {return !leaf;}
  inline bool is_nogood() const {return !good;}
  inline bool is_good() const {return good;}
  inline unsigned index() const {return idx;}

  inline st_index &set(const st_index &st) {*this = st; return *this;}

  inline st_index &set_good() {good = 1; return *this;}
  inline st_index &set_nogood()  {good = 0; return *this;}
  inline st_index &set_index(unsigned i) {idx = i; return *this;}
  inline st_index &set_leaf() {leaf = 1; return *this;}
  inline st_index &set_noleaf() {leaf = 0; return *this;}

  void print(ostream &o=cout) {
    o<<idx<< (leaf?"L":"N") << (good?"G":"g");
  }
};

class basic_suffix {
  template <class Node,class Suffix> friend class suffix_tree;

 protected:
  st_index sib;

 public:
  basic_suffix():sib() {}

};

class basic_node: public basic_suffix {
  template <class Node,class Suffix> friend class suffix_tree;

 protected:
  unsigned head, len;
  st_index child;

 public:
  basic_node(): basic_suffix(),head(0),len(0),child() {}

  inline unsigned length() const {return len;}
};

template <class Node=basic_node,class Suffix=basic_suffix>
class suffix_tree {
  protected:
  const char *S;
  unsigned Slen;
  Node *N;
  Suffix *L;
  unsigned Nsize,Lsize,nN,nL;
  char TERM;
  MapFile *map;

  suffix_tree<Node,Suffix> &operator=(const suffix_tree &);
  suffix_tree(const suffix_tree &);

  public:
  
  suffix_tree(const char *st,unsigned len,char t='\0'): 
  S(st),Slen(len),N(NULL),L(NULL),Nsize(0),Lsize(0),nN(0),nL(0),TERM(t),map(0) {};
  // note: if Lsize > 0 then it must be that Nsize > 0
  ~suffix_tree() {clear();}

  inline st_index root() const {
    st_index ret;
    return ret.set_index(0).set_good().set_noleaf();
  }
  inline bool is_root(st_index n) const {
    return n.index()==0 && n.is_noleaf();
  }
  inline void build() { build(Slen); }
  void build(unsigned suflen);
  void read(const char *f, bool mm=true) {
    if (!mm) {
      // get nN, read all Nn N's;
      // get nL, read all Ln L's;
      // use brute force, read and write unsigned chars.
      FILE *fp = fopen(f,"r");
      fread(&nN,1,sizeof(unsigned),fp);
      fread(&nL,1,sizeof(unsigned),fp);
      resize(nN,nL);
      fread(N,nN,sizeof(Node),fp);
      fread(L,nL,sizeof(Suffix),fp);
      fclose(fp);
    } else {
      map = new MapFile(f);
      char const * data = map->file_string();
      memcpy(&nN,data,sizeof(unsigned)*1);
      data += sizeof(unsigned)*1;
      memcpy(&nL,data,sizeof(unsigned)*1);
      data += sizeof(unsigned)*1;
      N = (Node*)data;
      data += sizeof(Node)*nN;
      L = (Suffix*)data;
    }
  }
  void write(const char *f) {
    FILE *fp = fopen(f,"w");
    fwrite(&nN,1,sizeof(unsigned),fp);
    fwrite(&nL,1,sizeof(unsigned),fp);
    fwrite(N,nN,sizeof(Node),fp);
    fwrite(L,nL,sizeof(Suffix),fp);
    fclose(fp);
  }

  void add(const char *,unsigned,unsigned);

  inline void clear() {resize(0,0);delete map;}
  inline bool empty() const {return !nN && !nL;}

  inline unsigned num_nodes() const {return nN;}
  inline unsigned num_suffix() const {return nL;}
  inline bool is_suffix(st_index si) const {return si.is_leaf();}
  inline bool is_node(st_index si) const {return si.is_noleaf();}

  inline Node &node(st_index si) {return N[si.index()];}
  inline Node &node(unsigned  i) {return N[i];}
  inline Suffix &suffix(st_index si) {return L[si.index()];}
  inline Suffix &suffix(unsigned  i) {return L[i];}
  inline const Node &node(st_index si) const {return N[si.index()];}
  inline const Suffix &suffix(st_index si) const {return L[si.index()];}
  inline const Node &node(unsigned  i) const {return N[i];}
  inline const Suffix &suffix(unsigned  i) const {return L[i];}

  inline const char *str() const {return S;}
  inline const char *str(const st_index n) const {
    if (n.is_leaf())
      return S+n.index();
    else
      return S+N[n.index()].head;
  }
  inline const unsigned position(const st_index n) const {
    if (n.is_leaf())
      return n.index();
    else
      return N[n.index()].head;
  }
  inline st_index first_child(const st_index n) const {
    return N[n.index()].child;
  }
  inline st_index next_child(const st_index n) const {
    if (n.is_leaf())
      return L[n.index()].sib;
    else
      return N[n.index()].sib;
  }
  void resize(unsigned,unsigned);
  inline char term() const {return TERM;}

  void dump(ostream &o=cerr) const;
  inline void print() const {print(cout);}
  inline void print(ostream &o) const {print(o,root());}
  inline void print(ostream &o,st_index n) const {print(o,n,0,0);}
  void print(ostream &,st_index,unsigned,unsigned) const;
  inline void suffix_print() const {suffix_print(cout);}
  inline void suffix_print(ostream &o) const {suffix_print(o,root());}
  void suffix_print(ostream &,st_index) const;

  unsigned allpos(st_index,unsigned *,unsigned,unsigned=0) const;

  inline unsigned exists(const char *query) const {
    st_index nd = locate(query);
    if (!is_root(nd)) {
      return true;
    } 
    return false;
  }

  inline unsigned count(const char *query) const {
    st_index nd = locate(query);
    if (!is_root(nd)) {
      return allpos(nd,NULL,0);
    }
    return 0;
  }

  inline unsigned find(const char *query, unsigned * const pos, unsigned const poslen) const {
    st_index nd = locate(query);
    if (!is_root(nd)) {
      return allpos(nd,pos,poslen);
    }
    return 0;
  }

  inline st_index locate(const char *query) const {
    return locate(query,0,root());
  }

  st_index locate(const char *query, unsigned const offset, st_index n) const;

  bool check(ostream &) const;
  bool check_order(ostream &,st_index,unsigned &) const;

  private:

  enum extend_result {branch_left,branch_right,
		      branch_terminal,
		      new_child,new_terminal};

  struct extend_state {
    const char *S;
    unsigned depth;  // output: how many characters of extension match
    unsigned fastlen;// input : how many characters of extension will match
    unsigned len;    // input : maximum number of characters to extend
                     //         (implicitly S[i] != TERM for i < len)
    st_index current;// leaf or branch were 1) first mismatch occurred
                     //  of 2) match terminated
    st_index parent; // parent of current
    st_index brother;// older sibling of current

    inline extend_state(const char *st): S(st),depth(0),fastlen(0),len(0),
      current(st_index()),parent(st_index()),brother(st_index()) {}
  };

  extend_result extend(extend_state &) const;

  public:

};

template <class Node,class Suffix>
bool suffix_tree<Node,Suffix>::check_order(ostream &o,st_index n,
					   unsigned &last_leaf) const {
  if (n.is_leaf()) {
    unsigned h = n.index();
    if (last_leaf < Slen) {
      unsigned ll = last_leaf;
      last_leaf = h;
      for(unsigned i=0; S[ll+i] != TERM || S[h+i] != TERM; ++i)
	if (S[ll+i] < S[h+i])
	  return true;
	else if (S[ll+i] > S[h+i]) {
	  o << "leaves: " << last_leaf <<" and " 
	    << h << " out of order" << endl;
	  return false;
	}
    }
    else
      last_leaf = h;
    return true;
  }
  else {
    n = N[n.index()].child;
    n.set_good();
    while(n.is_good()) {
      if (!check_order(o,n,last_leaf)) return false;
      if (n.is_leaf())
	n = L[n.index()].sib;
      else
	n = N[n.index()].sib;
    }
    return true;
  }
}
template <class Node,class Suffix>
bool suffix_tree<Node,Suffix>::check(ostream &o) const {
  bool pass = true;
  unsigned i;
  o << "there are " << nN << " Nodes" << endl;
  o << "there are " << nL << " Leaves" << endl;
  for (unsigned x=0; x < nN; ++x) {
    unsigned h = N[x].head;
    unsigned l = N[x].len;
    st_index n = N[x].child;
    if (n.is_nogood()) {
      n.set_good();
      while(n.is_good()) {
	for(i=0; S[n.index()+i] != TERM; ++i)
	  if (S[h+i] != S[n.index()+i]) {
	    o << "mismatch of leaf with terminal node" << endl;
	    o << n.index() <<" and "<<h<<"/"<<l<<endl;
	    pass = false;
	  }
	if (i != l) {
	  o << "terminal node different length as child" << endl;
	  pass = false;
	}
	if (n.is_noleaf()) {
	  o << "child of terminal node " << x << " is non-leaf" << endl;
	  pass = false;
	}
	n = L[n.index()].sib;
      }
      
    }
    else {
      if (n.is_leaf()) {
	for(i=0; i < l; ++i)
	  if (S[h+i] != S[n.index()+i]) {
	    o << "mismatch of leaf with parent node" << endl;
	    pass = false;
	  }
      }
      else {
	for(i=0; i < l; ++i)
	  if (S[h+i] != S[N[n.index()].head+i]) {
	    o << "mismatch of node with parent node" << endl;
	    pass = false;
	  }
      }
      while(n.is_good()) {
	if (n.is_leaf())
	  n = L[n.index()].sib;
	else
	  n = N[n.index()].sib;
      }
      if (n.is_leaf()) {
	o << "suffix link should be nonleaf but it isn't" << endl;
	pass = false;
      }
      else if (N[n.index()].child.is_nogood()) {
	o << "suffix link "; n.print(o); o << " of " << x;
	o << " is not a non-terminal node" << endl;
	pass = false;
      }
      else { 
	if (l>1) {
	  if (l != N[n.index()].len+1) {
	    o << "suffix link is not parent length-1" << endl;
	    o << x << " ";
	    o << N[n.index()].head <<"/"<< N[n.index()].len << " by ";
	    o << h<<"/"<<l << endl;
	    pass = false;
	  }
	  for(i=1; i < l; ++i)
	    if (S[h+i]!=S[N[n.index()].head+i-1]) {
	      o << "mismatch between suffix link and referrent" << endl;
	      o << N[n.index()].head << "/" << N[n.index()].len << " by ";
	      o << h << "/" << l << endl;
	      pass = false;
	    }
	}
	
      }
    }
  }
  if (pass) {
    o << "All sanity checks passed" << endl;
    o << "Checking order" << endl;
    unsigned ll=Slen;
    if (check_order(o,root(),ll)) {
      o << "Order check passed" << endl;
      return true;
    }
    else {
      o << "Order check failed" << endl;
      return false;
    }
  }
  else
    return false;
}


template <class Node,class Suffix>
void suffix_tree<Node,Suffix>::print(ostream &o,st_index n,unsigned len,unsigned dep) const {
  if (n.is_leaf()) {
    for (unsigned i=len; S[n.index()+i]!=TERM; ++i)
      o << char(S[n.index()+i]);
    o << TERM << n.index() << endl;
  }
  else {
    for (unsigned i = len; i < N[n.index()].len ; ++i)
      o << char(S[N[n.index()].head+i]);
    len = N[n.index()].len;
    n = N[n.index()].child;
    if (n.is_good())
      o << '+';
    else
      o << '-';
    n.set_good();

    print(o,n,len,dep+1);

    if (n.is_leaf())
      n = L[n.index()].sib;
    else
      n = N[n.index()].sib;

    while(n.is_good()) {
      for (unsigned i=0; i < len+dep+1; ++i) o << " ";
      print(o,n,len,dep+1);
      if (n.is_leaf())
	n = L[n.index()].sib;
      else
	n = N[n.index()].sib;
    }
  }
}

template <class Node,class Suffix>
void suffix_tree<Node,Suffix>::suffix_print(ostream &o,st_index n) const {
  if (n.is_leaf()) {
    o << n.index() << "\t";
    for (unsigned i=0; S[n.index()+i]!=TERM && i < 10; ++i)
      o << char(S[n.index()+i]);
    o << endl;
  }
  else {
    n = N[n.index()].child;
    n.set_good();
    while(n.is_good()) {
      suffix_print(o,n);
      if (n.is_leaf())
	n = L[n.index()].sib;
      else
	n = N[n.index()].sib;
    }
  }
}

template <class Node,class Suffix>
unsigned suffix_tree<Node,Suffix>::allpos(st_index n, unsigned *pos, unsigned poslen, unsigned npos) const {
  if (n.is_leaf()) {
    if (npos < poslen) {
      pos[npos] = n.index();
    }
    ++npos;
    return npos;
  }
  n = N[n.index()].child;
  n.set_good();
  while(n.is_good()) {
    npos = allpos(n,pos,poslen,npos);
    if (n.is_leaf())
      n = L[n.index()].sib;
    else
      n = N[n.index()].sib;
  }
  return npos;
}

template <class Node,class Suffix>
st_index suffix_tree<Node,Suffix>::locate(const char * query, unsigned offset, st_index n) const {
  // fprintf(stderr, "%s:%d: In find w/ query = \"%s\"\n",__FILE__,__LINE__,query+offset);
  if (n.is_leaf()) {
    if (strncmp(str(n),query,strlen(query)) == 0) {
      // fprintf(stderr, "%s:%d: node is leaf, returning allpos of leaf\n",__FILE__,__LINE__);
      return n;
    } else {
      return root();
    }
  }
  if (query[offset] == '\0') {
    // fprintf(stderr, "%s:%d: No more query, returning allpos of non-leaf node.\n",__FILE__,__LINE__);
    return n;
  }
  unsigned depth0 = node(n).length();
  n = first_child(n);
  for(n.set_good(); n.is_good(); n = next_child(n)) {
    char edge_char = str(n)[depth0];
    // char buffer[21];
    // strncpy(buffer,str(n),20);
    // buffer[20] = '\0';
    // fprintf(stderr, "%s:%d: Child %u %s... %d string starts with %c, query starts with %c\n",
    // __FILE__,__LINE__,n.index(),buffer,n.is_leaf()?1:0,edge_char,query[offset]);
    if (edge_char == query[offset]) {
      if (n.is_leaf()) {
	return locate(query,offset,n);
      }
      unsigned depth1 = node(n).length();
      unsigned delta = depth1-depth0;
      // Want this child.
      // fprintf(stderr, "%s:%d: Child edge length is %d, query length is %d\n",__FILE__,__LINE__,delta,strlen(query)-offset);
      unsigned step = delta;
      if (strlen(query)-offset < delta) {
	step = strlen(query)-offset;
      }
      // fprintf(stderr, "node %u\n",n.index());
      if (strncmp(str(n)+depth0,query+offset,step)==0) {
	// fprintf(stderr, "%s:%d: Match to step characters.\n",__FILE__,__LINE__);
	return locate(query,offset+step,n);
      } else {
	// fprintf(stderr, "%s:%d: No match to step characters.\n",__FILE__,__LINE__);
	return root();
      }
    } else if (edge_char > query[offset]) {
      break;
    }
  }
  return root();
}


template <class Node,class Suffix>
void suffix_tree<Node,Suffix>::dump(ostream &o) const {
  o << "Nodes:" << endl;
  for(unsigned i=0; i < nN; ++i) {
    o << i << " " << N[i].head << "/" << N[i].len << " ";
    N[i].child.print(o);
    o << " ";
    N[i].sib.print(o);
    o << endl;
  }
  o << "Leaves:" << endl;
  for(unsigned i=0; i < nL; ++i) {
    o << i << " ";
    L[i].sib.print(o);
    o << endl;
  }
}

template <class Node,class Suffix>
void suffix_tree<Node,Suffix>::resize(unsigned nn,unsigned nl) {
  if (!Lsize) {
    if (nl)
      L = (Suffix*)malloc(nl*sizeof(Suffix));
  }
  else {
    if (nl)
      L = (Suffix*)realloc(L,nl*sizeof(Suffix));
    else
      free(L);
    //{printf("freeing L,%d\n",Lsize); free(L);}
  }
  if (!Nsize) {
    if (nn)
      N = (Node*)malloc(nn*sizeof(Node));
  }
  else {
    if (nn)
      N = (Node*)realloc(N,nn*sizeof(Node));
    else
      free(N);
    //{printf("freeing N,%d\n",Nsize); free(N);}
  }
  Nsize = nn;
  Lsize = nl;
}

template <class Node,class Suffix>
typename suffix_tree<Node,Suffix>::extend_result 
suffix_tree<Node,Suffix>::extend(extend_state &es) const {
  for(;;) {
#ifdef PRINT_CRAP
    es.current.print(cerr);
    cerr << es.depth << "." << N[es.current.index()].len << endl;
#endif
    if (es.current.is_noleaf()) {
      if (es.depth < es.fastlen) {
	es.depth = N[es.current.index()].len;
	if (es.depth > es.fastlen) es.depth = es.fastlen;
      }
      // match our way down the branch
      for( ; es.depth < es.len && es.depth < N[es.current.index()].len 
	    && es.S[es.depth]==S[N[es.current.index()].head+es.depth];
	  es.depth++) ;

      if (es.depth == N[es.current.index()].len) { // we match all chars in node
#ifdef PRINT_CRAP
	cerr << "match succeeds: looking at children" << endl;
#endif
	// child of this node
	st_index chd = N[es.current.index()].child;
	if (chd.is_nogood()) { // this node is terminal
	  if (es.depth == es.len) // the extension also terminal
	    return new_terminal;
	  else if (es.S[es.depth] < TERM) return branch_left;
	  else if (es.S[es.depth] > TERM) return branch_right;
	  else { // es.S[es.depth] == TERM, but es.depth < es.len
	    cerr << "suffix_tree::extend: invalid extension requested" << endl;
	    cerr << "  TERM encountered at " << es.depth;
	    cerr << "  when len is " << es.len << endl;
	    exit(51);
	  }
	}
	else { // this node is non-terminal, we can explore children
	  es.parent = es.current;
	  es.current.set_nogood(); // causes first brother to be no good
	  char head;

	  // since we haven't checked es.depth < es.len, it may be that
	  // es.S[es.depth] == TERM
	  char c;
	  if (es.depth < es.len)
	    c = es.S[es.depth];
	  else if (es.depth == es.len)
	    c = TERM;
	  else {
	    cerr << "suffix_tree::extend:somehow depth has exceeded len" << endl;
	    cerr << "  this should not happen." << endl;
	    exit(52);
	  }
	  do {
	    es.brother = es.current;
	    es.current = chd;
	    if (es.current.is_leaf()) {
	      head = S[es.current.index()+es.depth];
	      chd = L[es.current.index()].sib;
	    }
	    else {
	      head = S[N[es.current.index()].head+es.depth];
	      chd = N[es.current.index()].sib;
	    }
	  } while( head < c && chd.is_good() );
	  if (head == c) // follow child down
	    continue;
	  else if (head > c) // in between branches
	    return new_child;
	  else { // ran out of
	    // ok to take that extra step
	    es.brother = es.current;
	    es.current = chd;
	    return new_child;
	  }
	}
      }
      else { // either extension termination or mismatch
	const char head = S[N[es.current.index()].head+es.depth];
	if (es.depth == es.len) { // extension is terminated
	  if (TERM < head)  return branch_left;
	  else if (TERM > head) return branch_right;
	  else { // node has a terminal symbol on an internal branch
	    cerr << "suffix_tree::extend: Terminal symbol found on interior";
	    cerr << "of internal node " << es.current.index() << endl;
	    exit(53);
	  }
	}
	else { // mismatch
	  if (es.S[es.depth] < head) return branch_left;
	  else if (es.S[es.depth] > head) return branch_right;
	  else { // I guess loops don't work
	    cerr << "suffix_tree::extend: Ross is a bonehead #1" << endl;
	    cerr << "   should never get to here" << endl;
	    exit(54);
	  }
	}
      }
    }
    else { // is leaf, so just match down it
      if (es.depth < es.fastlen)
	es.depth = es.fastlen;
      
      for( ; es.depth < es.len
	     && es.S[es.depth]==S[es.current.index()+es.depth]; es.depth++) ;

      // extension termination or mismatch
      const char head = S[es.current.index()+es.depth];
      if (es.depth == es.len) { // extension is terminated
	if (TERM < head) return branch_left;
	else if (TERM > head) return branch_right;
	else  // leaf is terminal too, make new terminal node
	  return branch_terminal;
      }
      else { // we had just a simple mismatch
	if (es.S[es.depth] < head) return branch_left;
	else if (es.S[es.depth] > head) return branch_right;
	else { // I guess loops don't work
	  cerr << "suffix_tree::extend: Ross is a bonehead #2" << endl;
	  cerr << "   should never get to here" << endl;
	  exit(55);
	}
      }
    }
  } // for(;;)
}


// This big ugly function
template <class Node,class Suffix>
void suffix_tree<Node,Suffix>::add(const char *start,
				   unsigned length,
				   unsigned nsuf) {
  st_index needlink,link;
  unsigned scr_len;
  extend_result scr_ext;

  extend_state es(start);
  es.len = length;
  es.current = root();
  es.depth = 0;
  es.fastlen = 0;

  if (empty() && nsuf) {
    L[nL].sib.set(root()).set_nogood();

    // current tree is empty so we initialize a root
    N[nN].child.set_index(nL).set_good().set_leaf();
    N[nN].sib.set(root()).set_nogood();
    N[nN].len = 0; // forces a branch at the root
    N[nN].head = 0; // meaningless if len = 0;

    // this creates a new node, and leaf
    ++nL;
    ++nN;
    ++es.S;
    --nsuf;
    --es.len;
#ifdef PRINT_CRAP
    cerr << "Root and first leaf created" << endl;
#endif
  }

  // ok now we are ready for the rest...
  extend_result ext;
  for(; nsuf ; --nsuf, ++nL ) {
#ifdef PRINT_CRAP
    cerr << "here we go:" << nsuf << endl;
#endif
    L[nL].sib.set(root()).set_nogood().set_index(Slen);
    
#ifdef PRINT_CRAP
    dump(cerr);
    cerr << "(extend " << es.S - S << " from ";
#endif
    ext = extend(es);
#ifdef PRINT_CRAP
    cerr << ")" << endl;
#endif

    // add to tree
    switch(ext) {
    case new_child: case new_terminal:
      switch(ext) {
      case new_child:
#ifdef PRINT_CRAP
	cerr << "Missed branch" << endl;
#endif
	// es.current == place before which to insert
	// es.brother == place after which to insert
	L[nL].sib = es.current;
	if (es.brother.is_good()) {
	  // reroute from older brother
	  if (es.brother.is_leaf())
	    L[es.brother.index()].sib.set_index(nL).set_good().set_leaf();
	  else
	    N[es.brother.index()].sib.set_index(nL).set_good().set_leaf();
	}
	else { 
	  // reroute from parent
	  N[es.parent.index()].child.set_index(nL).set_good().set_leaf();
	}
	break;
      case new_terminal:
#ifdef PRINT_CRAP
	cerr << "Terminal child" << endl;
#endif
	// es.current must be a Node
	L[nL].sib.set(N[es.current.index()].child).set_good();
	N[es.current.index()].child.set_index(nL).set_nogood().set_leaf();
	break;
      default: // for -Wall
	cerr << "Should not reach here.  Ross doesn't know how switch";
	cerr << "statements work." << endl;
	exit(56);
      }
      // still in case new_child: case new_terminal:
      // follow suffix link
      if (es.depth) {
	while(es.current.is_good())
	  if (es.current.is_leaf())
	    es.current = L[es.current.index()].sib;
	  else
	    es.current = N[es.current.index()].sib;
	es.current.set_good();
      }
      else
	es.current = root();

      // standard post-conditions
      ++es.S;
      --es.len;
      es.fastlen = es.depth? es.depth-1 : 0;
      es.depth = N[es.current.index()].len;
      break;
    case branch_left: case branch_right: case branch_terminal:
      switch (ext) {
      case branch_left:
#ifdef PRINT_CRAP
	cerr << "branch to the left" << endl;
#endif
	N[nN].child.set_index(nL).set_good().set_leaf();
	if (es.current.is_leaf())
	  N[nN].sib = L[es.current.index()].sib;
	else
	  N[nN].sib = N[es.current.index()].sib;
	L[nL].sib = es.current;
	needlink = es.current;
	break;
      case branch_terminal:
#ifdef PRINT_CRAP
	cerr << "branch to terminal" << endl;
#endif
	N[nN].child.set_index(nL).set_nogood().set_leaf();
	if (es.current.is_leaf())
	  N[nN].sib = L[es.current.index()].sib;
	else
	  N[nN].sib = N[es.current.index()].sib;
	L[nL].sib = es.current;
	needlink = es.current;
	break;
      case branch_right:
#ifdef PRINT_CRAP
	cerr << "branch to the right" << endl;
	es.current.print(cerr); cerr << " is es current\n";
#endif
	N[nN].child = es.current;
	if (es.current.is_leaf()) {
	  N[nN].sib = L[es.current.index()].sib;
	  L[es.current.index()].sib.set_index(nL).set_good().set_leaf();
	}
	else {
	  N[nN].sib = N[es.current.index()].sib;
	  N[es.current.index()].sib.set_index(nL).set_good().set_leaf();
	}
#ifdef PRINT_CRAP
	cerr << "make sib: " ; N[nN].sib.print(cerr); cerr << endl;
#endif
	needlink.set_index(nL).set_good().set_leaf();
	break;
      default: // for -Wall
	cerr << "Should not reach here.  Ross doesn't know how switch";
	cerr << "statements work." << endl;
	exit(57);
      }
      // still in case branch_left: case branch_right: case branch_terminal:
      N[nN].head = es.S-S;
      N[nN].len  = es.depth;
      if (es.brother.is_good()) {
	// reroute from older brother
	if (es.brother.is_leaf())
	  L[es.brother.index()].sib.set_index(nN).set_good().set_noleaf();
	else
	  N[es.brother.index()].sib.set_index(nN).set_good().set_noleaf();
      }
      else
	// reroute from parent to new eldest
	N[es.parent.index()].child.set_index(nN).set_good().set_noleaf();

      // find suffix link for new node
      // short (depth==1,0) branches link to root, and branch_terminals
      // aren't really branches at all.
      if (es.depth<=1 || ext == branch_terminal) { 
#ifdef PRINT_CRAP
	cerr << "root suffix link" << endl;
#endif
	es.current = root();
	++es.S;
	es.fastlen = 0;
	--es.len;
	es.depth = N[es.current.index()].len; // ok, this is just 0
	link.set(root()).set_nogood();
      }
      else {
#ifdef PRINT_CRAP
	cerr << "find suffix link" << endl;
#endif
	// current node leads to parents slink
	es.current = N[nN].sib; 
	while(es.current.is_good())
	  if (es.current.is_leaf())
	    es.current = L[es.current.index()].sib;
	  else
	    es.current = N[es.current.index()].sib;
	es.current.set_good();
	++es.S;
	// setup hacks for fastmatch to suffix link
	scr_len = es.len-1; es.len = es.depth-1;
	es.fastlen = es.len;
	es.depth = N[es.current.index()].len;

	if (es.depth > es.len) {
	  cerr << "problem #1" << endl;
	  exit(58);
	}
	// if next round will cause a true branching then we will
	// predict the suffix as nN+1, else we will use the parent
	// of the match found here
	scr_ext = extend(es);

#ifdef PRINT_CRAP
	cerr << "extension returns " << scr_ext <<endl;
#endif
	// because we twiddled es.len, we are in cases where we
	// may be in *_terminal in which we want to really be in
	// the other one.

	// no branching will occur if parent is matched entirely
	// and branching will occur if parent not entirely matched
	if (scr_ext == branch_terminal || scr_ext == new_terminal) {
	  if (N[es.parent.index()].len < es.depth)
	    scr_ext = branch_terminal;
	  else
	    scr_ext = new_terminal;
	}
	switch(scr_ext) {
	case new_terminal: case new_child:
#ifdef PRINT_CRAP
	  cerr << "old suffix link works (restart find)" << endl;
#endif
	  // since we may be in an iteration over children which
	  // leaves es in an inconsistent state, we should back up
	  // to the parent
	  link.set(es.parent).set_nogood();
	  es.current = es.parent;
	  es.depth = N[es.current.index()].len;
	  break;
	case branch_left: case branch_right: case branch_terminal:
#ifdef PRINT_CRAP
	  cerr << "new suffix link needed" << endl;
#endif
	  if (nsuf>1) // if this is the last suffix to be added then don't
	    link.set_index(nN+1).set_nogood().set_noleaf();
	  else
	    link.set(root()).set_nogood();
	  break;
	}
	// undo hacks for suffix link
	es.fastlen = es.len;
	es.len = scr_len;
      }
      // needlink is the guy who needs a new sib
      if (needlink.is_leaf())
	L[needlink.index()].sib = link;
      else
	N[needlink.index()].sib = link;

      if (es.depth > es.len) {
	cerr << "problem #2" << endl;
	exit(59);
      }
      ++nN;
      break;
    }
    if (es.depth > es.len) {
      cerr << "problem #3" << endl;
      exit(60);
    }
#ifdef PRINT_CRAP
    print(cerr);
#endif
  }
}


template <class Node,class Suffix>
void suffix_tree<Node,Suffix>::build(unsigned part) {

  if (S[Slen-1] != TERM) {
    cerr << "suffix_tree::build: suffix tree string is not properly terminated" << endl;
    cerr << "final character at position " << Slen << " is " << S[Slen-1] <<" not "<< TERM << endl;
    exit(61);
  }

  resize(part,part);
  unsigned i=0,j=0;
  for (; i < part; i = j+1) {
    for(j=i; j < part-1 && S[j] != TERM; ++j);

    add(S+i,j-i,j-i+1);
  }
}

#endif
