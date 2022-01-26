
#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "memory_debug.h"

#include <iostream>
#include <list>
#include <map>
#include <vector>
#include "util.h"
#include "tinylist.t"
#include "sortedvector.t"

template<class E>
class basenode {
  std::list<E*> out_;
  std::list<E*> in_;
  bool mark_;
 public:
  basenode() : mark_(false) {};
  virtual ~basenode() {};
  std::list<E*> const & out() const { return out_; };
  std::list<E*> const & in() const { return in_; };
  std::list<E*> & out() { return out_; };
  std::list<E*> & in() { return in_; };
  long unsigned int nout() const {
    return out_.size();
  };
  long unsigned int nout_mark(bool mark) const {
    long unsigned int count=0;
    typename std::list<E*>::const_iterator eit(out_.begin());
    while (eit != out_.end()) {
      if ((*eit)->mark() == mark) {
	++count;
      }
      ++eit;
    }
    return count;
  }
  long unsigned int nin() const {
    return in_.size();
  }
  long unsigned int nin_mark(bool mark) const {
    long unsigned int count=0;
    typename std::list<E*>::const_iterator eit(in_.begin());
    while (eit != in_.end()) {
      if ((*eit)->mark() == mark) {
	++count;
      }
      ++eit;
    }
    return count;
  }
  void add_out(E* e) {
    out_.push_back(e);
  }
  void add_in(E* e) {
    in_.push_back(e);
  }
  void remove_out(E* e) {
    out_.remove(e);
  }
  void remove_in(E* e) {
    in_.remove(e);
  }
  typename std::list<E*>::iterator erase_out(typename std::list<E*>::iterator eit) {
    return out_.erase(eit);
  }
  typename std::list<E*>::iterator  erase_in(typename std::list<E*>::iterator eit) {
    return in_.erase(eit);
  }
  bool mark() const { return mark_; }
  void set_mark() { mark_ = true; }
  void unset_mark() { mark_ = false; }
  MEMORY_DEBUG(basenode)
};

template<class N>
class baseedge {
  N* from_;
  N* to_;
  bool mark_;
 public:
  baseedge(N* f=0, N* t=0) : from_(f), to_(t), mark_(false) {};
  virtual ~baseedge() {};
  N* const from() const { return from_; };
  N* const to() const { return to_; };
  void from(N* f) { from_ = f; };
  void to(N* t) { to_ = t; };
  bool mark() const { return mark_; }
  void set_mark() { mark_ = true; }
  void unset_mark() { mark_ = false; }
  MEMORY_DEBUG(baseedge)
};

template<class N>
class dfsconf {
  bool dir_;
 public:
  dfsconf() : dir_(false) {};
  virtual ~dfsconf() {};
  virtual void root_node(N*)=0;
  virtual bool process_node(N*)=0;
  virtual void init()=0;
  virtual void fini()=0;
  bool undirected() { return !dir_; }
  bool directed() { return dir_; }
  void set_undirected() { dir_ = false; }
  void set_directed() { dir_ = true; }
  MEMORY_DEBUG(dfsconf)
};

template<class N>
class components : public dfsconf<N> {
  sortedvector<long unsigned int,N*> c2n_;
  long unsigned int components_;
  bool first_lookup;
  typename sortedvector<long unsigned int,N*>::const_iterator it;
 public:
  class iterator {
    typename sortedvector<long unsigned int,N*>::const_iterator it_;
  public:
    iterator(typename sortedvector<long unsigned int,N*>::const_iterator const & it) : it_(it) {};
    iterator(iterator const & it) : it_(it.it_) {};
    N* operator*() {
      return (it_->value());
    }
    bool operator!=(const iterator & it) const {
      return (it_ != it.it_);
    }
    bool operator==(const iterator & it) const {
      return (it_ == it.it_);
    }
    iterator & operator++() {
      ++it_;
      return (*this); 
    }
  };
  components(long unsigned int nds) : components_(0) {
    dfsconf<N>::set_undirected();
    c2n_.reserve(nds);
    first_lookup = true;
  }
  void root_node(N*) {
    components_++;
  }
  bool process_node(N* n) {
    c2n_.push_back(components_,n);
    return true;
  }
  iterator nodes_begin(long unsigned int c) {
    try {
      if (first_lookup) {
	return (it = c2n_.locate_first_at_least(c));
	first_lookup = false;
      } else {
	return (it = c2n_.finger_locate_first_at_least(it,c));
      }
    } 
    catch (typename sortedvector<long unsigned int,N*>::KeyOutOfRange &) {
      return ((sortedvector<long unsigned int,N*> const &)c2n_).end();
    }
  }
  iterator nodes_end(long unsigned int c) {
    try {
      if (first_lookup) {
	return (it = c2n_.locate_first_at_least(c+1));
	first_lookup = false;
      } else {
	return (it = c2n_.finger_locate_first_at_least(it,c+1));
      }
    }
    catch (typename sortedvector<long unsigned int,N*>::KeyOutOfRange &) {
      return ((sortedvector<long unsigned int,N*> const &)c2n_).end();
    }
  }
  long unsigned int number() const {
    return components_;
  }
  long unsigned int number(long unsigned int c) {
    long unsigned int begin,end;
    try {
      if (first_lookup) {
	it = c2n_.locate_first_at_least(c);
	first_lookup = false;
      } else {
	it = c2n_.finger_locate_first_at_least(it,c);
      }
      begin = c2n_.index(it);
    } 
    catch (typename sortedvector<long unsigned int,N*>::KeyOutOfRange &) {
      begin = c2n_.index(c2n_.end());
    }
    try {
      it = c2n_.finger_locate_first_at_least(it,c+1);
      end = c2n_.index(it);
    } 
    catch (typename sortedvector<long unsigned int,N*>::KeyOutOfRange &) {
      end = c2n_.index(c2n_.end());
    }
    return (end-begin);
  }
  void init() {}
  void fini() {
    c2n_.normalize();
    // for (long unsigned int i=1;i<=number();i++) {
    // iterator itb(nodes_begin(i));
    // iterator ite(nodes_end(i));
    // std::cerr << i << ": ";
    // long unsigned int count=0;
    // while (itb!=ite) {
    // std::cerr << (*itb)->name() << " ";
    // count ++;
    // ++itb;
    // }
    // std::cerr << ": " << count << endl;
    // }
  }
  MEMORY_DEBUG(components)
};

template <class N>
class componentsl : public dfsconf<N> {
  std::vector<std::list<N*> > c2n_;
  long unsigned int components_;
 public:
  componentsl(long unsigned int) : components_(0) {
    dfsconf<N>::set_undirected();
  }
  void root_node(N*) {
    components_++;
  }
  bool process_node(N* n) {
    c2n_.resize(components_);
    c2n_[components_-1].push_back(n);
    return true;
  }
  typedef typename std::list<N*>::const_iterator iterator;
  iterator nodes_begin(long unsigned int c) const {
    return c2n_[c-1].begin();
  }
  iterator nodes_end(long unsigned int c) const {
    return c2n_[c-1].end();
  }
  long unsigned int number() const {
    return c2n_.size();
  }
  void init() {}
  void fini() {
    // for (long unsigned int i=1;i<=number();i++) {
    //   iterator itb(nodes_begin(i));
    //   iterator ite(nodes_end(i));
    //   std::cerr << i << ": ";
    //   long unsigned int count=0;
    //   while (itb!=ite) {
    //     std::cerr << (*itb)->name() << " ";
    //     count++;
    //     ++itb;
    //   }
    //   std::cerr << ": " << count << endl;;
    // }
  }
};

template<class N, class E>
class basegraph {
  std::list<N*> nodes_;
  std::list<E*> edges_;
 protected:
  void delete_marked_nodes() {
    typename std::list<N*>::iterator nit;
    nit = nodes_.begin(); 
    while (nit != nodes_.end()) {
      if ((*nit)->mark()) {
	delete (*nit);
	nit = nodes_.erase(nit);
      } else {
	++nit;
      }
    }
  }
  void delete_marked_edges() {
    typename std::list<E*>::iterator eit;
    eit = edges_.begin(); 
    while (eit != edges_.end()) {
      if ((*eit)->mark()) {
	delete (*eit);
	eit = edges_.erase(eit);
      } else {
	++eit;
      }
    }
  }  
 public:
  basegraph() {};
  ~basegraph() {
    clear();
  }
  void clear() {
    set_node_marks();
    set_edge_marks();
    delete_marked_nodes();
    delete_marked_edges();
  }
  void new_node(N * n) {
    nodes_.push_back(n);
  }
  void new_edge(E * e) {
    edges_.push_back(e);
    e->to()->in().push_back(e);
    e->from()->out().push_back(e);
  }
  std::list<N*> const & nodes() const {
    return nodes_;
  }
  std::list<E*> const & edges() const {
    return edges_;
  }
  std::list<N*> & nodes() {
    return nodes_;
  }
  std::list<E*> & edges() {
    return edges_;
  }
  long unsigned int nnode() const {
    return nodes_.size();
  }
  long unsigned int nedge() const {
    return edges_.size();
  }
  void clear_node_marks() {
    typename std::list<N*>::iterator nit;
    nit = nodes_.begin(); 
    while (nit != nodes_.end()) {
      (*nit)->unset_mark();
      ++nit;
    }
  }
  void clear_edge_marks() {
    typename std::list<E*>::iterator eit;
    eit = edges_.begin(); 
    while (eit != edges_.end()) {
      (*eit)->unset_mark();
      ++eit;
    }
  }
  void set_node_marks() {
    typename std::list<N*>::iterator nit;
    nit = nodes_.begin(); 
    while (nit != nodes_.end()) {
      (*nit)->set_mark();
      ++nit;
    }
  }
  void set_edge_marks() {
    typename std::list<E*>::iterator eit;
    eit = edges_.begin(); 
    while (eit != edges_.end()) {
      (*eit)->set_mark();
      ++eit;
    }
  }
  void remove(N* const & n) {
    clear_edge_marks();
    bool remove_edges = false;
    typename std::list<E*>::iterator eit;
    eit = n->out().begin(); 
    while (eit != n->out().end()) {
      remove_edges = true;
      (*eit)->to()->remove_in(*eit);
      (*eit)->set_mark();
      ++eit;
    }
    eit = n->in().begin(); 
    while (eit != n->in().end()) {
      remove_edges = true;
      (*eit)->from()->remove_out(*eit);
      (*eit)->set_mark();
      ++eit;
    }
    if (remove_edges) {
      delete_marked_edges();
    }
    n->set_mark();
    delete_marked_nodes();
  }
  /* void remove(E* const & e) {
    delete (*eit);
    edges_.remove(*eit);
    }*/
  void remove_out_edges(N* const & n) {
    typename std::list<E*>::iterator eit;
    eit = n->out().begin(); 
    while (eit != n->out().end()) {
      (*eit)->to()->remove_in(*eit);
      n->remove_out(*eit);
      delete (*eit);
      eit = n->out().erase(eit);
    }
  }
  void mark_in_edges(N* const & n) {
    typename std::list<E*>::iterator eit;
    eit = n->in().begin(); 
    while (eit != n->in().end()) {
      (*eit)->set_mark();
      ++eit;
    }
  }
  void mark_out_edges(N* const & n) {
    typename std::list<E*>::iterator eit;
    eit = n->out().begin(); 
    while (eit != n->out().end()) {
      (*eit)->set_mark();
      ++eit;
    }
  }
  void remove_in_edges(N* const & n) {
    typename std::list<E*>::iterator eit;
    eit = n->in().begin(); 
    while (eit != n->in().end()) {
      (*eit)->from()->remove_out(*eit);
      n->remove_in(*eit);
      edges_.remove(*eit);
      delete (*eit);
      eit = n->in().erase(eit);
    }
  }
  void remove(std::list<N*> const & nl) {
    clear_edge_marks();
    clear_node_marks();
    bool remove_edges = false;
    typename std::list<N*>::const_iterator nit;
    nit = nl.begin(); 
    while (nit != nl.end()) {
      typename std::list<E*>::iterator eit;
      eit = (*nit)->out().begin(); 
      while (eit != (*nit)->out().end()) {
	remove_edges = true;
	(*eit)->to()->remove_in(*eit);
	(*eit)->set_mark();
	++eit;
      }
      eit = (*nit)->in().begin(); 
      while (eit != (*nit)->in().end()) {
	remove_edges = true;
	(*eit)->from()->remove_out(*eit);
	(*eit)->set_mark();
	++eit;
      }
      (*nit)->set_mark();
      ++nit;
    }
    if (remove_edges) {
      delete_marked_edges();
    }
    delete_marked_nodes();
  }
  void remove_marked_nodes() {
    clear_edge_marks();
    bool remove_edges = false;
    typename std::list<N*>::const_iterator nit;
    nit = nodes().begin(); 
    while (nit != nodes().end()) {
      typename std::list<E*>::iterator eit;
      eit = (*nit)->out().begin(); 
      while (eit != (*nit)->out().end()) {
	remove_edges = true;
	(*eit)->to()->remove_in(*eit);
	(*eit)->set_mark();
	++eit;
      }
      eit = (*nit)->in().begin(); 
      while (eit != (*nit)->in().end()) {
	remove_edges = true;
	(*eit)->from()->remove_out(*eit);
	(*eit)->set_mark();
	++eit;
      }
      ++nit;
    }
    if (remove_edges) {
      delete_marked_edges();
    }
    delete_marked_nodes();
  }
  void remove_marked_edges() {
    typename std::list<E*>::const_iterator eit;
    eit = edges().begin(); 
    while (eit != edges().end()) {
      if ((*eit)->mark()) {
	(*eit)->to()->remove_in(*eit);
	(*eit)->from()->remove_out(*eit);
      }
      ++eit;
    }
    delete_marked_edges();
  }
  void remove(std::list<E*> const & el) {
    clear_edge_marks();
    typename std::list<E*>::iterator eit;
    eit = el.begin(); 
    while (eit != el.end()) {
      (*eit)->to()->remove_in(*eit);
      (*eit)->from()->remove_out(*eit);
      (*eit)->set_mark();
    }
    delete_marked_edges();
  }
  typename std::list<N*>::iterator erase(typename std::list<N*>::iterator const & nit) {
    clear_edge_marks();
    bool remove_edges = false;
    typename std::list<E*>::iterator eit;
    eit = (*nit)->out().begin(); 
    while (eit != (*nit)->out().end()) {
      remove_edges = true;
      (*eit)->to()->remove_in(*eit);
      (*eit)->set_mark();
      ++eit;
    }
    eit = (*nit)->in().begin(); 
    while (eit != (*nit)->in().end()) {
      remove_edges = true;
      (*eit)->from()->remove_out(*eit);
      (*eit)->set_mark();
      ++eit;
    }
    if (remove_edges) {
      delete_marked_edges();
    }
    delete (*nit);
    return nodes_.erase(nit);
  }
  typename std::list<E*>::iterator erase(typename std::list<E*>::iterator const & eit) {
    (*eit)->to()->remove_in(*eit);
    (*eit)->from()->remove_out(*eit);
    delete (*eit);
    return edges_.erase(eit);
  }
  void move_edges(N* from, N* to) {
    typename std::list<E*>::iterator eit = from->in().begin();
    while (eit !=  from->in().end()) {
      (*eit)->to(to);
      to->add_in(*eit);
      eit = from->erase_in(eit);
    }
    eit = from->out().begin();
    while (eit !=  from->out().end()) {
      (*eit)->from(to);
      to->add_out(*eit);
      eit = from->erase_out(eit);
    }
  }
  E* find_one(N* f, N* t) const {
    if (f->nout() < t->nin()) {
      typename std::list<E*>::iterator eit = f->out().begin();
      typename std::list<E*>::iterator eit1 = f->out().end();
      while (eit != eit1) {
	if ((*eit)->to() == t) {
	  return (*eit);
	}
	++eit;
      }
      return 0;
    } else {
      typename std::list<E*>::iterator eit = t->in().begin();
      typename std::list<E*>::iterator eit1 = t->in().end();
      while (eit != eit1) {
	if ((*eit)->from() == f) {
	  return (*eit);
	}
	++eit;
      }
      return 0;
    }
  }
  template <class cmp>
  void sort_nodes(cmp f) {
    nodes_.sort(f);
  }
  void dfs(N* root, dfsconf<N> & dfsc) {
    if (this->nnodes() == 0 || root == 0) return;
    clear_node_marks();
    std::list<N*> stack;
    dfsc.root_node(root);
    stack.push_back(root);
    root->set_mark();
    while (!stack.empty()) {
      N* n = stack.front();
      stack.pop_front();
      if (!dfsc.process_node(n)) continue;
      typename std::list<E*>::iterator eit = n->out().begin();
      while (eit != n->out().end()) {
	if (!(*eit)->to()->mark()) {
	  stack.push_back((*eit)->to());
	  (*eit)->to()->set_mark();
	}
	++eit;
      }
      if (dfsc.undirected()) {
	eit = n->in().begin();
	while (eit != n->in().end()) {
	  if (!(*eit)->from()->mark()) {
	    stack.push_back((*eit)->from());
	    (*eit)->from()->set_mark();	    
	  }
	  ++eit;
	}
      }
    }
  }
  void dfs(dfsconf<N> & dfsc) {
    dfsc.init();
    if (nnode() == 0) {
      dfsc.fini();
      return;
    }
    clear_node_marks();
    std::list<N*> stack;
    typename std::list<N*>::iterator nit = nodes().begin();
    while (nit != nodes().end()) {
      N* n = *nit;
      if (!n->mark()) {
	dfsc.root_node(n);
	stack.push_back(n);
	n->set_mark();
	while (!stack.empty()) {
	  n = stack.front();
	  stack.pop_front();
	  if (!dfsc.process_node(n)) continue;
	  typename std::list<E*>::iterator eit = n->out().begin();
	  while (eit != n->out().end()) {
	    if (!(*eit)->to()->mark()) {
	      stack.push_back((*eit)->to());
	      (*eit)->to()->set_mark();
	    }
	    ++eit;
	  }
	  if (dfsc.undirected()) {
	    eit = n->in().begin();
	    while (eit != n->in().end()) {
	      if (!(*eit)->from()->mark()) {
		stack.push_back((*eit)->from());
		(*eit)->from()->set_mark();	    
	      }
	      ++eit;
	    }
	  }
	}
      }
      ++nit;
    }
    dfsc.fini();
  }
  MEMORY_DEBUG(basegraph)
};

template <class E>
class labelnode : public basenode<E> {
  long unsigned int name_;
 public:
  labelnode() : name_(0) {};
  virtual ~labelnode() {};
  long unsigned int name() const {
    return name_;
  }
  void name(long unsigned int n) {
    name_ = n;
  }
  MEMORY_DEBUG(labelnode)
};

template <class N>
class labeledge : public baseedge<N> {
 public:
  virtual ~labeledge() {};
  MEMORY_DEBUG(labeledge)
};

template <class N, class E>
class labelgraph : public basegraph<N,E> {
  typedef std::map<long unsigned int, N*> nodemap;
  nodemap nn2node_;
  long unsigned int maxlabel_;
 public:
  labelgraph() : maxlabel_(0) {};
  void new_node(N * n) {
    typename std::pair<typename nodemap::iterator,bool> res;
    res = nn2node_.insert(std::make_pair(n->name(),n));
    if (!res.second) {
      checkpoint;
      std::cerr << "labelgraph::new_node: Duplicate node name " 
		<< n->name() << "!" << endl;
      abort();
    } 
    if (maxlabel_ <= n->name()) {
      maxlabel_ = n->name()+1;
    }
    basegraph<N,E>::new_node(n);
  }
  N* find(long unsigned int name) {
    typename std::pair<typename nodemap::iterator,bool> res;
    typename nodemap::iterator fit;
    if ((fit=nn2node_.find(name))!=nn2node_.end()) {
      return (*fit).second;
    } else {
      return 0;
    }
  }
  void free_labelmap() {
    nn2node_.clear();
  }
  long unsigned int new_label() {
    return maxlabel_++;
  }
  void relabel() {
    nn2node_.clear();
    maxlabel_=0;
    typename std::list<N*>::iterator nit;
    nit = this->nodes().begin();
    while (nit != this->nodes().end()) {
      (*nit)->name(new_label());
      nn2node_.insert(std::make_pair((*nit)->name(),(*nit)));
      ++nit;
    }
  }
  MEMORY_DEBUG(labelgraph)
};

template <class N, class E>
class intlabelgraph : public basegraph<N,E> {
  typedef std::vector<N*> nodemap;
  nodemap *nn2node_;
  long unsigned int nn2node_size_;
  long unsigned int maxlabel_;
 public:
  intlabelgraph(long unsigned int maxlabelhint=10000) : maxlabel_(0), nn2node_(0) {
    nn2node_ = new std::vector<N*>(maxlabelhint,(N*)0);
    nn2node_size_ = nn2node_->size();
  };
  void new_node(N * n) {
    if (nn2node_) {
      while (n->name() >= nn2node_->size()) {
	nn2node_->resize(nn2node_->size()*2+1000,(N*)0);
	nn2node_size_ = nn2node_->size();
      }
      if ((*nn2node_)[n->name()]) {
	checkpoint;
	std::cerr << "intlabelgraph::new_node: Duplicate node name " 
		  << n->name() << "!" << endl;
	abort();
      } else {
	(*nn2node_)[n->name()] = n;
      }
    }
    if (maxlabel_ <= n->name()) {
      maxlabel_ = n->name()+1;
    }
    basegraph<N,E>::new_node(n);
  }
  N* find(long unsigned int name) {
    if (nn2node_ && name < nn2node_size_) {
      return (*nn2node_)[name];
    } else {
      return 0;
    }
  }
  long unsigned int new_label() {
    return maxlabel_++;
  }
  void relabel() {
    if (!nn2node_) {
      nn2node_ = new std::vector<N*>();
    }
    nn2node_->resize(this->nnodes());
    nn2node_->init((N*)0);
    maxlabel_=0;
    typename std::list<N*>::iterator nit;
    nit = this->nodes().begin();
    while (nit != this->nodes().end()) {
      (*nit)->name(new_label());
      nn2node_[(*nit)->name()] = (*nit);
      ++nit;
    }
  }
  void free_labelmap() {
    delete nn2node_;
    nn2node_ = 0;
  }
  long unsigned int labelmapsize() {
    if (nn2node_) {
      return nn2node_->capacity();
    } else {
      return 0;
    }
  }
  MEMORY_DEBUG(intlabelgraph)
};

#endif

