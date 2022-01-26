
// -*- C++ -*- //

#ifndef _KEYWORD_TREE_T
#define _KEYWORD_TREE_T

#include "keyword_tree.h"

template <class KTND>
ktnode<KTND>::ktnode() : fail_(0), output_(0), patid_(0) {};

template <class KTND>
ktnode<KTND>::~ktnode() {
  clear_patid();
}

template <class KTND>
void ktnode<KTND>::print_patid() {
  if (patid_) {
    tinylist<pattern_list::const_iterator>::iterator pid;
    for (pid=patid_->begin();pid!=patid_->end();++pid) {
      cerr << (*pid)->id() << ":" << (*pid)->pattern() << " ";
    }
    cerr << endl;
  } else {
    cerr << "pattern list is empty" << endl;
  }
}

template <class KTND>
void ktnode<KTND>::clear_patid() {
  //checkpoint;
  if (patid_) delete patid_;
}

template <class KTND>
KTND *  ktnode<KTND>::fail() const {
  return fail_;
}

template <class KTND>
void ktnode<KTND>::fail(KTND * const v) {
  fail_ = v;
}

template <class KTND>
KTND *  ktnode<KTND>::output() const {
  return output_;
}

template <class KTND>
void ktnode<KTND>::output(KTND * const v) {
  output_= v;
}

template <class KTND>
tinylist<pattern_list::const_iterator> * const & ktnode<KTND>::patid() const {
  return patid_;
}

template <class KTND>
void ktnode<KTND>::add_patid(pattern_list::const_iterator const & it) {
  if (!patid_) patid_ = new tinylist<pattern_list::const_iterator>;
  patid_->push_front(it);
}

template <class KTND>
void ktnode<KTND>::del_patid(pattern_list::const_iterator const & it) {
  // checkpoint;
  if (patid_) patid_->remove(it);
  if (patid_->empty()) {
    delete patid_;
    patid_ = 0;
  }
}

template<class KTND>
void keyword_tree<KTND>::add_keyword_(KTND * v, char const * const pat, unsigned int len,
                                      pattern_list::const_iterator  const & it) {
  assert(v!=0);
  assert(pat!=0);
  
  KTND *c=0;
  // cerr << __FILE__ << ":" << __LINE__ << " " << (void*) v << " " << (void*)pat << " " << pat << " " << (void*) li << endl;  
  if (len == 0) {
    // current node is end of keyword...
    // cerr << __FILE__ << ":" << __LINE__ << " " << (void*) v->patid() << endl;  
    // assert(!(v->patid().empty()));
    v->add_patid(it);
    // cerr << __FILE__ << ":" << __LINE__ << endl;  
  } else if ((c=v->lchild(*pat))) {
    relchars_[(unsigned char)(*pat)] = true;
    // cerr << __FILE__ << ":" << __LINE__ << endl;  
    add_keyword_(c,(pat+1),len-1,it);
    // cerr << __FILE__ << ":" << __LINE__ << endl;  
  } else {
    relchars_[(unsigned char)(*pat)] = true;
    // cerr << __FILE__ << ":" << __LINE__ << endl;  
    c = v->add_child(*this,*pat);
    // cerr << __FILE__ << ":" << __LINE__ << " " << (void*) c << endl;  
    add_keyword_(c,(pat+1),len-1,it);
    // cerr << __FILE__ << ":" << __LINE__ << endl;  
  }
  // cerr << __FILE__ << ":" << __LINE__ << endl;  
}

template<class KTND>
void keyword_tree<KTND>::del_keyword_(KTND * v, char const * const pat, int len,
				      long unsigned int id) {
  // checkpoint;
  assert(v!=0);
  assert(pat!=0);
  
  KTND *c=0;
  if (len == 0) {
    // checkpoint;
    if (id == 0) {
      v->clear_patid();
    } else {
      tinylist<pattern_list::const_iterator>::iterator it,eit;
      if (v->patid()) {
	it = v->patid()->begin();
	eit = v->patid()->end();
	while (it != eit) {
	  if ((*it)->id() == id) {
	    v->del_patid((*it));
	  }
	  ++it;
	}
      }
    }
  } else if ((c=v->child(*pat))) {
    del_keyword_(c,(pat+1),len-1,id);
  } else {
    // checkpoint;
  }
}

template<class KTND>
tinylist<pattern_list::const_iterator> * 
keyword_tree<KTND>::lookup_keyword_(KTND * v, char const * const pat, int len) {
  // checkpoint;
  assert(v!=0);
  assert(pat!=0);
  
  KTND *c=0;
  if (len == 0) {
    // checkpoint;
    return v->patid();
  } else if ((c=v->child(*pat))) {
    return lookup_keyword_(c,(pat+1),len-1);
  } else {
    return 0;
  }
}

template <class KTND>
keyword_tree<KTND>::keyword_tree() { 
  root_ = newktnode();
  w_ = root_;
  first_char_ = true;
  relchars_ = new bool[256];
  for (int i=0;i<256;i++) {
    relchars_[i] = false;
  }
};

template <class KTND>
void keyword_tree<KTND>::reset() {
  w_ = root_;
  first_char_ = true;
}

template <class KTND>
keyword_tree<KTND>::~keyword_tree() {
  kid_jump_table<KTND> *kjt;
  typename std::list<kid_jump_table<KTND>*>::iterator it1;
  for (it1=jump_table_delete_list.begin(); 
       it1!=jump_table_delete_list.end(); ++it1) {
    delete [] *it1;
  }
  typename std::list<KTND*>::iterator it2;
  for (it2=ktnode_delete_list.begin();it2!=ktnode_delete_list.end();++it2) {
    delete [] *it2;
  }
  delete [] relchars_;
}

template <class KTND>
void keyword_tree<KTND>::init(CharacterProducer & cp) {

  pattern_list::const_iterator it;
  char *nchpat=0;
  int nchpatlen = 0;
  for (it=patterns().begin();it!=patterns().end();++it) {
    const string & pat = it->pattern();
    int pl = pat.length();
    if (!nchpat || ((pl+1) > nchpatlen)) {
      if (nchpat) {
	delete [] nchpat;
      }
      nchpatlen = pl*2+1;
      nchpat = new char[nchpatlen];
    }
    for (int i=0;i<pl;i++) {
      nchpat[i] = cp.nch(pat[i]);
    }
    nchpat[pl] = '\0';
    add_keyword_(root_,nchpat,pl,it);
  }
  delete [] nchpat;

  compute_failure_links();
  // optimize_failure_links();
  optimize_nodes(cp);
  compress_relchars(cp);
}


template <class KTND>
unsigned long keyword_tree<KTND>::add_pattern(std::string const & pat, 
					      long unsigned int id,
					      int esb, int eeb) {
  // checkpoint;
  // cerr << "add pattern: " << pat << endl;
  pattern_list::const_iterator it;
  it = add_pattern_(pat,id,esb,eeb);
  // add_keyword_(root_,pat.c_str(),pat.length(),it);
  return id;
}

template <class KTND>
void keyword_tree<KTND>::lookup_pattern(CharacterProducer & cp, 
					std::string const & pat,
					std::list<long unsigned int> & l0) {
  // checkpoint;
  // cerr << "lookup pattern: " << pat << endl;
  char *nchpat = new char[pat.length()+1];
  for (int i=0;i<pat.length();i++) {
    nchpat[i] = cp.nch(pat[i]);
  }
  nchpat[pat.length()] = '\0';
  tinylist<pattern_list::const_iterator> * l = 
    lookup_keyword_(root_,nchpat,pat.length());
  tinylist<pattern_list::const_iterator>::iterator it;
  if (l) {
    // checkpoint;
    it = l->begin();
    while (it != l->end()) {
      // checkpoint;
      // (*it)->write(cerr);
      // cerr << endl;
      l0.push_back((*it)->id());
      ++it;
    }
  }
  delete [] nchpat;
}

template <class KTND>
void keyword_tree<KTND>::del_pattern(CharacterProducer & cp, std::string const & pat,
				     long unsigned int id) {
  // checkpoint;
  char *nchpat = new char[pat.length()+1];
  for (int i=0;i<pat.length();i++) {
    nchpat[i] = cp.nch(pat[i]);
  }
  nchpat[pat.length()] = '\0';
  del_keyword_(root_,nchpat,pat.length(),id);
  delete [] nchpat;
}

template <class KTND>
KTND *keyword_tree<KTND>::newktnode() {
  if (ktnodes.empty()) {
    int size=(GETPAGESIZE-128)/sizeof(KTND);
    if (size < 1) {
      size = 1;
    }
    KTND *buffer = new KTND[size];
    for (int i=0;i<size;i++) {
      ktnodes.push_front(buffer+i);
    }
    ktnode_delete_list.push_front(buffer);
  }
  KTND *retval = ktnodes.front();
  ktnodes.pop_front();
  return retval;
}

template <class KTND>
kid_jump_table<KTND> *keyword_tree<KTND>::newktjumptable(unsigned int size) {
  unsigned int number=(GETPAGESIZE-128)/(size*sizeof(kid_jump_table<KTND>));
  if (number < 1) {
    number = 1;
  }
  kid_jump_table<KTND> *buffer= new kid_jump_table<KTND>[size*number];
  for (unsigned int i=0;i<size*number;i++) {
    buffer[i].which=0;
  }
  for (unsigned int i=0;i<number;i++) {
    jump_tables_.push_front(buffer+i*size);
  }
  jump_table_delete_list.push_front(buffer);
  return buffer;
}

template <class KTND>
void keyword_tree<KTND>::compute_failure_links() {
  root_->fail(root_);
  std::list<KTND*> q;
  q.push_back(root_);
  while (!q.empty()) {
    // checkpoint;
    failure_links_(q);
  }
}

template <class KTND>
void keyword_tree<KTND>::optimize_failure_links() {
  root_->fail(root_);
  std::list<KTND*> q;
  q.push_back(root_);
  while (!q.empty()) {
    optimize_failure_links_(q);
  }
}

template <class KTND>
void keyword_tree<KTND>::optimize_nodes(CharacterProducer & cp) {
  // checkpoint;
  root_->optimize_node(*this,cp);
  // checkpoint;
}

template <class KTND>
void keyword_tree<KTND>::compress_relchars(CharacterProducer & cp) {
  bool *tmp;
  int size=cp.size();
  tmp = new bool[size];
  for (int i=0;i<size;i++) {
    if (relchars_[i]) {
      tmp[i] = true;
      // checkpoint;
      // std::cerr << "Relchars: " << i << ": true" << endl;
    } else {
      tmp[i] = false;
    }
  }
  delete [] relchars_;
  relchars_ = tmp;
}

template<class KTND>
bool keyword_tree<KTND>::has_novel_lchild(KTND *n1, KTND *n0) const {
  tinylist<kid_list_element<KTND> > *l = n1->kids();
  typename tinylist<kid_list_element<KTND> >::iterator it;
  if (l) {
    for (it=l->begin();it!=l->end();++it) {
      if (!n0 || !n0->lchild(it->ch)) {
	return true;
      }
    }
  }
  return false;
}

template<class KTND>
void keyword_tree<KTND>::failure_links_(std::list<KTND*> & q) {
  KTND *vp;
  vp=q.front();
  q.pop_front();
  tinylist<kid_list_element<KTND> > *l = vp->kids();
  typename tinylist<kid_list_element<KTND> >::iterator it;
  if (l) {
    for (it=l->begin();it!=l->end();++it) {
      KTND *v = it->nd;
      assert(v!=0);
      unsigned char x = it->ch;
      KTND* w = vp->fail();
      assert(w!=0);
      while (!w->lchild(x) && 
	     w != root_) {
	w = w->fail();
	if (!w) break;
      }
      KTND* u;
      if (w && (u=w->lchild(x)) && vp != root_) {
	v->fail(u);
	if (u->patid()) {
	  v->output(u);
	} else if (u->output()){
	  v->output(u->output());
	}
      } else {
	v->fail(root_);
      }
    }
    for (it=l->begin();it!=l->end();++it) {
      KTND *v = it->nd;
      q.push_back(v);
    } 
  }
}

template<class KTND>
void keyword_tree<KTND>::optimize_failure_links_(std::list<KTND*> & q) {
  KTND *v;
  v=q.front();
  q.pop_front();
  KTND *u = v->fail();
  while (u != root_ && u && !has_novel_lchild(u,v)) {
    u = u->fail();
  }
  v->fail(u);
  tinylist<kid_list_element<KTND> > *l = v->kids();
  typename tinylist<kid_list_element<KTND> >::iterator it;
  if (l) {
    for (it=l->begin();it!=l->end();++it) {
      KTND *v = it->nd;
      q.push_back(v);
    } 
  }
}

template <class KTND>
bool keyword_tree<KTND>::find_patterns(CharacterProducer & cp, 
				       pattern_hit_vector & kas,
				       long unsigned minka) {
  register KTND *wp, *wpp;
  long unsigned kacount=0;
  register bool eof;
  tinylist<pattern_list::const_iterator>::iterator pid,pide;
  if ((eof=cp.eof())) return false;
  if (first_char_) {
    ch_ = cp.getnch();
    first_char_ = false;
  }
  // checkpoint; 
  // std::cerr << (void*)w_ << std::endl;
  while (!eof) {
    // checkpoint; 
    while (relchars_[ch_] && (wp=w_->child(ch_))) {
      // checkpoint; 
      if (wp->patid()) {
	pide = wp->patid()->end();
	for (pid=wp->patid()->begin();pid!=pide;++pid) {
	  // checkpoint;
	  // std::cerr << cp.pos() << " " << (*pid)->id() << std::endl;
	  kas.push_back(cp.pos(),make_pair(*pid,0));
	  kacount++;
	} 
      }
      wpp = wp->output();
      while(wpp && wpp->patid()) {
	pide = wp->patid()->end();
	for (pid=wpp->patid()->begin();pid!=pide;++pid) {
	  // checkpoint;
	  // std::cerr << cp.pos() << " " << (*pid)->id() << std::endl;
	  kas.push_back(cp.pos(),make_pair(*pid,0));
	  kacount++;
	} 
	wpp = wpp->output();
      }
      w_ = wp;
      if ((eof=cp.eof())) break;
      ch_ = cp.getnch();
    }
    if (eof) break;
    if (w_ == root_) {
      if ((eof=cp.eof())) break;
      ch_ = cp.getnch();
    } else {
      w_ = w_->fail();
      // checkpoint; 
      // std::cerr << (void*)w_ << std::endl;
    }
    if (kacount >= minka) {
      report_progress(cp);
      return true;
    }
  } 
  report_progress(cp);
  if (kacount > 0) return true;
  return false;
}

template <class KTND>
bool keyword_tree<KTND>::find_suffixes(CharacterProducer & cp, 
				       pattern_hit_vector & kas,
				       long unsigned int minlevel) {
  // checkpoint;
  register KTND *wp, *wpp;
  register bool eof;
  tinylist<pattern_list::const_iterator>::iterator pid,pide;
  if ((eof=cp.eof())) return false;
  if (first_char_) {
    ch_ = cp.getnch();
    // cerr << "Got " << cp.pos() << " " << ch_ << endl;
    first_char_ = false;
  }
  while (!eof) {
    while (relchars_[ch_] && (wp=w_->child(ch_))) {
      // cerr << "Match " << ch_ << endl;
      w_ = wp;
      // cerr << w_->level() << endl;
      if ((eof=cp.eof())) break;
      ch_ = cp.getnch();
      // cerr << "Got " << cp.pos() << " " << ch_ << endl;
    }
    if (eof) break;
    if (w_ == root_) {
      if ((eof=cp.eof())) break;
      // cerr << "At root" << endl;
      ch_ = cp.getnch();
      // cerr << "Got " << cp.pos() << " " << ch_ << endl;
    } else {
      // cerr << "Fail " << endl;
      w_ = w_->fail();
      // cerr << w_->level() << endl;
    }
  } 
  if (w_ == root_ || w_->level() < minlevel) {
    // checkpoint;
    return false;    
  } else {
    // We want a list of all the suffixes supported by this input. 
    // the current node, w_, represents the longest suffix of our input 
    // checkpoint;
    // cerr << w_->level() << endl;    
    while (w_ != root_ && w_->level() >= minlevel) {
      tinylist<KTND *> nodestack;
      nodestack.push_front(w_);
      // checkpoint;
      while (!nodestack.empty()) {
	// checkpoint;
	wp = *(nodestack.begin());
	nodestack.pop(); 
	// cerr << wp->level() << endl;
	if (wp->patid()) {
	  pide = wp->patid()->end();
	  for (pid=wp->patid()->begin();pid!=pide;++pid) {
	    kas.push_back(w_->level(),*pid);
	    // checkpoint;
	    // cerr << w_->level() << " " << (*pid)->pattern() << endl;
	  } 
	}
	// checkpoint;
	if (wp->kids() && !wp->kids()->empty()) {
	  // checkpoint;
	  tinylist<kid_list_element<KTND> > *l = wp->kids();
	  typename tinylist<kid_list_element<KTND> >::iterator it;
	  it = l->begin();
	  while (it != l->end()) {
	    // cerr << "Added " << it->ch << endl;
	    nodestack.push_front(it->nd);
	    ++it;
	  }
	}
      }
      w_ = w_->fail();
    }
    kas.normalize_strict_byvalue();

    pattern_hit_vector::iterator phlit,phlit0;
    long unsigned int removed=0;
    phlit = kas.begin();
    phlit0 = phlit;
    if (phlit != kas.end()) ++phlit;
    while (phlit != kas.end()) {
      if (phlit0->value() == phlit->value()) {
	phlit0->key() = MAXINT;
	removed++;
      }
      phlit0 = phlit;
      ++phlit;
    }
    kas.normalize();
    kas.resize(kas.size()-removed);
    // checkpoint;
    return true;
  }
}

#endif
