
#include "keyword_tree.h"
#include "keyword_tree.t"
#include <string.h>
#include <assert.h>
#include "char_io.h"
#include "util.h"
#include <unistd.h>

ktnode_list::ktnode_list() : kidlist_(0) {};

ktnode_dna_list::ktnode_dna_list() : kidlist_(0) {
  acgt_[0] = 0;  acgt_[1] = 0;  acgt_[2] = 0;  acgt_[3] = 0;
}

ktnode_jtable::ktnode_jtable() : kidptr_(0) {};

ktnode_suffix::ktnode_suffix() : kidlist_(0), kidjtable_(0), level_(0) {};

ktnode_lsuffix::ktnode_lsuffix() : kidlist_(0), level_(0) {};

ktnode_list::~ktnode_list() {
  if (kidlist_) {
    delete kidlist_;
  }
};
ktnode_dna_list::~ktnode_dna_list() {
  if (kidlist_) {
    delete kidlist_;
  }
};

ktnode_jtable::~ktnode_jtable() {
  // kidptr points to things that don't need to be deleted, as long as
  // optimize_nodes has been called....
};

ktnode_suffix::~ktnode_suffix() {
  if (kidlist_) {
    delete kidlist_;
  }
};

ktnode_lsuffix::~ktnode_lsuffix() {
  if (kidlist_) {
    delete kidlist_;
  }
};

ktnode_list * const ktnode_list::lchild(unsigned char ch) const {
  // checkpoint;
  if (!kidlist_) return 0;
  // checkpoint;
  tinylist<kid_list_element<ktnode_list> >::iterator it;
  // checkpoint;
  for (it=kidlist_->begin(); it!=kidlist_->end(); ++it) {
    // checkpoint;
    // cerr << ch << " " << it->ch << " " << (void*)it->nd << endl;
    if (ch == it->ch) {
      return it->nd;
    } else if (ch < it->ch) {
      // checkpoint;
      return 0;
    }
  }
  // checkpoint;
  return 0;
}

ktnode_list * const ktnode_list::child(unsigned char ch) const {
  return this->lchild(ch);
}

ktnode_dna_list * const ktnode_dna_list::lchild(unsigned char ch) const {
  if (!kidlist_) return 0;
  tinylist<kid_list_element<ktnode_dna_list> >::iterator it;
  for (it=kidlist_->begin(); it!=kidlist_->end(); ++it) {
    if (ch == it->ch) {
      return it->nd;
    } else if (ch < it->ch) {
      return 0;
    }
  }
  return 0;
}

ktnode_dna_list * const ktnode_dna_list::child(unsigned char ch) const {
  if (ch < 4) {
    return acgt_[ch];
  }
  if (!kidlist_) return 0;
  tinylist<kid_list_element<ktnode_dna_list> >::iterator it;
  for (it=kidlist_->begin(); it!=kidlist_->end(); ++it) {
    if (ch == it->ch) {
      return it->nd;
    } else if (ch < it->ch) {
      return 0;
    }
  }
  return 0;
}

ktnode_jtable * const ktnode_jtable::lchild(unsigned char ch) const {
  if (!kidptr_) 
    return 0;
  
  tinylist<kid_list_element<ktnode_jtable> >::iterator it;
  for (it=((tinylist<kid_list_element<ktnode_jtable> >*)kidptr_)->begin(); 
       it!=((tinylist<kid_list_element<ktnode_jtable> >*)kidptr_)->end(); 
       ++it) {
    if (ch == it->ch) {
      return it->nd;
    } else if (ch < it->ch) {
      return 0;
    }
  } 
  return 0;
}

ktnode_suffix * const ktnode_suffix::lchild(unsigned char ch) const {
  if (!kidlist_) 
    return 0;
  
  tinylist<kid_list_element<ktnode_suffix> >::iterator it;
  for (it=kidlist_->begin(); 
       it!=kidlist_->end(); 
       ++it) {
    if (ch == it->ch) {
      return it->nd;
    } else if (ch < it->ch) {
      return 0;
    }
  } 
  return 0;
}

ktnode_lsuffix * const ktnode_lsuffix::lchild(unsigned char ch) const {
  if (!kidlist_) 
    return 0;
  
  tinylist<kid_list_element<ktnode_lsuffix> >::iterator it;
  for (it=kidlist_->begin(); 
       it!=kidlist_->end(); 
       ++it) {
    if (ch == it->ch) {
      return it->nd;
    } else if (ch < it->ch) {
      return 0;
    }
  } 
  return 0;
}

ktnode_list * const ktnode_list::add_child(keyword_tree<ktnode_list> & kt, 
					   unsigned char ch) {
  ktnode_list *n = kt.newktnode();
  if (!kidlist_) {
    kidlist_ = new tinylist<kid_list_element<ktnode_list> >;
    kidlist_->push_front(kid_list_element<ktnode_list>(ch,n));
    return n;
  }
  tinylist<kid_list_element<ktnode_list> >::iterator it;
  tinylist<kid_list_element<ktnode_list> >::iterator it1=kidlist_->end();
  for (it=kidlist_->begin(); it!=kidlist_->end(); ++it) {
    if (ch < it->ch) {
      break;
    }
    it1=it;
  }
  if (it1==kidlist_->end()) {
    kidlist_->push_front(kid_list_element<ktnode_list>(ch,n));
  } else {
    kidlist_->insert_after(it1,kid_list_element<ktnode_list>(ch,n));
  }
  return n;
}

ktnode_dna_list * const ktnode_dna_list::add_child(keyword_tree<ktnode_dna_list> & kt, 
						   unsigned char ch) {
  ktnode_dna_list *n = kt.newktnode();
  if (!kidlist_) {
    kidlist_ = new tinylist<kid_list_element<ktnode_dna_list> >;
    kidlist_->push_front(kid_list_element<ktnode_dna_list>(ch,n));
    return n;
  }
  tinylist<kid_list_element<ktnode_dna_list> >::iterator it;
  tinylist<kid_list_element<ktnode_dna_list> >::iterator it1=kidlist_->end();
  for (it=kidlist_->begin(); it!=kidlist_->end(); ++it) {
    if (ch < it->ch) {
      break;
    }
    it1=it;
  }
  if (it1==kidlist_->end()) {
    kidlist_->push_front(kid_list_element<ktnode_dna_list>(ch,n));
  } else {
    kidlist_->insert_after(it1,kid_list_element<ktnode_dna_list>(ch,n));
  }
  return n;
}

ktnode_jtable * const ktnode_jtable::add_child(keyword_tree<ktnode_jtable> & kt, 
					       unsigned char ch) {
  ktnode_jtable *n = kt.newktnode();
  if (!kidptr_) {
    kidptr_ = (void*)(new tinylist<kid_list_element<ktnode_jtable> >);
    ((tinylist<kid_list_element<ktnode_jtable> >*)kidptr_)->push_front(kid_list_element<ktnode_jtable>(ch,n));
    return n;
  }
  tinylist<kid_list_element<ktnode_jtable> >::iterator it;
  tinylist<kid_list_element<ktnode_jtable> >::iterator it1;
  it1=((tinylist<kid_list_element<ktnode_jtable> >*)kidptr_)->end();
  for (it=((tinylist<kid_list_element<ktnode_jtable> >*)kidptr_)->begin(); 
       it!=((tinylist<kid_list_element<ktnode_jtable> >*)kidptr_)->end(); 
       ++it) {
    if (ch < it->ch) {
      break;
    }
    it1=it;
  }
  if (it1==((tinylist<kid_list_element<ktnode_jtable> >*)kidptr_)->end()) {
    ((tinylist<kid_list_element<ktnode_jtable> >*)kidptr_)->push_front(kid_list_element<ktnode_jtable>(ch,n));
  } else {
    ((tinylist<kid_list_element<ktnode_jtable> >*)kidptr_)->insert_after(it1,kid_list_element<ktnode_jtable>(ch,n));
  }
  return n;
}

ktnode_suffix*const ktnode_suffix::add_child(keyword_tree<ktnode_suffix> & kt, 
					       unsigned char ch) {
  ktnode_suffix *n = kt.newktnode();
  n->level_ = level_+1;
  if (!kidlist_) {
    kidlist_ = new tinylist<kid_list_element<ktnode_suffix> >;
    kidlist_->push_front(kid_list_element<ktnode_suffix>(ch,n));
    return n;
  }
  tinylist<kid_list_element<ktnode_suffix> >::iterator it;
  tinylist<kid_list_element<ktnode_suffix> >::iterator it1;
  it1=kidlist_->end();
  for (it=kidlist_->begin(); 
       it!=kidlist_->end(); 
       ++it) {
    if (ch < it->ch) {
      break;
    }
    it1=it;
  }
  if (it1==kidlist_->end()) {
    kidlist_->push_front(kid_list_element<ktnode_suffix>(ch,n));
  } else {
    kidlist_->insert_after(it1,kid_list_element<ktnode_suffix>(ch,n));
  }
  return n;
}

ktnode_lsuffix*const ktnode_lsuffix::add_child(keyword_tree<ktnode_lsuffix> & kt, 
					       unsigned char ch) {
  ktnode_lsuffix *n = kt.newktnode();
  n->level_ = level_+1;
  if (!kidlist_) {
    kidlist_ = new tinylist<kid_list_element<ktnode_lsuffix> >;
    kidlist_->push_front(kid_list_element<ktnode_lsuffix>(ch,n));
    return n;
  }
  tinylist<kid_list_element<ktnode_lsuffix> >::iterator it;
  tinylist<kid_list_element<ktnode_lsuffix> >::iterator it1;
  it1=kidlist_->end();
  for (it=kidlist_->begin(); 
       it!=kidlist_->end(); 
       ++it) {
    if (ch < it->ch) {
      break;
    }
    it1=it;
  }
  if (it1==kidlist_->end()) {
    kidlist_->push_front(kid_list_element<ktnode_lsuffix>(ch,n));
  } else {
    kidlist_->insert_after(it1,kid_list_element<ktnode_lsuffix>(ch,n));
  }
  return n;
}

// std::list<kid_list_element> const & ktnode::children() const {
//   return *kids_;
// }
/* 

void ktnode::writeref(ostream & os) const {
  os << "{";
  os << tag_;
  os << "}"; 
}

void ktnode::write(ostream & os, int indent, bool recurse) const {
  char brack1='[', brack2=']';
  if (patid_) {
    brack1 = '(';
    brack2 = ')';
  }
  os << brack1 
     << std::string("%*d",(int)(floor(log10(kt_.nodetag_)+1)),tag_) 
     << brack2
     << " ";
  indent += 3;
  indent += (int)(floor(log10(kt_.nodetag_)+1));
  // os << "(" << tag_ << " " << tmp << ")";
  list_item li;
  forall_items(li,kids_) {
    os << kids_[li].first() << " => ";
    if (recurse) {
      kids_[li].second()->write(os,indent+5);
    } else {
      kids_[li].second()->writeref(os);
    }
    // if (li != kids_.last()) {
    os << endl;
    for(int i=0;i<indent;i++) {
      os << " ";
    }
    // }
  }
  os << "! => ";
  if (fail()) {
    fail()->writeref(os);
  } else {
    os << "NULL";
  }
  os << endl;
  for(int i=0;i<indent;i++) {
    os << " ";
  }
  os << "> => ";
  if (output()) {
    output()->writeref(os);
  } else {
    os << "NULL";
  }
  
}

*/


/* 

void keyword_tree::write(ostream & os) const {
  pattern_list_element ple;
  forall(ple,patterns()) {
    ple.write(os);
    os << endl;
  }
  if (root_) {
    root_->write(os,0);
  }
}

ostream & operator<<(ostream & os, keyword_tree const & kt) {
  kt.write(os);
  return os;
}

*/

void ktnode_list::optimize_node(keyword_tree<ktnode_list> & kt,
				   CharacterProducer & cp) {
  //noop
}

void ktnode_lsuffix::optimize_node(keyword_tree<ktnode_lsuffix> & kt,
				   CharacterProducer & cp) {
  //noop
}

void ktnode_dna_list::optimize_node(keyword_tree<ktnode_dna_list> & kt,
				    CharacterProducer & cp) {
  if (kidlist_) {
    tinylist<kid_list_element<ktnode_dna_list> >::iterator it;
    tinylist<kid_list_element<ktnode_dna_list> >::iterator it1;
    it=kidlist_->begin();
    it1=kidlist_->end();
    while (it!=kidlist_->end()) {
      if (it->nd) {
	it->nd->optimize_node(kt,cp);
      }
      unsigned char nch = it->ch;
      if (/*nch >= 0 &&*/ nch < 4) {
	acgt_[nch] = it->nd;
	if (it1 != kidlist_->end()) {
	  it = kidlist_->erase_after(it1);
	} else {
	  it = kidlist_->pop();
	}
      } else {
	// it->ch = nch;
	it1 = it;
	++it;
      }
    }
    if (kidlist_->empty()) {
      delete kidlist_;
      kidlist_=0; 
    }
  }
}

void ktnode_jtable::optimize_node(keyword_tree<ktnode_jtable> & kt, 
				  CharacterProducer & cp) {
  kid_jump_table<ktnode_jtable> *compatible_kjt=0;
  tinylist<kid_list_element<ktnode_jtable> > *kids=
    ((tinylist<kid_list_element<ktnode_jtable> > *)kidptr_);
  tinylist<kid_list_element<ktnode_jtable> >::iterator kit;
  std::list<kid_jump_table<ktnode_jtable> *>::iterator kjtit;
  for (kjtit=kt.jump_tables().begin();kjtit!=kt.jump_tables().end();++kjtit) {
    kid_jump_table<ktnode_jtable> *kjt = *kjtit;
    if (kids) {
      bool conflict=false;
      for (kit=kids->begin();kit!=kids->end();++kit) {
	int nch = kit->ch;
	if (nch >= 0 && kjt[nch].which) {
	  conflict=true;
	  break;
	}
      }
      if (!conflict) {
	compatible_kjt = kjt;
	break;
      }
    }
  }
  if (!compatible_kjt) {
    compatible_kjt = kt.newktjumptable(cp.size());
  }
  kidptr_ = (void*)compatible_kjt;
  if (kids) {
    for (kit=kids->begin();kit!=kids->end();++kit) {
      int nch = kit->ch;
      if (nch >= 0) {
	((kid_jump_table<ktnode_jtable>*)kidptr_)[nch].child = kit->nd;
	((kid_jump_table<ktnode_jtable>*)kidptr_)[nch].which = this;
      }
    }
    for (kit=kids->begin();kit!=kids->end();++kit) {
      kit->nd->optimize_node(kt,cp);
    }
  }
  delete kids;
}

void ktnode_suffix::optimize_node(keyword_tree<ktnode_suffix> & kt,
				  CharacterProducer & cp) {
  // checkpoint;
  // cerr << level_ << endl;
  kid_jump_table<ktnode_suffix> *compatible_kjt=0;
  tinylist<kid_list_element<ktnode_suffix> > *kids=kidlist_;
  tinylist<kid_list_element<ktnode_suffix> >::iterator kit;
  std::list<kid_jump_table<ktnode_suffix> *>::iterator kjtit;
  for (kjtit=kt.jump_tables().begin();kjtit!=kt.jump_tables().end();++kjtit) {
    kid_jump_table<ktnode_suffix> *kjt = *kjtit;
    if (kids) {
      bool conflict=false;
      for (kit=kids->begin();kit!=kids->end();++kit) {
	int nch = kit->ch;
	if (nch >= 0 && kjt[nch].which) {
	  conflict=true;
	  break;
	}
      }
      if (!conflict) {
	compatible_kjt = kjt;
	break;
      }
    }
  }
  if (!compatible_kjt) {
    compatible_kjt = kt.newktjumptable(cp.size());
  }
  kidjtable_ = compatible_kjt;
  if (kids) {
    for (kit=kids->begin();kit!=kids->end();++kit) {
      int nch = kit->ch;
      if (nch >= 0) {
	kidjtable_[nch].child = kit->nd;
	kidjtable_[nch].which = this;
      }
    }
    for (kit=kids->begin();kit!=kids->end();++kit) {
      kit->nd->optimize_node(kt,cp);
    }
  }
  // checkpoint;
}

