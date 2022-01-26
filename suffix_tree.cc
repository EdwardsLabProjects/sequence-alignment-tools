
#include "suffix_tree.h"
#include <string.h>
#include <assert.h>
#include "char_io.h"
#include "util.h"
#include <unistd.h>

stnode::stnode() : kidptr_(0), pos_(0), suffix_link_(0) {
#ifdef ASSUME_DNA_CHARS
  acgt_[0] = 0;  acgt_[1] = 0;  acgt_[2] = 0;  acgt_[3] = 0;
#endif /*ASSUME_DNA_CHARS*/
  tag_ = 0;
  level_ = 0;
}

stnode::~stnode() {
}

st_kid_list_element & const stnode::lchild(unsigned char ch) const {
#ifdef ASSUME_DNA_CHARS
  if (ch < 4) {
    return acgt_[ch];
  }
#endif
  slist_item li;
  if (!kidptr_) return 0;
  li = ((std::list<kid_list_element> *)kidptr_)->first();
  while (li) {
    kid_list_element const & kle = ((std::list<kid_list_element> *)kidptr_)->inf(li);
    if (ch == kle.ch) {
      return kle.nd;
    } else if (ch > kle.ch) {
      return 0;
    }
    li = ((std::list<kid_list_element> *)kidptr_)->succ(li);
  }
  return 0;
}

stnode * const stnode::add_child(suffix_tree & st, unsigned char ch, 
				 FILE_POSITION_TYPE s, FILE_POSITION_TYPE e) {
  stnode *n = st.newstnode();
  n->level_ = level_+1;
  st_kid_list_element stkle(ch,s,e,n);
  if (!kidptr_) {
    kidptr_ = (void*)(new std::list<st_kid_list_element>);
  }
  ((std::list<st_kid_list_element> *)kidptr_)->append(stkle);
  return n;
}

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

suffix_tree::suffix_tree() { 
  nodetag_ = 0;
  root_ = newstnode();
};

suffix_tree::~suffix_tree() {
  // cerr << __FILE__ << ":" << __LINE__ << endl;  
#ifdef USE_JUMP_TABLES
  st_kid_jump_table *kjt;
  forall(kjt,jump_table_delete_list) {
    delete [] kjt;
  }
#endif /*USE_JUMP_TABLES*/
  stnode *stn;
  forall(stn,stnode_delete_list) {
    delete [] stn;
  }
}

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

stnode *suffix_tree::newstnode() {
  if (stnodes.empty()) {
    int size=getpagesize()/sizeof(stnode);
    stnode *buffer = new stnode[size];
    for (int i=0;i<size;i++) {
      buffer[i].tag_ = node_tag_++;
      stnodes.append(buffer+i);
    }
    stnode_delete_list.append(buffer);
  }
  return stnodes.pop();
}

#ifdef ASSUME_DNA_CHARS
void suffix_tree::pack_dna_chars_array(CharacterProducer & cp) {
  root_->pack_dna_chars_array(cp);
}

void stnode::pack_dna_chars_array() {
  slist_item li,li0=0,li1=0;
  std::list<st_kid_list_element> *kids = 
    (std::list<st_kid_list_element> *)kidptr_;
  if (kids) {
    li = kids->first();
    while (li) {
      if (kids->inf(li).nd) {
	kids->inf(li).nd->pack_dna_chars_array(cp);
      }
      int nch = kids->inf(li).ch;
      if (nch >= 0 && nch < 4) {
	acgt_[nch] = kids->inf(li).nd;
      } 
      li0 = li;
      li = kids->succ(li);
      if (nch >= 0 && nch < 4) {
	if (li1) {
	  kids->del_succ_item(li1);
	  li1 = li0;
	} else {
	  kids->pop();
	  li1 = 0;
	}
      }
    }
    if (kids->empty()) {
      delete kids;
      kidptr_=0;
    }
  }
}
#endif

#ifdef USE_JUMP_TABLES

void suffix_tree::optimize_jump_tables(CharacterProducer & cp) {
  root_->optimize_jump_table(*this,cp);
}

void stnode::optimize_jump_table(suffix_tree & st, CharacterProducer & cp) {
  st_kid_jump_table *compatible_kjt=0;
  std::list<st_kid_list_element> *kids=
    ((std::list<st_kid_list_element> *)kidptr_);
  kid_jump_table *kjt;
  forall(kjt,kt.jump_tables) {
    slist_item li=0;
    if (kids) {
      li = kids->first();
    }
    while (li) {
      int nch = (*kids)[li].ch;
      if (nch >= 0 && kjt[nch].which) {
	break;
      }
      li = kids->succ(li);
    }
    if (!li) {
      compatible_kjt = kjt;
      break;
    }
  }
  if (!compatible_kjt) {
    int size=cp.size();
    int number=getpagesize()/(size*sizeof(kid_jump_table));
    kid_jump_table *buffer= new kid_jump_table[size*number];
    for (int i=0;i<size*number;i++) {
      buffer[i].which=0;
    }
    for (int i=0;i<number;i++) {
      kt.jump_tables.append(buffer+i*size);
    }
    kt.jump_table_delete_list.append(buffer);
    compatible_kjt = buffer;
  }
  kidptr_ = (void*)compatible_kjt;
  if (kids) {
    slist_item li;
    li = kids->first();
    while (li) {
      int nch = cp.nch(kids->inf(li).ch);
      if (nch >= 0) {
	((kid_jump_table*)kidptr_)[nch].child = kids->inf(li).nd;
	((kid_jump_table*)kidptr_)[nch].which = this;
      }
      li = kids->succ(li);
    }
    li = kids->first();
    while (li) {
      kids->inf(li).nd->optimize_jump_table(kt,cp);
      li = kids->succ(li);
    }
  }
  delete kids;
}

#endif /*USE_JUMP_TABLES*/

void suffix_tree::build_tree() {
  while (!cp.eof) {
    
  }
}
