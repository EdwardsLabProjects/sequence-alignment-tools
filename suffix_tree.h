
#ifndef _IBPEP_suffix_tree_h
#define _IBPEP_suffix_tree_h

#include <iostream>
#include "util.h"
#include "char_io.h"
#include "pattern_alignment.h"
#include "pattern_match.h"
#include "types.h"

#ifdef ASSUME_DNA_CHARS
#undef USE_JUMP_TABLES
#endif

class stnode;
class suffix_tree;

struct st_kid_list_element {
  char ch;
  FILE_POSITION_TYPE st;
  FILE_POSITION_TYPE ed;
  stnode *nd;
  st_kid_list_element(char c=0,FILE_POSITION_TYPE s=0, 
		      FILE_POSITION_TYPE e=0, stnode* n=0) 
    : ch(c),st(s),ed(e),nd(n) {};
};

#ifdef USE_JUMP_TABLES
struct st_kid_jump_table {
  stnode *which;
  stnode *child;
};
#endif /*USE_JUMP_TABLES*/

class suffix_tree  {
 private:
  unsigned char *string_;
  FILE_SIZE_TYPE string_len_;
  unsigned long nodetag_;
#ifdef USE_JUMP_TABLES
  std::list<st_kid_jump_table*> jump_tables;
  std::list<st_kid_jump_table*> jump_table_delete_list;
#endif /*USE_JUMP_TABLES*/
  stnode *root_;
  std::list<stnode*> stnodes;
  std::list<stnode*> stnode_delete_list;
#ifdef USE_JUMP_TABLES
  void optimize_jump_tables(CharacterProducer & cp);
#endif /*USE_JUMP_TABLES*/
  stnode *newstnode();
 public:
  suffix_tree();
  ~suffix_tree();
  void build_tree(Character_Producer & cp);
  void write(ostream & os) const;
  void optimize_tree(CharacterProducer & cp) {
#ifdef USE_JUMP_TABLES
    optimize_jump_tables(cp);
#endif /*USE_JUMP_TABLES*/
  }
  friend class stnode;
};

ostream & operator<<(ostream & os, suffix_tree const & kt);

class stnode {
 private:
  void* kidptr_;
  stnode* suffix_link_;
  FILE_POSITION_TYPE pos_;
  unsigned long tag_;
  unsigned long level_;
#ifdef ASSUME_DNA_CHARS
  stnode *acgt_[4];
#endif
public:
  stnode();
  ~stnode();
  stnode * const lchild(unsigned char ch) const;
#ifdef USE_JUMP_TABLES
  inline ktnode * const child(unsigned char ch) const {
    if (((kid_jump_table*)kidptr_)[ch].which == this) {
      return ((kid_jump_table*)kidptr_)[ch].child;
    } else {
      return 0;
    }
  } 
#else
  ktnode * const child(unsigned char ch) const {
    return lchild(ch);
  }
#endif /*USE_JUMP_TABLES*/
  stnode * const add_child(suffix_tree & st, unsigned char ch, 
			   FILE_POSITION_TYPE s, FILE_POSITION_TYPE e);
#ifdef USE_JUMP_TABLES
  void optimize_jump_table(suffix_tree & st, CharacterProducer & cp);
#endif /*USE_JUMP_TABLES*/
#ifdef ASSUME_DNA_CHARS
  void pack_dna_chars_array(Character_Producer & cp); 
#endif /*ASSUME_DNA_CHARS*/
  inline stnode * const suffix() const {
    return suffix_link_;
  }
  void suffix(stnode * const v) {
    suffix_link_= v;
  }
  inline FILE_POSITION_TYPE pos() const {
    return pos_;
  }
  void pos(FILE_POSITION_TYPE const & p) {
    pos_ = p;
  }
  inline long unsigned const & level() const {
    return level_;
  }
  void write(ostream & os, int indent=0, bool recurse=true) const;
  void writeref(ostream & os) const;
  friend class suffix_tree;
};

#endif
