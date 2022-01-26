
#ifndef _IBPEP_keyword_tree_h
#define _IBPEP_keyword_tree_h

#include <iostream>
#include <list>
#include "util.h"
#include "char_io.h"
#include "pattern_alignment.h"
#include "pattern_match.h"
#include "memory_debug.h"
#include "tinylist.t"

class ktnode_list;
class ktnode_dna_list;
class ktnode_jtable;

template<class KTND>
class kid_list_element {
 public:
  unsigned char ch;
  KTND *nd;
  kid_list_element(unsigned char c=255,KTND* n=0) : ch(c),nd(n) {};
  MEMORY_DEBUG(kid_list_element<KTND>)
};

template <class KTND>
struct kid_jump_table {
  KTND *which;
  KTND *child;
  MEMORY_DEBUG(kid_jump_table<KTND>)
};

template <class KTND>
class keyword_tree : public PatternMatch {
 private:
  bool first_char_;
  unsigned char ch_;
  KTND * w_;
  std::list<kid_jump_table<KTND>*> jump_tables_;
  std::list<kid_jump_table<KTND>*> jump_table_delete_list;
  KTND *root_;
  std::list<KTND*> ktnodes;
  std::list<KTND*> ktnode_delete_list;
  bool *relchars_;
  void add_keyword_(KTND * v, char const * const pat,  unsigned int len,
		    pattern_list::const_iterator const & it);
  tinylist<pattern_list::const_iterator> * 
    lookup_keyword_(KTND * v, char const * const pat, int len);
  void del_keyword_(KTND * v, char const * const pat, int len, long unsigned int);
  bool has_novel_lchild(KTND *u, KTND *v) const;
  void failure_links_(std::list<KTND*> & q);
  void compute_failure_links();
  void optimize_failure_links_(std::list<KTND*> & q);
  void optimize_failure_links();
  void optimize_nodes(CharacterProducer & cp);
  void compress_relchars(CharacterProducer & cp);
 public:
  keyword_tree();
  ~keyword_tree();
  unsigned long int add_pattern(std::string const & pat, 
				unsigned long int id=0,
				int esb=0, int eeb=0);
  void del_pattern(CharacterProducer & cp, std::string const & pat, 
		   long unsigned int=0);
  void lookup_pattern(CharacterProducer & cp, std::string const & pat,
		      std::list<long unsigned int> & l);
  void init(CharacterProducer & cp);
  bool find_patterns(CharacterProducer & cp, 
		     pattern_hit_vector & pas,
		     long unsigned minka=1);
  bool find_suffixes(CharacterProducer & cp, 
		     pattern_hit_vector & pas,
		     long unsigned int minlevel=1);
  void reset();
  KTND *newktnode();
  std::list<kid_jump_table<KTND>*> & jump_tables() {
    return jump_tables_;
  }
  kid_jump_table<KTND> *newktjumptable(unsigned int size);
  MEMORY_DEBUG(keyword_tree<KTND>)
};

template<class KTND>
class ktnode {
private:
  KTND *fail_;
  KTND *output_;
  tinylist<pattern_list::const_iterator> *patid_;
public:
  ktnode();
  ~ktnode();
  inline KTND * fail() const;
  void fail(KTND * const v);
  inline KTND * output() const;
  void output(KTND * const v);
  tinylist<pattern_list::const_iterator> * const & patid() const;
  void add_patid(pattern_list::const_iterator const & it);
  void del_patid(pattern_list::const_iterator const & it);
  void print_patid();
  void clear_patid();
  MEMORY_DEBUG(ktnode<KTND>)
};

class ktnode_list : public ktnode<ktnode_list> {
 private:
  tinylist<kid_list_element<ktnode_list> > *kidlist_;
 public:
  ktnode_list();
  ~ktnode_list();
  ktnode_list * const lchild(unsigned char ch) const;
  ktnode_list * const child(unsigned char ch) const;
  ktnode_list * const add_child(keyword_tree<ktnode_list> & kt, 
				unsigned char ch);
  tinylist<kid_list_element<ktnode_list> > * kids() const {
    return kidlist_;
  }
  void optimize_node(keyword_tree<ktnode_list> & kt, 
		     CharacterProducer & cp);
  MEMORY_DEBUG(ktnode_list)
};

class ktnode_dna_list : public ktnode<ktnode_dna_list> {
 private:
  tinylist<kid_list_element<ktnode_dna_list> > *kidlist_;
  ktnode_dna_list *acgt_[4];
 public:
  ktnode_dna_list();
  ~ktnode_dna_list();
  ktnode_dna_list * const lchild(unsigned char ch) const;
  ktnode_dna_list * const child(unsigned char ch) const;
  ktnode_dna_list * const add_child(keyword_tree<ktnode_dna_list> & kt, 
				    unsigned char ch);
  tinylist<kid_list_element<ktnode_dna_list> > * kids() const {
    return kidlist_;
  }
  void optimize_node(keyword_tree<ktnode_dna_list> & kt, 
		     CharacterProducer & cp);
  MEMORY_DEBUG(ktnode_dna_list)
};

class ktnode_jtable : public ktnode<ktnode_jtable> {
 private:
  void* kidptr_;
 public:
  ktnode_jtable();
  ~ktnode_jtable();
  ktnode_jtable * const lchild(unsigned char ch) const;
  inline ktnode_jtable * const child(unsigned char ch) const {
    if (((kid_jump_table<ktnode_jtable>*)kidptr_)[ch].which == this) {
      return ((kid_jump_table<ktnode_jtable>*)kidptr_)[ch].child;
    } else {
      return 0;
    }
  } 
  ktnode_jtable * const add_child(keyword_tree<ktnode_jtable> & kt, 
				  unsigned char ch);
  tinylist<kid_list_element<ktnode_jtable> > * kids() const {
    return ((tinylist<kid_list_element<ktnode_jtable> > *)kidptr_);
  }
  void optimize_node(keyword_tree<ktnode_jtable> & kt, 
		     CharacterProducer & cp);
  MEMORY_DEBUG(ktnode_jtable)
};

class ktnode_suffix : public ktnode<ktnode_suffix> {
 private:
  tinylist<kid_list_element<ktnode_suffix> > * kidlist_;
  kid_jump_table<ktnode_suffix>* kidjtable_;
  unsigned int level_;
 public:
  ktnode_suffix();
  ~ktnode_suffix();
  ktnode_suffix * const lchild(unsigned char ch) const;
  inline ktnode_suffix * const child(unsigned char ch) const {
    if (kidjtable_[ch].which == this) {
      return kidjtable_[ch].child;
    } else {
      return 0;
    }
  } 
  ktnode_suffix * const add_child(keyword_tree<ktnode_suffix> & kt, 
				  unsigned char ch);
  tinylist<kid_list_element<ktnode_suffix> > * kids() const {
    return kidlist_;
  }
  void optimize_node(keyword_tree<ktnode_suffix> & kt, 
		     CharacterProducer & cp);
  unsigned int level() const {
    return level_;
  }
  MEMORY_DEBUG(ktnode_suffix)
};

class ktnode_lsuffix : public ktnode<ktnode_lsuffix> {
 private:
  tinylist<kid_list_element<ktnode_lsuffix> > * kidlist_;
  unsigned int level_;
 public:
  ktnode_lsuffix();
  ~ktnode_lsuffix();
  ktnode_lsuffix * const lchild(unsigned char ch) const;
  inline ktnode_lsuffix * const child(unsigned char ch) const {
    return lchild(ch);
  } 
  ktnode_lsuffix * const add_child(keyword_tree<ktnode_lsuffix> & kt, 
				  unsigned char ch);
  tinylist<kid_list_element<ktnode_lsuffix> > * kids() const {
    return kidlist_;
  }
  void optimize_node(keyword_tree<ktnode_lsuffix> & kt, 
		     CharacterProducer & cp);
  unsigned int level() const {
    return level_;
  }
  MEMORY_DEBUG(ktnode_lsuffix)
};

#endif
