
#ifndef _IBPEP_SUFTREE_H
#define _IBPEP_SUFTREE_H

#include <iostream>
#include "char_io.h"
#include "pattern_match.h"
#include "types.h"

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

#include "rl_suffix_tree.h"

class RLSufTree: public suffix_tree<basic_node,basic_suffix> {
public:
  RLSufTree(const char *st,unsigned len,char t='\n'):
    suffix_tree<basic_node,basic_suffix>(st,len,t) {}
};

class suftree : public PatternMatch {
  unsigned char eos_;
  tinylist<pattern_list_element>::const_iterator curit;
  RLSufTree *stree;
  unsigned *pos;
  unsigned poslen;

public:
  suftree(unsigned char eos='\n');
  ~suftree();
  long unsigned int add_pattern(std::string const & pat, unsigned long id=0,
				int esb=0, int eeb=0);
  void write(ostream & os) const;
  bool find_patterns(CharacterProducer & cp, 
		     pattern_hit_vector & pas,
		     long unsigned minka=1);
  void init(CharacterProducer & cp);
  void reset();
  static std::string filename(char const * filename) {
    return (std::string(filename)+".st");
  }
};

#endif
