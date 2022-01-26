
#ifndef _IBPEP_FILTER_BITVEC_H
#define _IBPEP_FILTER_BITVEC_H

#include <iostream>
#include <utility>
#include <vector>
#include <map>
#include "char_io.h"
#include "keyword_tree.h"
#include "types.h"
#include "shift_and_inexact.h"

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

class filter_bitvec : public PatternMatch {
private:
  shift_and_inexact *pm_;
  long unsigned int num_patterns_;
  long unsigned int max_patlen_;
  pattern_hit_vector l;
  std::vector<pattern_list::const_iterator> plit_;
 public:
  filter_bitvec(unsigned int k=0, char eos='\n', 
		bool wc=false, bool tn=false, bool indels=true, bool dna_mut=false);
  filter_bitvec();
  ~filter_bitvec();
  long unsigned int add_pattern(std::string const & pat, unsigned long id=0,
				int esb=0, int eeb=0);
  bool find_patterns(CharacterProducer & cp, 
		     pattern_hit_vector & kas,
		     long unsigned minka=1);
  void init(CharacterProducer & cp);
  void reset();
  unsigned int mismatches() const;
  void mismatches(unsigned int k);
  bool wildcards() const;
  void wildcards(bool wc);
  bool wildcard_text_N() const;
  void wildcard_text_N(bool tn);
  bool indels() const;
  void indels(bool id);
  bool dna_mut() const;
  void dna_mut(bool dm);
  char eos_char() const;
  void eos_char(char c);
};

#endif


