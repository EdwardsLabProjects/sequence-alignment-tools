
#ifndef _IBPEP_EXACT_BASES_H
#define _IBPEP_EXACT_BASES_H

#include <iostream>
#include <utility>
#include <vector>
#include <map>
#include "char_io.h"
#include "keyword_tree.h"
#include "types.h"

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

class exact_bases : public PatternMatch {
private:
  PatternMatch *pm_;
  long unsigned int num_patterns_;
  unsigned int k_;
  unsigned char eos_;
  bool _wc;
  bool _textn;
  bool _indels;
  bool _dna_mut;
  long unsigned int _mpl;
  std::vector<std::string> rempat_;
  std::vector<bool> prefix_;
  std::vector<pattern_list::const_iterator> plit_;
public:
  exact_bases(PatternMatch *pm,
	      unsigned int k=0, unsigned char eos='\n', 
	      bool wc=false, bool tn=false, bool indels=true, bool dna_mut=false);
  exact_bases(PatternMatch *pm);
  ~exact_bases();
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
  void indels(bool tn);
  bool dna_mut() const;
  void dna_mut(bool dm);
  unsigned char eos_char() const;
  void eos_char(unsigned char c);
};

#endif


