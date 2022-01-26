
#ifndef _IBPEP_EXACT_HALVES_H
#define _IBPEP_EXACT_HALVES_H

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

class exact_halves : public PatternMatch {
private:
  PatternMatch *pm_;
  unsigned long int num_patterns_;
  unsigned int k_;
  unsigned char eos_;
  bool _wc;
  bool _textn;
  bool _indels;
  bool _dna_mut;
  long unsigned int _mpl;
  std::vector<FILE_POSITION_TYPE> lasthit_;
  std::vector<tinylist<pattern_list_element>::const_iterator> plit_;
  std::vector<std::string> pattern_halves_;
 private:
  static bool hit_lessthan(pattern_hit_vector::element const &, 
			   pattern_hit_vector::element const &);
 public:
  exact_halves(PatternMatch *pm,
	       unsigned int k=0, unsigned char eos='\n', 
	       bool wc=false, bool tn=false, bool indels=true, bool dna_mut=false);
  exact_halves(PatternMatch *pm);
  ~exact_halves();
  long unsigned int add_pattern(std::string const & pat, unsigned long id=0,
				int esb=0, int eeb=0);
  void del_pattern(CharacterProducer & cp, std::string const & pat,
		   unsigned long int id=0);
  void lookup_pattern(CharacterProducer & cp, std::string const & pat,
		      std::list<long unsigned int> & l);
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
  unsigned char eos_char() const;
  void eos_char(unsigned char c);
};

#endif


