
#ifndef _IBPEP_SHIFT_AND_INEXACT_H
#define _IBPEP_SHIFT_AND_INEXACT_H

#include <iostream>
#include "char_io.h"
#include "keyword_tree.h"
#include "types.h"

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

class shift_and_inexact : public PatternMatch {
  struct patbit {
    unsigned int bit;
    pattern_list::const_iterator it;
  };
  bigword **m_;
  bigword **u_;
  bigword *mask_;
  bigword *s_;
  unsigned int k_;
  char eos_;
  bool _wc;
  bool _textn;
  bool _indels;
  bool _dna_mut;
  std::vector<std::pair<int,int> > _eb;
  unsigned int _wordcount;
  bigword _highbit;
  patbit *_patbits;
  unsigned long int *_patbitind;
  void computeu(CharacterProducer & cp);
  void clearu();
public:
  shift_and_inexact(unsigned int k=0, unsigned char eos='\n', 
		    bool wc=false, bool tn=false, bool id=true, bool dna_mut=false);
  ~shift_and_inexact();
  long unsigned int add_pattern(std::string const & pat, unsigned long id=0,
				int esb=0, int eeb=0);
  void write(ostream & os) const;
  bool find_patterns(CharacterProducer & cp, 
		     pattern_hit_vector & kas,
		     long unsigned minka=1);
  void init(CharacterProducer & cp) { computeu(cp); }
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
