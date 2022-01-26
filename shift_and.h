
#ifndef _IBPEP_SHIFT_AND_H
#define _IBPEP_SHIFT_AND_H

#include <iostream>
#include "char_io.h"
#include "pattern_match.h"
#include "types.h"

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

class shift_and : public PatternMatch {
  struct patbit {
    unsigned int bit;
    tinylist<pattern_list_element>::const_iterator it;
  };
  bigword *m_;
  bigword **u_;
  bigword *mask_;
  bigword *s_;
  unsigned int k_;
  unsigned char eos_;
  bool _wc;
  bool _textn;
  bool _re;
  unsigned int _wordcount;
  bigword _highbit;
  unsigned int _wordbits_1;
  patbit *_patbits;
  unsigned int *_patbitind;
  void computeu(CharacterProducer & cp);
  void clearu();
public:
  shift_and(bool wc=false, bool tn=false, bool re=false, unsigned char eos='\n');
  ~shift_and();
  long unsigned int add_pattern(std::string const & pat, unsigned long id=0,
				int esb=0, int eeb=0);
  void write(ostream & os) const;
  bool find_patterns(CharacterProducer & cp, 
		     pattern_hit_vector & pas,
		     long unsigned minka=1);
  void init(CharacterProducer & cp) { computeu(cp); }
  void reset();
  bool wildcards() const;
  void wildcards(bool wc);
  bool wildcard_text_N() const;
  void wildcard_text_N(bool tn);
};

#endif
