
#ifndef _IBPEP_HASH_TABLE_H
#define _IBPEP_HASH_TABLE_H

#include "types.h"
#include "pattern_match.h"

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

class htele {
private:
  pattern_list::const_iterator patid_;
  unsigned int position_;
public:
  htele(pattern_list::const_iterator const & it,
	unsigned int p=0) : patid_(it), position_(p) {};
  ~htele() {};
  pattern_list::const_iterator const & pattern_list_it() const {
    return patid_;
  }
  void position(unsigned int p) {
    position_ = p;
  }
  unsigned int position() const {
    return position_;
  }
};

class hash_table_elt {
public:
  hash_table_elt();
  ~hash_table_elt();
  tinylist<htele> * const & patids() const;
  void add_patid(pattern_list::const_iterator const & it,
		 unsigned int p=0);
private:
  tinylist<htele> *patid_;
};

class hash_table : public PatternMatch {
  unsigned int ws_;
  unsigned int alphasize_;
  unsigned int alphalog_;
  int p_;
  bigword h_;
  bigword wsmask_;
  bool *relchars_;
  unsigned char *relcharmap_;
  std::vector<hash_table_elt> table_;
  std::vector<FILE_POSITION_TYPE> lastpos_;
  long unsigned int max_id_;
  unsigned int k_;
  unsigned char eos_;
  bool _wc;
  bool _textn;
  bool _indels;
  bool _dna_mut;
  long unsigned int _mpl;
public:
  hash_table(unsigned int ws=0, unsigned int k=0, unsigned char eos='\n', 
	     bool wc=false, bool tn=false, bool indels=true, bool dna_mut=false);
  ~hash_table();
  long unsigned int add_pattern(std::string const & pat, unsigned long id=0,
				int esb=0, int eeb=0);
  bool find_patterns(CharacterProducer & cp, 
		     pattern_hit_vector & pas,
		     long unsigned minka=1);
  void init(CharacterProducer & cp);
  void reset();
};

#endif
