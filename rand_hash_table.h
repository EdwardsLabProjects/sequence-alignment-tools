
#ifndef _IBPEP_RAND_HASH_TABLE_H
#define _IBPEP_RAND_HASH_TABLE_H

#include "types.h"
#include "pattern_match.h"
#include "hash_table.h"

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

class rhtele {
private:
  pattern_list::const_iterator patid_;
  unsigned int position_;
  unsigned int index_;
public:
  rhtele(pattern_list::const_iterator const & it,
	unsigned int p=0, unsigned int i=0) : patid_(it), position_(p), index_(i) {};
  ~rhtele() {};
  pattern_list::const_iterator const & pattern_list_it() const {
    return patid_;
  }
  void position(unsigned int p) {
    position_ = p;
  }
  unsigned int position() const {
    return position_;
  }
  void index(unsigned int i) {
    index_ = i;
  }
  unsigned int index() const {
    return index_;
  }
};

class rand_hash_table_elt {
public:
  rand_hash_table_elt();
  ~rand_hash_table_elt();
  tinylist<rhtele> * const & patids() const;
  void add_patid(pattern_list::const_iterator const & it,
		 unsigned int p=0, unsigned int i=0);
private:
  tinylist<rhtele> *patid_;
};

class rand_hash_table : public PatternMatch {
  unsigned int ws_;
  unsigned int alphasize_;
  int p_;
  bool *relchars_;
  unsigned char *relcharmap_;
  int nprime;
  std::vector<rand_hash_table_elt> table_;
  std::vector<bigword> h_;
  std::vector<unsigned long int> prime_;
  std::vector<unsigned long int> basetonmodp_;
  std::vector<unsigned char> lastnch_;
  unsigned int lastnchptr_;
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
  rand_hash_table(unsigned int nprimes = 4, 
		  unsigned int ws=0, unsigned int k=0, unsigned char eos='\n', 
	     bool wc=false, bool tn=false, bool indels=true, bool dna_mut=false);
  ~rand_hash_table();
  long unsigned int add_pattern(std::string const & pat, unsigned long id=0,
				int esb=0, int eeb=0);
  bool find_patterns(CharacterProducer & cp, 
		     pattern_hit_vector & pas,
		     long unsigned minka=1);
  void init(CharacterProducer & cp);
  void reset();
  static void random_primes_lt(unsigned long int m, std::vector<long unsigned int> &p);
};

#endif
