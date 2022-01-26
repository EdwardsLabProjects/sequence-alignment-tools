
#ifndef _IBPEP_GS_HASH_TABLE_H
#define _IBPEP_GS_HASH_TABLE_H

#include "types.h"
#include "pattern_match.h"
#include "hash_table.h"

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

class gapped_seed_set {
 private:
  void init(int*,int* =0);
 public:
  bool shift;
  int m;
  int n;
  int l;
  int k;
  bool indels;
  std::vector<std::vector<unsigned int> > patpos;
  std::vector<std::vector<unsigned int> > txtpos;
  int tmax;
  typedef enum {

    m5K1l3S,
    m6K1l4S,
    m7K1l5S,
    m8K1l6S,

    m9K2l5S,
    m10K2l5S,

    m19k2l9S,
    m19k2l10S,

    m19K2l10S,
    m19K2l11S,
    m19K2l12S,
    m19K2l13S,
    m19K2l14S,

    m19K3l9S,
    m19K3l10S,

    m19k3l6S,
    m19k3l7S,

    m20k2l9S,
    m20k2l10S,

    m20K2l12S, 
    
    m20K3l9S,
    m20K3l10S,

    m20k3l6S,
    m20k3l7S,

    m28K2l16S,

    m30K2l16S,

    no_scheme

  } scheme;
  gapped_seed_set() {};
  ~gapped_seed_set() {};
  bool initialize(scheme s);
  std::string str() const;
  static scheme select(int m, int k, bool indels);
};

class gsrhtele {
private:
  pattern_list::const_iterator patid_;
  unsigned int position_;
  unsigned int template_;
  unsigned int index_;
public:
  gsrhtele(pattern_list::const_iterator const & it,
	unsigned int p=0, unsigned int t=0, unsigned int i=0) : patid_(it), position_(p), template_(t), index_(i) {};
  ~gsrhtele() {};
  bool operator==(gsrhtele const & e) {
    return ((patid_ == e.patid_) && (position_ == e.position_) && 
	    (template_ == e.template_) && (index_ == e.index_));
  }
  pattern_list::const_iterator const & pattern_list_it() const {
    return patid_;
  }
  ostream & write(ostream & os) const {
    os << patid_->id() << " " 
       << patid_->pattern() << " "
       << template_ << " " 
       << position_ << " " 
       << index_;
    return os;
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
  void templ(unsigned int i) {
    template_ = i;
  }
  unsigned int templ() const {
    return template_;
  }
};

class gs_hash_table_elt {
public:
  gs_hash_table_elt();
  ~gs_hash_table_elt();
  inline tinylist<gsrhtele> * const & patids() { return patid_; }
  void add_patid(pattern_list::const_iterator const & it,
		 unsigned int p=0, unsigned int t=0, unsigned int i=0);
  void clear_patids();
 private:
  tinylist<gsrhtele> *patid_;
};

class gs_hash_table : public PatternMatch {
  gapped_seed_set gss;
  unsigned int asize_;
  unsigned long int must_mod_;
  bool *relchars_;
  unsigned char *relcharmap_;
  std::vector<int> buffer_;
  int bufp_;
  int nprime;
  std::vector<gs_hash_table_elt> table_;
  std::vector<unsigned long int> prime_;
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
  gs_hash_table(gapped_seed_set::scheme s, unsigned int nprimes = 4, 
		unsigned int k=0, unsigned char eos='\n', 
		bool wc=false, bool tn=false, bool indels=true, bool dna_mut=false);
  ~gs_hash_table();
  gapped_seed_set const & scheme() const {
    return gss;
  }
  long unsigned int add_pattern(std::string const & pat, unsigned long id=0,
				int esb=0, int eeb=0);
  void del_pattern(CharacterProducer & cp, std::string const & pat, long unsigned int =0);
  bool find_patterns(CharacterProducer & cp, 
		     pattern_hit_vector & pas,
		     long unsigned minka=1);
  void init(CharacterProducer & cp);
  void reset();
  unsigned long int hash(vector<unsigned char> const & d, 
			 unsigned long int p);
};

#endif
