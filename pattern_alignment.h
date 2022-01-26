
#ifndef _IBPEP_pattern_alignment_h
#define _IBPEP_pattern_alignment_h

#include <assert.h>
#include <iostream>
#include <vector>
#include "char_io.h"
#include "alignment_code.h"
#include "memory_debug.h"

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

class pattern_hit {
 private:
  unsigned long pattern_;
  FILE_POSITION_TYPE pos_;
 public:
  pattern_hit(unsigned long k=0, FILE_POSITION_TYPE p=0) 
    : pattern_(k), pos_(p) {}
  virtual ~pattern_hit() {};  
  unsigned long const & pattern_id() const {
    return pattern_;
  }
  void pattern_iter(unsigned long const & id) {
    pattern_ = id;
  }
  FILE_POSITION_TYPE const & pos() const {
    return pos_;
  }
  void pos(FILE_POSITION_TYPE const & e) {
    pos_ = e;
  }
  MEMORY_DEBUG(pattern_hit);
};

class pattern_hit_w_dist : public pattern_hit {
  int dist_;
 public:
  pattern_hit_w_dist(unsigned long k=0, FILE_POSITION_TYPE p=0, int d=-1) 
    : pattern_hit(k,p), dist_(d) {};
  pattern_hit_w_dist() {};
  void editdist(int const & d) {
    dist_ = d;
  }
  int const & editdist() const {
    return dist_;
  }
  MEMORY_DEBUG(pattern_hit_w_dist);
};

class pattern_alignment {  
public:
  pattern_alignment(FILE_POSITION_TYPE e=0, bool yn=false)
    : end_(e), alignment_done_(false), yesno_(yn) { 
    if (!yesno_) {
      stats_.resize(alignment_codes); 
      fill(stats_.begin(),stats_.end(),0); 
    }
  }
  pattern_alignment(pattern_hit const & ph, bool yn=false) 
    : end_(ph.pos()), alignment_done_(false) {
    if (!yesno_) {
      stats_.resize(alignment_codes); 
      fill(stats_.begin(),stats_.end(),0);
    }
  }
  virtual ~pattern_alignment();
  FILE_POSITION_TYPE const & end() const {
    return end_;
  }
  void end(FILE_POSITION_TYPE const & e) {
    end_ = e;
  }
  virtual bool align(CharacterProducer &, std::string const & pattern) =0;
  void reset() {
    alignment_done_ = false;
    if (!yesno_) {
      fill(stats_.begin(),stats_.end(),0);
    }
  }
  int size() const {
    assert(alignment_done_);
    return alignment_.size();
  };
  alignment_code operator[](int i) const {
    assert(alignment_done_);
    return alignment_[i];
  }
  FILE_POSITION_TYPE start() const {
    assert(alignment_done_);
    return start_;
  }
  FILE_POSITION_TYPE length() const {
    assert(alignment_done_);
    return end_-start_+1;
  }
  std::string const & matching_text() const {
    assert(alignment_done_);
    return matching_text_;
  }
  unsigned int stats(alignment_code ac) const {
    assert(alignment_done_);
    return stats_[ac];
  }
  unsigned int editdist() const {
    assert(alignment_done_);
    if (stats_[alignment_constraint_violation] > 0) {
      // checkpoint;
      return MAXINT;
    } else {
      return stats_[alignment_substitution] + 
	stats_[alignment_substitution_1] + 
	2*stats_[alignment_substitution_2] + 
	3*stats_[alignment_substitution_3] + 
	stats_[alignment_insertion] +
	3*stats_[alignment_insertion_3] +
	stats_[alignment_deletion] +
	3*stats_[alignment_deletion_3];
    }
  }
  std::string alignment_string() const {
    assert(alignment_done_);
    std::string r("");
    for (int i=0; i<size(); i++) {
      char ch;
      switch ((*this)[i]) {
      case alignment_equal :
	ch = '|';
	break;
      case alignment_wildcard_equal : 
	ch = '+';
	break;
      case alignment_substitution :
	ch = '*';
	break;
      case alignment_substitution_1 :
	ch = '.';
	break;
      case alignment_substitution_2 :
	ch = ':';
	break;
      case alignment_substitution_3 :
	ch = 'x';
	break;
      case alignment_insertion :
      case alignment_insertion_3 :
	ch = '^';
	break;
      case alignment_deletion :
      case alignment_deletion_3 :
	ch = 'v';
	break;
      case alignment_constraint_violation :
	ch = '!';
	break;
      default:
	ch = ' ';
	break;
      }
      r += ch;
    }
    return r;
  }
  std::string alignment_text() const {
    assert(alignment_done_);
    std::string r("");
    std::string const & mt = matching_text();
    int p=0;
    for (int i=0; i<size(); i++) {
      if ((*this)[i] != alignment_deletion &&
	  (*this)[i] != alignment_deletion_3) {
	r += mt[p];
	p++;
      } else {
	r += "-";
      }
    }
    return r;
  }
  std::string alignment_pattern(std::string const & pat) const {
    assert(alignment_done_);
    std::string r("");
    int p=0;
    for (int i=0; i<size(); i++) {
      if ((*this)[i] != alignment_insertion &&
	  (*this)[i] != alignment_insertion_3) {
	r += pat[p];
	p++;
      } else {
	r += "-";
      }
    }
    return r;
  }
  void yesno(bool yn) {
    yesno_ = yn;
  }
  void value(int v) {
    val_ = v;
  }
  int value() const {
    return val_;
  }
  // void write(ostream & os, 
  // FILE_POSITION_TYPE const seq_pos,
  // std::string const & pattern, 
  // long unsigned int id, bool revcomp);
  virtual void write(ostream & os) const;
  virtual void read(istream & is) {
    is >> pattern_ >> end_;
  }
  MEMORY_DEBUG(pattern_alignment);
private:
  unsigned long pattern_;
  FILE_POSITION_TYPE end_;
  int val_;

protected:
  std::vector<alignment_code> alignment_;
  std::vector<int> stats_;
  std::string matching_text_;
  FILE_POSITION_TYPE start_;
  bool alignment_done_;
  bool yesno_;
};

istream & operator>>(istream & is, pattern_alignment & ka);
ostream & operator<<(ostream & os, pattern_alignment const & ka);

class exact_alignment : public pattern_alignment {
public:
  exact_alignment(FILE_POSITION_TYPE e=0)
    : pattern_alignment(e) {};
  exact_alignment(pattern_hit const & ph)
    : pattern_alignment(ph) {};
  ~exact_alignment() {};
  bool align(CharacterProducer &, std::string const & pattern);
  MEMORY_DEBUG(exact_alignment)
};

class exact_peptide_alignment : public pattern_alignment {
 private:
  std::string lcontext_;
  std::string rcontext_;
  int context_;
 public:
  exact_peptide_alignment(FILE_POSITION_TYPE e=0)
    : pattern_alignment(e), context_(1) {};
  exact_peptide_alignment(pattern_hit const & ph)
    : pattern_alignment(ph), context_(1) {};
  ~exact_peptide_alignment() {};
  bool align(CharacterProducer &, std::string const & pattern);
  std::string lcontext() const {
    assert(alignment_done_);
    return lcontext_;
  }
  std::string rcontext() const {
    assert(alignment_done_);
    return rcontext_;
  }
  void set_context(int c) {
    context_ = c;
  }
  MEMORY_DEBUG(exact_peptide_alignment)
};

class exact_wc_alignment : public pattern_alignment {
 private:
  bool textn_;
public:
  exact_wc_alignment(FILE_POSITION_TYPE e=0, bool tn=false)
    : pattern_alignment(e), textn_(tn) {};
  exact_wc_alignment(pattern_hit const & ph, bool tn=false)
    : pattern_alignment(ph), textn_(tn) {};
  ~exact_wc_alignment() {};
  bool align(CharacterProducer &, std::string const & pattern);
  MEMORY_DEBUG(exact_wc_alignment)
};

class mismatch_alignment : public pattern_alignment {
public:
  mismatch_alignment(FILE_POSITION_TYPE e=0)
    : pattern_alignment(e) {};
  mismatch_alignment(pattern_hit const & ph)
    : pattern_alignment(ph) {};
  ~mismatch_alignment() {};
  bool align(CharacterProducer &, std::string const & pattern);
  MEMORY_DEBUG(mismatch_alignment)
};

class editdist_alignment : public pattern_alignment {
  FILE_POSITION_TYPE end2_;
  unsigned int k_;
  char eos_;
  bool wc_;
  bool textn_;
  bool indels_;
  bool dna_mut_;
  bool trans_;
  int lconst_;
  int rconst_;
  char *buffer_;
  char *buffer1_;
  FILE_POSITION_TYPE bufstart_;
  FILE_POSITION_TYPE bufend_;
  long unsigned int bufsize_;
  unsigned int maxpatlen_;
  unsigned int matsize_;
  unsigned int *dp_;
  int *best_;
public:
  editdist_alignment(FILE_POSITION_TYPE e=0, 
		     FILE_POSITION_TYPE e2=0,
		     unsigned int k=0, char eos='\n', 
		     bool wc=false, bool tn=false, bool id=true, bool dm=false,
		     int lconst=0, int rconst=0, bool yn=false, bool trans=false)
    : pattern_alignment(e,yn), end2_(e2), k_(k), eos_(eos), 
    wc_(wc), textn_(tn), indels_(id), dna_mut_(dm), trans_(trans), lconst_(lconst), rconst_(rconst),
    buffer_(0), buffer1_(0), bufstart_(0), bufend_(0), maxpatlen_(0), matsize_(0), dp_(0), best_(0) {};
  editdist_alignment(FILE_POSITION_TYPE e=0, 
		     unsigned int k=0, char eos='\n', 
		     bool wc=false, bool tn=false, bool id=true, bool dm=false,
		     int lconst=0, int rconst=0, bool yn=false, bool trans=false)
    : pattern_alignment(e,yn), end2_(e), k_(k), eos_(eos), wc_(wc), textn_(tn),
      indels_(id), dna_mut_(dm), trans_(trans), lconst_(lconst), rconst_(rconst),
      buffer_(0), buffer1_(0), bufstart_(0), bufend_(0), maxpatlen_(0), matsize_(0), dp_(0), best_(0) {};
  editdist_alignment(pattern_hit const & ph, unsigned int k=0, char eos='\n', 
		     bool wc=false, bool tn=false, bool id=true, bool dm=false,
		     int lconst=0, int rconst=0, bool yn=false, bool trans=false)
    : pattern_alignment(ph,yn), end2_(ph.pos()), k_(k), eos_(eos), 
    wc_(wc), textn_(tn), indels_(id), dna_mut_(dm), trans_(trans), lconst_(lconst), rconst_(rconst),
    buffer_(0), buffer1_(0), bufstart_(0), bufend_(0), maxpatlen_(0), matsize_(0), dp_(0), best_(0) {};
  ~editdist_alignment() {
    delete [] buffer_;
    delete [] buffer1_;
    delete [] dp_;
    delete [] best_;
  };
  void poslb(FILE_POSITION_TYPE p) {
    pattern_alignment::end(p);
  }
  void posub(FILE_POSITION_TYPE p) {
    end2_ = p;
  }
  void pos(FILE_POSITION_TYPE p) {
    pattern_alignment::end(p);
    end2_ = p;
  }
  void exact_start_bases(unsigned int esb) {
    lconst_ = esb;
  }
  void exact_end_bases(unsigned int eeb) {
    rconst_ = eeb;
  }
  void eos(char ch) {
    eos_ = ch;
  }
  char eos() {
    return eos_;
  }
  void kmax(int k) {
    k_ = k;
  }
  void wc(bool wc) {
    wc_ = wc;
  }
  void tn(bool tn) {
    textn_ = tn;
  }
  void indels(bool id) {
    indels_ = id;
  }
  void dna_mut(bool dm) {
    dna_mut_ = dm;
  }
  void trans(bool t) {
    trans_ = t;
  }
  void maxpatlen(long unsigned int mpl) {
    maxpatlen_ = mpl;
  }
  bool align(CharacterProducer &, std::string const & pattern);
  MEMORY_DEBUG(editdist_alignment)
};

class editdist_peptide_alignment: public editdist_alignment {
 private:
  std::string lcontext_;
  std::string rcontext_;
  int context_;
 public:
  editdist_peptide_alignment(FILE_POSITION_TYPE e=0, 
			     FILE_POSITION_TYPE e2=0,
			     unsigned int k=0, char eos='\n', 
			     bool wc=false, bool tn=false, bool id=true, bool dm=false,
			     int lconst=0, int rconst=0, bool yn=false, bool trans=false)
    : editdist_alignment(e,e2,k,eos,wc,tn,id,dm,lconst,rconst,yn,trans), context_(1) {};
  std::string lcontext() const {
    assert(alignment_done_);
    return lcontext_;
  }
  std::string rcontext() const {
    assert(alignment_done_);
    return rcontext_;
  }
  void set_context(int c) {
    context_ = c;
  }
  bool align(CharacterProducer & cp, std::string const & pattern) {
    bool retval = editdist_alignment::align(cp,pattern);
    cp.pos(max(start()-context_,(FILE_POSITION_TYPE)0));
    unsigned int len=context_;
    if (start() < context_) {
      len=start();
    }
    lcontext_ = cp.getstr(len);
    cp.pos(end());
    rcontext_ = cp.getstr((unsigned)context_);
    return retval;
  }
};
#endif
