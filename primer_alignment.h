
#ifndef _IBPEP_primer_alignment_h
#define _IBPEP_primer_alignment_h

#include <assert.h>
#include <iostream>
#include <string>
#include <vector>
#include "char_io.h"
#include "types.h"
#include "alignment_code.h"
#include "memory_debug.h"

class primer_alignment {
 public:
  primer_alignment() : 
    eos_('\n'), maxdist_(1), 
    wc_(false), tn_(false), 
    indels_(true), dna_mut_(false), 
    yesno_(false), 
    alignment_done_(false), maxpatlen_(0), 
    end_defined_(false), matsize_(0), dp_(0), best_(0)
{
    if (!yesno_) {
      stats_.resize(alignment_codes); 
      fill(stats_.begin(),stats_.end(),0); 
    }
  }
  virtual ~primer_alignment() {};
  virtual bool align(CharacterProducer &, std::string const & pattern1, std::string const & pattern2) =0;
  int size() const {
    assert(alignment_done_);
    return alignment_.size();
  };
  void reset() {
    alignment_done_ = false;
    end_defined_ = false;
    if (!yesno_) {
      fill(stats_.begin(),stats_.end(),0);
    }
  }
  alignment_code operator[](int i) const {
    assert(alignment_done_);
    return alignment_[i];
  }
  FILE_POSITION_TYPE start() const {
    assert(alignment_done_);
    return start_;
  }
  FILE_POSITION_TYPE const & end() const {
    assert(end_defined_ || alignment_done_);
    return end_;
  }
  int const & value() const {
    assert(end_defined_ || alignment_done_);
    return val_;
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
    // checkpoint;
    // cerr << stats_[alignment_constraint_violation] << endl;
    if (stats_[alignment_constraint_violation] > 0) {
      // checkpoint;
      return MAXINT;
    } else {
      // checkpoint;
      return stats_[alignment_substitution] + 
	stats_[alignment_substitution_1] + 
	2*stats_[alignment_substitution_2] + 
	3*stats_[alignment_substitution_3] + 
	stats_[alignment_insertion] +
	stats_[alignment_deletion] + 
	3*stats_[alignment_insertion_3] +
	3*stats_[alignment_deletion_3];
    }
  }
  std::string alignment_string() const {
    assert(alignment_done_);
    std::string r("");
    // checkpoint;
    for (int i=0; i<size(); i++) {
      // cerr << i << ": " << (int) (*this)[i] << endl;
      if ((*this)[i] == alignment_equal) {
	r += "|";
      } else if ((*this)[i] == alignment_wildcard_equal) {
	r += "+";
      } else if ((*this)[i] == alignment_substitution) {
	r += "*";
      } else if ((*this)[i] == alignment_substitution_1) {
	r += ".";
      } else if ((*this)[i] == alignment_substitution_2) {
	r += ":";
      } else if ((*this)[i] == alignment_substitution_3) {
	r += "x";
      } else if ((*this)[i] == alignment_insertion) {
	r += "^";
      } else if ((*this)[i] == alignment_deletion) {
	r += "v";
      } else if ((*this)[i] == alignment_constraint_violation) {
	r += "!";
      } else {
	checkpoint;
	abort();
      }
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
  // void write(ostream & os, 
  // FILE_POSITION_TYPE const seq_pos,
  //	     std::string const & pattern1, 
  // 	     std::string const & pattern2, 
  //	     long unsigned int, bool);
  void eos(char ch) {
    eos_ = ch;
  }
  void kmax(int k) {
    maxdist_ = k;
  }
  void wc(bool wc) {
    wc_ = wc;
  }
  void tn(bool tn) {
    tn_ = tn;
  }
  void indels(bool id) {
    indels_ = id;
  }
  void dna_mut(bool dm) {
    dna_mut_ = dm;
  }
  void yesno(bool yn) {
    yesno_ = yn;
  }
  void maxpatlen(long unsigned int mpl) {
    maxpatlen_ = mpl;
  }
  MEMORY_DEBUG(primer_alignment)
 protected:
  bool global_align(char* const & text, unsigned int textlen, 
		    std::string const & pattern,
		    int dirn, unsigned int lmatch, unsigned int rmatch, 
		    int & matchlen, int & value);
  
  std::vector<alignment_code> alignment_;
  std::vector<int> stats_;
  std::string matching_text_;
  FILE_POSITION_TYPE start_;
  FILE_POSITION_TYPE end_;
  int val_;
  bool alignment_done_;
  bool end_defined_;

 protected:
  char eos_;
  int maxdist_;
  bool wc_;
  bool tn_;
  bool indels_;
  bool dna_mut_;
  bool yesno_;
  long unsigned int maxpatlen_;
  unsigned long int matsize_;
  unsigned int *dp_;
  unsigned int *best_;
};

class primer_alignment_2match : public primer_alignment {
private:
  FILE_POSITION_TYPE end1_;
  FILE_POSITION_TYPE end2_;
  int lmatch_;
  int rmatch_;
public:
  primer_alignment_2match(FILE_POSITION_TYPE e1, FILE_POSITION_TYPE e2,
			  int lmatch,  int rmatch) 
    : end1_(e1), end2_(e2), lmatch_(lmatch), rmatch_(rmatch)  {}
  ~primer_alignment_2match() {}
  bool align(CharacterProducer &, 
	     std::string const & pattern1, std::string const & pattern2);
  MEMORY_DEBUG(primer_alignment_2match)
};

class primer_alignment_lmatch : public primer_alignment {
 private:
  char *buffer0_;
  char *buffer1_;
  FILE_POSITION_TYPE bufstart_;
  FILE_POSITION_TYPE bufend_;
  long unsigned int bufsize_;
  FILE_POSITION_TYPE end1_;
  unsigned int lmatch_;
  unsigned int rmatch_;
 public:
  primer_alignment_lmatch(FILE_POSITION_TYPE e1=0, 
			  unsigned int lmatch=0, 
			  unsigned int rmatch=0) 
    : buffer0_(0), buffer1_(0), bufstart_(0), bufend_(0), bufsize_(0),
    end1_(e1), lmatch_(lmatch), rmatch_(rmatch) {}
  ~primer_alignment_lmatch() {
    delete [] buffer0_;
    delete [] buffer1_;
  }
  FILE_POSITION_TYPE pos() const {
    return end1_;
  }
  void pos(FILE_POSITION_TYPE p) {
    end1_ = p;
  }
  unsigned int exact_start_bases() const {
    return lmatch_;
  }
  void exact_start_bases(unsigned int esb) {
    lmatch_ = esb;
  }
  unsigned int exact_end_bases() const {
    return rmatch_;
  }
  void exact_end_bases(unsigned int eeb) {
    rmatch_ = eeb;
  }
  bool align(CharacterProducer &, 
	     std::string const & pattern1, std::string const & pattern2);
  MEMORY_DEBUG(primer_alignment_lmatch)
};

class primer_alignment_rmatch : public primer_alignment {
 private:
  char *buffer0_;
  char *buffer1_;
  FILE_POSITION_TYPE bufstart_;
  FILE_POSITION_TYPE bufend_;
  long unsigned int bufsize_;
  FILE_POSITION_TYPE end2_;
  unsigned int lmatch_;
  unsigned int rmatch_;
 public:
  primer_alignment_rmatch(FILE_POSITION_TYPE e2=0, 
			  unsigned int lmatch=0, 
			  unsigned int rmatch=0) 
    :  buffer0_(0), buffer1_(0), bufstart_(0), bufend_(0), bufsize_(0),
    end2_(e2), lmatch_(lmatch), rmatch_(rmatch) {}
  ~primer_alignment_rmatch() {
    delete [] buffer0_;
    delete [] buffer1_;
  }
  FILE_POSITION_TYPE pos() const {
    return end2_;
  }
  void pos(FILE_POSITION_TYPE p) {
    end2_ = p;
  }
  unsigned int exact_start_bases() const {
    return lmatch_;
  }
  void exact_start_bases(unsigned int esb) {
    lmatch_ = esb;
  }
  unsigned int exact_end_bases() const {
    return rmatch_;
  }
  void exact_end_bases(unsigned int eeb) {
    rmatch_ = eeb;
  }
  bool align(CharacterProducer &, 
	     std::string const & pattern1, std::string const & pattern2);
  MEMORY_DEBUG(primer_alignment_rmatch)
};

#endif
