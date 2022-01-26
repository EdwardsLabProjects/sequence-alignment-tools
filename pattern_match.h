
#ifndef _IBPEP_pattern_match_h
#define _IBPEP_pattern_match_h

#include <list>
#include <utility>
#include <iostream>
#include <iomanip>
// #include <ios>
#include "char_io.h"
#include "types.h"
#include "tinylist.t"
#include "sortedvector.t"
#include "pattern_alignment.h"

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

class pattern_list_element {
private:
  std::string pattern_;
  unsigned long length_;
  unsigned long id_;
  int exact_start_bases_;
  int exact_end_bases_;
public:
  pattern_list_element() : pattern_(""), length_(0), id_(0), 
    exact_start_bases_(0), exact_end_bases_(0) {}
  pattern_list_element(std::string const & p, unsigned long const & i,
		       int esb=0, int eeb=0) 
    : pattern_(p), length_(0), id_(i), exact_start_bases_(esb), exact_end_bases_(eeb)  {
    length_ = pattern_.length();
  }
  ~pattern_list_element() {}
  std::string const & pattern() const {
    return pattern_;
  }
  void pattern(std::string const & p) {
    pattern_ = p;
    length_ = p.length();
  }
  unsigned long int const & id() const {
    return id_;
  }
  void id(unsigned long int i) {
    id_ = i;
  }
  unsigned long int const & length() const {
    return length_;
  }
  int exact_start_bases() const {
    return exact_start_bases_;
  }
  int exact_end_bases() const {
    return exact_end_bases_;
  }
  void exact_start_bases(int esb) {
    exact_start_bases_ = esb;
  }
  void exact_end_bases(int eeb) {
    exact_end_bases_ = eeb;
  }
  void read(istream & is) {
    is >> id_ >> pattern_;
  }
  void write(ostream & os) const {
    os << id_ << " " << pattern_;
  }
  bool operator==(pattern_list_element const & a) const {
    return ((a.id_ == id_) && (a.pattern_ == pattern_));
  }
};

ostream & operator<<(ostream & os, pattern_list_element const & ple);
istream & operator>>(ostream & is, pattern_list_element & ple);

typedef tinylist<pattern_list_element> pattern_list;

class pattern_hit;

typedef sortedvector<FILE_POSITION_TYPE,std::pair<pattern_list::const_iterator,unsigned char> > pattern_hit_vector;

class PatternMatch {
private:
  long unsigned int next_pattern_id_;
  pattern_list patterns_;
  pattern_list::iterator pit_;
  float pint_;
  float pcur_;
protected:
  pattern_list::const_iterator
    add_pattern_(std::string const & pat, long unsigned int & id,
		 int esb=0, int eeb=0) {
    if (id == 0) {
      id= ++next_pattern_id_;
    }
    if (pit_ == patterns_.end()) {
      pit_ = patterns_.push_front(pattern_list_element(pat,id,esb,eeb));
    } else {
      pit_ = patterns_.insert_after(pit_,pattern_list_element(pat,id,esb,eeb));
    }
    return pit_;
  }
  pattern_list const & patterns() const {
    return patterns_;
  }
 public:
  PatternMatch() { 
    next_pattern_id_=0; 
    pit_ = patterns_.end(); 
    pcur_=0.0; 
    pint_=-1.0; 
  };
  virtual ~PatternMatch();
  virtual long unsigned int add_pattern(std::string const & pat, 
					unsigned long id=0,
					int exact_start_bases=0,
					int exact_end_bases=0)=0; 
  virtual void del_pattern(CharacterProducer & cp, std::string const & pat, 
			   long unsigned int=0) {
    timestamp("Fatal error: del_pattern not implemented.");
    exit(1);
  };
  virtual void lookup_pattern(CharacterProducer & cp, std::string const & pat, 
			      std::list<long unsigned int> &l) {
    timestamp("Fatal error: lookup_pattern not implemented.");
    exit(1);
  };
  virtual void init(CharacterProducer & cp)=0;
  virtual bool find_patterns(CharacterProducer & cp,
			     pattern_hit_vector & pas,
			     long unsigned minka=1)=0;
  virtual void reset()=0;
  virtual void progress_interval(CharacterProducer & cp, float p=1.0) {
    if (p >= 0.0 && p <= 100) {
      pcur_ = cp.progress();
      pint_ = p/100;
    } else {
      pcur_ = 0.0;
      pint_ = -1.0;
    }
  }
  void report_progress(CharacterProducer const & cp) {
    if (pint_ >= 0 && cp.progress() > pcur_) {
      ostrstream ss;
      ss.setf(ios::fixed);
      ss << "Progress:";
      ss << setprecision(1) << setw(5) << cp.progress()*100;
      ss << "%" << ends;
      std::string v(ss.str());
      timestamp(v.c_str());
      pcur_ = cp.progress() + pint_;
    }
  }
};

#endif
