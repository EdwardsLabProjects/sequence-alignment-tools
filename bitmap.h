
#ifndef _BITMAP_H_
#define _BITMAP_H_

#include "bits.h"
#include <vector>
#include <sstream>
#include "util.h"
#include <string>

class bitmap {
 private:
  std::vector<bool> x_;
  uint32 size_;
  uint32 set_;
  uint32 unset_;
  uint32 runs_;
  int32 last_false_;
  int32 last_true_ed_;
  int32 last_true_st_;
  bool nocheck;
 public:
  static std::ostream & encode(std::ostream & os, uint32 x, bool ascii) {
    if (ascii) {
      os << x << std::endl;
      return os;
    } else {
      /* std::ostringstream ss;
      encode(ss,x);
      std::istringstream ss1(ss.str());
      uint32 x1=0;
      decode(ss1,x1);
      assert(x == x1); */
      return encode(os,x);
    }
  }
  static std::ostream & encode(std::ostream & os, uint32 x) {
    char b;
    while (true) {
      b = x & 127;
      x >>= 7;
      if (x == 0) {
	os.put(b);
	// os << b;
	return os;
      }
      b |= 128;
      os.put(b);
      // os << (b | 128);
      // x -= 1;
    }
  }
  static std::istream & decode(std::istream & is, uint32 & x, bool ascii) {
    std::string str;
    if (ascii) {
      is >> x; 
      getline(is,str);
      return is;
    } else {
      return decode(is,x);
    }
  }
  static std::istream & decode(std::istream & is, uint32 & x) {
    char b;
    unsigned int s=1;
    x = 0;
    while (true) {
      is.get(b);
      // fprintf(stderr," %u\n",b);
      // is >> b;
      x += (b & 127) * s;
      if ((b&128) == 0) {
	return is;
      }
      s <<= 7;
      // x += s;
    }
  }
 public:
  bitmap(uint32 size=1, bool def=false) : x_(size,def) {
    if (def) {
      set_=size;
    } else {
      set_=0;
    }
    unset_ = size-set_;
    runs_ = 1;
    size_ = size;
    last_false_ = -1;
    last_true_st_ = -1;
    last_true_ed_ = -1;
    nocheck = false;
  }
  ~bitmap() {};
  inline uint32 nset() const {
    return set_;
  }
  inline uint32 nunset() const {
    return unset_;
  }
  void set_nocheck(bool v=true) {
    nocheck = v;
  }
  void unset_nocheck() {
    nocheck = false;
  }
  inline uint32 size() const {
    return size_;
  }
  inline bool operator[](uint32 i) const {
    return x_[i];
  }
  inline bool get(uint32 i) const {
    return (*this)[i];
  }
  void set(uint32 i, bool v=true) {
    assert(i >= 0 && i < size_);
    if (x_[i] == v) {
      return;
    }
    x_[i] = v;
    if (v) {
      set_++;
      unset_--;
    } else {
      unset_++;
      set_--;
    }
    if (i > 0 && i < size_-1) {
      if (v == x_[i-1] && v == x_[i+1]) {
	runs_ -= 2;
      } else if (v != x_[i-1] && v != x_[i+1]) {
	runs_ += 2;
      } 
    } else if (i == 0) {
      if (v == x_[1]) {
	runs_ -= 1;
      } else /* if (v != x[1]) */ {
	runs_ += 1;
      }
    } else /* if (i == size_-1) */ {
      if (v == x_[size_-2]) {
	runs_ -= 1;
      } else /* if (v != x[size_-2]) */ {
	runs_ += 1;
      }
    }
  }
  inline void unset(uint32 i) {
    set(i,false);
  }
  inline bool toggle(uint32 i) {
    if (get(i)) {
      unset(i);
    } else {
      set(i);
    }
    return get(i);
  }
  bool all(uint32 i, uint32 w) {
    if (last_false_ >= i && last_false_ < i+w) {
      return false;
    } else if (last_true_st_ <= i && last_true_ed_+1 == i+w && get(i+w-1)) {
      last_true_ed_ = i+w;
      return true;
    }
    for (int j=w-1;j>=0;j--) {
      if (!get(i+j)) {
	last_false_ = i+j;
	return false;
      }
    }
    last_false_ = -1;
    last_true_ed_ = i+w;
    last_true_st_ = i;
    return true;
  }
  bool all0(uint32 i, uint32 w) {
    for (int j=w-1;j>=0;j--) {
      if (!get(i+j)) {
	return false;
      }
    }
    return true;
  }
  bool all1(int32 i, uint32 w) {
    bool a0 = all0(i,w);
    bool a =  all(i,w);
    assert(a0 == a);
    return a0;
  }
  bitmap & operator~() {
    for (int i=0;i<size_;i++) {
      x_[i] = !x_[i];
    }
    set_ = size_-set_;
    unset_ = size_-unset_;
    // runs_ unchanged
    assert(check());
    return *this;
  }
  bitmap & operator &=(bitmap const & y_) {
    assert(y_.size() == size_);
    for (int i=0;i<size_;i++) {
      set(i,(x_[i] && y_[i]));
    }
    assert(check());
    return *this;
  }
  bitmap & operator |=(bitmap const & y_) {
    assert(y_.size() == size_);
    for (int i=0;i<size_;i++) {
      set(i,(x_[i] || y_[i]));
    }
    assert(check());
    return *this;
  }
  bool check() {
    if (! nocheck) {
      // Check that set_, unset_, runs_ are correct.
      uint32 s=0;
      uint32 u=0;
      uint32 r=0;
      bool rstate=!get(0);
      for (int i=0;i<size_;i++) {
	if (get(i)) {
	  s++;
	} else {
	  u++;
	}
	if (get(i) != rstate) {
	  rstate=get(i);
	  r++;
	}
      }
      assert(s == set_);
      assert(u == unset_);
      assert(r == runs_);
    }
    return true;
  }
  void runs_internal(std::vector<uint32> & runlength) const {
    runlength.clear();
    runlength.reserve(runs_);
    std::vector<bool>::const_iterator curpos=x_.begin();
    std::vector<bool>::const_iterator theend=x_.end();
    std::vector<bool>::const_iterator nextpos;
    bool lookfor=true;
    while (curpos != theend) {
      nextpos = find(curpos,theend,lookfor);
      runlength.push_back(nextpos-curpos);
      lookfor = !lookfor;
      curpos = nextpos;
    }
  }
  typedef std::vector<std::pair<FILE_POSITION_TYPE,uint32> > runs_t;
  void runs(runs_t & runs, bool set=true) const {
    runs.clear();
    std::vector<uint32> runlength;
    runs_internal(runlength);
    runs.reserve(runlength.size()/2);
    FILE_POSITION_TYPE pos=0;
    bool sense=false;
    for (uint32 i=0;i<runlength.size();i++) {
      uint32 l=runlength[i];
      if (sense==set && l > 0) {
	runs.push_back(make_pair(pos,l));
      }
      sense = !sense;
      pos += l;
    }
  }
  std::ostream & write(std::ostream & os, bool ascii=false) const {
    std::vector<uint32> runlength;
    runs_internal(runlength);
    if (ascii) {
      os << "ASCII RUN LENGTHS START\n";
    } else {
      os << "BINARY RUN LENGTHS START\n";
    }
    encode(os,size_,/*ascii*/true);
    encode(os,set_,/*ascii*/true);
    encode(os,unset_,/*ascii*/true);
    encode(os,runs_,/*ascii*/true);
    encode(os,runlength.size(),ascii);
    for (int i=0;i<runlength.size();i++) {
      encode(os,runlength[i],ascii);
    }
    if (ascii) {
      os << "ASCII RUN LENGTHS END\n";
    } else {
      os << "BINARY RUN LENGTHS END\n";
    }
    return os;
  }
  std::istream & read(std::istream & is) {
    bool ascii=false;
    std::string line;
    assert(getline(is,line));
    if (line == "ASCII RUN LENGTHS START") {
      ascii = true;
    } else if (line == "BINARY RUN LENGTHS START") {
      ascii = false;
    } else {
      assert(0);
    }
    
    uint32 p,s,l,r;
    bool v = false;
    p = 0;
    r = 0;

    decode(is,size_,/*ascii*/true);
    x_.clear();
    x_.resize(size_,false);

    decode(is,set_,/*ascii*/true);
    decode(is,unset_,/*ascii*/true);
    decode(is,runs_,/*ascii*/true);

    decode(is,s,ascii);
    for (uint32 i=0;i<s;i++) {
      decode(is,l,ascii);
      if (l != 0) {
	r ++;
      }
      if (v) {
	for (uint32 j=0;j<l;j++) {
	  x_[p+j] = v;
	}
      }
      v = !v;
      p += l;
    }
    
    // cerr << "p " << " " << size_ << endl;
    assert(p == size_);
    assert(set_ + unset_ == size_);
    assert(r == runs_);
    
    assert(getline(is,line));
    if (ascii) {
      assert(line == "ASCII RUN LENGTHS END");
    } else {
      assert(line == "BINARY RUN LENGTHS END");
    }
    check();
    return is;
  }
};

#endif
