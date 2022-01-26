// -*- c++ -*-

#ifndef _CHAR_IO_T_
#define _CHAR_IO_T_

#include <iostream>
#include <fstream>
#include "util.h"
#include "types.h"
#include "math.h"

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

#define BUFTYPE bigword

template <class T>
class Compressed : public T {
private:
  char *chmap_;
  int *invchmap_;
  unsigned int chmapsize_;
  BUFTYPE *buffer_;
  BUFTYPE *bufp_;
  unsigned int ucharsperbuftype_;
  unsigned int bufsize_;
  unsigned int bufbits_;
  unsigned int bits_;
  unsigned int buftypebits_;
  unsigned int bufpos_;
  unsigned int bitpos_;
  unsigned int bitsperbuf_;
  BUFTYPE *mask0_;
  int *shift0_;
  BUFTYPE *mask1_;
  int *shift1_;
  bool eof_;
  FILE_POSITION_TYPE pos_;
public:
  Compressed(std::string const & filename, char eos_char='\n', int frame=0) 
    : T(change_extn(filename,"sqz").c_str(),eos_char,frame) {
    // timestamp("Begining of CompressedChars constructor...");
    eof_ = false;
    pos_ = 0;
    chmap_ = new char[256];
    std::string chmapfn = change_extn(filename,"tbz");
    ifstream chmapfile(chmapfn.c_str());
    chmapfile.read(chmap_,256);
    chmapsize_=chmapfile.gcount();
    chmapfile.close();
    invchmap_ = new int[256];
    for (int i=0;i<256;i++) {
      invchmap_[i] = -1;
    }
    for (unsigned int i=0;i<chmapsize_;i++) {
      invchmap_[chmap_[i]] = i;
    }
    bits_=1;
    while ( (((unsigned)1) << bits_) < chmapsize_ ) {
      bits_++;
    } 
    buftypebits_ = 8*sizeof(BUFTYPE);
    bufsize_ = least_common_multiple(bits_,buftypebits_)/buftypebits_;
    bufbits_ = buftypebits_*bufsize_;
    bitsperbuf_ = bufbits_/bits_;
    ucharsperbuftype_ = sizeof(BUFTYPE)/sizeof(unsigned char);
    buffer_ = new BUFTYPE[bufsize_];
    bitpos_ = 0;
    mask0_ = new BUFTYPE[buftypebits_*bufsize_];
    mask1_ = new BUFTYPE[buftypebits_*bufsize_];
    shift0_ = new int[buftypebits_*bufsize_];
    shift1_ = new int[buftypebits_*bufsize_];
    for (unsigned int i=0;i<buftypebits_*bufsize_;i++) {
      mask0_[i] = 0;
      mask1_[i] = 0;
      shift0_[i] = 0;
      shift1_[i] = 0;
    }
    for (unsigned int i=0;i<buftypebits_*bufsize_;i+=bits_) {
      int shift = (buftypebits_-bits_-(i%buftypebits_));
      BUFTYPE mask = ((1 << bits_)-1);
      if (shift>0) {
	mask0_[i] = mask << shift;
      } else if (shift<0) {
	mask0_[i] = mask >> -shift;
      } else {
	mask0_[i] = mask;
      }
      shift0_[i] = shift;
      if (i/buftypebits_ < (i+bits_-1)/buftypebits_) {
	shift = (buftypebits_ - (bits_+(i%buftypebits_)-buftypebits_));
	if (shift>0) {
	  mask1_[i] = mask << shift;
	} else if (shift<0) {
	  mask1_[i] = mask >> -shift;
	} else {
	  mask1_[i] = mask;
	}
	shift1_[i] = shift;
      }
    }
  }
  inline bool eof() const {
    return eof_;
  }
  unsigned char getnch() {
    // checkpoint;
    if (bitpos_ == 0) {	
      if (T::eof()) {
	eof_= true;
	return (unsigned char)EOF;
      } else {
	BUFTYPE *p(buffer_);
	for (unsigned int i=0;i<bufsize_;i++) {
	  for (unsigned int j=0;j<ucharsperbuftype_;j++) {
	    (*p) <<= 8;
	    (*p) |= T::getnch();
	  }
	  p++;
	}
      }
      bufp_ = buffer_;
    }
    BUFTYPE val;
    int shift0(shift0_[bitpos_]);
    if (shift0>0) {
      val = (((*bufp_)&mask0_[bitpos_]) >> shift0);
    } else if (shift0<0) {
      val = (((*bufp_)&mask0_[bitpos_]) << -shift0);      
    } else {
      val = ((*bufp_)&mask0_[bitpos_]);            
    }
    if (mask1_[bitpos_]) {
      bufp_++;
      val |= (((*bufp_)&mask1_[bitpos_]) >> shift1_[bitpos_]);
    }
    bitpos_ += bits_;
    if (bitpos_ == bufbits_) {
      bitpos_ = 0;
      if (T::eof()) {
	eof_ = true;
      }
    }
    pos_++;
    // cerr << val << endl;
    return val;
  }
  inline char getch() {
    return (char) chmap_[getnch()];
  }
  inline char ch(unsigned char nch) {
    return chmap_[nch];
  }
  inline int nch(char ch) {
    return invchmap_[ch];
  }
  inline FILE_POSITION_TYPE pos() const {
    // checkpoint;
    return pos_;
  }
  void pos(FILE_POSITION_TYPE p) {
    // checkpoint;
    // cerr << ((eof_)?"eof_ true":"eof_ false") << endl;
    // cerr << ((T::eof())?"T::eof() true":"T::eof() false") << endl;
    // cerr << p << endl; 
    // cerr << bitsperbuf_ << endl;
    // cerr << (p/bitsperbuf_)*bufsize_*ucharsperbuftype_ << endl;
    // cerr << T::pos() << endl;
    T::pos((p/bitsperbuf_)*bufsize_*ucharsperbuftype_);
    // cerr << T::pos() << endl;
    pos_ = p - (p%bitsperbuf_);
    // cerr << pos_ << endl;
    // cerr << ((eof_)?"eof_ true":"eof_ false") << endl;
    // cerr << ((T::eof())?"T::eof() true":"T::eof() false") << endl;
    bufpos_ = 0;
    bitpos_ = 0;
    // checkpoint;
    while (pos_ < p) {
      // cerr << pos_ << endl;
      getnch();
    }
    // cerr << ((eof_)?"eof_ true":"eof_ false") << endl;
    // cerr << ((T::eof())?"T::eof() true":"T::eof() false") << endl;
    eof_ = T::eof();
    // cerr << ((eof_)?"eof_ true":"eof_ false") << endl;
    // cerr << ((T::eof())?"T::eof() true":"T::eof() false") << endl;
  }
  inline void reset() {
    T::reset();
    pos_ = 0;
    bufpos_ = 0;
    bitpos_ = 0;
    eof_ = false;
  }
  inline unsigned int size() const {
    return chmapsize_;
  }
  FILE_POSITION_TYPE length() const {
    return (FILE_POSITION_TYPE)ceil((8.0/bits_)*T::length());
  }
  float progress() const {
    return T::progress();
  }
  ~Compressed() {
    delete [] chmap_;
    delete [] invchmap_;
    delete [] buffer_;
    delete [] mask0_;
    delete [] mask1_;
    delete [] shift0_;
    delete [] shift1_;
  }
};

template <class T>
class Normalized : public T {
private:
  unsigned int chmapsize_;
  int *invchmap_;
  char *chmap_;
public:
  Normalized(std::string const & filename, char eos_char='\n', int frame=0,int sizehint=0) 
    : T(change_extn(filename,"sqn").c_str(),eos_char,frame,sizehint) {
    chmap_ = new char[256];
    std::string mapfn = change_extn(filename,"tbl");
    ifstream chmapfile(mapfn.c_str());
    chmapfile.read(chmap_,256);
    chmapsize_=chmapfile.gcount();
    chmapfile.close();
    invchmap_ = new int[256];
    for (int i=0;i<256;i++) {
      invchmap_[i] = -1;
    }
    for (unsigned int i=0;i<chmapsize_;i++) {
      invchmap_[chmap_[i]] = i;
    }
  }
  inline bool eof() const {
    return T::eof();
  }
  inline char ch(unsigned char nch) {
    return chmap_[nch];
  }
  inline int nch(char ch) {
    return invchmap_[ch];
  }
  inline char getch() {
    return chmap_[T::getch()];
  }
  inline unsigned char getnch() {
    return T::getch();
  }
  inline FILE_POSITION_TYPE pos() const {
    return T::pos();
  }
  inline void pos(FILE_POSITION_TYPE p) {
    T::pos(p);
  }
  inline void reset() {
    T::reset();
  }
  inline char const * c_str() const { return T::c_str(); }
  inline char const * filename() const { return T::filename(); }
  ~Normalized() {
    delete [] chmap_;
    delete [] invchmap_;
  }
  inline unsigned int size() const {
    return chmapsize_;
  }
  FILE_POSITION_TYPE length() const {
    return T::length();
  }
  float progress() const {
    return T::progress();
  }
};

template <class T>
class Mapped : public T {
private:
  char *chmap_;
public:
  Mapped(std::string const & filename, char eos_char='\n', int frame=0,int sizehint=0) 
    : T(filename,eos_char,frame) {
    chmap_ = new char[128];
    for (int i=0;i<128;i++) {
      chmap_[i] = i;
    }
  }
  inline void mapto(char from, char to) {
    if (T::nch(from) >= 0) {
      chmap_[T::nch(from)] = T::nch(to);
    }
  }
  inline bool eof() const {
    return T::eof();
  }
  inline char ch(unsigned char nch) {
    return T::ch(nch);
  }
  inline int nch(char ch) {
    return chmap_[T::nch(ch)];
  }
  inline char getch() {
    return T::getch();
  }
  inline unsigned char getnch() {
    return chmap_[T::getnch()];
  }
  inline FILE_POSITION_TYPE pos() const {
    return T::pos();
  }
  inline void pos(FILE_POSITION_TYPE p) {
    T::pos(p);
  }
  inline void reset() {
    T::reset();
  }
  inline char const * c_str() const { return T::c_str(); }
  inline char const * filename() const { return T::filename(); }
  ~Mapped() {
    delete [] chmap_;
  }
  inline unsigned int size() const {
    return T::size();
  }
  FILE_POSITION_TYPE length() const {
    return T::length();
  }
  float progress() const {
    return T::progress();
  }
};

template <class T>
class Translated : public T {
private:
  int frame_in_;
  int frame_;
  bool eof_;
  FILE_POSITION_TYPE frame_end_pos_[6];
  FILE_POSITION_TYPE pos_;
  unsigned int chmapsize_;
  int *invchmap_;
  char *chmap_;
  char eos_char_;
  char codonid_;
public:
  Translated(std::string const & filename, char eos_char='\n', int frame = 0, int sizehint=0) 
    : T(filename.c_str(), eos_char, frame) {
    char locmap[23] = { 'A','C','D','E','F','G','H','I',
			'K','L','M','N','P','Q','R','S',
			'T','V','W','X','Y','*','$' };
    frame_in_ = frame;
    frame_ = 0;
    
    eof_ = false;
    pos_ = 0;
    for (int i=0;i<6;i++) {
      frame_end_pos_[i] = 0;
    }
    chmap_ = new char[23];
    chmapsize_ = 23;
    for (int i=0;i<23;i++) {
      chmap_[i] = locmap[i];
    }
    chmap_[22] = eos_char;
    eos_char_ = eos_char;
    invchmap_ = new int[256];
    for (int i=0;i<256;i++) {
      invchmap_[i] = -1;
    }
    for (unsigned int i=0;i<chmapsize_;i++) {
      invchmap_[chmap_[i]] = i;
    }
  }
  inline bool eof() const {
    return eof_;
  }
  inline char ch(unsigned char nch) {
    return chmap_[nch];
  }
  inline int nch(char ch) {
    return invchmap_[ch];
  }
  inline char getbasech() {
    return T::getch();
  }
  inline char getch() {
    if (T::eof()) {
      frame_end_pos_[frame_] = pos_;
      if ((frame_in_ == 0 && frame_ == 5) ||
	  (frame_in_ == 4 && frame_ == 2)) {
	eof_ = true;
      } else {
	T::reset();
	frame_++;
	for (int i=0;i<(frame_%3);i++) {
	  T::getnch();
	}
      }
    }
    char codon[3] = {0};
    for (int i=0;i<3;i++) {
      if (T::eof()) 
	break;
      codon[i]=T::getch();
    }
    // cout << " " << codon << " ";
    pos_++;
    if (codon[2] == eos_char_ || codon[2] == 0) {
      /* if (codon[0] != eos_char_) {
	T::getnch();
      }
      if (codon[1] != eos_char_) {
	T::getnch();
      }
      for (int i=0;i<(frame_%3);i++) {
	T::getnch();
	}*/
      return eos_char_;
    } else if (codon[0] == eos_char_) {
      return eos_char_;
    } else {
      return trans_codon(frame_,codon,&codonid_);
    }
  }
  inline char getcodonid() {
    return codonid_;
  }
  inline unsigned char getnch() {
    return invchmap_[this->getch()];
  }
  FILE_POSITION_TYPE pos() const {
    return pos_;
  }
  FILE_POSITION_TYPE basepos() const {
    return T::pos();
  }
  void getbasepos(FILE_POSITION_TYPE p, FILE_POSITION_TYPE & p1, int & f) const {
    // checkpoint;
    // cerr << p << " [";
    f = -1;
    for (int i=0;i<(frame_in_==0?6:3);i++) {
      // cerr << frame_end_pos_[i] << " ";
      if (p < frame_end_pos_[i] || frame_end_pos_[i] == 0 || (i == (frame_in_==0?6:3)-1 && p == frame_end_pos_[i]+1)) {
	f = i;
	break;
      }
    }
    // cerr << "]: " << f;
    assert(f >= 0);
    p1 = ((p-((f==0)?0:frame_end_pos_[f-1]))*3 + f%3);
    // cerr << " " << p1 << endl;
  }
  void pos(FILE_POSITION_TYPE p) {
    FILE_POSITION_TYPE bpos;
    getbasepos(p,bpos,frame_);
    pos_ = p;
    T::pos(bpos);
    eof_ = false;
  }
  inline void reset() {
    T::reset();
    eof_ = false;
  }
  ~Translated() {
    delete [] chmap_;
    delete [] invchmap_;
  }
  inline unsigned int size() const {
    return 23;
  }
  FILE_POSITION_TYPE length() const {
    return 2*T::length();
  }
  float progress() const {
    return (frame_ + T::progress())/(frame_in_==0?6.0:3.0);
  }
};

#endif
