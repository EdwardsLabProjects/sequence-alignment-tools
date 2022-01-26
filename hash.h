
#ifndef _HASH_H_
#define _HASH_H_

#include "char_io.h"
#include "bits.h"
#include <string>


typedef uint32 hash_t;
// typedef uint64 hash_t;

class hash {
 protected:
  CharacterProducer & cp_;
 public:
  hash(CharacterProducer & cp) :
    cp_(cp) {};
  virtual ~hash() {};
  virtual void update(unsigned char c)=0;
  virtual uint32 ns() const {
    return 0;
  }
  virtual hash_t size() const=0;
  virtual hash_t value() const=0;
  virtual hash_t rcvalue() const=0;
  virtual hash_t cvalue() const=0;
  virtual uint64 pos() const=0;
  virtual bool next()=0;
  virtual void reset(uint64 p=0)=0;
  virtual string str(hash_t v) const=0;
  virtual unsigned span() const=0;
  virtual unsigned weight() const=0;
  virtual unsigned minspan() const {
    return span();
  }
  virtual bool symmetric() const {
    return true;
  }
};

class contigmod : public hash {
private:
  hash_t h_;
  unsigned a_;
  unsigned w_;
  hash_t mask_;
public:
  contigmod(CharacterProducer & cp, int asize, int weight);
  hash_t size() const;
  hash_t value() const;
  hash_t rcvalue() const;
  hash_t cvalue() const;
  void update(unsigned char c);
  bool next();
  void reset(uint64 p=0);
  uint64 pos() const;
  string str(hash_t v) const;
  unsigned span() const;
  unsigned weight() const;
};

class contigshift : public hash {
private:
  hash_t h_;
  unsigned a_;
  unsigned w_;
  hash_t mask_;
  uint32 ns_;
public:
  contigshift(CharacterProducer & cp, int asize, int weight);
  hash_t size() const;
  hash_t value() const;
  hash_t rcvalue() const;
  hash_t cvalue() const;
  void update(unsigned char c);
  uint32 ns() const;
  bool next();
  void reset(uint64 p=0);
  uint64 pos() const;
  string str(hash_t v) const;
  unsigned span() const;
  unsigned weight() const;
  virtual bool symmetric() const {
    return true;
  }
};

class spaced : public hash {
private:
  hash_t* h_;
  unsigned p_;
  unsigned a_;
  unsigned w_;
  unsigned s_;
  unsigned pw_;
  unsigned ptr_;
  hash_t mask_;
  unsigned *hind_;
  unsigned *hptr_;
public:
  spaced(CharacterProducer & cp, int asize, const char *tstr, signed char period=0);
  ~spaced();
  hash_t size() const;
  hash_t value() const;
  hash_t rcvalue() const;
  hash_t cvalue() const;
  void update(unsigned char c);
  bool next();
  void reset(uint64 p=0);
  uint64 pos() const;
  string str(hash_t v) const;
  unsigned span() const;
  unsigned weight() const;
};

class shiftspaced : public hash {
protected:
  uint64 h0_;
  hash_t h_;
  unsigned a_;
  unsigned w_;
  unsigned s_;
  unsigned size_;
  uint64 *mask_;
  unsigned *shift_;
  unsigned nshift_;
public:
  shiftspaced(CharacterProducer & cp, int asize, const char *tstr,const char* type=NULL);
  ~shiftspaced();
  hash_t size() const;
  hash_t value() const;
  hash_t rcvalue() const;
  hash_t cvalue() const;
  void update(unsigned char c);
  bool next();
  void reset(uint64 p=0);
  uint64 pos() const;
  string str(hash_t v) const;
  string str(hash_t v, unsigned short w0) const;
  unsigned span() const;
  unsigned weight() const;
  virtual bool symmetric() const {
    return true;
  }
};

class asymmetric_shiftspaced : public shiftspaced {
 protected:
  uint64 h0rc_;
  hash_t hrc_;
 public:
  asymmetric_shiftspaced(CharacterProducer & cp, int asize, const char *tstr);
  virtual hash_t rcvalue() const;
  virtual void update(unsigned char c);  
  virtual bool symmetric() const {
    return false;
  }
};

class hashset : public hash {
  unsigned n_;
  unsigned p_;
  unsigned s_;
  hash** h_;
 public:
  hashset(CharacterProducer & cp, int asize, const char *tstr);
  ~hashset();
  hash_t size() const;
  hash_t value() const;
  hash_t rcvalue() const;
  hash_t cvalue() const;
  void update(unsigned char c);
  bool next();
  void reset(uint64 p=0);
  uint64 pos() const;
  string str(hash_t v) const;
  unsigned span() const;
  unsigned weight() const;
  unsigned minspan() const {
    return s_;
  }
  bool symmetric() const;
};

class taghashset : public hash {
  unsigned ab_;
  unsigned tb_;
  unsigned hb_;
  unsigned ts_;
  unsigned n_;
  unsigned tn_;
  unsigned tp_;
  unsigned s_;
  hash** h_;
  hash** t_;
  hash_t mask_;
  hash_t hm_;
  hash_t *tpm_;
 public:
  taghashset(CharacterProducer & cp, int asize, const char *tstr, int tsize=0);
  ~taghashset();
  hash_t size() const;
  hash_t value() const;
  hash_t rcvalue() const;
  hash_t cvalue() const;
  void update(unsigned char c);
  bool next();
  void reset(uint64 p=0);
  uint64 pos() const;
  string str(hash_t v) const;
  unsigned span() const;
  unsigned weight() const;
  unsigned minspan() const {
    return s_;
  }
  bool symmetric() const;
};

hash *hashselect(CharacterProducer & cp, int asize, const char *tstr);
hash *spacedselect(CharacterProducer & cp, int asize, const char *tstr);

#endif

