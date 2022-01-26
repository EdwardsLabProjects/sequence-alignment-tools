
#ifndef _PERFECTHASHTABLE_H_
#define _PERFECTHASHTABLE_H_

#include <vector>
#include "hash.h"
#include "bitmap.h"

class perfposht {
 private:
  int32 ls_;
  int32 le_;
  vector<uint32> ptr_;
  vector<int32>  pos_;
 public:
  perfposht();
  perfposht(hash & h, bool rc, bool cannon);
  perfposht(hash & h, bool rc, bool cannon, FILE_POSITION_TYPE, int, bitmap &);
  void init(hash & h, bool rc, bool cannon, FILE_POSITION_TYPE=0, int=0, bitmap* =0);
  ~perfposht();
  void set_invalid(hash_t v);
  bool is_valid(hash_t v);
  bool lookup(hash_t v);
  bool next();
  int32 value();
  void set_invalid();
  // void check(CharacterProducer & cp);
};

#endif
