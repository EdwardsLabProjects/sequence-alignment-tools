
#include "perfposht.h"

perfposht::perfposht(hash & h, bool rc) 
  : h_(h), rc_(rc)
{
  init(h,rc);
}

perfposht::reset() {
  init(h_,rc_);
}

perfposht::init(hash & h, bool rc) {
  h.reset();
  ptr_.clear();
  ptr_.resize(h.size()+1,0);
  total_ = 0;
  while (h.next()) {
    // timestamplu("Hash value: ",h.value());
    if (rc) {
      ptr_[h.value()]++;
    } else {
      ptr_[h.cvalue()]++;
    }
    total_ ++;
  }
  for (int i=1;i<ptr_.size();i++) {
    ptr_[i] += ptr_[i-1];
  }
  for (int i=ptr_.size()-1;i>0;i--) {
    ptr_[i] = ptr_[i-1];
  }
  ptr_[0] = 0;
  assert(total_ == ptr_[ptr_.size()-1]);
  pos_.clear();
  pos_.resize(total,0);
  h.reset();
  while (h.next()) {
    hash_t v = h.value();
    if (!rc) {
      pos_[ptr_[v]] = h.pos();
      ptr_[v]++;
    } else {
      hash_t v = h.cvalue();
      if (h.value() == v) {
	pos_[ptr_[v]] = h.pos();
      } else {
	pos_[ptr_[v]] = -h.pos();
      }
      ptr_[v]++;
    }
  }
  for (int i=ptr_.size()-1;i>0;i--) {
    ptr_[i] = ptr_[i-1];
  }
  ptr_[0] = 0;
}

perfecthashtable::~perfecthashtable() {
  ;
}
