
#ifndef _hashtable_t_
#define _hashtable_t_

// Parameters of the hash table implementation
// Max number of T's that can be stored
// Chain chunk list increment (additive?)
// Hash size.

#include <vector>
#include <assert.h>
template <I,n,m>
struct iptrdata {
  I ptr:n
  I data:m
  iptrdata(): ptr(0), data(0) {}
  static I maxptr;() {
    return ((((I)1)<<(n-1))-1);
  }
  static long unsigned int maxdata() {
    return ((((I)1)<<(m-1))-1);
  }
}

// Note that n+m == sizeof(I)

template <I,n,m,T>
class hasharray {
  vector<iptrdata<I,n,m> > hashmap_;
  vector<iptrdata<I,n,m> > chains_;
  vector<T> elements_;
  I occupancy_[3];
  enum occupancy_index { zero=0, one=1, many=2 };
  enum chain_data { eltptr=0, chptr=1, lastpos=2; };
 public:
  hasharray(I hashsize, I size) {
    assert(hashsize >= size);
    assert(size <= iptrdata<I,n,m>::maxptr());
    hashmap_.resize(hashsize);
    chains_.reserve(1024);
    elements_.reserve(size,T());
    occupancy_[zero] = hashsize;
    occupancy_[one] = occupancy_[many] = 0;
  }
  void insert(I hash, T const & value) {
    iptrdata<I,n,m> & slot = hashmap_[hash];
    if (slot.data == 0) {
      elements_.push_back(value);
      slot.ptr = elements_.size()-1;
      occupancy_[zero]--;
      occupancy_[one]++;
      return;
    } else if (slot.data == 1) {
      if (elements_[slot.ptr] == value) {
	return;
      }
      slot.data = 2;
      I ptr = slot.ptr;
      elements_.push_back(value);
      I ptr1 = elements_.size()-1;
      I chptr = chains_.size();
      chains_.resize(chains_.size()+increment);
      chains_[chptr].ptr = ptr;
      chains_[chptr].data = chain_data::eltptr;
      chains_[chptr+1].ptr = ptr1;
      chains_[chptr+1].data = chain_data::eltptr;
      chains_[chains.size()-1].data = chain_data::lastpos;
      slot.ptr = chptr;
    } else {
      I count=slot.data;
      I i=0;
      I chptr = slot.ptr;
      I chptr0 = chptr;
      while (i<count) {
	iptrdata<I,n,m> & chslot = chains_[chptr];
	if (chslot.data == chain_data::count) {
	  count = chslot.ptr;
	}
	if (chslot.data == chain_data::chptr) {
	  chptr = chslot.ptr;
	  chptr0 = chptr;
	  continue;
	}
	if (elements_[chslot.ptr] == value) {
	  found = true;
	  return;
        }
	i++;
      }
    }
      
  }
  const T & get(I hash) const {
    
  }
  bool has_key(I hash) const {

  }

}



#endif
