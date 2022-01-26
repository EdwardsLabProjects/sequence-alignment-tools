
#ifndef _IBPEP_sortedvector_h_
#define _IBPEP_sortedvector_h_

#include <vector>
#include <utility>
#include <algorithm>
#include "memory_debug.h"
#include "util.h"

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

template <class K, class V> 
class sortedvector {
 public:
  class iterator;
  class const_iterator;
 public:
  class element {
  private:
    std::pair<K,V> p_;
  public:
    element(K const & k=K(), V const & v=V()) : p_(k,v) {};
    K const & key() const{
      return p_.first;
    };
    K & key(){
      return p_.first;
    };
    V const & value() const {
      return p_.second;
    };
    V & value() {
      return p_.second;
    };
    bool operator<(element const & a) const {
      if (p_.first < a.p_.first) {
	return true;
      } else {
	return false;
      }     
    };
    bool operator==(element const & a) const {
      if (p_.first != a.p_.first) {
	return false;
      } else {
	return (p_.second == a.p_.second);
      }
    };
    static bool lt(element const & a,element const & b) {
      if (a.p_.first < b.p_.first) {
	return true;
      } else if (a.p_.first > b.p_.first) {
	return false;
      } else {
	return (a.p_.second < b.p_.second);
      }    
    };
    static bool valstrictlt(element const & a,element const & b) {
      if (a.p_.second < b.p_.second) {
	return true;
      } else if (a.p_.second > b.p_.second) {
	return false;
      } else {
	return (a.p_.first < b.p_.first);
      }    
    };
    static bool vallt(element const & a,element const & b) {
      if (a.p_.second < b.p_.second) {
	return true;
      } else {
	return false;
      }   
    };
    friend class iterator;
    friend class const_iterator;
    friend class sortedvector;
    MEMORY_DEBUG(sortedvector::element)
  };
 public:
  class iterator {
  private:
    element *elt;
  public:
    iterator() : elt(0) {};
    iterator(element* e) : elt(e) {};
    iterator(iterator const & it) : elt(it.elt) {};
    iterator & operator=(iterator const & it) { 
      elt = it.elt; 
      return (*this);
    };
    element & operator*() const {
      return (*elt);
    }
    element * operator->() const {
      return elt;
    }
    iterator & operator++() {
      elt++;
      return (*this);
    }
    iterator & operator--() {
      elt--;
      return (*this);
    }
    iterator & operator+=(long int inc) {
      elt+=inc;
      return (*this);
    }
    iterator & operator-=(long int inc) {
      elt-=inc;
      return (*this);
    }
    iterator operator+(long int inc) const {
      iterator retval(*this); 
      retval+=inc;
      return retval;
    }
    long int operator-(iterator const & it) const {
      return (elt - it.elt);
    }
    long int operator-(const_iterator const & it) const {
      return (elt - it.elt);
    }
    iterator operator-(long int inc) const {
      iterator retval(*this);
      retval+=(-inc);
      return retval;
    }
    bool operator==(const iterator & it) const {
      return (elt == it.elt);
    }
    bool operator!=(const iterator & it) const {
      return (elt != it.elt); 
    }
    bool operator<(const iterator & it) const {
      return ((void*)elt < (void*)it.elt);
    }
    bool operator>(const iterator & it) const {
      return ((void*)elt > (void*)it.elt);
    }
    bool operator>=(const iterator & it) const {
      return ((void*)elt >= (void*)it.elt);
    }
    friend class const_iterator;
    friend class sortedvector;
  };
 public:
  class const_iterator {
  private:
    element *elt;
  public:
    const_iterator() : elt(0) {};
    const_iterator(element *e) : elt(e) {};
    const_iterator(const_iterator const & it) : elt(it.elt) {};
    const_iterator(iterator const & it) : elt(it.elt) {};
    const_iterator & operator=(const_iterator const & it) { 
      elt = it.elt; 
      return (*this);
    };
    element const & operator*() const {
      return *elt;
    }
    element const * operator->() const {
      return elt;
    }
    const_iterator & operator++() {
      elt++;
      return (*this); 
    }
    const_iterator & operator--() {
      elt--;
      return (*this); 
    }
    const_iterator & operator+=(long int inc) {
      elt+=inc;
      return (*this);
    }
    const_iterator & operator-=(long int inc) {
      elt-=inc;
      return (*this);
    }
    const_iterator operator+(long int inc) const {
      const_iterator retval(*this); 
      retval+=inc;
      return retval;
    }
    long int operator-(const_iterator const & it) const {
      return (elt - it.elt);
    }
    const_iterator operator-(long int inc) const {
      const_iterator retval(*this);
      retval+=(-inc);
      return retval;
    }
    bool operator==(const const_iterator & it) const {
      return (elt == it.elt); 
    }
    bool operator!=(const const_iterator & it) const {
      return (elt != it.elt);
    }
    bool operator<(const const_iterator & it) const {
      return ((void*)elt < (void*)it.elt);
    }
    bool operator>(const const_iterator & it) const {
      return ((void*)elt > (void*)it.elt);
    }
    bool operator>=(const const_iterator & it) const {
      return ((void*)elt >= (void*)it.elt);
    }
    friend class iterator;
    friend class sortedvector;
  };
private:
  long unsigned int size_;
  long unsigned int capacity_;
  element *v_;
  const_iterator locate_between(const_iterator const & l, 
				const_iterator const & r,
				const K & k, 
				bool leftmost) const {
    // checkpoint;
    // cerr << k << " " << index(l) << " " << index(r) << endl;
    // Precondition:
    // l->key() < k || l->key == k && l == begin()
    // r->key() > k || r->key == end()
    if ((l != r) && ((l->key() >= k && (l->key() != k || l != begin())) ||
		     (r != end() && r->key() <= k))) {
      // checkpoint;
      // cerr << "located between precondition violated!\n";
      return end();
    } else {
      const_iterator l1=l;
      const_iterator r1=r;
      int length=(r-l);
      const_iterator m1=l1 + (length)/2;
      while ((l1+1) < r1) {
	if (m1->key() < k) {
	  l1 = m1;
	} else if (m1->key() > k) {
	  r1 = m1;
	} else if (leftmost) {
	  r1 = m1;
	} else {
	  l1 = m1;
	}
	length = r1-l1;
	m1=l1 + length/2;
	// cerr << k << " " << index(l1) << " " << index(r1) << endl;
      }
      if (leftmost && l1->key() == k) r1 = l1;
      // if (!leftmost && r1->key() == k) l1 = r1;
      return (leftmost?r1:l1);
    }
  };
  iterator locate_between(iterator const & l, 
			  iterator const & r,
			  const K & k, 
			  bool leftmost) {
    // checkpoint;
    // cerr << k << " " << l->key() << " (" << index(l) << ") " << r->key() << " (" << index(r) << ")" << endl;
    // Precondition:
    // l->key() < k || l->key == k && l == begin()
    // r->key() > k || r->key == end()
    if ((l != r) && ((l->key() >= k && (l->key() != k || l != begin())) ||
		     (r != end() && r->key() <= k))) {
      // checkpoint;
      // cerr << "located between precondition violated!\n";
      return end();
    } else {
      iterator l1=l;
      iterator r1=r;
      int length=(r-l);
      iterator m1=l1 + (length)/2;
      while ((l1+1) < r1) {
	if (m1->key() < k) {
	  l1 = m1;
	} else if (m1->key() > k) {
	  r1 = m1;
	} else if (leftmost) {
	  r1 = m1;
	} else {
	  l1 = m1;
	}
	length = r1-l1;
	m1=l1 + length/2;
	// cerr << k << " " << l1->key() << " (" << index(l1) << ") " << r1->key() << " (" << index(r1) << ")" << endl;
      }
      if (leftmost && l1->key() == k) r1 = l1;
      // if (!leftmost && r1->key() == k) l1 = r1;
      return (leftmost?r1:l1);
    }
  };

  void find_bracket(const_iterator const & c,
		    const_iterator & l,
		    const_iterator & r,
		    const K & k) const {
    // checkpoint;
    unsigned int jump=2;
    if (k > c->key()) {
      l = c;
      if (c != end()) {
	r = c+1;
      } else {
	r = c;
      }
      while (r->key() <= k && r < end()) {
      	if (r->key() < k) l = r;
	if (end() - r < jump) {
	  r = end();
	} else {
	  r += jump;
	  jump *= 2;
	}
      }
    } else if (k < c->key()) {
      if (c != begin()) {
	l = c-1;
      } else {
	l = c; 
      }
      r = c;
      while (l->key() >= k && l > begin()) {
	if (l->key() > k) r = l;
	if (l - begin() < jump) {
	  l = begin();
	} else {
	  l -= jump;
	  jump *= 2;
	}
      }
    } else { /* k == c->key() */
      if (c !=  begin()) {
	l = c-1;
      } else {
	l = c;
      }
      if (c != end()) {
	r = c+1;
      } else {
	r = c;
      }
      while (r->key() <= k && r < end()) {
	if (end() - r < jump) {
	  r = end();
	} else {
	  r += jump;
	  jump *= 2;
	}
      }
      while (l->key() >= k && l > begin()) {
	if (l - begin() < jump) {
	  l = begin();
	} else {
	  l -= jump;
	  jump *= 2;
	}
      }
    }
    //     checkpoint;
    //     cerr << l->key() << " (" << index(l) << ") " 
    // 	 << c->key() << " (" << index(c) << ") "
    // 	 << r->key() << " (" << index(r) << ") "
    // 	 << endl;
  };

  void find_bracket(iterator const & c,
		    iterator & l,
		    iterator & r,
		    const K & k) {
    unsigned int jump=2;
    // checkpoint;
    // cerr << k << " " << c->key() << endl;
    if (k > c->key()) {
      l = c;
      if (c != end()) {
	r = c+1;
      } else {
	r = c;
      }
      while (r->key() <= k && r < end()) {
      	if (r->key() < k) l = r;
	if (end() - r < jump) {
	  r = end();
	} else {
	  r += jump;
	  jump *= 2;
	}
      }
    } else if (k < c->key()) {
      r = c;
      if (c != begin()) {
	l = c-1;
      } else {
	l = c;
      }
      while (l->key() >= k && l > begin()) {
	if (l->key() > k) r = l;
	if (l - begin() < jump) {
	  l = begin();
	} else {
	  l -= jump;
	  jump *= 2;
	}
      }
    } else { /* k == c->key() */
      if (c != begin()) {
	l = c-1;
      } else {
	l = c;
      }
      if (c != end()) {
	r = c+1;
      } else {
	r = c;
      }
      while (r->key() <= k && r < end()) {
	if (end() - r < jump) {
	  r = end();
	} else {
	  r += jump;
	  jump *= 2;
	}
      }
      while (l->key() >= k && l > begin()) {
	if (l - begin() < jump) {
	  l = begin();
	} else {
	  l -= jump;
	  jump *= 2;
	}
      }
    }
    //     checkpoint;
    //     cerr << l->key() << " (" << index(l) << ") " 
    // 	 << c->key() << " (" << index(c) << ") "
    // 	 << r->key() << " (" << index(r) << ") "
    // 	 << endl;
  };
public:
  sortedvector(long unsigned int r=10) : size_(0), capacity_(0), v_(0) {
    reserve(r);
  };
  sortedvector(sortedvector const & sv) : size_(0), capacity_(0), v_(0) {
    capacity_ = sv.capacity_;
    v_ = new element[capacity_];
    size_ = sv.size_;
    for (int i=0;i<size_;i++) {
      v_[i] = sv.v_[i];
    }
  }
  sortedvector & operator=(sortedvector const & sv) {
    capacity_ = sv.capacity_;
    v_ = new element[capacity_];
    size_ = sv.size_;
    for (int i=0;i<size_;i++) {
      v_[i] = sv.v_[i];
    }
    return (*this);
  }
  ~sortedvector() {
    delete [] v_;
  }
  void append(K const & k, V const & v) {
    if (size_ >= capacity_) {
      reserve(capacity_*2);
    }
    if (size_ > 0 && v_[size_-1].key() > k) {
      // checkpoint;
      throw BadAppend();
    } else {
      v_[size_++] = element(k,v);
    }
  };
  void insert(K const & k, V const & v) {
    if (size_ >= capacity_) {
      reserve(capacity_*2);
    }
    if (size_ > 0 && v_[size_-1].key() > k) {
      // checkpoint;
      v_[size_++] = element(k,v);
      normalize();
    } else {
      v_[size_++] = element(k,v);
    }
  };
  void push_back(K const & k, V const & v) {
    if (size_ >= capacity_) {
      reserve(capacity_*2);
    }
    v_[size_++] = element(k,v);
  };
  void normalize_byvalue() {
    std::sort(&v_[0],&v_[size_],element::vallt);
  };
  void normalize_strict_byvalue() {
    std::sort(&v_[0],&v_[size_],element::valstrictlt);
  };
  void normalize() {
    std::sort(&v_[0],&v_[size_]);
  };
  void normalize_strict() {
    std::sort(&v_[0],&v_[size_],element::lt);
  };
  void normalize(bool (*cmp)(element const &,element const &)) {
    std::sort(&v_[0],&v_[size_],cmp);	
  }		
  void uniqueify() {
    if (size() < 2) return;
    normalize_strict();	
    long unsigned int newsize=0;
    iterator it = begin();
    if (it == end()) {
      return;
    }
    iterator it1 = it;
    ++it;
    while (it != end()) {
      if (!((*it) == (*it1))) {
	++it1;
	it1->key() = it->key();
	it1->value() = it->value();
	newsize++;
      } 
      ++it;
    }
    resize(newsize);
  }
  const_iterator locate_first_at_least(K const & k) const {
    // checkpoint;
    if (empty()) throw KeyOutOfRange();
    if (begin()->key() > k) {
      return begin();
    } else {
      const_iterator it = locate_between(begin(),end(),k,true);
      if (it == end()) {
	throw KeyOutOfRange();
      } else {
	return it;
      } 
    }
  };
  iterator locate_first_at_least(K const & k) {
    // checkpoint;
    if (empty()) throw KeyOutOfRange();
    if (begin()->key() > k) {
      return begin();
    } else {
      // checkpoint;
      iterator it = locate_between(begin(),end(),k,true);
      // checkpoint;
      // cerr << index(it) << endl;
      if (it == end()) {
	throw KeyOutOfRange();
      } else {
	return it;
      } 
    }
  };
  iterator finger_locate_first_at_least(iterator const & it0, K const & k) {
    if (empty()) throw KeyOutOfRange();
    iterator l,r;
    if (it0 < begin() || it0 >= end()) {
      throw InvalidFinger();
    }
    // checkpoint;
    // cerr << index(it0) << " " << it0->key() << " "  
    // << k << endl;
    find_bracket(it0,l,r,k);
//     cerr << index(it0) << " " << it0->key() << " "  
// 	 << index(l) << " " << l->key() << " "
// 	 << index(r) << " " << r->key() << " "
// 	 << k << endl;
    iterator it = locate_between(l,r,k,true);
    // checkpoint;
    if (it == end()) {
      // checkpoint;
      throw KeyOutOfRange();
    } else {
      // checkpoint;
      // cerr << index(it) << " " << it->key() << endl;
      return it;
    } 
  }
  ;
  const_iterator finger_locate_first_at_least(const_iterator const & it0, K const & k) const {
    if (empty()) throw KeyOutOfRange();
    const_iterator l,r;
    if (it0 < begin() || it0 >= end()) {
      throw InvalidFinger();
    }
    // checkpoint;
    find_bracket(it0,l,r,k);
    const_iterator it = locate_between(l,r,k,true);
    if (it == end()) {
      throw KeyOutOfRange();
    } else {
      return it;
    } 
  }
  ;
  iterator locate_last_at_most(K const & k) {
    // checkpoint;
    if (empty()) throw KeyOutOfRange();
    // checkpoint;
    // cerr << k << endl;
    // cerr << index(begin()) << endl;
    // cerr << index(end()) << endl;
    iterator it = locate_between(begin(),end(),k,false);
    // cerr << it << endl;
    if (it == end()) {
      throw KeyOutOfRange();
    } else {
      return it;
    } 
  };
  const_iterator locate_last_at_most(K const & k) const {
    // checkpoint;
    if (empty()) throw KeyOutOfRange();
    // checkpoint;
    // cerr << k << endl;
    // cerr << index(begin()) << endl;
    // cerr << index(end()) << endl;
    const_iterator it = locate_between(begin(),end(),k,false);
    // cerr << it << endl;
    if (it == end()) {
      throw KeyOutOfRange();
    } else {
      return it;
    } 
  };
  iterator finger_locate_last_at_most(iterator const & it0, K const & k) {
    if (empty()) throw KeyOutOfRange();
    // checkpoint;
    iterator l,r;
    if (it0 < begin() || it0 >= end()) {
      // checkpoint;
      throw InvalidFinger();
    }
    // checkpoint;
    find_bracket(it0,l,r,k);
    // checkpoint;
    iterator it = locate_between(l,r,k,false);
    // checkpoint;
    if (it == end()) {
      // checkpoint;
      throw KeyOutOfRange();
    } else {
      // checkpoint;
      return it;
    } 
    // checkpoint;
  };
  const_iterator finger_locate_last_at_most(const_iterator const & it0, K const & k) const {
    if (empty()) throw KeyOutOfRange();
    // checkpoint;
    const_iterator l,r;
    if (it0 < begin() || it0 >= end()) {
      // checkpoint;
      throw InvalidFinger();
    }
    // checkpoint;
    find_bracket(it0,l,r,k);
    // checkpoint;
    const_iterator it = locate_between(l,r,k,false);
    // checkpoint;
    if (it == end()) {
      // checkpoint;
      throw KeyOutOfRange();
    } else {
      // checkpoint;
      return it;
    } 
    // checkpoint;
  };
  element & operator[](long unsigned int const & i) {
    return v_[i];
  };
  element const & operator[](long unsigned int const & i) const {
    return v_[i];
  };
  iterator iter(long unsigned int const & i) {
    return &(v_[i]);
  };
  const_iterator iter(long unsigned int const & i) const {
    return &(v_[i]);
  };
  long int index(iterator const & it) const {
    return (it-begin());
  }
  long int index(const_iterator const & it) const {
    return (it-begin());
  }
  iterator begin(){
    return &(v_[0]);
  };
  iterator rbegin(){
    return &(v_[size_-1]);
  };
  iterator end() {
    return &(v_[size_]);
  };
  iterator rend() {
    return &(v_[-1]);
  };
  const_iterator begin() const {
    return &(v_[0]);
  };
  const_iterator end() const{
    return &(v_[size_]);
  };
  const_iterator rbegin() const {
    return &(v_[size_-1]);
  };
  const_iterator rend() const{
    return &(v_[-1]);
  };
  bool empty() const {
    return (size_==0);
  };
  long unsigned int size() const {
    return size_;
  };
  void resize(long unsigned int rs) {
    if (rs >= size_) {
      if (rs > capacity_) {
	reserve(rs);
      }
    } else {
      for (long unsigned int i=rs;i<size_;i++) {
	v_[i].~element();
      }
    }
    size_ = rs;
  };
  void clear() {
    resize(0);
  };
  long unsigned int capacity() const {
    return capacity_;
  };
  void reserve(long unsigned int rs){
    if (rs > capacity_) {
      element *tmp = new element[rs];
      for (long unsigned int i=0;i<size_;i++) {
	tmp[i] = v_[i];
      }
      delete [] v_;
      v_ = tmp;
      capacity_ = rs;
    }
  };
  void write(ostream & os) {
    for (const_iterator i=begin();i!=end(); ++i) {
      os << i->key() << " " << i->value() << endl;
    }
  }
  void read(istream & is) {
    K k;
    V v;
    while (is >> k >> v) {
      push_back(k,v);
    }
    normalize();
  }
  void bwrite(ostream & os) {
    long unsigned int s=size();
    os.write((char*)&s,sizeof(long unsigned int));
    os.write((char*)v_,sizeof(element)*s);
  }
  void bread(istream & is) {
    long unsigned int s;
    is.read((char*)&s,sizeof(long unsigned int));
    reserve(s);
    is.read((char*)v_,sizeof(element)*s);
    size_ = s;
    normalize();
  }
  // Exceptions...
  class BadAppend {};
  class KeyOutOfRange {};
  class InvalidFinger {};
};  

template <class K, class V>
std::ostream & operator<<(std::ostream & os, sortedvector<K,V> const & sv) {
  for (typename sortedvector<K,V>::const_iterator i=sv.begin();i<sv.end(); ++i) {
    os << sv.index(i) << ": (" << i->key() << "," << i->value() << ")" << endl;
  }
  return os;
}

#endif
