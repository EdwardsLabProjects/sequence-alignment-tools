
#ifndef _tiny_list_h_
#define _tiny_list_h_

#include "memory_debug.h"
#include <iostream>

template <class T>
class tinylist {
 public:
  class iterator;
  class const_iterator;

 public:
  class tinylist_elt {
  private:
    T data;
    tinylist_elt* next;
  public:
    tinylist_elt(const T & t, tinylist_elt * const & n)
      : data(t),next(n) {};
    friend class iterator;
    friend class const_iterator;
    friend class tinylist;
    MEMORY_DEBUG(tinylist::tinylist_elt)
  };

 public:
  class iterator {
  private:
    tinylist_elt* elt;
  public:
    iterator() : elt(0) {};
    iterator(tinylist_elt* e) : elt(e) {};
    iterator(iterator const & it) : elt(it.elt) {};
    iterator & operator=(iterator const & it) { 
      elt = it.elt; 
      return (*this);
    };
    std::ostream & write(std::ostream & os) const {
	os << (void*)elt << std::endl;
	return os;
    }
    T & operator*() const {
      return elt->data;
    }
    T * operator->() const {
      return (&(elt->data));
    }
    iterator & operator++() {
      elt = elt->next;
      return (*this); 
    }
    bool operator==(const iterator & it) const {
      return (elt == it.elt);
    }
    bool operator!=(const iterator & it) const {
      return (elt != it.elt); 
    }
    friend class tinylist;
    friend class const_iterator;
    MEMORY_DEBUG(tinylist::iterator)
  };
  class const_iterator {
  private:
    tinylist_elt* elt;
  public:
    const_iterator() : elt(0) {};
    const_iterator(tinylist_elt* e) : elt(e) {};
    const_iterator(const_iterator const & it) : elt(it.elt) {};
    const_iterator(iterator const & it) : elt(it.elt) {};
    std::ostream & write(std::ostream & os) const {
	os << (void*)elt << std::endl;
	return os;
    }
    const_iterator & operator=(const_iterator const & it) { 
      elt = it.elt; 
      return (*this);
    };
    T const & operator*() const {
      return elt->data;
    }
    T const * operator->() const {
      return (&(elt->data));
    }
    const_iterator & operator++() {
      elt = elt->next;
      return (*this); 
    }
    bool operator==(const const_iterator & it) const {
      return (elt == it.elt); 
    }
    bool operator!=(const const_iterator & it) const {
      return (elt != it.elt);
    }
    bool operator<(const const_iterator & it) const {
      return (elt < it.elt);
    }   
    bool operator>(const const_iterator & it) const {
      return (elt > it.elt);
    }   
    friend class tinylist;
    MEMORY_DEBUG(tinylist::const_iterator)
  };
 private:
  tinylist_elt *head;
 public:
  tinylist() : head(0) {};
  ~tinylist() {
    if (head) {
      tinylist_elt *tle = head;
      while (tle) {
	tinylist_elt *tle1 = tle->next;
	delete tle;
	tle = tle1;
      }
    }
  }
  iterator push_front(const T & t) {
    head = new tinylist_elt(t,head);
    return head;
  }
  iterator insert_after(const iterator & it, const T & t) {
    it.elt->next = new tinylist_elt(t,it.elt->next);
    return it.elt->next;
  }
  iterator erase_after(const iterator & it) {
    if (it.elt->next) {
      tinylist_elt *tle = it.elt->next;
      it.elt->next = tle->next;
      delete tle;
    } else {
      it.elt->next = 0;
    }
    return it.elt->next;
  }
  void remove(const T & t) {
    iterator lit=begin();
    iterator lit0=end();
    while (lit!=end()) {
      if ((*lit) == t) {
	if (lit0 == end()) {
	  pop();
	} else {
	  erase_after(lit0);
	}
      }
      lit0=lit;
      ++lit;
    }
  }	
  iterator pop() {
    if (head) {
      tinylist_elt *tle = head;
      head = tle->next;
      delete tle;
    }
    return head;
  }
  const_iterator begin() const {
    return head;
  }
  const_iterator end() const {
    return 0;
  }
  iterator begin() {
    return head;
  }
  iterator end() {
    return 0;
  }
  bool empty() {
    return (head==0);
  }
  bool single() {
    return !empty()&&!head->next;
  }
  MEMORY_DEBUG(tinylist<T>)
};

#endif
