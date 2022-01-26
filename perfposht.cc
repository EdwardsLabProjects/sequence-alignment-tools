
#include "perfposht.h"

perfposht::perfposht()
{
  ;
}

perfposht::perfposht(hash & h, bool rc, bool cannon) 
{
  init(h,rc,cannon);
}

perfposht::perfposht(hash & h, bool rc, bool cannon, FILE_POSITION_TYPE dboffset, int mersize, bitmap & bm) {
  init(h,rc,cannon,mersize,dboffset,&bm);
}

void
perfposht::init(hash & h, bool rc, bool cannon, FILE_POSITION_TYPE dboffset, int mersize, bitmap *ignore) {
  // uint32 window=mersize-h.span()+1;
  h.reset();
  ptr_.clear();
  timestampi("Size of hash index: ",sizeof(uint32)*(h.size()+1));
  ptr_.resize(h.size()+1,0);
  int total=0;
  while (h.next()) {
//     cerr << h.pos() << "\t" 
// 	 << h.value() << "\t" 
// 	 << h.rcvalue() << "\t"
// 	 << h.cvalue() << endl;
    // if (h.pos() % 10000 == 0) {
    // timestamplu("Hash pos: ",h.pos());
    // timestamplu("Hash value: ",h.value());
    // timestamplu("Hash rcvalue: ",h.rcvalue());
    // timestamplu("Hash cvalue: ",h.cvalue());
    // }
    if (mersize && ignore->all(h.pos()-dboffset,mersize-h.span()+1)) {
      continue;
    }
    if (!rc) {
//       if (h.value() == 43080) { 
// 	 checkpoint;
// 	 cerr << h.value() << '\t' << h.pos() << '\t' << h.str(h.value()) << endl;
//       }
      ptr_[h.value()]++;
      total++;
    } else {
      if (!cannon) {
	ptr_[h.value()]++;
	ptr_[h.rcvalue()]++;
	total += 2;
      } else {
	hash_t v0 = h.value();
	hash_t v1 = h.rcvalue();
	if (v0 < v1) {
	  ptr_[v0]++;
	  total++;
	} else if (v0 > v1) {
	  ptr_[v1]++;
	  total++;
	} else {
	  /* v0 == v1, rare but does happen! */
	  ptr_[v0]++;
	  ptr_[v1]++;
	  total+=2;
	}
      }
    }
  }
  for (int i=1;i<ptr_.size();i++) {
    ptr_[i] += ptr_[i-1];
  }
  for (int i=ptr_.size()-1;i>0;i--) {
    ptr_[i] = ptr_[i-1];
  }
  ptr_[0] = 0;
  assert(total == ptr_[ptr_.size()-1]);
  pos_.clear();
  timestampi("Size of hash positions: ",sizeof(int32)*(total));
  pos_.resize(total,0);
  h.reset();
  while (h.next()) {
    if (mersize && ignore->all(h.pos()-dboffset,mersize-h.span()+1)) {
      continue;
    }
    hash_t v = h.value();
    if (!rc) {
//       if (h.value() == 43080) { 
// 	checkpoint;
// 	cerr << h.value() << '\t' << h.pos() << '\t' << h.str(h.value()) << endl;
//       }
      pos_[ptr_[v]] = h.pos()-dboffset;
      ptr_[v]++;
    } else {
      hash_t v1 = h.rcvalue();
      if (!cannon) {
	pos_[ptr_[v]] = h.pos()-dboffset;
	ptr_[v]++;
	pos_[ptr_[v1]] = -(h.pos()-dboffset);
	ptr_[v1]++;
      } else {
	hash_t v0 = h.value();
	hash_t v1 = h.rcvalue();
	if (v0 < v1) {
	  pos_[ptr_[v]] = h.pos()-dboffset;
	  ptr_[v]++;
	} else if (v0 > v1) {
	  pos_[ptr_[v1]] = -(h.pos()-dboffset);
	  ptr_[v1]++;
	} else {
	  /* v1 == v, rare but does happen! */
	  pos_[ptr_[v]] = h.pos()-dboffset;
	  ptr_[v]++;
	  pos_[ptr_[v1]] = -(h.pos()-dboffset);
	  ptr_[v1]++;
	}
      }
    }
  }
  for (int i=ptr_.size()-1;i>0;i--) {
    ptr_[i] = ptr_[i-1];
  }
  ptr_[0] = 0;

//   for (int i=0;i<ptr_.size();i++) {
//     for (int j=ptr_[i];j<ptr_[i+1];j++) {
//       if (j == ptr_[i]) {
// 	cerr << i << ":";
//       }
//       cerr << " " << pos_[j];
//       if (j == (ptr_[i+1]-1)) {
// 	cerr << endl;
//       }
//     }
//   }
}

perfposht::~perfposht() {
  ;
}

void
perfposht::set_invalid(hash_t v) {
  if (ptr_[v]<ptr_[v+1]) {
    pos_[ptr_[v]] = -1;
  }
}

bool
perfposht::is_valid(hash_t v) {
  return ((ptr_[v]<ptr_[v+1]) && (pos_[ptr_[v]] != -1));
}

bool
perfposht::lookup(hash_t v) {
  if (!is_valid(v)) {
    return false;
  }
  ls_ = ptr_[v]-1;
  le_ = ptr_[v+1];
  return true;
}

bool 
perfposht::next() {
  do {
    ls_++;
  } while (ls_<le_ && pos_[ls_] == 0);
  return (ls_<le_);
}

int32
perfposht::value() {
  return pos_[ls_];
}
  
void
perfposht::set_invalid() {
  pos_[ls_] = 0;
}

/* 
void
perfposht::check(CharacterProducer & cp) {
  for (hash_t i=0;i<h_.size();i++) {
    lookup(i);
    while (next()) {
      int32 p = value();
      if (p > 0) {
	cp.pos(p-h_.span());
	string s = cp.getstr((unsigned int)h_.span());
	string s1 = h_.str(i);
	cerr << i << " " << p << " " << s << " " << s1 << " " << ((s==s1)?"TRUE":"FALSE") << endl;
      } else {
	cp.pos(-p-h_.span());
	string s = cp.getstr((unsigned int)h_.span());
	string s1 = reverse_comp(h_.str(i));
	cerr << i << " " << p << " " << s << " " << s1 << " " << ((s==s1)?"TRUE":"FALSE") << endl;
      }
    }
  }
}
*/
