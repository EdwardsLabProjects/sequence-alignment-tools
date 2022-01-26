
#include "bits.h"
#include "util.h"
#include <assert.h>

unsigned short bits(uint32 x) {
  return 32;
}

unsigned short bits(uint64 x) {
  return 64;
}

unsigned int bits(uint64 *x, unsigned short n) {
  return (64*n);
}

void zero(uint32 & x) {
  x = 0;
}

void zero(uint64 & x) {
  x = 0;
}

void zero(uint64 *x, unsigned short n) {
  for (int i=0;i<n;i++) {
    x[i] = 0;
  }
}

void set(uint32 & x, unsigned char p) {
  assert(p<bits(x));
  x |= ((1UL)<<p);
}

void set(uint64 & x, unsigned char p) {
  assert(p<bits(x));
  x |= ((1UL)<<p);
}

void set(uint64 *x, unsigned short n, unsigned int p) {
  assert(p<bits(x,n));
  x[p/64] |= ((1UL)<<(p%64));
}

bool get(uint32 x, unsigned char p) {
  assert(p<bits(x));
  return (x & ((1UL)<<p));
}

bool get(uint64 x, unsigned char p) {
  assert(p<bits(x));
  return (x & ((1UL)<<p));
}

bool get(uint64 *x, unsigned short n, unsigned int p) {
  assert(p<bits(x,n));
  return (x[p/64] & ((1UL)<<(p%64)));
}

char *
tostr(uint32 x) {
  char *s=new char[33];
  for (int i=31;i>=0;i--) {
    s[i] = (get(x,i)?'1':'0');
  }
  s[32] = 0;
  return s;
}

char *
tostr(uint64 x) {
  char *s=new char[65];
  for (int i=63;i>=0;i--) {
    s[i] = (get(x,i)?'1':'0');
  }
  s[64] = 0;
  return s;
}

char *tostr(uint64* x, unsigned short n) {
  char *s = new char[64*n+n+1];
  int j=0;
  for (int k=n-1;k>=0;k--) {
    for (int i=63;i>=0;i--) {
      s[j++] = (get(x[k],i)?'1':'0');
    }
    if (k != 0) {
      s[j++] = ',';
    }
  }
  s[j++] = 0;
  return s;
}

uint32
tobv32(const char * tmplstr) {
  uint32 templ;
  assert(strlen(tmplstr)<=32);
  zero(templ);
  for (int i=0;i<strlen(tmplstr);i++) {
    assert(tmplstr[i] == '1' || tmplstr[i] == '0');
    if (tmplstr[i] != '0') {
      set(templ,i);
    }
  }
  return templ;
}

uint64
tobv64(const char * tmplstr) {
  uint64 templ;
  assert(strlen(tmplstr)<=64);
  zero(templ);
  for (int i=0;i<strlen(tmplstr);i++) {
    assert(tmplstr[i] == '1' || tmplstr[i] == '0');
    if (tmplstr[i] != '0') {
      set(templ,i);
    }
  }
  return templ;
}


unsigned short
nlz(uint32 x) {
  uint32 y;
  unsigned short n;
  
  n=32;
  y = x >>16; if (y != 0) { n -= 16; x = y; }
  y = x >> 8; if (y != 0) { n -=  8; x = y; }
  y = x >> 4; if (y != 0) { n -=  4; x = y; }
  y = x >> 2; if (y != 0) { n -=  2; x = y; }
  y = x >> 1; if (y != 0) return n - 2;
  return n - x;
}

unsigned short
nlz(uint64 x) {
  uint64 y;
  unsigned short n;
  
  n=64;
  y = x >>32; if (y != 0) { n -= 32; x = y; }
  y = x >>16; if (y != 0) { n -= 16; x = y; }
  y = x >> 8; if (y != 0) { n -=  8; x = y; }
  y = x >> 4; if (y != 0) { n -=  4; x = y; }
  y = x >> 2; if (y != 0) { n -=  2; x = y; }
  y = x >> 1; if (y != 0) return n - 2;
  return n - x;
}

unsigned int
nlz(uint64* x, unsigned short n) {
  unsigned int v=0;
  for (int i=n-1;i>=0;i--) {
    if (x[i] == 0) {
      v += bits(x[i]);
    } else {
      v += nlz(x[i]);
      return v;
    }
  }
  return v;
}

unsigned short
hbi(uint32 x) {
  return bits(x)-nlz(x);
}

unsigned short
hbi(uint64 x) {
  return bits(x)-nlz(x);
}

unsigned int
hbi(uint64 *x, unsigned short n) {
  return bits(x,n)-nlz(x,n);
}

unsigned short
pop(uint32 x) {
  x -=  (x>>1) & 0x55555555UL;
  x  = ((x>>2) & 0x33333333UL) + (x & 0x33333333UL);
  x  = ((x>>4) + x) & 0x0f0f0f0fUL;
  x *= 0x01010101UL;
  return  x>>24;
}

#define TWO(c)     (0x1ul << (c))
#define MASK(c)    (((unsigned long int)(-1)) / (TWO(TWO(c)) + 1u))
#define COUNT(x,c) ((x) & MASK(c)) + (((x) >> (TWO(c))) & MASK(c))
     
unsigned short
pop(uint64 x) {
  x = COUNT(x, 0) ;
  x = COUNT(x, 1) ;
  x = COUNT(x, 2) ;
  x = COUNT(x, 3) ;
  x = COUNT(x, 4) ;
  x = COUNT(x, 5) ;
  return  x;
}

unsigned int
pop(uint64 *x, unsigned short n) {
  unsigned int v=0;
  for (int i=0;i<n;i++) {
    v += pop((uint64)x[i]);
  }
  return v;
}

unsigned short
flg2(uint32 x) {
  return bits(x) - 1 - nlz(x);
}

unsigned short
flg2(uint64 x) {
  return bits(x) - 1 - nlz(x);
}

unsigned int
flg2(uint64 *x, unsigned short n) {
  return bits(x,n) - 1 - nlz(x,n);
}

void
minus1(uint64 *x1, uint64 *x, unsigned short n) {
  for (int i=0;i<n;i++) {
    x1[i] = x[i];
  }
  for (int i=0;i<n;i++) {
    if (x1[i] != 0) {
      x1[i]--;
      return;
    } else {
      x1[i] = (~(0UL));
    }
  }
}

void
shiftr(uint64 * x1, uint64 * x, unsigned short n, unsigned short s) {
  for (int i=0;i<n;i++) {
    x1[i] = 0;
  }
  for (int i=0;i<n;i++) {
    if (s/64+i < n) {
      x1[i] |= x[s/64+i]>>(s%64);
      if (s/64+1+i < n) {
	x1[i] |= x[s/64+1+i]<<(64-s%64);
      }
    }
  }  
}

void
shiftl(uint64 * x1, uint64 * x, unsigned short n, unsigned short s) {
  for (int i=0;i<n;i++) {
    x1[i] = 0;
  }
  // timestampi("  s: ",s);
  // timestamps("  x: ",tostr(x,n));
  for (int i=0;i<n;i++) {
    if (s/64+i < n) {
      x1[s/64+i] |= x[i]<<(s%64);
      // timestamps("ax1: ",tostr(x1,n));
      if (s/64+1+i < n) {
	x1[s/64+1+i] |= x[i]>>(64-s%64);
	// timestamps("bx1: ",tostr(x1,n));
      }
    }
  }  
  // timestamps("cx1: ",tostr(x1,n));
}

bool
equal(uint64 *x1, uint64 *x2, unsigned short n) {
  for (int i=0;i<n;i++) {
    if (x1[i] != x2[i]) {
      return false;
    }
  }
  return true;
}

unsigned short
clg2(uint32 x) {
  return bits(x) - nlz(x-1);
}

unsigned short
clg2(uint64 x) {
  return bits(x) - nlz(x-1);
}

unsigned int
clg2(uint64 *x, unsigned short n) {
  uint64 *x1 = new uint64[n];
  minus1(x1,x,n);
  unsigned int v = bits(x,n) - nlz(x1,n);
  delete [] x1;
  return v;
}

uint64 rc(uint64 x, unsigned short w) {
  // reverse
  uint64 y = x;
  y = ((y >>  2) & 0x3333333333333333llu) | ((y <<  2) & 0xccccccccccccccccllu);
  y = ((y >>  4) & 0x0f0f0f0f0f0f0f0fllu) | ((y <<  4) & 0xf0f0f0f0f0f0f0f0llu);
  y = ((y >>  8) & 0x00ff00ff00ff00ffllu) | ((y <<  8) & 0xff00ff00ff00ff00llu);
  y = ((y >> 16) & 0x0000ffff0000ffffllu) | ((y << 16) & 0xffff0000ffff0000llu);
  y = ((y >> 32) & 0x00000000ffffffffllu) | ((y << 32) & 0xffffffff00000000llu);
  // complement
  y ^= 0xffffffffffffffffllu;
  //  Shift and mask out the bases not in the mer
  y >>= 64 - w*2;
  return y;
}

uint32 rc(uint32 x, unsigned short w) {
  // reverse
  uint32 y = x;
  y = ((y >>  2) & 0x33333333lu) | ((y <<  2) & 0xcccccccclu);
  y = ((y >>  4) & 0x0f0f0f0flu) | ((y <<  4) & 0xf0f0f0f0lu);
  y = ((y >>  8) & 0x00ff00fflu) | ((y <<  8) & 0xff00ff00lu);
  y = ((y >> 16) & 0x0000fffflu) | ((y << 16) & 0xffff0000lu);
  // complement
  y ^= 0xfffffffflu;
  //  Shift and mask out the bases not in the mer
  y >>= 32 - w*2;
  return y;
}

uint64 reverse(uint64 x, unsigned short w) {
  // reverse
  uint64 y = x;
  y = ((y >>  1) & 0x5555555555555555llu) | ((y <<  1) & 0xaaaaaaaaaaaaaaaallu);
  y = ((y >>  2) & 0x3333333333333333llu) | ((y <<  2) & 0xccccccccccccccccllu);
  y = ((y >>  4) & 0x0f0f0f0f0f0f0f0fllu) | ((y <<  4) & 0xf0f0f0f0f0f0f0f0llu);
  y = ((y >>  8) & 0x00ff00ff00ff00ffllu) | ((y <<  8) & 0xff00ff00ff00ff00llu);
  y = ((y >> 16) & 0x0000ffff0000ffffllu) | ((y << 16) & 0xffff0000ffff0000llu);
  y = ((y >> 32) & 0x00000000ffffffffllu) | ((y << 32) & 0xffffffff00000000llu);
  //  Shift and mask out the bases not in the mer
  y >>= 64 - w;
  return y;
}

uint32 reverse(uint32 x, unsigned short w) {
  // reverse
  uint32 y = x;
  y = ((y >>  1) & 0x55555555lu) | ((y <<  1) & 0xaaaaaaaalu);
  y = ((y >>  2) & 0x33333333lu) | ((y <<  2) & 0xcccccccclu);
  y = ((y >>  4) & 0x0f0f0f0flu) | ((y <<  4) & 0xf0f0f0f0lu);
  y = ((y >>  8) & 0x00ff00fflu) | ((y <<  8) & 0xff00ff00lu);
  y = ((y >> 16) & 0x0000fffflu) | ((y << 16) & 0xffff0000lu);
  //  Shift and mask out the bases not in the mer
  y >>= 32 - w;
  return y;
}

// Support for figuring out the "period" of a spaced seed, and its periodic weight.

void 
period(uint32 x, unsigned short & p, unsigned short & pw) {
  int hbi = ::hbi(x);
  pw = bits(x);
  p = 0;
  for (int pi=1;pi<hbi;pi++) {
    if (periodIs(pi,x)) {
      int uc = periodWeight(pi,x);
      if (uc < pw) {
        pw = uc;
        p = pi;
      }
    }
  }
}

void 
period(uint64 x, unsigned short & p, unsigned short & pw) {
  int hbi = ::hbi(x);
  pw = bits(x);
  p = 0;
  for (int pi=1;pi<hbi;pi++) {
    if (periodIs(pi,x)) {
      int uc = periodWeight(pi,x);
      if (uc < pw) {
        pw = uc;
        p = pi;
      }
    }
  }
}

void 
period(uint64 *x, unsigned short n, unsigned short & p, unsigned short & pw) {
  int hbi = ::hbi(x,n);
  pw = bits(x,n);
  p = 0;
  for (int pi=1;pi<hbi;pi++) {
    if (periodIs(pi,x,n)) {
      int uc = periodWeight(pi,x,n);
      if (uc < pw) {
        pw = uc;
        p = pi;
      }
    }
  }
}

bool
periodIs(unsigned short p, uint32 bv) {
  uint32 sbv = (bv >> p);
  int hbi = ::hbi(bv);
  int lshift = bits(bv) - hbi + p;
  if (lshift >= (int)bits(bv)) {
    return true;
  }
  return ( (bv << lshift) == (sbv << lshift) );
}

bool
periodIs(unsigned short p, uint64 bv) {
  uint64 sbv = (bv >> p);
  int hbi = ::hbi(bv);
  int lshift = bits(bv) - hbi + p;
  if (lshift >= (int)bits(bv)) {
    return true;
  }
  return ( (bv << lshift) == (sbv << lshift) );
}

bool
periodIs(unsigned short p, uint64 * bv, unsigned short n) {
  // timestampi("n: ",n);
  // timestampi("Checking for period: ",p);
  uint64 *sbv = new uint64[n];
  shiftr(sbv,bv,n,p);
  // timestamps(" bv: ",tostr(bv,n));
  // timestamps("sbv: ",tostr(sbv,n));
  int hbi = ::hbi(bv,n);
  // timestampi("hbi: ",hbi);
  int lshift = bits(bv,n) - hbi + p;
  // timestampi("lsh: ",lshift);
  if (lshift >= (int)bits(bv,n)) {
    return true;
  }
  uint64 *x1 = new uint64[n];
  shiftl(x1,bv,n,lshift);
  uint64 *x2 = new uint64[n];
  shiftl(x2,sbv,n,lshift);
  bool v = equal(x1,x2,n);
  delete [] sbv;
  delete [] x1;
  delete [] x2;
  // timestampi("Returning ",v);
  return v;
}

unsigned short
periodWeight(unsigned short p, uint32 bv) {
  int hbi = ::hbi(bv);
  uint32 sbv = (bv >> (hbi-p));
  return pop(sbv);
}

unsigned short
periodWeight(unsigned short p, uint64 bv) {
  int hbi = ::hbi(bv);
  uint64 sbv = (bv >> (hbi-p));
  return pop(sbv);
}

unsigned int
periodWeight(unsigned short p, uint64 * x, unsigned short n) {
  int hbi = ::hbi(x,n);
  uint64 *x1 = new uint64[n];
  shiftr(x1,x,n,hbi-p);
  unsigned int v = pop(x1,n);
  delete [] x1;
  return v;
}

unsigned short
runs(uint32 bv, unsigned short *r) {
  int p = 0;
  for (int i=1;i<32;i++) {
    if ((!get(bv,i-1) && get(bv,i)) ||
	(get(bv,i-1) && !get(bv,i))) {
      r[p] = i;
      p++;
    }
  }
  if (get(bv,r[p-1])) {
    r[p] = 32;
    p++;
  }
  for (int i=p-1;i>0;i--) {
    r[i] -= r[i-1];
  }
  return p;
}

unsigned short
runs(uint64 bv, unsigned short *r) {
  int p = 0;
  for (int i=1;i<64;i++) {
    if ((!get(bv,i-1) && get(bv,i)) ||
	(get(bv,i-1) && !get(bv,i))) {
      r[p] = i;
      p++;
    }
  }
  if (get(bv,r[p-1])) {
    r[p] = 64;
    p++;
  }
  for (int i=p-1;i>0;i--) {
    r[i] -= r[i-1];
  }
  return p;
}

unsigned short
nshift(uint32 bv) {
  int p=0;
  for (int i=1;i<32;i++) {
    if (!get(bv,i-1) && get(bv,i)) {
      p++;
    }
  }
  return p;
}

unsigned short
nshift(uint64 bv) {
  int p=0;
  for (int i=1;i<64;i++) {
    if (!get(bv,i-1) && get(bv,i)) {
      p++;
    }
  }
  return p;
}

uint32 
mask32(unsigned short st, unsigned short ed) {
  return (((((uint32)1)<<ed)-1) - ((((uint32)1)<<st)-1));
}

uint64
mask64(unsigned short st, unsigned short ed) {
  return (((((uint64)1)<<ed)-1) - ((((uint64)1)<<st)-1));
}
