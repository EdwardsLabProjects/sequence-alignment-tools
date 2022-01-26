
#ifndef _BITS_H_
#define _BITS_H_

#include <sys/types.h>
typedef u_int32_t uint32;
typedef u_int64_t uint64;
typedef int32_t   int32;
typedef int64_t   int64;

unsigned short bits(uint32);
unsigned short bits(uint64);
unsigned int   bits(uint64*,unsigned short);

void zero(uint32 &);
void zero(uint64 &);
void zero(uint64 *, unsigned short);
  
// set/get
void set(uint32 &, unsigned char);
void set(uint64 &, unsigned char);
void set(uint64 *, unsigned short, unsigned int);
bool get(uint32, unsigned char);
bool get(uint64, unsigned char);
bool get(uint64 *, unsigned short, unsigned int); 

char *tostr(uint32);
char *tostr(uint64);
char *tostr(uint64*, unsigned short);

uint32 tobv32(const char *);
uint64 tobv64(const char *);

// Number of leading zeros
unsigned short nlz(uint32);
unsigned short nlz(uint64);
unsigned int   nlz(uint64*, unsigned short);

// log_2 round down
unsigned short flg2(uint32);
unsigned short flg2(uint64);
unsigned int   flg2(uint64*, unsigned short);

// log_2 round up
unsigned short clg2(uint32);
unsigned short clg2(uint64);
unsigned int   clg2(uint64*, unsigned short);

// population (# of ones)
unsigned short pop(uint32);
unsigned short pop(uint64);
unsigned int   pop(uint64 *, unsigned short);

// reverse complement of 2-bit DNA hash
// complements of Brian Walenz
uint32 rc(uint32 hash_value, unsigned short hash_bits);
uint64 rc(uint64 hash_value, unsigned short hash_bits);
void rc(uint64 *h1, uint64 * hash_value, unsigned short n, unsigned int hash_bits);

uint32 reverse(uint32 hash_value, unsigned short hash_bits);
uint64 reverse(uint64 hash_value, unsigned short hash_bits);

void period(uint32 x, unsigned short & p, unsigned short & pw);
bool periodIs(unsigned short p, uint32 x);
unsigned short periodWeight(unsigned short p, uint32 x);

void period(uint64 x, unsigned short & p, unsigned short & pw);
bool periodIs(unsigned short p, uint64 x);
unsigned short periodWeight(unsigned short p, uint64 x);

void period(uint64 *x, unsigned short n, unsigned short & p, unsigned short & pw);
bool periodIs(unsigned short p, uint64 *x, unsigned short n);
unsigned int periodWeight(unsigned short p, uint64 *x, unsigned short n);

unsigned short nshift(uint32);
unsigned short nshift(uint64);

unsigned short runs(uint32, unsigned short *);
unsigned short runs(uint64, unsigned short *);

uint32 mask32(unsigned short, unsigned short);
uint64 mask64(unsigned short, unsigned short);

#endif
