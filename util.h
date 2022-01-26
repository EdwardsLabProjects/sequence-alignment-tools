
#ifndef _IBPEP_UTIL_H_
#define _IBPEP_UTIL_H_

#include "types.h"

unsigned int least_common_multiple(unsigned int a, unsigned int b);
bool is_true(char const * const s);
bool is_false(char const * const s);

#include <string>

void uppercase(std::string & s);
std::string reverse_comp(std::string const & sequence);
std::string reverse(std::string const & sequence);

char iupac_revcomp(char);
bool iupac_contains(char,char);
bool iupac_contained(char,char);
bool iupac_compatible(char,char);
char *iupac_contains(char);
char *iupac_contained(char);
char *iupac_compatible(char);

char charmap(int,char);

char trans_codon(int,char[],char*);

double monomolwt(char);
double avemolwt(char);
int aasubdist(char,char);
int aacodonsubdist(char,char,char);

#define bigword long unsigned int

std::string binarystr(bigword val, int length);
std::string binary(long unsigned int val);
std::string binary(long long unsigned int val);
std::string binary(unsigned char val);

bool exist(std::string const & filename);
time_t modtime(std::string const & filename);
FILE_POSITION_TYPE filesize(std::string const & filename);
std::string change_extn(std::string const & filename, std::string extn);

int anypos(std::string const & s0, std::string const & s1, int pos=0);

#ifdef PROFILE
void set_profile_signal_handler();
#endif

#include <stdio.h>
#include <time.h>
#include <string.h>
#define timestamp(m) { time_t t(time(NULL)); ASCTIME(s,localtime(&t)); fprintf(stderr,"[%24.*s] %s\n",(int)strlen(s)-1,s,(m));}
#define timestampli(m,v) { time_t t(time(NULL)); ASCTIME(s,localtime(&t)); fprintf(stderr,"[%24.*s] %s%ld\n",(int)strlen(s)-1,s,(m),(v));}
#define timestamplu(m,v) { time_t t(time(NULL)); ASCTIME(s,localtime(&t)); fprintf(stderr,"[%24.*s] %s%lu\n",(int)strlen(s)-1,s,(m),(v));}
#define timestampi(m,v) { time_t t(time(NULL)); ASCTIME(_s_,localtime(&t)); fprintf(stderr,"[%24.*s] %s%d\n",(int)strlen(_s_)-1,_s_,(m),(v));}
#define timestampd(m,v) { time_t t(time(NULL)); ASCTIME(s,localtime(&t)); fprintf(stderr,"[%24.*s] %s%g\n",(int)strlen(s)-1,s,(m),(v));}
#define timestamps(m,v) { time_t t(time(NULL)); ASCTIME(_s_,localtime(&t)); fprintf(stderr,"[%24.*s] %s%s\n",(int)strlen(_s_)-1,_s_,(m),(v));}
#define checkpoint { time_t t(time(NULL)); ASCTIME(s,localtime(&t)); fprintf(stderr,"[%24.*s] Checkpoint %s:%d.\n",(int)strlen(s)-1,s,__FILE__,__LINE__); }

extern time_t tictoctime_;
#define tic {tictoctime_ = time(NULL);}
#define toc { time_t t; time(&t); t-=tictoctime_; fprintf(stderr,"Elapsed time: %d:%02d:%02d\n",t/3600,(t%3600)/60,t%60);} 
#define tocassgn(tsec) { time_t t; time(&t); t-=tictoctime_; fprintf(stderr,"Elapsed time: %ld:%02ld:%02ld\n",t/3600,(t%3600)/60,t%60); tsec = t;} 

#ifdef PROFILE
void instrument1();
void instrument2();
void instrument3();
void instrument4();
void instrument5();
void instrument6();
void instrument7();
void instrument8();
void instrument9();
void instrument10();
void instrument11();
#else
#define instrument1()
#define instrument2()
#define instrument3()
#define instrument4()
#define instrument5()
#define instrument6()
#define instrument7()
#define instrument8()
#define instrument9()
#define instrument10()
#define instrument11()
#endif

#endif

