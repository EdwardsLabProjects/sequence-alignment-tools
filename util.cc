
#include <ctype.h>
#include <string.h>
#include <signal.h>
#include <unistd.h>
#include <iostream>
#include "util.h"
#include "types.h"

#ifdef PROFILE 

void profile_signal_handler(int sig) {
  switch (sig) {
  case SIGUSR1:
    exit(0);
    break;
  default:
    ; // do nothing, should not happen!
  }
}

void set_profile_signal_handler() {
  struct sigaction sa;
  sa.sa_handler = &profile_signal_handler;
  // SIGINITSET(sa.sa_mask);
  sigemptyset(&sa.sa_mask);
  sa.sa_flags = 0;
  sigaction(SIGUSR1,&sa,0);
}

#ifdef PROFILE
void instrument1() {1;};
void instrument2() {1;};
void instrument3() {1;};
void instrument4() {1;};
void instrument5() {1;};
void instrument6() {1;};
void instrument7() {1;};
void instrument8() {1;};
void instrument9() {1;};
void instrument10() {1;};
void instrument11() {1;};
#endif

#endif

time_t tictoctime_=0;

unsigned int least_common_multiple(unsigned int a, unsigned int b) {
  unsigned int mult = a*b;
  while (b != 0) {
    unsigned int tmp;
    tmp = a; a = b; b = tmp;
    b = b % a;
  }
  return mult/a;
};

static const char *trues[] = {"true", "True", "TRUE",
			 "T", "t", "1",
			 "Yes", "YES", "yes", ""};
static const char *falses[] = {"false", "False", "FALSE",
			 "F", "f", "0",
			 "No", "NO", "no", ""};


bool find_str(char const * const s, const char * values[]) {
  const char **p = values;
  while (p && strlen(*p) > 0) {
    if (*p && !strcmp(*p,s)) {
      return true;
    }
    p++;
  }
  return false;
}

bool is_false(char const * const s) {
  return find_str(s,falses);
}

bool is_true(char const * const s) {
  return find_str(s,trues);
}

std::string binarystr(bigword val, int length) {
  std::string str="";
  bigword mask = 1;
  for (int i=0;i<length;i++) {    
    if (val & mask) {
      str = "1" + str;
    } else {
      str = "0" + str;
    }
    if (i%8 == 7) {
      str = " " + str;
    }
    mask <<= 1;
  }
  return str;
}

std::string binary(unsigned char val) {
  return binarystr(val,sizeof(unsigned char)*8);
}

std::string binary(long unsigned int val) {
  return binarystr(val,sizeof(long unsigned int)*8);
}

std::string binary(long long unsigned int val) {
  return binarystr(val,sizeof(long long unsigned int)*8);
}

void uppercase(std::string & s) {
  for (unsigned int i=0;i<s.length();i++) {
    s[i] = toupper(s[i]);
  }
}

char *iupac_compatible(char w) {
  static char **compatiblestr=0;
  if (!compatiblestr) {
    compatiblestr = new char*[128];
    memset(compatiblestr,0,128*sizeof(char*));
    compatiblestr['A'] = strdup("ARMWDHVN");        
    compatiblestr['B'] = strdup("GTUCYKSBN");
    compatiblestr['C'] = strdup("CYMSBHVN");        
    compatiblestr['D'] = strdup("GATURWKDN");
    compatiblestr['G'] = strdup("GRKSBDVN");        
    compatiblestr['H'] = strdup("ACTUMYWHN");
    compatiblestr['K'] = strdup("GTKBDN");          
    compatiblestr['M'] = strdup("ACMHVN");
    compatiblestr['N'] = strdup("ACGTURYKMSWVDHVN"); 
    compatiblestr['R'] = strdup("GARDVN");
    compatiblestr['S'] = strdup("GCSBVN");          
    compatiblestr['T'] = strdup("TUYKWVDHN");
    compatiblestr['U'] = strdup("UTYKWVDHN");       
    compatiblestr['V'] = strdup("GCARSMVN");
    compatiblestr['W'] = strdup("ATUWDHN");         
    compatiblestr['Y'] = strdup("TUCYBHN");
    compatiblestr['X'] = strdup("MRWSYKVHDBXN");
    compatiblestr['a'] = strdup("armwdhvn");        
    compatiblestr['b'] = strdup("gtucyksbn");
    compatiblestr['c'] = strdup("cymsbhvn");        
    compatiblestr['d'] = strdup("gaturwkdn");
    compatiblestr['g'] = strdup("grksbdvn");        
    compatiblestr['h'] = strdup("actumywhn");
    compatiblestr['k'] = strdup("gtkbdn");          
    compatiblestr['m'] = strdup("acmhvn");
    compatiblestr['n'] = strdup("acgturykmswvdhvn"); 
    compatiblestr['r'] = strdup("gardvn");
    compatiblestr['s'] = strdup("gcsbvn");          
    compatiblestr['t'] = strdup("tuykwvdhn");
    compatiblestr['u'] = strdup("utykwvdhn");       
    compatiblestr['v'] = strdup("gcarsmvn");
    compatiblestr['w'] = strdup("atuwdhn");         
    compatiblestr['x'] = strdup("mrwsykvhdbxn");
    compatiblestr['y'] = strdup("tucybhn");
  }
  return compatiblestr[w];
}

bool iupac_compatible(char w, char c) {
  static bool **compatiblemap=0;
  if (!compatiblemap) {
    bool * buffer = new bool[128*128];
    memset(buffer,0,128*128*sizeof(bool));
    compatiblemap = new bool*[128];
    memset(compatiblemap,0,128*sizeof(bool*));
    char *compatible=0;
    int j;
    for (int w=0;w<128;w++) {
      compatiblemap[w] = buffer+w*128;
      if ((compatible=iupac_compatible(w))!=0) {
	j=0;
	while (compatible[j]) {
	  compatiblemap[w][compatible[j]] = true;
	  j++;
	}
      }
    }
  }
  return compatiblemap[w][c];
}

char *iupac_contains(char w) {
  static char **containsstr=0;
  if (!containsstr) {
    containsstr = new char*[256];
    memset(containsstr,0,256*sizeof(char));
    containsstr['A'] = strdup("A");        
    containsstr['B'] = strdup("GTUCYKSB");
    containsstr['C'] = strdup("C");        
    containsstr['D'] = strdup("GATURWKD");
    containsstr['G'] = strdup("G");        
    containsstr['H'] = strdup("ACTUMYWH");
    containsstr['K'] = strdup("GTK");          
    containsstr['M'] = strdup("ACM");
    containsstr['N'] = strdup("ACGTURYKMSWVDHVN"); 
    containsstr['R'] = strdup("GAR");
    containsstr['S'] = strdup("GCS");          
    containsstr['T'] = strdup("TU");
    containsstr['U'] = strdup("UT");       
    containsstr['V'] = strdup("GCARSMV");
    containsstr['W'] = strdup("ATUW");         
    containsstr['Y'] = strdup("TUCY");
    containsstr['X'] = strdup("MRWSYKVHDBXN");
    containsstr['a'] = strdup("a");        
    containsstr['b'] = strdup("gtucyksb");
    containsstr['c'] = strdup("c");        
    containsstr['d'] = strdup("gaturwkd");
    containsstr['g'] = strdup("g");        
    containsstr['h'] = strdup("actumywh");
    containsstr['k'] = strdup("gtk");          
    containsstr['m'] = strdup("acm");
    containsstr['n'] = strdup("acgturykmswvdhvn"); 
    containsstr['r'] = strdup("gar");
    containsstr['s'] = strdup("gcs");          
    containsstr['t'] = strdup("tu");
    containsstr['u'] = strdup("ut");       
    containsstr['v'] = strdup("gcarsmv");
    containsstr['w'] = strdup("atuw");         
    containsstr['y'] = strdup("tucy");
    containsstr['x'] = strdup("mrwsykvhdbxn");
  }
  return containsstr[w];
}

bool iupac_contains(char w, char c) {
  static bool **containsmap=0;
  if (!containsmap) {
    bool * buffer = new bool[256*256];
    memset(buffer,0,256*256*sizeof(bool));
    containsmap = new bool*[256];
    memset(containsmap,0,256*sizeof(bool*));
    char *contained;
    int j;
    for (int w=0;w<256;w++) {
      containsmap[w] = buffer+w*256;
      if ((contained=iupac_contains(w))!=0) {
	j=0;
	while (contained[j]) {
	  containsmap[w][contained[j]] = true;
	  j++;
	}
      }
    }
  }
  return containsmap[w][c];
}

char *iupac_contained(char w) {
  static char **containedstr=0;
  if (!containedstr) {
    containedstr = new char*[256];
    memset(containedstr,0,256*sizeof(char));
    containedstr['A'] = strdup("ARMWDHVN");        
    containedstr['B'] = strdup("BNX");
    containedstr['C'] = strdup("CYMSBHVN");        
    containedstr['D'] = strdup("DNX");
    containedstr['G'] = strdup("GRKSBDVN");        
    containedstr['H'] = strdup("HNX");
    containedstr['K'] = strdup("KBDNX");          
    containedstr['M'] = strdup("MHVNX");
    containedstr['N'] = strdup("NX"); 
    containedstr['R'] = strdup("RDVNX");
    containedstr['S'] = strdup("SBVNX");          
    containedstr['T'] = strdup("TUYKWVDHN");
    containedstr['U'] = strdup("UTYKWVDHN");       
    containedstr['V'] = strdup("VNX");
    containedstr['W'] = strdup("WDHNX");         
    containedstr['Y'] = strdup("YBHNX");
    containedstr['X'] = strdup("X");
    containedstr['a'] = strdup("armwdhvn");        
    containedstr['b'] = strdup("bnx");
    containedstr['c'] = strdup("cymsbhvn");        
    containedstr['d'] = strdup("dnx");
    containedstr['g'] = strdup("grksbdvn");        
    containedstr['h'] = strdup("hnx");
    containedstr['k'] = strdup("kbdnx");          
    containedstr['m'] = strdup("mhvnx");
    containedstr['n'] = strdup("nx"); 
    containedstr['r'] = strdup("rdvnx");
    containedstr['s'] = strdup("sbvnx");          
    containedstr['t'] = strdup("tuykwvdhn");
    containedstr['u'] = strdup("utykwvdhn");       
    containedstr['v'] = strdup("vnx");
    containedstr['w'] = strdup("wdhnx");         
    containedstr['y'] = strdup("ybhnx");
    containedstr['x'] = strdup("x");
  }
  return containedstr[w];
}

bool iupac_contained(char w, char c) {
  static bool **containedmap=0;
  if (!containedmap) {
    bool * buffer = new bool[256*256];
    memset(buffer,0,256*256*sizeof(bool));
    containedmap = new bool*[256];
    memset(containedmap,0,256*sizeof(bool*));
    char *contains;
    int j;
    for (int w=0;w<256;w++) {
      containedmap[w] = buffer+w*256;
      if ((contains=iupac_contained(w))!=0) {
	j=0;
	while (contains[j]) {
	  containedmap[w][contains[j]] = true;
	  j++;
	}
      }
    }
  }
  return containedmap[w][c];
}

char iupac_revcomp(char c) {
  static char *rcmap = 0;
  if (!rcmap) {
    rcmap = new char[256];
    for (int i=0;i<256;i++) {
      rcmap[i] = i;
    }

    rcmap['a'] = 't';          rcmap['A'] = 'T';
    rcmap['c'] = 'g';	       rcmap['C'] = 'G';
    rcmap['g'] = 'c';	       rcmap['G'] = 'C';
    rcmap['t'] = 'a';	       rcmap['T'] = 'A';
    rcmap['u'] = 'a';	       rcmap['U'] = 'A';
			                        
    rcmap['m'] = 'k';	       rcmap['M'] = 'K';
    rcmap['r'] = 'y';	       rcmap['R'] = 'Y';
    rcmap['w'] = 'w';	       rcmap['W'] = 'W';
    rcmap['s'] = 's';	       rcmap['S'] = 'S';
    rcmap['y'] = 'r';	       rcmap['Y'] = 'R';
    rcmap['k'] = 'm';	       rcmap['K'] = 'M';
			                        
    rcmap['v'] = 'b';	       rcmap['V'] = 'B';
    rcmap['h'] = 'd';	       rcmap['H'] = 'D';
    rcmap['d'] = 'h';	       rcmap['D'] = 'H';
    rcmap['b'] = 'v';	       rcmap['B'] = 'V';
  }
  return rcmap[c];
}

char charmap(int mapindex, char c) {
  if (mapindex == 2) {
    if (c == 'i') {
      return 'l';
    } else if (c == 'I') {
      return 'L';
    } else {
      return c;
    }
  } else if (mapindex == 3) {
    if (c == 'i') {
      return 'l';
    } else if (c == 'I') {
      return 'L';
    } else if (c == 'k') {
      return 'q';
    } else if (c == 'K') {
      return 'Q';
    } else {
      return c;
    }
  } else {
    return c;
  }
}

std::string reverse_comp(std::string const & sequence) {
  std::string str(sequence);
  int i,n = sequence.length();
  for(i=0; i<n; i++) {
    str[i] = iupac_revcomp(sequence[n-1-i]);
  }
  return str;
}

std::string reverse(std::string const & sequence) {
  std::string str(sequence);
  int i,n = sequence.length();
  for(i=0; i<n; i++) {
    str[i] = sequence[n-1-i];
  }
  return str;
}

double monomolwt(char c) {
  static double *monomolwt = 0;
  c = toupper(c);
  if (!monomolwt) {
    monomolwt = new double[256];
    for (int i=0;i<256;i++) {
      monomolwt[i] = -1;
    }
    monomolwt['A']=71.037113848;
    monomolwt['C']=103.009185648; 
    monomolwt['D']=115.026943128; 
    monomolwt['E']=129.042593208; 
    monomolwt['F']=147.068414008; 
    monomolwt['G']=57.021463768; 
    monomolwt['H']=137.058911944; 
    monomolwt['I']=113.084064088; 
    monomolwt['K']=128.094963136; 
    monomolwt['L']=113.084064088; 
    monomolwt['M']=131.040485808; 
    monomolwt['N']=114.042927536; 
    monomolwt['P']=97.052763928; 
    monomolwt['Q']=128.058577616; 
    monomolwt['R']=156.101111152; 
    monomolwt['S']=87.032028488; 
    monomolwt['T']=101.047678568; 
    monomolwt['V']=99.068414008;
    monomolwt['W']=186.079313056; 
    monomolwt['Y']=163.063328648; 
  }
  return monomolwt[c];
}

double avemolwt(char c) {
  static double *avemolwt = 0;
  c = toupper(c);
  if (!avemolwt) {
    avemolwt = new double[256];
    for (int i=0;i<256;i++) {
      avemolwt[i] = -1;
    }
    avemolwt['A']=71.078826901;
    avemolwt['C']=103.143216117;
    avemolwt['D']=115.088513436;
    avemolwt['E']=129.115401675;
    avemolwt['F']=147.176750991;
    avemolwt['G']=57.051938663;
    avemolwt['H']=137.141315021;
    avemolwt['I']=113.159491617;
    avemolwt['K']=128.174180322;
    avemolwt['L']=113.159491617;
    avemolwt['M']=131.196992594;
    avemolwt['N']=114.103877326;
    avemolwt['P']=97.116752043;
    avemolwt['Q']=128.130765564;
    avemolwt['R']=156.187706397;
    avemolwt['S']=87.078151717;
    avemolwt['T']=101.105039956;
    avemolwt['V']=99.132603378;
    avemolwt['W']=186.213513503;
    avemolwt['Y']=163.176075807;
  }
  return avemolwt[c];
}

int aacodonsubdist(char f, char c, char t) {
  static signed char *aasubdist = 0;
  static signed char *aachars =0;
  f = toupper(f);
  t = toupper(t);
  if (!aasubdist) {
    aasubdist = new signed char[20*6*20];
    for (int i=0;i<2400;i++) {
      aasubdist[i] = -1;
    }
    aachars = new signed char[256];
    for (int i=0;i<256;i++) {
      aachars[i] = -1;
    }
    char aas[] = "ARNDCQEGHILKMFPSTWYV";
    for (int i=0;i<20;i++) {
      aachars[aas[i]] = i;
    }
    int subs[] = {
0,2,2,1,2,3,2,1,2,2,2,3,3,2,1,1,1,3,2,1,
0,2,2,1,2,3,2,1,2,2,2,3,3,2,1,1,1,3,2,1,
0,2,3,2,3,2,1,1,3,2,2,2,3,3,1,1,1,3,3,1,
0,2,3,2,3,2,1,1,3,3,2,2,2,3,1,1,1,2,3,1,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
2,0,2,2,1,2,3,1,1,2,1,3,3,2,1,1,2,2,2,2,
2,0,2,2,1,2,3,1,1,2,1,3,3,2,1,1,2,2,2,2,
2,0,3,3,2,1,2,1,2,2,1,2,3,3,1,2,2,2,3,2,
2,0,3,3,2,1,2,1,2,3,1,2,2,3,1,2,2,1,3,2,
2,0,2,3,2,2,2,1,3,1,2,1,2,3,2,1,1,2,3,2,
2,0,2,3,2,2,2,1,3,2,2,1,1,3,2,1,1,1,3,2,
2,2,0,1,2,2,2,2,1,1,2,1,2,2,2,1,1,3,1,2,
2,2,0,1,2,2,2,2,1,1,2,1,2,2,2,1,1,3,1,2,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
1,2,1,0,2,2,1,1,1,2,2,2,3,2,2,2,2,3,1,1,
1,2,1,0,2,2,1,1,1,2,2,2,3,2,2,2,2,3,1,1,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
2,1,2,2,0,3,3,1,2,2,2,3,3,1,2,1,2,1,1,2,
2,1,2,2,0,3,3,1,2,2,2,3,3,1,2,1,2,1,1,2,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
2,1,2,2,3,0,1,2,1,2,1,1,3,3,1,2,2,3,2,2,
2,1,2,2,3,0,1,2,1,3,1,1,2,3,1,2,2,2,2,2,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
1,2,2,1,3,1,0,1,2,2,2,1,3,3,2,2,2,3,2,1,
1,2,2,1,3,1,0,1,2,3,2,1,2,3,2,2,2,2,2,1,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
1,1,2,1,1,3,2,0,2,2,2,3,3,2,2,1,2,2,2,1,
1,1,2,1,1,3,2,0,2,2,2,3,3,2,2,1,2,2,2,1,
1,1,3,2,2,2,1,0,3,2,2,2,3,3,2,2,2,2,3,1,
1,1,3,2,2,2,1,0,3,3,2,2,2,3,2,2,2,1,3,1,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
2,1,1,1,2,1,2,2,0,2,1,2,3,2,1,2,2,3,1,2,
2,1,1,1,2,1,2,2,0,2,1,2,3,2,1,2,2,3,1,2,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
2,2,1,2,2,3,3,2,2,0,1,2,1,1,2,1,1,3,2,1,
2,2,1,2,2,3,3,2,2,0,1,2,1,1,2,1,1,3,2,1,
2,1,2,3,3,2,2,2,3,0,1,1,1,2,2,2,1,3,3,1,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
2,2,3,3,2,2,2,2,3,1,0,2,2,1,2,1,2,2,2,1,
2,2,3,3,2,2,2,2,3,2,0,2,1,1,2,1,2,1,2,1,
2,1,2,2,2,2,3,2,1,1,0,3,2,1,1,2,2,3,2,1,
2,1,2,2,2,2,3,2,1,1,0,3,2,1,1,2,2,3,2,1,
2,1,3,3,3,1,2,2,2,1,0,2,2,2,1,2,2,3,3,1,
2,1,3,3,3,1,2,2,2,2,0,2,1,2,1,2,2,2,3,1,
2,1,1,2,3,1,1,2,2,1,2,0,2,3,2,2,1,3,2,2,
2,1,1,2,3,1,1,2,2,2,2,0,1,3,2,2,1,2,2,2,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
2,1,2,3,3,2,2,2,3,1,1,1,0,2,2,2,1,2,3,1,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
2,2,2,2,1,3,3,2,2,1,1,3,2,0,2,1,2,2,1,1,
2,2,2,2,1,3,3,2,2,1,1,3,2,0,2,1,2,2,1,1,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
1,1,2,2,2,2,3,2,1,2,1,3,3,2,0,1,1,3,2,2,
1,1,2,2,2,2,3,2,1,2,1,3,3,2,0,1,1,3,2,2,
1,1,3,3,3,1,2,2,2,2,1,2,3,3,0,1,1,3,3,2,
1,1,3,3,3,1,2,2,2,3,1,2,2,3,0,1,1,2,3,2,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
1,2,2,2,1,3,3,2,2,2,2,3,3,1,1,0,1,2,1,2,
1,2,2,2,1,3,3,2,2,2,2,3,3,1,1,0,1,2,1,2,
1,2,3,3,2,2,2,2,3,2,1,2,3,2,1,0,1,2,2,2,
1,2,3,3,2,2,2,2,3,3,1,2,2,2,1,0,1,1,2,2,
2,1,1,2,1,3,3,1,2,1,2,2,2,2,2,0,1,2,2,2,
2,1,1,2,1,3,3,1,2,1,2,2,2,2,2,0,1,2,2,2,
1,2,1,2,2,3,3,2,2,1,2,2,2,2,1,1,0,3,2,2,
1,2,1,2,2,3,3,2,2,1,2,2,2,2,1,1,0,3,2,2,
1,1,2,3,3,2,2,2,3,1,2,1,2,3,1,1,0,3,3,2,
1,1,2,3,3,2,2,2,3,2,2,1,1,3,1,1,0,2,3,2,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
2,1,3,3,1,2,2,1,3,3,1,2,2,2,2,1,2,0,2,2,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
2,2,1,1,1,2,2,2,1,2,2,2,3,1,2,1,2,2,0,2,
2,2,1,1,1,2,2,2,1,2,2,2,3,1,2,1,2,2,0,2,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
1,2,2,1,2,3,2,1,2,1,1,3,2,1,2,2,2,3,2,0,
1,2,2,1,2,3,2,1,2,1,1,3,2,1,2,2,2,3,2,0,
1,2,3,2,3,2,1,1,3,1,1,2,2,2,2,2,2,3,3,0,
1,2,3,2,3,2,1,1,3,2,1,2,1,2,2,2,2,2,3,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int k=0;
    for (int i=0;i<20;i++) {
     for (int l=0;l<6;l++) {
      for (int j=0;j<20;j++) {
	aasubdist[120*aachars[aas[i]]+20*l+aachars[aas[j]]] = subs[k];
	k++;
      }
     }
    }
  }
  if (aachars[f] < 0 || aachars[t] < 0) {
    return -1;
  } else {
    return aasubdist[120*aachars[f]+20*(c-'0')+aachars[t]];
  }
}

int aasubdist(char f, char t) {
  static signed char *aasubdist = 0;
  static signed char *aachars =0;
  f = toupper(f);
  t = toupper(t);
  if (!aasubdist) {
    aasubdist = new signed char[400];
    for (int i=0;i<400;i++) {
      aasubdist[i] = -1;
    }
    aachars = new signed char[256];
    for (int i=0;i<256;i++) {
      aachars[i] = -1;
    }
    char aas[] = "ARNDCQEGHILKMFPSTWYV";
    for (int i=0;i<20;i++) {
      aachars[aas[i]] = i;
    }
    int subs[] = {0,
		  2,0,
                  2,2,0,
                  1,2,1,0,
                  2,1,2,2,0,
                  2,1,2,2,3,0,
                  1,2,2,1,3,1,0,
                  1,1,2,1,1,2,1,0,
                  2,1,1,1,2,1,2,2,0,
                  2,1,1,2,2,2,2,2,2,0,
                  2,1,2,2,2,1,2,2,1,1,0,
                  2,1,1,2,3,1,1,2,2,1,2,0,
                  2,1,2,3,3,2,2,2,3,1,1,1,0,
                  2,2,2,2,1,3,3,2,2,1,1,3,2,0,
                  1,1,2,2,2,1,2,2,1,2,1,2,2,2,0,
                  1,1,1,2,1,2,2,1,2,1,1,2,2,1,1,0,
                  1,1,1,2,2,2,2,2,2,1,2,1,1,2,1,1,0,
                  2,1,3,3,1,2,2,1,3,3,1,2,2,2,2,1,2,0,
                  2,2,1,1,1,2,2,2,1,2,2,2,3,1,2,1,2,2,0,
                  1,2,2,1,2,2,1,1,2,1,1,2,1,1,2,2,2,2,2,0};
    int k=0;
    for (int i=0;i<20;i++) {
      for (int j=0;j<=i;j++) {
	aasubdist[20*aachars[aas[i]]+aachars[aas[j]]] = subs[k];
	aasubdist[20*aachars[aas[j]]+aachars[aas[i]]] = subs[k];
	k++;
      }
    }
  }
  if (aachars[f] < 0 || aachars[t] < 0) {
    return -1;
  } else {
    return aasubdist[20*aachars[f]+aachars[t]];
  }
}







#include <sys/stat.h>
// #include <errno.h>

bool exist(std::string const & filename)
{ 
  struct stat stat_buf;
  if (stat(filename.c_str(),&stat_buf) != 0) {
    // std::cerr << "Errno for stat of " << filename << " is " << errno << "(" << sys_errlist[errno] << ")" << std::endl;
    return false;
  }
  return (stat_buf.st_mode & S_IFMT) == S_IFREG;
}

time_t modtime(std::string const & filename) {
  struct stat stat_buf;
  if (stat(filename.c_str(),&stat_buf) != 0) return ((time_t)0);
  if (stat_buf.st_size == 0) return ((time_t)0);
  time_t lastmod(stat_buf.st_mtime);
  if (lastmod < stat_buf.st_ctime) {
    lastmod = stat_buf.st_ctime;
  }
  return lastmod;
}

FILE_POSITION_TYPE filesize(std::string const & filename) 
{
  struct stat st;
  if (stat(filename.c_str(),&st)<0) {
    return -1;
  }
  return st.st_size;
}

int anypos(std::string const & s0, std::string const & s1, int pos) {
  std::string::size_type retval;
  if ((retval=s0.find_first_of(s1,pos))==std::string::npos) {
    return -1;
  } else {
    return retval;
  }
} 

char trans_codon(int frame, char codon[3], char *codonid=NULL) {
  static char *codontable = 0;
  static char *maptable = 0;
  static char *rcmaptable = 0;
  static char *codonidtable = 0;
  if (!codontable) {
    unsigned int i;
    maptable = new char[256];
    rcmaptable = new char[256];
    for (i=0;i<256;i++) {
      rcmaptable[i] = maptable[i] = 4;
    }
    rcmaptable['T'] = maptable['A'] = 0;
    rcmaptable['G'] = maptable['C'] = 1;
    rcmaptable['C'] = maptable['G'] = 2;
    rcmaptable['A'] = maptable['T'] = 3;

    const char *aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    const char *b1 = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
    const char *b2 = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
    const char *b3 = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";
    
    char *aacodonid = new char[256];
    for (int i=0;i<256;i++) {
      aacodonid[i] = '0';
    }
    codontable = new char[5*5*5];
    codonidtable = new char[5*5*5];
    for (i=0;i<5*5*5;i++) {
      codontable[i] = 'X';
      codonidtable[i] = 'X';
    }
    for (i=0;i<strlen(aa);i++) {
      codontable[25*maptable[b1[i]]+5*maptable[b2[i]]+maptable[b3[i]]] = aa[i];
      codonidtable[25*maptable[b1[i]]+5*maptable[b2[i]]+maptable[b3[i]]] = aacodonid[aa[i]];
      aacodonid[aa[i]]++;
    }
  }
  if (frame < 3) {
    if (codonid) {
      (*codonid) = codonidtable[25*maptable[codon[0]]+5*maptable[codon[1]]+maptable[codon[2]]];
    }
    return codontable[25*maptable[codon[0]]+5*maptable[codon[1]]+maptable[codon[2]]];
  } else {
    if (codonid) {
      (*codonid) = codonidtable[25*rcmaptable[codon[2]]+5*rcmaptable[codon[1]]+rcmaptable[codon[0]]];
    }
    return codontable[25*rcmaptable[codon[2]]+5*rcmaptable[codon[1]]+rcmaptable[codon[0]]];
  }
}


std::string change_extn(std::string const & filename, std::string extn) {
  // std::cerr << filename << " " << extn << std::endl;
  if (exist(filename+'.'+extn)) {
    // std::cerr << "Returning: " << filename+'.'+extn << std::endl;
    return (filename+'.'+extn);
  }
  int p = filename.length();
  while ((p=filename.rfind('.',p-1))!=std::string::npos) {
    std::string testname = filename.substr(0,p+1)+extn;
    if (exist(testname)) {
      // std::cerr << "Returning: " << testname << std::endl;
      return testname;
    }
  }
  // std::cerr << "Returning: " << (filename.substr(0,filename.rfind('.')+1)+extn) << " *" << std::endl;
  return (filename.substr(0,filename.rfind('.')+1)+extn);
}
