#ifndef CHARMAP_H
#define CHARMAP_H

#include <ctype.h>

class charmap {
  char comp[256];
  char canon[256];
  char amino[256];
  char upper[256];

 public:
  inline char term1() const {return '$';}
  inline char term2() const {return '%';}
  inline char term3() const {return '!';}
  inline const char *canonical() const {return canon;}
  inline const char *complement() const {return comp;}
  inline const char *aminoacid() const {return amino;}
  inline const char *uppercase() const {return upper;}

  charmap() {
    for(int i=0; i < 256; ++i) {
      canon[i] = 'N';
      comp[i] = 'N';
      amino[i] = 'X';
      upper[i] = term3();
    }

    canon[term1()] = term1();
    canon[term2()] = term2();

    canon['a'] = canon['A'] = 'A';
    canon['t'] = canon['T'] = 'T';
    canon['c'] = canon['C'] = 'C';
    canon['g'] = canon['G'] = 'G';
    canon['u'] = canon['U'] = 'T';

    // reverse complement mapping
    comp[term1()] = term1();
    comp[term2()] = term2();

    comp['a'] = comp['A'] = 'T';
    comp['t'] = comp['T'] = 'A';
    comp['u'] = comp['U'] = 'A';
    comp['c'] = comp['C'] = 'G';
    comp['g'] = comp['G'] = 'C';

    // amino acid map
    amino[term1()] = term1();
    amino[term2()] = term2();

    amino['a']=amino['A']='A';
    amino['c']=amino['C']='C';
    amino['d']=amino['D']='D';
    amino['e']=amino['E']='E';
    amino['f']=amino['F']='F';
    amino['g']=amino['G']='G';
    amino['h']=amino['H']='H';
    amino['i']=amino['I']='I';
    amino['k']=amino['K']='K';
    amino['l']=amino['L']='L';
    amino['m']=amino['M']='M';
    amino['n']=amino['N']='N';
    amino['p']=amino['P']='P';
    amino['q']=amino['Q']='Q';
    amino['r']=amino['R']='R';
    amino['s']=amino['S']='S';
    amino['t']=amino['T']='T';
    amino['v']=amino['V']='V';
    amino['w']=amino['W']='W';
    amino['y']=amino['Y']='Y';

    // stop codons
    amino['.']=amino['*']=amino['@']='@';

    for (int i='A';i<='Z';i++) {
      upper[i]=upper[tolower(i)]=i;
    }

  }

  void map(const char *m,char *b,char *e) const {
    for( ; b < e ; ++b) {
      *b = m[(unsigned char)(*b)];
    }
  }
  void revmap(const char *m,char *b,char *e) const {
    char tmp;
    for( --e; b < e ; ++b,--e) {
      tmp = m[(unsigned char)(*b)];
      *b = m[(unsigned char)(*e)];
      *e = tmp;
    }
    if (b == e) *b = m[(unsigned char)(*b)];
  }
};

#endif
