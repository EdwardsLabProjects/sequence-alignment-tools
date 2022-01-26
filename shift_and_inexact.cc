
#include <assert.h>
#include <string.h>
#include "shift_and_inexact.h"
#include "util.h"

shift_and_inexact::shift_and_inexact(unsigned int k, unsigned char eos, bool wc, bool tn, bool id, bool dm) {
  m_ = 0;
  u_ = 0;
  mask_=0;
  k_ = k;
  eos_ = eos;
  _wc = wc;
  _textn = tn;
  _indels = id;
  _dna_mut = dm;
}

shift_and_inexact::~shift_and_inexact() {
  clearu();
}

unsigned int shift_and_inexact::mismatches() const {
  return k_;
}

void shift_and_inexact::mismatches(unsigned int k) {
  k_ = k;
}

bool shift_and_inexact::wildcards() const {
  return _wc;
}

void shift_and_inexact::wildcards(bool wc) {
  _wc = wc;
}

bool shift_and_inexact::wildcard_text_N() const {
  return _textn;
}

void shift_and_inexact::wildcard_text_N(bool tn) {
  _textn = tn;
}

bool shift_and_inexact::indels() const {
  return _indels;
}

void shift_and_inexact::indels(bool id) {
  _indels = id;
}

bool shift_and_inexact::dna_mut() const {
  return _dna_mut;
}

void shift_and_inexact::dna_mut(bool dm) {
  _dna_mut = dm;
}

char shift_and_inexact::eos_char() const {
  return eos_;
}

void shift_and_inexact::eos_char(char c) {
  eos_ = c;
}

void shift_and_inexact::clearu() {
  if (u_) {
    if (u_[0]) delete [] u_[0];
    delete [] u_;
  } 
  if (m_) {
    if (m_[0]) delete [] m_[0];
    delete [] m_;
  } 
  if (mask_) delete [] mask_;
}

unsigned long shift_and_inexact::add_pattern(std::string const & pat, 
					     long unsigned id, 
					     int esb, int eeb) {
  add_pattern_(pat,id,esb,eeb);
  return id;
}

void shift_and_inexact::computeu(CharacterProducer & cp) {
  clearu();
  long unsigned int patterns_length=0;
  long unsigned int patterns_count=0;
  tinylist<pattern_list_element>::const_iterator it;
  for (it=patterns().begin();it!=patterns().end();++it) {
    patterns_length+=it->pattern().length();
    patterns_count++;
  }
  unsigned int wordbits = (sizeof(bigword)*8);
  long unsigned int patterns_wordcount;
  patterns_wordcount = patterns_length/wordbits + 
                       ((patterns_length%wordbits)?1:0);
  // assert(patterns_wordcount*wordbits >= patterns_length);
  _wordcount = patterns_wordcount;
  // cerr << "# " << patterns_count << " L " << patterns_length
  // << " W " << _wordcount << " w " << wordbits << endl;
  _highbit = (((bigword)1)<<(wordbits-1));
  // cerr << "_highbit " << binary(_highbit) << endl;
  
  bigword *buffer = new bigword[_wordcount*cp.size()];
  memset(buffer,0,_wordcount*cp.size()*sizeof(bigword));
  u_ = new bigword*[cp.size()];
  for (unsigned int i=0;i<cp.size();i++) {
    u_[i] = buffer+_wordcount*i;
  }
  mask_=new bigword[_wordcount];
  memset(mask_,0,_wordcount*sizeof(bigword));
  s_=new bigword[_wordcount];
  memset(s_,0,_wordcount*sizeof(bigword));
  buffer=new bigword[_wordcount*(k_+1)];
  memset(buffer,0,_wordcount*(k_+1)*sizeof(bigword));
  m_ = new bigword*[k_+1];
  for (unsigned int i=0;i<=k_;i++) {
    m_[i] = buffer+_wordcount*i;
  }
  _patbits = new patbit[patterns_count];
  memset(_patbits,0,patterns_count*sizeof(patbit));
  _patbitind = new long unsigned int[_wordcount+1];
  memset(_patbitind,0,(_wordcount+1)*sizeof(long unsigned int));
  
  eos_ = cp.nch(eos_);

  unsigned int bitposition=0;
  unsigned int wordposition=0;
  unsigned int wordbitposition=0;
  unsigned int patbitsposition=0;
  _patbitind[0] = 0;

  for (it=patterns().begin();it!=patterns().end();++it) {
    std::string const & keyword = it->pattern();
    unsigned int pbits=keyword.length();
    
    for (unsigned int i=0;i<pbits;i++) {
      char *wccompat;
      if (_wc && ((wccompat=iupac_compatible(keyword[i])) != 0)) {
	unsigned int j=0;
	while (wccompat[j]) {
	  int nch1 = cp.nch(wccompat[j]);
	  if (nch1 >= 0 && (wccompat[j]!='N' || _textn)) {
	    u_[nch1][wordposition] 
	      |= ( ((bigword)1) << wordbitposition);
	  }
	  j++;
	}
      } else {
	int nch = cp.nch(keyword[i]);
	if (nch >=0 ) {
	  u_[nch][wordposition] 
	    |= ( ((bigword)1) << wordbitposition);
	}
      }
      for (unsigned int k=i+1;k<=k_;k++) {
	m_[k][wordposition] |= (((bigword)1) << wordbitposition);
      }
      if (i==0) { /* first position of this pattern */
	s_[wordposition] |= (((bigword)1) << wordbitposition);
	// cerr << wordposition << " " << binary(s_[wordposition]) << endl;
      }
      if (i==(pbits-1)) { /* last position of this pattern */
	mask_[wordposition] |= (((bigword)1) << wordbitposition);
	_patbits[patbitsposition].bit = wordbitposition;
	_patbits[patbitsposition].it = it;
	patbitsposition++;
      }
      bitposition++;
      wordposition = bitposition/wordbits;
      wordbitposition = bitposition%wordbits;
      if (wordbitposition%wordbits==0) {
	_patbitind[wordposition] = patbitsposition;
      }
    }
  }
  _patbitind[_wordcount] = patbitsposition;
  /* 
  cerr << "   s_ ";
  for (int i=_wordcount-1;i>=0;i--) {
    cerr << binary(s_[i]);
  }
  cerr << endl;
  cerr << "mask_ ";
  for (int i=_wordcount-1;i>=0;i--) {
    cerr << binary(mask_[i]);
  }
  cerr << endl;
  int p=0;
  for (int i=0;i<_wordcount;i++) {
    for (int j=_patbitind[i];j<_patbitind[i+1];j++) {
      cerr << "Word " << i << " bit " << j << " pattern " << p << endl;
      p++;
    }
  }
  for (int ch=0;ch<cp.size();ch++) {
    cerr << "u_[";
    cerr << cp.ch(ch) << "] ";
    for (int i=_wordcount-1;i>=0;i--) {
      cerr << binary(u_[ch][i]);
    }
    cerr << endl;
  }

  for (unsigned int k=0;k<=k_;k++) {
    cerr << "m_[";
    cerr << k << "] ";
    for (int i=_wordcount-1;i>=0;i--) {
      cerr << binary(m_[k][i]);
    }
    cerr << endl;
  }
  */
}


void shift_and_inexact::reset() {
  memset(m_,0,_wordcount*sizeof(bigword));
  tinylist<pattern_list_element>::const_iterator it;
  unsigned int bitposition=0;
  unsigned int wordposition=0;
  unsigned int wordbitposition=0;
  unsigned int patbitsposition=0;
  unsigned int wordbits = (sizeof(bigword)*8);
  for (it=patterns().begin();it!=patterns().end();++it) {
    std::string const & keyword = it->pattern();
    unsigned int pbits=keyword.length();
    
    for (unsigned int i=0;i<pbits;i++) {
      for (unsigned int k=i+1;k<=k_;k++) {
	m_[k][wordposition] |= (((bigword)1) << wordbitposition);
      }
      bitposition++;
      wordposition = bitposition/wordbits;
      wordbitposition = bitposition%wordbits;
      if (wordbitposition%wordbits==0) {
	_patbitind[wordposition] = patbitsposition;
      }
    }
  }
}

bool shift_and_inexact::find_patterns(CharacterProducer & cp,
				      pattern_hit_vector & pas,
				      long unsigned minpa) {
  // checkpoint;
  unsigned ccount=0;
  long unsigned pacount=0;
  FILE_POSITION_TYPE lastpapos=0;
  if (cp.eof()) return false;
  static bigword *m0=0,*m1=0,*m2=0,*m3=0; 
  if (!m0) {
    m0 = new bigword[4*_wordcount];
    m1 = m0+1*_wordcount;
    m2 = m0+2*_wordcount;
    m3 = m0+3*_wordcount;
  }
  unsigned char ch=cp.getnch();
  while (1) {
    bigword * const & m0_ = m_[0];
    if (_indels) {
      memcpy(m0,m0_,sizeof(bigword)*_wordcount);
    }
    for (int i=_wordcount-1;i>=1;i--) {
      m1[i] = ((m0_[i] << 1) | ((m0_[i-1]&_highbit)?1:0)) | s_[i] ;
    }
    m1[0] = ((m0_[0] << 1) | s_[0]);
    
    bigword * const & uch = u_[ch];
    for (unsigned int i=0;i<_wordcount;i++) {
      m0_[i] = m1[i] & uch[i];
      if (_indels) {
	m1[i] |= m0[i];
      }
    }
    for (unsigned int l=1;l<=k_;l++) {
      bigword * const & ml_ = m_[l];
      if (_indels) {
	memcpy(m2,ml_,sizeof(bigword)*_wordcount);
      }
      for (int i=_wordcount-1;i>=1;i--) {
	m3[i] = ((ml_[i] << 1) | ((ml_[i-1]&_highbit)?1:0) | s_[i]);
	ml_[i] = m3[i] & uch[i];
      }
      m3[0] = ((ml_[0] << 1) | s_[0]);
      ml_[0] = m3[0] & uch[0];
      if (ch != eos_) {
	for (int i=_wordcount-1;i>=1;i--) {
	  ml_[i] |= m1[i];
	  if (_indels) {
	    ml_[i] |= ((m_[l-1][i] << 1)|((m_[l-1][i-1]&_highbit)?1:0) | s_[i]);
	    ml_[i] |= m_[l-1][i];
	  }
	  m1[i] = m3[i];
	  if (_indels) {
	    m1[i] |= m2[i];
	  }
	} 
	ml_[0] |= m1[0];
	if (_indels) {
	  ml_[0] |= ((m_[l-1][0] << 1) | s_[0]);
	  ml_[0] |= (m_[l-1][0] | m1[0]);
	}
	m1[0] = m3[0];
	if (_indels) {
	  m1[0] |= m2[0];
	}
      }
    }
    for (unsigned int i=0;i<_wordcount;i++) {
      if (m_[k_][i] & mask_[i]) {
	/* There is a hit somewhere in this word... */
	for (unsigned int j=_patbitind[i];j<_patbitind[i+1];j++) {
	  if (m_[k_][i] & (((bigword)1)<<(_patbits[j].bit))) {
	    // checkpoint;
	    // cerr << _patbits[j].it->id() << " " << cp.pos() << " " << k_ << endl;
	    int k=k_-1;
	    while (k >= 0 && (m_[k][i] & (((bigword)1)<<(_patbits[j].bit)))) {
	      // checkpoint;
	      k--;
	    }
	    k++;
	    // cerr << _patbits[j].it->id() << " " << cp.pos() << " " << k << endl;
	    pas.push_back(cp.pos(),make_pair(_patbits[j].it,k));
	    pacount++;
	    lastpapos=cp.pos();
	  }
	}
      }
    }
    if (pacount >= minpa && cp.pos() > lastpapos+1) {
      return true;
    }
    if (ccount>1000) {
      report_progress(cp);
      ccount=0;
    } 
    ccount++;
    if (cp.eof()) break;
    ch=cp.getnch();
  }
  if (pacount > 0) {
    return true;
  }
  return false;
}


