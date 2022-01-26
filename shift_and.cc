
#include <assert.h>
#include <string.h>
#include <strings.h>
#include "shift_and.h"
#include "util.h"

shift_and::shift_and(bool wc, bool tn, bool re, unsigned char eos) {
  m_ = 0;
  u_ = 0;
  mask_=0;
  _wc = wc;
  _textn = tn;
  _re = re;
  eos_ = eos;
}

shift_and::~shift_and() {
  clearu();
}

bool shift_and::wildcards() const {
  return _wc;
}

void shift_and::wildcards(bool wc) {
  _wc = wc;
}

bool shift_and::wildcard_text_N() const {
  return _textn;
}

void shift_and::wildcard_text_N(bool tn) {
  _textn = tn;
}

void shift_and::clearu() {
  if (u_) {
    if (u_[0]) delete [] u_[0];
    delete [] u_;
  } 
  if (mask_) delete [] mask_;
  if (m_) delete [] m_;
}

unsigned long shift_and::add_pattern(std::string const & pat, long unsigned id,
				     int esb, int eeb) {
  add_pattern_(pat,id,esb,eeb);
  return id;
}

void shift_and::computeu(CharacterProducer & cp) {
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
  _highbit = (((bigword)1)<<(wordbits-1));
  _wordbits_1 = wordbits-1;
  // cerr << "# " << patterns_count << " L " << patterns_length
  // << " W " << _wordcount << " w " << wordbits << endl;
  // cerr << "_highbit " << binary(_highbit) << endl;
 
  bigword *buffer = new bigword[_wordcount*cp.size()];
  memset(buffer,0,_wordcount*cp.size()*sizeof(bigword));
  u_ = new bigword*[cp.size()];
  for (unsigned int i=0;i<cp.size();i++) {
    u_[i] = buffer+_wordcount*i;
  }
  // cerr << cp.size() << endl;
  mask_=new bigword[_wordcount];
  memset(mask_,0,_wordcount*sizeof(bigword));

  s_=new bigword[_wordcount];
  memset(s_,0,_wordcount*sizeof(bigword));

  m_=new bigword[_wordcount];
  memset(m_,0,_wordcount*sizeof(bigword));

  _patbits = new patbit[patterns_count];
  memset(_patbits,0,patterns_count*sizeof(patbit));

  _patbitind = new unsigned int[_wordcount+1];
  memset(_patbitind,0,(_wordcount+1)*sizeof(unsigned int));

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
      } else if (_re && keyword[i] == '.') {
	for (unsigned int nch1=0;nch1<cp.size();nch1++) {
	  if (cp.ch(nch1) == eos_) {
	    continue;
	  }
	  u_[nch1][wordposition] 
	    |= ( ((bigword)1) << wordbitposition);
	}
      } else if (_re && keyword[i] == ':') {
	for (unsigned int nch1=0;nch1<cp.size();nch1++) {
	  if (strchr("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy",cp.ch(nch1)) || cp.ch(nch1) == eos_) {
	    continue;
	  }
	  u_[nch1][wordposition] 
	    |= ( ((bigword)1) << wordbitposition);
	}
      } else if (_re && keyword[i] == ';') {
	for (unsigned int nch1=0;nch1<cp.size();nch1++) {
	  if (strchr("ACGTacgt",cp.ch(nch1)) || cp.ch(nch1) == eos_) {
	    continue;
	  }
	  u_[nch1][wordposition] 
	    |= ( ((bigword)1) << wordbitposition);
	}
      } else {
	int nch = cp.nch(keyword[i]);
	if (nch >=0 ) {
	  u_[nch][wordposition] 
	    |= ( ((bigword)1) << wordbitposition);
	}
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

  /* cerr << "   s_ ";
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
      cerr << "Word " << i << " bit " << _patbits[j].bit << " pattern " << p << endl;
      p++;
    }
  }
  checkpoint;
  for (int ch=0;ch<cp.size();ch++) {
    cerr << "u_[";
    cerr << ch << "] ";
    for (int i=_wordcount-1;i>=0;i--) {
      cerr << binary(u_[ch][i]);
    }
    cerr << endl;
  }
  
  cerr << "m_ ";
  for (int i=_wordcount-1;i>=0;i--) {
    cerr << binary(m_[i]);
  }
  cerr << endl;
  */
}

void shift_and::reset() {
  memset(m_,0,_wordcount*sizeof(bigword));
}

bool shift_and::find_patterns(CharacterProducer & cp,
			      pattern_hit_vector & pas,
			      long unsigned minpa) {
  unsigned ccount=0;
  long unsigned pacount=0;
  if (cp.eof()) return false;
  unsigned char ch=cp.getnch();
  // cerr << "m: " << binary<unsigned long>(m_) << endl;
  // cerr << "$: " << binary<unsigned long>(mask_) << endl;
  // cerr << "Read: " << ch << endl;
  while (1) {
    for (int i=_wordcount-1;i>=1;i--) {
      m_[i] = ((m_[i] << 1) | (m_[i-1] >> _wordbits_1) | s_[i]) & u_[ch][i];
    }
    m_[0] = ((m_[0] << 1) | s_[0]) & u_[ch][0];

    for (unsigned int i=0;i<_wordcount;i++) {
      // checkpoint;
      // cerr << "i: " << i << endl;
      // cerr << "_patbitind[i]: " << _patbitind[i] << endl;
      // cerr << "_patbitind[i+1]-1: " << _patbitind[i+1]-1 << endl;
      if (m_[i] & mask_[i]) {
	/* There is a hit somewhere in this word... */
	for (unsigned int j=_patbitind[i];j<_patbitind[i+1];j++) {
	  // checkpoint;
	  // cerr << "j: " << j << endl;
	  if (m_[i] & (((bigword)1)<<(_patbits[j].bit))) {
	    pas.push_back(cp.pos(),make_pair(_patbits[j].it,0));
	    pacount++;
	  }
	}
      }
    }
    if (pacount >= minpa) {
      return true;
    }
    if (ccount>1000) {
      report_progress(cp);
      ccount=0;
    } 
    ccount++;
    if (cp.eof()) break;
    ch=cp.getnch();
    // cerr << "Read: " << ch << endl;
  }
  if (pacount > 0) return true;
  return false;
}

