
#include "hash.h"
#include "bits.h"
#include "math.h"
#include "util.h"

contigmod::contigmod(CharacterProducer & cp, int asize, int weight) :
  hash(cp)
{
  a_ = asize;
  w_ = weight;
  mask_ = w_*a_;
  reset();
}

void
contigmod::reset(uint64 p) {
  cp_.reset();
  int n=w_-1;
  if (p > w_) {
    cp_.pos(p-w_);
  } else if (p > 0) {
    n = p-1;
  } 
  h_=0;
  for (int i=0;i<n;i++) {
    if (!cp_.eof()) {
      update(cp_.getnch());
    }
  }
}

hash_t
contigmod::size() const {
  return mask_;
}

hash_t
contigmod::value() const {
  return h_;
}

hash_t
contigmod::rcvalue() const {
  return rc(h_,w_);
}

hash_t
contigmod::cvalue() const {
  hash_t f = value();
  hash_t rc = rcvalue();
  if (f <= rc) {
    return f;
  } 
  return rc;
}

void 
contigmod::update(unsigned char c) {
  h_ *= a_;
  h_ += c;
  h_ %= mask_;
}

bool
contigmod::next() {
  if (cp_.eof()) {
    return false;
  }
  update(cp_.getnch());
  return true;
}

uint64
contigmod::pos() const {
  return cp_.pos();
}

string 
contigmod::str(hash_t v) const {
  string s;
  for (int w=w_-1;w>=0;w--) {
    s += cp_.ch(((int)(v/pow((double)a_,w)))%(int)a_);
  }
  return s;
}

unsigned
contigmod::span() const {
  return w_;
}

unsigned 
contigmod::weight() const {
  return w_;
}

contigshift::contigshift(CharacterProducer & cp, int asize, int weight) :
  hash(cp)
{
  a_ = clg2((uint32)asize);
  w_ = weight;
  mask_ = (((hash_t)1)<<w_*a_)-1;
  ns_ = 0;
  reset();
  timestampi("New contigshift hash with weight ",w_);
}

void
contigshift::reset(uint64 p) {
  cp_.reset();
  int n=w_-1;
  if (p > w_) {
    cp_.pos(p-w_);
  } else if (p > 0) {
    n = p-1;
  } 
  h_=0;
  for (int i=0;i<n;i++) {
    if (!cp_.eof()) {
      update(cp_.getnch());
    }
  }
}

hash_t
contigshift::size() const {
  return mask_+1;
}

hash_t
contigshift::value() const {
  return h_;
}

hash_t
contigshift::rcvalue() const {
  // assert(str(rc(h_,w_)) == reverse_comp(str(h_)));
  return rc(h_,w_);
}

hash_t
contigshift::cvalue() const {
  hash_t f = value();
  hash_t rc = rcvalue();
  if (f <= rc) {
    return f;
  } 
  return rc;
}

void 
contigshift::update(unsigned char c) {
  ns_ = ((c>>a_)&&!(ns_>>15))?(ns_+1):0;
  h_ = (((h_<<a_)|c)&mask_);
}

uint32
contigshift::ns() const {
  return ns_;
}

bool
contigshift::next() {
  if (cp_.eof()) {
    return false;
  }
  update(cp_.getnch());
  return true;
}

uint64
contigshift::pos() const {
  return cp_.pos();
}

string 
contigshift::str(hash_t v) const {
  string s;
  for (int w=w_-1;w>=0;w--) {
    s += cp_.ch((v>>(a_*w))&((((hash_t)1)<<a_)-1));
  }
  return s;
}

unsigned 
contigshift::span() const {
  return w_;
}

unsigned 
contigshift::weight() const {
  return w_;
}

spaced::spaced(CharacterProducer & cp, int asize, 
	       const char* tmplstr, signed char pin) :
  hash(cp)
{
  uint64* templ=0;
  int len=strlen(tmplstr);
  int blocks = len/64;
  if (len%64 != 0) {
    blocks += 1;
  }
  templ = new uint64[blocks];
  // timestampi("blocks  = ",blocks);
  zero(templ,blocks);
  for (int i=0;i<len;i++) {
    assert(tmplstr[i] == '1' || tmplstr[i] == '0');
    if (tmplstr[i] != '0') {
      set(templ,blocks,i);
    }
  }
  assert(get(templ,blocks,len-1));
  assert(get(templ,blocks,0));
  a_ = clg2((uint32)asize);
  // timestampi("a_  = ",a_);
  w_ = pop(templ,blocks);
  // timestampi("w_  = ",w_);
  assert(w_);
  //hash_t x;
  // assert(w_*a_ <= bits(x));
  s_ = bits(templ,blocks)-nlz(templ,blocks);
  // timestampi("s_  = ",s_);
  if (pin == -1) {
    p_ = s_;
    pw_ = w_;
  } else if (pin == 0) {
    unsigned short pin,pwin;
    period(templ,blocks,pin,pwin);
    p_ = pin;
    pw_ = pwin;
  } else {
    assert(periodIs(pin,templ,blocks));
    p_ = pin;
    pw_ = periodWeight(pin,templ,blocks);
  }
  // timestampi("p_  = ",p_);
  // timestampi("pw_ = ",pw_);
  h_ = new hash_t[p_];
  hind_ = new unsigned[p_*pw_];
  hptr_ = new unsigned[p_+1];
  for (int p=0;p<p_+1;p++) {
    hptr_[p]=p*pw_;
  }
  for (int pos=0;pos<p_;pos++) {
    for (int h=0;h<p_;h++) {
      uint32 p=(s_-1+pos-h)%p_;
      if (get(templ,blocks,p)) {
	hind_[hptr_[pos]] = h;
	hptr_[pos]++;
      }
    }
  }
  for (int p=0;p<p_+1;p++) {
    hptr_[p]=p*pw_;
  }
  /* for (int s=0;s<p_;s++) {
    cerr << s << ": ";
    for (int i=hptr_[s];i<hptr_[s+1];i++) {
      cerr << i << ":"<< (int)hind_[i] << " ";
    }
    cerr << endl;
  }
  // timestampi("init to: ",(p_-(s_%p_))%p_); */
  mask_ = (((hash_t)1)<<w_*a_)-1;
  reset();
  timestamps("New periodic spaced hash with template ",tmplstr)
}

spaced::~spaced() {
  delete [] h_;
  delete [] hind_;
  delete [] hptr_;
}

void
spaced::reset(uint64 p) {
  cp_.reset();
  int n=s_-1;
  if (p > s_) {
    cp_.pos(p-s_);
  } else if (p > 0) {
    n = p-1;
  } 
  fill(h_,h_+p_,0);
  ptr_ = (p_-(s_%p_))%p_;
  for (int i=0;i<n;i++) {
    if (!cp_.eof()) {
      update(cp_.getnch());
    }
  }
}

hash_t
spaced::size() const {
  return mask_+1;
}

hash_t
spaced::value() const {
  return h_[ptr_];
}

hash_t
spaced::rcvalue() const {
  // assert(str(rc(h_,w_)) == reverse_comp(str(h_)));
  return rc(h_[ptr_],w_);
}

hash_t
spaced::cvalue() const {
  hash_t f = value();
  hash_t rc = rcvalue();
  if (f <= rc) {
    return f;
  } 
  return rc;
}

void 
spaced::update(unsigned char c) {
  if (++ptr_==p_) ptr_ = 0;
  int stop=hptr_[ptr_+1];
  for (int i=hptr_[ptr_];i<stop;i++) {
    hash_t *h = h_+hind_[i];
    (*h) = (((*h)<<a_)|c);
  }
  h_[ptr_] &= mask_;
}

bool
spaced::next() {
  if (cp_.eof()) {
    return false;
  }
  update(cp_.getnch());
  return true;
}

uint64
spaced::pos() const {
  return cp_.pos();
}

string 
spaced::str(hash_t v) const {
  string s;
  for (int w=w_-1;w>=0;w--) {
    s += cp_.ch((v>>(a_*w))&((((hash_t)1)<<a_)-1));
  }
  return s;
}

unsigned 
spaced::span() const {
  return s_;
}

unsigned 
spaced::weight() const {
  return w_;
}


shiftspaced::shiftspaced(CharacterProducer & cp, int asize, const char* tmplstr,const char* type) :
  hash(cp)
{
  uint64 templ;
  zero(templ);
  for (int i=0;i<strlen(tmplstr);i++) {
    assert(tmplstr[i] == '1' || tmplstr[i] == '0');
    if (tmplstr[i] != '0') {
      set(templ,i);
    }
  }
  a_ = clg2((uint32)asize);
  // timestampi("a_  = ",a_);
  w_ = pop(templ);
  // timestampi("w_  = ",w_);
  assert(w_);
  s_ = bits(templ)-nlz(templ);
  // timestampi("s_  = ",s_);
  size_ = (((hash_t)1)<<w_*a_);
  unsigned short *run = new unsigned short[w_];
  unsigned short nrun = runs(templ,run);
  // timestampi("nrun = ",nrun);
  assert(nrun<=2*w_);
  unsigned short *pos = new unsigned short[nrun];
  pos[nrun-1] = run[nrun-1];
  for (int i=nrun-2;i>=0;i--) {
    pos[i] = run[i]+pos[i+1];
  }
//   for (int i=0;i<nrun;i++) {
//    cerr << run[i] << " ";
//   }
//   cerr << endl;
//   for (int i=0;i<nrun;i++) {
//    cerr << pos[i] << " ";
//   }
//   cerr << endl;

  nshift_ = (nrun+1)/2;
  mask_ = new uint64[nshift_];
  shift_ = new unsigned[nshift_];
  mask_[0] = (((uint64)1)<<(pos[nrun-1]*a_))-1;
  shift_[0] = 0;
  for (int i=1;i<nshift_;i++) {
    mask_[i] = mask64(pos[nrun-1-(2*i-1)]*a_,pos[nrun-1-(2*i)]*a_);
    shift_[i] = run[nrun-1-(2*i-1)]*a_+shift_[i-1];
  }
  
//   for (int i=0;i<nshift_;i++) {
//    cerr << tostr(mask_[i]) << "\n";
//   }
//   cerr << endl;
//   for (int i=0;i<nshift_;i++) {
//    cerr << shift_[i] << " ";
//   }
//   cerr << endl;

  delete [] run;
  delete [] pos;
  reset();
  if (type != NULL) {
    char buffer[1024];
    sprintf(buffer,"New %s shifted spaced hash with template ",type);
    timestamps(buffer,tmplstr);
  } else {
    timestamps("New shifted spaced hash with template ",tmplstr);
  }
}

asymmetric_shiftspaced::asymmetric_shiftspaced(CharacterProducer & cp, int asize, const char* tmplstr):
  shiftspaced(cp,asize,tmplstr,"asymmetric") {
}

shiftspaced::~shiftspaced() {
  delete [] mask_;
  delete [] shift_;
}

void
shiftspaced::reset(uint64 p) {
  cp_.reset();
  int n=s_-1;
  if (p > s_) {
    cp_.pos(p-s_);
  } else if (p > 0) {
    n = p-1;
  } 
  h0_ = 0;
  for (int i=0;i<n;i++) {
    if (!cp_.eof()) {
      update(cp_.getnch());
    }
  }
}

hash_t
shiftspaced::size() const {
  return size_;
}

hash_t
shiftspaced::value() const {
  return h_;
}

hash_t
shiftspaced::rcvalue() const {
  return rc(h_,w_);
}

hash_t
asymmetric_shiftspaced::rcvalue() const {
  return hrc_;
}

hash_t
shiftspaced::cvalue() const {
  hash_t f = value();
  hash_t rc = rcvalue();
  if (f <= rc) {
    return f;
  } 
  return rc;
}

void 
shiftspaced::update(unsigned char c) {
  h0_   = ((h0_<<a_)|c);
  h_   = h0_&mask_[0];
  for (int i=1;i<nshift_;i++) {
    h_   |= ((h0_&mask_[i])>>shift_[i]);
  }
}

void 
asymmetric_shiftspaced::update(unsigned char c) {
  h0_   = ((h0_<<a_)|c);
  h0rc_ = rc(h0_,s_);
  h_   = h0_&mask_[0];
  hrc_ = h0rc_&mask_[0];
  for (int i=1;i<nshift_;i++) {
    h_   |= ((h0_&mask_[i])>>shift_[i]);
    hrc_ |= ((h0rc_&mask_[i])>>shift_[i]);
  }
}

bool
shiftspaced::next() {
  if (cp_.eof()) {
    return false;
  }
  update(cp_.getnch());
  return true;
}

uint64
shiftspaced::pos() const {
  return cp_.pos();
}

string
shiftspaced::str(hash_t v) const {
  return str(v,w_);
}

string 
shiftspaced::str(hash_t v, unsigned short w0) const {
  string s;
  for (int w=w0-1;w>=0;w--) {
    s += cp_.ch((v>>(a_*w))&((((hash_t)1)<<a_)-1));
  }
  return s;
}

unsigned 
shiftspaced::span() const {
  return s_;
}

unsigned 
shiftspaced::weight() const {
  return w_;
}


hashset::hashset(CharacterProducer & cp, int asize, 
		 const char* tmplstr) :
  hash(cp)
{
  char *tstr;
  tstr=new char[strlen(tmplstr)+1];
  strcpy(tstr, tmplstr);

  char *p;
  
  int len = strlen(tstr);
  n_ = 0;
  for (int i=0;i<len;i++) {
    if (tstr[i]==';') {
      n_++;
    }
  }
  n_++;
  assert(n_ > 0);
  h_ = new hash*[n_];
  p = strtok(tstr,";");
  int i = 0;
  while (p != NULL) {
    if (!strchr(p,'0')) {
      h_[i] = new contigshift(cp,asize,strlen(p));      
    } else {
      h_[i] = spacedselect(cp,asize,p);
    }
    p = strtok(NULL,";");
    i++;
  }
  assert(i == n_);

  // span is minimum span.
  s_ = h_[0]->span();
  for (int i=1;i<n_;i++) {
    assert(h_[i]->span() >= s_);
  }
  delete [] tstr;
} 

hashset::~hashset() {
  for (int i=0;i<n_;i++) {
    delete h_[i];
  }
  delete [] h_;
}

void
hashset::reset(uint64 p) {
  for (int i=0;i<n_;i++) {
    if (p == 0) {
      h_[i]->reset(s_);
    } else {
      h_[i]->reset(p);
    }
  }
  p_ = n_-1;
}

hash_t
hashset::size() const {
  return h_[0]->size();
}

hash_t
hashset::value() const {
  return h_[p_]->value();
}

hash_t
hashset::rcvalue() const {
  return h_[p_]->rcvalue();
}

hash_t
hashset::cvalue() const {
  return h_[p_]->cvalue();
}

void 
hashset::update(unsigned char c) {
  for (int i=0;i<n_;i++) {
    h_[i]->update(c);
  }
}

bool
hashset::next() {
  if (++p_!=n_) {
    return true;
  } 
  p_ = 0;
  if (cp_.eof()) {
    return false;
  }
  update(cp_.getnch());
  return true;
}

uint64
hashset::pos() const {
  return cp_.pos();
}

string 
hashset::str(hash_t v) const {
  return h_[0]->str(v);
}

unsigned 
hashset::span() const {
  return h_[p_]->span();
}

unsigned 
hashset::weight() const {
  return h_[0]->weight();
}

bool
hashset::symmetric() const {
  for (int i=0;i<n_;i++) {
    if (!h_[i]->symmetric()) {
      return false;
    }
  }
  return true;
}

taghashset::taghashset(CharacterProducer & cp, int asize,
		       const char* tmplstr, int tsize) :
  hash(cp)
{
  ab_ = clg2((uint32)asize);
  
  string tstr(tmplstr);
  int len = tstr.length();
  n_ = 0;
  tn_ = 0;
  for (int i=0;i<len;i++) {
    if (tstr[i]==';') {
      n_++;
    }
    if (tstr[i]==':' || tstr[i]==',') {
      tn_++;
    }
  }
  n_++;

  assert(n_ > 0);
  assert(tn_ > 0);
  
  h_ = new hash*[n_];
  t_ = new hash*[tn_];
  fill(h_,h_+n_,(hash*)0);  
  fill(t_,t_+tn_,(hash*)0);

  int i = 0;
  int p0 = 0;
  int p1,p2,p3;
  p1 = tstr.find_first_of(":",p0);
  while (p1 != string::npos) {
    if (tstr.substr(p0,p1-p0).find_first_of("0") == string::npos) {
      h_[i] = new contigshift(cp,asize,p1-p0);      
    } else {
      h_[i] = spacedselect(cp,asize,tstr.substr(p0,p1-p0).c_str());
    }
    p2 = p1+1;
    while (((p3 = tstr.find_first_of(",;",p2))!=string::npos)&&tstr[p3]!=';') {
      ostringstream oss;
      oss << "Tag " << atoi(tstr.substr(p2,p3-p2).c_str()) << " associated with hash " << i << ends;
      timestamps("",oss.str().c_str());
      t_[atoi(tstr.substr(p2,p3-p2).c_str())] = h_[i];
      p2 = p3+1;
    }
    if (p3 != string::npos) {
      ostringstream oss;
      oss << "Tag " << atoi(tstr.substr(p2,p3-p2).c_str()) << " associated with hash " << i << ends;
      timestamps("",oss.str().c_str());
      t_[atoi(tstr.substr(p2,p3-p2).c_str())] = h_[i];
    } else {
      ostringstream oss;
      oss << "Tag " << atoi(tstr.substr(p2).c_str()) << " associated with hash " << i << ends;
      timestamps("",oss.str().c_str());
      t_[atoi(tstr.substr(p2).c_str())] = h_[i];
      break;
    }
    p0 = p3 + 1;
    p1 = tstr.find(":",p0);
    i++;
  }

  assert(i == n_-1);
  

  if (tsize <= 0) {
    tb_ = clg2((uint32)tn_);
    ts_ = tn_;
  } else {
    tb_ = clg2((uint32)tsize);
    ts_ = tsize;
  }

  hb_ = weight()*ab_;
  assert(hb_ <= bits((hash_t)0));
  mask_ = (((hash_t)1)<<(hb_+tb_))-1;
  hm_ = (((hash_t)1)<<hb_)-1;

  tpm_ = new hash_t[tn_];
  fill(tpm_,tpm_+tn_,0);
  for (int i=0;i<tn_;i++) {
    assert(t_[i]);
    tpm_[i] = (((hash_t)i)<<hb_)&mask_;;
  }

  // span is minimum span.
  s_ = h_[0]->span();
  for (int i=1;i<n_;i++) {
    assert(h_[i]->span() >= s_);
  }

} 

taghashset::~taghashset() {
  for (int i=0;i<n_;i++) {
    delete h_[i];
  }
  delete [] h_;
  delete [] t_;
}

void
taghashset::reset(uint64 p) {
  for (int i=0;i<n_;i++) {
    if (p == 0) {
      h_[i]->reset(s_);
    } else {
      h_[i]->reset(p);
    }
  }
  tp_ = tn_-1;
}

hash_t
taghashset::size() const {
  return ts_*h_[0]->size();
}

hash_t
taghashset::value() const {
  return (tpm_[tp_]|t_[tp_]->value());
}

hash_t
taghashset::rcvalue() const {
  return (tpm_[tp_]|t_[tp_]->rcvalue());
}

hash_t
taghashset::cvalue() const {
  return  (tpm_[tp_]|t_[tp_]->cvalue());
}

void 
taghashset::update(unsigned char c) {
  for (int i=0;i<n_;i++) {
    h_[i]->update(c);
  }
}

bool
taghashset::next() {
  if (++tp_!=tn_) {
    return true;
  } 
  tp_ = 0;
  if (cp_.eof()) {
    return false;
  }
  update(cp_.getnch());
  return true;
}

uint64
taghashset::pos() const {
  return cp_.pos();
}

string 
taghashset::str(hash_t v) const {
  ostringstream oss;
  oss << (v>>hb_) << ":" << h_[0]->str(v&hm_) << ends;
  return oss.str();
}

unsigned 
taghashset::span() const {
  return t_[tp_]->span();
}

unsigned 
taghashset::weight() const {
  return h_[0]->weight();
}

bool
taghashset::symmetric() const {
  for (int i=0;i<n_;i++) {
    if (!h_[i]->symmetric()) {
      return false;
    }
  }
  return true;
}

hash*
hashselect(CharacterProducer & cp, int asize, const char *tstr) {
  timestamps("Hash select with string: ",tstr);
  if (index(tstr,':')) {
    return new taghashset(cp,asize,tstr);
  } else if (index(tstr,';')) {
    return new hashset(cp,asize,tstr);
  } else if (strlen(tstr) > 2 && strspn(tstr,"01") == strlen(tstr) && index(tstr,'0')) {
    return spacedselect(cp,asize,tstr);
  } else if (strlen(tstr) > 2 && strspn(tstr,"1") == strlen(tstr)) {
    return new contigshift(cp,asize,strlen(tstr));
  }
  return new contigshift(cp,asize,atoi(tstr));
}

hash*
spacedselect(CharacterProducer & cp, int asize, const char *tstr) {
  uint64 templ = tobv64(tstr);
  unsigned short ns = nshift(templ);  
  unsigned short pd,pdwt;
  period(templ,pd,pdwt);
  float shiftcost = 9.4*ns+63.4;
  float periodcost = 10.5*pdwt+68.6;
  //timestampd("Running time estimate for shifted spaced hash: ",shiftcost);
  //timestampd("Running time estimate for periodic spaced hash: ",periodcost);
  if (shiftcost <= periodcost) {
    unsigned short s  = 64-nlz(templ);
    uint64 revtempl = reverse(templ,s);
    // cerr << tostr(templ) << endl;
    // cerr << tostr(revtempl) << endl;
    if (templ == revtempl) {
      return new shiftspaced(cp,asize,tstr);
    } else {
      return new asymmetric_shiftspaced(cp,asize,tstr);
    }
  } else {
    return new spaced(cp,asize,tstr);
  }
}
