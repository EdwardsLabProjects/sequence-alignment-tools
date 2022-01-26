
#include "exact_halves.h"
#include "primer_alignment.h"
// #include "sortedvector.t"

exact_halves::exact_halves(PatternMatch *pm,
			   unsigned int k, unsigned char eos, 
			   bool wc, bool tn, bool id, bool dm)
  : pm_(pm), num_patterns_(0), k_(k), eos_(eos), _wc(wc), _textn(tn), _indels(id), _dna_mut(dm)
{
}

exact_halves::exact_halves(PatternMatch *pm)
  : pm_(pm), num_patterns_(0), k_(0), eos_('\n'), _wc(false), _textn(false), _indels(true), _dna_mut(false)
{
}

exact_halves::~exact_halves() {
  delete pm_;
}

unsigned int exact_halves::mismatches() const {
  return k_;
}

void exact_halves::mismatches(unsigned int k) {
  k_ = k;
}

bool exact_halves::wildcards() const {
  return _wc;
}

void exact_halves::wildcards(bool wc) {
  _wc = wc;
}

bool exact_halves::wildcard_text_N() const {
  return _textn;
}

void exact_halves::wildcard_text_N(bool tn) {
  _textn = tn;
}

bool exact_halves::indels() const {
  return _indels;
}

void exact_halves::indels(bool id) {
  _indels = id;
}

unsigned char exact_halves::eos_char() const {
  return eos_;
}

void exact_halves::eos_char(unsigned char c) {
  eos_ = c;
}

long unsigned int 
exact_halves::add_pattern(std::string const & pat, unsigned long id,
			  int esb, int eeb) {
  add_pattern_(pat,id,esb,eeb);
  num_patterns_++;
  return id;
}

void 
exact_halves::del_pattern(CharacterProducer & cp, std::string const & pat, 
			  long unsigned int id) {
  // checkpoint;
  int patlen = pat.length();
  std::string const & patl = pat.substr(0,patlen/2);
  std::list<long unsigned int> l0;
  pm_->lookup_pattern(cp,patl,l0);
  std::list<long unsigned int>::iterator it;
  it = l0.begin();
  while (it != l0.end()) {
    if (plit_[(*it)]->pattern() == pat && (id == 0 || plit_[(*it)]->id() == id)) {
      pm_->del_pattern(cp,patl,(*it));
    }
    ++it;
  }
  std::string const & patr = pat.substr(patlen/2);
  l0.clear();
  pm_->lookup_pattern(cp,patr,l0);
  it = l0.begin();
  while (it != l0.end()) {
    if (plit_[(*it)]->pattern() == pat && (id == 0 || plit_[(*it)]->id() == id)) {
      pm_->del_pattern(cp,patr,(*it));
    }
    ++it;
  }
}

void exact_halves::lookup_pattern(CharacterProducer & cp, std::string const & pat,
				  std::list<long unsigned int> & l) {
  int patlen = pat.length();
  std::string const & patl = pat.substr(0,patlen/2);
  std::list<long unsigned int> l0;
  pm_->lookup_pattern(cp,pat,l0);
  std::list<long unsigned int>::iterator it;
  it = l0.begin();
  while (it != l0.end()) {
    if (plit_[(*it)]->pattern() == pat) {
      l.push_back(plit_[(*it)]->id());
    }
    ++it;
  }
}

bool exact_halves::hit_lessthan(pattern_hit_vector::element const & a, 
			   pattern_hit_vector::element const & b) {
  if (a.key() != b.key()) return (a.key()<b.key());
  return (a.value().first->id() > b.value().first->id());
}

bool 
exact_halves::find_patterns(CharacterProducer & cp, 
			    pattern_hit_vector & phs,
			    long unsigned minka) {


  pattern_hit_vector l(minka);
  pattern_hit_vector::iterator it;
  primer_alignment_lmatch pal;
  primer_alignment_rmatch par;
  pal.eos(eos_); pal.kmax(k_); pal.wc(_wc); 
  pal.tn(_textn); pal.indels(_indels); pal.dna_mut(_dna_mut);
  pal.maxpatlen(_mpl); pal.yesno(true);
  par.eos(eos_); par.kmax(k_); par.wc(_wc); 
  par.tn(_textn); par.indels(_indels); par.dna_mut(_dna_mut);
  par.maxpatlen(_mpl); par.yesno(true);
  bool more;
  // checkpoint;
  while ((more=pm_->find_patterns(cp,l,minka))||!l.empty()) {
    // checkpoint;
    FILE_POSITION_TYPE oldcharspos;
    oldcharspos = cp.pos();
    l.normalize(hit_lessthan);
    // checkpoint;
    it = l.begin();
    while (it != l.end()) {
      long unsigned int pid(it->value().first->id());
      FILE_POSITION_TYPE pos(it->key());
      // checkpoint;
      // cerr << it->value()->id() << " " << it->value()->pattern() << " " << it->key() << endl;
      tinylist<pattern_list_element>::const_iterator const & pit(plit_[pid]);
      int esb = pit->exact_start_bases();
      int eeb =  pit->exact_end_bases();
      if (pid%2==1) {
	// Exact match to first half
	// checkpoint;
	pal.reset();
	pal.pos(pos);
	pal.exact_start_bases(esb);
	pal.exact_end_bases(eeb);
	if (pal.align(cp,pattern_halves_[pid],pattern_halves_[pid+1])) {
	  // checkpoint;
	  // cerr << pal.end() << " " << pattern_halves_[pid] << " " << pattern_halves_[pid+1] << endl;
	  if (pal.end() > lasthit_[(pid+1)/2]+((_indels)?(2*k_):0)) {
	    phs.push_back(pal.end(),make_pair(pit,pal.value()));
	    lasthit_[(pid+1)/2] = pal.end();
	  }
	}
      } else {
	// Exact match to second half 
	// checkpoint;
	par.reset();
	par.pos(pos);
	par.exact_start_bases(esb);
	par.exact_end_bases(eeb);
	if (par.align(cp,pattern_halves_[pid-1],pattern_halves_[pid])) {
	  // checkpoint;
	  // cerr << pattern_halves_[pid-1] << " " << pattern_halves_[pid] << endl;
	  if (par.end() > lasthit_[pid/2]+((_indels)?(2*k_):0)) {
	    phs.push_back(par.end(),make_pair(pit,par.value()));
	    lasthit_[pid/2] = par.end();
	  }
	}
      }
      // checkpoint;
      ++it;
      // checkpoint;
    }
    // checkpoint;
    l.clear();
    // checkpoint;
    cp.pos(oldcharspos);
    report_progress(cp);
    if (phs.size() >= minka || 
	(more==false && phs.size() > 0)) return true;
  }
  return false;
}

void 
exact_halves::init(CharacterProducer & cp) {
  assert(pm_!=((void*)0));
  plit_.resize(num_patterns_*2+1);
  pattern_halves_.resize(num_patterns_*2+1);
  lasthit_.resize(num_patterns_+1);
  tinylist<pattern_list_element>::const_iterator it;
  long unsigned int id=0;
  _mpl=0;
  for (it=patterns().begin();it!=patterns().end();++it) {
    int patlen = it->pattern().length();
    if (patlen > _mpl) {
      _mpl = patlen;
    }
    std::string const & patl = it->pattern().substr(0,patlen/2);
    std::string const & patr = it->pattern().substr(patlen/2);
    id = pm_->add_pattern(patl);
    plit_[id] = it;
    pattern_halves_[id] = patl;
    id = pm_->add_pattern(patr);
    plit_[id] = it;
    pattern_halves_[id] = patr;
    lasthit_[id/2] = 0;
  } 
  pm_->init(cp);
}

void exact_halves::reset() {
  pm_->reset();
}




