
#include "exact_bases.h"
#include "primer_alignment.h"
// #include "sortedvector.t"

exact_bases::exact_bases(PatternMatch *pm,
			 unsigned int k, unsigned char eos, bool wc, bool tn, bool id, bool dm) 
  : pm_(pm), num_patterns_(0), k_(k), eos_(eos), _wc(wc), _textn(tn), _indels(id), _dna_mut(dm)
{
}

exact_bases::exact_bases(PatternMatch *pm)
  : pm_(pm), num_patterns_(0), k_(0), eos_('\n'), _wc(false), _textn(false), _indels(true), _dna_mut(false)
{
}

exact_bases::~exact_bases() {
  delete pm_;
}

unsigned int exact_bases::mismatches() const {
  return k_;
}

void exact_bases::mismatches(unsigned int k) {
  k_ = k;
}

bool exact_bases::wildcards() const {
  return _wc;
}

void exact_bases::wildcards(bool wc) {
  _wc = wc;
}

bool exact_bases::wildcard_text_N() const {
  return _textn;
}

void exact_bases::wildcard_text_N(bool tn) {
  _textn = tn;
}

bool exact_bases::indels() const {
  return _indels;
}

void exact_bases::indels(bool id) {
  _indels = id;
}

unsigned char exact_bases::eos_char() const {
  return eos_;
}

void exact_bases::eos_char(unsigned char c) {
  eos_ = c;
}

long unsigned int 
exact_bases::add_pattern(std::string const & pat, unsigned long id, 
			 int esb, int eeb) {
  add_pattern_(pat,id,esb,eeb);
  num_patterns_++;
  return id;
}

bool 
exact_bases::find_patterns(CharacterProducer & cp, 
			   pattern_hit_vector & phs,
			   long unsigned minka) {
  // checkpoint;
  pattern_hit_vector l;
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
    FILE_POSITION_TYPE oldcharspos;
    oldcharspos = cp.pos();
    it = l.begin();
    // checkpoint;
    while (it != l.end()) {
      // checkpoint;
      long unsigned int pid0(it->value().first->id());
      FILE_POSITION_TYPE pos(it->key());
      pattern_list::const_iterator const & plit(plit_[pid0]);
      long unsigned int pid = plit->id();
      int esb=plit->exact_start_bases();
      int eeb=plit->exact_end_bases();
      if (prefix_[pid]) {
	// Exact match to first part
	// checkpoint;
	pal.reset();
	pal.pos(it->key());
	pal.exact_start_bases(esb);
	pal.exact_end_bases(eeb);
	if (pal.align(cp,it->value().first->pattern(),rempat_[pid])) {
	  phs.push_back(pal.end(),make_pair(plit,pal.value()));
	}
      } else {
	// Exact match to second part
	// checkpoint;
	par.reset();
	par.pos(it->key());
	par.exact_start_bases(esb);
	par.exact_end_bases(eeb);
	if (par.align(cp,rempat_[pid],it->value().first->pattern())) {
	  phs.push_back(par.end(),make_pair(plit,par.value()));
	}
      }
      ++it;
    }
    l.clear();
    cp.pos(oldcharspos);
    report_progress(cp);
    if (phs.size() >= minka || 
	(more==false && phs.size() > 0)) return true;
  }
  return false;
}

void 
exact_bases::init(CharacterProducer & cp) {
  assert(pm_!=((void*)0));
  long unsigned int id=0;
  pattern_list::const_iterator it;
  plit_.resize(num_patterns_+1);
  prefix_.resize(num_patterns_+1);
  rempat_.resize(num_patterns_+1);
  _mpl=0;
  for (it=patterns().begin();it!=patterns().end();++it) {
    int esb=it->exact_start_bases();
    int eeb=it->exact_end_bases();
    if (it->pattern().length() > _mpl) {
      _mpl = it->pattern().length();
    }
    if (esb >= eeb) {
      id = pm_->add_pattern(it->pattern().substr(0,esb));
      prefix_[id] = true;
      rempat_[id] = it->pattern().substr(esb);
    } else {
      int patlen=it->pattern().length();
      id = pm_->add_pattern(it->pattern().substr(patlen-eeb));
      prefix_[id] = false;
      rempat_[id] = it->pattern().substr(0,patlen-eeb);
    }
    plit_[id] = it;
  } 
  pm_->init(cp);
  // checkpoint;
}

void exact_bases::reset() {
  pm_->reset();
}




