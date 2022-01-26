
#include "hash_table.h"

hash_table_elt::hash_table_elt() : patid_(0) {};

hash_table_elt::~hash_table_elt() {
  if (patid_) {
    delete patid_;
    patid_ = 0;
  }
}

void hash_table_elt::add_patid(pattern_list::const_iterator const & it,
			       unsigned int p) {
  // checkpoint;
  // cerr << (void*)patid_ << endl;
  if (!patid_) patid_ = new tinylist<htele>;
  // checkpoint;
  patid_->push_front(htele(it,p));
  // checkpoint;
}

tinylist<htele> * const & hash_table_elt::patids() const {
  return patid_;
}

hash_table::hash_table(unsigned int ws, unsigned int k, 
		       unsigned char eos, bool wc, bool tn, 
		       bool indels, bool dm) {
  ws_ = ws;
  relchars_ = new bool[256];
  for (int i=0;i<256;i++) {
    relchars_[i] = false;
  }
  max_id_=0;
  _mpl=0;
  k_ = k;
  eos_ = eos;
  _wc = wc;
  _textn = tn;
  _indels = indels;
  _dna_mut = dm;
}

hash_table::~hash_table() {
  delete [] relchars_;
  relchars_ = 0;
  delete [] relcharmap_;
  relcharmap_ = 0;
}

unsigned long hash_table::add_pattern(std::string const & pat, 
				      long unsigned id,
				      int esb, int eeb) {
  add_pattern_(pat,id,esb,eeb);
  for (int i=0;i<pat.length();i++) {
    relchars_[pat[i]] = true;
  }
  if (id > max_id_) {
    max_id_ = id;
  }
  return id;
}

void hash_table::init(CharacterProducer & cp) {

  bool *tmp;
  int size=cp.size();
  // cerr << size << endl;
  tmp = new bool[size];
  relcharmap_ = new unsigned char[size];
  int j=0;
  for (int i=0;i<size;i++) {
    if (relchars_[(unsigned char)cp.ch(i)]) {
      tmp[i] = true;
      relcharmap_[i] = j;
      // cerr << i << " " << cp.ch(i) << ": " 
      // << ((int)cp.ch(i)) << ": " << j << endl;
      j++;
    } else {
      tmp[i] = false;
    }
  }
  delete [] relchars_;
  relchars_ = tmp;

  alphasize_=j;
  // cerr << alphasize_ << endl;
  int alphausize=1;
  alphalog_ = 0;
  while (alphasize_ > alphausize) {
    alphausize <<= ((unsigned int)1);
    alphalog_++;
  }
  
  // cerr << ws_ << " " << alphalog_ << " " << sizeof(bigword)*8 << endl;

  assert(ws_ * alphalog_ <= sizeof(bigword)*8);
  bigword tablesize = (1 << ((bigword)(alphalog_*ws_)));
  wsmask_ = (tablesize-1);
  table_.resize(tablesize);
  // cerr << tablesize << endl;

  // checkpoint;
  // cerr << binary(wsmask_) << endl;

  _mpl=0;
  lastpos_.resize(max_id_+1);
  fill(lastpos_.begin(),lastpos_.end(),0);
  tinylist<pattern_list_element>::const_iterator it;
  for (it=patterns().begin();it!=patterns().end();++it) {
    std::string const & pat = it->pattern();
    int nch;
    h_ = 0;
    if (pat.length() > _mpl) {
      _mpl = pat.length();
    }
    // checkpoint;
    // cerr << binary(h_) << endl;
    for (int j=0, p=-ws_+1;j<pat.length();j++,p++) {
      // cerr << j << " " << p << " " << pat[j] << endl;
      if ((nch=cp.nch(pat[j])) == -1) {
	// checkpoint;
	p = -ws_;
	nch = 0;
      }
      // cerr << binary(h_) << endl;
      h_ = ((h_ << alphalog_) | relcharmap_[nch]) & wsmask_;
      // cerr << binary(h_) << endl;
      if (p >= 0) {
	// checkpoint;
	// cerr << h_ << " " << binary(h_) << ": " << j << endl;
	// cerr << h_ << " " << j << " " << pat.substr(j+1-ws_,ws_) << endl;
	table_[h_].add_patid(it,j);
	// checkpoint;
      }
    }
  }
  // checkpoint;
  h_ = 0;
  p_ = -ws_+1;
}

bool hash_table::find_patterns(CharacterProducer & cp, 
			       pattern_hit_vector & pas,
			       long unsigned minka) {
  // checkpoint;
  editdist_alignment pa(0,0,k_,eos_,_wc,_textn,_indels,_dna_mut,0,0,true);
  pa.maxpatlen(_mpl);

  long unsigned kacount=0;
  unsigned char nch;
  while (!cp.eof()) {
    nch = cp.getnch();
    // cerr << cp.ch(nch) << endl;
    if (!relchars_[nch]) {
      p_ = -ws_+1;
      continue;
    } else {
      h_ = ((h_ << alphalog_) | relcharmap_[nch]) & wsmask_;
    }
    // cerr << binary(h_);
    if (p_ >= 0) {
      // cerr << endl;
      tinylist<htele> * ptids;
      if ((ptids=table_[h_].patids())) {
	FILE_POSITION_TYPE oldpos=cp.pos();
	tinylist<htele>::iterator it=ptids->begin();
	while (it!=ptids->end()) {
	  pattern_list::const_iterator const & plit = it->pattern_list_it();
	  if (k_==0) {
	    pas.push_back(cp.pos(),make_pair(plit,0));
	    ++it;
	    continue;
	  }
	  std::string const & pat = plit->pattern();
	  unsigned int patlen = pat.length();
	  FILE_POSITION_TYPE patend = oldpos + patlen-(it->position())-1;
	  if (lastpos_[plit->id()]+((_indels)?k_:0) < patend) {
	    // checkpoint;
	    // cerr << "lastpos: " << lastpos_[plit->id()] << endl;
	    unsigned int esb = plit->exact_start_bases();
	    unsigned int eeb = plit->exact_end_bases();
	    pa.reset();
	    pa.poslb(patend-((_indels)?k_:0));
	    pa.posub(patend+((_indels)?k_:0));
	    pa.exact_start_bases(esb);
	    pa.exact_end_bases(eeb);
	    // checkpoint;
	    /// cerr << pat << " " << it->position() << " " << oldpos << endl;
	    if (pa.align(cp,pat)) {
	      // checkpoint;
	      // cerr << pa.end() << endl;
	      if (lastpos_[plit->id()]+((_indels)?k_:0) < pa.end()) {
		lastpos_[plit->id()] = pa.end();
		// cerr << pa.end() << " " << plit->id() << endl;
		pas.push_back(pa.end(),make_pair(plit,pa.value()));
		++kacount;
	      } else {
		// checkpoint;
		lastpos_[plit->id()] = patend;
	      }
	    } else {
	      // checkpoint;
	      lastpos_[plit->id()] = patend;
	    }
	  }
	  ++it;
	}
	cp.pos(oldpos);
      }
    } else {
      // cerr << " skip..." << endl;
    }
    if (p_ < 0) p_++;
    if (kacount>minka) {
      report_progress(cp);
      return true;
    }
  }
  if (kacount>0) {
    report_progress(cp);
    return true;
  } 
  return false;
}

void
hash_table::reset() {
  
}
