
#include "rand_hash_table.h"
// #include "primes.h"
#include "math.h"
#include <map>

rand_hash_table_elt::rand_hash_table_elt() : patid_(0) {};

rand_hash_table_elt::~rand_hash_table_elt() {
  if (patid_) {
    delete patid_;
    patid_ = 0;
  }
}

void rand_hash_table_elt::add_patid(pattern_list::const_iterator const & it,
				    unsigned int p, unsigned int i) {
  // checkpoint;
  // cerr << (void*)patid_ << endl;
  if (!patid_) patid_ = new tinylist<rhtele>;
  // checkpoint;
  patid_->push_front(rhtele(it,p,i));
  // checkpoint;
}

tinylist<rhtele> * const & rand_hash_table_elt::patids() const {
  return patid_;
}

rand_hash_table::rand_hash_table(unsigned int np, unsigned int ws, unsigned int k, 
				 unsigned char eos, bool wc, bool tn, 
				 bool indels, bool dm) {
  ws_ = ws;
  lastnch_.resize(ws_);
  lastnchptr_ = 0;
  relchars_ = new bool[256];
  for (int i=0;i<256;i++) {
    relchars_[i] = false;
  }
  nprime=np;
  max_id_=0;
  _mpl=0;
  k_ = k;
  eos_ = eos;
  _wc = wc;
  _textn = tn;
  _indels = indels;
  _dna_mut = dm;
}

rand_hash_table::~rand_hash_table() {
  delete [] relchars_;
  relchars_ = 0;
  delete [] relcharmap_;
  relcharmap_ = 0;
}

unsigned long rand_hash_table::add_pattern(std::string const & pat, 
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

bigword inchash(bigword h, 
		unsigned long int base, 
		unsigned long int d0, 
		unsigned long int dn, 
		unsigned long int powbasenmodp,
		unsigned long int p) {
  // cerr << h << " " << base << " " << d0 << " " << dn << " " << powbasenmodp << " " << p << endl;
  long unsigned int h1 = ((h*base+d0)%p + p - (dn*powbasenmodp)%p)%p;
  return h1;
}

bigword inchash(bigword h, 
		unsigned long int base, 
		unsigned long int d0, 
		unsigned long int p) {
  // cerr << h << " " << base << " " << d0 << " " << dn << " " << powbasenmodp << " " << p << endl;
  long unsigned int h1 = (h*base+d0)%p;
  return h1;
}

bigword hash(unsigned long int base, 
	     vector<unsigned long int> d, 
	     unsigned long int p) {
  // checkpoint;
  long unsigned h1 = 0;
  // for (int i=d.size()-1;i>=0;i--) {
  for (int i=0;i<d.size();i++) {
    h1 = (h1*base+d[i])%p;
  }
  return h1;
}

extern "C"
{
#include "primegen.h"
}

void rand_hash_table::random_primes_lt(unsigned long int m, std::vector<long unsigned int> &p) {
  primegen pg;
  primegen_init(&pg);
  uint64 np = primegen_count(&pg,(uint64)m);
  uint64 pr;
  
  // checkpoint;
  // cerr << m << " " << np << endl;

  char* ptr;
  if ((ptr=getenv("RAND48_SEED"))) {
    srand48(atol(ptr));
    timestampli("Using random seed:",atol(ptr));
  } else {
    srand48(time(NULL));
  }
  std::vector<long unsigned int> pi(p.size());
  for (int i=0;i<p.size();i++) {
    pi[i] = (long unsigned int) floor(np*DRAND48);
    // checkpoint;
    // cerr << pi[i] << endl;
  }
  sort(pi.begin(),pi.end());
  primegen_init(&pg);
  int j = 0,k = 0;
  while (j<=pi[p.size()-1]) {
    // checkpoint;
    pr = primegen_next(&pg);
    // cerr << j << " " << pr << endl;
    while (j == pi[k] && k < p.size()) {
      p[k] = pr;
      k++;
    }
    j++;
  }

}

void rand_hash_table::init(CharacterProducer & cp) {
  bool *tmp;
  int size=cp.size();
  tmp = new bool[size];
  relcharmap_ = new unsigned char[size];
  int j=0;
  for (int i=0;i<size;i++) {
    if (relchars_[(unsigned char)cp.ch(i)]) {
      tmp[i] = true;
      relcharmap_[i] = j;
      // cerr << i << " " << j << endl;
      j++;
    } else {
      tmp[i] = false;
    }
  }
  delete [] relchars_;
  relchars_ = tmp;

  alphasize_=j;
  // cerr << alphasize_ << endl;

  long unsigned int maxp = (((unsigned long int)1)<<27);
  long unsigned int astows=1,i=0;
  while (i<ws_ && astows < maxp) {
    astows *= alphasize_;
    i++;
  }
  if (astows < maxp) {
    timestamp("Exact hash table fits in memory - no false positives.");
    maxp = astows;
    nprime=1;
    h_.resize(nprime);
    prime_.resize(nprime);
    basetonmodp_.resize(nprime);
    prime_[0]=astows;
    basetonmodp_[0]=0;
  } else {
    h_.resize(nprime);
    prime_.resize(nprime);
    basetonmodp_.resize(nprime);
    random_primes_lt(maxp,prime_);
    for (int i=0;i<nprime;i++) {
      // cerr << i << " " << prime_[i] << endl;
      basetonmodp_[i] = 1;
      for (int j=0;j<ws_;j++) {
	basetonmodp_[i] = (basetonmodp_[i]*alphasize_)%prime_[i];
      }
    }
  }
  table_.resize(prime_[nprime-1]);

  _mpl=0;
  lastpos_.resize(max_id_+1);
  fill(lastpos_.begin(),lastpos_.end(),0);
  tinylist<pattern_list_element>::const_iterator it;
  for (it=patterns().begin();it!=patterns().end();++it) {
    std::string const & pat = it->pattern();
    int nch;
    for (int i=0;i<nprime;i++) {
      h_[i] = 0;
    }
    fill(lastnch_.begin(),lastnch_.end(),0);
    lastnchptr_=0;
    if (pat.length() > _mpl) {
      _mpl = pat.length();
    }
    for (int j=0, p=-ws_+1;j<pat.length();j++,p++) {
      if ((nch=cp.nch(pat[j])) == -1) {
	for (int i=0;i<nprime;i++) {
	  h_[i] = 0;
	}
	fill(lastnch_.begin(),lastnch_.end(),0);
	lastnchptr_=0;
	p = -ws_+1;
	continue;
      }
      for (int i=0;i<nprime;i++) {
	h_[i] = inchash(h_[i],alphasize_,
			relcharmap_[nch],lastnch_[lastnchptr_],
			basetonmodp_[i],prime_[i]);
      }
      if (p >= 0) {
	for (int i=0;i<nprime;i++) {
	  table_[h_[i]].add_patid(it,j,i);
	  // cerr << i << " " << h_[i] << " "  << (void*)table_[h_[i]].patids() << endl;
	}
      }
      lastnch_[lastnchptr_] = relcharmap_[nch];
      lastnchptr_ = (lastnchptr_ + 1)%ws_;
      // int i = 0;
      // while (i<ws_) {
      // cerr << (int)lastnch_[(i+lastnchptr_)%ws_] << " ";
      // i++;
      // }
      // cerr << " " << h_[0] << endl;
    }
  }
  for (int i=0;i<nprime;i++) {
    h_[i] = 0;
  }
  fill(lastnch_.begin(),lastnch_.end(),0);
  lastnchptr_=0;
  p_ = -ws_+1;
}

bool rand_hash_table::find_patterns(CharacterProducer & cp, 
				    pattern_hit_vector & pas,
				    long unsigned minka) {
  // checkpoint;
  // cerr << endl;
  editdist_alignment pa(0,0,k_,eos_,_wc,_textn,_indels,_dna_mut,0,0,true);
  pa.maxpatlen(_mpl);

  long unsigned kacount=0;
  unsigned char nch;
  while (!cp.eof()) {
    nch = cp.getnch();
    // cerr << cp.ch(nch) << endl;
    if (!relchars_[nch]) {
      for (int i=0;i<nprime;i++) {
	h_[i] = 0;
      }
      fill(lastnch_.begin(),lastnch_.end(),0);
      lastnchptr_=0;
      p_ = -ws_;
      continue;
    } else {
      for (int i=0;i<nprime;i++) {
	h_[i] = inchash(h_[i],alphasize_,
			relcharmap_[nch],lastnch_[lastnchptr_],
			basetonmodp_[i],prime_[i]);
      }
    }
    // cerr << binary(h_);
    if (p_ >= 0) {
      // cerr << endl;
      // cerr << h_[0] << " ";
      std::vector<tinylist<rhtele> *> ptids;
      ptids.resize(nprime);
      bool potential_hit = true;
      // checkpoint;
      for (int i=0;i<nprime;i++) {
	if (!(ptids[i]=table_[h_[i]].patids())) {
	  potential_hit = false;
	}
	// cerr << i << ": " << h_[i] << ": " << (void*)ptids[i] << endl;
      }
      if (potential_hit) {
	std::map<pattern_list::const_iterator,int> m0;
	std::map<pattern_list::const_iterator,int>::iterator m0it; 
	tinylist<rhtele>::iterator it=ptids[0]->begin();
	// checkpoint;
	while (it!=ptids[0]->end()) {
	  // cerr << "0: " << it->pattern_list_it()->id() << "," << it->pattern_list_it()->pattern() << "," << it->position() << "," << it->index() <<endl;
	  if (it->index() == 0) {
	    m0[it->pattern_list_it()] = 1;
	  }
	  ++it;
	}	    
	// checkpoint;
	for (int i=1;i<nprime;i++) {
	  tinylist<rhtele>::iterator it=ptids[i]->begin();
	  while (it!=ptids[i]->end()) {
	    m0it = m0.find(it->pattern_list_it());
	    if (m0it != m0.end() && it->index() == i) {
	      m0it->second++;
	      // cerr << i << ": " << it->pattern_list_it()->id() << "," << it->pattern_list_it()->pattern() << "," << it->position() << "," << it->index() << endl;
	    }
	    ++it;
	  }	    
	}
	// checkpoint;
	FILE_POSITION_TYPE oldpos=cp.pos();
	it=ptids[0]->begin();
	while (it!=ptids[0]->end()) {
	  pattern_list::const_iterator const & plit = it->pattern_list_it();
	  m0it = m0.find(it->pattern_list_it());
	  if (m0it == m0.end() || m0it->second != nprime) {
	    ++it;
	    continue;
	  }
	  std::string const & pat = plit->pattern();
	  unsigned int patlen = pat.length();
	  report_progress(cp);
	  /* mismatch_alignment ma(oldpos);
	  // cerr << pat << " " << pat.substr(it->position()-ws_+1,ws_) << " " << it->position() << " " << ws_ << endl;
	  ma.align(cp,pat.substr(it->position()-ws_+1,ws_));
	  if (ma.editdist()>0) {
	    checkpoint;
	    std::string w=ma.matching_text();
	    std::vector<unsigned long int> d(ws_);
	    for (int i=0;i<ws_;i++) {
	      d[i] = relcharmap_[cp.nch(w[i])];
	    }
	    std::string w1=pat.substr(it->position()-ws_+1,ws_);
	    std::vector<unsigned long int> d1(ws_);
	    for (int i=0;i<ws_;i++) {
	      d1[i] = relcharmap_[cp.nch(w1[i])];
	    }
	    cerr << "False positive hit: "
		 << w
		 << " to " 
		 << w1
		 << endl;
	    for (int i=0;i<nprime;i++) {
	      cerr << "hash(" << w << ") = " << hash(alphasize_,d,prime_[i]) 
		   << ", hash(" << w1 << ") = " << hash(alphasize_,d1,prime_[i])
		   << endl;
	    }
	  } else {
	    checkpoint;
	    cerr << "True positive hit: "
		 << ma.matching_text()
		 << " to " 
		 << pat.substr(it->position()-ws_+1,ws_) 
		 << endl;
		 } */
	  FILE_POSITION_TYPE patend = oldpos + patlen-(it->position())-1;
	  // checkpoint;
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
	    // cerr << pat << " " << it->position() << " " << oldpos << endl;
	    if (pa.align(cp,pat)) {
	      // checkpoint;
	      // cerr << pa.end() << endl;
	      if (lastpos_[plit->id()]+((_indels)?k_:0) < pa.end()) {
		lastpos_[plit->id()] = pa.end();
		// cerr << pa.end() << " " << plit->id() << endl;
		pas.push_back(pa.end(),make_pair(plit,pa.editdist()));
		++kacount;
	      } else {
		// checkpoint;
		lastpos_[plit->id()] = patend;
	      }
	    } else {
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
    lastnch_[lastnchptr_] = relcharmap_[nch];
    lastnchptr_ = (lastnchptr_ + 1)%ws_;
  }
  if (kacount>0) {
    report_progress(cp);
    return true;
  } 
  return false;
}

void
rand_hash_table::reset() {
  
}
