
#include "gs_hash_table.h"
#include "rand_hash_table.h"
// #include "primes.h"
#include "math.h"
#include <map>

gs_hash_table_elt::gs_hash_table_elt() : patid_(0) {};

gs_hash_table_elt::~gs_hash_table_elt() {
  clear_patids();
}

void gs_hash_table_elt::add_patid(pattern_list::const_iterator const & it,
				  unsigned int p, unsigned int l, unsigned int i) {
  // checkpoint;
  // cerr << (void*)patid_ << endl;
  if (!patid_) patid_ = new tinylist<gsrhtele>;
  // checkpoint;
  patid_->push_front(gsrhtele(it,p,l,i));
  // checkpoint;
}

void gs_hash_table_elt::clear_patids() {
  if (patid_) {
    delete patid_;
    patid_ = 0;
  }
}


gs_hash_table::gs_hash_table(gapped_seed_set::scheme s, unsigned int np, 
			     unsigned int k,unsigned char eos, bool wc, bool tn, 
			     bool indels, bool dm) {
  gss.initialize(s);
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

gs_hash_table::~gs_hash_table() {
  delete [] relchars_;
  relchars_ = 0;
  delete [] relcharmap_;
  relcharmap_ = 0;
}

unsigned long gs_hash_table::add_pattern(std::string const & pat, 
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

#include <limits.h>

bigword gs_hash_table::hash(vector<unsigned char> const & d, 
			    unsigned long int p) {

  // checkpoint;
  long unsigned h1 = 0;
  // for (int i=d.size()-1;i>=0;i--) {
  int i=0,j=d.size();
  for (;j>1;i++,j--) {
    if (h1 >= must_mod_) {
      h1 = (h1*asize_+d[i])%p;
    } else {
      h1 = h1*asize_+d[i];
    }
  }
  h1 = (h1*asize_+d[i])%p;
  return h1;
} 

/* bigword gs_hash_table::hash(vector<unsigned char> const & d, 
			    unsigned long int p) {

  long unsigned h1 = 0;
  for (int i=0;i<d.size();i++) {
    h1 = (h1*asize_+d[i])%p;
  }
  return h1;
} */

// void random_primes_lt(unsigned long int m, std::vector<long unsigned int> &p);

void gs_hash_table::init(CharacterProducer & cp) {
  // checkpoint;
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

  asize_=j;
  if (asize_ > 0) {
    must_mod_ = (ULONG_MAX - asize_)/asize_;
  }
  long unsigned int maxp = (((unsigned long int)1)<<27);
  long unsigned int astol=1,i=0;
  while (i<gss.l && astol < maxp) {
    astol *= asize_;
    i++;
  }
  if (astol < maxp) {
    // timestamp("Exact gapped seed hash table fits in memory - no false positives.");
    nprime=1;
    prime_.resize(nprime);
    prime_[0]=astol;
  } else {
    prime_.resize(nprime);
    rand_hash_table::random_primes_lt(maxp,prime_);
  }
  table_.resize(prime_[nprime-1]);

  // checkpoint;
  // std::cerr << gss.tmax << endl;

  _mpl=0;
  buffer_.resize(gss.m);
  lastpos_.resize(max_id_+1);
  fill(lastpos_.begin(),lastpos_.end(),0);
  std::vector<unsigned char> nch(gss.l);
  tinylist<pattern_list_element>::const_iterator it;
  for (it=patterns().begin();it!=patterns().end();++it) {
    std::string const & pat = it->pattern();
    // checkpoint;
    // cerr << pat << endl;
    if (pat.length() > _mpl) {
      _mpl = pat.length();
    }
    for (int l=0;l<gss.n;l++) {
      int lastpos;
      if (gss.shift) {
	lastpos = gss.patpos[l][gss.l-1];
      } else {
	lastpos = gss.m-1;
      }
      for (int j=0;j<gss.m-lastpos;j+=1) {
	bool nobad = true;
	for (int i=0;i<gss.l;i++) {
	  int nchi = cp.nch(pat[j+gss.patpos[l][i]]);
	  if (nchi == -1) {
	    nobad = false;
	    break;
	  } 
	  nch[i] = relcharmap_[nchi];
	}
	if (!nobad) {
	  continue;
	}
	/* checkpoint;
	cerr << l << ": ";
	for (int i=0;i<gss.l;i++) {
	  cerr << cp.ch(nch[i]);
	}
	cerr << " ";
	for (int i=0;i<gss.l;i++) {
	  cerr << j + gss.patpos[l][i] << ",";
	  }
	*/
	for (int p=0;p<nprime;p++) {
	  long unsigned int h = hash(nch,prime_[p]);
	  // cerr << h << " ";
	  table_[h].add_patid(it,j,l,p);
	  // checkpoint;
	  /* tinylist<gsrhtele>::iterator it=table_[h].patids()->begin();
	  while (it != table_[h].patids()->end()) {
	    it.write(std::cerr);
	    it->write(std::cerr) << endl;
	    ++it;
          } */
	}
	// cerr << endl;
      }
    }
  }
  bufp_ = 0;
  // checkpoint;
  // std::cerr << gss.tmax << endl;
  buffer_.resize(gss.tmax);
  // checkpoint;
  fill(buffer_.begin(),buffer_.end(),(char)-1);
  // checkpoint;
}

void gs_hash_table::del_pattern(CharacterProducer & cp, std::string const & pat,
				long unsigned int id) {
  std::vector<unsigned char> nch(gss.l);
  for (int l=0;l<gss.n;l++) {
    int lastpos;
    if (gss.shift) {
      lastpos = gss.patpos[l][gss.l-1]; 
    } else {
      lastpos = gss.m;
    }
    for (int j=0;j<pat.length()-lastpos;j+=1) {
      bool nobad = true;
      for (int i=0;i<gss.l;i++) {
	int nchi = cp.nch(pat[j+gss.patpos[l][i]]);
	if (nchi == -1) {
	  nobad = false;
	  break;
	} 
	nch[i] = relcharmap_[nchi];
      }
      if (!nobad) {
	continue;
      }
      for (int p=0;p<nprime;p++) {
	long unsigned int h = hash(nch,prime_[p]);
	// assert(table_[h].patids()!=0);
	if (!table_[h].patids()) {
	  continue;
	}
	tinylist<gsrhtele>::iterator eit=table_[h].patids()->begin();
	tinylist<gsrhtele>::iterator eit0=table_[h].patids()->end();
	while (eit != table_[h].patids()->end()) {
	  if (eit->pattern_list_it()->pattern() == pat &&
	      (id == 0 || eit->pattern_list_it()->id() == id)) {
	    if (eit0==table_[h].patids()->end()) {
	      table_[h].patids()->pop();
	      eit = table_[h].patids()->begin();
	    } else {
	      eit = table_[h].patids()->erase_after(eit0);
	    }
	  } else {
	    eit0 = eit;
	    ++eit;
	  }
	}
	if (table_[h].patids()->empty()) {
	  table_[h].clear_patids();
	}
      }
    }
  }
}

bool gs_hash_table::find_patterns(CharacterProducer & cp, 
				  pattern_hit_vector & pas,
				  long unsigned minka) {
  // checkpoint;
  // cerr << k_ << endl;
  // cerr << eos_ << endl;
  editdist_alignment pa(0,0,k_,eos_,_wc,_textn,_indels,_dna_mut,0,0,true);
  pa.maxpatlen(_mpl);

  int eofcnt=0;
  std::vector<unsigned char> nch(gss.l);
  long unsigned int h_;
  std::vector<tinylist<gsrhtele> *> ptids(nprime);
  long unsigned kacount=0;
  while (1) {
    if (!cp.eof()) {
      buffer_[bufp_] = cp.getnch();
    } else {
      buffer_[bufp_] = -1;
      eofcnt++;
    }
    if (eofcnt == gss.tmax) {
      break;
    }
    bufp_ = (bufp_+1)%gss.tmax;
    // checkpoint;
    /*cerr << bufp_ << ": ";
    for (int i=0;i<buffer_.size();i++) {
    cerr << buffer_[i] << " ";
    }
    cerr << endl;*/
    for (int l=0;l<gss.n;l++) { 
      bool nobad = true;
      for (int i=0;i<gss.l;i++) {
	int nchi = buffer_[(bufp_+gss.txtpos[l][i])%gss.tmax];
	if (nchi == -1 || !relchars_[nchi]) {
	  nobad = false;
	  break;
	}
	nch[i] = relcharmap_[nchi];
      }
      if (!nobad) {
	continue;
      }
      instrument1();
      //checkpoint;
      /*cerr << l << ": ";
      for (int i=0;i<gss.l;i++) {
	cerr << cp.ch(nch[i]);
      }
      cerr << " "; 
      for (int i=0;i<gss.l;i++) {
	cerr << gss.txtpos[l][i] << ",";
      }
      cerr << endl;*/
      bool potential_hit = true;
      for (int p=0;p<nprime;p++) {
	h_ = hash(nch,prime_[p]);
	if (!(ptids[p]=table_[h_].patids())) {
	  // checkpoint;
	  potential_hit = false;
	  break;
	}
      }
      if (potential_hit) {
	instrument2();
	// checkpoint;
	std::map<std::pair<pattern_list::const_iterator,int>,int> m0;
	std::map<std::pair<pattern_list::const_iterator,int>,int>::iterator m0it; 
	tinylist<gsrhtele>::iterator it=ptids[0]->begin();
	// ptids[0]->end().write(std::cerr) << endl;
	bool found=false;
	while (it!=ptids[0]->end()) {
	  // checkpoint;
	  // it.write(std::cerr) << endl;
	  // it->write(std::cerr) << endl;
	  // checkpoint;
	  if (it->templ() == l && it->index() == 0) {
	    m0[make_pair(it->pattern_list_it(),it->position())] = 1;
	    found = true;
	  }
	  // checkpoint;
	  ++it;
	}	    
	if (!found) {
	  continue;
	}
	instrument3();
	// checkpoint;
	for (int i=1;i<nprime;i++) {
	  found = false;
	  tinylist<gsrhtele>::iterator it=ptids[i]->begin();
	  while (it!=ptids[i]->end()) {
	    m0it = m0.find(make_pair(it->pattern_list_it(),it->position()));
	    if (m0it != m0.end() && it->templ() == l && it->index() == i) {
	      found = true;
	      m0it->second++;
	    }
	    ++it;
	  }
	  if (!found) {
	    break;
	  }
	}
	if (!found) {
	  continue;
	}
	instrument4();
	FILE_POSITION_TYPE oldpos=cp.pos();
	// checkpoint;
	// cerr << "oldpos: " << oldpos << endl;
	it=ptids[0]->begin();
	while (it!=ptids[0]->end()) {
	  pattern_list::const_iterator const & plit = it->pattern_list_it();
	  m0it = m0.find(make_pair(it->pattern_list_it(),it->position()));
	  if (m0it == m0.end()) {
	    ++it;
	    continue;
	  }
	  instrument5();
	  assert(m0it->second <= nprime);
	  if (m0it->second != nprime) {
	    ++it;
	    continue;
	  }
	  instrument6();
	  std::string const & pat = plit->pattern();
	  unsigned int patlen = pat.length();
	  /* 
	  report_progress(cp);
	  mismatch_alignment ma(oldpos);
	  cerr << pat << " " << pat.substr(it->position(),gss.m) << " " << it->position() << " " << gss.m << endl;
	  ma.align(cp,pat.substr(it->position(),gss.m));
	  checkpoint;
	  std::string w=ma.matching_text();
	  std::vector<unsigned long int> d(gss.l);
	  for (int i=0;i<gss.l;i++) {
	    d[i] = relcharmap_[cp.nch(w[gss.txtpos[l][i]])];
            w[gss.txtpos[l][i]] = w[gss.txtpos[l][i]]-'A'+'a';
	  }
	  std::string w1=pat.substr(it->position(),gss.m);
	  std::vector<unsigned long int> d1(gss.l);
	  for (int i=0;i<gss.l;i++) {
	    d1[i] = relcharmap_[cp.nch(w1[gss.patpos[l][i]])];
            w1[gss.patpos[l][i]] = w1[gss.patpos[l][i]]-'A'+'a';
	  }
	  cerr << "Hit: "
		 << w
		 << " to " 
		 << w1
		 << endl;
	  for (int i=0;i<nprime;i++) {
	    cerr << "hash(" << w << ") = " << hash(asize_,d,prime_[i]) 
	         << ", hash(" << w1 << ") = " << hash(asize_,d1,prime_[i])
		 << endl;
		 } */
	  // FILE_POSITION_TYPE patend = oldpos + patlen-(it->position())-1;
	  FILE_POSITION_TYPE patend = patlen;
	  if (oldpos > patlen) {
	    patend = oldpos + patlen - gss.tmax - it->position();
	  } 
	  // checkpoint; 
	  // cerr << "plit->id(): " << plit->id() << endl;
	  // cerr << "lastpos: " << lastpos_[plit->id()]+((_indels)?k_:0) << endl;
          // cerr << "patend: " << patend << endl;
	  if (lastpos_[plit->id()]+((_indels)?k_:0) < patend) {
	    // checkpoint;
            // cerr << "patend: " << patend << endl;
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
	    report_progress(cp);
	    instrument8();
	    if (pa.align(cp,pat)) {
	      // checkpoint;
	      // cerr << pat << endl;
	      // cerr << pa.alignment_string() << endl;
	      // cerr << pa.matching_text();
	      // cerr << " " << pa.end() << " " << endl;
	      if (lastpos_[plit->id()]+((_indels)?k_:0) < pa.end()) {
		instrument9();
		lastpos_[plit->id()] = pa.end();
		// cerr << pa.end() << " " << plit->id() << endl;
		pas.push_back(pa.end(),make_pair(plit,pa.value()));
		++kacount;
	      } else {
		// checkpoint;
		instrument10();
		lastpos_[plit->id()] = patend;
	      }
	    } else {
	      instrument11();
	      lastpos_[plit->id()] = patend;
	    }
	  }
	  ++it;
	}
	cp.pos(oldpos);
      } else {
	// cerr << " skip..." << endl;
      }
    }
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
gs_hash_table::reset() {
  
}

void gapped_seed_set::init(int *pos,int *tpos) {
  tmax = 0;
  txtpos.resize(n);
  patpos.resize(n);
  for (int i=0;i<n;i++) {
    txtpos[i].resize(l);
    patpos[i].resize(l);
    for (int j=0;j<l;j++) {
      patpos[i][j] = pos[i*l+j];
      if (!tpos) {
	txtpos[i][j] = pos[i*l+j];
      } else {
	txtpos[i][j] = tpos[i*l+j];	
      }
      if (txtpos[i][j] >= tmax) {
	tmax = txtpos[i][j]+1;
      }
    }
  }
}

std::string gapped_seed_set::str() const {
  ostrstream ss;
  ss << "m" << m << (indels?"k":"K") << k << "l" << l << (shift?"S":"") << "(n" << n << ")" << std::ends;
  return string(ss.str());
}

gapped_seed_set::scheme
gapped_seed_set::select(int m, int k, bool indels) {
  if (!indels) {
    if (k <=1) {
      if (m >= 8) {
	return m8K1l6S;	
      } else if (m >= 7) {
	return m7K1l5S;	
      } else if (m >= 6) {
	return m6K1l4S;
      } else if (m >= 5) {
	return m5K1l3S;
      }
    } else if (k <=2) {
      if (m >= 30) {
	return m30K2l16S;
      } else if (m >= 28) {
	return m28K2l16S;	
      } else if (m >= 20) {
	return m20K2l12S;
      } else if (m >= 19) {
	return m19K2l12S;
      } else if (m >= 10) {
	return m10K2l5S;
      } else if (m >= 9) {
	return m9K2l5S;
      }
    } else if (k <= 3) {
      if (m >= 20) {
	return m20K3l9S;
      } else if (m >= 19) {
	return m19K3l9S;
      } 
    }
  } else {
    if (k <= 2) {
      if (m >= 20) {
	return m20k2l10S;
      } else if (m >= 19) {
	return m19k2l10S;
      }
    } else if (k <= 3) {
      if (m >= 20) {
	return m20k3l7S;
      } else if (m >= 19) {
	return m19k3l7S;
      } 
    }
  }
  /*  checkpoint;
  cerr << "Returning no_scheme for m = " << m 
       << ", k = " << k 
       << ", indels = " << ((indels)?"TRUE":"FALSE")
       << endl; */
  return no_scheme;
}

bool gapped_seed_set::initialize(scheme s) {
  switch (s) {

  case m5K1l3S:
    {
      int pos[] = {
	0, 1, 2,
	0, 1, 3,
      };
      m = 5;
      l = 3;
      n = 2;
      k = 1; indels = false;
      shift = true;
      init(pos);
    }
    break;

  case m6K1l4S:
    {
      int pos[] = {
	0, 1, 2, 4,
	0, 2, 3, 4,
      };
      m = 6;
      l = 4;
      n = 2;
      k = 1; indels = false;
      shift = true;
      init(pos);
    }
    break;

  case m7K1l5S:
    {
      int pos[] = {
	0, 1, 2, 3, 4,
	0, 1, 4, 5, 6,
	0, 1, 2, 3, 5,
      };
      m = 7;
      l = 5;
      n = 3;
      k = 1; indels = false;
      shift = true;
      init(pos);
    }
    break;

  case m8K1l6S:
    {
      int pos[] = {
	0, 1, 2, 3, 4, 5,
	0, 1, 4, 5, 6, 7,
	0, 1, 2, 3, 6, 7,
      };
      m = 8;
      l = 6;
      n = 3;
      k = 1; indels = false;
      shift = true;
      init(pos);
    }
    break;

  case m9K2l5S:
    {
      int pos[] = {
	0, 1, 3, 4, 5,
	0, 1, 2, 3, 6,
	0, 1, 5, 7, 8,
	0, 2, 4, 5, 6,
      };
      m = 9;
      l = 5;
      n = 4;
      k = 2; indels = false;
      shift = true;
      init(pos);
    }
    break;

  case m10K2l5S:
    {
      int pos[] = {
	0, 1, 3, 4, 5,
	0, 1, 2, 6, 7,
	0, 1, 5, 6, 8,
      };
      m = 10;
      l = 5;
      n = 3;
      k = 2; indels = false;
      shift = true;
      init(pos);
    }
    break;

  case m20K2l12S:
    {
      int pos[] = {
	0, 1, 2, 3, 4, 6, 7, 8, 10, 11, 12, 13, 
	0, 1, 2, 3, 4, 5, 6, 8, 9, 12, 13, 14, 
	0, 1, 2, 3, 5, 6, 7, 8, 10, 13, 14, 15, 
	0, 1, 2, 4, 8, 9, 10, 11, 13, 14, 15, 16 
      };
      m = 20;
      l = 12;
      n = 4;
      k = 2; indels = false;
      shift = true;
      init(pos);
    }
    break;
  case m19K2l10S:
    {
      int pos[] = {
	0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 
	0, 1, 2, 3, 5, 6, 8, 9, 10, 11, 
	0, 1, 2, 3, 5, 9, 10, 12, 13, 14 
      };
      m = 19;
      l = 10;
      n = 3;
      k = 2; indels = false;
      shift = true;
      init(pos);
    }
    break;
  case m19K2l11S:
    {
      int pos[] = {
	0, 1, 2, 3, 5, 6, 7, 9, 10, 11, 12,
	0, 1, 2, 4, 5, 6, 7, 9, 12, 13, 14, 
	0, 1, 2, 4, 7, 8, 9, 10, 12, 13, 14
      };
      m = 19;
      l = 11;
      n = 3;
      k = 2; indels = false;
      shift = true;
      init(pos);
    }
    break;
  case m19K2l12S:
    {
      int pos[] = {
	0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 13, 
	0, 1, 2, 3, 4, 6, 7, 8, 10, 11, 12, 13, 
	0, 1, 2, 4, 5, 6, 7, 9, 10, 12, 13, 14, 
	0, 1, 3, 6, 7, 8, 9, 10, 12, 13, 14, 15, 
	0, 1, 2, 3, 4, 5, 7, 10, 13, 14, 15, 16
      };
      m = 19;
      l = 12;
      n = 5;
      k = 2; indels = false;
      shift = true;
      init(pos);
    }
    break;
  case m19K2l13S:
    {
      int pos[] = {
	0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 
	0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 12, 13, 14, 
	0, 1, 2, 3, 6, 7, 8, 9, 10, 11, 12, 14, 15, 
	0, 1, 2, 4, 6, 7, 8, 10, 11, 12, 13, 14, 15, 
	0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 16, 
	0, 1, 2, 4, 5, 6, 7, 9, 12, 13, 14, 15, 16, 
	0, 1, 2, 3, 4, 6, 7, 10, 11, 12, 14, 16, 17, 
	0, 1, 2, 4, 5, 8, 10, 11, 12, 14, 15, 16, 17
      };
      m = 19;
      l = 13;
      n = 8;
      k = 2; indels = false;
      shift = true;
      init(pos);
    }
    break;
  case m19K2l14S:
    {
      int pos[] = {
	0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 
	0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 
	0, 1, 2, 3, 4, 6, 7, 9, 10, 11, 12, 13, 14, 15, 
	0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 13, 15, 16, 
	0, 1, 3, 4, 5, 7, 8, 9, 11, 12, 13, 14, 15, 16, 
	0, 1, 2, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 
	0, 1, 3, 4, 5, 6, 7, 9, 10, 11, 12, 14, 15, 17, 
	0, 1, 2, 3, 5, 6, 9, 10, 11, 12, 13, 14, 15, 17, 
	0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 13, 14, 16, 17, 
	0, 1, 3, 4, 6, 8, 9, 10, 11, 12, 13, 14, 16, 17, 
	0, 1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 15, 16, 17, 
	0, 1, 2, 3, 4, 7, 8, 9, 11, 12, 14, 15, 16, 17, 
	0, 1, 2, 3, 5, 6, 8, 10, 11, 12, 14, 15, 16, 17, 
	0, 1, 3, 4, 5, 6, 7, 8, 10, 13, 14, 15, 16, 17 
      };
      m = 19;
      l = 14;
      n = 14;
      k = 2; indels = false;
      shift = true;
      init(pos);
    }
    break;

  case m30K2l16S:
    {
      int pos[] = {
	0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
	0,1,2,3,4,5,8,9,10,11,14,15,16,17,18,19,
	0,1,2,3,4,5,8,9,14,15,20,21,22,23,24,25
      };
      m = 30;
      l = 16;
      n = 3;
      k = 2; indels = false;
      shift = true;
      init(pos);
    }
    break;
  case m28K2l16S:
    {
      int pos[] = {
	0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,
	0,1,2,3,6,7,10,11,12,13,14,15,16,17,18,19,
	0,1,2,3,4,5,10,11,12,13,16,17,20,21,22,23
      };
      m = 28;
      l = 16;
      n = 3;
      k = 2; indels = false;
      shift = true;
      init(pos);
    }
    break;
  case m19K3l9S: 
    {
      int pos[] = {
	0, 1, 2, 3, 5, 6, 7, 8, 9,
	0, 1, 2, 4, 8, 10, 11, 12, 13,
	0, 1, 3, 4, 5, 8, 10, 13, 14,
	0, 1, 2, 5, 7, 10, 11, 14, 16,
	0, 1, 5, 8, 9, 10, 15, 16, 17
      };
      m = 19;
      l = 9;
      n = 5;
      k = 3; indels = false;
      shift = true;
      init(pos);
    }
    break;
  case m19K3l10S:
    {
      int pos[] = {
	0, 1, 2, 3, 4, 5, 6, 7, 8, 10,
	0, 1, 2, 3, 5, 6, 8, 9, 10, 11,
	0, 1, 2, 3, 5, 7, 8, 10, 11, 12,
	0, 1, 2, 3, 5, 6, 7, 9, 11, 13,
	0, 1, 2, 4, 6, 7, 8, 11, 12, 13,
	0, 1, 3, 7, 8, 9, 10, 11, 12, 13,
	0, 1, 2, 4, 5, 9, 11, 13, 14, 15,
	0, 1, 2, 3, 7, 10, 11, 13, 14, 15,
	0, 1, 3, 6, 8, 9, 13, 15, 16, 17
      };
      m = 19;
      l = 10;
      n = 9;
      k = 3; indels = false;
      shift = true;
      init(pos);
    }
    break;
  case m20K3l9S:
    {
      int pos[] = {
	0, 1, 2, 3, 5, 6, 7, 8, 9,
	0, 1, 2, 3, 5, 9, 12, 13, 14,
	0, 1, 3, 4, 8, 9, 11, 13, 14,
	0, 1, 2, 5, 7, 10, 11, 12, 13,
      };
      m = 20;
      l = 9;
      n = 4;
      k = 3; indels = false;
      shift = true;
      init(pos);
    }
    break;
  case m20K3l10S:
    {
      int pos[] = {
	0, 1, 2, 3, 5, 6, 7, 9, 10, 11,
	0, 1, 2, 3, 5, 6, 9, 12, 13, 14,
	0, 2, 3, 8, 10, 11, 12, 13, 14, 15,
	0, 1, 3, 4, 6, 7, 8, 12, 13, 14,
	0, 2, 4, 6, 8, 9, 11, 14, 16, 17,
	0, 1, 2, 3, 6, 10, 12, 15, 16, 17,
	0, 1, 5, 6, 9, 10, 12, 13, 15, 16,
      };
      m = 20;
      l = 10;
      n = 7;
      k = 3; indels = false;
      shift = true;
      init(pos);
    }
    break;
  case m20k2l10S:
    {
      int pos[] = {
	0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
	0, 1, 2, 3, 5, 6, 8, 9, 10, 11,
	0, 1, 2, 3, 4, 8, 9, 10, 11, 12,
	0, 1, 2, 3, 4, 8, 9, 10, 11, 12,
	0, 2, 3, 4, 14, 15, 16, 17, 18, 19,
	0, 1, 2, 3, 4, 15, 16, 17, 18, 19,
	0, 1, 2, 3, 4, 13, 14, 15, 16, 17,
	0, 1, 2, 3, 4, 5, 7, 8, 9, 10,
	0, 1, 2, 3, 4, 6, 7, 8, 9, 10
      };
      int tpos[] = {
	0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
	0, 1, 2, 3, 5, 6, 8, 9, 10, 11,
	0, 1, 2, 3, 4, 7, 8, 9, 10, 11,
	0, 1, 2, 3, 4, 9, 10, 11, 12, 13,
	0, 2, 3, 4, 14, 15, 16, 17, 18, 19,
	0, 1, 2, 3, 4, 14, 15, 16, 17, 18,
	0, 1, 2, 3, 4, 11, 12, 13, 14, 15,
	0, 1, 2, 4, 5, 6, 8, 9, 10, 11,
	0, 1, 2, 3, 4, 6, 7, 8, 10, 11,
      };
      m = 20;
      l = 10;
      n = 9;
      k = 2; indels = true;
      shift = true;
      init(pos,tpos);
    }
    break;
  case m20k2l9S:
    {
      int pos[] = {
	0, 1, 2, 3, 4, 5, 6, 7, 8,
	0, 1, 2, 4, 5, 7, 8, 9, 10,
	0, 1, 2, 3, 4, 6, 7, 8, 9,
	0, 1, 3, 4, 5, 7, 8, 9, 10,
	0, 1, 3, 4, 5, 6, 7, 8, 9,
	0, 1, 2, 3, 4, 5, 6, 7, 9,
      };
      int tpos[] = {
	0, 1, 2, 3, 4, 5, 6, 7, 8,
	0, 1, 2, 4, 5, 7, 8, 9, 10,
	0, 1, 2, 3, 4, 5, 6, 7, 8,
	0, 1, 3, 4, 5, 8, 9, 10, 11,
	0, 1, 2, 3, 4, 5, 6, 7, 8,
	0, 1, 2, 3, 4, 5, 6, 7, 8,
      };
      m = 20;
      l = 9;
      n = 6;
      k = 2; indels = true;
      shift = true;
      init(pos,tpos);
    }
    break;
  case m19k2l9S:
    {
      int pos[] = {
	0, 1, 2, 3, 4, 5, 6, 7, 8,
	0, 1, 2, 4, 5, 7, 8, 9, 10,
	0, 1, 2, 3, 4, 6, 7, 8, 9,
	0, 1, 2, 3, 4, 5, 6, 7, 8,
	0, 1, 2, 3, 13, 15, 16, 17, 18,
	0, 1, 2, 3, 4, 5, 6, 8, 9,
	0, 1, 2, 3, 12, 15, 16, 17, 18,
	0, 1, 2, 3, 4, 5, 6, 7, 18,
	0, 1, 2, 3, 13, 14, 15, 16, 17,
      };
      int tpos[] = {
	0, 1, 2, 3, 4, 5, 6, 7, 8,
	0, 1, 2, 4, 5, 7, 8, 9, 10,
	0, 1, 2, 3, 4, 5, 6, 7, 8,
	0, 1, 2, 3, 4, 5, 7, 8, 9,
	0, 1, 2, 3, 12, 14, 15, 16, 17,
	0, 1, 2, 3, 4, 5, 6, 7, 8,
	0, 1, 2, 3, 13, 16, 17, 18, 19,
	0, 1, 2, 3, 4, 5, 6, 7, 18,
	0, 1, 2, 3, 11, 12, 13, 14, 15,
      };
      m = 19;
      l = 9;
      n = 9;
      k = 2; indels = true;
      shift = true;
      init(pos,tpos);
    }
    break;
  case m19k2l10S:
    {
      int pos[] = {
	0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
	0, 1, 2, 3, 5, 6, 8, 9, 10, 11,
	0, 1, 2, 3, 4, 7, 8, 9, 10, 11,
	0, 1, 2, 3, 4, 8, 9, 10, 11, 12,
	0, 1, 2, 3, 4, 13, 14, 15, 16, 17,
	0, 1, 2, 3, 4, 14, 15, 16, 17, 18,
	0, 1, 2, 3, 4, 12, 13, 14, 15, 16,
	0, 1, 2, 3, 4, 14, 15, 16, 17, 18,
	0, 1, 2, 3, 13, 14, 15, 16, 17, 18,
	0, 1, 2, 4, 5, 6, 7, 9, 11, 12,
	0, 1, 2, 3, 11, 12, 13, 14, 15, 16,
      };
      int tpos[] = {
	0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
	0, 1, 2, 3, 5, 6, 8, 9, 10, 11,
	0, 1, 2, 3, 4, 6, 7, 8, 9, 10,
	0, 1, 2, 3, 4, 9, 10, 11, 12, 13,
	0, 1, 2, 3, 4, 13, 14, 15, 16, 17,
	0, 1, 2, 3, 4, 13, 14, 15, 16, 17,
	0, 1, 2, 3, 4, 10, 11, 12, 13, 14,
	0, 1, 2, 3, 4, 15, 16, 17, 18, 19,
	0, 1, 2, 3, 15, 16, 17, 18, 19, 20,
	0, 1, 2, 4, 5, 6, 7, 8, 10, 11,
	0, 1, 2, 3, 9, 10, 11, 12, 13, 14,
      };
      m = 19;
      l = 10;
      n = 11;
      k = 2; indels = true;
      shift = true;
      init(pos,tpos);
    }
    break;

  case m19k3l6S:
    {
      int pos[] = {
	0,  1,  2,  3,  4,  5,
	0,  1,  2,  4,  6,  7,
	0,  1,  2,  3,  4,  6,
	0,  1,  2,  3,  4,  5,
	0,  1,  2,  3,  5, 10,
	0,  1,  2,  3,  5,  6,
      };
      int tpos[] = {
	0,  1,  2,  3,  4,  5,
	0,  1,  2,  4,  6,  7,
	0,  1,  2,  3,  4,  5,
	0,  1,  2,  3,  5,  6,
	0,  1,  2,  3,  4,  8,
	0,  1,  2,  3,  4,  5,
      };
      m = 19;
      l = 6;
      n = 6;
      k = 3; indels = true;
      shift = true;
      init(pos,tpos);
    }
    break;
  case m19k3l7S:
    {
      int pos[] = {
	0,  1,  2,  3,  4,  5,  6,
	0,  1,  2,  3,  4,  5,  6,
	0,  1,  2,  3,  5,  6,  8,
	0,  1,  2,  3,  4,  6,  8,
	0,  1,  2,  3,  4,  5,  7,
	0,  1,  2,  3,  5,  6,  8,
	0,  1,  2,  3,  5,  6,  7,
	0,  1,  2,  3,  4,  5,  7,
	0,  1,  2,  3,  9, 10, 11,
	0,  1,  2,  3,  4,  6,  7,
	0,  1,  2,  3,  5,  6,  7,
      };
      int tpos[] = {
	0,  1,  2,  3,  4,  5,  6,
	0,  1,  2,  3,  4,  6,  7,
	0,  1,  2,  3,  5,  6,  8,
	0,  1,  2,  3,  4,  5,  7,
	0,  1,  2,  3,  4,  5,  7,
	0,  1,  2,  3,  4,  5,  7,
	0,  1,  2,  3,  5,  6,  7,
	0,  1,  2,  3,  4,  5,  6,
	0,  1,  2,  3,  8,  9, 10,
	0,  1,  2,  3,  4,  5,  6,
	0,  1,  2,  3,  4,  5,  6,
      };
      m = 19;
      l = 7;
      n = 11;
      k = 3; indels = true;
      shift = true;
      init(pos,tpos);
    }
    break;
  case m20k3l6S:
    {
      int pos[] = {
	0,  1,  2,  3,  4,  5,
	0,  1,  2,  3,  5,  6,
	0,  1,  2,  3,  5,  6,
	0,  1,  2,  3,  4,  6,
      };
      int tpos[] = {
	0,  1,  2,  3,  4,  5,
	0,  1,  2,  3,  4,  5,
	0,  1,  2,  3,  5,  6,
	0,  1,  2,  3,  4,  7,
      };
      m = 20;
      l = 6;
      n = 4;
      k = 3; indels = true;
      shift = true;
      init(pos,tpos);
    }
    break;
  case m20k3l7S:
    {
      int pos[] = {
	0,  1,  2,  3,  4,  5,  6,
	0,  1,  2,  3,  4,  5,  6,
	0,  1,  2,  3,  4,  5,  7,
	0,  1,  2,  3,  4,  7,  8,
	0,  1,  2,  3,  5,  7,  8,
	0,  1,  2,  3,  5,  6,  7,
	0,  1,  2,  3,  4,  6,  7,
      };
      int tpos[] = {
	0,  1,  2,  3,  4,  5,  6,
	0,  1,  2,  3,  4,  6,  7,
	0,  1,  2,  3,  4,  5,  7,
	0,  1,  2,  3,  4,  6,  7,
	0,  1,  2,  3,  5,  7,  8,
	0,  1,  2,  3,  4,  5,  6,
	0,  1,  2,  3,  4,  6,  7,
      };
      m = 20;
      l = 7;
      n = 7;
      k = 3; indels = true;
      shift = true;
      init(pos,tpos);
    }
    break;

  default:
    return false;
  }

  return true;
  
}

