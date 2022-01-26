
#include "filter_bitvec.h"
#include "primer_alignment.h"
// #include "sortedvector.t"

filter_bitvec::filter_bitvec(unsigned int k, char eos, 
			     bool wc, bool tn, bool id, bool dm) :
  num_patterns_(0)
{
  pm_ = new shift_and_inexact(k,eos,wc,tn,id,dm);
}

filter_bitvec::~filter_bitvec() {
  delete pm_;
}

unsigned int filter_bitvec::mismatches() const {
  return pm_->mismatches();
}

void filter_bitvec::mismatches(unsigned int k) {
  pm_->mismatches(k);
}

bool filter_bitvec::wildcards() const {
  return pm_->wildcards();
}

void filter_bitvec::wildcards(bool wc) {
  pm_->wildcards(wc);
}

bool filter_bitvec::wildcard_text_N() const {
  return pm_->wildcard_text_N();
}

void filter_bitvec::wildcard_text_N(bool tn) {
  pm_->wildcard_text_N(tn);
}

bool filter_bitvec::indels() const {
  return pm_->indels();
}

void filter_bitvec::indels(bool id) {
  pm_->indels(id);
}

bool filter_bitvec::dna_mut() const {
  return pm_->dna_mut();
}

void filter_bitvec::dna_mut(bool dm) {
  pm_->dna_mut(dm);
}

char filter_bitvec::eos_char() const {
  return pm_->eos_char();
}

void filter_bitvec::eos_char(char c) {
  pm_->eos_char(c);
}

long unsigned int 
filter_bitvec::add_pattern(std::string const & pat, unsigned long id,
			  int esb, int eeb) {
  add_pattern_(pat,id,esb,eeb);
  num_patterns_++;
  return id;
}

bool 
filter_bitvec::find_patterns(CharacterProducer & cp, 
			     pattern_hit_vector & phs,
			     long unsigned minka) {
  // checkpoint;

  // the eos_char we store in the shift_and_inexact is normalized, so
  // we much unnormalize it before passing it on to the
  // editdist_alignment class.
  editdist_alignment pa(0,0,mismatches(),cp.ch(eos_char()),
			wildcards(),wildcard_text_N(),
			indels(),dna_mut(),0,0,true);
  pattern_hit_vector::iterator it,it1,it2;
  l.reserve(minka*2);
  bool more;
  while ((more=pm_->find_patterns(cp,l,minka))||!l.empty()) {
    // checkpoint;
    FILE_POSITION_TYPE oldcharspos;
    oldcharspos = cp.pos();
    l.normalize();
    // checkpoint;
    it = l.begin();
    while (it != l.end()) {
      // checkpoint;
      FILE_POSITION_TYPE firstpos = it->key();
      if (firstpos > 0) {
	long unsigned int pid = it->value().first->id();
	FILE_POSITION_TYPE pos = firstpos;
	// checkpoint;
	// cerr << pid << " " << firstpos << endl;
	tinylist<pattern_hit_vector::iterator> adj;
	adj.push_front(it);
	it1 = it+1;
	while (it1 != l.end() && it1->key() <= pos+2*mismatches()+1) {
	  // checkpoint;
	  if (it1->value().first->id() == pid) {
	    // checkpoint;
	    // cerr << it1->value()->id() << " " << it1->key() << endl;
	    pos = it1->key();
	    adj.push_front(it1);
	  }
	  ++it1;
	  // checkpoint;
	}
	// checkpoint;
	if (oldcharspos < pos+2*mismatches()+1 && more) {
	  // end of returned hit condition...
	  // checkpoint;
	  break;
	} else {
	  // checkpoint;
	  // cerr << pid << " " << firstpos << " " << pos << endl;
	  pattern_list::const_iterator plit = plit_[pid];
	  pa.reset();
	  pa.poslb(firstpos);
	  pa.posub(pos);
	  pa.exact_start_bases(plit->exact_start_bases());
	  pa.exact_end_bases(plit->exact_end_bases());
	  // cerr << plit->pattern() << endl;
	  // checkpoint;
	  if (pa.align(cp,plit->pattern())) {
	    // checkpoint;
	    phs.push_back(pa.end(),make_pair(plit,pa.value()));
	  }
	  tinylist<pattern_hit_vector::iterator>::iterator adjit;
	  adjit = adj.begin();
	  while (adjit != adj.end()) {
	    (*adjit)->key() = 0;
	    ++adjit;
	  }
	}
      }
      ++it;
    }
    // Move everything with valid key to beginning...
    // int i,i1;
    it=l.begin();
    it1=l.begin();
    // i=0;
    // i1=0;
    long unsigned int newsize=0;
    while (it!=l.end()) {
      // cerr << i << " " << i1 << endl;
      // cerr << l.key(it) << " " << l.value(it)->id() << endl;
      if (it->key() != 0) {
	if (it != it1) {
	  it1->key() = it->key();
	  it1->value() = it->value();
	}
	++newsize;
	++it1;
	// i1++;
      }
      ++it;
      // i++;
    }
    // cerr << newsize << endl;
    l.resize(newsize);
    // checkpoint;
    // it=l.begin();
    // while (it!=l.end()) {
    // cerr << it->key() << " " << it->value()->id() << endl;
    // ++it;
    // }
    cp.pos(oldcharspos);
    report_progress(cp);
    if (phs.size() >= minka || 
	(more==false && phs.size() > 0)) return true;
  }
  return false;
}

void 
filter_bitvec::init(CharacterProducer & cp) {
  assert(pm_!=((void*)0));
  plit_.resize(num_patterns_+1);
  long unsigned int id;
  tinylist<pattern_list_element>::const_iterator it;
  for (it=patterns().begin();it!=patterns().end();++it) {
    id = pm_->add_pattern(it->pattern());
    plit_[id] = it;
  }
  pm_->init(cp);
}

void filter_bitvec::reset() {
  pm_->reset();
}



/*


      // checkpoint;
      std::list<pattern_hit*>::iterator it,it1,nextit;
      unsigned long int mindist=MAXINT;
      bool assigned_nextit=false;
      std::list<std::list<pattern_hit*>::iterator> adj;
      std::list<std::list<pattern_hit*>::iterator>::iterator it2;
      it = pas.begin();
      // cerr << "it: " << (*it)->pattern_id() << " " << (*it)->pos() << " " << ((pattern_hit_w_dist*)(*it))->editdist() << endl;
      long unsigned int pos;
      while (it!=pas.end()) {
	// checkpoint;
	it1=it;
	mindist=((pattern_hit_w_dist*)(*it1))->editdist();
	// cerr << "it1: " << (*it1)->pattern_id() << " " << (*it1)->pos() << " " << ((pattern_hit_w_dist*)(*it1))->editdist() << endl;
	// cerr << "mindist: " << mindist << endl;
	pos=(*it)->pos();
	// checkpoint;
	adj.clear();
	// checkpoint;
	assigned_nextit=false;
	while (it1 != pas.end() && 
	       (*it1)->pos() <= pos+2*k_+1) {
	  // checkpoint;
	  if ((*it)->pattern_id() == (*it1)->pattern_id()) {
	    if (mindist > ((pattern_hit_w_dist*)(*it1))->editdist()) {
	      // cerr << "mindist: " << mindist << endl;
	      mindist = ((pattern_hit_w_dist*)(*it1))->editdist();
	    }
	    adj.push_back(it1);
	  } else {
	    if (!assigned_nextit) {
	      nextit=it1;
	      assigned_nextit=true;
	    }
	  }
	  // checkpoint;
	  pos = (*it1)->pos();
	  it1++;
	  // if (it1 != pas.end()) {
	  // cerr << "it1: " << (*it1)->pattern_id() << " " << (*it1)->pos() << " " << ((pattern_hit_w_dist*)(*it1))->editdist() << endl;
	  // }
	}
	// checkpoint;
	// cerr << adj.size() << " " << mindist << endl;
	if (adj.size() > 1) {
	  // checkpoint;
	  it2  = adj.begin();
	  bool assigned_keptit=false;
	  std::list<pattern_hit*>::iterator keptit;
	  while (it2 != adj.end()) {
	    if (((pattern_hit_w_dist*)(**it2))->editdist() > mindist || 
		assigned_keptit) {
	      // cerr << "erase *it2: " << (**it2)->pattern_id() << " " << (**it2)->pos() << " " << ((pattern_hit_w_dist*)(**it2))->editdist() << endl;
	      delete **it2;
	      pas.erase(*it2);
	    } else {
	      // cerr << "keep  *it2: " << (**it2)->pattern_id() << " " << (**it2)->pos() << " " << ((pattern_hit_w_dist*)(**it2))->editdist() << endl;	      
	      if (!assigned_keptit) {
		keptit=*it2;
		assigned_keptit=true;
	      }
	    }
	    it2++;
	  }
	  if (assigned_nextit) {
	    it = nextit;
	  } else {
	    it = keptit;
	    it++;
	  }
	} else {
	  it++;
	}
	// if (it != pas.end()) {
	// cerr << "it: " << (*it)->pattern_id() << " " << (*it)->pos() << " " << ((pattern_hit_w_dist*)(*it))->editdist() << endl;
	// }
      }
      // checkpoint;






*/
