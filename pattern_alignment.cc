
#include "pattern_alignment.h"
#include "util.h"

#include <list>
#include <iomanip>

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

istream & operator>>(istream & is, pattern_alignment & pa) {
  pa.read(is);
  return is;
}

ostream & operator<<(ostream & os, pattern_alignment const & pa) {
  pa.write(os);
  return os;
}

pattern_alignment::~pattern_alignment() {}

void pattern_alignment::write(ostream & os) const {
  os << pattern_ << " " << end_;
}


bool exact_alignment::align(CharacterProducer & /*cp*/, std::string const & pat) {
  if (yesno_) {
    return true;
  }
  matching_text_ = pat;
  start_ = end() - matching_text_.length();
  alignment_.resize(matching_text_.length());
  for (int i=0;i<matching_text_.length();i++) {
    alignment_[i] = alignment_equal;
  }
  stats_[alignment_equal] = matching_text_.length();
  value(0);
  alignment_done_ = true;
  return true;
}

bool exact_peptide_alignment::align(CharacterProducer & cp, std::string const & pat) {
  if (yesno_) {
    return true;
  }
  matching_text_ = pat;
  start_ = end() - matching_text_.length();
  alignment_.resize(matching_text_.length());
  for (int i=0;i<matching_text_.length();i++) {
    alignment_[i] = alignment_equal;
  }
  stats_[alignment_equal] = matching_text_.length();
  cp.pos(max(start_-context_,(FILE_POSITION_TYPE)0));
  unsigned int len=context_;
  if (start_ < context_) {
    len=start_;
  }
  lcontext_ = cp.getstr(len);
  matching_text_ = cp.getstr((unsigned)pat.length());
  // cp.pos(end());
  rcontext_ = cp.getstr((unsigned)context_);
  alignment_done_ = true;
  value(0);
  return true;
}

bool exact_wc_alignment::align(CharacterProducer & cp, std::string const & pat) {
  start_ = end() - pat.length();
  alignment_.resize(pat.length());
  char *buffer = new char[pat.length()+1];
  cp.pos(start_);
  for (int i=0;i<pat.length();i++) {
    buffer[i] = cp.getch();
    if (buffer[i]==pat[i]) {
      alignment_[i] = alignment_equal;
    } else if (iupac_compatible(buffer[i],pat[i]) && 
	       (textn_ || buffer[i] != 'N')) {
      alignment_[i] = alignment_wildcard_equal;
    } else {
      alignment_[i] = alignment_substitution;
    }
    stats_[alignment_[i]]++;
  }
  buffer[pat.length()] = '\0';
  matching_text_ = buffer;
  delete [] buffer;
  alignment_done_ = true;
  if (editdist()<=0) return true;
  return false;
}

bool mismatch_alignment::align(CharacterProducer & cp, std::string const & pat) {
  start_ = end() - pat.length();
  alignment_.resize(pat.length());
  char *buffer = new char[pat.length()+1];
  cp.pos(start_);
  for (int i=0;i<pat.length();i++) {
    buffer[i] = cp.getch();
    if (buffer[i]==pat[i]) {
      alignment_[i] = alignment_equal;
    } else {
      alignment_[i] = alignment_substitution;
    }
    stats_[alignment_[i]]++;
  }
  buffer[pat.length()] = '\0';
  matching_text_ = buffer;
  delete [] buffer;
  alignment_done_ = true;
  value(editdist());
  return true;
}

bool editdist_alignment::align(CharacterProducer & cp, std::string const & pat) {
  bool problem_case=false;
  // const bool problem_case(end() <= 217834338 && 217834338 <= end2_);
  if (problem_case) {
    checkpoint;
    cerr << "end():" << end() << endl;
    cerr << "end2_:" << end2_ << endl;
    cerr << "pat.length():" << pat.length() << endl;
    cerr << "k_:" << k_ << endl;
    cerr << "wc_:" << ((wc_)?"true":"false") << endl;
    cerr << "indels_:" << ((indels_)?"true":"false") << endl;
    cerr << "yesno_: " << ((yesno_)?"true":"false") << endl;
    cerr << "eos_: " << (int)eos_ << endl;
  } 
  int const_viol_penalty = 5*k_+1;
  FILE_POSITION_TYPE textstart=0;
  // checkpoint;
  // cerr << end() << endl;
  // cerr << pat.length() << endl;
  // cerr << k_ << endl;
  if (end() > pat.length() + k_) {
    textstart = end() - pat.length() - k_;
  }
  // if (problem_case) checkpoint;
  // if (problem_case) cerr << textstart << endl;
  // checkpoint;
  // cerr << end2_ << endl;
  // cerr << textstart << endl;
  // cerr << k_ << endl;
  // cerr << end() << endl;
  // cerr << pat.length() << endl;
  FILE_SIZE_TYPE buflen = end2_ - textstart;
  if (maxpatlen_ < pat.length()) {
    maxpatlen_ = pat.length();
  }
  // if (problem_case) cerr << buflen << endl;
  bool realloc=false;
  if (!buffer_ || (buflen+1)>bufsize_) {
    // if (problem_case) checkpoint;
    if (buffer_) {
      delete [] buffer_;
      if (dna_mut_ && trans_) {
        delete [] buffer1_;
      }
    }
    // size should always be larger than buflen+1;
    bufsize_ = 2*maxpatlen_+((indels_)?(dna_mut_?1:k_):0)+1;
    if (bufsize_<(buflen+1)) {
      bufsize_ = (buflen+1);
    }
    // checkpoint;
    // cerr << maxpatlen_ << " " << bufsize_ << endl;
    buffer_ = new char[bufsize_+1];
    if (dna_mut_ && trans_) {
      buffer1_ = new char[bufsize_+1];
    }
    realloc=true;
    // if (problem_case) checkpoint;
  }
  // if (problem_case) checkpoint;
  char *buffer;
  char *buffer1;
  if (textstart < bufstart_ || textstart + buflen >= bufend_ || realloc) {
    // if (problem_case) checkpoint;
    cp.pos(textstart);
    for (int i=0;i<bufsize_;i++) {
      buffer_[i] = cp.getch();
      if (dna_mut_ && trans_) {
        buffer1_[i] = cp.getcodonid();
      }
    }
    buffer_[bufsize_]='\0';
    // if (problem_case) cerr << bufsize_ << " " << buffer_ << endl;
    bufstart_ = textstart;
    bufend_ = textstart+bufsize_;
    buffer = buffer_;
    buffer1 = buffer1_;
  } else {
    // if (problem_case) checkpoint;
    buffer = buffer_+(textstart-bufstart_);
    buffer1 = buffer1_+(textstart-bufstart_);
  }
  if (problem_case) checkpoint;
  if (problem_case) cerr << pat.length() << " " << pat << endl;
  if (problem_case) cerr << buflen << " ";
  if (problem_case) {
    for (int i=0;i<buflen;i++) {
      cerr << buffer[i];
    }
  }
  if (problem_case) cerr << endl;

   // if (problem_case) checkpoint;
  if (!dp_ || 
      matsize_ < (buflen+1)*(pat.length()+1)) {
    if (dp_) {
      delete [] dp_;
    }
    if (best_) {
      delete [] best_;
    }
    // if (problem_case) checkpoint;
    matsize_ = (bufsize_+1)*(maxpatlen_+1);
    if (matsize_ < (buflen+1)*(pat.length()+1)) {
      matsize_ = (buflen+1)*(pat.length()+1);
    }
    // if (problem_case) cerr << matsize_ << endl;
    dp_ = new unsigned int[matsize_];
    // if (problem_case) checkpoint;
    best_ = new int[matsize_];
    // if (problem_case) checkpoint;
  }
  // if (problem_case) checkpoint;
  int lbexact=0;
  int rbexact=pat.length()+1;
  if (lconst_ > 0) rbexact = pat.length()+1-lconst_;
  if (rconst_ > 0) lbexact = rconst_;

  // if (problem_case) checkpoint;
  // if (problem_case) cerr << lbexact << " " << rbexact << endl;

  unsigned int *dp = dp_;
  int *best = best_;
#define index(i,j) ((i)*(buflen+1)+(j))
  int p,t;
  dp[index(0,0)] = 0; 
  best[index(0,0)] = alignment_mask_alignment_end; 
  // if (problem_case) cerr << "dp " << 0 << " " << 0 << " " << dp[index(0,0)] << endl;
  // if (problem_case) cerr << "bt " << 0 << " " << 0 << " " << best[index(0,0)] << endl;
  int lb,ub;
  ub = (indels_?(dna_mut_?1:k_):0);
  if (ub > pat.length()) {
    ub = pat.length();
  }
  // if (problem_case) checkpoint;
  // if (problem_case) cerr << "ub: " << ub << endl;
  for (p=1;p<=ub;p++) {
    if (!indels_ || p < lbexact || p >= rbexact || pat[pat.length()-p] == eos_) {
      dp[index(p,0)] = const_viol_penalty; 
      best[index(p,0)] = alignment_mask_constraint_violation;
    } else {
      if (!dna_mut_) {
	dp[index(p,0)] = dp[index(p-1,0)] + 1; /* cost for a deletion */
	best[index(p,0)] = alignment_mask_deletion;
      } else {
	dp[index(p,0)] = dp[index(p-1,0)] + 3; /* cost for a deletion */
	best[index(p,0)] = alignment_mask_deletion_3;
      }
    }
    // if (problem_case) cerr << "dp " << p << " " << 0 << " " << dp[index(p,0)] << endl;
    // if (problem_case) cerr << "bt " << p << " " << 0 << " " << best[index(p,0)] << endl;
  }
  ub = end2_-end()+(indels_?(dna_mut_?1:k_):0);
  if (buflen < ub) {
    ub = buflen;
  }
  // if (problem_case) checkpoint;
  // if (problem_case) cerr << "ub: " << ub << endl;
  // if (problem_case) cerr << "end2_ - end(): " << end2_ - end() << endl;
  for (t=1;t<=ub;t++) {
    if (t <= end2_ - end()) {
      dp[index(0,t)] = 0;
      best[index(0,t)] = alignment_mask_alignment_end;
    } else if (!indels_ || lbexact > 0) {
      dp[index(0,t)] = const_viol_penalty; 
      best[index(0,t)] = alignment_mask_constraint_violation;
    } else {
      if (!dna_mut_) {
	dp[index(0,t)] = dp[index(0,t-1)] + 1; /* cost for a insertion */
	best[index(0,t)] = alignment_mask_insertion;
      } else {
	dp[index(0,t)] = dp[index(0,t-1)] + 3; /* cost for a insertion */
	best[index(0,t)] = alignment_mask_insertion_3;	
      }
    }
    // if (problem_case) cerr << "dp " << 0 << " " << t << " " << index(0,t) << " " << dp[index(0,t)] << endl;
    // if (problem_case) cerr << "bt " << 0 << " " << t << " " << index(0,t) << " " << best[index(0,t)] << endl;
  }
  // if (problem_case) checkpoint;
  for (p=1;p<=pat.length();p++) {
    lb = p-(indels_?(dna_mut_?1:k_):0);
    if (lb < 1) {
      lb = 1;
    }
    ub = p+end2_-end()+(indels_?(dna_mut_?1:k_):0);
    if (buflen < ub) {
      ub = buflen;
    }
    // if (problem_case) checkpoint;
    // cerr << 0 << " " << 27 << " " << index(0,27) << " " << dp[index(0,27)] << endl;
    // cerr << 0 << " " << 27 << " " << index(0,27) << " " << best[index(0,27)] << endl;
    // if (problem_case) cerr << "lb: " << lb << endl;
    // if (problem_case) cerr << "ub: " << ub << endl;
    int bestscorerow=const_viol_penalty;
    for (t=lb;t<=ub;t++) {
      unsigned int v,v1;
      int ac=0;
      if (buffer[buflen-t] == pat[pat.length()-p]) {
	v = dp[index(p-1,t-1)]; /* characters match, no charge */
	ac = alignment_mask_equal;
      } else if (wc_ && iupac_compatible(pat[pat.length()-p],buffer[buflen-t]) && (buffer[buflen-t]!='N' || textn_)) {
	v = dp[index(p-1,t-1)]; /* characters match, no charge */
	ac = alignment_mask_wildcard_equal;	
      } else if (buffer[buflen-t] == eos_ || pat[pat.length()-p] == eos_ || 
		 p <= lbexact || p >=rbexact) {
	v = const_viol_penalty; 
	ac = alignment_mask_constraint_violation;	
      } else {
	int mut=-1;
	if (!dna_mut_) {
	  v = dp[index(p-1,t-1)]+1; /* cost for a substitution */
	  ac = alignment_mask_substitution;
	} else if (trans_ && ((mut=aacodonsubdist(buffer[buflen-t],buffer1[buflen-t],pat[pat.length()-p]))>=0)) {
	  v = dp[index(p-1,t-1)]+mut; /* cost for subst*/
	  switch(mut) {
	  case 1:
	    ac = alignment_mask_substitution_1;
	    break;
	  case 2:
	    ac = alignment_mask_substitution_2;
	    break;
	  case 3:
	    ac = alignment_mask_substitution_3;
	    break;
	  default:
	    checkpoint;
	    assert(0);
	  }
	} else if (!trans_ && ((mut=aasubdist(buffer[buflen-t],pat[pat.length()-p]))>=0)) {
	  v = dp[index(p-1,t-1)]+mut; /* cost for subst*/
	  switch(mut) {
	  case 1:
	    ac = alignment_mask_substitution_1;
	    break;
	  case 2:
	    ac = alignment_mask_substitution_2;
	    break;
	  case 3:
	    ac = alignment_mask_substitution_3;
	    break;
	  default:
	    checkpoint;
	    assert(0);
	  }
	} else {
	  v = const_viol_penalty;
	  ac = alignment_mask_constraint_violation;
	}
      }
      if (buffer[buflen-t] == eos_ || pat[pat.length()-p] == eos_ || !indels_ || t <= lb ||
	  p < lbexact || p >= rbexact) {
	v1 = const_viol_penalty; /* cost for an insertion */
	if (v1 < v) {
	  v = v1;
	  ac = alignment_mask_constraint_violation;
	}
      } else {
	if (!dna_mut_) {
	  v1 = dp[index(p,t-1)]+1; /* cost for an insertion */
	  if (v1 < v) { 
	    v = v1;
	    ac = alignment_mask_insertion;
	  } else if (v1 == v) { 
	    ac |= alignment_mask_insertion;
	  }
	} else {
	  v1 = dp[index(p,t-1)]+3; /* cost for an insertion */
	  if (v1 < v) { 
	    v = v1;
	    ac = alignment_mask_insertion_3;
	  } else if (v1 == v) { 
	    ac |= alignment_mask_insertion_3;
	  }
	}
      }
      if (!indels_ || pat[pat.length()-p] == eos_ || t >= ub || p <= lbexact || p >= rbexact) {
	v1 = const_viol_penalty; /* cost for a deletion */
	if (v1 < v) {
	  v = v1;
	  ac = alignment_mask_constraint_violation;
	}
      } else {
	if (!dna_mut_) {
	  v1 = dp[index(p-1,t)]+1; /* cost for a deletion */
	  if (v1 < v) {
	    v = v1;
	    ac = alignment_mask_deletion;
	  } else if (v1 == v) {
	    ac |= alignment_mask_deletion;
	  }
	} else {
	  v1 = dp[index(p-1,t)]+3; /* cost for a deletion */
	  if (v1 < v) {
	    v = v1;
	    ac = alignment_mask_deletion_3;
	  } else if (v1 == v) {
	    ac |= alignment_mask_deletion_3;
	  }
	}
      }
      dp[index(p,t)] = v;
      best[index(p,t)] = ac;
      // if (problem_case) cerr << "dp " << p << " " << t << " " << index(p,t) << " " << dp[index(p,t)] << endl;
      // if (problem_case) cerr << "bt " << p << " " << t << " " << index(p,t) << " " << best[index(p,t)] << endl;
      if (v < bestscorerow) {
	bestscorerow = v;
      }
    }
    if (bestscorerow > k_) {
      // if (problem_case) checkpoint;
      // if (problem_case) cerr << bestscorerow << " " << k_ << endl;
      if (!yesno_) {
	// checkpoint;
	alignment_.push_back(alignment_constraint_violation);
	stats_[alignment_constraint_violation]++;
      }
      alignment_done_ = true;
      // if (problem_case) checkpoint;
      return false;
    }
  }

  // cerr << 0 << " " << 27 << " " << index(0,27) << " " << dp[index(0,27)] << endl;
  // cerr << 0 << " " << 27 << " " << index(0,27) << " " << best[index(0,27)] << endl;
  // if (problem_case) checkpoint;

  int beststart=pat.length()-(indels_?(dna_mut_?1:k_):0);
  if (beststart>buflen) {
    beststart = buflen;
  }
  if (beststart < 0) {
    beststart = 0;
  }
//	cerr << "beststart: " << beststart << endl;
  int bestval=dp[index(pat.length(),beststart)]; 
  ub = pat.length()+end2_-end()+(indels_?(dna_mut_?1:k_):0);
  if (buflen < ub) {
    ub = buflen;
  }
  // checkpoint;
//	cerr << "ub: " << ub << endl;
//   if (problem_case) checkpoint;
//   if (problem_case) cerr << bestval << " " << beststart << endl;
//   if (problem_case) cerr << dp[index(pat.length(),beststart)] << " " << beststart << endl;
  for (int k=beststart+1;k<=ub;k++) {
    // if (problem_case) cerr << dp[index(pat.length(),k)] << " " << k << endl;
    if ((dp[index(pat.length(),k)]) < bestval ||
	((dp[index(pat.length(),k)]) <= bestval &&
	 (best[index(pat.length(),k)]&(alignment_mask_equal|
				       alignment_mask_wildcard_equal|
				       alignment_mask_substitution|
				       alignment_mask_substitution_1|
				       alignment_mask_substitution_2|
				       alignment_mask_substitution_3)
				       ))) {
      bestval = dp[index(pat.length(),k)];
      beststart = k;
    }
  }

  // if (problem_case) checkpoint;
  // if (problem_case) cerr << bestval << " " << beststart << endl;
 
  p=pat.length(); t=beststart;
  // if (problem_case) cerr << p << " " << t << endl;
  if (t < p - (indels_?(dna_mut_?1:k_):0) || t > p + (indels_?(dna_mut_?1:k_):0) + (end2_ - end())) {
    if (!yesno_) {
      alignment_.push_back(alignment_constraint_violation);
      stats_[alignment_constraint_violation]++;
    }
    // if (problem_case) checkpoint;
    alignment_done_ = true;
    return false;
  }

  // if (problem_case) checkpoint;

  std::list<alignment_code> alignment;
  std::list<alignment_code>::iterator alit;

  // if (problem_case) checkpoint;

  int ac;
  bool match,ins,del,wc,sub;
  // checkpoint;

  // cerr << 0 << " " << 27 << " " << index(0,27) << " " << dp[index(0,27)] << endl;
  // cerr << 0 << " " << 27 << " " << index(0,27) << " " << best[index(0,27)] << endl;

  assert(dp[index(0,0)] == 0); 
  assert(best[index(0,0)] == alignment_mask_alignment_end); 

  alignment_code lastac=alignment_none;
  // checkpoint;
  // cerr << p << " " << t << endl;
  // cerr << p << " " << t << " " << dp[index(p,t)] << endl;
  // cerr << p << " " << t << " " << best[index(p,t)] << endl;
  while (!(best[index(p,t)] & alignment_mask_alignment_end)) {
    // if (problem_case) cerr << p << " " << t << " " << index(p,t) << " " << dp[index(p,t)] << endl;
    // if (problem_case) cerr << p << " " << t << " " << index(p,t) << " " << best[index(p,t)] << endl;
    assert(p >= 0 && t >= 0);
    ac = best[index(p,t)];
    match = (ac & (alignment_mask_equal|
		   alignment_mask_wildcard_equal|
		   alignment_mask_substitution|
		   alignment_mask_substitution_1|
		   alignment_mask_substitution_2|
		   alignment_mask_substitution_3)
		   );
    wc = (ac & alignment_mask_wildcard_equal);
    sub = (ac & (alignment_mask_substitution|
		 alignment_mask_substitution_1|
		 alignment_mask_substitution_2|
		 alignment_mask_substitution_3)
		 );
    ins = (ac & (alignment_mask_insertion|
		 alignment_mask_insertion_3));
    del = (ac & (alignment_mask_deletion|
		 alignment_mask_deletion_3));
    // cerr << "eodp: " << ((ac & alignment_mask_alignment_end)?"true":"false") << " ";
    // cerr << "match: " << (match?"true":"false") << " ";
    // cerr << "wc: " << (wc?"true":"false") << " ";
    // cerr << "sub: " << (sub?"true":"false") << " ";
    // cerr << "ins: " << (ins?"true":"false") << " ";
    // cerr << "del: " << (del?"true":"false") << " ";
    // cerr << endl;
    if (match && !(((lastac == alignment_insertion ||
		     lastac == alignment_insertion_3) && ins) ||
		   ((lastac == alignment_deletion ||
		     lastac == alignment_deletion_3) && del) ||
		   (lastac == alignment_wildcard_equal && !wc && (ins || del)))) {
      p--; t--;
      if ((ac & alignment_mask_equal) && 
	  !((lastac == alignment_wildcard_equal && wc) || 
	    (lastac == alignment_substitution && sub))) {
	lastac = alignment_equal;
      } else if (wc) {
	lastac = alignment_wildcard_equal;
      } else if (sub) {
	if (ac & alignment_mask_substitution) {
	  lastac = alignment_substitution;
	} else if (ac & alignment_mask_substitution_1) {
	  lastac = alignment_substitution_1;
	} else if (ac & alignment_mask_substitution_2) {
	  lastac = alignment_substitution_2;
	} else if (ac & alignment_mask_substitution_3) {
	  lastac = alignment_substitution_3;
	}
      }
    } else if (del) {
      p--;
      if (ac & alignment_mask_deletion) {
	lastac = alignment_deletion;
      } else if (ac & alignment_mask_deletion_3) {
	lastac = alignment_deletion_3;
      }
    } else if (ins) {
      t--;
      if (ac & alignment_mask_insertion) {
	lastac = alignment_insertion;
      } else if (ac & alignment_mask_insertion_3) {
	lastac = alignment_insertion_3;
      }
    } else if (ac & alignment_mask_constraint_violation) {
      p = 0; t = 0;
      lastac = alignment_constraint_violation;
    } else {
      assert(0);
    }
    if (!yesno_) {
      stats_[lastac]++;
      alignment.push_back(lastac);
    }
  }
  if (!yesno_) {
    alignment_.resize(alignment.size());
    for (int i=0;i<alignment_.size();i++) {
      alignment_[i] = alignment.front();
      alignment.pop_front();
    }
  }
  int endt = t;
  int endp = p;
  
  // if (problem_case) cerr << buflen << " " << beststart << " " << endt << endl;

  matching_text_ = "";
  for (int i=0;i<(buflen-(buflen-beststart)-endt);i++) {
    matching_text_ += *(buffer+(buflen-beststart)+i);
  }
  // if (problem_case) cerr << matching_text_ << endl;
  start_ = end2_ - beststart;
  end(start_ + matching_text_.length());
  value(bestval);
  // if (problem_case) cerr << start_ << endl;
  // if (problem_case) cerr << end() << endl;

  // if (problem_case)   checkpoint;

  if (/*problem_case || */ endp < 0 || endt < 0) {
    cerr << "Something is wrong with the dp..." << endl;
    // Dump everything we know about the DP!
    checkpoint;
    cerr << "Abort DP: Dumping everything!\n";
    cerr << "Text: ";
    for (t=0;t<buflen;t++) {
      cerr << buffer[t];
    } 
    cerr << endl;
    cerr << "Pattern: " << pat << endl;
    cerr << "DP Matrix: " << endl;
    for (t=0;t<=buflen;t++) {
      if (t>0) {
	cerr << setw(3) << (char)buffer[buflen-t] << " ";
      } else {
	cerr << "         ";
      }
    }
    cerr << endl;
    // cerr << end2_ << " " << end() << endl;
    for (p=0;p<=pat.length();p++) {
      lb = p-(indels_?(dna_mut_?1:k_):0);
      if (lb < 0) {
	lb = 0;
      }
      ub = p+end2_-end()+(indels_?(dna_mut_?1:k_):0);
      if (buflen < ub) {
	ub = buflen;
      }
      // cerr << "p = " << p << ", t lb = " << lb << " ub = " << ub << endl;
      if (p > 0) {
	cerr << pat[pat.length()-p] << ": ";
      } else {
	cerr << " : ";
      }
      for (t=0;t<lb;t++) {
	cerr << "    ";
      }
      for (t=lb;t<=ub;t++) {
	cerr << setw(3) << dp[index(p,t)] << " ";
      }
      cerr << endl;
    }
    cerr << "Best Matrix: " << endl;
    for (t=0;t<=buflen;t++) {
      if (t>0) {
	cerr << setw(3) << (char)buffer[buflen-t] << " ";
      } else {
	cerr << "         ";
      }
    }
    cerr << endl;
    for (p=0;p<=pat.length();p++) {
      lb = p-(indels_?(dna_mut_?1:k_):0);
      if (lb < 0) {
	lb = 0;
      }
      ub = p+end2_-end()+(indels_?(dna_mut_?1:k_):0);
      if (buflen < ub) {
	ub = buflen;
      }
      // cerr << "p = " << p << ", t lb = " << lb << " ub = " << ub << endl;
      if (p > 0) {
	cerr << pat[pat.length()-p] << ": ";
      } else {
	cerr << " : ";
      }
      for (t=0;t<lb;t++) {
	cerr << "    ";
      }
      for (t=lb;t<=ub;t++) {
	cerr << setw(3) << best[index(p,t)] << " ";      }
      cerr << endl;
    }
    cerr << "Bestval: " << bestval << endl;
    cerr << "Bestpos: " << beststart << endl;
    cerr << "Matching Text: " << matching_text_ << endl;
    cerr << "Alignment chars: ";
    
    for (int i=0;i<alignment_.size();i++) {
      cerr << alignment_[i] << " ";
    }
    cerr << endl;
    cerr << "Current p,t: " << endp << "," << endt << endl;
    // abort();
  }
  alignment_done_ = true;
  return (bestval<=k_);
}

//  void alignment_write(ostream & os, 
//  			      FILE_POSITION_TYPE const seq_pos,
//  			      std::string const & pattern, 
//  			      long unsigned int id, bool revcomp) {
//    assert(alignment_done_);
//    std::string const & mt = matching_text();
//    int p=0;
//    os << " ";
//    for (int i=0; i<size(); i++) {
//      if ((*this)[i] != alignment_deletion) {
//        os << mt[p];
//        p++;
//      } else {
//        os << "-";
//      }
//    }
//    os /* << " " << start() << " " << end() */
//       << " " << seq_pos - length() + 1 << " " << seq_pos
//       << " " << editdist();
//    os << endl;
//    os << " ";
//    for (int i=0; i<size(); i++) {
//      if ((*this)[i] == alignment_equal) {
//        os << "|";
//        p++;
//      } else if ((*this)[i] == alignment_substitution) {
//        os << "*";
//      } else {
//        os << " ";
//      }
//    }
//    os << endl;
//    std::string kw = pattern;
//    p = 0;
//    os << " ";
//    for (int i=0; i<size(); i++) {
//      if ((*this)[i] != alignment_insertion) {
//        os << kw[p];
//        p++;
//      } else {
//        os << "-";
//      }
//    }
//    os << " " << id;
//    if (revcomp) {
//      os << " REVCOMP";
//    }
//    os << endl;
//  }
