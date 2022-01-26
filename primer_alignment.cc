
#include "primer_alignment.h"
#include "util.h"
#include <list>

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

bool primer_alignment::global_align(char * const & text, unsigned int textlen,
				    std::string const & pattern, 
				    int dirn, 
				    unsigned int lmatch, 
				    unsigned int rmatch,
				    int & matchlen, int & value) {
  unsigned int textlenp1 = textlen+1;
  unsigned int patlen = pattern.length();
  unsigned int patlenp1 = patlen+1;
  
  if (!dp_ || matsize_ < textlenp1*patlenp1) {
    if (dp_) {
      delete [] dp_;
      delete [] best_;
    }
    // size should always be larger than (textlen+1)*(patlen+1)
    if ((maxpatlen_+1) < patlenp1) {
      maxpatlen_ = patlen;
    }
    matsize_ = (maxpatlen_+maxdist_+1)*(maxpatlen_+1);
    if (matsize_ < textlenp1*patlenp1) {
      matsize_ = textlenp1*patlenp1;
    }
    dp_ = new unsigned int[matsize_];
    best_ = new unsigned int[matsize_];
  }
#define index(i,j) ((i)*(textlenp1)+(j))
  // checkpoint;   
  // cerr << text << endl;
  // cerr << pattern << endl;
  // cerr << dirn << endl;
  // cerr << lmatch << endl;
  // cerr << rmatch << endl;
  // cerr << (int)eos_ << endl;
  int lbexact=0;
  int rbexact=patlenp1;
  int const_violation = 5*maxdist_+1;
  if (dirn < 0) {
    if (lmatch > 0) rbexact = patlenp1-lmatch;
    if (rmatch > 0) lbexact = rmatch;
  } else {
    if (lmatch > 0) lbexact = lmatch;
    if (rmatch > 0) rbexact = patlenp1-rmatch;
  }
  // cerr << dirn << " " << lbexact << " " << rbexact << endl;
  int p,t;
  long unsigned int ind,ind10,ind01,ind11;
  dp_[/*index(0,0)=*/0] = 0;
  // cerr << "dp_ " << 0 << " " << 0 << " " << dp_[index(0,0)] << endl;
  int lb,ub;
  ub = (indels_?(dna_mut_?1:maxdist_):0);
  if (patlen < ub) {
    ub = patlen;
  }
  for (p=1,ind=textlenp1,ind10=0;p<=ub;p++,ind10=ind,ind+=textlenp1) {
    // assert(index(p,0) == ind);
    // assert(index(p-1,0) == ind10);
    if (!indels_ || p < lbexact || p >= rbexact) {
      dp_[/*index(p,0)=*/ind] = const_violation; /* price away a possible
						   deletion of this pattern
						   character... */
      best_[/*index(p,0)=*/ind] = alignment_mask_constraint_violation;
    } else {
      if (!dna_mut_) {
	dp_[/*index(p,0)=*/ind] = dp_[/*index(p-1,0)=*/ind10]+1; /* cost for a deletion */
	best_[/*index(p,0)=*/ind] = alignment_mask_deletion;
      } else {
	dp_[/*index(p,0)=*/ind] = dp_[/*index(p-1,0)=*/ind10]+3; /* cost for a deletion */
	best_[/*index(p,0)=*/ind] = alignment_mask_deletion_3;
      }
    }
    // cerr << "dp_ " << p << " " << 0 << " " << dp_[index(p,0)] << endl;
  }
  ub = (indels_?(dna_mut_?1:maxdist_):0);
  if (textlen < ub) {
    ub = textlen;
  }
  char textch,patch;
  for (t=1,ind01=0;t<=ub;ind01=t++) {
    // assert(index(0,t) == t && index(0,t-1) == ind01);
    if (dirn>0) {
      textch = text[t-1];
    } else {
      textch = text[textlen-t];
    }
    if (!indels_ || 0 < lbexact || 0 >= rbexact || textch == eos_ || patch == eos_) {
      dp_[/*index(0,t)=*/t] = const_violation; /* price away a
						  possible
						  insertion at
						  this pattern
						  character... */
      best_[/*index(0,t)=*/t] = alignment_mask_constraint_violation;
    } else {
      if (!dna_mut_) {
	dp_[/*index(0,t)=*/t] = dp_[/*index(0,t-1)=*/ind01]+1; /* cost for a insertion */
	best_[/*index(0,t)=*/t] = alignment_mask_insertion;
      } else {
	dp_[/*index(0,t)=*/t] = dp_[/*index(0,t-1)=*/ind01]+3; /* cost for a insertion */
	best_[/*index(0,t)=*/t] = alignment_mask_insertion_3;	
      }
    }
    // cerr << "dp_ " << 0 << " " << t << " " << dp_[index(0,t)] << endl;
  }
  int v,v1;
  unsigned int ac,ac1;
  // checkpoint;
  for (p=1;p<=patlen;p++) {
    lb = p-(indels_?(dna_mut_?1:maxdist_):0);
    if (lb < 1) {
      lb = 1;
    }
    ub = p+(indels_?(dna_mut_?1:maxdist_):0);
    if (textlen < ub) {
      ub = textlen;
    }
    int bestscorerow=const_violation;
    // checkpoint;
    ind = index(p,lb);
    ind01 = ind-1;
    ind11 = ind01-textlenp1;
    ind10 = ind11+1;
    for (t=lb;t<=ub;t++) {
//       cerr << index(p,t) << " " << ind << endl;
//       cerr << index(p-1,t) << " " << ind10 << endl;
//       cerr << index(p-1,t-1) << " " << ind11 << endl;
//       cerr << index(p,t-1) << " " << ind01 << endl;
//       assert(index(p,t) == ind &&
// 	     index(p-1,t) == ind10 &&
// 	     index(p,t-1) == ind01 &&
// 	     index(p-1,t-1) == ind11);
      int ac=0;
      if (dirn>0) {
	textch = text[t-1];
	patch = pattern[p-1];
      } else {
	textch = text[textlen-t];
	patch = pattern[patlen-p];
      }
      if (textch == patch){
	v = dp_[/*index(p-1,t-1)=*/ind11]; /* characters match, get profit*/
	ac = alignment_mask_equal;
      } else if (wc_ && iupac_compatible(textch,patch) && 
		 (tn_ || textch != 'N')) {
	v = dp_[/*index(p-1,t-1)=*/ind11]; /* characters wildcard match, get profit*/
	ac = alignment_mask_wildcard_equal;
      } else if (textch == eos_ || patch == eos_ ||  p <= lbexact || p >=rbexact) {
	/* cannot tolerate a substitution error with eos char or if it
           is in the exact range... */
	v = const_violation;
	ac = alignment_mask_constraint_violation;
      } else {
	int mut=-1;
	if (!dna_mut_) {
	  v = dp_[/*index(p-1,t-1)=*/ind11]+1; /* cost for subst*/
	  ac = alignment_mask_substitution;
	} else if ((mut=aasubdist(textch,patch))>=0) {
	  v = dp_[/*index(p-1,t-1)=*/ind11]+mut; /* cost for subst*/
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
	  v = const_violation;
	  ac = alignment_mask_constraint_violation;
	}
      }
      // checkpoint;
      if (textch == eos_ || patch == eos_ || !indels_ || t <= lb ||
	  p < lbexact || p >= rbexact) {
	// checkpoint;
	v1 = const_violation; 
	                   /* cannot tolerate an insertion in the text
			      if that insertion is an eos or if the
			      pattern is in the exact range */
	ac1 = alignment_mask_constraint_violation;
      } else {
	// checkpoint;
	if (!dna_mut_) {
	  v1 = dp_[/*index(p,t-1)=*/ind01]+1; /* cost for an insertion */
	  ac1 = alignment_mask_insertion;
	} else {
	  v1 = dp_[/*index(p,t-1)=*/ind01]+3; /* cost for an insertion */
	  ac1 = alignment_mask_insertion_3;	  
	}
      }
      if (v1 < v) {
	// checkpoint;
	v = v1;
	ac = ac1;
      } else if (v1 == v) {
	ac |= ac1;
      }
      // checkpoint;
      if (!indels_ || t >= ub || p <= lbexact || p >= rbexact) {
	v1 = const_violation; 
	ac1 = alignment_mask_constraint_violation;
      } else {
	if (!dna_mut_) {
	  v1 = dp_[/*index(p-1,t)=*/ind10]+1; /* cost for a deletion. Note that we
						 permit a deletion in the text, even
						 if the text charcter is an eos... */
	  ac1 = alignment_mask_deletion;
	} else {
	  v1 = dp_[/*index(p-1,t)=*/ind10]+3; 
	  ac1 = alignment_mask_deletion_3;
	}
      }
      if (v1 < v) {
	v = v1;
	ac = ac1;
      } else if (v1 == v) {
	ac |= ac1;
      }
      dp_[/*index(p,t)=*/ind] = v;
      // cerr << "dp_ " << p << " " << t << " " << dp_[index(p,t)] << endl;
      best_[/*index(p,t)=*/ind] = ac;
      if (bestscorerow > v) {
	bestscorerow = v;
      }
      ind01 = ind++;
      ind11 = ind10++;
    }
    // checkpoint;
    if (bestscorerow > maxdist_) {
      // delete [] dp_;
      // delete [] best_;
      return false;
    } 
    // checkpoint;
  }
  // checkpoint;
  
  int val;
  int bestpos = patlen-(indels_?(dna_mut_?1:maxdist_):0);
  if (textlen < bestpos) {
    bestpos = textlen;
  }
  if (bestpos < 0) {
    bestpos = 0;
  }
  ind = index(patlen,bestpos);
  int k,bestval=dp_[ind]; 
  ub = patlen+(indels_?(dna_mut_?1:maxdist_):0);
  if (textlen < ub) {
    ub = textlen;
  }
  for (k=bestpos+1,ind++;k<=ub;k++,ind++) {
    // assert(index(patlen,k)==ind);
    if ((val=dp_[/*index(patlen,k)=*/ind]) < bestval ||
	((val=dp_[/*index(patlen,k)=*/ind]) <= bestval &&
	 (best_[/*index(patlen,k)=*/ind]&
	  (alignment_mask_equal|
	   alignment_mask_wildcard_equal|
	   alignment_mask_substitution|
	   alignment_mask_substitution_1|
	   alignment_mask_substitution_2|
	   alignment_mask_substitution_3)))) {
      bestval = val;
      bestpos = k;
    }
  }

  // checkpoint;

  p=patlen; t=bestpos;
  if (t < p - (indels_?(dna_mut_?1:maxdist_):0) || t > p + (indels_?(dna_mut_?1:maxdist_):0)) {
    // delete [] dp_;
    // delete [] best_;
    return false;
  } 
  if (yesno_) {
    // checkpoint;
    // cerr << bestval << endl;
    // cerr << bestpos << endl;
    // delete [] dp_;
    // delete [] best_;    
    matchlen = bestpos;
    value = bestval;
    return true;
  }
  // checkpoint;
  matching_text_.resize(bestpos);
  if (dirn>0) {
    for (int i=0;i<bestpos;i++) {
      matching_text_[i] = text[i];
    }
  } else {
    for (int i=0,j=textlen-bestpos;i<bestpos;i++,j++) {
      matching_text_[i] = text[j];
    }
  }
  std::list<alignment_code> alignment;
  std::list<alignment_code>::iterator alit;
  bool match,ins,del,wc,sub;
  alignment_code lastac=alignment_none;
  while (p!=0||t!=0) { 
    ac = best_[index(p,t)];
    match = (ac & (alignment_mask_equal |
		   alignment_mask_wildcard_equal |
		   alignment_mask_substitution |
		   alignment_mask_substitution_1 |
		   alignment_mask_substitution_2 |
		   alignment_mask_substitution_3)
	     );
    wc = (ac & alignment_mask_wildcard_equal);
    sub = (ac & (alignment_mask_substitution|
		 alignment_mask_substitution_1|
		 alignment_mask_substitution_2|
		 alignment_mask_substitution_3)
	   );
    ins = (ac & (alignment_mask_insertion|
		 alignment_mask_insertion_3)
	   );
    del = (ac & (alignment_mask_deletion|
		 alignment_mask_deletion_3)
	   );
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
      // Dump everything we know about the DP!
      // checkpoint;
      cerr << "Abort DP: Dumping everything!\n";
      cerr << "Text: " << text << endl;
      cerr << "Pattern: " << pattern << endl;
      cerr << "Dirn: " << dirn << endl;
      cerr << "LMatch: " << lmatch << endl;
      cerr << "RMatch: " << rmatch << endl;
      cerr << "DP Matrix: " << endl;
      for (p=0;p<=patlen;p++) {
	for (t=0;t<=textlen;t++) {
	  cerr << dp_[index(p,t)] << " ";
	}
	cerr << endl;
      }
      cerr << "Best Matrix: " << endl;
      for (p=0;p<=patlen;p++) {
	for (t=0;t<=textlen;t++) {
	  cerr << best_[index(p,t)] << " ";
	}
	cerr << endl;
      }
      cerr << "Bestval: " << bestval << endl;
      cerr << "Bestpos: " << bestpos << endl;
      cerr << "Matching Text: " << matching_text_ << endl;
      cerr << "Alignment chars: ";
      
      for (alit=alignment.begin();alit!=alignment.end();++alit) {
	cerr << *alit << " ";
      }
      cerr << endl;
      cerr << "Current p,t: " << p << "," << t << endl;
      assert(0);
    }
    if (dirn > 0) {
      alignment.push_front(lastac);
    } else if (dirn < 0) {
      alignment.push_back(lastac);
    }
    stats_[lastac]++;
  }
  // checkpoint;
  if (p != 0 || t!=0) {
    // Dump everything we know about the DP!
    cerr << "Abort DP: Dumping everything!\n";
    cerr << "Text: " << text << endl;
    cerr << "Pattern: " << pattern << endl;
    cerr << "Dirn: " << dirn << endl;
    cerr << "LMatch: " << lmatch << endl;
    cerr << "RMatch: " << rmatch << endl;
    cerr << "DP Matrix: " << endl;
    for (p=0;p<=patlen;p++) {
      for (t=0;t<=textlen;t++) {
	cerr << dp_[index(p,t)] << " ";
      }
      cerr << endl;
    }
    cerr << "Best Matrix: " << endl;
    for (p=0;p<=patlen;p++) {
      for (t=0;t<=textlen;t++) {
	cerr << best_[index(p,t)] << " ";
      }
      cerr << endl;
    }
    cerr << "Bestval: " << bestval << endl;
    cerr << "Bestpos: " << bestpos << endl;
    cerr << "Matching Text: " << matching_text_ << endl;
    cerr << "Alignment chars: ";
    for (alit=alignment.begin();alit!=alignment.end();++alit) {
      cerr << *alit << " ";
    }
    cerr << endl;
    cerr << "Current p,t: " << p << "," << t << endl;
    assert(0);
    // cerr << "Something is wrong with the dp..." << endl;
  }
  alignment_.resize(alignment.size());
  for (int i=0;i<alignment_.size();i++) {
    alignment_[i] = alignment.front();
    alignment.pop_front();
  }
  // delete [] dp_;
  // delete [] best_;
  matchlen = bestpos;
  value = bestval;
  return true;
}

bool primer_alignment_2match::align(CharacterProducer & cp, 
				    std::string const & pattern1, 
				    std::string const & pattern2) {
  bool retval = (((end2_ - pattern2.length()) >= end1_-1) &&
		 ((end2_ - pattern2.length()) <= end1_+1));
  if (yesno_) {
    end_ = end2_;
    end_defined_ = true;
    return retval;
  }
  // checkpoint;
  start_ = end1_ - pattern1.length();
  end_ = end2_;
  int txtlen = end_-start_;
  char *buffer = new char[txtlen+1];
  cp.pos(start_);
  for (int i=0;i<txtlen;i++) {
    buffer[i] = cp.getch();
  }
  buffer[txtlen] = '\0';
  // checkpoint;
  alignment_.resize(std::max(txtlen,(int)(pattern1.length() + pattern2.length())));
  // cerr << pattern1.length() << " " << pattern2.length() << " " << " " 
  // << txtlen << " " << alignment_.size() << endl;
  for (int i=0;i<pattern1.length();i++) {
    if (buffer[i] == pattern1[i]) {
      alignment_[i] = alignment_equal;
    } else if (wc_ && iupac_compatible(buffer[i],pattern1[i]) && 
	       (tn_ || buffer[i] != 'N')) {
      alignment_[i] = alignment_wildcard_equal;
    } else {
      // checkpoint;
      alignment_[i] = alignment_substitution;
    }
    // cerr << i << ": " << (int)alignment_[i] << endl;
  }
  for (int i=0;i<pattern2.length();i++) {
    if (buffer[txtlen-i-1] == pattern2[pattern2.length()-i-1]) {
      alignment_[alignment_.size()-i-1] = alignment_equal;
    } else if (wc_ && 
	       iupac_compatible(buffer[txtlen-i-1],pattern2[pattern2.length()-i-1]) && 
	       (tn_ || buffer[txtlen-i-1] != 'N')) {
      alignment_[alignment_.size()-i-1] = alignment_wildcard_equal;
    } else {
      alignment_[i] = alignment_substitution;
    }
    // cerr << alignment_.size()-i-1 << ": " << (int)alignment_[alignment_.size()-i-1] << endl;
  }
  if (pattern1.length() + pattern2.length() < txtlen) {
    if (buffer[pattern1.length()] != '\n' 
	&& lmatch_ < pattern1.length() 
	&& rmatch_ < pattern2.length()) {
      alignment_[pattern1.length()] = alignment_insertion;
    } else {
      alignment_[pattern1.length()] = alignment_constraint_violation;
    } 
  } else if (pattern1.length() + pattern2.length() > txtlen) {
    if (lmatch_ >= pattern1.length() 
	|| rmatch_ >= pattern2.length()) {
      alignment_[pattern1.length()] = alignment_constraint_violation;
    } else if (buffer[pattern1.length()-1] == pattern1[pattern1.length()-1]) {
      // checkpoint;
      // Position may match last char of first piece...
      // so put the deletion on the first char of the second piece...
      alignment_[pattern1.length()] = alignment_deletion;    
      // cerr << pattern1.length() << ": " << (int)alignment_[pattern1.length()] << endl;
    } else if (buffer[pattern1.length()-1] == pattern2[0]) {
      // checkpoint;
      // Position may match first char of second piece...
      // so put the deletion on the last char of the first piece...
      alignment_[pattern1.length()-1] = alignment_deletion;    
      // cerr << pattern1.length()-1 << ": " << (int)alignment_[pattern1.length()-1] << endl;
    } else if (wc_ && iupac_compatible(buffer[pattern1.length()-1],pattern1[pattern1.length()-1]) &&
	       (tn_ || buffer[pattern1.length()-1]!='N')) {
      // checkpoint;
      // Position may wildcard match last char of first piece...
      // so put the deletion on the first char of the second piece...
      alignment_[pattern1.length()] = alignment_deletion;    
      // cerr << pattern1.length() << ": " << (int)alignment_[pattern1.length()] << endl;
    } else if (wc_ && iupac_compatible(buffer[pattern1.length()-1],pattern2[0]) &&
	       (tn_ || buffer[pattern1.length()-1]!='N')) {
      // checkpoint;
      // Position may wildcard match first char of second piece...
      // so put the deletion on the last char of the first piece...
      alignment_[pattern1.length()-1] = alignment_deletion;    
      // cerr << pattern1.length()-1 << ": " << (int)alignment_[pattern1.length()-1] << endl;
    } else {
      // checkpoint;
      // No match, so put the deletion on the first char of the second piece...
      alignment_[pattern1.length()] = alignment_deletion;
      // cerr << pattern1.length() << ": " << (int)alignment_[pattern1.length()] << endl;
    }
    // checkpoint;
  } 
  for (int i=0;i<alignment_.size();i++) {
    stats_[alignment_[i]]++;
  }
  matching_text_ = buffer;
  delete [] buffer;
  alignment_done_ = true;
  return retval;
}

bool primer_alignment_lmatch::align(CharacterProducer & cp, 
				    std::string const & pattern1, 
				    std::string const & pattern2) {
  // checkpoint;
  // cerr << cp.pos() << " " << pattern1 << " " << pattern2 << endl;
  FILE_POSITION_TYPE textstart = end1_;
  int buflen=pattern2.length()+maxdist_;
  if (!buffer1_ || bufsize_ < buflen+1) {
    // checkpoint;
    if (buffer1_) {
      delete [] buffer0_;
      delete [] buffer1_;
    }
    // size should always be larger than buflen+1;
    bufsize_ = maxpatlen_+((indels_)?(dna_mut_?1:maxdist_):0);
    if (bufsize_ < buflen+1) {
      bufsize_ = buflen+1;
    }
    buffer1_ = new char[bufsize_];
    buffer0_ = new char[bufsize_];
    // checkpoint;
  }
  // checkpoint;
  char *buffer;
  if (textstart < bufstart_ || textstart + buflen >= bufend_) {
    // cerr << "LMatch reading: " << textstart << "-" << textstart+bufsize_ << endl;
    cp.pos(textstart);
    for (int i=0;i<bufsize_;i++) {
      buffer1_[i] = cp.getch();
    }
    bufstart_ = textstart;
    bufend_ = textstart+bufsize_;
    buffer = buffer1_;
  } else {
    buffer = buffer1_+(textstart-bufstart_);
  }
  // checkpoint;
  // cerr << buffer << endl;
  int matchlen;
  int value;
  bool retval = global_align(buffer,buflen,pattern2,1,lmatch_-pattern1.length(),rmatch_,matchlen,value);
  // checkpoint;
  if (yesno_) {
    // checkpoint;
    // delete buffer;
    end_ = end1_+matchlen;
    val_ = value;
    end_defined_ = true;
    return retval;
  }
  textstart = end1_-pattern1.length();
  cp.pos(textstart);
  int buflen0 = pattern1.length();
  for (int i=0;i<buflen0;i++) {
    buffer0_[i] = cp.getch();
  }
  buffer0_[buflen0] = '\0';
  matching_text_ = std::string(buffer0_)+matching_text_;
  start_ = end1_-pattern1.length();
  end_ = start_ + matching_text_.length();
  int oldsize=alignment_.size();
  alignment_.resize(oldsize+pattern1.length());
  for (int i=oldsize-1;i>=0;i--) {
    alignment_[i+pattern1.length()] = alignment_[i];
  }
  for (int i=0;i<pattern1.length();i++) {
    if (buffer0_[i] == pattern1[i]) {
      alignment_[i] = alignment_equal;
    } else if (wc_ && iupac_compatible(buffer0_[i],pattern1[i]) &&
	       (tn_ || buffer0_[i]!='N')) {
      alignment_[i] = alignment_wildcard_equal;
    } else {
      alignment_[i] = alignment_substitution;
    }
    stats_[alignment_[i]]++;
  }
  alignment_done_=true;
  val_ = value;
  // delete [] buffer;
  // delete [] buffer0;
  return retval;
}

bool primer_alignment_rmatch::align(CharacterProducer & cp, 
				    std::string const & pattern1, 
				    std::string const & pattern2) {
  // checkpoint;
  // cerr << end2_ << endl;
  // cerr << pattern1 << " " << pattern2 << endl;
  FILE_POSITION_TYPE textstart=0;
  int patlen = pattern1.length()+pattern2.length()+maxdist_;
  if (end2_>(FILE_POSITION_TYPE)patlen) {
    textstart = end2_-patlen;
  }
  int buflen = end2_-pattern2.length()-textstart;
  // cerr << textstart << " " << buflen << endl;
  if (!buffer1_ || bufsize_ < buflen+1) {
    // checkpoint;
    if (buffer1_) {
      delete [] buffer0_;
      delete [] buffer1_;
    }
    // size should always be larger than buflen+1;
    bufsize_ = maxpatlen_+((indels_)?(dna_mut_?1:maxdist_):0);
    if (bufsize_ < buflen+1) {
      bufsize_ = buflen+1;
    }
    buffer1_ = new char[bufsize_];
    buffer0_ = new char[bufsize_];
    // checkpoint;
  }
  // checkpoint;
  char *buffer;
  if (textstart < bufstart_ || textstart + buflen >= bufend_) {
    // cerr << "RMatch reading: " << textstart << "-" << textstart+bufsize_ << endl;
    cp.pos(textstart);
    for (int i=0;i<bufsize_;i++) {
      buffer1_[i] = cp.getch();
    }
    bufstart_ = textstart;
    bufend_ = textstart+bufsize_;
    buffer = buffer1_;
  } else {
    buffer = buffer1_+(textstart - bufstart_);
  }
  int matchlen;
  int value;
  bool retval = global_align(buffer,buflen,pattern1,-1,lmatch_,rmatch_-pattern2.length(),matchlen,value);
  // checkpoint;
  if (yesno_) {
    // checkpoint;
    // delete buffer;
    end_ = end2_;
    val_ = value;
    end_defined_ = true;
    return retval;
  }
  for (int i=0;i<pattern2.length();i++) {
    buffer0_[i] = cp.getch();
  }
  buffer0_[pattern2.length()] = '\0';
  matching_text_ = matching_text_+buffer0_;
  end_ = end2_;
  start_ = end_-matching_text_.length();
  int oldsize=alignment_.size();
  alignment_.resize(alignment_.size()+pattern2.length());
  for (int i=0;i<pattern2.length();i++) {
    if (buffer0_[i] == pattern2[i]) {
      alignment_[oldsize+i] = alignment_equal;
    } else if (wc_ && iupac_compatible(buffer0_[i],pattern2[i]) &&
	       (tn_ || buffer0_[i]!='N')) {
      alignment_[oldsize+i] = alignment_wildcard_equal;
    }
    stats_[alignment_[oldsize+i]] ++;
  }
  alignment_done_=true;
  // delete [] buffer;
  // delete [] buffer0;
  val_ = value;
  return retval;
}

//  void primer_alignment::write(ostream & os, FILE_POSITION_TYPE const seq_pos, 
//  			     std::string const & pattern1, 
//  			     std::string const & pattern2, 
//  			     long unsigned int id, bool revcomp) {
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
//    std::string kw = pattern1 + pattern2;
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
