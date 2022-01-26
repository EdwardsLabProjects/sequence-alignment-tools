
#include "select.h"

#include "keyword_tree.h"
#include "keyword_tree.t"
#include "shift_and.h"
#include "shift_and_inexact.h"
#include "exact_bases.h"
#include "exact_halves.h"
#include "filter_bitvec.h"
#include "hash_table.h"
#include "rlst.h"
#ifndef NOPRIMEGEN
#include "rand_hash_table.h"
#include "gs_hash_table.h"
#endif
#include "math.h"

PatternMatch * pick_pattern_index(CharacterProducer * const & ff,
				  int pmselect,
				  int nmismatch,
				  std::vector<std::pair<int,int> > *exact_const,
				  std::vector<int> *patlen,
				  int seedlen,
				  bool wildcard,
				  bool textn,
				  bool indels,
				  bool dna_mut,
				  char eos,
				  bool verbose) {
  long int min_exact_const=0;
  long int cumbooldiff=0;
  long int cumdiff=0;
  long int cum_exact_const=0;
  long int min_inexact_bases=MAXINT;
  double avexact=0;
  double avbool=0;
  double avdiff=0;
  if (exact_const) {
    min_exact_const=MAXINT;
    for (int i=1;i<(*exact_const).size();i++) {
      if ((*exact_const)[i].first >= (*exact_const)[i].second) {
	if (min_exact_const > (*exact_const)[i].first) {
	  min_exact_const = (*exact_const)[i].first;
	}
	cum_exact_const += (*exact_const)[i].first;
	cumdiff += ((*exact_const)[i].first - ((*patlen)[i]/2));
	cumbooldiff += ((((*exact_const)[i].first - ((*patlen)[i]/2))>=0)?1:0);
        if (min_inexact_bases > -((*exact_const)[i].first - (*patlen)[i])) {
          min_inexact_bases = -((*exact_const)[i].first - (*patlen)[i]);
        }
      } else /* (*exact_const)[i].first < (*exact_const)[i].second */ {
	if (min_exact_const > (*exact_const)[i].second) {
	  min_exact_const = (*exact_const)[i].second;
	}
	cum_exact_const += (*exact_const)[i].second;
	cumdiff += ((*exact_const)[i].second - ((*patlen)[i]/2));
	cumbooldiff += ((((*exact_const)[i].second - ((*patlen)[i]/2))>=0)?1:0);
        if (min_inexact_bases > -((*exact_const)[i].second - (*patlen)[i])) {
          min_inexact_bases = -((*exact_const)[i].second - (*patlen)[i]);
        }
      }
    }
    avexact=((double)cum_exact_const)/((*exact_const).size()-1);
    avbool=((double)cumbooldiff)/((*exact_const).size()-1);
    avdiff=((double)cumdiff)/((*exact_const).size()-1);
  } 

  long int min_length=0;
  long int cum_len=0;
  double avlen=0;
  if (patlen) {
    min_length=MAXINT;
    for (int i=1;i<(*patlen).size();i++) {
      if (min_length > (*patlen)[i]) {
	min_length = (*patlen)[i];
      }
      cum_len += (*patlen)[i];
    }
    avlen=((double)cum_len)/((*patlen).size()-1);
  }

  if (min_inexact_bases > min_length) {
    min_inexact_bases = min_length;
  }

  if (exact_const && nmismatch >= min_inexact_bases && nmismatch > 0) {
    timestampli("Fatal error: Number of edits >= Minimum number of inexact bases: ",min_inexact_bases);
    exit(1);
  }

  /*
  if (dna_mut && ff->size() < 20) {
    timestamp("Fatal error: DNA Mutations for non amino-acid alphabet (size < 20)");
    exit(1);
  }
  */

#ifndef NOPRIMEGEN
  gapped_seed_set::scheme gss = gapped_seed_set::no_scheme;
#endif
  PatternMatch * pm=0;
  if (pmselect == 0) {
    if (wildcard) {
      pmselect = 4;
    } else if (ff->size() < 255) {
      if (ff->nch('A') == 0 && 
	  ff->nch('C') == 1 &&
	  ff->nch('G') == 2 &&
	  ff->nch('T') == 3) {
	pmselect = 2;
      } else {
	pmselect = 3;
      }
    } else {
      pmselect = 3;
    }
    // pmselect is now the offset of the best exact pattern lookup.
    // Now determine if we can use any of our tricks for inexact lookup.
    if (nmismatch>0) {
      if (nmismatch == 1 && 
		 ((min_length >= 12 && ff->size() < 10) || 
		  (min_length >= 8 && ff->size() >=10)) &&
		 (cumbooldiff <= 0 || cumdiff <= 0)) {
	// Use exact_halves...
	pmselect = 11 + pmselect - 1;
#ifndef NOPRIMEGEN
      } else if ((gss = gapped_seed_set::select(min_length,nmismatch,indels)) != gapped_seed_set::no_scheme) {
	pmselect = 15;
#endif
      } else if (min_exact_const >= 6) {
	// Use exact_bases...
	pmselect = 7 + pmselect - 1;
      } else if (seedlen > 0) {
	// use hash table
	pmselect = 6;
      } else {
	// Use inexact bitvector...
	pmselect = 5;
      }
    }
  }

  if (dna_mut && pmselect == 5) {
    timestamp("Fatal error: DNA Mutations incompatible with inexact bitvector iimplementation");
    exit(1);
  }

  if (verbose && patlen) {
    timestampli("Primer stats: min length: ",min_length);
    timestampd("              average len: ",floor(avlen*10+.5)/10);
    if (nmismatch>0) {
      timestampli("              min exact bases: ",min_exact_const);
      timestampd("              average exact: ",floor(avexact*10+.5)/10);
      timestampd("              average (exact - len/2): ",floor(avdiff*10+.5)/10);
      timestampli("              count (exact >= len/2): ",cumbooldiff);
      timestampli("              seed length: ",(long)seedlen);
    } 
    timestampli("              number of primers: ",((*patlen).size()-1));
  }

  if (verbose) {
    if (indels) {
      timestampi("Options summary: string edits: ",nmismatch);
    } else {
      timestampi("Options summary: mismatches: ",nmismatch);      
    }
    if (dna_mut) {
      timestamp("                 DNA mutation scoring");
    } 
    if (wildcard) {
      if (textn) {
	timestamp("                 wildcard, w/ text N");
      } else {
	timestamp("                 wildcard, no text N");
      }
    } else {
      timestamp("                 no wildcard");      
    }
  }
    
  switch (pmselect) {
  case 1:
    if (ff->has_filename() && exist(suftree::filename(ff->filename()))) {
      pm = new suftree(eos);
      if (verbose) timestamp("Using suffix tree...");
    } else {
      pm = new keyword_tree<ktnode_list>;
      if (verbose) timestamp("Using keyword tree with list nodes...");
    }
    break;
  case 2:
    if (ff->has_filename() && exist(suftree::filename(ff->filename()))) {
      pm = new suftree(eos);
      if (verbose) timestamp("Using suffix tree...");
    } else {
      pm = new keyword_tree<ktnode_dna_list>;
      if (verbose) timestamp("Using keyword tree with nodes optimized for DNA...");
    }
    break;
  case 3:
    if (ff->has_filename() && exist(suftree::filename(ff->filename()))) {
      pm = new suftree(eos);
      if (verbose) timestamp("Using suffix tree...");
    } else {
      pm = new keyword_tree<ktnode_jtable>;
      if (verbose) timestamp("Using keyword tree with jump table nodes...");
    }
    break;
  case 4:
    pm = new shift_and(wildcard,textn);
    if (verbose) timestamp("Using bitvector...");
    break;
  case 5:
    pm = new filter_bitvec(nmismatch,eos,wildcard,textn,indels,dna_mut);
    if (verbose) timestamp("Using inexact bitvector...");
    break;
  case 6:
    if ( (log((double)ff->size())/log((double)2))*seedlen <= 25) {
      pm = new hash_table(seedlen,nmismatch,eos,wildcard,textn,indels,dna_mut);
      if (verbose) timestamp("Using exact seed with hash table...");
    } 
#ifndef NOPRIMEGEN
    else { 
      pm = new rand_hash_table(2,seedlen,nmismatch,eos,wildcard,textn,indels,dna_mut);
      if (verbose) timestamp("Using (large) exact seed with randomized hash table...");
    }
#endif
    break;
  case 7:
    pm = new exact_bases(new keyword_tree<ktnode_list>,nmismatch,eos,wildcard,textn,indels,dna_mut);
    if (verbose) timestamp("Using keyword tree with list nodes for exact portion...");
    break;
  case 8:
    pm = new exact_bases(new keyword_tree<ktnode_dna_list>,nmismatch,eos,wildcard,textn,indels,dna_mut);
    if (verbose) timestamp("Using keyword tree with nodes optimized for DNA for exact portion...");
    break;
  case 9:
    pm = new exact_bases(new keyword_tree<ktnode_jtable>,nmismatch,eos,wildcard,textn,indels,dna_mut);
    if (verbose) timestamp("Using keyword tree with jump table nodes for exact portion...");
    break;
  case 10:
    pm = new exact_bases(new shift_and(wildcard,textn),nmismatch,eos,wildcard,textn,indels,dna_mut);
    if (verbose) timestamp("Using bitvector for exact portion...");
    break;
  case 11:
    pm = new exact_halves(new keyword_tree<ktnode_list>,nmismatch,eos,wildcard,textn,indels,dna_mut);
    if (verbose) timestamp("Using keyword tree with list nodes for exact halves...");
    break;
  case 12:
    pm = new exact_halves(new keyword_tree<ktnode_dna_list>,nmismatch,eos,wildcard,textn,indels,dna_mut);
    if (verbose) timestamp("Using keyword tree with nodes optimized for DNA for exact halves...");
    break;
  case 13:
    pm = new exact_halves(new keyword_tree<ktnode_jtable>,nmismatch,eos,wildcard,textn,indels,dna_mut);
    if (verbose) timestamp("Using keyword tree with jump table nodes for exact halves...");
    break;
  case 14:
    pm = new exact_halves(new shift_and(wildcard,textn),nmismatch,eos,wildcard,textn,indels,dna_mut);
    if (verbose) timestamp("Using bitvector for exact halves...");
    break;
#ifndef NOPRIMEGEN
  case 15:
    if (gss != gapped_seed_set::no_scheme) {
      pm = new gs_hash_table(gss,2,nmismatch,eos,wildcard,textn,indels,dna_mut);
      if (verbose) timestamps("Using gapped seed set, scheme ",
			      ((gs_hash_table*)pm)->scheme().str().c_str());
    }  
    break;
#endif
  default:
    timestamp("Bad pattern index selected...");
    assert(0);
    break;
  }
  assert(pm != NULL);
  return pm;
}
