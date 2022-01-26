#include <stdio.h>
#include <iostream>
using namespace std;

#include "rl_time.h"
#include "rl_charmap.h"
static charmap CharMap;

#include "rl_index.h"
#include "xspacefsm.h"

#include "types.h"

static unsigned VERBOSE=0;
static char *oname = NULL;
static char *iname = NULL;
static char *fname = NULL;
static char *rname = NULL;
static char *whoami;
static unsigned MB=4000;
static bool planonly = false;
static unsigned mersize = 0;
static bool amino_acid_map = false;
static bool uppercase_map = false;

void usage(FILE *f) {
  fprintf(f,
	  "%s: [-v] [-h] [-p] [-A|-U] [-M MB limit] -k mersize -o output -i index_file -f fwd_file [ -r revc_file ]\n",
	  whoami);
}

static int get_args(int argc, const char **argv) {
  int errflg = 0;
  whoami = strdup(*argv);
  argv++;
#define OPT_ARG ( (argv[0][1])?(++(*argv)):(*(++argv)) )
  while(!errflg && *argv) {
    if (**argv == '-') {
      (*argv)++;
      while (!errflg && **argv) {
        switch (**argv) {
        case 'v':
          VERBOSE++;
          goto loopin;
        case 'p':
          planonly = true;
          goto loopin;
        case 'A':
          amino_acid_map = true;
          goto loopin;
        case 'U':
          uppercase_map = true;
          goto loopin;
        case 'M':
          MB = atoi(OPT_ARG);
          goto loopout;
        case 'o':
          oname = strdup(OPT_ARG);
          goto loopout;
        case 'i':
          iname = strdup(OPT_ARG);
          goto loopout;
        case 'f':
          fname = strdup(OPT_ARG);
          goto loopout;
        case 'r':
          rname = strdup(OPT_ARG);
          goto loopout;
        case 'k':
          mersize = atoi(OPT_ARG);
          goto loopout;
        case 'h':
          usage(stdout);
          exit(0);
          goto loopin;
        default:
          cerr << whoami << ": unknown flag '-" << *argv <<"'"<< endl;
          errflg++;
        }
      loopin:
        (*argv)++;
      }
    loopout:
      argv++;
    }
    else {
      errflg++;
    }
  }
  if (!oname) ++errflg;
  if (!iname) ++errflg;
  // if (!rname) ++errflg;
  if (!fname) ++errflg;
  if (!mersize) ++errflg;
  return errflg;
}

// the tree can only go up to 2^30 bases of sequence
// each node of the tree and fsm cost 20+1+2 bytes
static const unsigned long MaxTreeMB = ((1*23)<<10);
inline unsigned long Bytes(const min_index_elt &e) {
  return 1+((e.stop - e.start+1)*23);
}

int form_Plan(list<min_index_list> &P,const min_index_list &L) {
  unsigned long part_B;
  min_index_list part;

  // cerr << "Input MB: " << MB << endl;
  // cerr << "MaxTreeMB: " << MaxTreeMB << endl;

  unsigned long MaxMB = MB;
  if (MaxMB > MaxTreeMB) MaxMB = MaxTreeMB;

  P.clear();

  min_index_list::const_iterator i;
  for( part.clear(), part_B = 0, i = L.begin();
       i != L.end(); ) {
    if ( Bytes(*i) + part_B < (MaxMB<<20)) {
      part.push_back(*i);
      part_B += Bytes(*i);
      i++;
    }
    else {
      if (part.empty()) return 1;
      P.push_back(part);
      part.clear();
      part_B = 0;
    }
  }
  if (!part.empty()) {
    P.push_back(part);
  }
  return 0;
}

int main(int argc,const char **argv) {
  char MAP1[256],MAP2[256];
  TIC

  if (get_args(argc,argv)) {
    usage(stderr);
    exit(1);
  }
  FILE *stdverb = stdout;

  min_index_list li;

  FILE *o;
  if (!(o = fopen(oname,"w"))) {
    fprintf(stderr,"%s:main: unable to open %s\n",whoami,oname);
    exit(101);
  }
  FILE *fwd,*rev=0,*idx;
  if (!(idx = fopen(iname,"r"))) {
    fprintf(stderr,"%s:main: unable to open %s\n",whoami,iname);
    exit(102);
  }
  if (!(fwd = fopen(fname,"r"))) {
    fprintf(stderr,"%s:main: unable to open %s\n",whoami,fname);
    exit(103);
  }
  if (rname && !(rev = fopen(rname,"r"))) {
    fprintf(stderr,"%s:main: unable to open %s\n",whoami,rname);
    exit(104);
  }

  if (li.iload(idx)) {
    fprintf(stderr,"%s:main: trouble loading index from %s\n",whoami,iname);
    exit(105);
  }
  if (VERBOSE) {
    TOC(stdverb) fprintf(stdverb,": loaded the index\n");
  }

  list<min_index_list> Plan;
  if (form_Plan(Plan,li)) {
    fprintf(stderr,"%s:main: unable to form a decent plan\n",whoami);
    exit(106);
  }

  if (VERBOSE>0) {
    TOC(stdverb) fprintf(stdverb,": The plan has %lu blocks\n",Plan.size());
    if(VERBOSE>1) {
      for(list<min_index_list>::const_iterator i = Plan.begin();
	  i != Plan.end(); ++i) {
	long unsigned total = 0;
	for(min_index_list::const_iterator j = i->begin();
	    j != i->end(); ++j) {
	  if (VERBOSE>2) fprintf(stdverb,"%llu ",j->stop - j->start-1);
	  total += j->stop - j->start;
	}
	if (VERBOSE>2)fprintf(stdverb,"\n",total);
	fprintf(stdverb,"start: %12llu stop: %12llu total: %12llu\n",
		i->begin()->start,i->rbegin()->stop,total);
      }
    }
    fflush(stdverb);
  }

  if (planonly) exit(0);

  if (amino_acid_map) {
    for(int i=0; i < 256; ++i) {
      MAP1[i] = MAP2[i] = CharMap.aminoacid()[i];
    }
    // turn X's into mismatches
    MAP1[(unsigned char)CharMap.term1()] = CharMap.term1();
    MAP1[(unsigned char)CharMap.term2()] = CharMap.term1();
    for(int i=0; i < 256; ++i) if (MAP1[i] == 'X') MAP1[i] = CharMap.term1();
    
    MAP2[(unsigned char)CharMap.term1()] = CharMap.term2();
    MAP2[(unsigned char)CharMap.term2()] = CharMap.term2();
    for(int i=0; i < 256; ++i) if (MAP2[i] == 'X') MAP2[i] = CharMap.term2();
  } else if (uppercase_map) {
    for(int i=0; i < 256; ++i) {
      MAP1[i] = MAP2[i] = CharMap.uppercase()[i];
    }
    // turn non-letter's into mismatches
    MAP1[(unsigned char)CharMap.term1()] = CharMap.term1();
    MAP1[(unsigned char)CharMap.term2()] = CharMap.term1();
    for(int i=0; i < 256; ++i) if (MAP1[i] == CharMap.term3()) MAP1[i] = CharMap.term1();
    
    MAP2[(unsigned char)CharMap.term1()] = CharMap.term2();
    MAP2[(unsigned char)CharMap.term2()] = CharMap.term2();
    for(int i=0; i < 256; ++i) if (MAP2[i] == CharMap.term3()) MAP2[i] = CharMap.term2();
  } else /* DNA */ {
    for(int i=0; i < 256; ++i) {
      MAP1[i] = MAP2[i] = CharMap.canonical()[i];
    }
    MAP1[(unsigned char)CharMap.term1()] = CharMap.term1();
    MAP1[(unsigned char)CharMap.term2()] = CharMap.term1();
    for(int i=0; i < 256; ++i) if (MAP1[i] == 'N') MAP1[i] = CharMap.term1();
    
    MAP2[(unsigned char)CharMap.term1()] = CharMap.term2();
    MAP2[(unsigned char)CharMap.term2()] = CharMap.term2();
    for(int i=0; i < 256; ++i) if (MAP2[i] == 'N') MAP2[i] = CharMap.term2();
  }

  long unsigned Slen = 1;
  unsigned Snum = 0;
  unsigned maxlen = 0;
  for(min_index_list::const_iterator i = li.begin() ;
      i != li.end(); ++i) {
    Slen += i->stop - i->start;
    Snum ++;
    if (maxlen < i->stop - i->start) maxlen = i->stop - i->start;
  }
  

  char *S = new char[maxlen+1];
  FILE_POSITION_TYPE pos=1;
  for(min_index_list::const_iterator i = li.begin() ; i != li.end(); ++i) {
      min_sequence seq(*i);
      seq.chars = S;
      seq.sload(fwd);

      CharMap.map(MAP1,S,S+maxlen);

      if (pos == 1) {
	// fake node to give the term character
	fprintf(o," %llu.%c %llu.%c\n",0,S[0],Slen-1,S[0]);
	// fake nodes to mark sequence boundaries
	fprintf(o," %llu.%c %llu.%c\n",1,S[1],Slen,S[1]);
	fprintf(o," %llu.%c\n",1+mersize-1,S[1+mersize-1]);
	fprintf(o," %llu.%c\n",1+mersize,S[1+mersize]);
      }
      else {
	fprintf(o," %llu.%c\n",pos-1,S[0]);
	fprintf(o," %llu.%c\n",pos,S[1]);
	fprintf(o," %llu.%c\n",pos+mersize-1,S[1+mersize-1]);
	fprintf(o," %llu.%c\n",pos+mersize,S[1+mersize]);
      }
      
      pos += i->stop - i->start;
  }
  fflush(o);
  delete [] S;

  int parti = 0;
  for(list<min_index_list>::const_iterator i = Plan.begin();
      i != Plan.end(); ++i) {

    // setup the FSM for the sequences in this part
    long unsigned Slen = 1;
    for(min_index_list::const_iterator s = i->begin();
	s != i->end(); ++s) {
      Slen += s->stop - s->start;
    }

    char *S = new char[Slen];

    FILE_POSITION_TYPE pos = 0;
    for(min_index_list::const_iterator s = i->begin();
	s != i->end(); ++s) {
      min_sequence tmp(*s);
      tmp.chars = S+pos;
      tmp.sload(fwd);
      pos += tmp.stop - tmp.start;
    }

    CharMap.map(MAP1,S,S+Slen);

    ++parti;
    if (VERBOSE>1) {
      TOC(stdverb) fprintf(stdverb,": building fsm for part %d\n",parti);
      fflush(stdverb);
    }
    XSpaceFSM FSM(S,Slen,mersize);

    if (VERBOSE>1) {
      TOC(stdverb) fprintf(stdverb,": initializing fsm for part %d\n",parti);
      fflush(stdverb);
    }
    FSM.initialize();
    FSM.selfstream();

    if (VERBOSE>1) {
      TOC(stdverb) fprintf(stdverb,": streaming for part %d\n",parti);
      fflush(stdverb);
    }
    int partj = 0;
    // stream the FWD/REV versions of the other parts through
    for(list<min_index_list>::const_iterator j = Plan.begin();
	j != Plan.end(); ++j) {

      for(min_index_list::const_iterator s = j->begin();
	  s != j->end(); ++s) {

	// stream the RC part through
	if (rev) {
	  if (s->seek(rev) <0) {
	    fprintf(stderr,"main: unable to seek in rev\n");
	  }
	  FSM.stream(MAP2,rev,s->stop - s->start+1);
	}	  

	if (j != i) {
	  
	  // stream this part through
	  if (s->seek(fwd) <0) {
	    fprintf(stderr,"main: unable to seek in fwd\n");
	  }
	  FSM.stream(MAP2,fwd,s->stop - s->start+1);
	}
      }
      ++partj;
      if (VERBOSE>1 && parti != partj) {
	TOC(stdverb) fprintf(stdverb,": done %d vs %d\n",parti,partj);
	fflush(stdverb);
      }
    }

    if (VERBOSE>1) {
      TOC(stdverb) fprintf(stdverb,": finalizing scores\n");
      fflush(stdverb);
    }
    // FSM.finalize();

    // the Zscores for this part have now been computed
    // Write them
    if (VERBOSE>1) {
      TOC(stdverb) fprintf(stdverb,": written out the scores for part %d\n",parti);
      fflush(stdverb);
    }
    FSM.output(o,(FILE_POSITION_TYPE)(i->begin()->start));
    delete[] S;
    fflush(o);
  }
  
  if (VERBOSE) {
    TOC(stdverb) fprintf(stdverb,": written out the scores\n");
  }

  fputc('\n',o);

  fclose(o);
  fclose(idx);
  fclose(fwd);
  if (rev) fclose(rev);
}
