#include <cstdio>
#include <iostream>
using namespace std;

#include "rl_suffix_tree.h"
#include "rl_time.h"
#include "rl_charmap.h"
static charmap CharMap;

#include "rl_index.h"

static unsigned VERBOSE=0;
static char *oname = NULL;
static char *iname = NULL;
static char *fname = NULL;
static unsigned MerSize=10;
static char *whoami;
static bool Amino=false;
static bool UC=false;
static bool AllSeqs=false;

void usage(FILE *f) {
  fprintf(f,
	  "%s: [-v] [-h] [-m mer_size] [-a] [-A] [-U] -o output -i index_file -f fwd_file\n",
	  whoami);
  fprintf(f, "\t-m is the mersize for the word graph\n");
  fprintf(f, "\t-a means to do the analysis on all sequences simultaneously\n");
  fprintf(f, "\t-A assume amino acid files\n");
  fprintf(f, "\t-U assume unrestricted (any letter) amino acid files\n");
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
        case 'a':
          AllSeqs=true;
          goto loopin;
        case 'A':
	  Amino = true;
          goto loopin;
        case 'U':
	  UC = true;
          goto loopin;
        case 'm':
          MerSize = atoi(OPT_ARG);
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
  if (!iname) ++errflg;
  if (!fname) ++errflg;
  return errflg;
}

class SufTree: public suffix_tree<basic_node,basic_suffix> {
public:
  SufTree(const char *st,unsigned len,char t='\0'):
    suffix_tree<basic_node,basic_suffix>(st,len,t) {}
  void process(FILE *o,st_index n) const;
  bool interesting(st_index n,char cl,char cr,st_index par) const;
  void output(FILE *o,st_index n,st_index par) const;
};

/** debug **
static bool debug=false;
**/
void SufTree::process(FILE *o,st_index n) const {
  if (!n.is_leaf()) {
    if (node(n).length() >= MerSize) {
      if (interesting(n,str(n)[-1],str(n)[MerSize],n)) {
	output(o,n,n);
	fputc('\n',o);
      }
    }
    else {
      n = first_child(n);
      for(n.set_good(); n.is_good(); n = next_child(n)) process(o,n);
    }
  }
}

bool SufTree::interesting(st_index n,char cl,char cr,st_index par) const {

  if (str(n)[-1]!=cl || str(n)[-1]==str()[0]
      ||
      str(n)[MerSize]!=cr || str(n)[MerSize]==str()[0])
    return true;

  if (!n.is_leaf()) {
    n = first_child(n);
    for(n.set_good(); n.is_good(); n = next_child(n)) {
      if (interesting(n,cl,cr,par)) return true;
    }
  }

  return false;
}

void SufTree::output(FILE *o,st_index n,st_index par) const {
  if (!n.is_leaf()) {
    n = first_child(n);
    for(n.set_good(); n.is_good(); n = next_child(n)) output(o,n,par);
  }
  else {
    /** debug **
    debug = (str(n)-str()+MerSize==13709759) ||(str(n)-str()+MerSize==23250365);
    if (debug) {
      unsigned offset = str(n)-str()+MerSize-1;
      fprintf(stderr,"here I am, printing under %d %u\n",par.is_leaf(),par.index());
      for(unsigned i=0; i < MerSize; ++i) fputc(str()[offset-MerSize+1+i],stderr);
      fputc('\n',stderr);
    }
    **/
    fprintf(o," %u.%c",str(n)+MerSize-str(),str(n)[MerSize]);
  }
}

int main(int argc,const char **argv) {
  char MAP[256];
  TIC

  if (get_args(argc,argv)) {
    usage(stderr);
    exit(1);
  }
  FILE *stdverb = stdout;

  index_list li;

  FILE *o=stdout;
  if (oname) {
    if (!(o = fopen(oname,"w"))) {
      fprintf(stderr,"%s:main: unable to open %s\n",whoami,oname);
      exit(101);
    }
    stdverb = stderr;
  }
  FILE *fwd,*idx;
  if (!(idx = fopen(iname,"r"))) {
    fprintf(stderr,"%s:main: unable to open %s\n",whoami,iname);
    exit(102);
  }
  if (!(fwd = fopen(fname,"r"))) {
    fprintf(stderr,"%s:main: unable to open %s\n",whoami,fname);
    exit(103);
  }
  if (li.iload(idx)) {
    fprintf(stderr,"%s:main: trouble loading index from %s\n",whoami,iname);
    exit(105);
  }
  if (VERBOSE) {
    TOC(stdverb) fprintf(stdverb,": loaded the index\n");
  }

  if (Amino) {
    for(int i=0; i < 256; ++i) {
      MAP[i] = CharMap.aminoacid()[i];
    }
    // turn X's into mismatches
    for(int i=0; i < 256; ++i) if (MAP[i] == 'X') MAP[i] = CharMap.term1();
  } 
  else if (UC) {
    for(int i=0; i < 256; ++i) {
      MAP[i] = CharMap.uppercase()[i];
    }
    // turn non-letter's into mismatches
    for(int i=0; i < 256; ++i) if (MAP[i] == CharMap.term3()) MAP[i] = CharMap.term1();
  }
  else {
    for(int i=0; i < 256; ++i) {
      MAP[i] = CharMap.canonical()[i];
    }
    // turn n's into mismatches
    for(int i=0; i < 256; ++i) if (MAP[i] == 'N') MAP[i] = CharMap.term1();
  }

  if (AllSeqs) {
    if (VERBOSE>1) {
      TOC(stdverb) fprintf(stdverb,": doing all sequences in file\n");
      fflush(stdverb);
    }
    unsigned Slen = 1;
    for(index_list::const_iterator i = li.begin() ;
	i != li.end(); ++i) {
      Slen += i->stop - i->start;
    }
    char *S = new char[Slen];
    unsigned pos=0;
    for(index_list::const_iterator i = li.begin() ; i != li.end(); ++i) {
      sequence seq(*i);
      seq.chars = S+pos;
      seq.sload(fwd);
      pos += i->stop - i->start;
    }

    CharMap.map(MAP,S,S+Slen);

    SufTree T(S,Slen,S[0]);
    T.build();

    if (VERBOSE>1) {
      TOC(stdverb) fprintf(stdverb,": tree built\n");
      fflush(stdverb);
    }

    // fake node to give the term character
    fprintf(o," %llu.%c\n",0,S[0]);

    // fake nodes to mark sequence boundaries
    fprintf(o," %llu.%c %llu.%c\n",1,S[1],Slen,S[1]);
    
    pos = 1;
    for(index_list::const_iterator i = li.begin() ; i != li.end(); ++i) {
      pos += i->stop - i->start;
      if (pos < Slen) fprintf(o," %llu.%c\n",pos,S[pos]);
    }

    T.process(o,T.root());

    // empty set terminates output
    fprintf(o,"\n");

    T.clear();

    if (VERBOSE>1) {
      TOC(stdverb) fprintf(stdverb,": xspaces dumped\n");
      fflush(stdverb);
    }

    delete[] S;
  }
  else {
    int parti=0;
    for(index_list::const_iterator i = li.begin() ;
	i != li.end(); ++i) {

      // setup the FSM for the sequences in this part
      unsigned Slen = 1+i->stop - i->start;

      char *S = new char[Slen];
      
      sequence seq(*i);
      seq.chars = S;
      seq.sload(fwd);

      CharMap.map(MAP,S,S+Slen);

      /** debug **
      {
	unsigned offset;
	offset=13709759;
	unsigned MerSize=128;
	for(unsigned i=0; i < MerSize; ++i) fputc(S[offset-MerSize+1+i],stderr);	fputc('\n',stderr);
	offset=23250365;
	for(unsigned i=0; i < MerSize; ++i) fputc(S[offset-MerSize+1+i],stderr);	fputc('\n',stderr);
	
	exit(0);
      }
      **/

      SufTree T(S,Slen,S[0]);
      T.build();

      if (VERBOSE>1) {
	TOC(stdverb) fprintf(stdverb,": tree built for part %d\n",++parti);
	fflush(stdverb);
      }

      // fake node to give the term character
      fprintf(o," %llu.%c\n",0,S[0]);

      // fake nodes to mark sequence boundaries
      fprintf(o," %llu.%c %u.%c\n",unsigned(1),S[1],Slen,S[1]);

      T.process(o,T.root());

      // empty set terminates output
      fprintf(o,"\n");

      T.clear();

      if (VERBOSE>1) {
	TOC(stdverb) fprintf(stdverb,": xspaces dumped for part %d\n",parti);
	fflush(stdverb);
      }

      delete[] S;
    }
  }
  
  if (VERBOSE) {
    TOC(stdverb) fprintf(stdverb,": written out the Xspaces\n");
  }

  if (oname) fclose(o);
  fclose(idx);
  fclose(fwd);
}
