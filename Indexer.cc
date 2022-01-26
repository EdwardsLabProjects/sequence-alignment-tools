#include <iostream>
using namespace std;

#include "rl_time.h"
#include "rl_charmap.h"
static charmap CharMap;

#include "rl_index.h"


static unsigned VERBOSE=0;
static char *oname = NULL;
static char *iname = NULL;
static char *fname = NULL;
static char *rname = NULL;
static char *whoami;

void usage(FILE *f) {
  fprintf(f,
	  "%s: [-v] [-i index file] [-f fwd file] [-r revc file] file.fasta\n",
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
        case 'i':
          oname = strdup(OPT_ARG);
          goto loopout;
        case 'f':
          fname = strdup(OPT_ARG);
          goto loopout;
        case 'r':
          rname = strdup(OPT_ARG);
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
      if (!iname) {iname = strdup(*argv); ++argv;}
      else
	errflg++;
    }
  }
  if (!iname) ++errflg;
  return errflg;
}

int main(int argc,const char **argv) {
  if (get_args(argc,argv)) {
    usage(stderr);
    exit(1);
  }

  TIC

  index_list li;

  FILE *f = fopen(iname,"r");
  if (!f) {
    fprintf(stderr,"%s:main: unable to open %s\n",whoami,iname);
    exit(31);
  }
  FILE *o = stdout;
  if (oname) {
    if (!(o = fopen(oname,"w"))) {
      fprintf(stderr,"%s:main: unable to open %s\n",whoami,oname);
      exit(32);
    }
  }
  FILE *fwd = NULL;
  if (fname) {
    if (!(fwd = fopen(fname,"w"))) {
      fprintf(stderr,"%s:main: unable to open %s\n",whoami,fname);
      exit(33);
    }
  }
  FILE *rev = NULL;
  if (rname) {
    if (!(rev = fopen(rname,"w"))) {
      fprintf(stderr,"%s:main: unable to open %s\n",whoami,rname);
      exit(34);
    }
  }

  if (li.iload_fasta(f) < 0) {
    fprintf(stderr,"%s:main: trouble indexing %s\n",whoami,iname);
    exit(35);
  }
  if (li.isave(o) < 0) {
    fprintf(stderr,"%s:main: trouble writing %s\n",whoami,
	    oname ? oname : "stdout" );
    exit(36);
  }

  if (VERBOSE) {
    TOC(stderr) fprintf(stderr,"written out the index\n");
  }
    
    
  if (fname || rname) {
    for(index_list::const_iterator i=li.begin(); i != li.end(); ++i) {
      sequence tmp(*i);

      tmp.chars = new char[ tmp.stop - tmp.start + 1];
      tmp.chars[0] = CharMap.term1();
      tmp.chars[tmp.stop - tmp.start] = CharMap.term1();

      tmp.sload_fasta(f);
      CharMap.map(CharMap.canonical(),
		  tmp.chars+1 ,
		  tmp.chars+tmp.stop-tmp.start);

      if (fname) {
	tmp.ssave(fwd);
      }
      if (rname) {
	CharMap.revmap(CharMap.complement(),
		       tmp.chars+1 ,
		       tmp.chars+tmp.stop-tmp.start);
	tmp.ssave(rev);
      }

      delete[] tmp.chars;
    }
  }

  if (VERBOSE) {
    TOC(stderr) fprintf(stderr,"written out the sequences\n");
  }

  fclose(f);
  if (oname) fclose(o);
  if (fname) fclose(fwd);
  if (rname) fclose(rev);
}
