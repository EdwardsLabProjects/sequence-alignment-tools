
#include <iostream>
#include <sstream>
#include <list>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <assert.h>
#include <string.h>
#ifdef __sun
#include <strings.h>
#endif

#include "select.t"
#include "fasta_io.t"
#include "word_graph.h"

struct Options {
  std::string graphfile;
  std::string sequencefile;
  char eos_char;
  bool reuse;
  bool reuseopt;
  ostream *out;
  bool delout;
  bool verbose;
  bool peel;
  int mersize;
  int ctspec;
  int ctsign;
  long maxlabel;
  Options(int argc, char *argv[]);
  ~Options();
  void usage(const char *msg=NULL);
};

Options::Options(int argc, char *argv[]) : out(&std::cout) {
  signed char c;
  optarg = NULL;
  graphfile = "";
  sequencefile = "";
  eos_char = '$';
  delout = false;
  verbose = false;
  reuse = false;
  reuseopt = false;
  mersize = -1;
  ctsign = 0;
  ctspec = 0;
  maxlabel = 10000;
  peel = false;
  while ((c = getopt(argc, argv, "g:i:E:o:k:C:PrRvhM:")) != -1)
    switch (c) {
    case 'g':
      graphfile = optarg;
      break;
    case 'i':
      sequencefile = optarg;
      break;
    case 'o':
      if (delout) {
	delete out;
      }
      out = new ofstream(optarg);
      delout = true;
      break;
    case 'E': {
      int eos_chari;
      if (!sscanf(optarg,"%i",&eos_chari)) {
	usage("Invalid end-of-sequence specification.\n");
      }
      eos_char = (char)eos_chari;
    }
      break;
    case 'v':
      verbose = true;
      break; 
    case 'P':
      peel = true;
      break; 
    case 'k':
      mersize = atoi(optarg);
      break; 
    case 'M':
      maxlabel = atol(optarg);
      break; 
    case 'C':
      ctspec = atoi(optarg);
      if (ctspec < 0) {
	ctspec = -ctspec;
	ctsign = -1;
      } else {
	if (index(optarg,'+')!=NULL) {
	  ctsign = 1;
	} else {
	  ctsign = 0;
	}
      } 
      break; 
    case 'r':
      reuse = true;
      reuseopt = false;
      break; 
    case 'R':
      reuse = true;
      reuseopt = true;
      break; 
    case 'h':
    default :
      usage();
    }
    if ((graphfile == "" || sequencefile == "")) usage();
}

Options::~Options() {
  if (delout) {
    delete out;
  }
}

void 
Options::usage(const char *message) {

  if (message != NULL && strlen(message) > 0) {
    cerr << message << endl;
    cerr << endl;
  }
  cerr << "Usage: cannon_csbh_graph [options]\n\n";
  cerr << "Options: \n";
  cerr << "  -g <graph-file> Word graph file. Required.\n";
  cerr << "  -i <fasta-file> Fasta file. Required.\n";
  cerr << "  -k <int>        Length of node sequence, if fixed\n";
  cerr << "  -C (c|+c|-c)    Keep only those edges with count exactly c, more than c or\n";
  cerr << "                  less than c.\n";
  cerr << "  -r              Reuse edges. Heurtistic reuse.\n";
  cerr << "  -R              Reuse edges. Optimal reuse.\n";
  cerr << "  -P              Peel edges.\n";
  cerr << "  -E <int>        End-of-sequence character. Default: \'$\'\n";
  cerr << "  -o <out-file>   Sequence output file. Default: Output to stdout.\n";
  cerr << "  -v              Verbose.\n";
  cerr << "  -h              Help.\n";
  exit(1);
}

int
main(int argc, char **argv) {
  
#if defined(PROFILE)
  set_profile_signal_handler();
#endif

  Options opt(argc,argv);
  if (opt.verbose) {
    timestamp("Begin.");
  }
  
  word_graph g(opt.maxlabel);
  g.read(opt.graphfile,opt.mersize,opt.ctspec,opt.ctsign);
  if (opt.verbose) {
    timestamp("Finished reading word_graph");
    g.print_stats();
  }
  if (opt.verbose) {
    timestamp("Finished computing stats");
  }

  FastaFile<Lazy_Header_SI>* ff;
  fasta_file_seq_params ffp;
  ffp.upper_case = true;
  ffp.eos_char = opt.eos_char;
  ffp.eos_start = true;
  ff = pick_fasta_file<Lazy_Header_SI>(opt.sequencefile.c_str(),0,
				       /*memmap=*/true,/*load_index=*/true,
				       ffp,opt.verbose);

  // checkpoint;
  // g.dump(cerr,*ff);

  if (opt.peel) {
    timestamp("Starting edge peel.");
    int i=0;
    while (i < 50) {
      if (!g.peel_edges(*ff,opt.mersize,opt.eos_char)) {
	break;
      }
      i++;
    }
    timestamp("Edge peel done.");
  }

  // checkpoint;
  // g.dump(cerr,*ff);


  if (opt.reuse) {
    if (opt.verbose) {
      timestamp("Adding shortcut edges");
    }
    g.find_joiners(opt.reuseopt,opt.verbose);
    if (opt.verbose) {
      timestamp("Finished adding shortcut edges");
    }
  }
  if (opt.verbose) {
    timestamp("Balancing word_graph");
  }
  g.balance_nodes(opt.eos_char,opt.verbose);
  if (opt.verbose) {
    timestamp("Finished balancing word_graph");
  }

  // g.dump(std::cerr,*ff,false);
  if (opt.verbose) {
    timestamp("Writing word_graph sequence");
  }
  g.writeseq(*opt.out,*ff,opt.eos_char,opt.verbose);
  if (opt.verbose) {
    timestamp("Finished writing word_graph sequence");
  }
  if (opt.verbose) {
    timestamp("Done.");
  }
}

