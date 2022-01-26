
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
  std::string annotatefile;
  char eos_char;
  ostream *out;
  bool delout;
  bool verbose;
  int transform;
  int format;
  int mersize;
  Options(int argc, char *argv[]);
  ~Options();
  void usage(char *msg=NULL);
};

Options::Options(int argc, char *argv[]) : out(&std::cout) {
  signed char c;
  optarg = NULL;
  graphfile = "";
  sequencefile = "";
  annotatefile = "";
  eos_char = '$';
  delout = false;
  verbose = false;
  transform = 0;
  mersize = -1;
  format = 0;
  while ((c = getopt(argc, argv, "g:i:a:E:o:k:t:f:vh")) != -1)
    switch (c) {
    case 'g':
      graphfile = optarg;
      break;
    case 'i':
      sequencefile = optarg;
      break;
    case 'a':
      annotatefile = optarg;
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
    case 't':
      transform = atoi(optarg);
      break; 
    case 'f':
      format = atoi(optarg);
      break; 
    case 'k':
      mersize = atoi(optarg);
      break; 
    case 'h':
    default :
      usage();
    }
    if ((graphfile == "" || sequencefile == "")) usage();
    if (annotatefile == "") {
      annotatefile = sequencefile;
    }
}

Options::~Options() {
  if (delout) {
    delete out;
  }
}

void 
Options::usage(char *message) {

  if (message != NULL && strlen(message) > 0) {
    cerr << message << endl;
    cerr << endl;
  }
  cerr << "Usage: csbh_annotate [options]\n\n";
  cerr << "Options: \n";
  cerr << "  -g <graph-file> Word graph file. Required.\n";
  cerr << "  -i <fasta-file> Fasta file for graph. Required.\n";
  cerr << "  -a <fasta-file> Fasta file to annotate. Optional.\n";
  cerr << "  -k <int>        Length of node sequence, if fixed\n";
  cerr << "  -E <int>        End-of-sequence character. Default: \'$\'\n";
  cerr << "  -t <transform>  Output transform for counts. Default: 0;\n";
  cerr << "                  0: no transform; 1: log2(c)+1; 2: greater than 1\n";
  cerr << "  -f <format>     Output format. Default: 0;\n";
  cerr << "                  0: Fasta alpha; 1: UCSC WIG (non-unique only);\n";
  cerr << "                  2: Fasta sequence; 3: UCSC WIG (all counts)\n";
  cerr << "  -o <out-file>   Sequence output file. Default: Output to stdout.\n";
  cerr << "  -v              Verbose.\n";
  cerr << "  -h              Help.\n";
  exit(1);
}

int
main(int argc, char **argv) {
  
  Options opt(argc,argv);
  if (opt.verbose) {
    timestamp("Begin.");
  }
  
  word_graph g;
  g.read(opt.graphfile,opt.mersize,0,0,true);
  if (opt.verbose) {
    timestamp("Finished reading word_graph");
    g.print_stats();
    timestamp("Finished computing stats");
  }

  FastaFile<Lazy_Header_SI>* ffs,*ffa;
  fasta_file_seq_params ffp;
  ffp.upper_case = true;
  ffp.eos_char = opt.eos_char;
  ffp.eos_start = true;
  ffs = pick_fasta_file<Lazy_Header_SI>(opt.sequencefile.c_str(),0,
					/*memmap=*/true,
					/*load_index=*/true,
					ffp,opt.verbose);
  ffa = pick_fasta_file<Lazy_Header_SI>(opt.annotatefile.c_str(),0,
					/*memmap=*/true,
					/*load_index=*/true,
					ffp,opt.verbose);

  if (!g.check_out_edges(*opt.out,*ffs,opt.verbose)) {
    exit(1);
  }

  if (opt.verbose) {
    timestamp("Annotating sequence");
  }
  g.annotateseq(*opt.out,*ffs,*ffa,opt.eos_char,opt.transform,opt.format,opt.verbose);
  if (opt.verbose) {
    timestamp("Finished annotating sequence");
  }
  if (opt.verbose) {
    timestamp("Done.");
  }
}

