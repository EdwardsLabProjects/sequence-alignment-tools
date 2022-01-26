
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
  ostream *out;
  bool delout;
  bool verbose;
  int mersize;
  int iterations;
  Options(int argc, char *argv[]);
  ~Options();
  void usage(char *msg=NULL);
};

Options::Options(int argc, char *argv[]) : out(&std::cout) {
  signed char c;
  optarg = NULL;
  graphfile = "";
  sequencefile = "";
  eos_char = '$';
  delout = false;
  verbose = false;
  mersize = -1;
  iterations = 10000;
  while ((c = getopt(argc, argv, "g:i:E:o:k:c:vh")) != -1)
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
    case 'k':
      mersize = atoi(optarg);
      break; 
    case 'c':
      iterations = atoi(optarg);
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
Options::usage(char *message) {

  if (message != NULL && strlen(message) > 0) {
    cerr << message << endl;
    cerr << endl;
  }
  cerr << "Usage: solid_assembly [options]\n\n";
  cerr << "Options: \n";
  cerr << "  -g <graph-file> Word graph file. Required.\n";
  cerr << "  -i <fasta-file> Fasta file. Required.\n";
  cerr << "  -k <int>        Length of node sequence, if fixed\n";
  cerr << "  -c <int>        Widget removal cycles\n";
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
  
  word_graph g;
  g.read(opt.graphfile,opt.mersize,0,1,true);

  FastaFile<Lazy_Header_SI>* ff;
  fasta_file_seq_params ffp;
  ffp.upper_case = true;
  ffp.eos_char = opt.eos_char;
  ffp.eos_start = true;
  ff = pick_fasta_file<Lazy_Header_SI>(opt.sequencefile.c_str(),0,
				       /*memmap=*/true,/*load_index=*/true,
				       ffp,opt.verbose);

  g.print_stats();

  int i=0;
  while (i < opt.iterations) {
    if (!g.peel_edges(*ff,opt.mersize,opt.eos_char)) {
      break;
    }
    i++;
  }

  // g.remove_trivial_nodes(*ff,opt.verbose);
  // g.print_stats();

  g.writetrivialpaths(*opt.out,*ff,opt.eos_char,opt.verbose);

}

