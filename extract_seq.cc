
#include "types.h"
#include "char_io.h"
#include "char_io.t"
#include "fasta_io.h"
#include "fasta_io.t"
#include "select.t"
#include <assert.h>

#include <list>

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

char release_tag[] = "$Name:  $";

void usage(char *message=NULL) {
  if (message != NULL && strlen(message) > 0) {
    cerr << message << endl;
    cerr << endl;
  }
  cerr << "Usage: extract_seq [options] \n\n";
  cerr << "Options: \n";
  cerr << "  -i <sequence-file>\n";
  cerr << "  -o <output-fasta>\n";
  cerr << "  -A <pos-file>   Line based alignment records, format:\n";
  cerr << "                  <id-string> <fasta-index> <start-pos> <length>\n";
  cerr << "                  <id-string> need not be unique.\n";   
  cerr << "                  <fasta-index> is 0,1,2,...\n";   
  cerr << "                  <start-pos> is space-based.\n";   
  cerr << "                  Fields are white space separated.\n";   
  cerr << "                  \"-\" indicates standard input.\n";
  cerr << "  -n              Output sequence between alignments.\n";
  cerr << "  -I              Include \"ends\" of sequence for between alignments.\n";
  cerr << "                  Default: false.\n";
  cerr << "  -v              Verbose.\n";
  cerr << "  -h\n";
  cerr << "\n";
  exit(1);
}

struct Options {
  Options() : out(&cout) {};
  std::string atac_file;
  std::string seq_file;
  bool notin;
  bool includeends;
  bool verbose;
  ostream *out;
  bool delout;
  char eos_char;
};

void options(int argc, char *argv[], Options & opt) {
  signed char c;
  optarg = NULL;
  opt.notin = false;
  opt.includeends = false;
  opt.verbose = false;
  opt.delout = false;
  opt.eos_char = '\n';
  while ((c = getopt(argc, argv, "E:A:i:o:nIvh")) != -1) {
    switch (c) {
    case 'A':
      opt.atac_file = optarg;
      break;
    case 'i':
      opt.seq_file = optarg;
      break;
    case 'n':
      opt.notin = true;
      break;
    case 'I':
      opt.includeends = true;
      break;
    case 'E': {
      int eos;
      if (!sscanf(optarg,"%i",&eos)) {
	usage("Invalid end-of-sequence specification.\n");
      }
      opt.eos_char = (char)eos;
    }
    break;
    case 'v':
      opt.verbose = true;
      break;
    case 'o':
      if (opt.delout) {
	delete opt.out;
      }
      opt.out = new ofstream(optarg,(ios::out|ios::app|ios::ate));
      opt.delout = true;
      break;
    case 'h':
    default :
      usage();
    }
  }
  if (opt.verbose) return;
  if (opt.atac_file == "") usage();
  if (opt.seq_file == "") usage();
}

struct match {
  match(long signed int f, 
	FILE_POSITION_TYPE s, FILE_POSITION_TYPE l)
    : fasta_entry(f), start(s), length(l) {};
  long signed int fasta_entry;
  FILE_POSITION_TYPE start;
  FILE_POSITION_TYPE length;
};

bool operator<(match const & a, match const & b) {
  if (a.fasta_entry < b.fasta_entry) {
    return true;
  } else if (a.fasta_entry > b.fasta_entry) {
    return false;
  }
  if (a.start < b.start) {
    return true;
  } else if (a.start > b.start) {
    return false;
  }
  if (a.length < b.length) {
    return true;
  } else if (a.length > b.length) {
    return false;
  }
  return false;
}

int main(int argc,char *argv[]) {

  Options opt;
  options(argc,argv,opt);

  if (opt.verbose) {
    ostrstream ss;
    ss << "Release Tag: " << release_tag << ends;
    std::string v(ss.str());
    timestamp(v.c_str());
  }

  FastaFile<Lazy_Header_SI>* ff;
  fasta_file_seq_params ffp;
  ffp.upper_case = true;
  ffp.eos_char = opt.eos_char;
  ffp.check_params = true;
  ff = pick_fasta_file<Lazy_Header_SI>(opt.seq_file,
				       /*ffselect=*/0,/*memmap=*/true,
				       /*alignments=*/true,ffp,opt.verbose);
  
  bool delete_stream = false;
  istream * ifs(&cin);
  if (opt.atac_file != "-") {
    ifs = new ifstream(opt.atac_file.c_str());
    delete_stream = true;
  }

  const int bufsize=8192;
  char buffer[bufsize];
  
  long signed int fe;
  FILE_POSITION_TYPE s;
  FILE_POSITION_TYPE l;
  string toeol;
  match m(0,0,0);
  match *it = &m;

  if (!opt.notin) {
    
    FILE_POSITION_TYPE abspos;
    while ((*ifs)) {
      (*ifs) >> fe >> s >> l;
      getline(*ifs,toeol);
      if (l > 0) {
	if (fe >= 0) {
	  m.fasta_entry = fe;
	  m.start = s;
	  m.length = l;
	} else {
	  Lazy_Header_SI const & hi(ff->get_header_data(s));
	  s = ff->get_seq_pos(s+1)-1;
	  fe = hi.index()-1;
	  m.fasta_entry = fe;
	  m.start = s;
	  m.length = l;
	}
      }
      bool b=true;
      if (it->fasta_entry >= 0) {
	b = ff->fasta_pos(it->fasta_entry,it->start);
      } else {
	ff->pos(it->start);
      }
      if (b) {
	abspos = ff->pos();
	if (ff->is_subseq(abspos,abspos+it->length)) {
	  ff->pos(abspos);
	  int i,j;
	  for (i=0,j=0;i<it->length;i++) {
	    buffer[j++] = ff->getch();
	    if (j >= bufsize) {
	      opt.out->write(buffer,j);
	      j=0;
	    }
	  }
	  if (j>0) {
	    opt.out->write(buffer,j);
	  }
	  *opt.out << opt.eos_char;
	} else {
	  if (opt.verbose) 
	    timestamp("Warning: Sequence is not in a single fasta entry.");
	}
      } else {
	if (opt.verbose) 
	  timestamp("Warning: Can't set Fasta file postion.");
      }
      // delete it;
    }

  } else {
    /*
    std::vector<match>::iterator it=mymatches.begin();
    FILE_POSITION_TYPE abspos,endabspos;
    int current_fasta_entry=0;
    string current_run=it->run;
    FILE_POSITION_TYPE current_start_position=0;
    FILE_POSITION_TYPE current_end_position=0;
    FILE_POSITION_TYPE current_len;
    while (it != mymatches.end()) {
      current_start_position=0;
      while (it != mymatches.end() && it->fasta_entry==current_fasta_entry) {
	current_end_position=it->start; 
	if (current_start_position <= current_end_position &&  
	    !(!opt.includeends && 
	      ((it->run != current_run) || current_start_position==0))) {
	  current_len=current_end_position-current_start_position;
	  if (ff->fasta_pos(current_fasta_entry,current_start_position)) {
	    abspos = ff->pos();
	    if (ff->is_subseq(abspos,abspos+current_len)) {
	      Lazy_Header_SI const & hi(ff->get_header_data(abspos));
	      *opt.out << ">" << hi.header();
	      if (current_start_position>0) {
		*opt.out << " /" << runword << "_before=" << "{" << (it-1)->id << "}"
		     << (it-1)->fasta_entry << ":"
		     << (it-1)->start << "-"
		     << (it-1)->start + (it-1)->length << "("
		     << (it-1)->sense_seq << ","
		     << (it-1)->sense_comp << ")";
	      } 
	      *opt.out << " /" << runword << "_after=" << "{" << it->id << "}" 
		   << it->fasta_entry << ":"
		   << it->start << "-"
		   << it->start + it->length << "("
		   << it->sense_seq << ","
		   << it->sense_comp << ")";
	      if (current_start_position>0) {
		*opt.out << " /cannonical_pair_id={";
		if ((it-1)->id <= it->id) {
		  *opt.out << (it-1)->id << "," << it->id;
		} else {
		  *opt.out << it->id << "," << (it-1)->id;
		}
		*opt.out << "}";
	      }
	      *opt.out << " /between_" << runword << "="
		   << current_fasta_entry << ":"
		   << current_start_position << "-"
		   << current_end_position;
	      *opt.out << " /length=" 
		   << current_len;
	      long unsigned int ncount=0;
	      long unsigned int maxcontign=0,curcontign=0;
	      for (i=0;i<current_len;i++) {
		if (ff->eof()) break; 
		if (ff->getch() == 'N') {
		  ncount++;
		  curcontign++;
		} else {
		  if (curcontign > maxcontign) {
		    maxcontign = curcontign;
		  }
		  curcontign = 0;
		}
	      }
	      if (curcontign > maxcontign) {
		maxcontign = curcontign;
	      }
	      *opt.out << " /Ns=" << ncount  
		       << " /nonNs=" << current_len - ncount
		       << " /maxContigN=" << maxcontign << endl;
	      if (!opt.headersonly) {
		ff->pos(abspos);
		for (i=0,j=0;i<current_len;i++) {
		  buffer[j++] = ff->getch();
		  if (i%60==59) {
		    buffer[j++] = '\n';
		  }
		  if (j >= 8190) {
		    opt.out->write(buffer,j);
		    j=0;
		  }
		}
		if (j>0) {
		  opt.out->write(buffer,j);
		}
		if (j==0 || (j > 0 && buffer[j-1] != '\n')) {
		  *opt.out << endl;
		}
	      }
	    } else {
	      if (opt.verbose) 
		timestamp("Warning: Sequence is not in a single fasta entry.");
	    }
	  } else {
	    if (opt.verbose) 
	      timestamp("Warning: Can''t set Fasta file postion.");
	  }
	} else {
	  if (current_start_position > current_end_position &&
	      (opt.includeends || it->run == current_run)) {
	    if (opt.verbose && opt.nooverlap) 
	      timestamp("Overlapping alignments observed.");
	  }
	}
	current_start_position=(it->start+it->length);
	current_run = it->run;
	++it;
      }
      if (ff->fasta_pos(current_fasta_entry,current_start_position)) {
	abspos = ff->pos();
	ff->fasta_pos(current_fasta_entry+1,0);
	endabspos = ff->pos()-1 // eos char;
	ff->pos(abspos);
	if (endabspos > abspos && opt.includeends) {
	  current_len = endabspos-abspos;
	  current_end_position = current_start_position+current_len;
	  Lazy_Header_SI const & hi(ff->get_header_data(abspos));
	  *opt.out << ">" << hi.header();
	  if (current_start_position>0) {
	    *opt.out << " /" << runword << "_before=" << "{" << (it-1)->id << "}" 
		 << (it-1)->fasta_entry << ":"
		 << (it-1)->start << "-"
		 << (it-1)->start + (it-1)->length << "("
		 << (it-1)->sense_seq << ","
		 << (it-1)->sense_comp << ")";
	  } 
	  *opt.out << " /between_" << runword << "="
	       << current_fasta_entry << ":"
	       << current_start_position << "-"
	       << current_end_position;
	  *opt.out << " /length=" 
	       << current_len;
	  long unsigned int ncount=0;
	  long unsigned int maxcontign=0,curcontign=0;
	  i=0;
	  char ch=ff->getch();
	  while (ch!='\n' && i<current_len) {
	    if (ch == 'N') {
	      ncount++;
	      curcontign++;
	    } else {
	      if (curcontign > maxcontign) {
		maxcontign = curcontign;
	      }
	      curcontign = 0;
	    }
	    ch=ff->getch();
	    i++;
	  }
	  if (curcontign > maxcontign) {
	    maxcontign = curcontign;
	  }
	  *opt.out << " /Ns=" << ncount   
		   << " /nonNs=" << current_len - ncount
		   << " /maxContigN=" << maxcontign << endl;
	  if (!opt.headersonly) {
	    ff->pos(abspos);
	    i=0;j=0;
	    ch=ff->getch();
	    while (ch!='\n' && i<current_len) {
	      buffer[j++] = ch;
	      if (i%60==59) {
		buffer[j++] = '\n';
	      }
	      if (j >= 8190) {
		opt.out->write(buffer,j);
		j=0;
	      }
	      ch = ff->getch();
	      i++;
	    }
	    if (j>0) {
	      opt.out->write(buffer,j);
	    }
	    if (j==0 || (j > 0 && buffer[j-1] != '\n')) {
	      *opt.out << endl;
	    }
	  }
	}
	if (it != mymatches.end()) {
	  current_fasta_entry++;
	}
      } else {
	if (opt.verbose) 
	  timestamp("Warning: Can''t set Fasta file postion.");
      }
    }  
    */
  } 
  opt.out->flush();
  // *opt.out << "\nDone" << endl;
  if (delete_stream) {
    delete ifs;
  }
  if (opt.verbose) 
    timestamp("Pulling out sequence from fasta file...done.");
}
