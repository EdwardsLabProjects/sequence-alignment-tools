
#include "types.h"
#include "char_io.h"
#include "char_io.t"
#include "fasta_io.h"
#include "fasta_io.t"
#include "select.t"
#include <assert.h>

#include <list>
#include <vector>
#include <algorithm>
// #include "AtacMatch.h"
// #include "AtacComment.h"
// #include "AtacAttribute.h"
// #include "AtacDataset.h"

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

char release_tag[] = "$Name:  $";

void usage(char *message=NULL) {
  if (message != NULL && strlen(message) > 0) {
    cerr << message << endl;
    cerr << endl;
  }
  cerr << "Usage: atac_seq [options] \n\n";
  cerr << "Options: \n";
  cerr << "  -i <sequence-file>\n";
  cerr << "  -o <output-fasta>\n";
  /*   cerr << "\n";
  cerr << "  -a <ATAC-file>\n";
  cerr << "  -s (0|1)        Output ATAC sequence 0 or 1.\n";
  cerr << "  -m              Consider matches.\n";
  cerr << "  -r              Consider runs.\n";
  cerr << "  -t <char>       Consider type <char> alignments.\n";
  cerr << "                  Set to \'u\' by -m.\n";
  cerr << "                  Set to \'r\' by -r.\n";
  cerr << "\n"; */
  cerr << "  -A <pos-file>   Line based alignment records, format:\n";
  cerr << "                  <id-string> <fasta-index> <start-pos> <length>\n";
  cerr << "                  <id-string> need not be unique.\n";   
  cerr << "                  <fasta-index> is 0,1,2,...\n";   
  cerr << "                  <start-pos> is space-based.\n";   
  cerr << "                  Fields are white space separated.\n";   
  cerr << "                  \"-\" indicates standard input.\n";
  cerr << "\n";
  cerr << "  -n              Output sequence between matches or runs.\n";
  cerr << "  -I              Include \"ends\" of sequence for between runs.\n";
  cerr << "                  Default: false.\n";
  cerr << "                  Set to true by -r.\n"; 
  cerr << "                  Set to false by -m.\n";
  cerr << "  -O              Permit overlap in runs or matches.\n";
  cerr << "                  Default: false.\n";
  cerr << "  -e              Basic extract.\n";
  cerr << "  -H              Output headers only.\n";
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
  int seqno;
  std::string type;
  bool notin;
  bool nooverlap;
  bool includeends;
  bool verbose;
  bool atac_format;
  bool extract;
  bool headersonly;
  ostream *out;
  bool delout;
  char eos_char;
};

void options(int argc, char *argv[], Options & opt) {
  signed char c;
  optarg = NULL;
  opt.seqno = -1;
  opt.notin = false;
  opt.includeends = false;
  opt.nooverlap = true;
  opt.type = "";
  opt.verbose = false;
  opt.atac_format = true;
  opt.delout = false;
  opt.headersonly = false;
  opt.eos_char = '\n';
  opt.extract = false;
  while ((c = getopt(argc, argv, "E:A:i:s:o:Omrt:nHeIvh")) != -1) {
    switch (c) {
    case 'a':
      opt.atac_file = optarg;
      opt.atac_format = true;
      break;
    case 'A':
      opt.atac_file = optarg;
      opt.atac_format = false;
      break;
    case 'i':
      opt.seq_file = optarg;
      break;
    case 's':
      opt.seqno = atoi(optarg);
      break;
    case 'n':
      opt.notin = true;
      break;
    case 'I':
      opt.includeends = true;
    case 'm':
      opt.includeends = false;
      opt.type = "u";
      break;
    case 'r':
      opt.includeends = true;
      opt.type = "r";
      break;
    case 't':
      opt.type = std::string(optarg);
      break;
    case 'e':
      opt.extract = true;
      opt.nooverlap = false;
      break;
    case 'O':
      opt.nooverlap = false;
      break;
    case 'H':
      opt.headersonly = true;
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
  if (opt.atac_format && opt.seqno != 0 && opt.seqno != 1) usage();
  if (opt.atac_format && opt.type == "") usage();
}

struct match {
  match(std::string i, 
      long signed int f, std::string r, 
      FILE_POSITION_TYPE s, FILE_POSITION_TYPE l, 
      int ss, int sc)
    : id(i), 
    fasta_entry(f), run(r), 
    start(s), length(l), sense_seq(ss), sense_comp(sc) {};
  std::string id;
  long signed int fasta_entry;
  std::string run;
  FILE_POSITION_TYPE start;
  FILE_POSITION_TYPE length;
  int sense_seq;
  int sense_comp;
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
  if (a.sense_seq < b.sense_seq) {
    return true;
  } else if (a.sense_seq > b.sense_seq) {
    return false;
  }
  if (a.sense_comp < b.sense_comp) {
    return true;
  } else if (a.sense_comp > b.sense_comp) {
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

  if (opt.atac_file == "" || opt.seq_file == "") {
    exit(1);
  }

  FastaFile<Lazy_Header_SI>* ff;
  fasta_file_seq_params ffp;
  ffp.upper_case = true;
  ffp.eos_char = opt.eos_char;
  ffp.check_params = true;
  ff = pick_fasta_file<Lazy_Header_SI>(opt.seq_file,
				       /*ffselect=*/0,/*memmap=*/true,
				       /*alignments=*/true,ffp,opt.verbose);
  
  std::vector<match> mymatches;
    
  if (opt.atac_format) {
    /* 
    if (opt.verbose) 
      timestamp("Reading ATAC data...");

    AtacDataset *dataset = new AtacDataset();
    dataset->setInputFile(opt.atac_file);
    dataset->readRecords();
    
    if (opt.verbose) 
      timestamp("Reading ATAC data...done.");
    
    if (opt.verbose) 
      timestamp("Looking through runs & matches...");
    
    std::string typechar=opt.type;
    const vector<AtacMatch*> &matches = dataset->getMatches();
    vector<AtacMatch*>::const_iterator matit;
    // checkpoint;
    for (matit=matches.begin();matit!=matches.end();++matit) {
      AtacMatch *am;
      // checkpoint;
      if ((am=*matit) && am->isValid() && am->getTypeId() == typechar) {
	// checkpoint;
	AtacGenomicAxis *axis; 
	axis = dataset->getGenomicAxis(am->getGenomicAxisId(opt.seqno));
	// cerr << (void*) axis << endl;
	// cerr << axis->getOrdinal() << ":"
	// << am->getStart(opt.seqno) << "-"
	// << am->getStart(opt.seqno)+am->getLength(opt.seqno) << endl;
	int sense_seq = am->getOrientation(opt.seqno);
	int sense_comp = am->getOrientation(1-opt.seqno);
	// checkpoint;
	mymatches.push_back(match(am->getId(),axis->getOrdinal(),am->getParentId(),
				  am->getStart(opt.seqno),am->getLength(opt.seqno),
				  sense_seq,sense_comp));
      }
      // checkpoint;
    }

    if (opt.verbose) 
      timestamp("Looking through runs & matches...done.");
    */ 
  } else {

    if (opt.verbose) 
      timestamp("Reading alignment records...");

    bool delete_stream = false;
    istream * ifs(&cin);
    if (opt.atac_file != "-") {
      ifs = new ifstream(opt.atac_file.c_str());
      delete_stream = true;
    }
    while (true) {
      std::string id=""; 
      long signed int fe=-1;
      FILE_POSITION_TYPE s=0;
      FILE_POSITION_TYPE l=0;
      unsigned int ss=0;
      unsigned int sc=0;
      char buffer[1024];
      (*ifs) >> id >> fe >> s >> l >> ss >> sc;
      if (ifs->fail()) {
	ifs->clear();
      }
      ifs->ignore();
      if (id == "") {
	break;
      }
      if (fe >= 0) {
	if (id != "" && l != 0) {
	  mymatches.push_back(match(id,fe,"",s,l,ss,sc));
	}
      } else {
	Lazy_Header_SI const & hi(ff->get_header_data(s+1));
	s = ff->get_seq_pos(s+1)-1;
	fe = hi.index()-1;
	mymatches.push_back(match(id,fe,"",s,l,ss,sc));
      }
    }
    if (delete_stream) {
      delete ifs;
    }

    if (opt.verbose) 
      timestamp("Reading alignment records...done.");
  }
  
  if (opt.verbose) 
    cerr << " Found " << mymatches.size() << " intervals." << endl;
  if (!opt.extract) {
  if (opt.verbose) 
    timestamp("Sorting runs & matches...");
  std::sort(mymatches.begin(),mymatches.end());
  if (opt.verbose) 
    timestamp("Sorting runs & matches...done.");
  }

  if (opt.nooverlap) {
    if (opt.verbose) 
      timestamp("Ensuring no run/match overlap...");
    std::vector<match>::iterator mymatchit,mymatchit1;
    std::list<std::vector<match>::iterator> dlist;
    mymatchit = mymatches.begin();
    while (mymatchit!=mymatches.end()) {
      mymatchit1 = mymatchit;
      ++mymatchit1;
      while (mymatchit1 != mymatches.end() && 
	     mymatchit1->fasta_entry == mymatchit->fasta_entry &&
	     mymatchit1->start < mymatchit->start + mymatchit->length) {
	if (opt.verbose) 
	  timestamp("Overlap observed...");
	if (opt.verbose) 
	  cerr << mymatchit1->id << " to be suppressed.\n";
	// We have overlap...
	if (mymatchit1->start + mymatchit1->length >
	    mymatchit->start + mymatchit->length) {
	  mymatchit->length = mymatchit1->start + mymatchit1->length -
	    mymatchit->start;
	}
	dlist.push_back(mymatchit1);
	mymatchit->id += ":";
	mymatchit->id += mymatchit1->id;
	++mymatchit1;
      }
      mymatchit = mymatchit1;
    }
    std::list<std::vector<match>::iterator>::iterator dlit=dlist.begin();
    while (dlit!=dlist.end()) {
      mymatches.erase(*dlit);
    }
    std::sort(mymatches.begin(),mymatches.end());
    if (opt.verbose) 
      timestamp("Ensuring no run/match overlap...done");
  }

  const int bufsize=8192;
  char buffer[bufsize];
  
  if (opt.verbose) 
    timestamp("Pulling out sequence from fasta file...");

  FILE_POSITION_TYPE i,j;

  std::string runword;
  if (opt.type == "r") {
    runword = "run";
  } else if (opt.type == "u") {
    runword = "match";
  } else if (opt.type == "s") {
    runword = "signature";
  } else {
    runword = "alignment";
  }

  if (!opt.notin) {
    
    std::vector<match>::iterator it;
    FILE_POSITION_TYPE abspos;
    for (it=mymatches.begin();it!=mymatches.end();++it) {
      // checkpoint;
      // cerr << it->fasta_entry << " " << it->start;
      bool b=true;
      if (it->fasta_entry >= 0) {
	b = ff->fasta_pos(it->fasta_entry,it->start);
      } else {
	ff->pos(it->start);
      }
      if (b) {
	abspos = ff->pos();
	// cerr << " " << abspos << endl;
	if (ff->is_subseq(abspos,abspos+it->length)) {
          if (!opt.extract) {
	  Lazy_Header_SI const & hi(ff->get_header_data(abspos+1));
	  it->start = ff->get_seq_pos(abspos+1)-1;
	  *opt.out << ">" << hi.short_header()
		   << " /" << runword << "=" << "{" << it->id << "}";
	  // *opt.out << " /entry=" << it->fasta_entry;
	  *opt.out /* << it->sense_seq << ","
		   << it->sense_comp << ")"*/
		   << " /start=" 
		   << it->start
		   << " /end=" 
		   << it->start + it->length
		   << " /length=" 
		   << it->length;
	  long unsigned int ncount=0;
	  long unsigned int maxcontign=0,curcontign=0;
	  for (i=0;i<it->length;i++) {
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
	  if (ncount > 0) {
	    *opt.out << " /Ns=" << ncount 
		     << " /nonNs=" << it->length - ncount
		     << " /maxContigN=" << maxcontign << endl;
	  } else {
	    *opt.out << endl;
          }
          }
	  if (!opt.headersonly) {
	    ff->pos(abspos);
	    string thestr=ff->getstr((unsigned int)it->length);
	    if (it->sense_seq && it->sense_comp) {
	      thestr = reverse_comp(thestr);
	    }
	    for (i=0,j=0;i<it->length;i++) {
	      buffer[j++] = thestr[i];
	      if (i%60==59 && !opt.extract) {
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
              if (!opt.extract) {
	        *opt.out << endl;
              } else {
	        *opt.out << opt.eos_char;
              }
	    }
	  }
	} else {
	  if (opt.verbose) 
	    timestamp("Warning: Sequence is not in a single fasta entry.");
	}
      } else {
	if (opt.verbose) 
	  timestamp("Warning: Can't set Fasta file postion.");
      }
    }

  } else {

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
	if (/*current_start_position <= current_end_position && */ 
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
	endabspos = ff->pos()-1/*eos char*/;
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
  }
  opt.out->flush();
  // *opt.out << "\nDone" << endl;
  if (opt.verbose) 
    timestamp("Pulling out sequence from fasta file...done.");
}
