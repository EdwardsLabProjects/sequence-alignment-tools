

#include <unistd.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include "keyword_tree.h"
#include "keyword_tree.t"
#include "shift_and_inexact.h"
#include "rand_hash_table.h"
#include "shift_and.h"
#include "fasta_io.h"
#include "fasta_io.t"
#include "char_io.t"
#include "primer_alignment.h"
#include "util.h"
#include "types.h"
#include "select.h"
#include "select.t"

#include "hash.h"
#include "perfposht.h"
#include "bitmap.h"

#if defined(PRIMER3TM)
extern "C" {
#include "oligotm.h"
}
#endif

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

char release_tag[] = "$Name:  $";

#if defined(__alpha) 
// Hack to get around problems with mmap
const unsigned long __sbrk_override = 1;
#endif

void newfail_handler() {
  timestamp("Fatal Error: Attempt to allocate memory with new() failed.");
  exit(1);
  return;
}

struct Options {
  std::string database;
  std::string background;
  std::string output0;
  std::string output;
  int verbose;
  int mersize;
  bool indels;
  bool rc;
  bool ignore;
  int nmismatch;
  bool self;
  int exitthresh;
  bool inexonly;
  bool cannon;
  string qtemp;
  string ttemp;
  int dbchunksize;
  int bgchunksize;
  bool headerself;
  int threeprime;
  int fiveprime;
  float tmtarget;
  float tmdelta;
  bool posmatch;
  bool chkpt;
  bool ascout;
  Options(int argc, char *argv[]);
  ~Options();
  void usage(char *msg=NULL);
};
Options::Options(int argc, char *argv[]) {
  signed char c;
  optarg = NULL;
  verbose = 0;
  mersize = 0;
  indels = false;
  nmismatch = 0;
  self = false;
  rc = false;
  ignore = false;
  exitthresh = -1;
  inexonly = false;
  cannon = false;
  qtemp = "";
  ttemp = "";
  output = "";
  output0 = "";
  dbchunksize = 0;
  bgchunksize = 0;
  headerself = false;
  threeprime = 0;
  fiveprime = 0;
  tmtarget = 0;
  tmdelta = 3;
  posmatch = true;
  chkpt = false;
  ascout = false;
  while ((c = getopt(argc, argv, "i:b:SC:O:o:Ihrvcm:k:K:l:e:Xt:T:HP3:5:M:D:AZ")) != -1)
    switch (c) {
    case 'm':
      mersize = atoi(optarg);
      break;
    case 'k':
      nmismatch = atoi(optarg);
      indels = true;
      break;
    case 'K':
      nmismatch = atoi(optarg);
      indels = false;
      break;
    case 'C':
      char *p;
      if ((p=strchr(optarg,','))!=NULL) {
	*p = '\0';
	dbchunksize = atoi(optarg);
	bgchunksize = atoi(p+1);	
	*p = ',';
      } else {
	dbchunksize = atoi(optarg);
	bgchunksize = atoi(optarg);
      }
      break;
    case 'e':
      exitthresh = atoi(optarg);
      break;
    case '3':
      threeprime = atoi(optarg);
      break;
#ifdef PRIMER3TM
    case '5':
      fiveprime = atoi(optarg);
      break;
    case 'M':
      tmtarget = atof(optarg);
      break;
#endif
    case 'D':
      tmdelta = atof(optarg);
      break;
    case 'i':
      database = optarg;
      break;
    case 'b':
      background = optarg;
      break;
    case 'S':
      self = true;
      break;
    case 'I':
      ignore = true;
      break;
    case 'c':
      cannon = false;
      break;
    case 'r':
      rc = true;
      cannon = true;
      break;
    case 'o':
      output = optarg;
      break;
    case 'O':
      output0 = optarg;
      break;
    case 'X':
      inexonly = true;
      break;
    case 'H':
      headerself = true;
      break;
    case 'P':
      posmatch = false;
      break;
    case 'A':
      ascout = true;
      break;
    case 'Z':
      chkpt = true;
      break;
    case 'l':
      stderr = fopen(optarg,"a");
      setlinebuf(stderr);
      break;
    case 'v':
      verbose++;
      break; 
    case 't':
      qtemp = string(optarg);
      break;
    case 'T':
      ttemp = string(optarg);
      break;
    case 'h':
    default :
      usage();
    }
  if (database == "" || mersize == 0) usage();
  if (self) {
    background = database;
  } 
  if (background == "") usage();
  if (output == "-") {
    chkpt = false;
  }
  if (cannon && !rc) {
    usage("Cannonical mers are only used with reverse complement searching");
  }
  if (qtemp == "") {
    usage("Required option -t not specified.");
  }
  if (ttemp == "" && qtemp != "") {
    ttemp = qtemp;
  }
}


Options::~Options() {
  ;
}

void Options::usage(char *message) {
  if (message != NULL && strlen(message) > 0) {
    cerr << message << endl;
    cerr << endl;
  }
  cerr << "Usage: allvall [options] \n\n";
  cerr << "Options: \n";
  cerr << "  -m <int>               Mersize of mers.\n";
  cerr << "  -k <int>               Edit distance.\n";
  cerr << "  -K <int>               Hamming distance.\n";
  cerr << "  -5 <int>               Number of exact match bases at 5' end of mer.\n";
  cerr << "  -3 <int>               Number of exact match bases at 3' end of mer.\n";
#ifdef PRIMER3TM
  cerr << "  -M <float>             Melting temperature target. Default: No Tm constraint.\n";
  cerr << "  -D <float>             Melting temperature max delta. Default: at most 3 degrees.\n";
#endif
  cerr << "  -r                     Consider reverse complement matches too. Default: False.\n";
  cerr << "  -X                     Consider inexact matches only. Default: False.\n";
  cerr << "  -H                     Consider matches to different last header word only. Default: False.\n";
  cerr << "  -P                     Consider matches at same sequence position offset. Default: False.\n";
  cerr << "  -i <sequence-database> Input sequence database. Required.\n";
  cerr << "  -b <sequence-database> Background sequence database. -b or -S required.\n";
  cerr << "  -S                     Search self as background. -b or -S required.\n";
  cerr << "  -C <int>               Sequence database chunk size.\n";
  cerr << "  -t <seed-template>     Seed template for input sequence database. Required.\n";
  cerr << "  -T <seed-template>     Seed template for background sequence database. Default: Same as for -t.\n";
  cerr << "  -c                     Do not use cannonical mer for forward and reverse comp.\n";
  cerr << "  -o <output-file>       Output file name. The empty-string implies no matches will be saved, while - implies stdout.\n";
  cerr << "  -O <output-file>       Initialize match bitmap from output-file. Default: Same as for -o.\n";
  cerr << "  -Z                     Write periodic match-bitmap checkpoints.\n";
  cerr << "  -A                     Ascii match bitmap-format.\n";
  cerr << "  -I                     Ignore background sequence position in match bitmap checkpoint file.\n";
  cerr << "  -l <log-file>          Redirect stderr.\n";
  cerr << "  -e <int>               Exit status 2 if less than threshold.\n";
  cerr << "                         Default: Exit status 2 on error only.\n";
  cerr << "  -v                     Verbose (version & diagnostic) output.\n";
  cerr << "  -h                     Command line option help.\n";
  cerr << "\n";
  exit(1);
}

void
progress(CharacterProducer const & cp1, bitmap const & match, int interval, bool force=false) {
  static int nextprogress=0;
  if (nextprogress == 0 || force) {
    nextprogress = time(NULL);
  }
  if (time(NULL) >= nextprogress) {
    ostrstream ss;
    ss.setf(ios::fixed);
    ss << "Progress:";
    ss << setprecision(3) << setw(7) << cp1.progress()*100;
    ss << "%";
    ss << " Matched:"; 
    ss << setprecision(3) << setw(7) << 100.0*match.nset()/match.size();
    ss << "% ";
    ss << setw(10) << match.nunset();
    ss << " left.";
    ss << ends;
    std::string v(ss.str());
    timestamp(v.c_str());
    nextprogress += interval;
  }
}

void 
write_chkpnt(FILE_POSITION_TYPE pos, bitmap const & match, std::string & filename, int interval, bool ascii, bool force=false) {
  static int nextprogress=0;
  if (nextprogress == 0) {
    nextprogress = time(NULL)+interval;
  }
  if (time(NULL) >= nextprogress || force) {
    if (filename != "") {
      ostrstream ss;
      ss << "BEGIN" << endl;
      ss << 0 << " " << pos << endl;
      match.write(ss, ascii);
      ss << "END" << endl;
      
      string s(ss.str());
      
      // Atomic write of match state
      if (filename == "-") {
	std::cout.write(s.c_str(),s.length());
      } else {
        ofstream ofs(filename.c_str());
        ofs.write(s.c_str(),s.length());
        ofs.close();
      }
    }
    nextprogress += interval;
  }
}

void fpw(FILE* fp, ostrstream const & oss) {
  string s(oss.str());
  fwrite((void*)s.c_str(),sizeof(char),s.length(),fp);
}

int main(int argc,char *argv[]) {
  
  set_new_handler(newfail_handler);
#if defined(PROFILE)
  set_profile_signal_handler();
#endif
  assert(sizeof(FILE_POSITION_TYPE)==8);
  assert(sizeof(uint64)==8);
  assert(sizeof(uint32)==4);
  assert(sizeof(int64)==8);
  assert(sizeof(int32)==4);

  Options opt(argc,argv);

  int dbblock=0;
  int bgblock=0;
  int p = opt.database.rfind('.');
  if (p != std::string::npos) {
    dbblock = atoi(opt.database.substr(p+1).c_str());
  }
  p = opt.background.rfind('.');
  if (p != std::string::npos) {
    bgblock = atoi(opt.background.substr(p+1).c_str());
  }

  FILE_POSITION_TYPE dboffset = ((FILE_POSITION_TYPE)opt.dbchunksize)*dbblock;
  FILE_POSITION_TYPE bgoffset = ((FILE_POSITION_TYPE)opt.bgchunksize)*bgblock;
  timestamplu("database offset: ",dboffset);
  timestamplu("background offset: ",bgoffset);
  FastaFile<Lazy_Header_SI>* ff, *ff1;
  fasta_file_seq_params ffp;
  ffp.upper_case = true;
  ffp.eos_char = '$';
  ffp.check_params = false;
  ffp.translate = false;
  ffp.eos_start = true;
  ffp.offset = dboffset;
  ff  = pick_fasta_file<Lazy_Header_SI>(opt.database,0,false,true,ffp,opt.verbose);
  ffp.offset = bgoffset;
  if (opt.bgchunksize == 0 && opt.dbchunksize > 0) {
    // Guess when to memmap vs load entire file into memory...
    ff1 = pick_fasta_file<Lazy_Header_SI>(opt.background,0,true,true,ffp,opt.verbose);
  } else {
    ff1 = pick_fasta_file<Lazy_Header_SI>(opt.background,0,false,true,ffp,opt.verbose);
  }

  FastaFile<Lazy_Header_SI> & cp  = *ff;
  FastaFile<Lazy_Header_SI> & cp1 = *ff1;

  editdist_alignment pa(0,0,'$');
  pa.kmax(opt.nmismatch);
  pa.indels(opt.indels);
  pa.eos('$');
  pa.tn(false);
  pa.wc(false);
  pa.dna_mut(false);
  if (opt.verbose >= 3) {
    pa.yesno(false);
  } else {
    pa.yesno(true);
  }

  int maxdelta = opt.indels?opt.nmismatch:0;
  int mindist = opt.mersize;
  int interval = 5;
  int chkpnt_write_interval = 60;

  pa.maxpatlen(opt.mersize+maxdelta);
  
  bitmap match(cp.length()+1);

  // cerr << cp.length() << " " << match.size() << endl;

  ifstream *match_in=0;
  if (opt.output0 != "") {
    match_in = new ifstream(opt.output0.c_str());
  } else if (opt.output != "" && opt.output != "-") {
    match_in = new ifstream(opt.output.c_str());
  }
  bool newoutfile = true;
  FILE_POSITION_TYPE posin = 0;
  uint32 spanin = 0;
  if (match_in) {
    if (*match_in) {
      newoutfile = false;
      string l;
      assert(getline(*match_in,l) && l == "BEGIN");
      (*match_in) >> spanin >> posin;
      assert(getline(*match_in,l));
      match.read(*match_in);
      assert(getline(*match_in,l) && l == "END");
    }
    match_in->close();
    delete match_in;
  }

  if (opt.ignore) {
    posin = 0;
    spanin = 0;
  }

  

  // cerr << cp.length() << " " << match.size() << endl;

  for (int i=0;i<opt.mersize;i++) {
    // ostrstream oss;
    // oss << " " << i << " " << -1 << endl;
    // fpw(stderr,oss);
    match.set(i);
  }
  if (opt.dbchunksize > 0) {
    for (int i=opt.dbchunksize+opt.mersize;i<match.size();i++) {
      // ostrstream oss;
      // oss << " " << i << " " << -1 << endl;
      // fpw(stderr,oss);
      match.set(i);
    }
  }

  hash *h = hashselect(cp,4,opt.qtemp.c_str());
  hash *h1ptr = hashselect(cp1,4,opt.ttemp.c_str());

  if (!h->symmetric() || !h1ptr->symmetric()) {
    timestamp("Asymmetric hash => cannonical mers cannot be used");
    opt.cannon = false;
  }
  perfposht mers;
  if (newoutfile) {
    mers.init(*h,opt.rc,opt.cannon,dboffset);
  } else {
    mers.init(*h,opt.rc,opt.cannon,dboffset,opt.mersize,&match);
  }
  delete h;

  hash & h1 = *h1ptr;
  if (posin > 0) {
    h1.reset(posin+1);
  } else {
    h1.reset();
  }

  if (opt.verbose) {
    progress(cp1,match,interval);
  }

  unsigned int lastp1=0;
  vector<bool> pairseenvec(2*cp.length()+1,false);
  list<unsigned int> pairseenlist;

  hash_t v1;
  bool v1rc=false;
  while (h1.next()) {
    if (!opt.cannon) {
      v1 = h1.value();
    } else {
      v1 = h1.cvalue();
      v1rc = (v1 != h1.value());
    }
    if (!mers.lookup(v1)) {
      continue;
    }
    if (h1.ns() > 0) {
      // checkpoint;
      // cerr << h1.ns() << endl;
      continue;
    }

    FILE_POSITION_TYPE p1 = h1.pos();
    FILE_POSITION_TYPE oldpos = cp1.pos();

    // checkpoint;
    // cerr << p1 << endl;

    /* if (v1 == 43080) {
       checkpoint;
       cerr << p1 << '\t'
       << v1 << '\t'
       << h1.str(v1) << endl;
       }
    */

    if (p1 != lastp1) {
      // clear pairseen...
      list<unsigned int>::iterator it=pairseenlist.begin();
      while (it != pairseenlist.end()) {
	pairseenvec[*it] = false;
	++it;
      }
      pairseenlist.clear();
      lastp1 = p1;
    }

    bool anypos=false;
    while (mers.next()) {
      FILE_POSITION_TYPE p;
      bool rc;
      if (!v1rc) {
	if (mers.value() > 0) {
	  p = mers.value();
	  rc = false;
	} else {
	  p = -mers.value();
	  rc = true;
	} 
      } else {
	if (mers.value() > 0) {
	  p = mers.value();
	  rc = true;
	} else {
	  p = -mers.value();
	  rc = false;
	} 
      }

      anypos = true;

      ostrstream oss;
      if (opt.verbose >= 4) {
	oss << "Hash hit at q= " << p << " b= " << p1 << " (" << h1.span() << ")";
      }

      unsigned int psvi = 2*p+1*rc;
      if (pairseenvec[psvi]) { // Need more than position, need revcomp, monotonic span
	if (opt.verbose >= 4) {
	  oss << " suppressed" << endl;
	  fpw(stderr,oss);
	}
	continue;
      } 

      if (opt.verbose >= 4) {
	oss << endl;
	fpw(stderr,oss);
      }

      pairseenvec[psvi] = true;
      pairseenlist.push_back(psvi);

      p += dboffset;

      // checkpoint;
      // cerr << "Hash hit at: " << p << endl;

      // Check all of the query mers consistent with p
      // checkpoint;
      FILE_POSITION_TYPE startp = max(p-opt.mersize,dboffset);
      FILE_POSITION_TYPE endp = min(p+(opt.mersize-h1.span()),dboffset+cp.length());

      FILE_POSITION_TYPE offset = 0;
      if (startp == dboffset) {
	offset = dboffset+opt.mersize-p;
      }

      // if (v1 == 43080) {
      // checkpoint;
      // fprintf(stderr,"p %d startp %d endp %d h1.span() %d \n",p,startp,endp,h1.span());
      // fprintf(stderr,"p1 %d\n",p1);
      // }

      string qregion;
      FILE_POSITION_TYPE matchpos;
      FILE_POSITION_TYPE npos = endp-startp-opt.mersize;
      FILE_POSITION_TYPE npos1 = npos+(h1.span()-h1.minspan());
      bool anychecked=false;
      for (uint32 l=0;l<=npos;l++) {
	// checkpoint;
	// timestampi("l: ",l);
	if (!rc) {
	  matchpos = p+l+offset;
	} else {
	  matchpos = p+npos-l;
	}
	// timestampli("matchpos: ",matchpos);
	// timestampli("matchpos: ",matchpos-dboffset);
	if (!match[matchpos-dboffset]) {
	  if (!anychecked) {
	    anychecked = true;
	    cp.pos(startp);
	    qregion = cp.getstr((unsigned int)(endp-startp));
	    if (rc) {
	      qregion = reverse_comp(qregion);
	    }
	  }
	  string qstr = qregion.substr(l,opt.mersize);
	  // timestamplu("matchpos: ",matchpos);
	  // if (v1 == 43080) {
	  // fprintf(stderr,"qstr %s\n",qstr.c_str());
	  // } 
	  if (qstr.find_first_not_of("ACGT") != string::npos) {
	    // match[matchpos] = cp.length()+mindist;
	    // *opt.out << matchpos << " " << cp.length()+mindist << endl;
	    if (opt.verbose >= 2) {
	      if (opt.verbose >= 3) {
		// checkpoint;
		const Lazy_Header_SI & h  = cp.get_header_data(matchpos);
		ostrstream oss;
		oss << " " 
		    << opt.database << " "
		    << matchpos << " "
		    << h.short_header() << " "
		    << cp.get_seq_pos(matchpos) << " " 
		    << qstr.c_str() << " "
		    << "Non-ACGT" << endl;
		  fpw(stderr,oss);
	      } else {
		// cerr << matchpos << " " << cp.length() << " " << mindist << " " << cp.length()+mindist << endl;
		ostrstream oss;
		oss << " " << matchpos << " " << 0 << endl;
		fpw(stderr,oss);
		// fprintf(stderr," %lu %ld\n",matchpos,cp.length()+mindist);
	      }
	    }
	    match.set(matchpos-dboffset);
	    // total++;
	    continue;
	  }
#ifdef PRIMER3TM
	  if (opt.tmtarget > 0) {
	    /* Defaults from Primer3 source */
	    float tm = oligotm(qstr.c_str(),50.0,50.0,0.0,0.0,TM_METHOD_SANTALUCIA,SALT_CORRECTION_SANTALUCIA);
	    if (fabs(tm-opt.tmtarget) > opt.tmdelta) {
	      if (opt.verbose >= 2) {
		if (opt.verbose >= 3) {
		  // checkpoint;
		  const Lazy_Header_SI & h  = cp.get_header_data(matchpos);
		  ostrstream oss;
		  oss << " " 
		      << opt.database << " "
		      << matchpos << " "
		      << h.short_header() << " "
		      << cp.get_seq_pos(matchpos) << " " 
		      << qstr.c_str() << " "
		      << "Bad-Tm" << " "
		      << tm << endl;
		  fpw(stderr,oss);
		} else {
		  ostrstream oss;
		  oss << " " << matchpos << " " << 0 << endl;
		  fpw(stderr,oss);
		}
	      }
	      match.set(matchpos-dboffset);
	      continue;
	    }
	  }
#endif
	  // fprintf(stderr,"%lu %lu\n",matchpos,p1+l);
	  // fprintf(stderr,"%s %s %lu %lu %lu \n",((opt.self)?"TRUE":"FALSE"),
	  // ((rc)?"TRUE":"FALSE"),matchpos - mindist,p1 + l,matchpos + mindist);
	  // # Fix this!
	  if (opt.posmatch && !rc && 
	      p1 + l >= matchpos - mindist && 
	      p1 + l <= matchpos + mindist) {
	    continue;
	  }
	  if (opt.headerself) {
	    const Lazy_Header_SI & h1 = cp1.get_header_data(p1 + l);
	    const Lazy_Header_SI & h  = cp.get_header_data(matchpos);
	    // cerr << h1.index() << " " << h.index() << endl;
	    if (h1.index() == h.index()) {
	      // cerr << h1.index() << " " << h.index() << endl;
	      continue;
	    }
	    std::string header1 = h1.header();
	    std::string header = h.header();
	    int spacepos1 = header1.rfind(' ');
	    int spacepos = header.rfind(' ');
	    assert(spacepos1 != std::string::npos);
	    // cerr << header << " " << spacepos << endl;
	    assert(spacepos != std::string::npos);
	    if (header1.substr(spacepos1+1) == header.substr(spacepos+1)) {
	      // cerr << header1.substr(spacepos1+1) << " " << header.substr(spacepos+1) << endl;
	      continue;
	    }
	  }
	  pa.reset();
	  if (!rc) {
	    pa.exact_start_bases(opt.fiveprime);
	    pa.exact_end_bases(opt.threeprime);
	  } else {
	    pa.exact_start_bases(opt.threeprime);
	    pa.exact_end_bases(opt.fiveprime);
	  }

	  FILE_POSITION_TYPE lb=p1+l-maxdelta;
	  FILE_POSITION_TYPE ub=p1+l+maxdelta;
	  // checkpoint;
	  // cerr << lb << " " << ub << endl;
	  if (lb < bgoffset+opt.mersize-maxdelta) {
	    lb = bgoffset+opt.mersize-maxdelta;
	  }
	  // checkpoint;
	  // cerr << lb << " " << ub << endl;
	  if (ub >= bgoffset+cp1.length()) {
	    ub = bgoffset+cp1.length()-1;
	  }
	  // if (v1 != 43080) {
	  // checkpoint;
	  // cerr << lb << " " << ub << endl;
	  // }
	  /*
	    if (opt.self && !rc) {
	    if (matchpos <= lb + mindist && matchpos + mindist >= lb) {
	    lb = matchpos+mindist+1; 
	    } 
	    if (matchpos + mindist >= ub && matchpos <= ub + mindist) {
	    ub = matchpos-mindist-1;
	    }
	    }
	    if (ub < 0 || lb > ub) {
	    continue;
	    }
	  */
	  pa.poslb(lb);
	  pa.posub(ub);
	  if (pa.align(cp1,qstr) && (!opt.inexonly || pa.value() > 0)) {
	    // cerr << matchpos << endl;
	    match.set(matchpos-dboffset);
	    if (opt.verbose >= 2) {
	      if (!rc) {
		if (opt.verbose >= 3) {
		  // checkpoint;
		  const Lazy_Header_SI & h1 = cp1.get_header_data(pa.end());
		  const Lazy_Header_SI & h  = cp.get_header_data(matchpos);
		  ostrstream oss;
		  oss << " " 
		      << opt.database << " " 
		      << matchpos << " " 
		      << h.short_header() << " "
		      << cp.get_seq_pos(matchpos) << " "
		      << qstr << " " 
		      << opt.background << " "
		      << pa.end() << " "
		      << h1.short_header() << " "
		      << cp1.get_seq_pos(pa.end()) << " "
		      << pa.matching_text() << " "
		      << pa.alignment_text() << " " 
		      << pa.alignment_string() << " "
		      << pa.alignment_pattern(qstr) << " "
		      << "F" << endl;
		  fpw(stderr,oss);
		  // checkpoint;
		} else {
		  ostrstream oss;
		  oss << " " << matchpos << " " << pa.end() << endl;
		  fpw(stderr,oss);
		  // fprintf(stderr," %lu %ld\n",matchpos,pa.end());free(apstr);
		}
	      } else {
		if (opt.verbose >= 3) {
		  // checkpoint;
		  const Lazy_Header_SI & h1 = cp1.get_header_data(pa.end());
		  const Lazy_Header_SI & h  = cp.get_header_data(matchpos);
		  ostrstream oss;
		  oss << " " 
		       << opt.database << " " 
		       << matchpos << " " 
		       << h.short_header() << " "
		       << cp.get_seq_pos(matchpos) << " "
		       << reverse_comp(qstr) << " " 
		       << opt.background << " "
		       << pa.end() << " "
		       << h1.short_header() << " "
		       << cp1.get_seq_pos(pa.end()) << " "
		       << reverse_comp(pa.matching_text()) << " "
		       << reverse_comp(pa.alignment_text()) << " " 
		       << reverse(pa.alignment_string()) << " "
		       << reverse_comp(pa.alignment_pattern(qstr)) << " "
		       << "R" << endl;
		  fpw(stderr,oss);
		  // checkpoint;
		} else {
		  ostrstream oss;
		  oss << " " << matchpos << " " << -pa.end() << endl;
		  fpw(stderr,oss);
		  // fprintf(stderr," %lu %ld\n",matchpos,-pa.end());
		}
	      }
	    }

	    /* if (v1 == 43080) {
	       checkpoint;
	       cerr << pa.value() << endl;
	       } */
	    /*total ++;*/
	    /* 
	       << " " << h1.str(h1.value()) << " " << h1.str(h1.rcvalue()) << " " << ((rc)?" RC ":"    ") << " " << p << " " << p1 << endl;
	       *opt.out << "  " << pa.alignment_text() << " " << pa.editdist() << endl;
	       *opt.out << "  " << pa.alignment_string() << endl;
	       *opt.out << "  " << pa.alignment_pattern(qstr) 
	       << endl;
	    */
	  }
	}
      }
      if (!anychecked) {
	/* if (v1 == 43080) {
	   checkpoint;
	   cerr << npos << " " << npos1 << endl;
	   } */
	for (int l=npos+1;l<=npos1;l++) {
	  matchpos = p+l;
	  /* if (v1 == 43080) {
	     checkpoint;
	     cerr << "Checking " << matchpos << endl;
	     }*/
	  if (!match[matchpos-dboffset]) {
	    anychecked = true;
	  }
	}
	if (!anychecked) {
	  /* if (v1 == 43080) {
	     checkpoint;
	     std::cerr << "Set mer ending at " << p << " (" << (rc?"R":"F") << ") to invalid" << std::endl;
	     } */
	  mers.set_invalid();
	}
      }
    }
    if (!anypos) {
      mers.set_invalid(v1);
    }
    cp1.pos(oldpos);
    if (opt.verbose) {
      progress(cp1,match,interval);
    }
    if (opt.chkpt) {
      write_chkpnt(p1,match,opt.output,chkpnt_write_interval,opt.ascout);
    }
    if (match.nunset() == 0) {
      break;
    }
  }

  if (opt.verbose) {
    progress(cp1,match,interval,true);
  }

  delete h1ptr;

  write_chkpnt(0,match,opt.output,chkpnt_write_interval,opt.ascout,true);

  if (opt.exitthresh > 0 && match.nunset() < opt.exitthresh) {
    return 2;
  }
  return 0;
}

