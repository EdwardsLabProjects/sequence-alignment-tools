
#include <unistd.h>
#include <assert.h>
#include <ctype.h>
#include <iostream>
#include <vector>
#include "sortedvector.t"
#include "fasta_io.h"
#include "char_io.h"
#include "util.h"

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

char release_tag[] = "$Name:  $";

#if defined(__alpha) 
// Hack to get around problems with mmap
const unsigned long __sbrk_override = 1;
#endif

#define _NO_LARGEFILE_STREAMS

#include "rl_suffix_tree.h"

void make_suftree(std::string base, std::string ext1, std::string ext2) {
  string infile = base + ext1;
  string outfile = base + ext2;
  FILE_POSITION_TYPE size = filesize(infile);
  char *buffer=new char[size+1];
  FILE* in=fopen(infile.c_str(),"r");
  fread(buffer,sizeof(char),size,in);
  buffer[size]='\0';
  fclose(in);

  suffix_tree<basic_node,basic_suffix> T(buffer,size,buffer[0]);
  T.build();
  T.write(outfile.c_str());
}


void usage(const char *message=NULL) {
  if (message != NULL && strlen(message) > 0) {
    cerr << message << endl;
    cerr << endl;
  }
  cerr << "Usage: compress_seq [options] \n\n";
  cerr << "Options: \n";
  cerr << "  -i <sequence-database> Input sequence database.\n";
  cerr << "                         Required.\n";
  cerr << "  -e [true|false]        Insert end-of-sequence marker.\n";
  cerr << "                         Default: true. Optional.\n";
  cerr << "  -S [true|false]        Insert end-of-sequence marker before\n";
  cerr << "                         initial sequence. Default: true. Optional.\n";
  cerr << "  -E <char>              Use specified character as single character\n";
  cerr << "                         end-of-sequence marker.\n";
  cerr << "                         Specified as decimal, octal, or hex. integer.\n";
  cerr << "                         Default: \"\\n\". Optional.\n";
  cerr << "  -3 <char>              Use specified character as three character\n";
  cerr << "                         end-of-sequence marker.\n";
  cerr << "                         Specified as decimal, octal, or hex. integer.\n";
  cerr << "                         Default: \"\\n\". Optional.\n";
  cerr << "  -u [true|false]        Uppercase sequence characters.\n";
  cerr << "                         Default: true. Optional.\n";
  cerr << "  -n [true|false]        Normalize sequence information.\n";
  cerr << "                         Default: false. Optional.\n";
  cerr << "  -D [true|false]        Optimize normalized sequence for DNA.\n";
  cerr << "                         Default: true. Optional.\n";
  cerr << "  -R [true|false]        Ensure reverse complement characters are in\n";
  cerr << "                         normalized character table. Default: false.\n";
  cerr << "                         Optional.\n";
  cerr << "  -z [true|false]        Compress normalized sequence information.\n";
  cerr << "                         Default: false. Optional.\n";
  cerr << "  -t [true|false]        Compute suffix tree for sequence.\n";
  cerr << "                         Default: false. Optional.\n";
  cerr << "  -I [true|false]        Write sequence and header index in binary format.\n";
  cerr << "                         Default: true. Optional.\n";
  cerr << "  -T [true|false]        Output character table only. Default: false. Optional.\n";
  cerr << "  -c [true|false]        Non-zero exit status indicates that output files will\n";
  cerr << "                         be rebuilt. Default: false. Optional.\n";
  cerr << "  -F [true|false]        Force all output files to be re-built.\n";
  cerr << "                         Default: false. Optional.\n";
  cerr << "  -C [true|false]        Cleanup (delete) unnecessary files.\n";
  cerr << "                         Default: true. Optional.\n";
  cerr << "  -B                     Use buffered I/O rather than mmap.\n";
  cerr << "  -v                     Print version information.\n";
  cerr << "  -h                     Command line option help.\n";
  cerr << "\n";
  exit(1);
}

struct Options {
  bool quiet;
  std::string database;
  bool eos;
  bool init_eos;
  bool compress;
  bool uc;
  bool normalize;
  int eos_len;
  char eos_char;
  bool dnaopt;
  bool force;
  bool cleanup;
  bool bufferedio;
  bool binindex;
  bool verbose;
  bool tableonly;
  bool addrc;
  bool guard;
  bool checkonly;
  bool suftree;
  bool ingz;
};

void options(int argc, char *argv[], Options & opt) {
  signed char c;
  optarg = NULL;
  opt.eos = true;
  opt.compress = false;
  opt.uc = true;
  opt.normalize = false;
  opt.eos_char = '\n';
  opt.eos_len = 1;
  opt.dnaopt = true;
  opt.force = false;
  opt.cleanup = true;
  opt.bufferedio = false;
  opt.binindex = true;
  opt.init_eos = true;
  opt.verbose = false;
  opt.tableonly = false;
  opt.addrc = false;
  opt.guard = false;
  opt.checkonly = false;
  opt.suftree = false;
  opt.ingz = false;
  while ((c = getopt(argc, argv, "i:e:S:z:u:D:E:3:n:F:C:I:T:BR:hvG:c:t:")) != -1)
  switch (c) {
  case 'i':
    opt.database = std::string(optarg);
    break;
  case 'e':
    if (is_true(optarg)) {
      opt.eos = true;
    } else if (is_false(optarg)) {
      opt.eos = false;
    } else {
      usage("Invalid value for -e option.");
    }
    break;
  case 'S':
    if (is_true(optarg)) {
      opt.init_eos = true;
      opt.eos = true;
    } else if (is_false(optarg)) {
      opt.init_eos = false;
    } else {
      usage("Invalid value for -S option.");
    }
    break;
  case '3': {
    int ch;
    sscanf(optarg,"%i",&ch);
    opt.eos_char = ch;
    opt.eos_len = 3;
    }
    break;
  case 'E': {
    int ch;
    sscanf(optarg,"%i",&ch);
    opt.eos_char = ch;
    opt.eos_len = 1;
    }
    break;
  case 'u':
    if (is_true(optarg)) {
      opt.uc = true;
    } else if (is_false(optarg)) {
      opt.uc = false;
    } else {
      usage("Invalid value for -u option.");
    }
    break;
  case 'D':
    if (is_true(optarg)) {
      opt.dnaopt = true;
    } else if (is_false(optarg)) {
      opt.dnaopt = false;
    } else {
      usage("Invalid value for -D option.");
    }
    break;
  case 'R':
    if (is_true(optarg)) {
      opt.addrc = true;
    } else if (is_false(optarg)) {
      opt.addrc = false;
    } else {
      usage("Invalid value for -R option.");
    }
    break;
  case 'I':
    if (is_true(optarg)) {
      opt.binindex = true;
    } else if (is_false(optarg)) {
      opt.binindex = false;
    } else {
      usage("Invalid value for -I option.");
    }
    break;
  case 'T':
    if (is_true(optarg)) {
      opt.tableonly = true;
    } else if (is_false(optarg)) {
      opt.tableonly = false;
    } else {
      usage("Invalid value for -T option.");
    }
    break;
  case 'n':
    if (is_true(optarg)) {
      opt.normalize = true;
    } else if (is_false(optarg)) {
      opt.normalize = false;
    } else {
      usage("Invalid value for -n option.");
    }
    break;
  case 'z':
    if (is_true(optarg)) {
      opt.compress = true;
    } else if (is_false(optarg)) {
      opt.compress = false;
    } else {
      usage("Invalid value for -z option.");
    }
    break;
  case 'F':
    if (is_true(optarg)) {
      opt.force = true;
    } else if (is_false(optarg)) {
      opt.force = false;
    } else {
      usage("Invalid value for -F option.");
    }
    break;
  case 'C':
    if (is_true(optarg)) {
      opt.cleanup = true;
    } else if (is_false(optarg)) {
      opt.cleanup = false;
    } else {
      usage("Invalid value for -C option.");
    }
    break;
  case 't':
    if (is_true(optarg)) {
      opt.suftree = true;
    } else if (is_false(optarg)) {
      opt.suftree = false;
    } else {
      usage("Invalid value for -t option.");
    }
    break;
  case 'G':
    if (is_true(optarg)) {
      opt.guard = true;
    } else if (is_false(optarg)) {
      opt.guard = false;
    } else {
      usage("Invalid value for -G option.");
    }
    break;
  case 'c':
    if (is_true(optarg)) {
      opt.checkonly = true;
    } else if (is_false(optarg)) {
      opt.checkonly = false;
    } else {
      usage("Invalid value for -c option.");
    }
    break;
  case 'B':
    opt.bufferedio = true;
    break;
  case 'v':
    opt.verbose = true;
    break; 
  case 'h':
  default :
    usage();
  }
  if (opt.database == ""&&!opt.verbose) usage("Arguement -i missing from commandline.");
  if (opt.suftree && !opt.init_eos) usage("Option intial EOS (-S true) required to build suffix tree (-t true).");
#if !defined(NO_ZLIB)
  int p = opt.database.rfind('.');
  if (p > 0 && opt.database.substr(p) == string(".gz")) {
    opt.database = opt.database.substr(0,p);
    opt.ingz = true;
  } 
#endif
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

  if (opt.database == "") {
    exit(1);
  }
  time_t fasta_time = modtime(opt.database);
  time_t grd_time = modtime(opt.database+".grd");
  time_t seq_time = modtime(opt.database+".seq");
  time_t hdr_time = modtime(opt.database+".hdr");
  time_t idx_time = modtime(opt.database+".idx");
  time_t idb_time = modtime(opt.database+".idb");
  time_t tbl_time = modtime(opt.database+".tbl");
  time_t sqn_time = modtime(opt.database+".sqn");
  time_t tbz_time = modtime(opt.database+".tbz");
  time_t sqz_time = modtime(opt.database+".sqz");

  time_t seqst_time = modtime(opt.database+".seq.st");
  time_t sqnst_time = modtime(opt.database+".sqn.st");
  time_t sqzst_time = modtime(opt.database+".sqz.st");

  if (opt.guard) {
    ofstream grd;
#if   ! defined(__CYGWIN__)
    grd.open((opt.database+".grd").c_str());
#else 
    grd.open((opt.database+".grd").c_str(),ios::binary);
#endif
    grd << "Guard file for compress_seq" << endl;
    grd.close();
  }

  if (idx_time < idb_time) idx_time = idb_time;

  bool didinitscan=false;
  std::vector<bool> obs(256,false);
  if (opt.eos) {
    obs[(int)opt.eos_char] = true;
  }

  FILE_POSITION_TYPE outputi=0;
#define OUTPUTBUFSIZE 8192
  char outputbuffer[OUTPUTBUFSIZE];
  
  if (opt.force || 
      (opt.guard && grd_time > 0) ||
      ((!opt.compress && !opt.normalize &&
	fasta_time > seq_time) || 
       fasta_time > hdr_time ||
       fasta_time > idx_time) ||
      (opt.tableonly && fasta_time > tbl_time) ||
      (opt.compress && 
       (fasta_time > idx_time || fasta_time > tbz_time || fasta_time > sqz_time)) ||
      (opt.normalize && 
       (fasta_time > idx_time || fasta_time > tbl_time || fasta_time > sqn_time))
      ) {
    if (opt.checkonly) {
      return 1;
    }
    didinitscan=true;
    CharacterProducer* db;
    if (!opt.ingz) {
      if (!opt.bufferedio) {
	db = new MapFileChars(opt.database.c_str(),'\n');
      } else {
	db = new BufferedFileChars(opt.database.c_str(),'\n');
      }
#if !defined(NO_ZLIB)
    } else {
      db = new GZChars((opt.database+".gz").c_str(),'\n');
#endif
    }
#if defined(_NO_LARGEFILE_STREAMS)
    FILE *seq;
#else
    ofstream seq;
#endif
    if (!opt.tableonly) {
#if defined(_NO_LARGEFILE_STREAMS)
      seq = fopen((opt.database+".seq").c_str(),"wb");
#else
#if ! defined(__CYGWIN__) 
      seq.open((opt.database+".seq").c_str());
#else
      seq.open((opt.database+".seq").c_str(),ios::binary);
#endif
#endif
    }

#if defined(_NO_LARGEFILE_STREAMS)
    FILE *hdr;
#else
    ofstream hdr;
#endif
    if (!opt.tableonly) {
#if defined(_NO_LARGEFILE_STREAMS)
      hdr = fopen((opt.database+".hdr").c_str(),"wb");
#else
#if ! defined(__CYGWIN__) 
      hdr.open((opt.database+".hdr").c_str());
#else
      hdr.open((opt.database+".hdr").c_str(),ios::binary);
#endif 
#endif
    }

#if defined(_NO_LARGEFILE_STREAMS)
    FILE *idx;
#else
    ofstream idx;
#endif
    if (!opt.binindex && !opt.tableonly) {
#if defined(_NO_LARGEFILE_STREAMS)
    idx = fopen((opt.database+".idx").c_str(),"wb");
#else
#if ! defined(__CYGWIN__)
      idx.open((opt.database+".idx").c_str());
#else 
      idx.open((opt.database+".idx").c_str(),ios::binary);
#endif
#endif
    } 
    
    long unsigned int count=0;
    FILE_POSITION_TYPE headerpos=0;
    FILE_POSITION_TYPE seqpos=0;
    sortedvector<FILE_POSITION_TYPE,FILE_POSITION_TYPE> svi;
    if (opt.init_eos) {
      for (int ii=0;ii<opt.eos_len;ii++) {
        outputbuffer[outputi++] = opt.eos_char;
        seqpos++;
      }
    }
    if (!opt.tableonly) {
      if (!opt.binindex) {
#if defined(_NO_LARGEFILE_STREAMS)
	fprintf(idx,"%lu %llu %llu %llu\n",count,headerpos,seqpos,(long long unsigned int)0);
#else
	idx << count << " " << headerpos << " " << seqpos << " " << 0 << '\n';
#endif
      } else {
	svi.push_back(seqpos,headerpos);
      }
    }
    
    bool inseq=false;
    bool inheader=false;
    bool startofline=true;
    
    char ch;
    
    while (!db->eof()) {
      ch = db->getch();
      if (startofline && ch == '>') {
	if (inseq) {
	  if (opt.eos) {
            for (int ii=0;ii<opt.eos_len;ii++) {
              outputbuffer[outputi++] = opt.eos_char;
              seqpos++;
            }
	  }
	  if (outputi >= (OUTPUTBUFSIZE-4)) {
	    if (!opt.tableonly) {
#if defined(_NO_LARGEFILE_STREAMS)
	      fwrite(outputbuffer,sizeof(char),outputi,seq);
#else
	      seq.write(outputbuffer,outputi);
#endif
	    }
	    outputi=0;
	  }
	  if (!opt.tableonly) {
	    if (!opt.binindex) {
#if defined(_NO_LARGEFILE_STREAMS)
	      fprintf(idx,"%llu %llu\n",seqpos,(db->pos())-1);
#else
	      idx << seqpos << " " << (db->pos())-1 << '\n';
#endif
	    } else {
	      svi.push_back(seqpos,headerpos);
	    }
	  }
	}
	inheader = true;
	inseq = false;
	startofline = false;
	continue;
      } else if (inheader) {
	if (ch == '\n' || ch == '\r') {
	  if (ch == '\r') {
	    ch = db->getch();
	    assert(ch == '\n');
	  } 
	  if (!opt.tableonly) {
#if defined(_NO_LARGEFILE_STREAMS)
	    fputc(ch,hdr);
#else
	    hdr << ch;
#endif
	  }
	  headerpos++;
	  inheader=false;
	  inseq=true;
	  startofline=true;
	  count++;
	  if (!opt.binindex && !opt.tableonly) {
#if defined(_NO_LARGEFILE_STREAMS)
	    fprintf(idx,"%lu %llu ",count,headerpos);
#else
	    idx << count << " " << headerpos << " ";
#endif
	  }
	  continue;
	} else {
	  if (!opt.tableonly) {
#if defined(_NO_LARGEFILE_STREAMS)
	    fputc(ch,hdr);
#else
	    hdr << ch;
#endif
	  }
	  headerpos++;
	}
	if (startofline) startofline=false;
	continue;
      } else if (inseq) {
	if (ch == '\n' || ch == '\r') {
	  if (ch == '\r') {
	    ch = db->getch();
	    assert(ch == '\n');
	  }
	  startofline = true;
	  continue;
	} else if ((int)ch < 33 || (int)ch > 126) {
	  if (startofline) startofline=false;
	  continue;
	} else {
	  if (opt.uc) ch = toupper(ch);
	  outputbuffer[outputi++] = ch;
	  if (outputi >= (OUTPUTBUFSIZE-4)) {
	    if (!opt.tableonly) {
#if defined(_NO_LARGEFILE_STREAMS)
	      fwrite(outputbuffer,sizeof(char),outputi,seq);
#else
	      seq.write(outputbuffer,outputi);
#endif
	    }
	    outputi=0;
	  }
	  seqpos++;
	  if (opt.normalize || opt.compress || opt.tableonly) {
	    obs[(int)ch] = true;
	    if (opt.addrc) {
	      obs[(int)iupac_revcomp(ch)] = true;
	    }
	  }
	  if (startofline) startofline=false;
	  continue;
	}
      }
    }
    if (inheader && !opt.tableonly) {
#if defined(_NO_LARGEFILE_STREAMS)
      fputc('\n',hdr);
#else
      hdr << '\n';
#endif
      headerpos++;
      count++;
      if (!opt.binindex) {
#if defined(_NO_LARGEFILE_STREAMS)
	fprintf(idx,"%lu %llu %llu %llu\n",count,headerpos,seqpos,(db->pos())-1);
#else
	idx << count << " " << headerpos << " " << seqpos << " " << (db->pos())-1 << '\n';
#endif
      } else {
	svi.push_back(seqpos,headerpos);
      }
    } else if (inseq) {
      if (opt.eos) {
        for (int ii=0;ii<opt.eos_len;ii++) {
          outputbuffer[outputi++] = opt.eos_char;
          seqpos++;
        }
	if (outputi >= (OUTPUTBUFSIZE-4)) {
	  if (!opt.tableonly) {
#if defined(_NO_LARGEFILE_STREAMS)
	    fwrite(outputbuffer,sizeof(char),outputi,seq);
#else
	    seq.write(outputbuffer,outputi);
#endif
	  }
	  outputi=0;
	}
	seqpos++;
      }
      if (!opt.tableonly) {
	if (!opt.binindex) {
#if defined(_NO_LARGEFILE_STREAMS)
	  fprintf(idx,"%llu %llu\n",seqpos,(db->pos())-1);
#else
	  idx << seqpos << " " << (db->pos())-1 << '\n';
#endif
	} else {
	  svi.push_back(seqpos,headerpos);
	}
      }
    }
    if (outputi > 0) {
      if (!opt.tableonly) {
#if defined(_NO_LARGEFILE_STREAMS)
	fwrite(outputbuffer,sizeof(char),outputi,seq);
#else
	seq.write(outputbuffer,outputi);
#endif
      }
      outputi=0;
    }
    delete db;
    if (!opt.tableonly) {
#if defined(_NO_LARGEFILE_STREAMS)
      fclose(seq);
#else
      seq.close();
#endif
#if defined(_NO_LARGEFILE_STREAMS)
      fclose(hdr);
#else
      hdr.close();
#endif
      if (!opt.binindex) {
#if defined(_NO_LARGEFILE_STREAMS)
	fclose(idx);
#else
	idx.close();
#endif
      }
      
      if (opt.binindex) {
	ofstream idb;
#if ! defined(__CYGWIN__) and ! defined(__MINGW32__)
	idb.open((opt.database+".idb").c_str());
#else 
	idb.open((opt.database+".idb").c_str(),ios::binary);
#endif
	svi.bwrite(idb);
	idb.close();
      }

    }
  
  }

  seq_time = modtime(opt.database+".seq");
  hdr_time = modtime(opt.database+".hdr");
  idx_time = modtime(opt.database+".idx");
  idb_time = modtime(opt.database+".idb");
  if (idx_time < idb_time) idx_time = idb_time;
 
  if (!opt.normalize && !opt.compress && !opt.tableonly) {
    if (opt.guard) {
      unlink((opt.database+".grd").c_str());
    }
    return 0;
  }

  if (didinitscan) {
    if (opt.checkonly) {
      return 1;
    }
    std::vector<int> chmap(256);
    std::vector<int> order(256);
    unsigned int index=0;
    ofstream tbl;
    ofstream tbz;
    if (opt.normalize || opt.tableonly) {
#if   ! defined(__CYGWIN__)
      tbl.open((opt.database+".tbl").c_str());
#else 
      tbl.open((opt.database+".tbl").c_str(),ios::binary);
#endif
    }
    if (opt.compress) {
#if   ! defined(__CYGWIN__)
      tbz.open((opt.database+".tbz").c_str());
#else 
      tbz.open((opt.database+".tbz").c_str(),ios::binary);
#endif      
    }
    for (int i=0;i<256;i++) {
      order[i] = i;
    }
    if (opt.dnaopt) {
      order[0] = 'A'; order['A'] = 0;
      order[1] = 'C'; order['C'] = 1;
      order[2] = 'G'; order['G'] = 2;
      order[3] = 'T'; order['T'] = 3;
    }
    for (int i=0;i<256;i++) {
      if (obs[order[i]]) {
		tbl << (unsigned char)order[i];
		tbz << (unsigned char)order[i];
		chmap[order[i]] = index;
		index++;
      }
    }
    tbl.close();
    tbz.close();
  }

  seq_time = modtime(opt.database+".seq");

  if (!opt.normalize && !opt.compress) {
    if (opt.suftree && seqst_time < seq_time) {
      // Build stuffix tree...
      make_suftree(opt.database,".seq",".seq.st");
    }
    if (opt.guard) {
      unlink((opt.database+".grd").c_str());
    }
    return 0;
  }

  tbl_time = modtime(opt.database+".tbl");
  tbz_time = modtime(opt.database+".tbz");

  if (opt.compress &&
      (opt.force ||
       (opt.guard && grd_time > 0) || 
       (seq_time > sqz_time) ||
       (tbz_time > sqz_time))
      ) {

    if (opt.checkonly) {
      return 1;
    }

    unsigned char invchmap[256];
    long unsigned int index;
    {
      std::string mapfn(opt.database+".tbz");
      ifstream chmapfile(mapfn.c_str());
      char chmap[256];
      chmapfile.read(chmap,256);
      index=chmapfile.gcount();
      chmapfile.close();
      for (int i=0;i<256;i++) {
	invchmap[i] = 255;
      }
      for (unsigned int i=0;i<index;i++) {
	invchmap[chmap[i]] = i;
      }
    }
   
    unsigned int bits = 1;
    while ( (((unsigned int)1) << bits) < index ) {
      bits++;
    } 
    unsigned int ucharbits = 8*sizeof(unsigned char);
    assert(bits <= ucharbits);
    unsigned int bufsize = least_common_multiple(bits,ucharbits)/8;
    bufsize *= 8;//sizeof(bigword);
    // cerr << bits << " " << ucharbits << " " << bufsize << endl;
    
    CharacterProducer* seqin;
    if (!opt.bufferedio) {
      seqin = new MapFileChars((opt.database+".seq").c_str(),'\n');
    } else {
      seqin = new BufferedFileChars((opt.database+".seq").c_str(),'\n');
    }
#if defined(_NO_LARGEFILE_STREAMS) 
    FILE *sqz = fopen((opt.database+".sqz").c_str(),"wb");
#else
#if ! defined(__CYGWIN__)
    ofstream sqz((opt.database+".sqz").c_str());
#else
    ofstream sqz((opt.database+".sqz").c_str(),ios::binary);
#endif
#endif
    
    FILE_POSITION_TYPE zposition;
    unsigned int bufposition;
    unsigned int bitposition;
    unsigned int position;
    position = zposition = bufposition = bitposition = 0;
    
    unsigned char *buffer = new unsigned char[bufsize];
    for (unsigned int i=0;i<bufsize;i++) {
      buffer[i] = 0;
    }

    outputi=0;

    bool partial_buffer = false;
    bool eofcache;
    while (!(eofcache=seqin->eof()) || partial_buffer) {
      partial_buffer = true;
      char ch;
      // cerr << ((eofcache)?"eofcache true":"eofcache false") << endl;
      // cerr << ((seqin->eof())?"seqin->eof() true":"seqin->eof() false") << endl;
      if (eofcache) {
	// checkpoint;
	ch = opt.eos_char;
      } else {
	// checkpoint;
	ch = seqin->getch();
      }
      // cerr << ((seqin->eof())?"seqin->eof() true":"seqin->eof() false") << endl;
      // cerr << "\"" << (unsigned int)ch << "\" ";
      // cerr << bufposition << "." << bitposition << " " << endl;
      unsigned char mask = invchmap[ch];
      // cerr << (int)mask << " " << (int)ch << endl;
      assert(mask < 255);
      // binary_format(cerr,mask);
      // cerr << endl;
      mask = mask << (ucharbits - bits);
      // binary_format(cerr,mask);
      // cerr << endl;
      if (bits + bitposition <= ucharbits) {
	mask = mask >> bitposition;
	// binary_format(cerr,mask);
	// cerr << endl;
	buffer[bufposition] |= mask;
	// binary_format(cerr,buffer[bufposition]);
	// cerr << endl;
	bitposition += bits;
      } else {
	unsigned char mask1 = mask;
	mask = mask >> bitposition;
	mask1 = mask1 << (bits - (bitposition + bits - ucharbits));
	// binary_format(cerr,mask);
	// cerr << " " << (int) mask << endl;
	// binary_format(cerr,mask1);
	// cerr << " " << (int) mask1 << endl;
	// binary_format(cerr,buffer[bufposition]);
	// cerr << " " << (int) buffer[bufposition] << endl;
	buffer[bufposition] |= mask;
	// binary_format(cerr,buffer[bufposition]);
	// cerr << " " << (int) buffer[bufposition] << endl;
	buffer[bufposition+1] |= mask1;
	// binary_format(cerr,buffer[bufposition+1]);
	// cerr << " " << (int) buffer[bufposition+1] << endl;
	bitposition += (bits - ucharbits);
	bufposition ++;
      }
      if (bufposition == bufsize-1 && bitposition == ucharbits) {
	// int z=0;
	for (unsigned int i=0;i<bufsize;i++) {
	  // binary_format(cerr,buffer[i],bits,z);
	  outputbuffer[outputi++] = buffer[i];
	  buffer[i] = 0;
	}
	if (outputi >= (OUTPUTBUFSIZE-bufsize)) {
#if defined(_NO_LARGEFILE_STREAMS)
	  fwrite(outputbuffer,sizeof(char),outputi,sqz);
#else
	  sqz.write(outputbuffer,outputi);
#endif
	  outputi=0;
	}
	// cerr << endl;
	bufposition = 0;
	bitposition = 0;
	zposition += bufsize;
	partial_buffer = false;
      }
      position++;
    }
    // int z=0;
    // for (unsigned int i=0;i<bufsize;i++) {
    // binary_format(cerr,buffer[i],bits,z);
    // sqz << buffer[i];
    // }
    // cerr << endl;
    
    if (outputi > 0) {
#if defined(_NO_LARGEFILE_STREAMS)
      fwrite(outputbuffer,sizeof(char),outputi,sqz);
#else
      sqz.write(outputbuffer,outputi);
#endif
      outputi=0;
    }
    delete [] buffer;
    
    delete seqin;
#if defined(_NO_LARGEFILE_STREAMS)
    fclose(sqz);
#else
    sqz.close();
#endif

  } 

  sqz_time = modtime(opt.database+".sqz");

  if (opt.compress && opt.suftree && sqzst_time < sqz_time) {
    // Build stuffix tree...
    make_suftree(opt.database,".sqz",".sqz.st");
  }

  if (opt.normalize &&
      (opt.force ||
       (opt.guard && grd_time > 0) || 
       (seq_time > sqn_time) ||
       (tbl_time > sqn_time))
      ) {

    if (opt.checkonly) {
      return 1;
    }

    unsigned char invchmap[256];
    long unsigned int index;
    {
      std::string mapfn(opt.database+".tbl");
      ifstream chmapfile(mapfn.c_str());
      char chmap[256];
      chmapfile.read(chmap,256);
      index=chmapfile.gcount();
      chmapfile.close();
      for (int i=0;i<256;i++) {
	invchmap[i] = 255;
      }
      for (unsigned int i=0;i<index;i++) {
	invchmap[chmap[i]] = i;
      }
    }
    CharacterProducer* seqin;
    if (!opt.bufferedio) {
      seqin = new MapFileChars((opt.database+".seq").c_str(),'\n');
    } else {
      seqin = new BufferedFileChars((opt.database+".seq").c_str(),'\n');
    }
#if defined(_NO_LARGEFILE_STREAMS)
    FILE *sqn = fopen((opt.database+".sqn").c_str(),"wb");
#else
#if ! defined(__CYGWIN__)
    ofstream sqn((opt.database+".sqn").c_str());
#else
    ofstream sqn((opt.database+".sqn").c_str(),ios::binary);
#endif
#endif

    char ch;

    outputi=0;
    while (!seqin->eof()) {
      ch = seqin->getch();
      outputbuffer[outputi++] = (unsigned char)invchmap[ch];
      if (outputi >= (OUTPUTBUFSIZE-4)) {
#if defined(_NO_LARGEFILE_STREAMS)
	fwrite(outputbuffer,sizeof(char),outputi,sqn);
#else
	sqn.write(outputbuffer,outputi);
#endif
	outputi=0;
      }
    }
    if (outputi > 0) {
#if defined(_NO_LARGEFILE_STREAMS)
      fwrite(outputbuffer,sizeof(char),outputi,sqn);
#else
      sqn.write(outputbuffer,outputi);
#endif
      outputi=0;
    }
    delete seqin;
#if defined(_NO_LARGEFILE_STREAMS)
    fclose(sqn);
#else
    sqn.close();
#endif
  }

  sqn_time = modtime(opt.database+".sqn");

  if (opt.normalize && opt.suftree && sqnst_time < sqn_time) {
    // Build stuffix tree...
    make_suftree(opt.database,".sqn",".sqn.st");
  }

  if (opt.cleanup && (opt.compress || opt.normalize) ) {
    unlink((opt.database+".seq").c_str());
  }

  if (opt.guard) {
    unlink((opt.database+".grd").c_str());
  }

  // checkpoint;
  return 0;
}

