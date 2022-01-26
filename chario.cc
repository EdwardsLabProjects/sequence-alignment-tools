
#include "types.h"
#include "char_io.h"

#if defined(__linux)
# include <unistd.h>
#endif /* defined(__linux) */

#include <sys/stat.h>
#include "char_io.t"
#include "fasta_io.t"
#include <assert.h>

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

void usage(char *message=NULL) {
  if (message != NULL && strlen(message) > 0) {
    cerr << message << endl;
    cerr << endl;
  }
  cerr << "Usage: chario [options] \n\n";
  cerr << "Options: \n";
  cerr << "  -i <character-file>    Input file.\n";
  cerr << "  -h                     Command line option help.\n";
  cerr << "\n";
  exit(1);
}

struct Options {
  std::string char_file;
};

void options(int argc, char *argv[], Options & opt) {
  signed char c;
  optarg = NULL;
  while ((c = getopt(argc, argv, "i:h")) != -1) {
    switch (c) {
    case 'i':
      opt.char_file = optarg;
      break;
    case 'h':
    default :
      usage();
    }
  }
  if (opt.char_file == "") usage();
}

int main(int argc,char *argv[]) {

  assert(sizeof(FILE_POSITION_TYPE)>=8);
  
  Options opt;
  options(argc,argv,opt);

  timestamp("File open");
  // MapFileChars cp(opt.char_file.c_str()); timestamp("MapFileChars");
  // IStreamChars cp(opt.char_file.c_str()); timestamp("IStreamChars");
  // BufferedFileChars cp(opt.char_file.c_str()); timestamp("BufferedFileChars");
  // FileStarChars cp(opt.char_file.c_str()); timestamp("FileStarChars");
  // Compressed<MapFileChars> cp(opt.char_file.c_str()); timestamp("Compressed<MapFileChars>");
  // Normalized<MapFileChars> cp(opt.char_file.c_str()); timestamp("Normalized<MapFileChars>");
  // IndexedFastaFile<MapFileChars,Header> cp(opt.char_file.c_str(),true);
  // StreamedFastaFile<MapFileChars,Header_SI> cp(opt.char_file.c_str(),true,true,'\n');
  Translated<MapFileChars> cp(opt.char_file.c_str(),'$'); timestamp("Translated<MapFileChars>");

  char ch=0,ch1=0;

  timestamp("Start scan");
  tic;
  FILE_SIZE_TYPE size=0;
  while (!cp.eof()) {
    ch = cp.getch();
    // This attempts to ensure the compiler doesn't optimize out the getnch();
    if (ch > ch1) {
      ch1 = ch;
    }
    size++;
    cout << ch;
    if (size%50==0) {
      cout << " " << cp.pos() << " " << cp.basepos();
      FILE_POSITION_TYPE p = cp.pos();
      cp.pos(p);
      cout << " " << cp.pos() << " " << cp.basepos() << endl;
    }
  }
  cout << " " << cp.pos() << " " << cp.basepos();
  FILE_POSITION_TYPE p = cp.pos();
  cp.pos(p);
  cout << " " << cp.pos() << " " << cp.basepos() << endl;
  cout << '\n';
  timestamp("End scan");
  time_t tsec;
  tocassgn(tsec);
  double size1;
  if (size > 1024*1024*1024) {
    size1 = ((double)size)/(1024*1024*1024);
    fprintf(stderr,"File size: %.2f GB\n",size1);
  } else if (size > 1024*1024) {
    size1 = ((double)size)/(1024*1024);
    fprintf(stderr,"File size: %.2f MB\n",size1);
  } else if (size > 1024) {
    size1 = ((double)size)/(1024);
    fprintf(stderr,"File size: %.2f kB\n",size1);
  } else {
    size1 = ((double)size);
    fprintf(stderr,"File size: %.2f bytes\n",size1);
  }
  double cps = ((double)size)/tsec;
  if (cps > 1024*1024*1024) {
    cps /= 1024*1024*1024;
    fprintf(stderr,"Scan rate: %.2f GB/s, %.2f Gb/s\n",cps,cps*8);
  } else if (cps > 1024*1024) {
    cps /= 1024*1024;
    fprintf(stderr,"Scan rate: %.2f MB/s, %.2f Mb/s\n",cps,cps*8);
  } else if (cps > 1024) {
    cps /= 1024;
    fprintf(stderr,"Scan rate: %.2f kB/s, %.2f kb/s\n",cps,cps*8);
  } else {
    fprintf(stderr,"Scan rate: %.2f B/s, %.2f b/s\n",cps,cps*8);
  }

  timestamp("First 50 chars");
  int count=0;
  // cerr << cp.pos() << endl;
  cp.pos((FILE_POSITION_TYPE)0);
  // cerr << cp.pos() << endl;
  while (!cp.eof() && count < 50) {
    ch = cp.getch();
    fprintf(stderr,"%c",ch);
    count++;
  }
  fprintf(stderr,"\n");

  timestamp("Middle 50 chars");
  // cerr << cp.pos() << endl;
  FILE_POSITION_TYPE newpos=size/2;
  if (newpos > 25) {
    newpos -= 25;
  } else {
    newpos = 0;
  }
  cp.pos(newpos);
  // cerr << cp.pos() << endl;
  count=0;
  while (!cp.eof() && count < 50) {
    ch = cp.getch();
    fprintf(stderr,"%c",ch);
    count++;
  }
  fprintf(stderr,"\n");

  timestamp("Last 50 chars");
  if (size <= 50) {
    cp.pos(0);
  } else {
    cp.pos((FILE_POSITION_TYPE)(size-50));
  }
  while (!cp.eof()) {
    ch = cp.getch();
    fprintf(stderr,"%c",ch);
  }
  fprintf(stderr,"\n");

  srand48(time(NULL));
  int nseek=1000;
  double s0;
  FILE_POSITION_TYPE s;
  timestamp("Start random seeks");
  tic;
  for (int i=0; i<nseek; i++) {
    s0 = drand48();
    s = ((FILE_POSITION_TYPE)(s0*size));
    cp.pos(s);
    ch = cp.getch();
    // This attempts to ensure the compiler doesn't optimize out the getnch();
    if (ch > ch1) {
      ch1 = ch;
    }
  }
  timestamp("End random seeks")
  tocassgn(tsec);
  double sps = ((double)nseek)/tsec;
  fprintf(stderr,"Seeks per second: %.2f\n",sps);

  timestamp("50 random seeks")
  for (int i=0; i<50; i++) {
    s0 = drand48();
    s = ((FILE_POSITION_TYPE)(s0*size));
    cp.pos(s);
    ch = cp.getch();
    fprintf(stderr,"%c",ch);
  }
  fprintf(stderr,"\n");

}

