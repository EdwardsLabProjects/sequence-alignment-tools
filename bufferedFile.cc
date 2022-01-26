//================================================================
// Brian Walenz's code. Provides methods for buffered I/O.
//
//===============================================================

#include "bufferedFile.h"
#include "math.h"
#include "util.h"

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

BufferedFile::BufferedFile() {
  int  buffersize = 256*1024*1024;
  _buffer = new char[buffersize];
  if (_buffer == 0L) {
    fprintf(stderr, "Couldn't allocate the file buffer\n");
    exit(1);
  }

  _bufferMax = buffersize;
  _bufferLen = 0;
  _readPos   = 0;
  _size = _pos = 0;
}

BufferedFile::BufferedFile(std::string const &infile, int buffersize) {

  _size = _pos = 0;

  if (buffersize <= 0) {
    buffersize = filesize(infile);
    if (buffersize <= 0) {
      // fprintf(stderr, "WARNING:  Illegal buffer size;resetting to 1024*1024.\n");
      buffersize = 256*1024*1024;
    }
  }

  _buffer = new char[buffersize];
  if (_buffer == 0L) {
    fprintf(stderr, "Couldn't allocate the file buffer\n");
    exit(1);
  }
  _bufferMax = buffersize;
  openBuf(infile);
}

void BufferedFile::openBuf(std::string const &fname) {
  // cerr << "file: " << fname << endl;
  struct stat st;
  if (stat(fname.c_str(),&st)<0) {
    cerr << "stat failed for file: " << fname << endl;
    perror("BufferedFile::openBuf:stat");
    exit(1);
  }
  _size = st.st_size;

#if ! defined(__CYGWIN__) 
  bufF = fopen(fname.c_str(),"r");
#else
  bufF = fopen(fname.c_str(),"rb");
#endif
  if (!bufF) {
    cerr << "couldn't open input file: " << fname << endl;
    exit(1);
  }
  _bufferLen = 0;
  _readPos   = 0;
}


BufferedFile::~BufferedFile() {
  delete [] _buffer;
  fclose(bufF);
}

void BufferedFile::reset() {
  _bufferLen = 0;
  _readPos   = 0;
  rewind(bufF);
  _pos = 0;
}

FILE_POSITION_TYPE BufferedFile::getPos() const {
  return _pos;
}

void BufferedFile::setPos(FILE_POSITION_TYPE p) {
  _pos = p;
  FILE_POSITION_TYPE tg = FTELL(bufF);
  if (tg > p && tg <= p + _bufferLen) {
    _readPos = p - (tg - _bufferLen);
  } else {
    // bufF.clear();
    FSEEK(bufF,p,SEEK_SET);
    _bufferLen = fread(_buffer,sizeof(unsigned char),_bufferMax,bufF);
    _readPos   = 0;
  }
}

unsigned char BufferedFile::getCharacter(void) {
  if (_readPos >= _bufferLen) {
    long signed int lookback=(long signed int)floor(((double)_bufferMax)/20);
    FILE_POSITION_TYPE tg = FTELL(bufF);
    if (tg > lookback) {
      FSEEK(bufF,-lookback,SEEK_CUR);
    } else {
      FSEEK(bufF,0,SEEK_SET);
      lookback=0;
    }
    _bufferLen = fread(_buffer,sizeof(unsigned char),_bufferMax,bufF);
    _readPos   = lookback;
  }
  _pos++;
  return(_buffer[_readPos++]);
}



