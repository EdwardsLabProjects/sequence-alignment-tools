#include "gzFile.h"
#include <iostream>

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

GZFile::GZFile() {
  _pos = 0;
  bufF = 0;
}

GZFile::GZFile(std::string const &infile) {
  bufF = 0;
  _pos = 0;
  openBuf(infile);
}

void GZFile::openBuf(std::string const &fname) {

  bufF = gzopen(fname.c_str(),"rb");
  
  if (!(*bufF)) {
    cerr << "couldn't open input file: " << fname << endl;
    exit(1);
  }
}

GZFile::~GZFile() {
  gzclose(bufF);
  _pos=0;
}

void GZFile::reset() {
  gzseek(bufF,0,SEEK_SET);
  _pos=0;
}

FILE_POSITION_TYPE GZFile::getPos() const {
  return _pos;
}

void GZFile::setPos(FILE_POSITION_TYPE p) {
  gzseek(bufF,p,SEEK_SET);
  _pos = p;
}

unsigned char GZFile::getCharacter(void) {
  _pos++;
  return gzgetc(bufF);
}




