#include "fileStar.h"
#include <iostream>

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

FileStarFile::FileStarFile() {
  _size = _pos = 0;
}

FileStarFile::FileStarFile(std::string const &infile) {
  _size = _pos = 0;
  openBuf(infile);
}

void FileStarFile::openBuf(std::string const &fname) {

  struct stat st;
  if (stat(fname.c_str(),&st)<0) {
    perror("FileStarFile::openBuf:stat");
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
}

FileStarFile::~FileStarFile() {
  fclose(bufF);
  _size=_pos=0;
}

void FileStarFile::reset() {
  rewind(bufF);
  _pos=0;
}

FILE_POSITION_TYPE FileStarFile::getPos() const {
  return _pos;
}

void FileStarFile::setPos(FILE_POSITION_TYPE p) {
  FSEEK(bufF,p,SEEK_SET);
  _pos = p;
}

unsigned char FileStarFile::getCharacter(void) {
  _pos++;
  return fgetc(bufF);
}




