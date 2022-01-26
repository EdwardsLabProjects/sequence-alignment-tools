//===============================================================

#ifndef IBPEP_MS_NDX_gzstar_h
#define IBPEP_MS_NDX_gzstar_h

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "types.h"

#include "zlib.h"

class GZFile {
private:
  gzFile bufF;
  FILE_POSITION_TYPE _pos;
  FILE_SIZE_TYPE _size;
public:
  GZFile(std::string const &);
  GZFile();
  ~GZFile();
  
  unsigned char     getCharacter(void);
  void openBuf(std::string const &fname);
  inline bool isEndOfFile(void) const {
    return(_pos >= _size);
  }
  FILE_POSITION_TYPE getPos() const ;
  void setPos(FILE_POSITION_TYPE p);
  FILE_SIZE_TYPE getSize() const {
    return _size;
  }
  void reset();
};

#endif

