//================================================================
// Modified Brian Walenz's code. Provides methods for buffered I/O.
// VB 9/18/00
//===============================================================

#ifndef IBPEP_MS_NDX_BUFFEREDFILE_H
#define IBPEP_MS_NDX_BUFFEREDFILE_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include "types.h"

class BufferedFile {
private:
  FILE* bufF;
  char* _buffer;
  int               _bufferMax;
  int               _bufferLen;
  int               _readPos;
  FILE_POSITION_TYPE _pos;
  FILE_SIZE_TYPE     _size;
public:
  BufferedFile(std::string const &, int buffersize=1024*1024);
  BufferedFile();
  ~BufferedFile();
  
  unsigned char     getCharacter(void);
  void openBuf(std::string const &fname);
  bool              isEndOfFile(void) const {
    return(_pos>=_size);
  }
  FILE_POSITION_TYPE getPos() const ;
  void setPos(FILE_POSITION_TYPE p);
  FILE_SIZE_TYPE getSize() const {
    return _size;
  }
  void reset();
};

#endif  //  BUFFEREDFILE_H
