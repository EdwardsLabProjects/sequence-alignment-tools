//================================================================
// Modified Brian Walenz's code. Provides methods for buffered I/O.
// VB 9/18/00
//===============================================================

#ifndef IBPEP_MS_NDX_filestar_h
#define IBPEP_MS_NDX_filestar_h

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "types.h"

class FileStarFile {
private:
  FILE* bufF;
  FILE_POSITION_TYPE _pos;
  FILE_SIZE_TYPE _size;
public:
  FileStarFile(std::string const &);
  FileStarFile();
  ~FileStarFile();
  
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

#endif  //  BUFFEREDFILE_H

