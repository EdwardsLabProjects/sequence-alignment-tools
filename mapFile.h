#ifndef IBPEP_MMAPUTIL_H
#define IBPEP_MMAPUTIL_H

#include "types.h"

class MapFile {

 public:
  char *fileString;
  char *fileName;
  FILE_SIZE_TYPE mapsize;
  MapFile();
  MapFile(const char *fname);
  virtual ~MapFile();
  void createMap(const char *fname);
  char * const & file_string() const {
    return fileString;
  }
  char * const & filename() const {
    return fileName;
  }
  FILE_SIZE_TYPE map_size() const {
    return mapsize;
  }
  inline bool endOfFile() const {
    return(_pos>=mapsize);
  }
  inline bool endOfFile(const FILE_POSITION_TYPE pos) const {
    return(pos>=mapsize);
  }
  inline FILE_SIZE_TYPE size() const {
    return mapsize;
  }
  inline signed char readCharacter() const{
    if (_pos < mapsize)
      return fileString[_pos];
    else return(-1);
  }
  inline char readCharacter( FILE_POSITION_TYPE const indx) const {
    if (indx < mapsize-1) 
      return fileString[indx];
    else 
      return fileString[mapsize-1];
  }


  inline char getCharacter () {
#ifdef CHECK_GETCHARACTER
    if (_pos < mapsize)
      return fileString[_pos++];
    else
      return -1;
#else
    return  fileString[_pos++];
#endif
  }

  inline char getCharacter (FILE_POSITION_TYPE indx) const {
#ifdef CHECK_GETCHAR
    if (indx < mapsize - 1)
      return  fileString[indx];
    else
      return  fileString[mapsize-1];
#else
    return  fileString[indx];
#endif
  }

  int setPos(const FILE_POSITION_TYPE);
  FILE_POSITION_TYPE getPos() const{
    return _pos;
  }
  virtual void reset() {
    _pos=0;
  }
 protected:
  FILE_POSITION_TYPE _pos;
  
};

#endif  //  MMAPUTIL_H

