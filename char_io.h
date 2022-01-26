
#ifndef _IBPEP_CHAR_IO_H_
#define _IBPEP_CHAR_IO_H_

#include "mapFile.h"
#include "bufferedFile.h"
#include "fileStar.h"
#include "util.h"
#if !defined(NO_ZLIB)
#include <zlib.h>
#endif
#include <assert.h>

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

class CharacterProducer {
public:
  CharacterProducer() {};
  virtual ~CharacterProducer() {};
  virtual char getch()=0;
  virtual unsigned char getnch()=0;
  virtual char ch(unsigned char)=0;
  virtual int nch(char)=0;
  string getstr(char eos, bool keep=true) {
    char ch;
    string str = "";
    while (!eof() && ((ch=getch()) != eos)) {
      str += ch;
    }
    if (keep && ch == eos) {
      str += ch;
    }
    return str;
  }
  string getstr(unsigned int n) {
    char ch;
    string str = "";
    for (unsigned int i=0;i<n;i++) {
      if (eof()) break;
      ch = getch();
      str += ch;
    }
    return str;    
  }
  virtual unsigned int size() const =0;
  virtual FILE_POSITION_TYPE length() const =0;
  virtual bool eof() const =0;
  virtual FILE_POSITION_TYPE pos() const =0;
  virtual void pos(FILE_POSITION_TYPE p) =0;
  virtual void reset()=0;
  virtual float progress() const =0;
  virtual inline bool has_filename() const { return false; }
  virtual inline char const * c_str() const { assert(0); return NULL; }
  virtual inline char const * filename() const { assert(0); return NULL; }
  void getbasepos(FILE_POSITION_TYPE, FILE_POSITION_TYPE &, int &) const {
    assert(0);
  };
  char getbasech() {
    assert(0);
    return 0;
  }
  virtual inline char getcodonid() {
    assert(0);
    return 0;
  }
  void mapto(char f, char t) {
    assert(0);
  };
};

class IStreamChars : public CharacterProducer {
private:
  bool delstream;
  std::istream & is_;
  unsigned char cached_char;
public:
  IStreamChars(const std::string & filename,char,int=0) : is_(*(new std::ifstream(filename.c_str()))) {
    delstream = true;
    cached_char = is_.get();
  }
  IStreamChars(std::istream & is) : is_(is) {
    delstream = false;
    cached_char = is_.get();
  }
  ~IStreamChars() {
    if (delstream) {
      delete (&is_);
    }
  }
  char getch() { char tmp(cached_char); cached_char = is_.get(); return tmp; }
  unsigned char getnch() { return (unsigned char) getch(); }
  inline char ch(unsigned char nch) { return nch; }
  inline int nch(char ch) { return ch; }
  inline unsigned int size() const { return 256; }
  inline FILE_POSITION_TYPE length() const { return 0; }
  bool eof() const { return is_.eof(); }
  FILE_POSITION_TYPE pos() const { return ((FILE_POSITION_TYPE)(is_.tellg())-1); }
  void pos(FILE_POSITION_TYPE p) { is_.clear(); is_.seekg(p); cached_char = is_.get();}
  void reset() { is_.clear(); is_.seekg(0); }
  float progress() const { return 0.0; }
};

#if !defined(NO_ZLIB)
class GZChars : public CharacterProducer {
private:
  gzFile gzf;
  unsigned char cached_char;
public:
  GZChars(const std::string & filename,char,int=0) {
    gzf = gzopen(filename.c_str(),"rb");
    cached_char = gzgetc(gzf);
  }
  char getch() { char tmp(cached_char); cached_char = gzgetc(gzf); return tmp; }
  unsigned char getnch() { return (unsigned char) getch(); }
  inline char ch(unsigned char nch) { return nch; }
  inline int nch(char ch) { return ch; }
  inline unsigned int size() const { return 256; }
  inline FILE_POSITION_TYPE length() const { return 0; }
  bool eof() const { return gzeof(gzf); }
  FILE_POSITION_TYPE pos() const { return ((FILE_POSITION_TYPE)gztell(gzf)-1); }
  void pos(FILE_POSITION_TYPE p) { gzseek(gzf,p,SEEK_SET); cached_char = gzgetc(gzf);}
  void reset() { pos(0); }
  float progress() const { return 0.0; }
};
#endif

class CharStarChars : public CharacterProducer {
  const char *ptr_;
  const char *p_;
public:
  CharStarChars(const char *chars,char,int=0) : ptr_(chars), p_(chars) {}
  ~CharStarChars() {}
  inline char getch() { return *p_++; }
  inline unsigned char getnch() { return *p_++; }
  inline char ch(unsigned char nch) { return nch; }
  inline int nch(char ch) { return ch; }
  inline unsigned int size() const { return 256; }
  inline FILE_POSITION_TYPE length() const { return 0; }
  inline bool eof() const { return (*p_=='\0'); }
  FILE_POSITION_TYPE pos() const { return (p_-ptr_); }
  void pos(FILE_POSITION_TYPE p) { p_ = ptr_+p; }
  void reset() { p_ = ptr_; }
  float progress() const { 
    return 0; 
  }
};

class MapFileChars : public MapFile, public CharacterProducer {
public:
  MapFileChars(std::string const & filename,char,int=0,int=0) : MapFile(filename.c_str()) {}
  ~MapFileChars() {}
  inline char getch() { return (char) MapFile::getCharacter(); }
  inline unsigned char getnch() { return MapFile::getCharacter(); }
  inline char ch(unsigned char nch) { return nch; }
  inline int nch(char ch) { return ch; }
  inline unsigned int size() const { return 256; }
  inline FILE_POSITION_TYPE length() const { return MapFile::map_size(); }
  inline bool eof() const { return MapFile::endOfFile(); }
  FILE_POSITION_TYPE pos() const { return MapFile::getPos(); }
  void pos(FILE_POSITION_TYPE p) { MapFile::setPos(p); }
  void reset() { MapFile::reset(); }
  float progress() const { 
    return ((float)MapFile::getPos())/MapFile::map_size(); 
  }
  inline char const * c_str() const { return MapFile::file_string(); }
  inline char const * filename() const { return MapFile::filename(); }
  inline bool has_filename() const { return true; }
};

class BufferedFileChars : public BufferedFile, public CharacterProducer {
public:
  BufferedFileChars(std::string const & filename,char,int=0,int sizehint=0) : 
    BufferedFile(filename,sizehint) {}
  ~BufferedFileChars() {}
  inline char getch() { return (char) BufferedFile::getCharacter(); }
  inline unsigned char getnch() { return BufferedFile::getCharacter(); }
  inline char ch(unsigned char nch) { return nch; }
  inline int nch(char ch) { return ch; }
  inline unsigned int size() const { return 256; }
  inline FILE_POSITION_TYPE length() const { return BufferedFile::getSize(); }
  inline bool eof() const { return BufferedFile::isEndOfFile(); }
  FILE_POSITION_TYPE pos() const { return BufferedFile::getPos(); }
  void pos(FILE_POSITION_TYPE p) { BufferedFile::setPos(p); }
  void reset() { BufferedFile::reset(); }
  float progress() const { 
    return ((float)BufferedFile::getPos())/BufferedFile::getSize(); 
  }
};

class FileStarChars : public FileStarFile, public CharacterProducer {
public:
  FileStarChars(std::string const & filename,char,int=0) : FileStarFile(filename) {}
  ~FileStarChars() {}
  inline char getch() { return (char) FileStarFile::getCharacter(); }
  inline unsigned char getnch() { return FileStarFile::getCharacter(); }
  inline char ch(unsigned char nch) { return nch; }
  inline int nch(char ch) { return ch; }
  inline unsigned int size() const { return 256; }
  inline FILE_POSITION_TYPE length() const { return FileStarFile::getSize(); }
  inline bool eof() const { return FileStarFile::isEndOfFile(); }
  FILE_POSITION_TYPE pos() const { return FileStarFile::getPos(); }
  void pos(FILE_POSITION_TYPE p) { FileStarFile::setPos(p); }
  void reset() { FileStarFile::reset(); }
  float progress() const { 
    return ((float)FileStarFile::getPos())/FileStarFile::getSize(); 
  }

};

class CompressedCharUtil {
public:
  static unsigned int lcm(unsigned int, unsigned int);
};

// #include "char_io.t"

#endif

