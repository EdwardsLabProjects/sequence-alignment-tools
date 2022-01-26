
#include <stdio.h>

/* posix compliant mmap */
# include <fcntl.h>
# include <errno.h>
# include <sys/types.h>
# include <unistd.h>
# include <stdlib.h>
# include <string.h>
# include <sys/stat.h>
# ifndef __USE_MISC
#   define __USE_MISC
# endif
# if ! defined(__MINGW32__)
#   include <sys/mman.h>
#   if !defined(MAP_SHARED)
#     error MAP_SHARED
#   endif
#endif

#ifdef __alpha
#define MMAPFLAGS    (MAP_FILE | MAP_VARIABLE | MAP_SHARED)
#endif

#ifdef __linux
#define MMAPFLAGS    (MAP_FILE | MAP_SHARED)
#endif

#ifdef __FreeBSD__
#define MMAPFLAGS    (MAP_FILE | MAP_SHARED)
#endif

#ifdef _AIX
#define MMAPFLAGS    (MAP_FILE | MAP_SHARED)
#endif

#ifdef __CYGWIN__
#define MMAPFLAGS    (MAP_FILE | MAP_SHARED)
#endif

#ifdef __sun
#define MMAPFLAGS    (MAP_SHARED)
#endif

#ifdef __APPLE__
#define MMAPFLAGS    (MAP_FILE | MAP_SHARED)
#endif

#include "mapFile.h"

/* posix mmap */

MapFile::MapFile() {
  fileString=0;
  fileName = 0;
  _pos=0;
  mapsize=0;
}

MapFile::MapFile(const char *fname) {
  errno = 0;
  fileName = 0;

  int f = open(fname, O_RDONLY);
  if ((f < 0) || (errno)) {
    fprintf(stderr, "Couldn't open '%s'\n", fname);
    perror("open");
    exit(1);
  }

  struct stat  sb;
  fstat(f, &sb);
  if (errno) {
    fprintf(stderr, "Couldn't stat '%s'\n", fname);
    perror("fstat\n");
    exit(1);
  }

  /* if (sizeof(FILE_SIZE_TYPE) > sizeof(void*) &&
      ((((FILE_SIZE_TYPE)1)<<(8*sizeof(void*)))-1) < sb.st_size) {
    fprintf(stderr, "Couldn't map '%s'\n", fname);
    exit(1);
    }*/

  mapsize = sb.st_size;

  // fprintf(stderr, "%40.40s: Mapping %lu bytes.\n", fname, mapsize);
  // fflush(stderr);

# ifndef __MINGW32__
  fileString = (char *)mmap(0L, mapsize, PROT_READ, MMAPFLAGS, f, 0);
# else
  fileString = new char[mapsize];
  read(f,fileString,mapsize);
# endif

  if (errno) {
    fprintf(stderr, "Couldn't map '%s'\n", fname);
    perror("mmap");
    exit(1);
  }

  _pos=0;

  close(f);
  fileName = new char[strlen(fname)+1];
  strcpy(fileName,fname);
}

void MapFile::createMap(const char *fname) {
  errno = 0;

  int f = open(fname, O_RDONLY);
  if ((f < 0) || (errno)) {
    fprintf(stderr, "Couldn't open '%s'\n", fname);
    perror("open");
    exit(1);
  }

  struct stat  sb;
  fstat(f, &sb);
  if (errno) {
    fprintf(stderr, "Couldn't stat '%s'\n", fname);
    perror("fstat\n");
    exit(1);
  }

  /* if (sizeof(FILE_SIZE_TYPE) > sizeof(void*) &&
      ((((FILE_SIZE_TYPE)1)<<(8*sizeof(void*)))-1) < sb.st_size) {
    fprintf(stderr, "Couldn't map '%s'\n", fname);
    exit(1);
    } */

  mapsize = sb.st_size;

  // fprintf(stderr, "%40.40s: Mapping %lu bytes.\n", fname, mapsize);
  // fflush(stderr);

# ifndef __MINGW32__
  fileString = (char *)mmap(0L, mapsize, PROT_READ, MMAPFLAGS, f, 0);
# else
  fileString = new char[mapsize];
  read(f,fileString,mapsize);
# endif

  if (errno) {
    fprintf(stderr, "Couldn't map '%s'\n", fname);
    perror("mmap");
    exit(1);
  }

  _pos=0;
  
  close(f);
  fileName = new char[strlen(fname)+1];
  strcpy(fileName,fname);
}

MapFile::~MapFile() {
  if (fileName) {
    delete [] fileName;
  }
  if (fileString) {
#ifdef __CYGWIN__
    (void)munmap((char*)fileString, mapsize);
#elif defined(__MINGW32__) 
    delete [] fileString;
#else
    (void)munmap(fileString, mapsize);
#endif
  }
}

int MapFile::setPos(const FILE_POSITION_TYPE pos) {
  _pos=pos;
  return (pos>=mapsize);  
}









