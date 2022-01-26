#ifndef _IBPEP_TYPES_H_
#define _IBPEP_TYPES_H_

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>

#if defined(__linux) && defined(_FILE_OFFSET_BITS) && _FILE_OFFSET_BITS == 64
    typedef __off64_t FILE_POSITION_TYPE;
    typedef __off64_t FILE_SIZE_TYPE;
#   define FSEEK fseeko
#   define FTELL ftello
#   include <sstream>
#   define istrstream istringstream
#   define ostrstream ostringstream
#   include <values.h>
    typedef long long unsigned int bigword;
#   define ASCTIME(s,t) char s[1024]; asctime_r(t,s)
#   define GETPAGESIZE getpagesize()
#   define DRAND48 drand48()
#elif defined(__linux)
    typedef __off_t FILE_POSITION_TYPE;
    typedef __off_t FILE_SIZE_TYPE;
#   define FSEEK fseek
#   define FTELL ftell
#   include <sstream>
#   define istrstream istringstream
#   define ostrstream ostringstream
#   include <values.h>
    typedef long long unsigned int bigword;
#   define ASCTIME(s,t) char *s; s = asctime(t)
#   define GETPAGESIZE getpagesize()
#   define DRAND48 drand48()
#elif defined(__alpha)
    typedef off_t FILE_POSITION_TYPE;
    typedef off_t FILE_SIZE_TYPE;
#   define FSEEK fseek
#   define FTELL ftell
#   include <iostream>
#   include <strstream>
#   include <values.h>
    typedef long unsigned int bigword;
#   define ASCTIME(s,t) char s[1024]; s = asctime_r(t,s)
#   define GETPAGESIZE getpagesize()
#   define DRAND48 drand48()
#elif defined(_AIX)
    typedef off_t FILE_POSITION_TYPE;
    typedef off_t FILE_SIZE_TYPE;
#   define FSEEK fseek
#   define FTELL ftell
#   include <strstream>
#   include <values.h>
    typedef long unsigned int bigword;
#   define ASCTIME(s,t) char s[1024]; asctime_r(t,s)
#   define GETPAGESIZE getpagesize()
#   define DRAND48 drand48()
#elif defined(__CYGWIN__)
    typedef off_t FILE_POSITION_TYPE;
    typedef off_t FILE_SIZE_TYPE;
#   define FSEEK fseek
#   define FTELL ftell
#   include <sstream>
#   define istrstream istringstream
#   define ostrstream ostringstream
#   define MAXINT ((int)(~(1U<<31)))
    typedef long long unsigned int bigword;
#   define ASCTIME(s,t) char s[1024]; asctime_r(t,s)
#   define GETPAGESIZE getpagesize()
#   define DRAND48 drand48()
#elif defined(__MINGW32__)
    typedef off_t FILE_POSITION_TYPE;
    typedef off_t FILE_SIZE_TYPE;
#   define FSEEK fseek
#   define FTELL ftell
#   include <sstream>
#   define istrstream istringstream
#   define ostrstream ostringstream
#   define MAXINT ((int)(~(1U<<31)))
    typedef long long unsigned int bigword;
#   define ASCTIME(s,t) char *s; s = asctime(t)
#   define GETPAGESIZE 8192
#   define DRAND48 ((double)rand()/RAND_MAX)
#   define srand48(s) srand(s)
#elif defined(__sun)
    typedef off_t FILE_POSITION_TYPE;
    typedef off_t FILE_SIZE_TYPE;
#   define FSEEK fseek
#   define FTELL ftell
#   include <iostream>
#   include <strstream>
#   include <values.h>
    typedef long unsigned int bigword;
#   define ASCTIME(s,t) char *s; s = asctime(t)
#   define GETPAGESIZE getpagesize()
#   define DRAND48 drand48()
#elif defined(__APPLE__)
    typedef fpos_t FILE_POSITION_TYPE;
    typedef size_t FILE_SIZE_TYPE;
#   define FSEEK fseek
#   define FTELL ftell
#   include <sstream>
#   define istrstream istringstream
#   define ostrstream ostringstream
#   define MAXINT ((int)(~(1U<<31)))
    typedef long long unsigned int bigword;
#   define ASCTIME(s,t) char *s; s = asctime(t)
#   define GETPAGESIZE getpagesize()
#   define DRAND48 drand48()
#endif


int compare(FILE_POSITION_TYPE const &, FILE_POSITION_TYPE const &);
#endif
