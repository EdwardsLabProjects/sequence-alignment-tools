#ifndef INDEX_H
#define INDEX_H


#include <cstdio>
#include <iostream>
#include <string>
#include <list>

#include "types.h"

using namespace std;

/**
  Indexing conventions:
  start is the place in the file where there sequence character information
  starts.  All sequences start with a special terminal character and end at
  that terminal character.

  Thus, if I seek to start, the first thing I read will be a term,
  and if I seek to stop, the first thing I real will be a term
**/

class min_index_elt {
 public:
  FILE_POSITION_TYPE start,stop;
  min_index_elt():
    start(0),stop(0) {}
  inline int seek(FILE *f,bool ends = true) const {
    return FSEEK(f,start+(ends?0:1),SEEK_SET);
  }
};

class min_index_list: public list<min_index_elt> {
 public:
  min_index_list():list<min_index_elt>() {}

  int iload(FILE *,char flag=0);
};


class index_elt {
 public:
  FILE_POSITION_TYPE cstart,cstop,start,stop;
  string defline;

  unsigned seqno;
  char flags;

  index_elt():
    cstart(0),cstop(0), // read only fields
    start(0),stop(0),   // read only fields
    defline("NODEF"),   // read only fields
    seqno(unsigned(-1)),      // mutable field for the user
    flags(0)                  // mutable field for the user
    {}
  inline int seek_fasta(FILE *f) const {
    return FSEEK(f,cstart,SEEK_SET);
  }
  inline int seek(FILE *f,bool ends = true) const {
    return FSEEK(f,start+(ends?0:1),SEEK_SET);
  }
  void fprint(FILE *f) const;
};

class index_list : public list<index_elt> {
 public:
  index_list():list<index_elt>() {}

  int iload_fasta(FILE *,char flag=0);
  int iload(FILE *,char flag=0);
  int isave(FILE *) const;
};

class sequence : public index_elt {
 public:
  char *chars;

  sequence():index_elt(),chars(NULL){}
  sequence(const index_elt &e):index_elt(e),chars(NULL){}

  int sload_fasta(FILE *);
  int sload(FILE *,bool ends = true);
  int ssave(FILE *,bool ends = true) const;
};

class min_sequence : public min_index_elt {
 public:
  char *chars;

  min_sequence():min_index_elt(),chars(NULL){}
  min_sequence(const min_index_elt &e):min_index_elt(e),chars(NULL){}

  int sload(FILE *,bool ends = true);
};


#endif
