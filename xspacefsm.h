#ifndef ZSCOREFSM_H
#define ZSCOREFSM_H

#include <iostream>
using namespace std;

#include <stdio.h>
#include "rl_suffix_tree.h"
#include "types.h"

typedef suffix_tree<basic_node,basic_suffix> SufTree;

class XSpaceFSM {
  unsigned int mersize;
  char *S;
  unsigned Slen;
  SufTree *Tree;
  unsigned char *Z,*nZ;

 public:
  ~XSpaceFSM() {
    if (Z) delete[] Z;
    if (nZ) delete[] nZ;
    if (Tree) delete Tree;
    Tree = NULL; Z = NULL; nZ = NULL;
  }
  XSpaceFSM(char *s,unsigned len,unsigned int m): S(s),Slen(len),mersize(m) {
    Tree = new SufTree(S,Slen,S[0]);
    Tree->build();
  }

  void initialize();
  void selfstream();
  void stream(const char *map,FILE *f,long unsigned len);
  void output(FILE* o, FILE_POSITION_TYPE offset) const;

 private:
  bool interesting(st_index n,char cl,char cr) const;
  void selfprocess(st_index n) const;
  void printnode(st_index,bool,unsigned) const;
  void output_nodes(FILE* o, FILE_POSITION_TYPE offset, st_index n) const;
  void print_nodes(FILE *o, FILE_POSITION_TYPE offset, st_index n) const;
};

#endif
