
#include <cstdio>
#include <iostream>
using namespace std;

#include "rl_suffix_tree.h"
#include "util.h"

class SufTree: public suffix_tree<basic_node,basic_suffix> {
public:
  SufTree(const char *st,unsigned len,char t='$'):
    suffix_tree<basic_node,basic_suffix>(st,len,t) {}
};

int main(int argc,const char **argv) {

  if (argc < 2) {
    exit(1);
  }

  string infile(argv[1]);
  string stfile = (infile+".st");

  FILE_POSITION_TYPE size = filesize(argv[1]);
  char *buffer=new char[size+1];
  FILE* in=fopen(infile.c_str(),"r");
  fread(buffer,sizeof(char),size,in);
  buffer[size]='\0';
  fclose(in);

  SufTree T(buffer,size,buffer[0]);
  if (!exist(stfile)) {
    T.build();
    T.write(stfile.c_str());
  } else {
    T.read((const char *)stfile.c_str());
  }
  // T.suffix_print(); 
  // T.print();
  T.check(cout);
  char *query="LQGSA";
  //  char *query="LQGSATAAEAfdjal;fjaj;f";
  unsigned poslen = 100;
  unsigned *pos=new unsigned[poslen];
  unsigned npos;
  while ((npos = T.find(query,pos,poslen))>poslen) {
    delete [] pos;
    pos = new unsigned[npos];
    poslen = npos;
  }
  cout << "npos " << npos << endl;
  for (int i=0;i<npos;i++) {
    cout << i << ' ' << pos[i] << endl;
  }
}
