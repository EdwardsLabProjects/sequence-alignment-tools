
#include "fasta_io.h"
#include <ctype.h>
#include <iostream>
#include <cstdio>
#include <assert.h>

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

void fasta_entry::read(istream & is) {
  // 8K buffer sufficient for a single line? Hope so. 
  static unsigned long BUFLEN = 1024*8;
  static char *buffer = new char[BUFLEN];
  int peek;

  peek=is.peek();
  while (peek!=EOF && (peek=='#' || peek=='\n')) {
    // Read and discard line from the stream...
    is.getline(buffer,BUFLEN-1,'\n');
    // assert(is.gcount() < BUFLEN-1);
    peek=is.peek();
  }

  if (peek == EOF) {
    peek=is.get();
    defline_ = "";
    sequence_ = "";
    return;
  }

  // assert(peek == '>');
  is.getline(buffer,BUFLEN-1,'\n');
  // assert(is.gcount() < BUFLEN-1);
  // Ensure the string is null terminated
  buffer[is.gcount()] = '\0';
  defline_ = (buffer+1);

  // Read the sequence
  sequence_="";
  peek=is.peek();
  while (peek!=EOF && peek!='>' && peek!='#' && peek!='\n') {
    is.getline(buffer,BUFLEN-1,'\n');
    // assert(is.gcount() < BUFLEN-1);
    buffer[is.gcount()] = '\0';
    sequence_ += buffer;
    peek=is.peek();
  }
  
  while (peek!=EOF && (peek=='#' || peek=='\n')) {
    // Read and discard line from the stream...
    is.getline(buffer,BUFLEN-1,'\n');
    // assert(is.gcount() < BUFLEN-1);
    peek=is.peek();
  }
}

void fasta_entry::write(ostream & os) const {
  static int LINELEN=60;
  os << '>' << defline_ << endl;
  for (unsigned int i=0;i<sequence_.length();i+=LINELEN) {
    os << sequence_.substr(i,LINELEN) << endl;
  }
}

void fasta_entry::uppercase() {
  for (unsigned int i=0;i<sequence_.length();i++) {
    if (sequence_[i] >= 'a' && sequence_[i] <= 'z') {
      sequence_[i] = toupper(sequence_[i]);
    }
  }
}

void fasta_entry::lowercase() {
  for (unsigned int i=0;i<sequence_.length();i++) {
    if (sequence_[i] >= 'A' && sequence_[i] <= 'Z') {
      sequence_[i] = tolower(sequence_[i]);
    }
  }
}

ostream & operator<<(ostream & os, fasta_entry const & fe) {
  fe.write(os);
  return os;
}
istream & operator>>(istream & is, fasta_entry & fe) {
  fe.read(is);
  return is;
}


