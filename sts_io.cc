
#include "sts_io.h"
#include <ctype.h>
#include <iostream>
#include <assert.h>
#include "util.h"

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

void sts_entry::read(istream & is) {
  // 8K buffer sufficient for a single line? Hope so. 
  static unsigned long BUFLEN = 1024*8;
  static char *buffer = new char[BUFLEN];
  static char *buffer1 = new char[BUFLEN];
  std::string size;
  is.getline(buffer,BUFLEN-2,'\n');
  if (is.gcount() >= BUFLEN-2) {
    timestamp("STS line too long!\n");
    exit(1);
  }
  buffer[is.gcount()] = '\0';
  istrstream iss(buffer);
  iss >> id_
      >> forward_primer_
      >> reverse_primer_
      >> size
      >> acc_
      >> chrom_
      >> altacc_;
  std::string::size_type p;
  if ((p=size.find("-"))!=std::string::npos) {
    sizelb_ = atoi(size.substr(0,p).c_str());
    sizeub_ = atoi(size.substr(p+1).c_str());
  } else {
    sizelb_ = sizeub_ = atoi(size.c_str());
  }
  // scoop up rest of line...
  iss.getline(buffer1,BUFLEN-1,'\n');
  assert(iss.gcount() < BUFLEN-1);
  buffer1[iss.gcount()] = '\0';
  species_ = buffer1;
}

void sts_entry::write(ostream & os) const {
  os << id_ << '\t'
     << forward_primer_ << '\t'
     << reverse_primer_ << '\t';
  if (sizelb_ == sizeub_) {
    os << sizelb_ << '\t';
  } else {
    os << sizelb_ << "-" << sizeub_ << '\t';
  }
  os << acc_ << '\t'
     << chrom_ << '\t'
     << altacc_ << '\t'
     << species_ << endl;
}

void sts_entry::uppercase() {
  for (int i=0;i<forward_primer_.length();i++) {
    if (forward_primer_[i] >= 'a' && forward_primer_[i] <= 'z') {
      forward_primer_[i] = toupper(forward_primer_[i]);
    }
  }
  for (int i=0;i<reverse_primer_.length();i++) {
    if (reverse_primer_[i] >= 'a' && reverse_primer_[i] <= 'z') {
      reverse_primer_[i] = toupper(reverse_primer_[i]);
    }
  }
}

void sts_entry::lowercase() {
  for (int i=0;i<forward_primer_.length();i++) {
    if (forward_primer_[i] >= 'a' && forward_primer_[i] <= 'z') {
      forward_primer_[i] = toupper(forward_primer_[i]);
    }
  }
  for (int i=0;i<reverse_primer_.length();i++) {
    if (reverse_primer_[i] >= 'a' && reverse_primer_[i] <= 'z') {
      reverse_primer_[i] = toupper(reverse_primer_[i]);
    }
  }
}

ostream & operator<<(ostream & os, sts_entry const & fe) {
  fe.write(os);
  return os;
}
istream & operator>>(istream & is, sts_entry & fe) {
  fe.read(is);
  return is;
}


