#ifndef _IBPEP_FASTA_IO_H_
#define _IBPEP_FASTA_IO_H_

#include <assert.h>
#include <string>

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

class fasta_entry {
  std::string defline_;
  std::string sequence_;
public:
  fasta_entry(std::string const & seq="", std::string const & dl="") {
    defline_ = dl;
    sequence_ = seq;
  }
  ~fasta_entry() {}
  std::string const & sequence() const {
    return sequence_;
  }
  std::string const & defline() const {
    return defline_;
  }
  void sequence(std::string const & seq) {
    sequence_ = seq;
  }
  void defline(std::string const & dl) {
    defline_ = dl;
  }
  void uppercase();
  void lowercase();
  void read(istream & is);
  void write(ostream & os) const;
};

ostream & operator<<(ostream & os, fasta_entry const & fe);
istream & operator>>(istream & is, fasta_entry & fe);

// #include "fasta_io.t"

#endif

