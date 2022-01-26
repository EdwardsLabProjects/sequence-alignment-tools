#ifndef _IBPEP_STS_IO_H_
#define _IBPEP_STS_IO_H_

#include <assert.h>
#include <string>

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

class sts_entry {
  std::string id_;
  std::string forward_primer_;
  std::string reverse_primer_;
  long unsigned int sizelb_;
  long unsigned int sizeub_;
  std::string acc_;
  std::string chrom_;
  std::string altacc_;
  std::string species_;
public:
  sts_entry(std::string const & id="", 
	    std::string const & fp="", std::string const & rp="",
	    long unsigned int sizelb=0, long unsigned int sizeub=0, 
	    std::string const & accession="",std::string const & chrom="",
	    std::string const & altacc="", std::string const & species="") {
    id_ = id;
    forward_primer_ = fp;
    reverse_primer_ = rp;
    sizelb_ = sizelb;
    sizeub_ = sizeub;
    acc_ = accession;
    chrom_ = chrom;
    altacc_ = altacc;
    species_ = species;
  }
  ~sts_entry() {}
  std::string const & id() const {
    return id_;
  }
  void id(std::string const & id) {
    id_ = id;
  }
  std::string const & forward_primer() const {
    return forward_primer_;
  }
  void forward_primer(std::string const & fp) {
    forward_primer_ = fp;
  }
  std::string const & reverse_primer() const {
    return reverse_primer_;
  }
  void reverse_primer(std::string const & rp) {
    reverse_primer_ = rp;
  }
  long unsigned int const & sizelb() const {
    return sizelb_;
  }
  void sizelb(long unsigned int const & sizelb) {
    sizelb_ = sizelb;
  }
   long unsigned int const & sizeub() const {
    return sizeub_;
  }
  void sizeub(long unsigned int const & sizeub) {
    sizeub_ = sizeub;
  }
  std::string const & accession() const {
    return acc_;
  }
  void accession(std::string const & acc) {
    acc_ = acc;
  }
  std::string const & altacc() const {
    return altacc_;
  }
  void altacc(std::string const & altacc) {
    altacc_ = altacc;
  }
  std::string const & species() const {
    return species_;
  }
  void species(std::string const & species) {
    species_ = species;
  }
  std::string const & chrom() const {
    return chrom_;
  }
  void chrom(std::string const & chrom) {
    chrom_ = chrom;
  }
  void uppercase();
  void lowercase();
  void read(istream & is);
  void write(ostream & os) const;
};

ostream & operator<<(ostream & os, sts_entry const & fe);
istream & operator>>(istream & is, sts_entry & fe);

#endif

