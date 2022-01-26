
#ifndef WORD_GRAPH_H_
#define WORD_GRAPH_H_

#include "memory_debug.h"
#include "graph.h"
#include "char_io.h"
#include "fasta_io.h"
#include "fasta_io.t"
#include "types.h"

class word_graph_edge;
class word_graph_node : public labelnode<word_graph_edge> {
 public:
  virtual ~word_graph_node() {};
  std::string sequence(CharacterProducer & cp) const {
    std::string retval;
    retval.resize(length());
    cp.pos(seq_start());
    for (unsigned short i=0;i<length();i++) {
      retval[i] = cp.getch();
    }
    return retval;
  }
  FILE_POSITION_TYPE seq_start() const {
    return seq_end() - length();
  }
  virtual unsigned short const & length() const=0;
  FILE_POSITION_TYPE seq_end() const;
  void dump(ostream & os, CharacterProducer & cp, bool mark=false) const; 
  MEMORY_DEBUG(word_graph_node)
};

class varlen_word_graph_node : public word_graph_node {
  unsigned short length_;
public:
  varlen_word_graph_node() : length_(0) {};
  virtual ~varlen_word_graph_node() {};
  unsigned short const & length() const {
    return length_;
  }
  void length(unsigned short const & l) {
    length_ = l;
  }
  MEMORY_DEBUG(varlen_word_graph_node)
};

class fixedlen_word_graph_node : public word_graph_node {
  static unsigned short length_;
public:
  fixedlen_word_graph_node() {};
  virtual ~fixedlen_word_graph_node() {};
  unsigned short const & length() const {
    return length_;
  }
  static void length(unsigned short const & l) {
    length_ = l;
  }
  MEMORY_DEBUG(fixedlen_word_graph_node)
};

class word_graph_edge : public labeledge<word_graph_node> {
 public:
  word_graph_edge() {}
  virtual ~word_graph_edge() {};
  virtual unsigned length() const=0;
  virtual std::string sequence(CharacterProducer & cp)=0;
  virtual bool isreal() const=0;
  virtual unsigned int count() const=0;
  MEMORY_DEBUG(word_graph_edge)
};

class real_word_graph_edge : public word_graph_edge {
  unsigned length_;
  FILE_POSITION_TYPE seq_end_;
 public:
  real_word_graph_edge() : seq_end_(0), length_(0) {}
  virtual ~real_word_graph_edge() {};
  FILE_POSITION_TYPE seq_start() const {
    return seq_end_-length_;
  }
  FILE_POSITION_TYPE seq_end() const {
    return seq_end_;
  }
  void seq_end(FILE_POSITION_TYPE const & p) {
    seq_end_ = p;
  }
  bool isreal() const { return true; }
  std::string sequence(CharacterProducer & cp) {
    std::string retval;
    retval.resize(length_);
    FILE_POSITION_TYPE p = seq_end_ - length_;
    cp.pos(p);
    for (unsigned int i=0;i<length_;retval[i++]=cp.getch());
    return retval;
  }
  unsigned length() const {
    return length_;
  }
  void length(FILE_POSITION_TYPE const & l) {
    length_ = l;
  }
  unsigned int count() const {
    return 0;
  };
  MEMORY_DEBUG(real_word_graph_edge)
};

class real_count_word_graph_edge : public real_word_graph_edge {
  unsigned count_;
 public:
  real_count_word_graph_edge() : count_(0) {}
  virtual ~real_count_word_graph_edge() {};
  unsigned int count() const {
    return count_;
  };
  void count(unsigned int c) {
    count_ = c;
  }
  MEMORY_DEBUG(real_count_word_graph_edge)
};

class artificial_word_graph_edge : public word_graph_edge {
  std::string sequence_;
 public:
  artificial_word_graph_edge() {}
  virtual ~artificial_word_graph_edge() {};
  std::string sequence(CharacterProducer & cp) {
    return sequence_;
  }
  bool isreal() const { return false; }
  unsigned length() const {
    return sequence_.length();
  }
  void sequence(std::string const & s) {
    sequence_ = s;
  }
  unsigned int count() const {
    return 0;
  };
  MEMORY_DEBUG(artificial_word_graph_edge)
};

class restart_word_graph_edge : public word_graph_edge {
  static char eos_char_;
 public:
  restart_word_graph_edge() {}
  virtual ~restart_word_graph_edge() {};
  std::string sequence(CharacterProducer & cp) {
    std::string retval(1,eos_char_);
    retval += to()->sequence(cp);
    return retval;
  }
  unsigned length() const {
    return to()->length()+1;
  }
  static void eos_char(char ec) {
    eos_char_ = ec;
  }
  bool isreal() const { return false; }
  unsigned int count() const {
    return 0;
  };
  MEMORY_DEBUG(restart_word_graph_edge)
};

class sim_word_graph_edge : public word_graph_edge {
 public:
  sim_word_graph_edge() {}
  virtual ~sim_word_graph_edge() {};
  std::string sequence(CharacterProducer & cp) {
    return "";
  }
  unsigned length() const {
    return 0;
  }
  bool isreal() const { return false; }
  unsigned int count() const {
    return 0;
  };
  MEMORY_DEBUG(sim_word_graph_edge)
};

class word_graph : public intlabelgraph<word_graph_node,word_graph_edge> {
public:
  word_graph(long maxlabelhint=10000) : 
    intlabelgraph<word_graph_node,word_graph_edge>(maxlabelhint) {}
  void read(std::string const &, int=-1, int=0, int=0, bool=false);
  void dump(ostream &, CharacterProducer &, bool=false);
  void print_stats();
  long unsigned int balance_nodes(char, bool);
  long unsigned int find_joiners(bool, bool);
  void writeseq(ostream &,CharacterProducer &,char,bool);
  void writetrivialpaths(ostream &,CharacterProducer &,char,bool);
  void annotateseq(ostream &,CharacterProducer &,FastaFile<Lazy_Header_SI> &,
		   char,int,int,bool);
  bool check_out_edges(ostream &os, CharacterProducer & cp, bool);
  long unsigned int remove_trivial_nodes(CharacterProducer & cp, bool verbose);
  void sort_nodes(CharacterProducer & cp);
  long unsigned int peel_edges(FastaFile<Lazy_Header_SI> & cp, int, char);
  MEMORY_DEBUG(word_graph)
};

#endif
