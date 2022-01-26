
#include <iostream>
#include <sstream>
#include <list>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <assert.h>
#include "char_io.h"
#include "types.h"
#include <string.h>

#include "graph.h"
#include "trans_prob.h"

#include "keyword_tree.h"
#include "keyword_tree.t"

class word_graph_edge;
class word_graph_node : public labelnode<word_graph_edge> {
  std::string sequence_;
public:
  word_graph_node() {};
  std::string const & sequence() const {
    return sequence_;
  }
  void sequence(std::string const & s) {
    sequence_ = s;
  }
  void dump(ostream &, unsigned int);
};
class word_graph_edge : public labeledge<word_graph_node> {
  FILE_POSITION_TYPE seq_start_;
  FILE_POSITION_TYPE seq_end_;
  std::string sequence_;
public:
  std::string const & sequence() const {
    return sequence_;
  }
  void sequence(std::string const & s) {
    sequence_ = s;
  }
  FILE_POSITION_TYPE seq_start() const {
    return seq_start_;
  }
  void seq_start(FILE_POSITION_TYPE const & p) {
    seq_start_ = p;
  }
  FILE_POSITION_TYPE seq_end() const {
    return seq_end_;
  }
  void seq_end(FILE_POSITION_TYPE const & p) {
    seq_end_ = p;
  }
};

class word_graph : public labelgraph<word_graph_node,word_graph_edge> {
public:
  void read(std::string const &, std::string const &,unsigned int);
  void dump(ostream &, unsigned int);
  unsigned long int uniquify_nodes(bool);
  long unsigned int uniquify_edges(bool);
  void print_stats();
  long unsigned int remove_trivial_nodes(bool);
  long unsigned int remove_eos(char,unsigned int,bool);
  void rename_nodes(bool);
  long unsigned int balance_nodes(char,bool allownew, bool);
  void join_components(char eos, bool allownew, bool verbose);
  long unsigned int find_joiners(unsigned int,bool,bool);
  long unsigned int find_jumpers(unsigned int,bool,bool);
  void writeseq(ostream &,bool);
  void write(ostream &) const;
};

typedef std::list<word_graph_node*> nodelist;
typedef std::list<word_graph_edge*> edgelist;
typedef std::list<word_graph_node*>::iterator nodelist_it;
typedef std::list<word_graph_edge*>::iterator edgelist_it;
typedef std::list<word_graph_node*>::const_iterator nodelist_cit;
typedef std::list<word_graph_edge*>::const_iterator edgelist_cit;

void
word_graph::read(std::string const & gf, 
		 std::string const & sf, 
		 unsigned int mersize) {
  std::string label;
  long unsigned int count;
  long unsigned int nodefrom;
  long unsigned int nodeto;
  FILE_POSITION_TYPE seqst;
  FILE_POSITION_TYPE seqed;
  std::stringstream gfss(gf.c_str());
  std::stringstream sfss(sf.c_str());
  std::string gfs="";
  std::string sfs="";
  while (gfss.good()) {
    gfss >> gfs;
    sfss >> sfs;
    if (gfs == "" || sfs == "") break;
    MapFileChars seq(sfs.c_str());
    ifstream gis(gfs.c_str());
    long unsigned int label_base = new_label();
    std::string line;
    while (getline(gis,line)) {
      istringstream ss(line);
      count = 0;
      ss >> label >> nodefrom >> nodeto >> seqst >> seqed >> count;
      if (label != "E") {
	continue;
      }
      nodefrom += label_base;
      nodeto += label_base;
      word_graph_node *f;
      if (!(f=find(nodefrom))) {
	f = new word_graph_node;
	f->name(nodefrom);
	if (seqst >= mersize-1) {
	  char *buffer = new char[mersize];
	  seq.pos(seqst-(mersize-1));
	  long unsigned int i=0;
	  while (i < mersize-1) {
	    buffer[i] = seq.getch();
	    i++;
	  }
	  buffer[i]='\0';
	  f->sequence(buffer);
	  delete [] buffer;
	}
	f->unset_mark();
	new_node(f);
      }
      word_graph_node *t;
      if (!(t=find(nodeto))) {
	t = new word_graph_node;
	t->name(nodeto);
	if (seqed >= mersize-1) {
	  char *buffer = new char[mersize];
	  seq.pos(seqed-(mersize-1));
	  long unsigned int i=0;
	  while (i < mersize-1) {
	    buffer[i] = seq.getch();
	    i++;
	  }
	  buffer[i]='\0';
	  t->sequence(buffer);
	  delete [] buffer;
	}
	t->unset_mark();
	new_node(t);
      }
      word_graph_edge *e;
      e = new word_graph_edge;
      e->from(f);
      e->to(t);
      e->seq_start(seqst);
      e->seq_end(seqed);
      char *buffer = new char[seqed-seqst+1];
      seq.pos(seqst);
      long unsigned int i=0;
      while (i < seqed-seqst) {
	buffer[i] = seq.getch();
	i++;
      }
      buffer[i]='\0';
      e->sequence(buffer);
      delete [] buffer;
      e->unset_mark();
      new_edge(e);
    }
    gis.close();
  }

  checkpoint;

  nodelist_it nit=nodes().begin();
  while (nit!= nodes().end()) {
    if ((*nit)->sequence() == "") {
      checkpoint;
      std::string nodeseq="";
      edgelist_it eit = (*nit)->in().begin();
      while (nodeseq.length() < (mersize-1)) {
	long unsigned int edgelen=(*eit)->sequence().length();
	long unsigned int needed=((mersize-1)-nodeseq.length());
	if (needed>edgelen) {
	  needed=edgelen;
	}
	nodeseq = (*eit)->sequence().substr(edgelen-needed,needed) + nodeseq;
	eit = (*eit)->from()->in().begin();
      }
      (*nit)->sequence(nodeseq);
    }
    ++nit;
  }

}

void
word_graph_node::dump(ostream & os, unsigned int mersize) {
  edgelist_cit eit=in().begin();
  while (eit != in().end()) {
    os << "    <- " << ((*eit)->mark()?"(T)":"(F)") << " [" << (*eit)->from()->name() << "] " << (*eit)->sequence() << endl;
    (*eit)->set_mark();
    ++eit;
  }
  os << sequence() << " [" << name() << "] " 
     << (mark()?"(T)":"(F)") << endl;
  eit=out().begin();
  while (eit != out().end()) {
    os << "    -> "<< ((*eit)->mark()?"(T)":"(F)") << " [" << (*eit)->to()->name() << "] " << (*eit)->sequence() << endl;
    (*eit)->set_mark();
    ++eit;
  }
  os << endl;
}

void
word_graph::dump(ostream & os, unsigned int mersize) {
  os << "Number of nodes: " << nnode() << endl;
  os << "Number of edges: " << nedge() << endl;
  clear_edge_marks();
  nodelist_cit nit=nodes().begin();
  while (nit != nodes().end()) {
    (*nit)->dump(os,mersize);
    ++nit;
  }
  edgelist_cit eit=edges().begin();
  while (eit != edges().end()) {
    if (!(*eit)->mark()) {
      os << "Unattached edge: " << 
	" [" << (*eit)->from()->name() << "] " 
	 << ((*eit)->from()->mark()?"(T)":"(F)") << " -- "<< ((*eit)->mark()?"(T)":"(F)") << "->" <<  ((*eit)->to()->mark()?"(T)":"(F)") << " [" << (*eit)->to()->name() << "] " << (*eit)->sequence() << endl;
      }
      ++eit;
  }
  clear_edge_marks();
}

unsigned long int
word_graph::uniquify_nodes(bool verbose) {
  long unsigned int count=0;
  typedef std::multimap<std::string,nodelist::iterator> nodeseqmap;
  nodeseqmap nsm;
  nodelist_it nit=nodes().begin();
  while (nit != nodes().end()) {
    nsm.insert(make_pair((*nit)->sequence(),nit));
    ++nit;
  }

  clear_node_marks();

  nodeseqmap::iterator sit,sit0;
  long unsigned int dupcount=0;
  sit0=sit=nsm.begin();
  ++sit;
  while (sit!=nsm.end()) {
    if (sit->first != sit0->first) {
      if (dupcount>0) {
	if (verbose) {
	  cout << "Duplicate node: " << sit0->first 
	       << ", " << dupcount+1 << " copies" << endl;
	}
	dupcount=0;
      }
      sit0=sit;
      ++sit;
      continue;
    } 
    // OK, remove the sit node but re-attach all the edges
    dupcount++;
    move_edges(*(sit->second),*(sit0->second));
    (*(sit->second))->set_mark();
    count++;
    ++sit;
  }
  delete_marked_nodes();
  return count;
}

unsigned long int
word_graph::uniquify_edges(bool verbose) {
  long unsigned int count=0;
  typedef std::multimap<std::string,edgelist::iterator> edgeftmap;
  edgeftmap em;
  edgelist_it eit=edges().begin();
  while (eit != edges().end()) {
    ostrstream ss;
    ss << ":" << (*eit)->from()
       << ":" << (*eit)->to() 
       << ":" << (*eit)->sequence()
       << ":" << ends;
    std::string key(ss.str());
    em.insert(make_pair(key,eit));
    ++eit;
  }

  edgeftmap::iterator ftit,ftit0;
  long unsigned int dupcount=0;
  ftit0=ftit=em.begin();
  ++ftit;
  while (ftit!=em.end()) {
    if (ftit->first != ftit0->first) {
      if (dupcount>0) {
	if (verbose) {
	  cout << "Duplicate edge: " <<  "[" << (*(ftit0->second))->from()->name() << "]"
	       << " -> [" << (*(ftit0->second))->to()->name() << "], " << dupcount+1 << " copies" << endl;
	}
	dupcount=0;
      }
      ftit0=ftit;
      ++ftit;
      continue;
    }
    dupcount++;
    count++;
    erase(ftit->second);
    ++ftit;
  }
  return count;
}

long unsigned int 
word_graph::remove_trivial_nodes(bool verbose) {
  long unsigned int count=0;

  clear_node_marks();
  clear_edge_marks();

  count=0;
  nodelist_it nit;
  nit=nodes().begin();
  while (nit!=nodes().end()) {
    if ((*nit)->nin() == 1 && (*nit)->nout() == 1) {
      word_graph_node *n=(*nit);
      word_graph_edge *ein=n->in().front();
      word_graph_edge *eout=n->out().front();
      word_graph_node *nin=ein->from();
      word_graph_node *nout=eout->to();
      if (ein != eout) { //self loop!
	word_graph_edge *enew = new word_graph_edge;
	enew->from(nin);
	enew->to(nout);
	std::string seq = ein->sequence() + eout->sequence();
	enew->sequence(seq);
	enew->seq_start(ein->seq_start());
	enew->seq_end(eout->seq_end());
	enew->unset_mark();
	new_edge(enew);
	count++;
	nout->remove_in(eout);
	nin->remove_out(ein);
	n->remove_out(eout);
	n->remove_in(ein);
	eout->set_mark();
	ein->set_mark();
	n->set_mark();
      }
    } else if ((*nit)->nin() == 0 && (*nit)->nout() == 0) {
      // checkpoint;
      // cerr << "0 -> " << (*nit)->name() << " <- 0" << endl;
      (*nit)->set_mark();
      count++;
    }
    ++nit;
  }

  delete_marked_nodes();
  delete_marked_edges();
  checkpoint;
  return count;
}

void
word_graph::rename_nodes(bool verbose) {
  nodelist_it nit;
  nit = nodes().begin();
  long unsigned int index=1;
  while (nit!=nodes().end()) {
    (*nit)->name(index);
    ++index;
    ++nit;
  }
}

long unsigned int 
word_graph::remove_eos(char eos, unsigned int mersize, bool verbose) {
  if (verbose) { timestamp("graph::remove_eos"); }
  long unsigned int count=0;
  edgelist_it eit;
  std::string::size_type pos;
  eit = edges().begin();
  while (eit != edges().end()) {
    if (((pos=((*eit)->sequence().find(eos)))==std::string::npos)||
	((*eit)->sequence().length() == 1 && (*eit)->sequence()[0] == eos)) {
      ++eit;
      continue;
    }
    /* if (verbose) {
       cout << "EOS edge: " <<  "[" << (*eit)->from->name << "]"
       << " -> [" << (*eit)->to->name << "], " 
       // << (*eit)->sequence << " EOS position " << pos 
       << endl;
       }*/
    
    // Create a new node to represent the last valid (mersize-1)-mer.
    // Create a new node to represent the next (mersize-1)-mer (has
    // eos as last char)
    std::string nodeedgeseq = (*eit)->from()->sequence() + (*eit)->sequence();
    
    /*	if (verbose) 
	cerr << "N E: " << (*eit)->from->sequence << " " << (*eit)->sequence << endl; */
    word_graph_node *nvalid = (*eit)->from();
    if (pos > 0) {
      nvalid = new word_graph_node;
      nvalid->sequence(nodeedgeseq.substr(pos,mersize-1));
      nvalid->name(new_label());
      nvalid->unset_mark();
      new_node(nvalid);
      /* if (verbose) 
	 cerr << "Created node [" << nvalid->name << "] " << nvalid->sequence << endl; */
    }
    
    word_graph_node *ninvalid = (*eit)->to();
    if (pos < (*eit)->sequence().length()-1) {
      ninvalid = new word_graph_node;
      ninvalid->sequence(nodeedgeseq.substr(pos+1,mersize-1));
      ninvalid->name(new_label());
      ninvalid->unset_mark();
      new_node(ninvalid);
      /* if (verbose) 
	 cerr << "Created node [" << ninvalid->name << "] " << ninvalid->sequence << endl; */
    }
    
    if ((*eit)->from() != nvalid) {
      word_graph_edge *e1 = new word_graph_edge;
      e1->from((*eit)->from());
      e1->to(nvalid);
      e1->sequence((*eit)->sequence().substr(0,pos));
      e1->seq_start((*eit)->seq_start());
      e1->seq_end((*eit)->seq_start()+pos);
      e1->unset_mark();
      /* if (verbose)
	 cerr << "Created edge [" << e1->from->name << "] -> [" << e1->to->name << "] " << e1->sequence << endl; */
      new_edge(e1);
    }
    
    word_graph_edge *e2 = new word_graph_edge;
    e2->from(nvalid);
    e2->to(ninvalid);
    e2->sequence(std::string()+eos);
    e2->seq_start((*eit)->seq_start()+pos);
    e2->seq_end((*eit)->seq_start()+pos+1);
    e2->unset_mark();
    /* if (verbose) 
       cerr << "Created edge [" << e2->from->name << "] -> [" << e2->to->name << "] " << e2->sequence << endl; */
    new_edge(e2);
    
    if ((*eit)->to() != ninvalid) {
      word_graph_edge *e3 = new word_graph_edge;
      e3->from(ninvalid);
      e3->to((*eit)->to());
      e3->sequence((*eit)->sequence().substr(pos+1));
      e3->seq_start((*eit)->seq_start()+pos+1);
      e3->seq_end((*eit)->seq_end());
      e3->unset_mark();
      /* if (verbose) 
	 cerr << "Created edge [" << e3->from->name << "] -> [" << e3->to->name << "] " << e3->sequence << endl; */
      new_edge(e3);
    }
    eit = erase(eit);
    count++;
  }
  
  clear_edge_marks();

  nodelist_it nit;
  nit = nodes().begin();
  while (nit != nodes().end()) {
    if ((*nit)->sequence()[mersize-2] != eos) {
      ++nit;
      continue;
    }
    // This node represents the first invalid (mersize-1)-mer. We must
    // find all positions after this node and mark then with a new
    // pair of nodes, as above.
    /* if (verbose) 
       cerr << "[" << (*nit)->name << "] > " << (*nit)->nout << ": " << (*nit)->sequence << endl; */
    typedef std::pair<word_graph_node*,unsigned int> stackelt;
    std::list<stackelt> stack;
    stack.push_front(stackelt(*nit,0));
    while (!stack.empty()) {
      stackelt se = stack.front();
      stack.pop_front();
      // if (verbose) 
      // cerr << "Pop (" << se.second << ") [" << se.first->name << "] > " << se.first->nout << ": " << se.first->sequence << endl;
      eit = se.first->out().begin();
      while (eit != se.first->out().end()) {
	// if (verbose)
	// cerr << "Considering edge [" << (*eit)->from->name << "] -> [" << (*eit)->to->name << "] " << (*eit)->sequence << endl;	
	if (se.second+(*eit)->sequence().length()>=mersize-1) {
	  // The last invalid (mersize-1)-mer is on this edge. 
	  if ((*eit)->sequence().length()>1 && !(*eit)->mark()) {
	    word_graph_node * n = se.first;
	    word_graph_edge * e = (*eit);
	    long unsigned int pos = (mersize-1) - se.second - 1;
	    std::string neseq = n->sequence() + e->sequence();

	    // if (verbose)
	    // cerr << "N E: " << n->sequence << " " << e->sequence << endl;
	    
	    word_graph_node *ninvalid = e->from();
	    if (pos > 0) {
	      ninvalid = new word_graph_node;
	      ninvalid->sequence(neseq.substr(pos,mersize-1));
	      ninvalid->name(new_label());
	      ninvalid->unset_mark();
	      /* if (verbose) 
		 cerr << "Created node [" << ninvalid->name << "] " 
		 << ninvalid->sequence << endl; */
	      new_node(ninvalid);
	    }

	    word_graph_node *nvalid = e->to();
	    if (pos < e->sequence().length()-1) {
	      nvalid = new word_graph_node;
	      nvalid->sequence(neseq.substr(pos+1,mersize-1));
	      nvalid->name(new_label());
	      nvalid->unset_mark();
	       /* if (verbose) 
		 cerr << "Created node [" << nvalid->name << "] " 
		 << nvalid->sequence << endl; */
	      new_node(nvalid);
	    }

	    if (e->from() != ninvalid) {
	      word_graph_edge *e1 = new word_graph_edge;
	      e1->from(e->from());
	      e1->to(ninvalid);
	      e1->sequence(e->sequence().substr(0,pos));
	      e1->seq_start(e->seq_start());
	      e1->seq_end(e->seq_start()+pos);
	      e1->unset_mark();
	      // if (verbose)
	      // cerr << "Created edge [" << e1->from->name << "] -> [" << e1->to->name << "] " << e1->sequence << endl;
	      new_edge(e1);
	    }

	    word_graph_edge *e2 = new word_graph_edge;
	    e2->from(ninvalid);
	    e2->to(nvalid);
	    e2->sequence(e->sequence().substr(pos,1));
	    e2->seq_start(e->seq_start()+pos);
	    e2->seq_end(e->seq_start()+pos+1);
	    e2->unset_mark();
	    // if (verbose) 
	    // cerr << "Created edge [" << e2->from->name << "] -> [" << e2->to->name << "] " << e2->sequence << endl;
	    new_edge(e2);
	    
	    if (e->to() != nvalid) {
	      word_graph_edge *e3 = new word_graph_edge;
	      e3->from(nvalid);
	      e3->to(e->to());
	      e3->sequence(e->sequence().substr(pos+1));
	      e3->seq_start(e->seq_start()+pos+1);
	      e3->seq_end(e->seq_end());
	      e3->unset_mark();
	      // if (verbose) 
	      // cerr << "Created edge [" << e3->from->name << "] -> [" << e3->to->name << "] " << e3->sequence << endl;
	      new_edge(e3);
	    }
	    e->set_mark();
	    count++;
	  }
	} else {
	  // if (verbose) 
	    /* cerr << "Push ( " << se.second+(*eit)->sequence.length() 
	       << ") [" << (*eit)->to->name << "] > " << (*eit)->to->nout 
	       << ": " << (*eit)->to->sequence << endl; */
	  stack.push_front(stackelt((*eit)->to(),
				    se.second+(*eit)->sequence().length()));
	}
	++eit;
      }
    }
    ++nit;
  }
  
  // checkpoint;
  
  eit = edges().begin();
  while (eit!=edges().end()) {
    if ((*eit)->mark()) {
      eit = erase(eit);
    } else {
      ++eit;
    }
  }

  checkpoint;
   
  // If we do not distinguish mers that occur at the start and end of
  // fasta entires, we just remove all nodes whose sequence contains
  // an eos character. If we do distinguish mers at the start and end
  // of fasta entries, we leave nodes that start or end in a eos
  // character, but ensure that they have no in or out edges as
  // appropriate.

  bool distinguish_eos = /*false; */ true;

  nodelist del;

  if (!distinguish_eos) {
    nit = nodes().begin();
    while (nit!=nodes().end()) {
      if ((*nit)->sequence().find(eos)!=std::string::npos) {
	del.push_back(*nit);
      }
      ++nit;
    }
  } else {
    nit = nodes().begin();
    while (nit!=nodes().end()) {
      std::string::size_type pos=std::string::npos;
      pos=(*nit)->sequence().find(eos);
      if (pos != std::string::npos && 
	  pos != 0 && 
	  pos != mersize-2) {
	del.push_back(*nit);
      }
      ++nit;
    }
  }
  
  remove(del);

  if (distinguish_eos) {
  
    clear_edge_marks();

    nit = nodes().begin();
    while (nit!=nodes().end()) {
      std::string::size_type pos=std::string::npos;
      pos=(*nit)->sequence().find(eos);
      if (pos != std::string::npos) {
	if (pos == 0) {
	  mark_in_edges(*nit);
	} else if (pos == mersize-2) {
	  mark_out_edges(*nit);
	} else {
	  assert(0);
	}
      }
      ++nit;
    }

    remove_marked_edges();

  }

  

  checkpoint;

  return count;
}

/* 
void
word_graph::join_components(char eos, bool allownew, bool verbose) {
  typedef std::multimap<long unsigned int,long unsigned int> nodecompmap;
  nodecompmap n2c;
  typedef std::multimap<long unsigned int,word_graph_node*> compnodemap;
  compnodemap c2n;
  long unsigned int component_number = 0;
  long unsigned int elts = 0;
  std::list<word_graph_node*> stack;
  clear_node_marks();

  nodelist_it nit = nodes().begin();
  while (nit != nodes().end()) {
    if (!(*nit)->mark()) {
      elts = 0;
      stack.clear();
      word_graph_node *n = (*nit);
      stack.push_back(n);
      n->set_mark();
      while (!stack.empty()) {
	n = stack.front();
	stack.pop_front();
	n2c.insert(make_pair(n->name(),component_number));
	c2n.insert(make_pair(component_number,n));
	elts++;
	edgelist_it eit = n->out().begin();
	while (eit != n->out().end()) {
	  if (!(*eit)->to()->mark()) {
	    stack.push_back((*eit)->to());
	    (*eit)->to()->set_mark();
	  }
	  ++eit;
	}
	eit = n->in().begin();
	while (eit != n->in().end()) {
	  if (!(*eit)->from()->mark()) {
	    stack.push_back((*eit)->from());
	    (*eit)->from()->set_mark();	    
	  }
	  ++eit;
	}
      }
      // checkpoint;
      // cerr << "Component " << component_number << " has "
      // << elts << " elements.\n";
      component_number++;
      
    }
    ++nit;
  }
  
  // Now we have components, we need to find unbalanced nodes in each. 
  std::vector<word_graph_node*> comp2s(component_number,0);
  std::vector<word_graph_node*> comp2d(component_number,0);
  compnodemap::iterator it=c2n.begin();
  while (it != c2n.end()) {
    if (it->second->nin() > it->second->nout() &&
	comp2s[it->first] == 0) {
      comp2s[it->first] = it->second;
    } else if (it->second->nin() < it->second->nout() &&
	comp2d[it->first] == 0) {
      comp2d[it->first] = it->second;
    }
    ++it;
  }
  
  // check that all components have at least one in and one out...
  std::list<long unsigned int> snd,sonly,donly;
  long unsigned int nosnd=0;
  long unsigned int i=0;
  for (i=0;i<component_number;i++) {
    if (comp2s[i]==0 && comp2d[i]!=0) {
      donly.push_back(i);
    } else 
    if (comp2s[i]!=0 && comp2d[i]==0) {
      sonly.push_back(i);
    } else 
    if (comp2s[i]!=0 && comp2d[i]!=0) {
      snd.push_back(i);
    } else 
    if (comp2s[i]==0 && comp2d[i]==0) {
      nosnd++;
    } 
  }

  checkpoint;
  cerr << "Components with only supply nodes: " << sonly.size() << endl;
  cerr << "Components with only demand nodes: " << donly.size() << endl;
  cerr << "Components with supply and demand nodes: " << snd.size() << endl;
  cerr << "Components no supply and demand nodes: " << nosnd << endl;
  
  while (!sonly.empty() && !snd.empty()) {
    long unsigned int c1 = sonly.front();
    sonly.pop_front();
    long unsigned int c2 = snd.front();
    snd.pop_front();
    word_graph_edge * e = new word_graph_edge;
    e->from(comp2s[c1]);
    e->to(comp2d[c2]);
    std::string seq;
    if (!allownew) {
      seq += eos;
    }
    seq += e->to()->sequence();
    e->sequence(seq);
    e->unset_mark();
    new_edge(e);
    sonly.push_back(c2);
  }

  while (!donly.empty() && !snd.empty()) {
    long unsigned int c1 = snd.front();
    sonly.pop_front();
    long unsigned int c2 = donly.front(); 
    snd.pop_front();
    word_graph_edge * e = new word_graph_edge;
    e->from(comp2s[c1]);
    e->to(comp2d[c2]);
    std::string seq;
    if (!allownew) {
      seq += eos;
    }
    seq += e->to()->sequence();
    e->sequence(seq);
    e->unset_mark();
    new_edge(e);
    donly.push_back(c1);
  }
  
  while (!sonly.empty() && !donly.empty()) {
    long unsigned int c1 = sonly.front();
    sonly.pop_front();
    long unsigned int c2 = donly.front();
    snd.pop_front();
    word_graph_edge * e = new word_graph_edge;
    e->from(comp2s[c1]);
    e->to(comp2d[c2]);
    std::string seq;
    if (!allownew) {
      seq += eos;
    }
    seq += e->to()->sequence();
    e->sequence(seq);
    e->unset_mark();
    new_edge(e);
  }

  std::list<long unsigned int>::iterator lit,lit_1;
  
  lit = snd.begin();
  if (lit != snd.end()) {
    lit_1 = lit;
    lit++;
    while (lit != snd.end()) {
      word_graph_edge * e = new word_graph_edge;
      e->from(comp2s[*lit_1]);
      e->to(comp2d[*lit]);
      std::string seq;
      if (!allownew) {
	seq += eos;
      }
      seq += e->to()->sequence();
      e->sequence(seq);
      e->unset_mark();
      new_edge(e);
      lit_1 = lit;
      ++lit;
    }
  }
  
}

*/ 

void word_graph::print_stats() {

  cerr << "CSBH-graph nodes: " << nnode() << endl;
  cerr << "CSBH-graph edges: " << nedge() << endl;

  long unsigned int l=0;
  long unsigned int edge_length=0;
  long unsigned int suppressed_nodes=0;
  edgelist_it eit=edges().begin();
  while (eit != edges().end()) {
    l = (*eit)->sequence().length();
    edge_length += l;
    suppressed_nodes += l - 1;
    ++eit;
  }
  cerr << "CSBH-graph total edge length: " << edge_length << endl;
  cerr << "Suppressed SBH-graph nodes: " << suppressed_nodes << endl;

  components<word_graph_node> cmps;
  dfs(cmps);

  cerr << "CSBH-graph components: " << cmps.size() << endl;
  
  long unsigned int comp_bal_plus;
  long unsigned int comp_bal_minus;
  long unsigned int comp_n_plus;
  long unsigned int comp_n_minus;
  long unsigned int nbalanced=0;
  long unsigned int nunbalanced=0;
  long unsigned int total_bal_plus=0;
  long unsigned int total_bal_minus=0;  
  long unsigned int total_n_plus=0;
  long unsigned int total_n_minus=0;
  
  for (long unsigned int i=0;i<cmps.size();i++) {
    comp_bal_plus = 0;
    comp_bal_minus = 0;
    comp_n_plus = 0;
    comp_n_minus = 0;
    nodelist const & cmp = cmps.nodes(i);
    nodelist_cit cit=cmp.begin();
    while (cit != cmp.end()) {
      if ((*cit)->nin() > (*cit)->nout()) {
	comp_n_plus++;
	comp_bal_plus += (*cit)->nin() - (*cit)->nout();
      } else if ((*cit)->nin() < (*cit)->nout()) {
	comp_n_minus++;
	comp_bal_minus += (*cit)->nout() - (*cit)->nin();
      }
      ++cit;
    }
    if (comp_n_plus > 0 || comp_n_minus > 0) {
      nunbalanced++;
    } else {
      nbalanced++;
    }
    total_bal_plus += comp_bal_plus;
    total_bal_minus += comp_bal_minus;
    total_n_plus += comp_n_plus;
    total_n_minus += comp_n_minus;
  }

  cerr << "CSBH-graph unbalanced components: " << nunbalanced << endl;
  cerr << "CSBH-graph balanced components: " << nbalanced << endl;
  
  cerr << "Degree surplus nodes: " << total_n_plus << endl;
  cerr << "Degree deficit nodes: " << total_n_minus << endl;
  cerr << "Total degree surplus: " << total_bal_plus << endl;
  cerr << "Total degree deficit: " << total_bal_minus << endl;
} 


long unsigned int
word_graph::balance_nodes(char eos, bool allownew, bool verbose) {

  components<word_graph_node> cmps;
  dfs(cmps);

  for (long unsigned int i0=0;i0<cmps.size();i0++) {
    
    long unsigned int bplus=0;
    nodelist_cit nit = cmps.nodes(i0).begin();
    while (nit != cmps.nodes(i0).end()) {
      if ((*nit)->nin() > (*nit)->nout()) {
	bplus += (*nit)->nin() - (*nit)->nout();
      }
      ++nit;
    }

    long unsigned int edges_added=0;
    
    nodelist_cit nitneg;
    nodelist_cit nitpos;
    nitneg = cmps.nodes(i0).begin();
    nitpos = cmps.nodes(i0).begin();
    while (nitneg != cmps.nodes(i0).end() && (*nitneg)->nout() >= (*nitneg)->nin()) ++nitneg;
    while (nitpos != cmps.nodes(i0).end() && (*nitpos)->nout() <= (*nitpos)->nin()) ++nitpos;
    
    int negrem=0;
    if (nitneg != cmps.nodes(i0).end()) {
      negrem = (*nitneg)->nin() - (*nitneg)->nout();
    }
    int posrem=0;
    if (nitpos != cmps.nodes(i0).end()) {
      posrem = (*nitpos)->nout() - (*nitpos)->nin();
    }

    while (edges_added < bplus - 1 && 
	   nitneg != cmps.nodes(i0).end() && 
	   nitpos != cmps.nodes(i0).end()) {
      int nnew = negrem;
      if (posrem < negrem) {
	nnew = posrem;
      }
      // cerr << "[" << (*nitneg)->name << "] " << negrem << " [" << (*nitpos)->name << "] " << posrem << " " << nnew << endl;
      for (int i=0;i<nnew;i++) {
	edges_added++;
	if (edges_added >= bplus) break;
	word_graph_edge * e = new word_graph_edge;
	e->from(*nitneg);
	e->to(*nitpos);
	std::string seq;
	if (!allownew) {
	  seq += eos;
	}
	seq += e->to()->sequence();
	e->sequence(seq);
	e->unset_mark();
	new_edge(e);
      }
      negrem -= nnew;
      if (negrem == 0) {
	while (nitneg != cmps.nodes(i0).end() && (*nitneg)->nout() >= (*nitneg)->nin()) 
	  ++nitneg;
	if (nitneg != cmps.nodes(i0).end()) {
	  negrem = (*nitneg)->nin() - (*nitneg)->nout();
	}
      }
      posrem -= nnew;
      if (posrem == 0) {
	while (nitpos != cmps.nodes(i0).end() && (*nitpos)->nout() <= (*nitpos)->nin()) 
	  ++nitpos;
	if (nitpos != cmps.nodes(i0).end()) {
	  posrem = (*nitpos)->nout() - (*nitpos)->nin();
	}
      }
    }
    if (edges_added < bplus - 1 && bplus > 1) {
      checkpoint;
      assert(0);
    }
  }
  return 0;
}

void
word_graph::write(ostream& os) const {
  edgelist_cit eit;
  eit = edges().begin();
  while (eit != edges().end()) {
    os << "E " 
       << (*eit)->from()->name() << " "
       << (*eit)->to()->name() << " "
       << (*eit)->seq_start() << " "
       << (*eit)->seq_end() << endl;
    ++eit;
  }
}

void 
word_graph::writeseq(ostream &os, bool verbose) {

  clear_edge_marks();
  clear_node_marks();

  long unsigned int total_tour_length=0;
  long unsigned int subtourlen = 0;
  checkpoint;
  nodelist_cit nit;

  components<word_graph_node> cmps;
  dfs(cmps);

  checkpoint;

  for (long unsigned int i=0;i<cmps.size();i++) {
    nit = cmps.nodes(i).begin();
    while (nit != cmps.nodes(i).end()) {
      if ((*nit)->nin() < (*nit)->nout()) {
	break;
      } 
      ++nit;
    }
    if (nit == cmps.nodes(i).end()) {
      nit = cmps.nodes(i).begin();
    }
    // checkpoint;
    edgelist tour;
    edgelist_it tip=tour.end();
    
    word_graph_node *n;
    word_graph_edge *e;
    
    edgelist_it eit;

    n = (*nit);
    n->set_mark();

    while (1) {

      // checkpoint;
      // cerr << "Starting sub-tour at node " << "[" 
      // << n->name() << "]" << endl;
      subtourlen = 0;
      edgelist_it eit = n->out().begin();
      while (eit != n->out().end()) {
	if (!(*eit)->mark()) {
	  // Follow edge
	  e = (*eit);
	  tour.insert(tip,e);
	  e->set_mark();
	  n = e->to();
	  n->set_mark();
	  eit = n->out().begin();
	  // cerr << "[" << n->name() << "]" << endl;
	  subtourlen++;
	} else {
	  ++eit;
	}
      }

      // checkpoint;
      // cerr << "Ending sub-tour at node " << "[" 
      // << n->name() << "]: Length: " << subtourlen << endl;
      
      assert(tip == tour.end() || n == (*tip)->from());

      tip = tour.begin();
      while (tip != tour.end()) {
	n = (*tip)->from();
	bool edgeout = false;
	eit = n->out().begin();
	while (eit != n->out().end()) {
	  if (!(*eit)->mark()) {
	    edgeout = true;
	    break;
	  }
	  ++eit; 
	}
	if (edgeout) {
	  break;
	}
	++tip;
      }
      if (tip == tour.end()) {
	break;
      }
    }

    // checkpoint;
    
    // cerr << "Tour edges: " << tour.size() << endl;
    total_tour_length += tour.size();

    if (tour.size()>0) {
      eit = tour.begin();
      os << (*eit)->from()->sequence();
      while (eit != tour.end()) {
	os << (*eit)->sequence();
	++eit;
      }
      os << '|';
    }
  }  
  cerr << "Total number of edges: " << edges().size() << endl;
  cerr << "Total tour length: " << total_tour_length << endl;

  edgelist_it eit = edges().begin();
  while (eit != edges().end()) {
    if (!(*eit)->mark()) {
      checkpoint;
      cerr << "Unmarked edge: " 
	   << "[" << (*eit)->from()->name() << "]" << " --> " 
	   << "[" << (*eit)->to()->name() << "]" << endl;
	
    }
    ++eit;
  }
}

struct stackelt{
  word_graph_node * n;
  int dist;
  edgelist el;
  stackelt(word_graph_node *n0, int d0, edgelist const & el0) 
    : n(n0), dist(d0), el(el0) {};
  stackelt(word_graph_node *n0=0, int d0=0) 
    : n(n0), dist(d0){};
};

long unsigned int 
word_graph::find_joiners(unsigned int mersize, bool optimize, bool verbose) {
  
  trans_prob_graph<word_graph_node*,edgelist> tpg;

  nodelist_it nit;
  nit = nodes().begin();
  while (nit != nodes().end()) {
    // cerr << "Consider " << (*nit)->nin << " > ["
    // << (*nit)->name << "] < " << (*nit)->nout 
    // << endl;
    if ((*nit)->nin() > (*nit)->nout()) {
      // Node is unbalanced, more in than out..
      word_graph_node *n = (*nit);

      std::list<stackelt> stack;
      stack.push_front(stackelt(n,0));
      while (!stack.empty()) {
	stackelt se = stack.front();
	stack.pop_front();
	if (se.n->nout() > se.n->nin()) {
	  trans_prob_node<word_graph_node*,edgelist> *f;
	  if (!(f=tpg.find(n->name()))) {
	    f = new trans_prob_node<word_graph_node*,edgelist>(0);
	    f->name(n->name());
	    f->netflow(n->nin() - n->nout());
	    f->data(n);
	    tpg.new_node(f);
	  }

	  trans_prob_node<word_graph_node*,edgelist> *t;
	  if (!(t=tpg.find(se.n->name()))) {
	    t = new trans_prob_node<word_graph_node*,edgelist>(0);
	    t->name(se.n->name());
	    t->netflow(se.n->nin() - se.n->nout());
	    t->data(se.n);
	    tpg.new_node(t);
	  }
	  
	  trans_prob_edge<word_graph_node*,edgelist> *e;
	  e = new trans_prob_edge<word_graph_node*,edgelist>;
	  e->from(f);
	  e->to(t);
	  e->cost(se.dist);
	  e->data(se.el);
	  tpg.new_edge(e);

	} else {

	  edgelist_it eit = se.n->out().begin();
	  while (eit != se.n->out().end()) {
	    if (se.dist+(*eit)->sequence().length()<mersize) {
	      stackelt se1((*eit)->to(),
			   se.dist+(*eit)->sequence().length(),
			   se.el);
	      se1.el.push_back(*eit);
	      stack.push_front(se1);
	    }
	    ++eit;
	  }
	}
      }
      //     } else if ((*nit)->nin() < (*nit)->nout()) {
      //       word_graph_node *n = (*nit);
      //       trans_prob_node *t;
      //       if (!(t=tpg.find(n->name()))) {
      // 	t = new trans_prob_node;
      // 	t->name(n->name());
      // 	t->netflow(n->nin() - n->nout());
      // 	t->original_node(n);
      // 	tpg.new_node(t);
      //       }
    }
    ++nit;
  }

  checkpoint;

  long unsigned int supply_nodes=0;
  long unsigned int total_supply=0;
  long unsigned int demand_nodes=0;
  long unsigned int total_demand=0;
  std::list<trans_prob_node<word_graph_node*,edgelist>*>::iterator tn;
  tn = tpg.nodes().begin();
  while (tn != tpg.nodes().end()) {
    if ((*tn)->netflow() > 0) {
      supply_nodes ++;
      total_supply += (*tn)->netflow();
    } else {
      demand_nodes ++;
      total_demand += -((*tn)->netflow());
    }
    ++tn;
  }
  
  checkpoint;
  cerr << "Supply nodes: " << supply_nodes << endl;
  cerr << "Total supply: " << total_supply << endl;
  cerr << "Demand nodes: " << demand_nodes << endl;
  cerr << "Total demand: " << total_demand << endl;

  trans_prob_node<word_graph_node*,edgelist> *dummyf 
    = new trans_prob_node<word_graph_node*,edgelist>(0);
  dummyf->name(tpg.new_label());
  // dummyf->netflow(total_demand);
  dummyf->data(0);
  tpg.new_node(dummyf);

  trans_prob_node<word_graph_node*,edgelist> *dummyt 
    = new trans_prob_node<word_graph_node*,edgelist>(0);
  dummyt->name(tpg.new_label());
  // dummyt->netflow(-total_supply);
  dummyt->data(0);
  tpg.new_node(dummyt);
  
  if (total_demand > total_supply) { 
    dummyf->netflow(total_demand - total_supply);
  } else if (total_demand < total_supply) {
    dummyt->netflow(total_demand - total_supply);
  }
  
  trans_prob_edge<word_graph_node*,edgelist> *e 
    = new trans_prob_edge<word_graph_node*,edgelist>;
  e->from(dummyf);
  e->to(dummyt);
  e->cost(mersize);
  tpg.new_edge(e);

  checkpoint;
  tn = tpg.nodes().begin();
  while (tn != tpg.nodes().end()) {
    if ((*tn) == dummyf ||
	(*tn) == dummyt) {
      ++tn;
      continue;
    }
    e = new trans_prob_edge<word_graph_node*,edgelist>;
    if ((*tn)->netflow() > 0) {
      e->from(*tn);
      e->to(dummyf);
    } else {
      e->from(dummyt);
      e->to(*tn);
    }
    tpg.new_edge(e);
    ++tn;
  }
  
  checkpoint;
  cerr << "Number of nodes: " << tpg.nnode() << endl;
  cerr << "Number of edges: " << tpg.nedge() << endl;
   
  if (optimize) {
    checkpoint;
    tpg.solver(trans_prob_graph<word_graph_node*,edgelist>::netflo);
    tpg.solve();
    checkpoint; 
    cerr << "Solution: " << tpg.evaluate_solution() << endl;
  } else {
    checkpoint;
    tpg.solver(trans_prob_graph<word_graph_node*,edgelist>::vogels);
    tpg.solve();
    checkpoint;
    cerr << "Solution: " << tpg.evaluate_solution() << endl;
  }
  checkpoint;
 
  std::list<trans_prob_edge<word_graph_node*,edgelist>*>::iterator te;
  te = tpg.edges().begin();
  while (te!=tpg.edges().end()) {
    if ((*te)->flow() > 0 &&
	(*te)->from()->data() != 0 &&
	(*te)->to()->data() != 0) {
      std::string seq;
      edgelist const & el = (*te)->data();
      edgelist_cit eit = el.begin();
      while (eit!=el.end()) {
	seq += (*eit)->sequence();
	++eit;
      }
      for (int i=0;i<(*te)->flow();i++) {
	word_graph_edge * e = new word_graph_edge;
	e->from((*te)->from()->data());
	e->to((*te)->to()->data());
	e->sequence(seq);
	e->unset_mark();
	new_edge(e);
      }      
    }
    ++te;
  }

  return 0;
}

std::map<std::string,std::pair<trans_prob_node<word_graph_node*,std::string>*,trans_prob_node<word_graph_node*,std::string>* > >::iterator 
add_widget(trans_prob_graph<word_graph_node*,std::string> & tpg,
		       std::map<std::string,std::pair<trans_prob_node<word_graph_node*,std::string>*,trans_prob_node<word_graph_node*,std::string>* > >  & dcm, 
		       std::string const & key, 
		       long unsigned int cost) {

  trans_prob_node<word_graph_node*,std::string> *dummyf 
    = new trans_prob_node<word_graph_node*,std::string>(0);
  dummyf->name(tpg.new_label());
  // dummyf->netflow(total_demand);
  tpg.new_node(dummyf);

  trans_prob_node<word_graph_node*,std::string> *dummyt 
    = new trans_prob_node<word_graph_node*,std::string>(0);
  dummyt->name(tpg.new_label());
  tpg.new_node(dummyt);

  trans_prob_edge<word_graph_node*,std::string> *e 
    = new trans_prob_edge<word_graph_node*,std::string>;
  e->from(dummyf);
  e->to(dummyt);
  e->cost(cost);
  e->data(key);
  tpg.new_edge(e);

  std::pair<std::map<std::string,std::pair<trans_prob_node<word_graph_node*,std::string>*,trans_prob_node<word_graph_node*,std::string>* > >::iterator,bool> res;
  res = dcm.insert(make_pair(key,make_pair(dummyf,dummyt)));
  assert(res.second);
  return res.first;
};

/* 
long unsigned int 
word_graph::find_jumpers(unsigned int mersize, bool optimize, bool verbose) {
  
  long unsigned int supply_nodes=0;
  long unsigned int total_supply=0;
  long unsigned int demand_nodes=0;
  long unsigned int total_demand=0;
  nodelist_it nit = nodes().begin();
  while (nit != nodes().end()) {
    if ((*nit)->nin() > (*nit)->nout()) {
      supply_nodes++;
      total_supply += (*nit)->nin() - (*nit)->nout();
    } else if ((*nit)->nin() < (*nit)->nout()) {
      demand_nodes++;
      total_demand += (*nit)->nout() - (*nit)->nin();
    }
    ++nit;
  }
  
  checkpoint;
  cerr << "Total nodes: " << nnode() << endl;
  cerr << "Supply nodes: " << supply_nodes << endl;
  cerr << "Total supply: " << total_supply << endl;
  cerr << "Demand nodes: " << demand_nodes << endl;
  cerr << "Total demand: " << total_demand << endl;

  trans_prob_graph<word_graph_node*,std::string> tpg;

  char* buffer = new char[mersize+1];
  CharStarChars cp(buffer);

  keyword_tree<ktnode_lsuffix> kt;

  checkpoint;

  std::vector<word_graph_node*> ndarray(10000,0);
  nit=nodes().begin();
  while (nit != nodes().end()) {
    if ((*nit)->nin() < (*nit)->nout()) {
      ndarray.push_back(*nit);
      kt.add_pattern((*nit)->sequence(),ndarray.size()-1);
    }
    ++nit;
  }
  
  checkpoint;

  kt.init(cp);

  // For each node, if it has defficient degree, try all of its
  // suffixes to see if these are valid prefixes for some other nodes.

  checkpoint;

  nit = nodes().begin();
  while (nit != nodes().end()) {
    // cerr << "Consider " << (*nit)->nin << " > ["
    // << (*nit)->name << "] < " << (*nit)->nout 
    // << endl;
    if ((*nit)->nin() > (*nit)->nout()) {
      // Node is unbalanced, more in than out..
      // cerr << "Node " << (*nit)->name() << " " << (*nit)->sequence() << endl;
      word_graph_node *n = (*nit);

      strcpy(buffer,n->sequence().c_str());
      cp.reset();
      kt.reset();

      long unsigned int added_edges=0;
      pattern_hit_vector phl;
      pattern_hit_vector::iterator it;
      if (kt.find_suffixes(cp,phl,7)) {
	// checkpoint;
	it = phl.rbegin();
	while (it != phl.rend()) {
	  // checkpoint;
	  word_graph_node *n1 = ndarray[it->value()->id()];
	  unsigned int l = it->key();
	  
	  trans_prob_node<word_graph_node*,std::string> *f;
	  if (!(f=tpg.find(n->name()))) {
	    f = new trans_prob_node<word_graph_node*,std::string>(0);
	    f->name(n->name());
	    f->netflow(n->nin() - n->nout());
	    f->data(n);
	    tpg.new_node(f);
	  }
	  
	  trans_prob_node<word_graph_node*,std::string> *t;
	  if (!(t=tpg.find(n1->name()))) {
	    t = new trans_prob_node<word_graph_node*,std::string>(0);
	    t->name(n1->name());
	    t->netflow(n1->nin() - n1->nout());
	    t->data(n1);
	    tpg.new_node(t);
	  }
	  
	  trans_prob_edge<word_graph_node*,std::string> *e;
	  e = new trans_prob_edge<word_graph_node*,std::string>;
	  e->from(f);
	  e->to(t);
	  e->cost(mersize-1-l);
	  e->data(n1->sequence().substr(l,mersize-1-l));
	  tpg.new_edge(e);


	  ++added_edges;
	  if (added_edges > 10) break;
// 	  checkpoint;
//  	  cerr << "Creating edge between nodes ["
//  	       << n->name() << "] " << n->sequence() << " and [" 
//  	       << n1->name() << "] " 
//  	       << n1->sequence() << ", cost " << mersize-1-l << endl;
	  
	  --it;
	}
      } else {
	trans_prob_node<word_graph_node*,std::string> *f;
	if (!(f=tpg.find(n->name()))) {
	  f = new trans_prob_node<word_graph_node*,std::string>(0);
	  f->name(n->name());
	  f->netflow(n->nin() - n->nout());
	  f->data(n);
	  tpg.new_node(f);
	}
      }
    } else if ((*nit)->nin() < (*nit)->nout()) {
      trans_prob_node<word_graph_node*,std::string> *t;
      if (!(t=tpg.find((*nit)->name()))) {
	t = new trans_prob_node<word_graph_node*,std::string>(0);
	t->name((*nit)->name());
	t->netflow((*nit)->nin() - (*nit)->nout());
	t->data(*nit);
	tpg.new_node(t);
      }
    }
    ++nit;
  }

  checkpoint;

  supply_nodes=0;
  total_supply=0;
  demand_nodes=0;
  total_demand=0;
  std::list<trans_prob_node<word_graph_node*,std::string>*>::iterator tn;
  tn = tpg.nodes().begin();
  while (tn != tpg.nodes().end()) {
    if ((*tn)->netflow() > 0) {
      supply_nodes ++;
      total_supply += (*tn)->netflow();
    } else {
      demand_nodes ++;
      total_demand += -((*tn)->netflow());
    }
    ++tn;
  }
  
  checkpoint;
  cerr << "Supply nodes: " << supply_nodes << endl;
  cerr << "Total supply: " << total_supply << endl;
  cerr << "Demand nodes: " << demand_nodes << endl;
  cerr << "Total demand: " << total_demand << endl;

  typedef std::map<std::string,std::pair<trans_prob_node<word_graph_node*,std::string>*,trans_prob_node<word_graph_node*,std::string>* > > dupcostmap;
  dupcostmap dcm;
  dupcostmap::iterator it;

  trans_prob_node<word_graph_node*,std::string> *dummyf 
    = new trans_prob_node<word_graph_node*,std::string>(0);
  dummyf->name(tpg.new_label());
  // dummyf->netflow(total_demand);
  tpg.new_node(dummyf);

  trans_prob_node<word_graph_node*,std::string> *dummyt 
    = new trans_prob_node<word_graph_node*,std::string>(0);
  dummyt->name(tpg.new_label());
  // dummyt->netflow(-total_supply);
  tpg.new_node(dummyt);

  if (total_demand > total_supply) { 
    dummyf->netflow(total_demand - total_supply);
  } else if (total_demand < total_supply) {
    dummyt->netflow(total_demand - total_supply);
  }
  
  trans_prob_edge<word_graph_node*,std::string> *e 
    = new trans_prob_edge<word_graph_node*,std::string>;
  e->from(dummyf);
  e->to(dummyt);
  e->cost(mersize-1);
  tpg.new_edge(e);

//   cerr << "New edge [" 
//        << e->from()->name()
//        << "] -- " 
//        << e->data() 
//        << "("
//        << e->cost()
//        << ") -> ["
//        << e->to()->name()
//        << "]"
//        << endl;

  dcm.insert(make_pair("",make_pair(dummyf,dummyt)));

  checkpoint;
  tn = tpg.nodes().begin();
  while (tn != tpg.nodes().end()) {
    // Check to make sure this is a real node! 
    if ((*tn)->data()==0) {
      ++tn;
      continue;
    }
    std::string nodeseq = (*tn)->data()->sequence();
    std::string key="";
    long unsigned int bound=7;
    if ((*tn)->nin() + (*tn)->nout() > 3) {
      bound -= (*tn)->nin() + (*tn)->nout() - 3;
    }
    if (bound == 0) {
      bound = 1;
    }
    if ((*tn)->netflow() > 0) {
      for (int j=bound-1;j>=0;j--) {
	key = nodeseq.substr(nodeseq.length()-j,j);
	if ((it=dcm.find(key))==dcm.end()) {
	  it = add_widget(tpg,dcm,key,mersize-1-j);
	}
	e = new trans_prob_edge<word_graph_node*,std::string>;
	e->from(*tn);
	e->to((*it).second.first);
	tpg.new_edge(e);
      }
    } else {
      for (int j=bound-1;j>=0;j--) {
	key = nodeseq.substr(0,j);
	if ((it=dcm.find(key))==dcm.end()) {
	  it = add_widget(tpg,dcm,key,mersize-1-j);
	}
	e = new trans_prob_edge<word_graph_node*,std::string>;
	e->from((*it).second.second);
	e->to(*tn);
	tpg.new_edge(e);
      }
    }
    ++tn;
  }
  
  checkpoint;
  cerr << "Number of nodes: " << tpg.nnode() << endl;
  cerr << "Number of edges: " << tpg.nedge() << endl;

  std::list<trans_prob_edge<word_graph_node*,std::string>*>::iterator e3it;
  long unsigned int total_edges=0;
  tn = tpg.nodes().begin();
  while (tn != tpg.nodes().end()) {
    if ((*tn)->data() != 0 && (*tn)->nout() > 0) {
      // Real node!
      e3it = (*tn)->out().begin();
      while (e3it!= (*tn)->out().end()) {
	if ((*e3it)->to()->data() != 0) {
	  total_edges++;
	} else {
	  total_edges += (*e3it)->to()->out().front()->to()->nout();
	}
	++e3it;
      }
    }
    ++tn;
  }  

  checkpoint;
  cerr << "Total number of implicit edges: " << total_edges << endl;
   
  if (optimize) {
    checkpoint;
    tpg.solver(trans_prob_graph<word_graph_node*,std::string>::netflo);
    tpg.solve();
    checkpoint; 
    cerr << "Solution: " << tpg.evaluate_solution() << endl;
  } else {
    checkpoint;
    tpg.solver(trans_prob_graph<word_graph_node*,std::string>::vogels);
    tpg.solve();
    checkpoint;
    cerr << "Solution: " << tpg.evaluate_solution() << endl;
  }
  checkpoint;
 
  std::list<trans_prob_edge<word_graph_node*,std::string>*>::iterator te,eiti,eito;
  te = tpg.edges().begin();
  while (te!=tpg.edges().end()) {
    if ((*te)->flow() > 0 &&
	(*te)->from()->data() != 0 &&
	(*te)->to()->data() != 0) {
      for (int i=0;i<(*te)->flow();i++) {
	word_graph_edge * e = new word_graph_edge;
	e->from((*te)->from()->data());
	e->to((*te)->to()->data());
	e->sequence((*te)->data());
	e->unset_mark();
	new_edge(e);
      }      
    } else if ((*te)->flow() > 0 &&
	       (*te)->from()->data() == 0 &&
	       (*te)->to()->data() == 0 &&
	       (*te)->cost() < mersize-1) {
      // This is one of the dummy edges for one or two characters in common
      eiti = (*te)->from()->in().begin();
      while (eiti != (*te)->from()->in().end() && (*eiti)->flow() == 0) ++eiti;
      long unsigned int flowin = 0;
      eito = (*te)->to()->out().begin();
      while (eito != (*te)->to()->out().end() && (*eito)->flow() == 0) ++eito;

      long unsigned int flowout = 0;
      while (eito != (*te)->to()->out().end() &&
	     eiti != (*te)->from()->in().end()) {
	long unsigned int send;
	if ((*eiti)->flow() < (*eito)->flow()) {
	  send = (*eiti)->flow();
	} else {
	  send = (*eito)->flow();
	}
	flowin += send;
	flowout += send;
	for (int i=0;i<send;i++) {
	  word_graph_edge * e = new word_graph_edge;
	  e->from((*eiti)->from()->data());
	  e->to((*eito)->to()->data());
	  e->sequence((*eito)->to()->data()->sequence().substr(mersize-1-(*te)->cost()));
	  e->unset_mark();
	  new_edge(e);
	}      
	if (flowin >= (*eiti)->flow()) {
	  ++eiti;
	  flowin=0;
	  while (eiti != (*te)->from()->in().end() && (*eiti)->flow() == 0) ++eiti;	  
	}
	if (flowout >= (*eito)->flow()) {
	  ++eito;
	  flowout=0;
	  while (eito != (*te)->to()->out().end() && (*eito)->flow() == 0) ++eito;
	}
      }
    }
    ++te;
  }

  return 0;
}

*/

struct Options {
  std::string graphfile;
  std::string sequencefile;
  char eos_char;
  unsigned int mersize;
  ostream *out;
  bool delout;
  bool verbose;
  bool optimize;
  bool redundant;
  bool remove_eos;
  bool allownew;
  Options(int argc, char *argv[]);
  ~Options();
  void usage(char *msg=NULL);
};

Options::Options(int argc, char *argv[]) : out(&cout) {
  signed char c;
  optarg = NULL;
  graphfile = "";
  sequencefile = "";
  eos_char = '$';
  mersize = 30;
  delout = false;
  verbose = false;
  optimize = false;
  redundant = false;
  remove_eos = true;
  allownew = false;
  while ((c = getopt(argc, argv, "g:s:k:E:o:vORNeh")) != -1)
    switch (c) {
    case 'g':
      graphfile = optarg;
      break;
    case 's':
      sequencefile = optarg;
      break;
    case 'k':
      mersize = atoi(optarg);
      break;
    case 'O':
      optimize = true;
      break;
    case 'R':
      redundant = true;
      break;
    case 'N':
      allownew = true;
      break;
    case 'e':
      remove_eos = false;
      break;
   case 'o':
      if (delout) {
	delete out;
      }
      out = new ofstream(optarg);
      delout = true;
      break;
    case 'E': {
      int eos_chari;
      if (!sscanf(optarg,"%i",&eos_chari)) {
	usage("Invalid end-of-sequence specification.\n");
      }
      eos_char = (char)eos_chari;
    }
      break;
    case 'v':
      verbose = true;
      break; 
    case 'h':
    default :
      usage();
    }
    if ((graphfile == "" || sequencefile == "")) usage();
    if (allownew && redundant) usage();
}

Options::~Options() {
  if (delout) {
    delete out;
  }
}

void 
Options::usage(char *message) {

  if (message != NULL && strlen(message) > 0) {
    cerr << message << endl;
    cerr << endl;
  }
  cerr << "Usage: walk_graph [options]\n\n";
  cerr << "Options: \n";
  cerr << "  -g <graph-file> Word graph file. Required.\n";
  cerr << "  -s <seq-file>   Word graph sequence file. Required.\n";
  cerr << "  -k <mer-size>   Mersize of word graph. Default: 30.\n";
  cerr << "  -e              Retain end-of-seqence character in output. Default:false.\n";
  cerr << "  -E <int>        End-of-sequence character. Default: \'$\'\n";
  cerr << "  -o <out-file>   Output file. Default is standard out.\n";
  cerr << "  -R              Permit redundant k-mers to be output. Default: false.\n";
  cerr << "                  At most one of -R and -N can be specified.\n";
  cerr << "  -N              Permit new k-mers to be output. Default: false.\n";
  cerr << "                  At most one of -R and -N can be specified.\n";
  cerr << "  -O              Find optimal redundant or new compression.\n";
  cerr << "  -v              Verbose.\n";
  cerr << "  -h              Help.\n";
  exit(1);
}

int
main(int argc, char **argv) {
  
  Options opt(argc,argv);

  word_graph g;
  checkpoint;
  g.read(opt.graphfile,opt.sequencefile,opt.mersize);
  checkpoint;
  if (opt.remove_eos) {
    // g.remove_eos(opt.eos_char,opt.mersize,opt.verbose);
  }
  checkpoint;
  // g.uniquify_nodes(opt.verbose);
  checkpoint;
  // g.uniquify_edges(opt.verbose);
  checkpoint;
  // g.remove_trivial_nodes(opt.verbose);
  checkpoint;
  // g.uniquify_edges(opt.verbose);
  checkpoint;
  g.print_stats();
  checkpoint;
  if (opt.redundant) {
    g.find_joiners(opt.mersize,opt.optimize,opt.verbose);
    checkpoint;
  } 
  // checkpoint;
  // g.join_components(opt.eos_char,opt.allownew,opt.verbose);
  checkpoint;
  g.balance_nodes(/*opt.eos_char*/'|',opt.allownew,opt.verbose);
  if (opt.verbose) {
    checkpoint;
    g.dump(*opt.out,opt.mersize);
  }
  checkpoint;
  g.writeseq(*opt.out,opt.verbose);
}

