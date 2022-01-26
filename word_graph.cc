
#include "word_graph.h"
#include "fasta_io.h"
#include "fasta_io.t"
#include "select.h"
#include "select.t"
#include "trans_prob.h"
#include "math.h"

#include <malloc.h>
#include <map>
#include <iomanip>
#include <set>

typedef std::list<word_graph_node*> nodelist;
typedef std::list<word_graph_edge*> edgelist;
typedef std::list<word_graph_node*>::iterator nodelist_it;
typedef std::list<word_graph_edge*>::iterator edgelist_it;
typedef std::list<word_graph_node*>::const_iterator nodelist_cit;
typedef std::list<word_graph_edge*>::const_iterator edgelist_cit;
typedef std::map<word_graph_node*,word_graph_node*> nodenodemap;
typedef std::map<word_graph_edge*,word_graph_edge*> edgeedgemap;
typedef std::map<word_graph_node*,word_graph_node*>::iterator nodenodemap_it;
typedef std::map<word_graph_edge*,word_graph_edge*>::iterator edgeedgemap_it;
typedef std::map<word_graph_node*,word_graph_node*>::const_iterator nodenodemap_cit;
typedef std::map<word_graph_edge*,word_graph_edge*>::const_iterator edgeedgemap_cit;

unsigned short fixedlen_word_graph_node::length_(0);
char restart_word_graph_edge::eos_char_('\n');

FILE_POSITION_TYPE 
word_graph_node::seq_end() const {
  edgelist_cit ei;
  ei = out().begin();
  while (ei != out().end()) {
    if ((*ei)->isreal()) {
      break;
    }
    ++ei;
  }
  if (ei == out().end()) {
    ei = in().begin();
    while (ei != in().end()) {
      if ((*ei)->isreal()) {
	break;
      }
      ++ei;
    }
    assert(ei != in().end());
    return ((real_word_graph_edge*)(*ei))->seq_end();
  } else {
    return ((real_word_graph_edge*)(*ei))->seq_start();
  }
}

void
word_graph::read(std::string const & gf, int mersize, int ctin, int sign, bool readcounts) {
  long unsigned int nodefrom;
  long unsigned int nodeto;
  FILE_POSITION_TYPE fseqst,tseqst;
  FILE_POSITION_TYPE fseqed,tseqed;
  unsigned short flen,tlen;
  long unsigned int count;
  signed int edgeno;

  // checkpoint;
  // std::cerr << "CtIn: " << ctin << endl;
  // std::cerr << "Sign: " << sign << endl;

  if (mersize > 0) {
    fixedlen_word_graph_node::length(mersize);
  }

  ifstream gis(gf.c_str());
  if (!gis) {
    checkpoint;
    cerr << "Fatal error: Can't open graph file " << gf << "." << endl;
    exit(1);
  }
  
  std::string line;
  int lineno = 0;
  while (getline(gis,line)) {
    // if (lineno % (50*1000) == 0) {
    //   cerr << endl << setw(10) << lineno << " ";
    // }
    // if (lineno % 1000 == 0) {
    //   cerr << ".";
    // }
    lineno++;
    istringstream ss(line);
    if (mersize <= 0) {
      ss >> nodefrom >> nodeto >> fseqst >> fseqed >> tseqst >> tseqed >> count;
      if (ctin != 0 &&
	  ((sign == 0 && count == ctin) || 
	   (sign < 0 && count < ctin) || 
	   (sign > 0 && count > ctin))) {
	continue;
      }
      flen = fseqed - fseqst;
      tlen = tseqed - tseqst;
      varlen_word_graph_node *f;
      if (!(f=(varlen_word_graph_node*)find(nodefrom)) && count) {
	f = new varlen_word_graph_node;
	f->name(nodefrom);
	f->length(flen);
	new_node(f);
      }
      varlen_word_graph_node *t;
      if (!(t=(varlen_word_graph_node*)find(nodeto)) && count) {
	t = new varlen_word_graph_node;
	t->name(nodeto);
	t->length(tlen);
	new_node(t);
      }
      if (f && t) {
	if (readcounts) {
	  if (count) {
	    real_count_word_graph_edge *e = new real_count_word_graph_edge;
	    e->from(f);
	    e->to(t);
	    e->seq_end(tseqed);
	    e->length(tseqed-fseqed);
	    e->count(count);
	    new_edge(e);
	  } else {
	    sim_word_graph_edge *e = new sim_word_graph_edge;
	    e->from(f);
	    e->to(t);
	    new_edge(e);
	  }
	} else {
	  real_word_graph_edge *e = new real_word_graph_edge;
	  e->from(f);
	  e->to(t);
	  e->seq_end(tseqed);
	  e->length(tseqed-fseqed);
	  new_edge(e);
	}
      }
    } else {
      ss >> nodefrom >> nodeto >> fseqed >> tseqed >> count;
      if (ctin != 0 &&
	  ((sign == 0 && count != ctin) || 
	   (sign < 0 && count >= ctin) || 
	   (sign > 0 && count <= ctin))) {
	// checkpoint;
	// std::cerr << "Skip: " << nodefrom << " " <<  nodeto << " " <<  fseqed << " " << tseqed << " " << count << endl;
	continue;
      } /* else {
	checkpoint;
	std::cerr << "Keep: " << nodefrom << " " <<  nodeto << " " <<  fseqed << " " << tseqed << " " << count << endl;
	}*/
      fixedlen_word_graph_node *f;
      if (!(f=(fixedlen_word_graph_node*)find(nodefrom)) && count) {
	f = new fixedlen_word_graph_node;
	f->name(nodefrom);
	new_node(f);
      }
      fixedlen_word_graph_node *t;
      if (!(t=(fixedlen_word_graph_node*)find(nodeto)) && count) {
	t = new fixedlen_word_graph_node;
	t->name(nodeto);
	new_node(t);
      }
      if (f && t) {
	if (readcounts) {
	  if (count) {
	    real_count_word_graph_edge *e = new real_count_word_graph_edge;
	    e->from(f);
	    e->to(t);
	    e->seq_end(tseqed);
	    e->length(tseqed-fseqed);
	    e->count(count);
	    new_edge(e);
	  } else {
	    sim_word_graph_edge *e = new sim_word_graph_edge;
	    e->from(f);
	    e->to(t);
	    new_edge(e);
	  }
	} else {
	  real_word_graph_edge *e = new real_word_graph_edge;
	  e->from(f);
	  e->to(t);
	  e->seq_end(tseqed);
	  e->length(tseqed-fseqed);
	  new_edge(e);
	}
      }
    }
  }
  gis.close();
  free_labelmap();
  cerr << endl;
}

void
word_graph_node::dump(ostream & os, CharacterProducer & cp, bool mark) const {
  edgelist_cit eit=in().begin();
  while (eit != in().end()) {
    if ((*eit)->isreal()) {
      os << "    <- ";
    } else {
      os << "    <~ ";
    }
    os << " [" << (*eit)->from()->name();
    if (!mark) {
      if ((*eit)->from()->mark()) {
	os << "*";
      } 
    }
    os << "]" ;
    if (!mark) {
      if ((*eit)->mark()) {
	os << "*";
      } else {
	os << " ";
      }
    }
    os << (*eit)->sequence(cp);
    if ((*eit)->count() != 0) {
      os << "(" << (*eit)->count() << ")";
    }
    os << endl;
    if (mark) {
      (*eit)->set_mark();
    }
    ++eit;
  }
  os << sequence(cp) << " [" << name();
  if (!mark) {
    if (this->mark()) {
      os << "*";
    } 
  }
  os << "] " << endl;
  eit=out().begin();
  while (eit != out().end()) {
    if ((*eit)->isreal()) {
      os << "    -> ";
    } else {
      os << "    ~> ";
    }
    os << " [" << (*eit)->to()->name();
    if (!mark) {
      if ((*eit)->to()->mark()) {
	os << "*";
      } 
    }
    os << "]";
    if (!mark) {
      if ((*eit)->mark()) {
	os << "*";
      } else {
	os << " ";
      }
    }
    os << (*eit)->sequence(cp);
    if ((*eit)->count() != 0) {
      os << "(" << (*eit)->count() << ")";
    }
    os << endl;
    if (mark) {
      (*eit)->set_mark();
    }
    ++eit;
  }
  os << endl;
}

void
word_graph::dump(ostream & os, CharacterProducer & cp, bool mark) {
  os << "Number of nodes: " << nnode() << endl;
  os << "Number of edges: " << nedge() << endl;
  if (mark) { 
    clear_edge_marks();
  }
  
  nodelist_cit nit=nodes().begin();
  while (nit != nodes().end()) {
    (*nit)->dump(os,cp,mark);
    ++nit;
  }
  if (mark) {
    edgelist_cit eit=edges().begin();
    while (eit != edges().end()) {
      if (!(*eit)->mark()) {
	os << "Unattached edge: " 
	   << " [" << (*eit)->from()->name() << "] " 
	   << "->" 
	   << " [" << (*eit)->to()->name() << "] " 
	   << (*eit)->sequence(cp) << endl;
      }
      ++eit;
    }
    clear_edge_marks();
  }
}

void word_graph::print_stats() {

  std::cerr << "CSBH-graph nodes: " << nnode() << endl;
  std::cerr << "CSBH-graph edges: " << nedge() << endl;
  
  long unsigned int l=0;
  long unsigned int edge_length=0;
  long unsigned int suppressed_nodes=0;
  edgelist_it eit=edges().begin();
  while (eit != edges().end()) {
    l = (*eit)->length();
    edge_length += l;
    suppressed_nodes += l - 1;
    ++eit;
  }
  std::cerr << "CSBH-graph total edge length: " << edge_length << endl;
  std::cerr << "Suppressed SBH-graph nodes: " << suppressed_nodes << endl;

  components<word_graph_node> cmps(nnode());
  dfs(cmps);

  std::cerr << "CSBH-graph components: " << cmps.number() << endl;
  
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
  long unsigned int comp_line=0;
  long unsigned int min_node_length=MAXINT;
  long unsigned int total_bal_plus_length=0;
  long unsigned int nbaln=0;
  long unsigned int ntrivial=0;
  long unsigned int nterminal=0;
  long unsigned int comp_line_len=0;
  
  for (long unsigned int i=1;i<=cmps.number();i++) {
    comp_bal_plus = 0;
    comp_bal_minus = 0;
    comp_n_plus = 0;
    comp_n_minus = 0;
    bool line = true;
    min_node_length = MAXINT;
    components<word_graph_node>::iterator cit(cmps.nodes_begin(i));
    components<word_graph_node>::iterator cit1(cmps.nodes_end(i));
    while (cit != cit1) {
      if ((*cit)->nin() > (*cit)->nout()) {
	comp_n_plus++;
	comp_bal_plus += (*cit)->nin() - (*cit)->nout();
	if ((*cit)->nout() == 0) {
	  nterminal++;
	}
      } else if ((*cit)->nin() < (*cit)->nout()) {
	comp_n_minus++;
	comp_bal_minus += (*cit)->nout() - (*cit)->nin();
	total_bal_plus_length += 
	  ((*cit)->nout() - (*cit)->nin())*(((*cit)->length()+1));
	if ((*cit)->nin() == 0) {
	  nterminal++;
	}
      } else if ( (*cit)->nin() != 1 ) {
	nbaln++;
      } else if ( (*cit)->nin() == 1 ) {
	ntrivial++;
	nbaln++;
      }
      if ((*cit)->nin() > 1 ||
	  (*cit)->nout() > 1) {
	line = false;
      }
      if (min_node_length > ((*cit)->length()+1)) {
	min_node_length = (*cit)->length()+1;
      }
      ++cit;
    }
    if (comp_n_plus > 0 || comp_n_minus > 0) {
      nunbalanced++;
    } else {
      nbalanced++;
      total_bal_plus_length += min_node_length;
    }
    if (line) {
      comp_line++;
      comp_line_len += cmps.number(i)-1;
    }
    total_bal_plus += comp_bal_plus;
    total_bal_minus += comp_bal_minus;
    total_n_plus += comp_n_plus;
    total_n_minus += comp_n_minus;
  }

  std::cerr << "CSBH-graph unbalanced components: " << nunbalanced << endl;
  std::cerr << "CSBH-graph balanced components: " << nbalanced << endl;
  
  std::cerr << "Degree surplus nodes: " << total_n_plus << endl;
  std::cerr << "Degree deficit nodes: " << total_n_minus << endl;
  std::cerr << "Total degree surplus: " << total_bal_plus << endl;
  std::cerr << "Total degree deficit: " << total_bal_minus << endl;
  std::cerr << "Line components: " << comp_line << endl;
  std::cerr << " Average length: " << ((double)comp_line_len)/comp_line << endl;
  std::cerr << "Balanced nodes: " << nbaln << endl;
  std::cerr << "Trivial nodes: " << ntrivial << endl;
  std::cerr << "Terminal nodes: " << nterminal << endl;
  std::cerr << "Restart sequence: " << total_bal_plus_length << endl;
} 

bool
word_graph::check_out_edges(ostream &os, CharacterProducer & cp, bool verbose) {

  std::vector<bool> seen;
  seen.resize(255,false);

  nodelist_it nit=nodes().begin();
  while (nit != nodes().end()) {
    fill(seen.begin(),seen.end(),false);
    edgelist_it eit=(*nit)->out().begin();
    while (eit != (*nit)->out().end()) {
      if (seen[(*eit)->sequence(cp)[0]]) {
	checkpoint;
	cerr << "Node " << (*nit)->name() << ": " << (*nit)->sequence(cp) << " "
	     << "has two edges out that start with " << (*eit)->sequence(cp)[0] << ".";
	return false;
      }
      seen[(*eit)->sequence(cp)[0]] = true;
      ++eit;
    }
    ++nit;
  }
  return true;
}

void 
word_graph::writetrivialpaths(ostream &os, CharacterProducer & cp, char eos_char, bool verbose) {

  os << eos_char;
  edgelist_it eit = edges().begin();
  while (eit != edges().end()) {
    if ((*eit)->from()->nin()!=1 || (*eit)->from()->nout()!=1) {
      os << (*eit)->from()->sequence(cp);
      os << (*eit)->sequence(cp);
       word_graph_node *n=(*eit)->to();
      while (n->nin() == 1 && n->nout() == 1) {
	edgelist_cit eit1 = n->out().begin();
	os << (*eit1)->sequence(cp);
	n = (*eit1)->to();
      }
      os << eos_char;
    }
    ++eit;
  }
}


void 
word_graph::writeseq(ostream &os, CharacterProducer & cp, char eos_char, bool verbose) {
  clear_edge_marks();
  // clear_node_marks();

  // long unsigned int total_tour_length=0;
  // long unsigned int subtourlen = 0;
  components<word_graph_node> cmps(nnode());
  dfs(cmps);

  // checkpoint;

  // long unsigned int opos = 0;
  // Mark beginning of the sequence
  os << eos_char;
  for (long unsigned int i=1;i<=cmps.number();i++) {
    components<word_graph_node>::iterator nit(cmps.nodes_begin(i));
    components<word_graph_node>::iterator nit1(cmps.nodes_end(i));
    components<word_graph_node>::iterator minnit(nit);
    while (nit != nit1) {
      if ((*nit)->nin() < (*nit)->nout()) {
	break;
      } 
      if ((*nit)->length() < (*minnit)->length()) {
	minnit = nit;
      }
      ++nit;
    }
    if (nit == nit1) {
      nit = minnit;
    }

    // checkpoint;

    edgelist tour;
    edgelist_it tip=tour.end();
    edgelist_it sts=tour.end();
    
    word_graph_node *n;
    word_graph_edge *e;
    
    edgelist_it eit;

    n = (*nit);
    n->set_mark();

    bool tourstart = true;

    // checkpoint;

    while (1) {

      edgelist_it eit = n->out().begin();
      while (eit != n->out().end()) {
	if (!(*eit)->mark()) {
	  // Follow edge
	  e = (*eit);
	  if (sts == tour.end()) {
	    sts = tour.insert(tip,e);
	  } else {
	    tour.insert(tip,e);
	  } 
	  e->set_mark();
	  n = e->to();
	  n->set_mark();
	  eit = n->out().begin();
	} else {
	  ++eit;
	}
      }

      // checkpoint;

      tip = sts;
      sts = tour.end();
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

	// checkpoint;

	// output the tip edge and delete from tour...
	if (tourstart) {
	  // checkpoint;
	  os << (*tip)->from()->sequence(cp);
	  // opos += (*tip)->from()->length();
	  tourstart = false;
	}
	// checkpoint;
	os << (*tip)->sequence(cp);
	// opos += (*tip)->length();
	// (*tip)->seq_end(opos);
	// checkpoint;

	tip = tour.erase(tip);
      }
      // checkpoint;
      if (tip == tour.end()) {
	break;
      }
      // checkpoint;
    }

    // checkpoint;


    os << eos_char;
    // opos++;

  }

  // checkpoint;
 
  edgelist_it eit = edges().begin();
  while (eit != edges().end()) {
    if (!(*eit)->mark()) {
      checkpoint;
      std::cerr << "Unmarked edge: " 
	   << "[" << (*eit)->from()->name() << "]" << " --> " 
	   << "[" << (*eit)->to()->name() << "]" << endl;
	
    }
    ++eit;
  }

  // checkpoint;

}


unsigned hash(const char* s)
{
   unsigned int hash = 0;

   while (*s) {
     hash = (*s) + (hash << 6) + (hash << 16) - hash;
     s++;
   }

   return (hash & 0x7FFFFFFF);
}

struct ltuint
{
  bool operator()(const unsigned a, const unsigned b) const
  {
    return (a<b);
  }
};

void 
word_graph::annotateseq(ostream &os, 
			CharacterProducer &ffs, 
			FastaFile<Lazy_Header_SI> &ffa, 
			char eos_char, 
			int transform, 
			int format, 
			bool verbose) {

  typedef std::multimap<unsigned, word_graph_node*,ltuint> nodemap;
  typedef nodemap::iterator nodemap_it;
 
  nodemap nm;
  int minnodelen=MAXINT,maxnodelen=0;
  nodelist_it nit=nodes().begin();
  while (nit != nodes().end()) {
    if (true) {
      std::string nseq = (*nit)->sequence(ffs);
      nm.insert(nodemap::value_type(hash(nseq.c_str()),(*nit)));
      if (nseq.length() < minnodelen) {
	minnodelen = nseq.length();
      }
      if (nseq.length() > maxnodelen) {
	maxnodelen = nseq.length();
      }
    }
    ++nit;
  }

  if (format == 1 || format == 3) {
    os << "track type=wiggle_0 name=";
    if (minnodelen == maxnodelen) {
      os << minnodelen+1 << "-";
    }
    os << "mer-cnt priority=10 graphType=bar autoScale=off viewLimits=1:20" << endl;
  }

  // checkpoint;
  
  nodemap_it nmit,nmitb,nmite;
  std::pair<nodemap_it,nodemap_it> nmeqits;
  char *buffer = new char[maxnodelen+1];
  int j = 0, p = 0;
  while (ffa.fasta_pos(j,p)) {

    int lastval = 0;
    const Lazy_Header_SI & h = ffa.get_header_data(ffa.pos()+1);
    if (format == 0 || format == 2) {
      os << ">";
      os << h.header();
      os << "\n";
    } 

    int i=0;
    for (i=0;i<minnodelen;i++) {
      buffer[i] = ffa.getch();
    }
    i = minnodelen;
    buffer[i] = '\0';
    
    while (i <= maxnodelen) {
      // cerr << j << " " << i << " looking for " << buffer << endl;
      bool found = false;
      nmeqits = nm.equal_range(hash(buffer));
      nmit = nmeqits.first;
      while (nmit != nmeqits.second) {
	if (!strcmp(nmit->second->sequence(ffs).c_str(),buffer)) {
	  found = true;
	  break;
	}
	++nmit;
      }
      if (found) {
	break;
      }
      buffer[i] = ffa.getch();
      i++;
      buffer[i] = '\0';
    }

    if (i > maxnodelen) {
      checkpoint;
      buffer[i-1] = '\0';
      cerr << "Can't find node for beginning of sequence " << j << " that starts with " << buffer << endl;
    } else {
      word_graph_node* n = nmit->second;

      // checkpoint;
      // cerr << "Found node " << n->name() << " that has sequence " << j << ":" << p << ":" << buffer << endl;

      p = i;
      
      if (format == 0) {
	for (int k=0;k<i;k++) {
	  os << (unsigned char)(0+'a');
	}
      } else if (format == 1) {
	// os << h.short_header() << " " << 0 << " " << p << " " << 0 << endl;
      } else if (format == 2) {
	os << n->sequence(ffs);
	os.flush();
      } else if (format == 3) {
	os << h.short_header() << " " << 0 << " ";
	lastval = 0;
      }

      while (1) {
	ffa.fasta_pos(j,p);
	char first_edge_char = ffa.getch();
	
	// cerr << "First edge char: " << first_edge_char << endl;
	
	if (first_edge_char == eos_char) {
	  // checkpoint
	  if ((format == 1 && lastval > 1) || format == 3) {
	    os << p << " " << lastval << endl;
	  }
	  break;
	}

	bool found = false;
	edgelist_it eit=n->out().begin();
	while (eit != n->out().end()) {
	  if ((*eit)->sequence(ffs)[0] == first_edge_char) {
	    break;
	  }
	  ++eit;
	}
	
	if (eit == n->out().end()) {
	  checkpoint;
	  cerr << "Can't find edge out of node " << n->name() << " for character " << first_edge_char << endl;
	  break;
	} 
	
	// checkpoint;
	// cerr << "Found edge leaving node " << n->name() << ": " << n->sequence(ffs) << " starting with character " << first_edge_char  << " to " << (*eit)->to()->name() << ": " << (*eit)->sequence(ffs) << endl;

	int val;
	if (transform == 0) {
	  val = (*eit)->count();
	} else if (transform == 1) {
	  val = (int)floor((log(((float)(*eit)->count()))/log(2.0))+.001)+1;
	} else if (transform == 2) {
	  val = ((*eit)->count()>1)?2:1;
	} 
	if (format == 0) {
	  for (int k=0;k<(*eit)->length();k++) {
	    if (val > 'z'-'a') {
	      os << 'z';
	    } else {
	      os << (char)(val+'a');
	    }
	  }
	} else if (format == 1 || format == 3) {
	  if (val != lastval) {
	    if (lastval > 1 || format == 3) {
	      os << p << " " << lastval << endl;
            }
	    if (val > 1 || format == 3) {
	      os << h.short_header() << " " << p << " ";
	    }
	  }
	} else if (format == 2) {
	  os << (*eit)->sequence(ffs);
	}

	lastval = val;
	n = (*eit)->to();
	p += (*eit)->length();
      }
    }

    j++;
    p = 0;
    if (format == 0 || format == 2) {
      os << "\n";
    }
  }
}


long unsigned int
word_graph::balance_nodes(char eos_char, bool verbose) {

  components<word_graph_node> cmps(nnode());
  dfs(cmps);

  restart_word_graph_edge::eos_char(eos_char);

  for (long unsigned int i=1;i<=cmps.number();i++) {
    
    long unsigned int bplus=0;
    components<word_graph_node>::iterator nit(cmps.nodes_begin(i));
    components<word_graph_node>::iterator nit1(cmps.nodes_end(i));
    while (nit != nit1) {
      if ((*nit)->nin() > (*nit)->nout()) {
	bplus += (*nit)->nin() - (*nit)->nout();
      }
      ++nit;
    }

    long unsigned int edges_added=0;
    
    components<word_graph_node>::iterator nitneg(cmps.nodes_begin(i));
    components<word_graph_node>::iterator nitpos(cmps.nodes_begin(i));
    while (nitneg != nit1 && (*nitneg)->nout() >= (*nitneg)->nin()) ++nitneg;
    while (nitpos != nit1 && (*nitpos)->nout() <= (*nitpos)->nin()) ++nitpos;
    
    int negrem=0;
    if (nitneg != nit1) {
      negrem = (*nitneg)->nin() - (*nitneg)->nout();
    }
    int posrem=0;
    if (nitpos != nit1) {
      posrem = (*nitpos)->nout() - (*nitpos)->nin();
    }

    while (edges_added < bplus - 1 && 
	   nitneg != nit1 && 
	   nitpos != nit1 ) {
      int nnew = negrem;
      if (posrem < negrem) {
	nnew = posrem;
      }
      for (int i1=0;i1<nnew;i1++) {
	edges_added++;
	if (edges_added >= bplus) break;
	restart_word_graph_edge* e = new restart_word_graph_edge;
	e->from(*nitneg);
	e->to(*nitpos);
	new_edge(e);
      }
      negrem -= nnew;
      if (negrem == 0) {
	while (nitneg != nit1 && (*nitneg)->nout() >= (*nitneg)->nin()) 
	  ++nitneg;
	if (nitneg != nit1) {
	  negrem = (*nitneg)->nin() - (*nitneg)->nout();
	}
      }
      posrem -= nnew;
      if (posrem == 0) {
	while (nitpos != nit1 && (*nitpos)->nout() <= (*nitpos)->nin()) 
	  ++nitpos;
	if (nitpos != nit1) {
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
};

class queueelt {
public:
  word_graph_node *n;
  long unsigned int dist;
  queueelt(std::pair<const long unsigned int, word_graph_node*> const & e) {
    n = e.second;
    dist = e.first;
  }
};



long unsigned int 
word_graph::find_joiners(bool optimize, bool verbose) {

  unsigned long problem_count = 0;

  components<word_graph_node> cmps(nnode());
  dfs(cmps);

  for (long unsigned int i=1;i<=cmps.number();i++) {

    // checkpoint;
    // std::cerr << "Component " << i << ": size: " << cmps.number(i) << endl;

    trans_prob_graph<word_graph_node*,bool> tpg;
    long signed int total_demand = 0;
    long signed int total_supply = 0;
    components<word_graph_node>::iterator nit(cmps.nodes_begin(i));
    components<word_graph_node>::iterator nit1(cmps.nodes_end(i));
    while (nit != nit1) {
      
      if ((*nit)->nin() < (*nit)->nout()) {
	word_graph_node *n = (*nit);

	// checkpoint;
	// std::cerr << "Unbalanced node: " << endl;
	// n->dump(std::cerr);

	std::map<word_graph_node*,long unsigned int> seen;
	std::multimap<const long unsigned int, word_graph_node*> queue;
	std::multimap<const long unsigned int, word_graph_node*>::iterator qit;
	long unsigned int const maxqueuesize=100;
	long unsigned int const max_degree=20;
	
	long unsigned int tpg_edges_added=0;
	queue.insert(std::pair<const unsigned long,word_graph_node*>(0,n));
	while (!queue.empty() && tpg_edges_added < max_degree) {
	  // checkpoint;
	  // std::cerr << "queue size: " << queue.size() << endl;
	  qit = queue.begin();
	  queueelt se(*qit);
	  queue.erase(qit);
	  if (se.n->nout() < se.n->nin()) {
	    trans_prob_node<word_graph_node*,bool> *f;
	    if (!(f=tpg.find(se.n->name()))) {
	      f = new trans_prob_node<word_graph_node*,bool>(0);
	      f->name(se.n->name());
	      f->netflow(se.n->nin() - se.n->nout());
	      f->data(se.n);
	      tpg.new_node(f);
	      total_supply += f->netflow();
	    }

	    trans_prob_node<word_graph_node*,bool> *t;
	    if (!(t=tpg.find(n->name()))) {
	      t = new trans_prob_node<word_graph_node*,bool>(0);
	      t->name(n->name());
	      t->netflow(n->nin() - n->nout());
	      t->data(n);
	      tpg.new_node(t);
	      total_demand += - t->netflow();
	    }

	    trans_prob_edge<word_graph_node*,bool> *e0;
	    if (!(e0=tpg.find_one(f,t))) {
	      trans_prob_edge<word_graph_node*,bool> *e;
	      e = new trans_prob_edge<word_graph_node*,bool>;
	      e->from(f);
	      e->to(t);
	      e->cost(se.dist);
	      tpg.new_edge(e);
	      tpg_edges_added++;
	    } else if (se.dist < e0->cost()) {
	      e0->cost(se.dist);
	    }

	  } 

	  if (seen.count(se.n) == 0 ||
	      seen[se.n] > se.dist) {
	    if (seen.count(se.n) == 0) {
	      seen.insert(make_pair(se.n,se.dist));
	    } else {
	      seen[se.n] = se.dist;
	    }
	    edgelist_it eit = se.n->in().begin();
	    while (eit != se.n->in().end()) {
	      long unsigned int newdist = se.dist+(*eit)->length();
	      if (newdist < (n->length()+1) ) {
		if (queue.size() < maxqueuesize || queue.rbegin()->first > newdist) {
		  queue.insert(std::pair<const unsigned long, word_graph_node*>(newdist,(*eit)->from()));
		  while (queue.size() > maxqueuesize) {
		    qit = queue.end();
		    --qit;
		    queue.erase(qit);
		  }
		}
	      } 
	      ++eit;
	    }
	  }
	}
      } 
      ++nit;
    }
    
    if (tpg.nedge() == 0) {
      continue;
    }
    
    // checkpoint;
    
    trans_prob_node<word_graph_node*,bool> *dummyf 
      = new trans_prob_node<word_graph_node*,bool>(0);
    dummyf->name(tpg.new_label());
    if (total_demand > total_supply) {
      dummyf->netflow(total_demand - total_supply);
    } 
    dummyf->data(0);
    tpg.new_node(dummyf);
    
    trans_prob_node<word_graph_node*,bool> *dummyt 
      = new trans_prob_node<word_graph_node*,bool>(0);
    dummyt->name(tpg.new_label());
    if (total_demand < total_supply) {
      dummyt->netflow(total_demand - total_supply);
    }
    dummyt->data(0);
    tpg.new_node(dummyt);
    
    trans_prob_edge<word_graph_node*,bool> *e 
      = new trans_prob_edge<word_graph_node*,bool>;
    e->from(dummyf);
    e->to(dummyt);
    e->cost(0);
    tpg.new_edge(e);

    std::list<trans_prob_node<word_graph_node*,bool>*>::iterator tn;
    long signed int realfrom=0;
    long signed int realto=0;
    long signed int realfromdeg=0;
    long signed int realtodeg=0;
    long signed int costub=0;
    tn = tpg.nodes().begin();
    while (tn != tpg.nodes().end()) {
      if ((*tn) == dummyf ||
	  (*tn) == dummyt) {
	++tn;
	continue;
      }
      e = new trans_prob_edge<word_graph_node*,bool>;
      if ((*tn)->netflow() > 0) {
	e->from(*tn);
	e->to(dummyf);
	e->cost(0);
	realfrom++;
	realfromdeg += (*tn)->netflow();
      } else {
	e->from(dummyt);
	e->to(*tn);
	e->cost((*tn)->data()->length()+1);
	realto++;
	realtodeg += (*tn)->netflow();
	costub += (-(*tn)->netflow())*e->cost();
      }
      tpg.new_edge(e);
      ++tn;
    }


    if (verbose && tpg.nnode() > 10000) {

      checkpoint;
      
      std::cerr << "Number of nodes: " << tpg.nnode() << endl;
      std::cerr << "Real from nodes: " << realfrom << " (" << realfromdeg << ")" << endl;
      std::cerr << "Real to nodes: " << realto << " (" << realtodeg << ")" << endl;
      std::cerr << "Number of edges: " << tpg.nedge() << endl;
      std::cerr << "Real edges: " << tpg.nedge() - realfrom - realto - 1 << endl;
      std::cerr << "Cost Upperbound: " << costub << endl;
 
      // ostrstream ss;
      // ss << ".input" << problem_count++ << ".min" << ends;
      // tpg.write_dimacs_input(ss.str().c_str());

    }

    if (optimize) {
      tpg.solver(trans_prob_graph<word_graph_node*,bool>::cs2);
      tpg.solve();
      if (verbose && tpg.nnode() > 10000) {
        std::cerr << "Solution: " << tpg.evaluate_solution() << endl;
      }
      if (tpg.evaluate_solution() > costub) {
	// tpg.write_dimacs_input("problem_input.txt");
	checkpoint;
	std::cerr << "Bogus solution to optimization problem" << endl;
	exit(1);
      }
    } else {
      tpg.solver(trans_prob_graph<word_graph_node*,bool>::vogels);
      tpg.solve();
      if (verbose && tpg.nnode() > 10000) {
        std::cerr << "Solution: " << tpg.evaluate_solution() << endl;
      }
    }
 
    std::list<trans_prob_edge<word_graph_node*,bool>*>::iterator te;
    te = tpg.edges().begin();
    while (te!=tpg.edges().end()) {
      if ((*te)->flow() > 0 &&
	  (*te)->from()->data() != 0 &&
	  (*te)->to()->data() != 0) {
	for (int j=0;j<(*te)->flow();j++) {
	  real_word_graph_edge * e = new real_word_graph_edge;
	  e->from((*te)->from()->data());
	  e->to((*te)->to()->data());
	  e->seq_end(e->to()->seq_end());
	  e->length((*te)->cost());
	  new_edge(e);
	}      
      }
      ++te;
    }

    // checkpoint;

  }
  return 0;
}

long unsigned int 
word_graph::remove_trivial_nodes(CharacterProducer & cp, bool verbose) {
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
	artificial_word_graph_edge *enew = new artificial_word_graph_edge;
	enew->from(nin);
	enew->to(nout);
	std::string seq = ein->sequence(cp) + eout->sequence(cp);
	enew->sequence(seq);
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
};

class node_cmp {
  CharacterProducer & cp_;
public:
  node_cmp(CharacterProducer & cp): cp_(cp) {};
  bool operator()(word_graph_node* a, word_graph_node* b) {
    return (a->sequence(cp_) < b->sequence(cp_));
  }
};

class node_less {
private:
  CharacterProducer &cp_;
public:
  node_less(CharacterProducer &cp) : cp_(cp) {};
  bool operator()(word_graph_node * const & a, word_graph_node * const & b) const {
    // checkpoint;
    // cerr << a->sequence(cp_) << " < " << b->sequence(cp_) << " " << ((a->sequence(cp_) < b->sequence(cp_))?"TRUE":"FALSE") << endl;
    return (a->sequence(cp_) < b->sequence(cp_));
  }
};

void word_graph::sort_nodes(CharacterProducer & cp) {
  intlabelgraph<word_graph_node,word_graph_edge>::sort_nodes(node_less(cp));
}

class edge_ident {
public:
  bool operator()(word_graph_edge * const & a, word_graph_edge * const & b) const {
    if (a->from()->name() < b->from()->name()) {
      return true;
    } else if (a->from()->name() > b->from()->name()) {
      return false;
    } 
    return (a->to()->name() < b->to()->name());
  }
};

long unsigned int 
word_graph::peel_edges(FastaFile<Lazy_Header_SI> &cp, int mersize, char eos) {
  
  timestamp("peel_edges:: begin");

  clear_edge_marks();
  clear_node_marks();

  int twonodewidgets=0;
  int onenodewidgets=0;

  nodenodemap widgetmap;
  nodelist_cit nit=nodes().begin();
  while (nit != nodes().end()) {
    if ((*nit)->nin()>1) {
      if ((*nit)->nout()==1) {
	// possible widget start..
	word_graph_node *n = (*((*nit)->out().begin()))->to();
	while (n->nin()==1&&n->nout()==1) {
	  n = (*(n->out().begin()))->to();
	}
	if (n->nin()==1 && n->nout()>1) {
	  // we have a widget
	  // check neighbours
	  bool good=true;
	  edgelist_cit eit;
	  eit=(*nit)->in().begin();
	  while (eit!=(*nit)->in().end()) {
	    if (widgetmap.find((*eit)->from()) != widgetmap.end()) {
	      good=false;
	      break;
	    }
	    ++eit;
	  }
	  eit=n->out().begin();
	  while (good && eit!=n->out().end()) {
	    if (widgetmap.find((*eit)->to()) != widgetmap.end()) {
	      good=false;
	      break;
	    }
	    ++eit;
	  }
	  if (good) {
	    widgetmap.insert(make_pair(*nit,n));
	    widgetmap.insert(make_pair(n,(word_graph_node*)0));
	    twonodewidgets++;
	  }
	} 
      } else if ((*nit)->nout()>1) {
	// single node widget
	// check neighbours
	bool good=true;
	edgelist_cit eit;
	eit=(*nit)->in().begin();
	while (eit!=(*nit)->in().end()) {
	  if (widgetmap.find((*eit)->from()) != widgetmap.end()) {
	    good=false;
	    break;
	  }
	  ++eit;
	}
	eit=(*nit)->out().begin();
	while (good && eit!=(*nit)->out().end()) {
	  if (widgetmap.find((*eit)->to()) != widgetmap.end()) {
	    good=false;
	    break;
	  }
	  ++eit;
	}
	if (good) {
	  widgetmap.insert(make_pair(*nit,*nit));
	  onenodewidgets++;
	}
      }
    }
    ++nit;
  }

  long unsigned int nchanges = 0;

  checkpoint;
  cerr << "Found " << widgetmap.size() << "("<< onenodewidgets<<","<<twonodewidgets<<") widget nodes" << endl;

  timestamp("peel_edges:: begin sort");

  sort_nodes(cp);

  timestamp("peel_edges:: sort done");
  
  int usefulreads=0;

  int j = 0;
  nodelist_cit nit0=nodes().begin();
  while (1) {
    // Find node(s) corresponding to read # j
    bool exit = false;
    cp.fasta_pos(j,0);
    FILE_POSITION_TYPE pos = cp.pos();
    string buffer = cp.getstr(eos);
    // checkpoint;
    // cerr << j << " " << pos << " " << buffer << endl;
    string bufferprefix = buffer.substr(0,mersize);
    string node_buffer = (*nit0)->sequence(cp);
    // checkpoint;
    // cerr << j << " " << bufferprefix << " " << node_buffer << endl;
    while (node_buffer != bufferprefix) {
      // checkpoint;
      // cerr << j << " " << bufferprefix << " " << node_buffer << endl;
      if (node_buffer < bufferprefix) {
	++nit0;
	if (nit0 == nodes().end()) {
	  exit=true;
	  break;
	}
	node_buffer = (*nit0)->sequence(cp);
      }
      if (node_buffer > bufferprefix) {
	checkpoint;
	cerr << "No start node found for read: " << j << " " << buffer << endl;
	j++;
	if (!cp.fasta_pos(j,0)) {
	  exit=true;
	  break;
	}
	pos = cp.pos();
	buffer = cp.getstr(eos);
      }
    }
    if (exit) {
      break;
    }
    // checkpoint;
    // cerr << j << " " << bufferprefix << " " << node_buffer << endl;

    // nit0 is beginning of nodes with correct sequence
    // nit1 is first node with differnt sequence
    
    nodelist_cit nit1 = nit0;
    ++nit1;
    while (nit1 != nodes().end() && (*nit1)->sequence(cp) == node_buffer) {
      ++nit1;
    }

    // For each node that has the initial sequence corresponding to read j
    
    int npath=0;
    nodelist partialpathend;
    bool readisuseful=false;
    int p;
    nodelist_cit nit = nit0;
    while (nit != nit1) {

      p = mersize;
      word_graph_node *n = (*nit);
      ++nit;

      edgelist path;

      bool nopath=false;
      while (p < buffer.length() && buffer[p] != eos) {
	edgelist_cit eit=n->out().begin();
	while (eit != n->out().end()) {
	  if (!(*eit)->mark() && (*eit)->sequence(cp)[0] == buffer[p]) {
	    break;
	  }
	  ++eit;
	}
	if (eit == n->out().end()) {
	  nopath=true;
	  partialpathend.push_back(n);
	  break;
	} 
	path.push_back(*eit);
	p += (*eit)->length();
	n = (*eit)->to();
      }
      if (nopath) {
	continue;
      } else {
	npath++;
      }
      
      // Have the path corresponding to the read, now look for a
      // widget on the path.

//       checkpoint;
//       edgelist_cit peit=path.begin();
//       (*peit)->from()->dump(cerr,cp);
//       while (peit != path.end()) {
// 	(*peit)->to()->dump(cerr,cp);
// 	++peit;
//       }

      nodenodemap_cit wmit;
      
      int prefix_len=0;
      int edge_len=0;
      edgelist_cit eit;
      eit = path.begin();
      while (eit != path.end()) { //0
	edgelist_cit left_end,right_end;
	while (eit != path.end()) { //1
	  wmit = widgetmap.find((*eit)->to());
	  if (wmit != widgetmap.end() && wmit->second != 0) {
	    break;
	  }
	  prefix_len += (*eit)->length();
	  ++eit;
	}
	if (eit == path.end()) {
	  break;
	}
	left_end = eit;
	word_graph_node* xnd1 = wmit->second;
	edge_len += (*eit)->length();
	
	eit++;
	while (eit != path.end()) {
	  if ((*eit)->from() == xnd1) {
	    break;
	  }
	  edge_len += (*eit)->length();
	  ++eit;
	}
	if (eit == path.end()) {
	  break;
	}
	right_end = eit;
	edge_len += (*eit)->length();
	
	word_graph_node *f = (*left_end)->from();
	word_graph_node *s = (*left_end)->to();
	word_graph_node *t = (*right_end)->to();
	word_graph_node *e = (*right_end)->from();
	
	readisuseful = true;

	checkpoint;
	s->dump(cerr,cp);
	checkpoint;
	e->dump(cerr,cp);

	edgelist_cit eit1 = f->out().begin();
	while (eit1 != f->out().end()) {
	  if (((*eit1)->to() == t) && (*eit1)->mark()){ 
	    break;
	  }
	  ++eit1;
	}
	
	if (eit1 == f->out().end()) {
	  checkpoint;
	  real_count_word_graph_edge *e = new real_count_word_graph_edge;
	  e->from(f);
	  e->to(t);
	  e->length(edge_len);
	  e->seq_end(pos+mersize+prefix_len+edge_len);
	  e->count(1);
	  e->set_mark();
	  new_edge(e);

	  // checkpoint;
	  // f->dump(cerr,cp);
	  // checkpoint;
	  // t->dump(cerr,cp);

	}

	eit = left_end;
	++eit;
      }
    }

    if (npath == 0) {
      checkpoint;
      cerr << "Can't find path for read:" << j << " " << buffer << endl;
      checkpoint;
      nodelist_it nit = partialpathend.begin();
      while (nit != partialpathend.end()) {
	(*nit)->dump(cerr,cp);
	++nit;
      }
    }
    if (readisuseful) {
      usefulreads++;
    }
    
    j++;

  }

  checkpoint;
  cerr << "Useful reads: " << usefulreads << endl;

  nodenodemap_cit wit = widgetmap.begin();
  while (wit != widgetmap.end()) {
    if (wit->second == 0) {
      ++wit;
      continue;
    }
    checkpoint;
    word_graph_node *xnd = wit->first;
    edgelist_cit eit = xnd->in().begin();
    bool good = true;
    while (eit != xnd->in().end()) {
      (*eit)->from()->dump(cerr,cp);
      if ((*eit)->from()->nout_mark(true) < 1) {
	good = false;
      }
      ++eit;
    }
    xnd->dump(cerr,cp);
    word_graph_node *xnd1 = wit->second;
    if (xnd1 != xnd) {
      xnd1->dump(cerr,cp);
    }
    eit = xnd1->out().begin();
    while (eit != xnd1->out().end()) {
      (*eit)->to()->dump(cerr,cp);
      if ((*eit)->to()->nin_mark(true) < 1) {
	good = false;
      }
      ++eit;
    }
    checkpoint;
    cerr << ((good)?"GOOD":"NOT GOOD") << endl;

    set<word_graph_node*> sourcenodes;
    set<word_graph_node*>::iterator snit;

    edgelist newedges;
    
    checkpoint;
    eit = xnd->in().begin();
    while (eit != xnd->in().end()) {
      edgelist_it eit1 = (*eit)->from()->out().begin();
      while (eit1 != (*eit)->from()->out().end()) {
	if ((*eit1)->mark()) {
	  bool good = false;
	  edgelist_cit eit2=xnd1->out().begin();
	  while (eit2 != xnd1->out().end()) {
	    if ((*eit2)->to() == (*eit1)->to()) {
	      good=true;
	      break;
	    }
	    ++eit2;
	  }
	  if (good) {
	    newedges.push_back(*eit1);
	    sourcenodes.insert((*eit1)->from());
	  }
	} 
	++eit1;
      }
      ++eit;
    }

    checkpoint;
    
    eit=newedges.begin();
    while (eit != newedges.end()) {
      checkpoint;
      (*eit)->from()->dump(cerr,cp);
      (*eit)->to()->dump(cerr,cp);
      ++eit;
    }

    checkpoint;
    
    snit = sourcenodes.begin();
    while (snit != sourcenodes.end()) {
      checkpoint;
      (*snit)->dump(cerr,cp);
      ++snit;
    }

    checkpoint;
    typedef map<word_graph_node*,pair<word_graph_node*,word_graph_node*> > s2cm_t;
    s2cm_t source2clone;
    s2cm_t::iterator cmit;
    
    // for each source node, clone xnd -> xnd1 and save the begin and end.
    snit = sourcenodes.begin();
    if (good && snit != sourcenodes.end()) {
      // if good, then the first source gets the original xnd,xnd1 pair
      source2clone.insert(make_pair(*snit,make_pair(xnd,xnd1)));
      ++snit;
    }
    while (snit != sourcenodes.end()) {
      fixedlen_word_graph_node *cxnd = new fixedlen_word_graph_node;
      fixedlen_word_graph_node *c0,*c1;
      cxnd->name(new_label());
      cxnd->unset_mark();
      new_node(cxnd);
      word_graph_node *n0,*n1;
      c0 = cxnd;
      n0 = xnd;
      while (n0 != xnd1) {
	edgelist_cit eit2 = n0->out().begin();
	fixedlen_word_graph_node *c1 = new fixedlen_word_graph_node;
	c1->name(new_label());
	c1->unset_mark();
	new_node(c1);
	real_count_word_graph_edge *ce = new real_count_word_graph_edge;
	ce->from(c0);
	ce->to(c1);
	ce->unset_mark();
	ce->count(1);
	ce->length(((real_count_word_graph_edge*)(*eit2))->length());
	ce->seq_end(((real_count_word_graph_edge*)(*eit2))->seq_end());
	new_edge(ce);
	n0 = (*eit2)->to();
	c0 = c1;
      }
      word_graph_node *cxnd1 = c0;
      source2clone.insert(make_pair(*snit,make_pair(cxnd,cxnd1)));
      ++snit;
    }

    checkpoint;

    snit = sourcenodes.begin();
    while (snit != sourcenodes.end()) {
      cmit = source2clone.find((*snit));
      assert(cmit != source2clone.end());

      real_count_word_graph_edge *e;
      if (cmit->second.first != xnd) {
	real_count_word_graph_edge *cse = new real_count_word_graph_edge;
	cse->from((*snit));
	cse->to(cmit->second.first);
	cse->unset_mark();
	cse->count(1);
	e = (real_count_word_graph_edge*)find_one((*snit),xnd);
	e->set_mark();
	cse->seq_end(e->seq_end());
	cse->length(e->length());
	new_edge(cse);
      }
      ++snit;
    }

    checkpoint;
      
    // Now connect up the from and to nodes of the new edges to the
    // xnd->xnd1 clone paths...
    eit = newedges.begin();
    while (eit != newedges.end()) {
      
      cmit = source2clone.find((*eit)->from());
      assert(cmit != source2clone.end());

      real_count_word_graph_edge *e;
      if (cmit->second.first != xnd) {
	real_count_word_graph_edge *cee = new real_count_word_graph_edge;
	cee->from(cmit->second.second);
	cee->to((*eit)->to());
	cee->unset_mark();
	cee->count(1);
	e = (real_count_word_graph_edge*)find_one(xnd1,(*eit)->to());
	e->set_mark();
	cee->seq_end(e->seq_end());
	cee->length(e->length());
	new_edge(cee);
      }
      ++eit;
    }

    checkpoint;

    eit = newedges.begin();
    while (eit != newedges.end()) {
      cmit = source2clone.find((*eit)->from());
      if (cmit->second.first == xnd) {
	checkpoint;
	(*eit)->from()->dump(cerr,cp);
	checkpoint;
	xnd->dump(cerr,cp);
	if (xnd1 != xnd) {
	  checkpoint;
	  xnd1->dump(cerr,cp);
	}
	checkpoint;
	(*eit)->to()->dump(cerr,cp);
	word_graph_edge *e = (real_count_word_graph_edge*)find_one(xnd1,(*eit)->to());
	assert(e!=0);
	e->unset_mark();
      }
      ++eit;
    }
    
    checkpoint;

    if (!newedges.empty()) {
      nchanges++;
    }

    ++wit;
  }

  remove_marked_edges();
  // remove_marked_nodes();

  // sort_nodes(cp);

  checkpoint;
  dump(cerr,cp);

  return nchanges;
  
}
