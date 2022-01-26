
#ifndef _TRANS_PROB_SOVLER_H_
#define _TRANS_PROB_SOVLER_H_

#include "graph.h"
#include "netflo.h"
#include "CS2.h"
#include <sstream>

template <class ND, class ED>
class trans_prob_edge;

template <class ND, class ED>
class trans_prob_node : public labelnode<trans_prob_edge<ND,ED> > {
  long signed int netflow_;
  long signed int potential_;
  ND data_;
 public:
  trans_prob_node() : netflow_(0), potential_(0) {};
  trans_prob_node(ND const & t) : netflow_(0), data_(t),  potential_(0) {};
  void netflow(long signed int f) {
    netflow_ = f;
  }
  long signed int netflow() const {
    return netflow_;
  }
  void potential(long signed int i) {
    potential_ = i;
  }
  long signed int potential() const {
    return potential_;
  }
  ND const & data() const {
    return data_;
  }
  void data(ND const & t) {
    data_ = t;
  }
  void write(ostream &os) const {
    os << this->name();
  }
};

template <class ND, class ED>
class trans_prob_edge : public labeledge<trans_prob_node<ND,ED> > {
  long signed int cost_;
  long signed int flow_;
  long signed int reduced_cost_;
  bool basic_;
  ED data_;
public:
  trans_prob_edge() : cost_(0), flow_(0), reduced_cost_(0), basic_(false) {};
  trans_prob_edge(ED const & t) : cost_(0), flow_(0), data_(t), reduced_cost_(0), basic_(false) {};
  void cost(long signed int c) {
    cost_ = c;
  }
  long signed int cost() const {
    return cost_;
  }
  void flow(long signed int f) {
    flow_ = f;
  }
  long signed int flow() const {
    return flow_;
  }
  void basic(bool b) {
    basic_ = b;
  }
  bool basic() const {
    return basic_;
  }
  void reduced_cost(long signed int i) {
    reduced_cost_ = i;
  }
  long signed int reduced_cost() const {
    return reduced_cost_;
  }
  ED const & data() const {
    return data_;
  }
  void data(ED const & t) {
    data_ = t;
  }
};

template <class ND,class ED>
class trans_prob_graph : public labelgraph<trans_prob_node<ND,ED>,trans_prob_edge<ND,ED> > {
public:
  typedef enum {
    netflo,
    cs2,
    greedy,
    vogels,
    netsimplex
  } solver_id;
 private:
  static bool vogels_lt(trans_prob_edge<ND,ED>* const & a, 
			trans_prob_edge<ND,ED>* const & b);
  static bool min_cost_edge_lt(trans_prob_edge<ND,ED>* const & a, 
			       trans_prob_edge<ND,ED>* const & b);
  long signed int compute_reduced_cost();
  void compute_node_potentials();
  bool iterate() const;
  trans_prob_edge<ND,ED>* select_entering_edge() const;
  void add_edge_to_basis(trans_prob_edge<ND,ED>*);
  void solve_my(long signed int iter=-1);
  void solve_nf();
  void solve_cs2();
  solver_id solver_;
 public:
  trans_prob_graph() { };
  void heuristic(bool vogels);
  void check_solution() const;
  long signed int evaluate_solution() const;
  void solver(solver_id const & s) {
    solver_ = s;
  } 
  solver_id const & solver() const {
    return solver_;
  }
  void solve();
  void write_dimacs_input(const char *filename);
};

template <class ND,class ED>
long signed int trans_prob_graph<ND,ED>::evaluate_solution() const {
  long signed int val=0;
  typename std::list<trans_prob_edge<ND,ED>*>::const_iterator eit;
  eit = this->edges().begin();
  while (eit != this->edges().end()) {
    val += (*eit)->cost()*(*eit)->flow();
    ++eit;
  }
  return val;
}

template <class ND,class ED>
void trans_prob_graph<ND,ED>::check_solution() const {
  typename std::list<trans_prob_node<ND,ED>*>::const_iterator tn;
  typename std::list<trans_prob_edge<ND,ED>*>::const_iterator eit;
  tn = this->nodes().begin();
  while (tn != this->nodes().end()) {
    trans_prob_node<ND,ED> *n = *tn;
    long signed int flowin=0;
    eit=n->in().begin();
    while (eit!=n->in().end()) {
      assert((*eit)->flow() >= 0);
      flowin += (*eit)->flow();
      ++eit;
    }
    long signed int flowout=0;
    eit=n->out().begin();
    while (eit!=n->out().end()) {
      assert((*eit)->flow() >= 0);
      flowout += (*eit)->flow();
      ++eit;
    }
    if (n->netflow() != (flowout - flowin)) {
      n->write(std::cerr);
      std::cerr << " " << n->netflow()
	   << " " << flowin 
	   << " " << flowout << endl;
      assert(n->netflow() == (flowout - flowin));
    }
    ++tn;
  }
};

template <class ND,class ED>
void trans_prob_graph<ND,ED>::write_dimacs_input(const char *filename) {
  ofstream ofs(filename);
  this->relabel();
  ofs << "p min " << this->nnode() << " " << this->nedge() << endl;
  typename std::list<trans_prob_node<ND,ED>*>::const_iterator tn;
  typename std::list<trans_prob_edge<ND,ED>*>::const_iterator eit;
  tn = this->nodes().begin();
  while (tn != this->nodes().end()) {
    trans_prob_node<ND,ED> *n = *tn;
    if (n->netflow() != 0) {
      ofs << "n " << n->name()+1 << " " << n->netflow() << endl;
    }
    ++tn;
  }
  eit = this->edges().begin();
  while (eit != this->edges().end()) {
    trans_prob_edge<ND,ED> *e = *eit;
    if (e->from()->nin() == 0) {
      ofs << "a" 
	  << " " << e->from()->name()+1
	  << " " << e->to()->name()+1
	  << " " << 0
	  << " " << e->from()->netflow()
	  << " " << e->cost()
	  << endl;
    } else if (e->to()->nout() == 0) {
      ofs << "a" 
	  << " " << e->from()->name()+1
	  << " " << e->to()->name()+1
	  << " " << 0
	  << " " << - e->to()->netflow()
	  << " " << e->cost()
	  << endl;
    } else {
      ofs << "a" 
	  << " " << e->from()->name()+1
	  << " " << e->to()->name()+1
	  << " " << 0
	  << " " << 999999
	  << " " << e->cost()
	  << endl;
    }
    ++eit;
  }
  ofs.close();
};

template <class ND, class ED>
void trans_prob_graph<ND,ED>::solve() {
  switch (solver_) {
  case netflo:
    solve_nf();
    break;
  case cs2:
    solve_cs2();
    break;
  case greedy:
    heuristic(false);
    break;
  case vogels:
    heuristic(true);
    break;
  case netsimplex:
    solve_my();
    break;
  default:
    std::cerr << "trans_prob_graph::solve: Invalid solver setting\n";
    assert(0);
  }
  check_solution();
}

template <class ND, class ED>
void trans_prob_graph<ND,ED>::solve_my(long signed int maxit) {
  long signed int val;
  long unsigned int its=0;
  heuristic(true);
  check_solution();
  val = evaluate_solution();
  // checkpoint;
  // std::cerr << "Solution: " << val << endl;
  compute_reduced_cost();
  trans_prob_edge<ND,ED> *e;
  while (iterate() && (maxit < 0 || its <= maxit)) {
    e = select_entering_edge();
    add_edge_to_basis(e);
    check_solution();
    val = evaluate_solution();
    // checkpoint;
    // std::cerr << "Solution: " << val << endl;
    compute_reduced_cost();
    its++;
  }
}

template <class ND, class ED>
trans_prob_edge<ND,ED> *trans_prob_graph<ND,ED>::select_entering_edge() const {
  typename std::list<trans_prob_edge<ND,ED>*>::const_iterator eit;
  trans_prob_edge<ND,ED>* best_edge=0;
  long signed int best_score=0,score=0;
  long signed int potential_flow=0;
  eit = this->edges().begin();
  while (eit!=this->edges().end()) {
    if (!(*eit)->basic()) {
      if ((score=(*eit)->reduced_cost()) < 0) {
	potential_flow = (*eit)->from()->netflow();
	if (potential_flow > -(*eit)->to()->netflow()) {
	  potential_flow = -(*eit)->to()->netflow();
	}
	score *= potential_flow;
	if (score < best_score) {
	  best_edge = (*eit);
	  best_score = score;
	}
      }
    }
    ++eit;
  }
  assert(best_edge != 0);
  return best_edge;
}


template <class ND, class ED>
bool trans_prob_graph<ND,ED>::iterate() const {
  typename std::list<trans_prob_edge<ND,ED>*>::const_iterator eit;
  eit = this->edges().begin();
  while (eit != this->edges().end()) {
    if ((*eit)->reduced_cost() < 0) {
      return true;
    }
    ++eit;
  }
  return false;
}

template <class ND, class ED>
void trans_prob_graph<ND,ED>::compute_node_potentials() {

  typename std::list<trans_prob_node<ND,ED>*>::iterator nit;
  typename std::list<trans_prob_edge<ND,ED>*>::iterator eit;
  typename std::list<trans_prob_node<ND,ED>*> frontier;

  this->clear_node_marks();
  trans_prob_node<ND,ED>* n = this->nodes().front();
  n->potential(0);
  n->set_mark();
  frontier.push_back(n);

  while (!frontier.empty()) {
    n = frontier.front();
    frontier.pop_front();
    eit = n->out().begin();
    while (eit != n->out().end()) {
      trans_prob_edge<ND,ED> *e = (*eit);
      if (e->basic() && !e->to()->mark()) {
	e->to()->potential(n->potential() + e->cost());
	e->to()->set_mark();
	frontier.push_back(e->to());
      }
      ++eit;
    }
    eit = n->in().begin();
    while (eit != n->in().end()) {
      trans_prob_edge<ND,ED> *e = (*eit);
      if (e->basic() && !e->from()->mark()) {
	e->from()->potential(n->potential() - e->cost());
	e->from()->set_mark();
	frontier.push_back(e->from());
      }
      ++eit;
    }
  }
  
}

template <class ND, class ED>
long signed int trans_prob_graph<ND,ED>::compute_reduced_cost() {
  typename std::list<trans_prob_node<ND,ED>*>::iterator nit;
  typename std::list<trans_prob_node<ND,ED>*> frontier;

  long signed int minrc=0;

  compute_node_potentials();

  typename std::list<trans_prob_edge<ND,ED>*>::iterator eit = this->edges().begin();
  while (eit != this->edges().end()) {
    trans_prob_edge<ND,ED> *e = (*eit);
    e->reduced_cost(e->cost() - (e->to()->potential() - e->from()->potential()));
    assert(!e->basic() || e->reduced_cost() == 0);
    if (minrc > e->reduced_cost()) {
      minrc = e->reduced_cost();
    }
    ++eit;
  }
  // checkpoint;
  // std::cerr << "Minimum reduced cost: " << minrc << endl;
  return minrc;
}

template <class ND, class ED>
struct tpgpathsearchelt {
  trans_prob_node<ND,ED> * n;
  typename std::list<trans_prob_edge<ND,ED>*> el;
  tpgpathsearchelt(trans_prob_node<ND,ED> * const n0=0) : n(n0) {}
  tpgpathsearchelt(trans_prob_node<ND,ED> * const n0, 
		   typename std::list<trans_prob_edge<ND,ED>*> const & el0)
    : n(n0), el(el0) {};
};

template <class ND, class ED>
void trans_prob_graph<ND,ED>::add_edge_to_basis(trans_prob_edge<ND,ED> *entering_edge) {
  // We need to find a path from entering_edge's to node to its from
  // node, on baisic edges only. Directionality of edges can be ignored. 

  checkpoint;

  typename std::list<trans_prob_edge<ND,ED>*>::iterator eit;
  std::list<trans_prob_edge<ND,ED>*> cycle;

  this->clear_node_marks();

  //   checkpoint;
  //   std::cerr << "Entering edge: ";
  //   entering_edge->from()->write(std::cerr);
  //   entering_edge->write(std::cerr);
  //   entering_edge->to()->write(std::cerr);
  //   std::cerr << endl;
  
  std::list<tpgpathsearchelt<ND,ED> > queue;
  queue.push_back(tpgpathsearchelt<ND,ED>(entering_edge->to()));
  while (!queue.empty()) {
    // checkpoint;
    // std::cerr << "Queue length: " << queue.size() << endl;
    tpgpathsearchelt<ND,ED> se = queue.front();
    queue.pop_front();
    trans_prob_node<ND,ED> *n = se.n;
    n->set_mark();
    // checkpoint;
    // std::cerr << "Pop node: ";
    // n->write(std::cerr) << endl;
    // std::cerr << "Path length: " << se.el.size() << endl;
    if (n == entering_edge->from()) {
      cycle = se.el;
      break;
    } else {
      eit = n->out().begin();
      while (eit!=n->out().end()) {
	if ((*eit)->basic() && !(*eit)->to()->mark()) {
	  tpgpathsearchelt<ND,ED> se1((*eit)->to(),se.el);
	  se1.el.push_back(*eit);
	  queue.push_back(se1);
	}
	++eit;
      }
      eit = n->in().begin();
      while (eit!=n->in().end()) {
	if ((*eit)->basic() && !(*eit)->from()->mark()) {
	  tpgpathsearchelt<ND,ED> se1((*eit)->from(),se.el);
	  se1.el.push_back(*eit);
	  queue.push_back(se1);
	}
	++eit;
      }
    }
  }
  
  cycle.push_front(entering_edge);

//   checkpoint;
//   std::cerr << "Cycle length: " << cycle.size() << endl;
//   assert(cycle.size() >= 2);
  
  long unsigned int cyclepos=0;
  long signed int maxaug=MAXINT;
  trans_prob_edge<ND,ED> *leaving_edge=0;
  eit = cycle.begin();
  while (eit!=cycle.end()) {
    if (cyclepos%2==1) {
      if (maxaug > (*eit)->flow()) {
	maxaug = (*eit)->flow();
	leaving_edge = (*eit);
      }
    }
    ++cyclepos;
    ++eit;
  }
  
  
//   checkpoint;
//   std::cerr << "Leaving edge: ";
//   leaving_edge->from()->write(std::cerr);
//   leaving_edge->write(std::cerr);
//   leaving_edge->to()->write(std::cerr);
//   std::cerr << endl;

//   checkpoint;
//   std::cerr << "Cycle: " << endl;
//   eit = cycle.begin();
//   while (eit!=cycle.end()) {
//     std::cerr << "      ";
//     (*eit)->from()->write(std::cerr);
//     (*eit)->write(std::cerr);
//     (*eit)->to()->write(std::cerr);
//     std::cerr << endl;
//     ++eit;
//   }
  
//   checkpoint;
//   std::cerr << "Maximum augmenting flow: " << maxaug << endl;
  assert(maxaug>=0);
  assert(leaving_edge!=0);

  eit = cycle.begin();
  while (eit!=cycle.end()) {
    if (cyclepos%2==0) {
      (*eit)->flow((*eit)->flow()+maxaug);
    } else {
      (*eit)->flow((*eit)->flow()-maxaug);
    }
    ++cyclepos;
    ++eit;
  }
 
  entering_edge->basic(true);
  leaving_edge->basic(false);

  // update potentials
  
};

template <class ND,class ED>
void trans_prob_graph<ND,ED>::solve_nf() {

  class netflo *nf;
  this->relabel();

  // write_dimacs_input(".input0.min");

  typename std::list<trans_prob_node<ND,ED>*>::iterator nit;
  typename std::list<trans_prob_edge<ND,ED>*>::iterator eit;

  typedef std::map<std::pair<long unsigned int,long unsigned int>,trans_prob_edge<ND,ED>*> edgemap;
  edgemap nn2e;
  typename edgemap::iterator meit;
  
  // Insert best edge between two nodes into the map...
  eit = this->edges().begin();
  while (eit!=this->edges().end()) {
    if ((meit=nn2e.find(make_pair((*eit)->from()->name()+1,(*eit)->to()->name()+1)))!=nn2e.end()) {
      if ((*meit).second->cost() > (*eit)->cost()) {
	(*meit).second = (*eit);
      }
    } else {
      nn2e.insert(make_pair(make_pair((*eit)->from()->name()+1,(*eit)->to()->name()+1),(*eit)));
    }
    ++eit;
  }

  // Remove any edge that isn't in the map
  eit = this->edges().begin();
  while (eit!=this->edges().end()) {
    meit=nn2e.find(make_pair((*eit)->from()->name()+1,(*eit)->to()->name()+1));
    if (meit != nn2e.end() && (*meit).second != (*eit)) {
      eit = erase(eit);
    } else {
      eit++;
    }
  }

  nf = new class netflo;
  
  nf->nodes(this->nnode());
  nf->edges(this->nedge());
  nf->netflow_input_begin();
  nit = this->nodes().begin();
  while (nit != this->nodes().end()) {
    if ((*nit)->netflow() != 0) {
      nf->netflow((*nit)->name()+1,(*nit)->netflow());
    }
    ++nit;
  }
  nf->netflow_input_end();

  nf->arcsin_input_begin();
  nit = this->nodes().begin();
  while (nit != this->nodes().end()) {
    nf->arcsin((*nit)->name()+1,(*nit)->nin());
    ++nit;
  }
  nf->arcsin_input_end();

  nf->arc_input_begin();
  nit = this->nodes().begin();
  while (nit != this->nodes().end()) {
    eit = (*nit)->in().begin();
    while (eit != (*nit)->in().end()) {
      nf->arc((*eit)->from()->name()+1,
	     (*eit)->to()->name()+1,
	     (*eit)->cost());
      ++eit;
    }
    ++nit;
  }
  nf->arc_input_end();
  netflo::solution_code sc;
  sc = nf->solve();
  assert(sc == netflo::optimal);
  // checkpoint;
  // std::cerr << "Netflo:: Objective: " << nf->objective() << endl;
  // std::cerr << "Netflo:: Iterations: " << nf->iterations() << endl;
  // checkpoint;
  
  for (long unsigned int i=1;i<=nf->arcs();i++) {
    trans_prob_edge<ND,ED>* e=0;
    meit = nn2e.find(make_pair(nf->arcfrom(i),nf->arcto(i)));
    if (meit != nn2e.end()) {
      (*meit).second->flow(nf->arcflow(i));
    }
  }
  // checkpoint;
  check_solution();
};

template <class ND,class ED>
void trans_prob_graph<ND,ED>::solve_cs2() {

  this->relabel();

  CS2* cs2=0;
  
  // checkpoint;

  typename std::list<trans_prob_node<ND,ED>*>::iterator nit;
  typename std::list<trans_prob_edge<ND,ED>*>::iterator eit;

  typedef std::map<std::pair<long unsigned int,long unsigned int>,trans_prob_edge<ND,ED>*> edgemap;
  edgemap nn2e;
  typename edgemap::iterator meit;
  
  eit = this->edges().begin();
  while (eit!=this->edges().end()) {
    if ((meit=nn2e.find(make_pair((*eit)->from()->name()+1,(*eit)->to()->name()+1)))!=nn2e.end()) {
      if ((*meit).second->cost() > (*eit)->cost()) {
	(*meit).second = (*eit);
      }
    } else {
      nn2e.insert(make_pair(make_pair((*eit)->from()->name()+1,(*eit)->to()->name()+1),(*eit)));
    }
    ++eit;
  }

  eit = this->edges().begin();
  while (eit!=this->edges().end()) {
    meit=nn2e.find(make_pair((*eit)->from()->name()+1,(*eit)->to()->name()+1));
    if (meit != nn2e.end() && (*meit).second != (*eit)) {
      eit = erase(eit);
    } else {
      eit++;
    }
  }

  // checkpoint;

  cs2 = new CS2(this->nnode(),this->nedge());

  cs2->netflow_input_begin();
  nit = this->nodes().begin();
  while (nit != this->nodes().end()) {
    if ((*nit)->netflow() != 0) {
      cs2->netflow((*nit)->name()+1,(*nit)->netflow());
    }
    ++nit;
  }
  cs2->netflow_input_end();

  // checkpoint;
  
  cs2->arc_input_begin();
  
  // checkpoint;
  
  nit = this->nodes().begin();
  while (nit != this->nodes().end()) {
    eit = (*nit)->in().begin();
    while (eit != (*nit)->in().end()) {
      if ((*eit)->from()->nin() == 0) {
	cs2->arc((*eit)->from()->name()+1,
		 (*eit)->to()->name()+1,
		 0,
		 (*eit)->from()->netflow(), 
		 (*eit)->cost());	
      } else if ((*eit)->to()->nout() == 0) {
	cs2->arc((*eit)->from()->name()+1,
		 (*eit)->to()->name()+1,
		 0,
		 -(*eit)->to()->netflow(), 
		 (*eit)->cost());	
      } else {
	cs2->arc((*eit)->from()->name()+1,
		 (*eit)->to()->name()+1,
		 0,
		 -1, 
		 (*eit)->cost());
      }
      ++eit;
    }
    ++nit;
  }
  // checkpoint;
  cs2->arc_input_end();

  // checkpoint;

  CS2::solution_code sc;
  sc = cs2->solve();
  assert(sc == CS2::optimal);

  // checkpoint;
  // std::cerr << "CS2:: Objective: " << cs2->objective() << endl;
  // checkpoint;
  
  //  checkpoint;

  for (long unsigned int i=1;i<=cs2->narc();i++) {
    trans_prob_edge<ND,ED>* e=0;
    meit = nn2e.find(make_pair(cs2->arcfrom(i),cs2->arcto(i)));
    if (meit != nn2e.end()) {
      (*meit).second->flow(cs2->arcflow(i));
    }
  }

  // checkpoint;

  delete cs2;
};

template <class ND, class ED>
bool trans_prob_graph<ND,ED>::min_cost_edge_lt(trans_prob_edge<ND,ED>* const & a, 
		      trans_prob_edge<ND,ED>* const & b) {
  long signed int d(a->cost() - b->cost());
  return (d<0);
}

template <class ND, class ED>
bool trans_prob_graph<ND,ED>::vogels_lt(trans_prob_edge<ND,ED>* const & a, 
	       trans_prob_edge<ND,ED>* const & b) {
  if (a->reduced_cost() > b->reduced_cost()) {
    return true;
  } else if (a->reduced_cost() < b->reduced_cost()) {
    return false;
  } else if (a->cost() < b->cost()) {
    return true;
  } 
  return false;
}

template <class ND, class ED>
void trans_prob_graph<ND,ED>::heuristic(bool vogels) {

  std::vector<trans_prob_edge<ND,ED>*> sorted_edges;
  sorted_edges.reserve(this->edges().size());
  
  typename std::list<trans_prob_edge<ND,ED>*>::iterator eit;
  eit = this->edges().begin();
  while (eit!=this->edges().end()) {
    sorted_edges.push_back(*eit);
    ++eit;
  }

  typename std::list<trans_prob_node<ND,ED>*>::iterator nit;
  trans_prob_node<ND,ED> *n;
  trans_prob_edge<ND,ED> *e;
  
  if (vogels) {
    // checkpoint;
    // std::cerr << "Vogel's initial solution heuristic." << endl;
    nit = this->nodes().begin();
    while (nit!=this->nodes().end()) {
      n = (*nit);
      long signed int mincost=MAXINT,nextmincost=MAXINT;
      if (n->netflow() > 0) {
	eit = n->out().begin();
	while (eit != n->out().end()) {
	  e = (*eit);
	  if (e->cost() < mincost) {
	    mincost = e->cost();
	    nextmincost = mincost;
	  } else if (e->cost() < nextmincost) {
	    nextmincost = e->cost();
	  }
	  ++eit;
	}
      } else {
	eit = n->in().begin();
	while (eit != n->in().end()) {
	  e = (*eit);
	  if (e->cost() < mincost) {
	    mincost = e->cost();
	    nextmincost = mincost;
	  } else if (e->cost() < nextmincost) {
	    nextmincost = e->cost();
	  }
	  ++eit;
	}
      }
      (*nit)->potential(nextmincost - mincost);
      ++nit;
    }
    eit = this->edges().begin();
    while (eit != this->edges().end()) {
      (*eit)->reduced_cost((*eit)->from()->potential());
      if ((*eit)->reduced_cost() < (*eit)->to()->potential()) {
	(*eit)->reduced_cost((*eit)->to()->potential());	
      }
      ++eit;
    }
    sort(sorted_edges.begin(),sorted_edges.end(),vogels_lt);
  } else {
    // checkpoint;
    // std::cerr << "Min cost initial solution heuristic." << endl;
    sort(sorted_edges.begin(),sorted_edges.end(),min_cost_edge_lt);    
  }

  // checkpoint;
  
  this->clear_edge_marks();
  this->clear_node_marks();

  nit = this->nodes().begin();
  while (nit!=this->nodes().end()) {
    (*nit)->potential(0);
    ++nit;
  }
  
  long unsigned int maxcost=0;
  trans_prob_node<ND,ED> *tdummy,*fdummy;
  eit = this->edges().begin();
  while (eit != this->edges().end()) {
    if ((*eit)->from()->data() == 0 &&
	(*eit)->to()->data() == 0) {
      tdummy = (*eit)->to();
      fdummy = (*eit)->from();
    }
    ++eit;
  }

  // checkpoint;
  // sdummy->write(std::cerr);
  // std::cerr << endl;
  // ddummy->write(std::cerr);
  // std::cerr << endl;

  long unsigned int nbasic=0;
  
  typename std::vector<trans_prob_edge<ND,ED>*>::iterator seit;
  seit = sorted_edges.begin();
  while (seit != sorted_edges.end()) {
    if ((*seit)->from()->mark() ||
	(*seit)->to()->mark() ||
	(*seit)->from()->data() == 0 ||
	(*seit)->to()->data() == 0) {
      ++seit;
      continue;
    }
    long unsigned int potflow;
    potflow = (*seit)->from()->netflow() - (*seit)->from()->potential();
    if (potflow > (*seit)->to()->potential() - (*seit)->to()->netflow()) {
      potflow = (*seit)->to()->potential() - (*seit)->to()->netflow();
      (*seit)->to()->set_mark();
    } else {
      (*seit)->from()->set_mark();
    }
    (*seit)->basic(true);
    nbasic++;
    (*seit)->flow(potflow);
    (*seit)->to()->potential((*seit)->to()->potential()-potflow);
    (*seit)->from()->potential((*seit)->from()->potential()+potflow);
    ++seit;
  }

  // checkpoint;
  // std::cerr << "Number of basic edges: " << nbasic << endl;

  this->clear_node_marks();
  std::list<trans_prob_node<ND,ED>*> queue;
  trans_prob_node<ND,ED>* unsat_node=0;

  nit = this->nodes().begin();
  while (nit!=this->nodes().end()) {
    if (!(*nit)->mark() && 
	(*nit)->data()!=0 &&
	(*nit)->data()!=0) {
      // checkpoint;
      // std::cerr << "Component:" << endl;
      queue.clear();
      unsat_node=0;
      queue.push_back(*nit);
      (*nit)->set_mark();
      while (!queue.empty()) {
	trans_prob_node<ND,ED> *n = queue.front();
	queue.pop_front();
	// n->write(std::cerr) << endl;
	if (n->netflow()!=n->potential()) {
	  assert(unsat_node==0);
	  unsat_node = n;
	}
	eit=n->out().begin();
	while (eit!=n->out().end()) {
	  if ((*eit)->basic() && !(*eit)->to()->mark()) {
	    (*eit)->to()->set_mark();
	    queue.push_back((*eit)->to());
	  }
	  ++eit;
	}
	eit=n->in().begin();
	while (eit!=n->in().end()) {
	  if ((*eit)->basic() && !(*eit)->from()->mark()) {
	    (*eit)->from()->set_mark();
	    queue.push_back((*eit)->from());
	  }
	  ++eit;
	}
      }
      if (unsat_node==0) {
	// unsat_node = *nit;
	++nit;
	continue;
      }
      // checkpoint;
      // std::cerr << "Unsaturated node: " << endl;
      // unsat_node->write(std::cerr) << endl;
      long unsigned int potflow;
      if (unsat_node->netflow()>0) {
	eit=unsat_node->out().begin();
	while (eit!=unsat_node->out().end()) {
	  if ((*eit)->to() == fdummy) {
	    potflow = unsat_node->netflow()-unsat_node->potential();
	    (*eit)->basic(true);
	    nbasic++;
	    (*eit)->flow(potflow);
	    (*eit)->to()->potential((*eit)->to()->potential()-potflow);
	    (*eit)->from()->potential((*eit)->from()->potential()+potflow);
	    break;
	  }
	  ++eit;
	}
	assert(eit!=unsat_node->out().end());
      } else {
	eit=unsat_node->in().begin();
	while (eit!=unsat_node->in().end()) {
	  if ((*eit)->from() == tdummy) {
	    potflow = unsat_node->potential()-unsat_node->netflow();
	    (*eit)->basic(true);
	    nbasic++;
	    (*eit)->flow(potflow);
	    (*eit)->to()->potential((*eit)->to()->potential()-potflow);
	    (*eit)->from()->potential((*eit)->from()->potential()+potflow);
	    break;
	  }
	  ++eit;
	}
	assert(eit!=unsat_node->in().end());
      }
    }
    ++nit;
  }

  // checkpoint;
  // std::cerr << "Number of basic edges: " << nbasic << endl;

  long unsigned int potflow;
  potflow = fdummy->netflow() - fdummy->potential();
  assert(potflow == (tdummy->potential() - tdummy->netflow()));

  e = fdummy->out().front();
  assert(e->to() == tdummy);
  
  e->basic(true);
  nbasic++;
  e->flow(potflow);

  // checkpoint;
  // std::cerr << "Number of basic edges: " << nbasic << endl;

  check_solution();
}

#endif
