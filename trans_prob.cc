#include "trans_prob.h"

  /*   
  
  
  std::list<trans_prob_node*>::iterator tn,sit,dit;
  tn = nodes().begin();
  sit = nodes().end();
  while (tn != nodes().end()) {
    if ((*tn)->netflow() > 0) {
      sit = tn;
      break;
    } 
    ++tn;
  }
  if (sit == nodes().end()) {
    checkpoint;
    timestamp("trans_prob_graph::init_solution: No supply nodes!");
    abort();
  }

  long unsigned int nn = nnode();
  clear_node_marks();
  trans_prob_node* s;
  trans_prob_node* d;
  std::list<trans_prob_edge*>::iterator eit;
  trans_prob_edge* e;
  s = *sit;
  eit = s->out().begin();
  e = *eit;
  d = (*eit)->to();
  bool spivot = true;
  long unsigned int fl,pivotflow = 0;
  
  long unsigned int basicedges=0;

  while (1) {
    if (spivot) {
      fl = s->netflow()-pivotflow;
      if (fl > -d->netflow()) {
	fl = -d->netflow();
      }
      e->flow(fl);
      assert(!e->basic());
      e->basic(true);
      basicedges++;
      pivotflow += fl;
      checkpoint;
      cerr << "Basic edge " << basicedges << ": " << pivotflow << " ";
      s->write(cerr);
      e->write(cerr);
      d->write(cerr);
      cerr << endl;
      if (s->netflow() == pivotflow) {
	checkpoint;
	s->set_mark();
	nn--;
	spivot = false;
	pivotflow = fl;
	eit = d->in().begin();
	while (eit != d->in().end() && 
	       ((*eit)->basic() || (*eit)->from()->mark())) {
	  eit++;
	}
	assert(eit!=d->in().end());
	e = (*eit);
	s = e->from();
      } else {
	d->set_mark();
	nn--;
	eit++;
	while (eit != s->out().end()  && 
	       ((*eit)->basic() || (*eit)->to()->mark())) {
	  eit++;
	}
	assert(eit!=s->out().end());
	e = (*eit);
	d = e->to();
      }
    } else {
      fl = (-d->netflow())-pivotflow;
      if (fl > s->netflow()) {
	fl = s->netflow();
      }
      e->flow(fl);
      assert(!e->basic());
      e->basic(true);
      basicedges++;
      pivotflow += fl;
      checkpoint;
      cerr << "Basic edge " << basicedges << ":";
      s->write(cerr);
      e->write(cerr);
      d->write(cerr);
      cerr << " " << pivotflow << endl;
      if ((-d->netflow()) == pivotflow) {
	checkpoint;
	d->set_mark();
	nn--;
	spivot = true;
	pivotflow = fl;
	eit = s->out().begin();
	while (eit != s->out().end() && 
	       ((*eit)->basic() || (*eit)->to()->mark())) {
	  eit++;
	}
	assert(eit!=s->out().end());
	e = (*eit);
	d = e->to();
      } else {
	s->set_mark();
	nn--;
	eit++;
	while (eit != d->in().end() && 
	       ((*eit)->basic() || (*eit)->from()->mark())) {
	  eit++;
	}
	assert(eit!=d->in().end());
	e = (*eit);
	s = e->from();
      }      
    }
  }

  tn = nodes().begin();
  while (tn != nodes().end()) {
    if (!(*tn)->mark()) {
      checkpoint;
      cerr << "Node ";
      (*tn)->write(cerr);
      cerr << " not marked!" << endl;
    } 
    ++tn;
  }
  */
}

