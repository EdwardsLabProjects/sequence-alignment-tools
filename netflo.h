
#ifndef _NETFLO_H_
#define _NETFLO_H_

#include <iostream>
#include <fstream>

class netflo {
 private:
  /* number of nodes */
  long unsigned int nnnn;
  long unsigned int lnode;
  long unsigned int lnodp1;
  long unsigned int mnode;
  long unsigned int mnodp1;
  long unsigned int mnodp2;
  /* dimensioned to number of nodes */
  long signed int *down;
  long signed int *next;
  long signed int *level;
  long signed int *arcid;
  long signed int *flow;
  long signed int *equiv;
  long signed int *dual;
  long signed int *cat;

  /* number of edges */
  long unsigned int aaaa;
  long unsigned int larc;
  long unsigned int larcp1;
  long unsigned int marc;
  long unsigned int medge;
  /* dimensioned to number of arcs */
  long signed int *from;
  long signed int *cost;
  long signed int *capac;
  long signed int *floor;
  long signed int *arcflow_;
  long signed int *arcto_;
  long unsigned int arcnam;

  /* other variables */
  long unsigned int iter;
  long signed int mreg;
  long signed int mslk;
  long signed int mtree;
  long signed int kost0;
  long signed int kost;
  long signed int net;
  long signed int msorc;
  long signed int itdopt;
  long signed int to;
  long signed int slack;
  long signed int artif;
  long signed int dummy;
  long signed int xcess;
  long signed int price0;
  long signed int too;
  long signed int try_;
  long signed int price;
  long signed int newarc;
  long signed int newpr;
  long signed int newfrm;
  long signed int newto;
  long signed int arty;
  long signed int big;
  long signed int artyp1;
  long signed int thd;
  long signed int dw[3];
  long signed int ch[3];
  long signed int dwn;
  long signed int chg;
  long signed int theta;
  long signed int itheta;
  long signed int jtheta;
  long signed int ktheta;
  long signed int poss;
  long signed int jposs;
  long signed int frm;
  long signed int lvj;
  long signed int fm;
  long signed int lst;
  long signed int dwe;
  long signed int flw;
  long signed int aid;
  long signed int q1;
  long signed int q2;
  long signed int dir;
  long signed int ref;
  long signed int u1;
  long signed int u2;
  long signed int u3;
  long signed int u4;
  long signed int nxt;
  long signed int lsave;
  long signed int n;
  long signed int lstar;
  long signed int k440;
  long signed int ldiff;
  bool infeas;
  bool optim;
  bool dmp;
  bool ppr;
  bool pres;
  void init();
  void clear();
  long signed int isign(long signed int a,long signed int b);
  long signed int max0(long signed int a,long signed int b);
  long signed int ieor(long signed int a,long signed int b);
  long signed int iabs(long signed int a);

 public:
  netflo();
  void nodes(long unsigned int n);
  void edges(long unsigned int e);
  long unsigned int nodes();
  long unsigned int edges();
  long unsigned int arcs();
  long signed int infty();
 ~netflo();
  void netflow_input_begin();
  void netflow(long unsigned int node, long signed int netflow);
  void netflow_input_end();
  long signed int netflow(long unsigned int node);
  void arcsin_input_begin();
  void arcsin(long unsigned int node, long unsigned int narcs);
  void arcsin_input_end();
  void arc_input_begin();
  void arc(long unsigned int from, long unsigned int to, long signed int cost, 
	   long signed int ub, long signed int lb);
  void arc(long unsigned int from, long unsigned int to, long signed int cost, 
	   long signed int ub);
  void arc(long unsigned int from, long unsigned int to, long signed int cost);
  void arc_input_end();
  bool arcreal(long unsigned int arcid);
  long signed int arcflow(long unsigned int arcid);
  long unsigned int arcfrom(long unsigned int arcid);
  long unsigned int arcto(long unsigned int arcid);
  long signed int arccost(long unsigned int arcid);
  long signed int arclb(long unsigned int arcid);
  long signed int arcub(long unsigned int arcid);
  typedef enum  {
    optimal,
    infeasible,
    error
  } solution_code;
  solution_code solve();
  long signed int objective();
  long signed int iterations();
};

#endif /* _NETFLO_H_ */
