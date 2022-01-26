#include "netflo.h"
#include <assert.h>
#include <stdio.h>

#include <iostream>
#include <ios>

long signed int netflo::isign(long signed int a, long signed int b) {
  long signed int x;
  x = (a >= 0 ? a : - a);
  return ( b >= 0 ? x : -x);
}
long signed int netflo::max0(long signed int a, long signed int b) {
  if (a > b) {
    return a;
  } else {
    return b;
  }
}
long signed int netflo::iabs(long signed int a) {
  if (a > 0) {
    return a;
  } else {
    return -a;
  }
}
long signed int netflo::ieor(long signed int a, long signed int b) {
  return (a ^ b);
}
netflo::netflo() {};

void netflo::nodes(long unsigned int n) {
  /* if (outformat == 1) {
    ofs_.width(10);
    ofs_ << n << endl;
    } */
  mnode = n;
}

long unsigned int netflo::nodes() {
  return mnode;
}

void netflo::edges(long unsigned int e) {
  medge = e;
}

long unsigned int netflo::edges() {
  return medge;
}

long unsigned int netflo::arcs() {
  return marc;
}

netflo::~netflo() {
  clear();
}

void netflo::init() {
  // fprintf(stderr,"Checkpoint %s:%d\n",__FILE__,__LINE__);
  itdopt = 1;
  arty = 1;
  big = 100000000;
  u1 = u2 = u3 = u4 = 5;
  nnnn=mnode+3;
  aaaa=medge+2*mnode+2;

  equiv = new long signed int[nnnn+1];
  cat = dual = equiv;
  down = new long signed int[nnnn+1];
  flow = new long signed int[nnnn+1];
  next = new long signed int[nnnn+1];
  arcid = new long signed int[nnnn+1];
  level = new long signed int[nnnn+1];
  from = new long signed int[aaaa+1];
  cost = new long signed int[aaaa+1];
  capac = new long signed int[aaaa+1];
  floor = new long signed int[aaaa+1];
  arcflow_ = new long signed int[aaaa+1];
  arcto_ = new long signed int[aaaa+1];

  artyp1 = arty+1;
  lnode = nnnn-2;
  larc = aaaa-2;
  larcp1 = larc+1;
  lnodp1 = lnode+1;
  ppr = true;
  pres = false;
  mnodp1 = mnode+1;
  mnodp2 = mnode+2;
  assert(mnodp1<=lnode);
  dual[mnodp1]=0;
  for (long unsigned int j10=1;j10<=mnodp1;j10++) {
    down[j10] =0;
    next[j10] =0;
    level[j10] =0;
    arcid[j10] =0;
    flow[j10] =0;
  }
  from[arty] = mnodp1;
  cost[arty] = big;
  capac[arty] = 0;
  floor[arty] = 0;
  net = 0;
  msorc = 0;
  marc = 0;
  next[mnodp1] = mnodp1;
  down[mnodp1] = mnodp1;
  // fprintf(stderr,"Checkpoint %s:%d\n",__FILE__,__LINE__);
}
void 
netflo::clear() {
  delete [] equiv;
  delete [] down;
  delete [] flow;
  delete [] next;
  delete [] arcid;
  delete [] level;
  delete [] from;
  delete [] cost;
  delete [] capac;
  delete [] floor; 
  delete [] arcflow_;
}
void netflo::netflow_input_begin() {
  init();
  // if (outformat==2) {
  // ofs_ << "p min " << mnode << " " << medge << endl;
  // }
}
void netflo::netflow(long unsigned int i, long signed int j) {
  // fprintf(stderr,"Checkpoint %s:%d\n",__FILE__,__LINE__);
  // fprintf(stderr,"netflow: %ld %ld\n",i,j);
  assert(i >= 1 && i <= mnode);
  assert(flow[i] == 0);
  /*  if (outformat == 1) {
    if (j != 0) {
      ofs_.width(10);
      ofs_ << i;
      ofs_.width(10);
      ofs_ << j << endl;
    }
  } else if (outformat == 2) {
    if (j != 0) {
      ofs_ << "n " << i << " " << j << endl;
    }
  }
  */
  flow[i] = j;
  net = net+j;
  if (j > 0) {
    msorc = msorc+1;
    level[i] = j;
    next[i] = next[mnodp1];
    next[mnodp1] = i;
  }
  // fprintf(stderr,"Checkpoint %s:%d\n",__FILE__,__LINE__);
}
void netflo::netflow_input_end() {
  assert(net>=0);
  /* if (outformat == 1) {
    ofs_.width(10);
    ofs_ << 0;
    ofs_.width(10);
    ofs_ << 0 << endl;
    }*/
}
void netflo::arcsin_input_begin() {};
void netflo::arcsin(long unsigned int i, long unsigned int c) {
  // fprintf(stderr,"Checkpoint %s:%d\n",__FILE__,__LINE__);
  cat[i] = c;
  // fprintf(stderr,"Cat[%ld]=%ld\n",i,cat[i]);
  // fprintf(stderr,"Checkpoint %s:%d\n",__FILE__,__LINE__);
}
void netflo::arcsin_input_end() {
  // fprintf(stderr,"Checkpoint %s:%d\n",__FILE__,__LINE__);
  /* if (outformat == 1) {
    for (long unsigned int myi=1;myi<=((mnode/8)+1)*8;myi++) {
      if (myi > mnode) {
	ofs_.width(10);
	ofs_ << 0;
      } else {
	ofs_.width(10);
	ofs_ << cat[myi];
      }
      if (myi % 8 == 0) {
	ofs_ << endl;
      } 
    }
  } 
  */
  long signed int i = 1;
  long signed int j = arty;
  for (long unsigned int j80 = 1; j80<=mnode; j80++) {
    i = -i;
    long signed int k70 = max0(1,cat[j80]);
    // fprintf(stderr,"k70=%ld\n",k70);
    assert(j+k70 <= larc);
    cat[j80] = isign(j+1,i);
    // fprintf(stderr,"1: Cat[%ld]=%ld\n",j80,cat[j80]);
    for (long unsigned int i70=1;i70<=k70;i70++) {
      j = j + 1;
      from[j] = isign(j80,i);
      // fprintf(stderr,"from[%ld]=%ld\n",j,from[j]);
      cost[j] = 0;
      capac[j] = -big;
      floor[j] = 0;
    }
  }
  marc = j+1;
  assert(marc <= larc);
  from[marc] = isign(mnodp1,-i);
  // fprintf(stderr,"from[%ld]=%ld\n",marc,from[marc]);
  // fprintf(stderr,"Checkpoint %s:%d\n",__FILE__,__LINE__);
};

void netflo::arc_input_begin() {
  // fprintf(stderr,"Checkpoint %s:%d\n",__FILE__,__LINE__);
  kost0 = 0;
  // fprintf(stderr,"Checkpoint %s:%d\n",__FILE__,__LINE__);
}

long signed int netflo::infty() {
  return big-1;
}

void netflo::arc(long unsigned int i, long unsigned int j, long signed int k) {
  arc(i,j,k,infty(),0);
}

void netflo::arc(long unsigned int i, long unsigned int j, long signed int k,
		 long signed int l) {
  arc(i,j,k,l,0);
}

void netflo::arc(long unsigned int i, long unsigned int j, long signed int k, 
		 long signed int l, long signed int m) {
  // fprintf(stderr,"Checkpoint %s:%d\n",__FILE__,__LINE__);
  // fprintf(stderr,"%ld,%ld,%ld,%ld,%ld,%ld\n",i,j,k,l,m,big);
  assert(i>=1 && i<=mnode);
  assert(j>=1 && j<=mnode);
  assert(k <  big);
  assert(k > -big);
  assert(l <  big);
  assert(l > -big);
  assert(m <  big);
  assert(m > -big);
  // fprintf(stderr,"%ld\t%ld\n",m,l);
  assert(m <= l);
  long signed int ii = cat[j];
  long signed int jj = iabs(ii);
  /* if (outformat == 1) {
    ofs_.width(10);
    ofs_ << jj;
    ofs_.width(10);
    ofs_ << i;
    ofs_.width(10);
    ofs_ << j;
    ofs_.width(10);
    ofs_ << k;
    ofs_.width(10);
    ofs_ << l;
    ofs_.width(10);
    ofs_ << m << endl;
  } else if (outformat == 2) {
    ofs_ << "a " 
	 << i << " " 
	 << j << " " 
	 << m << " " 
	 << l << " " 
	 << k << endl;
	 }*/
  long signed int kk = isign(lnodp1,ii);
  // fprintf(stderr,"ii=%ld\n",ii);
  // fprintf(stderr,"jj=%ld\n",jj);
  // fprintf(stderr,"kk=%ld\n",kk);
  // fprintf(stderr,"from[jj]=%ld\n",from[jj]);
  // fprintf(stderr,"ieor(kk,from[jj])=%ld\n",ieor(kk,from[jj]));
  if (ieor(kk,from[jj]) <= 0) {
    // fprintf(stderr,"Checkpoint %s:%d\n",__FILE__,__LINE__);
    assert(marc != larc);
    marc++;
    long signed int k120 = marc-jj;
    long signed int m120 = marc;
    for (long unsigned int j120=1;j120<=k120;j120++) {
      long signed int l120 = m120-1;
      from[m120] = from[l120];
      // fprintf(stderr,"from[%ld]=%ld\n",m120,from[m120]);
      cost[m120] = cost[l120];
      capac[m120] = capac[l120];
      floor[m120] = floor[l120];
      m120 = l120;
    }
    // fprintf(stderr,"j=%ld\n",j);
    // fprintf(stderr,"mnode=%ld\n",mnode);
    for (long unsigned int j130=j;j130<=mnode;j130++) {
      cat[j130] = cat[j130]+isign(1,cat[j130]);
      // fprintf(stderr,"2: Cat[%ld]=%ld\n",j130,cat[j130]);
    }
  }
  // // fprintf(stderr,"Checkpoint %s:%d\n",__FILE__,__LINE__);
  from[jj] = isign(i,ii);
  // fprintf(stderr,"from[%ld]=%ld\n",jj,from[jj]);
  cost[jj] = k;
  kost0 = kost0+k*m;
  capac[jj] = l-m;
  floor[jj] = m;
  flow[i] = flow[i]-m;
  flow[j] = flow[j]+m;
  cat[j] = isign(jj+1,ii);
  // fprintf(stderr,"Edge has tag %ld\n",jj);
  arcid[i] = -1;
  // fprintf(stderr,"Checkpoint %s:%d\n",__FILE__,__LINE__);
}
void netflo::arc_input_end() {
  // fprintf(stderr,"Checkpoint %s:%d\n",__FILE__,__LINE__);
  /* if (outformat == 1) {
    ofs_.width(10);
    ofs_ << 0;
    ofs_.width(10);
    ofs_ << 0;
    ofs_.width(10);
    ofs_ << 0;
    ofs_.width(10);
    ofs_ << 0;
    ofs_.width(10);
    ofs_ << 0;
    ofs_.width(10);
    ofs_ << 0 << endl;
    }*/
  long signed int i = lnodp1;
  long signed int k = arty;
  long signed int l = 0;
  marc = marc-1;
  for (long unsigned int j190=artyp1;j190<=marc;j190++) {
    long signed int j = from[j190];
    if (ieor(i,j) <= 0) {
      i = -i;
      l++;
      cat[l] = k+1;
      // fprintf(stderr,"4: Cat[%ld]=%ld\n",l,cat[l]);
    } else {
      if (iabs(j) == l) continue;
    }
    k++;
    if (k == j190) continue;
    from[k] = from[j190];
    // fprintf(stderr,"from[%ld]=%ld\n",k,from[k]);
    cost[k] = cost[j190];
    capac[k] = capac[j190];
    floor[k] = floor[j190];
  }
  marc = k;
  mreg = k;
  assert(marc+max0(1,msorc)+1 <= larc);
  i = -from[marc];
  thd = next[mnodp1];
  next[mnodp1] = mnodp1;
  if (thd == mnodp1) {
    marc++;
    from[marc] = isign(mnodp1,i);
    // fprintf(stderr,"from[%ld]=%ld\n",marc,from[marc]);
    cost[marc] = 0;
    capac[marc] = -big;
    floor[marc] = 0;
  } else {
    do {
      marc++;
      from[marc] = isign(thd,i);
      // fprintf(stderr,"from[%ld]=%ld\n",marc,from[marc]);
      cost[marc] = 0;
      capac[marc] = level[thd];
      level[thd] = 0;
      floor[marc] = 0;
      long signed int nxt = next[thd];
      next[thd] = 0;
      thd = nxt;
    } while (thd != mnodp1);
  }
  mslk = marc;
  marc++;
  from[marc] = isign(mnodp2,-i);
  // fprintf(stderr,"from[%ld]=%ld\n",marc,from[marc]);
  cost[marc] = big;
  capac[marc] = 0;
  floor[marc] = 0;
  net = 0;
  mtree = 0;
  thd = mnodp1;
  for (long unsigned int i200 = 1;i200<=mnode;i200++) {
    long signed int j = flow[i200];
    net = net + j;
    if (j < 0) {
      flow[i200] = -j;
      dwn = mnodp1;
      long signed int nxt;
      while (1) {
	nxt = down[dwn];
	if (flow[nxt]+j <= 0) break;
	dwn = nxt;
      }
      down[dwn] = i200;
      down[i200] = nxt;
      level[i200] = -1;
    } else if (j > 0) {
      mtree = mtree+1;
      arcid[i200] = -marc;
      flow[i200] = j;
      next[thd] = i200;
      down[i200] = mnodp1;
      next[i200] = mnodp1;
      level[i200] = 1;
      dual[i200] = big;
      thd = i200;
    }
  }
  assert(net >= 0);
  // fprintf(stderr,"Checkpoint %s:%d\n",__FILE__,__LINE__);
}

netflo::solution_code 
netflo::solve() {
  // fprintf(stderr,"Checkpoint %s:%d\n",__FILE__,__LINE__);
  long signed int i,j,k,l,m;
  i=j=k=l=m=0;
 L1010:
      // fprintf(stderr,"L1010:\n");

  to = down[mnodp1];
  if (to == mnodp1) {
    goto L210;
  }
  
  // fprintf(stderr,"to=%ld\n",to);

 L1020:
  // fprintf(stderr,"L1020:\n");

  newarc = arty;
  newpr = big;
  if (flow[to] == 0) {
    goto L1110;
  }
  try_ = cat[to];
  // fprintf(stderr,"try_=%ld\n",try_);
  frm = from[try_];
  lst = isign(lnodp1,frm);
 L1030:
      // fprintf(stderr,"L1030:\n");
      
  // fprintf(stderr,"try_=%ld\n",try_);
  // fprintf(stderr,"capac[try_]=%ld\n",capac[try_]);
  // fprintf(stderr,"frm=%ld\n",frm);


  if (capac[try_] <= 0) {
    goto L1050;
  }
  fm = iabs(frm);
  // fprintf(stderr,"fm=%ld\n",fm);
  // fprintf(stderr,"level[fm]=%ld\n",level[fm]);
  if (level[fm] != 1) {
    goto L1050;
  }
  // fprintf(stderr,"arcid[fm]=%ld\n",arcid[fm]);
  // fprintf(stderr,"arty=%ld\n",arty);
  if (arcid[fm] == arty) {
    goto L1050;
  }
  price = cost[try_];
  // fprintf(stderr,"cost[try_]=%ld\n",cost[try_]);
  // fprintf(stderr,"price=%ld\n",price);
  if (price >= newpr) {
    goto L1050;
  }
  // fprintf(stderr,"flow[to]=%ld\n",flow[to]);
  // fprintf(stderr,"to=%ld\n",to);  
  if (capac[try_] > flow[to]) {
    goto L1040;
  }
  // fprintf(stderr,"flow[fm]=%ld\n",flow[fm]);
  if (flow[fm] < capac[try_]) {
    goto L1050;
  }
  newarc = -try_;
  newpr = price;
  if (newpr == 0) {
    goto L1055;
  }
  goto L1050;
 L1040:
      // fprintf(stderr,"L1040:\n");

  if (flow[fm] < flow[to]) {
    goto L1050;
  }
  newarc = try_;
  newpr = price;
  if (newpr == 0) {
    goto L1055;
  }
 L1050:
      // fprintf(stderr,"L1050:\n");

  try_++;
  frm = from[try_];
  // fprintf(stderr,"frm=%ld\n",frm);
  // fprintf(stderr,"lst=%ld\n",lst);
  // fprintf(stderr,"ieor(frm,lst)=%ld\n",ieor(frm,lst));
  if (ieor(frm,lst)>0) {
    goto L1030;
  }
  if (newarc == arty) {
    goto L1070;
  }
 L1055:
      // fprintf(stderr,"L1055:\n");

  if (newarc > 0) {
    goto L1060;
  }
  newarc = -newarc;
  fm = iabs(from[newarc]);
  flw = capac[newarc];
  capac[newarc] = -flw;
  flow[fm] = flow[fm]-flw;
  flow[to] = flow[to]-flw;
  goto L1020;
 L1060:
      // fprintf(stderr,"L1060:\n");

  capac[newarc] = -capac[newarc];
  fm = iabs(from[newarc]);
  flow[fm] = flow[fm]-flow[to];
  k = big;
  goto L1115;
 L1070:
      // fprintf(stderr,"L1070:\n");

  try_ = cat[to];
  frm = from[try_];
 L1080:
      // fprintf(stderr,"L1080:\n");

  if (capac[try_] <= 0) {
    goto L1090;
  }
  fm = iabs(frm);
  if (level[fm] != 0) {
    goto L1090;
  }
  price = cost[try_];
  if (price >= newpr) {
    goto L1090;
  }
  newarc = try_;
  newpr = price;
  if (newpr == 0) {
    goto L1095;
  }
 L1090:
      // fprintf(stderr,"L1090:\n");

  try_++;
  frm = from[try_];
  if (ieor(frm,lst)>0) {
    goto L1080;
  }
  if (newarc == arty) {
    goto L1110;
  }
 L1095:
      // fprintf(stderr,"L1095:\n");

  fm = iabs(from[newarc]);
  if (capac[newarc] > flow[to]) {
    goto L1100;
  }
  flw = capac[newarc];
  capac[newarc] = -flw;
  flow[fm] = flw;
  flow[to] = flow[to]-flw;
  down[fm]= to;
  down[mnodp1] = fm;
  level[fm] = -1;
  goto L1010;
 L1100:
      // fprintf(stderr,"L1100:\n");

  capac[newarc] = -capac[newarc];
  flow[fm] = flow[to];
  down[fm] = down[to];
  down[to] = fm;
  down[mnodp1] = fm;
  next[fm] = to;
  arcid[to] = newarc;
  level[fm] = level[to]-1;
  dual[to] = newpr;
  mtree++;
  goto L1010;
 L1110:
      // fprintf(stderr,"L1110:\n");

  k = 0;
 L1115:
      // fprintf(stderr,"L1115:\n");

  down[mnodp1] = down[to];
  fm = iabs(from[newarc]);
  arcid[to] = newarc;
  dual[to] = newpr;
  down[to] = fm;
  i = next[fm];
  next[fm] = to;
  j = level[fm]-level[to]+1;
  thd = fm;
 L1120:
      // fprintf(stderr,"L1120:\n");

  thd = next[thd];
  l = level[thd];
  level[thd] = l+j;
  k = k-dual[thd];
  dual[thd] = k;
  if (l != -1) {
    goto L1120;
  }
  next[thd] = i;
  mtree++;
  goto L1010;
 L210:
      // fprintf(stderr,"L210:\n");

  to = 1;
  try_ = artyp1;
  frm = from[try_];
 L220:
      // fprintf(stderr,"L220:\n");

  if (mtree == mnode) {
    goto L285;
  }
  too = to;
  newpr = big;
 L230:
      // fprintf(stderr,"L230:\n");

  lvj = level[to];
  lst = isign(lnodp1,frm);
 L235:
      // fprintf(stderr,"L235:\n");

  if (capac[try_] <= 0) {
    goto L260;
  }
  m = cost[try_];
  if (newpr < m) {
    goto L260;
  }
  fm = iabs(frm);
  if (level[fm] == 0) {
    goto L240;
  }
  if (lvj != 0) {
    goto L260;
  }
  i = fm;
  j = to;
  k = m;
  l = try_;
  goto L250;
 L240:
      // fprintf(stderr,"L240:\n");

  if (lvj == 0) {
    goto L260;
  }
  i = to;
  j = fm;
  k = -m;
  l = -try_;
 L250:
      // fprintf(stderr,"L250:\n");

  newpr = m;
 L260:
      // fprintf(stderr,"L260:\n");

  try_++;
  frm = from[try_];
  if (ieor(frm,lst)>0) {
    goto L235;
  }
  to++;
  if (to != mnodp1) {
    goto L270;
  }
  to = 1;
  try_ = artyp1;
  frm = from[try_];
 L270:
      // fprintf(stderr,"L270:\n");

  if (newpr != big) {
    goto L280;
  }
  if (to != too) {
    goto L230;
  }
  for (long unsigned int i275=1;i275<=mnode;i275++) {
    long signed int j275=0;
    if (level[i275] != 0) {
      continue;
    }
    if (arcid[i275] == -1) {
      goto L274;
    }
    j275 = cat[i275];
    if (iabs(from[j275]) != i275) {
      goto L274;
    }
    fprintf(stderr,"Node %ld is isolated.\n",i275);
  L274:
      // fprintf(stderr,"L274:\n");

    arcid[i275] = arty;
    flow[i275] = 0;
    next[i275] = next[mnodp1];
    next[mnodp1] = i275;
    down[i275] = mnodp1;
    level[i275] = 1;
    dual[i275] = -big;
  }
  goto L285;
 L280:
      // fprintf(stderr,"L280:\n");

  arcid[j] = l;
  down[j] = i;
  next[j] = next[i];
  next[i] = j;
  level[j] = level[i]+1;
  dual[j] = dual[i]-k;
  newarc = iabs(l);
  capac[newarc] = -capac[newarc];
  mtree++;
  goto L220;
 L285:
      // fprintf(stderr,"L285:\n");

  for (long unsigned int i290=1;i290<=mnode;i290++) {
    long signed int j290= iabs(arcid[i290]);
    capac[j290] = -capac[j290];
  }
      // fprintf(stderr,"L290:\n");
  for (long unsigned int i295=1;i295<=marc;i295++) {
    if (capac[i295]+big == 0) {
      capac[i295] = 0;
    }
  }
      // fprintf(stderr,"L295:\n");
  capac[arty] = big;
  capac[marc] = big;
  to = 1;
  try_ = artyp1;
  frm = from[try_];
  iter = 0;
  optim = false;
  dmp = true;
 L300:
      // fprintf(stderr,"L300:\n");

  // fprintf(stderr,"Checkpoint %s:%d\n",__FILE__,__LINE__);
  iter++;
 L305:
      // fprintf(stderr,"L305:\n");

  // fprintf(stderr,"Checkpoint %s:%d\n",__FILE__,__LINE__);
  too = to;
  newpr = 0;
 L310:
      // fprintf(stderr,"L310:\n");

  price0 = -dual[to];
  lst = isign(lnodp1,frm);
 L320:
      // fprintf(stderr,"L320:\n");

  fm = iabs(frm);
  price = dual[fm]+price0-cost[try_];
  if (capac[try_] < 0) goto L325;
  if (capac[try_] == 0) goto L330;
  if (capac[try_] > 0) goto L326;
 L325:
      // fprintf(stderr,"L325:\n");

  price = -price;
 L326:
      // fprintf(stderr,"L326:\n");

  if (price <= newpr) {
    goto L330;
  }
  newarc = try_;
  newpr = price;
  newto = to;
 L330:
      // fprintf(stderr,"L330:\n");

  try_++;
  frm = from[try_];
  if (ieor(frm,lst)>0) {
    goto L320;
  }
  to++;
  if (to != mnodp2) {
    goto L350;
  }
  to = 1;
  try_ = artyp1;
  frm = from[try_];
 L350:
      // fprintf(stderr,"L350:\n");

  if (newpr != 0) {
    goto L360;
  }
  if (to != too) {
    goto L310;
  }
  
  if (itdopt != 0) {
    dmp = true;
  }
  optim = true;
  goto L795;
 L360:
      // fprintf(stderr,"L360:\n");

  newfrm = iabs(from[newarc]);
  // fprintf(stderr,"newfrm=%ld\n",newfrm);
  // fprintf(stderr,"newto=%ld\n",newto);
  // fprintf(stderr,"newarc=%ld\n",newarc);  
  // fprintf(stderr,"from[newarc]=%ld\n",from[newarc]);
  
  theta = iabs(capac[newarc]);
  jtheta = 0;
  ch[2] = isign(larcp1,capac[newarc]);
  ch[1] = -ch[2];
  dw[1] = newfrm;
  dw[2] = newto;
  ldiff = level[newfrm]-level[newto];
  ktheta = 1;
  // fprintf(stderr,"ldiff=%ld\n",ldiff);  
  // fprintf(stderr,"level[newfrm]=%ld\n",level[newfrm]);  
  // fprintf(stderr,"level[newto]=%ld\n",level[newto]);  
  if (ldiff < 0) goto L380;
  if (ldiff == 0) goto L450;
  if (ldiff > 0) goto L390;
 L380:
      // fprintf(stderr,"L380:\n");

  ktheta = 2;
 L390:
      // fprintf(stderr,"L390:\n");

  dwn = dw[ktheta];
  chg = ch[ktheta];
  k440 = iabs(ldiff);
  for (long unsigned int i440=1;i440<=k440;i440++) {
    if (ieor(chg,arcid[dwn])>0) goto L410;
    i = iabs(arcid[dwn]);
    poss = capac[i]-flow[dwn];
    jposs = -dwn;
    goto L420;
  L410:
      // fprintf(stderr,"L410:\n");

    poss = flow[dwn];
    jposs = dwn;
  L420:
      // fprintf(stderr,"L420:\n");

    if (theta <= poss) goto L430;
    theta = poss;
    jtheta = jposs;
    if (theta == 0) goto L530;
  L430:
      // fprintf(stderr,"L430:\n");

    dwn = down[dwn];
  }
  // fprintf(stderr,"theta=%ld\n",theta);
  // fprintf(stderr,"L440:\n");
  dw[ktheta] = dwn;
 L450:
      // fprintf(stderr,"L450:\n");

 L460:
      // fprintf(stderr,"L460:\n");
      
  // fprintf(stderr,"dw[1]=%ld\n",dw[1]);
  // fprintf(stderr,"dw[2]=%ld\n",dw[2]);
  if (dw[1] == dw[2]) goto L520;
  for (long unsigned int l510=1;l510<=2;l510++) {
    dwn = dw[l510];
    // fprintf(stderr,"l510=%ld\n",l510);
    // fprintf(stderr,"dwn=%ld\n",dwn);
    // fprintf(stderr,"ch[l510]=%ld\n",ch[l510]);
    // fprintf(stderr,"arcid[dwn]=%ld\n",arcid[dwn]);
    // fprintf(stderr,"ieor(ch[l510],arcid[dwn])=%ld\n",ieor(ch[l510],arcid[dwn]));
    if (ieor(ch[l510],arcid[dwn])>0) goto L480;
    i = iabs(arcid[dwn]);
    // fprintf(stderr,"i=%ld\n",i);    
    // fprintf(stderr,"dwn=%ld\n",dwn);    
    // fprintf(stderr,"arcid[dwn]=%ld\n",arcid[dwn]);    
    // fprintf(stderr,"capac[i]=%ld\n",capac[i]);    
    // fprintf(stderr,"flow[dwn]=%ld\n",flow[dwn]);    
    poss = capac[i]-flow[dwn];
    jposs = -dwn;
    goto L490;
  L480:
      // fprintf(stderr,"L480:\n");

    poss = flow[dwn];
    jposs = dwn;
  L490:
      // fprintf(stderr,"L490:\n");
      // fprintf(stderr,"theta=%ld\n",theta);
      // fprintf(stderr,"poss=%ld\n",poss);
      
    if (theta <= poss) goto L500;
    theta = poss;
    jtheta = jposs;
    ktheta = l510;
    if (theta == 0) goto L530;
  L500:
    // fprintf(stderr,"L500:\n");
    
    dw[l510] = down[dwn];
    // fprintf(stderr,"L510:\n");
  }
  // fprintf(stderr,"L510:\n");

  goto L460;
 L520:
  // fprintf(stderr,"L520:\n");

  dwe = dw[1];
 L530:
      // fprintf(stderr,"L530:\n");
      // fprintf(stderr,"theta=%ld\n",theta);
  if (theta == 0) goto L560;
  dw[1] = newfrm;
  dw[2] = newto;
  // fprintf(stderr,"dwe=%ld\n",dwe);
  // fprintf(stderr,"dw[1]=%ld\n",dw[1]);
  // fprintf(stderr,"dw[2]=%ld\n",dw[2]);
  // fprintf(stderr,"jtheta=%ld\n",jtheta);
  // fprintf(stderr,"ktheta=%ld\n",ktheta);
  if (jtheta != 0) {
    dw[ktheta] = iabs(jtheta);
  }
  for (long unsigned int i550=1;i550<=2;i550++) {
    dwn = dw[i550];
    chg = isign(theta,ch[i550]);
  L540:
      // fprintf(stderr,"L540:\n");

    if (dwn == dwe) continue;
    flow[dwn] = flow[dwn]-chg*isign(1,arcid[dwn]);
    dwn = down[dwn];
    goto L540;
  }
      // fprintf(stderr,"L550:\n");
 L560:
      // fprintf(stderr,"L560:\n");

  // fprintf(stderr,"jtheta=%ld\n",jtheta);
  if (jtheta != 0) goto L570;
  capac[newarc] = -capac[newarc];
  goto L300;
 L570:
      // fprintf(stderr,"L570:\n");

  itheta = iabs(jtheta);
  if (jtheta>0) goto L590;
  j = iabs(arcid[itheta]);
  capac[j] = -capac[j];
 L590:
  // fprintf(stderr,"L590:\n");
  
  // fprintf(stderr,"theta=%ld\n",theta);
  // fprintf(stderr,"itheta=%ld\n",itheta);
  // fprintf(stderr,"jtheta=%ld\n",jtheta);
  // fprintf(stderr,"ktheta=%ld\n",ktheta);
  // fprintf(stderr,"newarc=%ld\n",newarc);
  
  flw = theta;
  if (capac[newarc]>0) goto L600;
  capac[newarc] = -capac[newarc];
  flw = capac[newarc]-flw;
  newpr = -newpr;
 L600:
      // fprintf(stderr,"L600:\n");

  if (ktheta == 2) goto L610;
  q1 = newfrm;
  q2 = newto;
  aid = -newarc;
  newpr = -newpr;
  goto L620;
 L610:
      // fprintf(stderr,"L610:\n");

  q1 = newto;
  q2 = newfrm;
  aid = newarc;
 L620:
      // fprintf(stderr,"L620:\n");

      // fprintf(stderr,"q1=%ld\n",q1);
      // fprintf(stderr,"q2=%ld\n",q2);
  i = q1;
  j = down[i];
  lstar = level[q2]+1;
  if (theta == 0) goto L730;
  chg = isign(theta,ch[ktheta]);
 L680:
      // fprintf(stderr,"L680:\n");

  dual[i] = dual[i]+newpr;
  n = flow[i];
  flow[i] = flw;
  dir = isign(1,arcid[i]);
  ref = iabs(arcid[i]);
  arcid[i] = aid;
  lsave = level[i];
  ldiff = lstar-lsave;
  level[i] = lstar;
  thd = i;
 L690:
      // fprintf(stderr,"L690:\n");

  // fprintf(stderr,"thd=%ld\n",thd);
  nxt = next[thd];
  // fprintf(stderr,"nxt=%ld\n",nxt);
  // fprintf(stderr,"leval[nxt]=%ld\n",level[nxt]);
  // fprintf(stderr,"lsave=%ld\n",lsave);
  if (level[nxt] <= lsave) goto L700;
  level[nxt] = level[nxt]+ldiff;
  dual[nxt] = dual[nxt]+newpr;
  thd = nxt;
  goto L690;
 L700:
      // fprintf(stderr,"L700:\n");

  k = j;
 L710:
      // fprintf(stderr,"L710:\n");

  l = next[k];
  if (l == i) goto L720;
  k = l;
  goto L710;
 L720:
      // fprintf(stderr,"L720:\n");

  if (i == itheta) goto L790;
  flw = n-chg*dir;
  aid = -isign(ref,dir);
  next[k] = nxt;
  next[thd] = j;
  k = i;
  i = j;
  j = down[j];
  down[i] = k;
  lstar++;
  goto L680;
 L730:
      // fprintf(stderr,"L730:\n");

 L740:
      // fprintf(stderr,"L740:\n");

  dual[i] = dual[i]+newpr;
  n = flow[i];
  flow[i] = flw;
  dir = isign(1,arcid[i]);
  ref = iabs(arcid[i]);
  arcid[i] = aid;
  lsave = level[i];
  ldiff = lstar-lsave;
  level[i] = lstar;
  thd = i;
 L750:
      // fprintf(stderr,"L750:\n");

  nxt = next[thd];
  if (level[nxt] <= lsave) goto L760;
  level[nxt] = level[nxt]+ldiff;
  dual[nxt] = dual[nxt]+newpr;
  thd = nxt;
  goto L750;
 L760:
      // fprintf(stderr,"L760:\n");

  k = j;
 L770:
      // fprintf(stderr,"L770:\n");

  l = next[k];
  if (l == i) goto L780;
  k = l;
  goto L770;
 L780:
      // fprintf(stderr,"L780:\n");

  if (i == itheta) goto L790;
  flw = n;
  aid = -isign(ref,dir);
  next[k] = nxt;
  next[thd] = j;
  k = i;
  i = j;
  j = down[j];
  down[i] = k;
  lstar++;
  goto L740;
 L790:
      // fprintf(stderr,"L790:\n");

  next[k] = nxt;
  next[thd] = next[q2];
  next[q2] = q1;
  down[q1] = q2;
  goto L300;
 L795:
      // fprintf(stderr,"L795:\n");

  infeas = false;
  kost = kost0;
  for (long unsigned int i830=1;i830<=mnode;i830++) {
    i = iabs(arcid[i830]);
    if (flow[i830]!=0 && cost[i] == big) infeas = true;
    kost = kost+cost[i]*flow[i830];
  }
  for (long unsigned int i840=1;i840<=mslk;i840++) {
    if (capac[i840] >= 0) continue;
    long signed int j840 = -capac[i840];
    kost = kost+cost[i840]*j840;
  }
  if (optim) {
    // fprintf(stderr,"Iterations: %ld\n",iter);
  }
  if (infeas) goto L850;
  if (optim) {
    // fprintf(stderr,"Cost: %ld\n",kost);
  }
  goto L860;
 L850:
      // fprintf(stderr,"L850:\n");

  if (optim) {
    // fprintf(stderr,"Infeasible.\n");
  }
 L860:
      // fprintf(stderr,"L860:\n");

  if (!optim) goto L305;
  // fprintf(stderr,"Iterations: %ld\n",iter);
  // fprintf(stderr,"Optimal indication.\n");
  
  pres = true;
  if (!pres) goto L885;
  for (long unsigned int i884=1;i884<=mnode;i884++) {
    long signed int j884 = iabs(arcid[i884]);
    capac[j884] = -flow[i884];
  }
 L885:
      // fprintf(stderr,"L885:\n");

  for (long unsigned int arcnam=1;arcnam<=marc;arcnam++) {
    arcflow_[arcnam] = 0;
    arcto_[arcnam] = 0;
  }

  // fprintf(stderr,"Arc Flows\n");
  // fprintf(stderr,"Arc\tFrom\tTo\tFlow\tCost\n");
  to = 1;
  try_ = artyp1;
  frm = from[try_];
 L8886:
      // fprintf(stderr,"L8886:\n");

  lst = isign(lnodp1,frm);
 L8888:
      // fprintf(stderr,"L8888:\n");

  flw = max0(0,-capac[try_])+floor[try_];
  arcto_[try_] = to;
  arcflow_[try_] = flw;
  if (flw == 0) goto L8889;
  fm = iabs(frm);
  // fprintf(stderr,"%ld\t%ld\t%ld\t%ld\t%ld\n",try_,fm,to,flw,cost[try_]);
 L8889:
      // fprintf(stderr,"L8889:\n");

  try_ = try_ + 1;
  frm = from[try_];
  if (ieor(frm,lst) > 0) goto L8888;
  to = to + 1;
  if (to != mnodp1) goto L8886;
 L892:
      // fprintf(stderr,"L892:\n");

  if (optim) {
    return optimal;
  } else if (infeas) {
    return infeasible;
  } else {
    return error;
  }
}

long signed int netflo::objective() {
  return kost;
}

long signed int netflo::iterations() {
  return iter;
}

long signed int netflo::netflow(long unsigned int nodeid) {
  assert(nodeid >= 1 && nodeid <= mnode);
  return flow[nodeid];
}

bool netflo::arcreal(long unsigned int arcid) {
  assert(arcid >= 1 && arcid<=marc);
  return (cost[arcid]==big);
}

long signed int netflo::arcflow(long unsigned int arcid) {
  assert(arcid >= 1 && arcid<=marc);
  return arcflow_[arcid];
}

long unsigned int netflo::arcfrom(long unsigned int arcid) {
  assert(arcid >= 1 && arcid<=marc);
  return iabs(from[arcid]);
}

long unsigned int netflo::arcto(long unsigned int arcid) {
  assert(arcid >= 1 && arcid<=marc);
  return arcto_[arcid];
}

long signed int netflo::arccost(long unsigned int arcid) {
  assert(arcid >= 1 && arcid<=marc);
  return cost[arcid];
}

long signed int netflo::arclb(long unsigned int arcid) {
  assert(arcid >= 1 && arcid<=marc);
  return floor[arcid];
}

long signed int netflo::arcub(long unsigned int arcid) {
  assert(arcid >= 1 && arcid<=marc);
  return floor[arcid]+capac[arcid];
}

#ifdef MAIN

#include <iostream> 
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

int
main(int argc,char *argv[]) {

  netflo nf;
  if (argc > 1) {
    ifstream ifs(argv[1]);
    assert(ifs);
    int nodes=0;
    char buffer[1024];
    while (ifs) {
      ifs.getline(buffer,1024);
      stringstream ss(buffer);
      ss >> nodes;
      break;
    }
    // cerr << "Nodes: " << nodes << endl;
    long signed int r,n;
    vector<long signed int> req(nodes,0);
    while (ifs) {
      ifs.getline(buffer,1024);
      stringstream ss(buffer);
      ss >> n >> r;
      if (n == 0) break;
      req[n-1] = r;
    }
    for (int i=0;i<nodes;i++) {
      // cerr << "Node " << i+1 << " requirement " << req[i] << endl;
    }
    vector<long signed int> arcsin(nodes,0);
    int i = 0;
    while (i < nodes && ifs) {
      ifs.getline(buffer,1024);
      stringstream ss(buffer);
      for (int j=0;j<8&&i<nodes;j++) {
	ss >> arcsin[i];
	i++;
      }
    }
    long unsigned int edges=0;
    for (int i=0;i<nodes;i++) {
      edges += arcsin[i];
    }
    // cerr << "Edges: " << edges << endl;
    for (int i=0;i<nodes;i++) {
      // cerr << "Node " << i+1 << " arcsin " << arcsin[i] << endl;
    }
    nf.nodes(nodes);
    nf.edges(edges);
    nf.netflow_input_begin();
    for (int i=0;i<nodes;i++) {
      if (req[i]!=0) {
	nf.netflow(i+1,req[i]);
      }
    }   
    nf.netflow_input_end();
    nf.arcsin_input_begin();
    for (int i=0;i<nodes;i++) {
      nf.arcsin(i+1,arcsin[i]);
    }   
    nf.arcsin_input_end();
     
    nf.arc_input_begin();
    long signed int nam,f,t,c,ub,lb;
    // char buffer[1024];
    while (ifs) {
      ifs.getline(buffer,1024);
      // cerr  << "B: " << buffer << endl;
      stringstream ss(buffer);
      nam = 0;
      t = 0;
      c = 0;
      ub = 0;
      lb = 0;
      ss >> nam >> f >> t >> c >> ub >> lb;
      // cerr << nam << endl;
      if (nam == 0) break;
//       cerr << "Arc: " 
// 	   << nam << " " 
// 	   << f << " " 
// 	   << t << " " 
// 	   << c << " " 
// 	   << ub << " " 
// 	   << lb << endl;
      if (lb != 0) {
	nf.arc(f,t,c,ub,lb);
      } else {
	nf.arc(f,t,c,ub);
      }
    }
    nf.arc_input_end();
  } else {
    nf.nodes(10);
    nf.edges(50);
    nf.netflow_input_begin();
    nf.netflow(8,98315);
    nf.netflow(9,-71754);
    nf.netflow_input_end();
    
    nf.arcsin_input_begin();
    nf.arcsin(1,3);
    nf.arcsin(2,8);
    nf.arcsin(3,4);
    nf.arcsin(4,3);
    nf.arcsin(5,6);
    nf.arcsin(6,6);
    nf.arcsin(7,6);
    nf.arcsin(8,6);
    nf.arcsin(9,5);
    nf.arcsin(10,3);
    nf.arcsin_input_end();
    
    nf.arc_input_begin();
    nf.arc(1,6,89643,45305);
    nf.arc(6,5,6564,8812);
    nf.arc(6,1,21387,73319);
    nf.arc(1,6,88766,49827);
    nf.arc(7,8,82939,23238);
    nf.arc(9,7,15298,14626);
    nf.arc(8,7,62907,5781);
    nf.arc(4,8,70292,56237);
    nf.arc(9,6,20322,48307);
    nf.arc(4,9,74205,83575);
    nf.arc(10,1,33008,63284);
    nf.arc(10,1,72746,1204);
    nf.arc(2,10,92988,46353);
    nf.arc(3,4,5730,67937);
    nf.arc(7,3,11796,5444);
    nf.arc(2,10,77810,61205);
    nf.arc(8,2,20910,75561);
    nf.arc(7,6,7321,36592);
    nf.arc(3,5,33186,17591);
    nf.arc(8,7,74729,64548);
    nf.arc(4,5,57803,6219);
    nf.arc(1,8,71630,31372);
    nf.arc(2,5,50033,71291);
    nf.arc(6,9,60744,66879);
    nf.arc(1,6,70905,13504);
    nf.arc(10,8,529,45529);
    nf.arc(4,3,18757,7315);
    nf.arc(3,2,93144,45602);
    nf.arc(6,5,60890,29893);
    nf.arc(2,7,45867,30878);
    nf.arc(7,10,12574,80777);
    nf.arc(8,4,60875,30490);
    nf.arc(6,2,31757,29600);
    nf.arc(9,2,75246,46560);
    nf.arc(9,2,95626,77814);
    nf.arc(7,3,88887,54290);
    nf.arc(9,7,11085,6138);
    nf.arc(1,7,79809,33792);
    nf.arc(2,5,1492,21364);
    nf.arc(3,9,67277,29512);
    nf.arc(9,2,40058,59982);
    nf.arc(6,9,3868,60435);
    nf.arc(6,2,5840,93849);
    nf.arc(3,9,1009,84544);
    nf.arc(3,8,20921,10224);
    nf.arc(5,2,54491,58076);
    nf.arc(5,3,56508,93927);
    nf.arc(4,8,28496,67105);
    nf.arc(5,4,23380,85587);
    nf.arc(5,6,95556,78970);
    nf.arc_input_end();
  }

  netflo::solution_code sc;
  sc = nf.solve();
 
  if (sc == netflo::optimal) {
    cout << "Nodes: " << nf.nodes() << endl;
    cout << "Edges: " << nf.edges() << endl;
    cout << "Arcs: " << nf.arcs() << endl;
    cout << "Objective: " << nf.objective() << endl;
    cout << "From\tTo\tFlow\tCost\n";
    for (int i=1;i<=nf.arcs();i++) {
      if (nf.arcflow(i)!=0) {
	cout << nf.arcfrom(i) << "\t"
	     << nf.arcto(i) << "\t"
	     << nf.arcflow(i) << "\t"
	     << nf.arccost(i) << "\t"
	     << endl;
      }
    }
  } else if (sc == netflo::infeasible) {
    cout << "Infeasible!" << endl;
  } else {
    cout << "Error!" << endl;
  }
}
#endif
