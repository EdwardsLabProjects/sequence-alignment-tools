
#include <assert.h>
#include <string.h>
#include <strings.h>
#include "rlst.h"
#include "util.h"

suftree::suftree(unsigned char eos) {
  eos_ = eos;
  poslen = 1000;
  pos = new unsigned[poslen];
}

suftree::~suftree() {
  delete [] pos;
  delete [] stree;
}

unsigned long suftree::add_pattern(std::string const & pat, long unsigned id,
				   int esb, int eeb) {
  add_pattern_(pat,id,esb,eeb);
  return id;
}

void suftree::init(CharacterProducer & cp) {
  stree = new RLSufTree(cp.c_str(),cp.length(),eos_);
  stree->read(suftree::filename(cp.filename()).c_str());
  reset();
}

void suftree::reset() {
  curit = patterns().begin();
}

bool suftree::find_patterns(CharacterProducer & cp,
			    pattern_hit_vector & pas,
			    long unsigned minpa) {
  unsigned pacount = 0;
  while (curit != patterns().end()) {
    char *buffer = new char[curit->pattern().length()+1];
    strcpy(buffer,curit->pattern().c_str());
    for (int i=0;i<curit->pattern().length();i++) {
      buffer[i] = cp.nch(buffer[i]);
    }
    unsigned npos = stree->find(buffer,pos,poslen);
    if (npos > poslen) {
      delete [] pos;
      poslen = npos;
      pos = new unsigned[poslen];
      npos = stree->find(buffer,pos,poslen);
      assert(npos <= poslen);
    }
    unsigned patlen = curit->pattern().length();
    for (int i=0;i<npos;i++) {
      pas.push_back(pos[i]+patlen,make_pair(curit,0));
      pacount++;
    }
    ++curit;
    if (pacount > minpa) {
      break;
    }
  }
  pas.normalize();
  if (pacount > 0) return true;
  return false;
}

