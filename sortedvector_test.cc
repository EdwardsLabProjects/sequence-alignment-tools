
#include <math.h>
#include <iostream>
#include "sortedvector.t"
#include "util.h"
#include "types.h"

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

int
main(int argc, char *argv[]) {

  sortedvector<FILE_POSITION_TYPE,FILE_POSITION_TYPE> v;

  sortedvector<FILE_POSITION_TYPE,FILE_POSITION_TYPE>::const_iterator it1,it2;

  srand48(time(NULL));

  for (int i=0;i<100;i++) {
    int count = 1+(int)floor(10*DRAND48);
    if (count == 1) {
      v.append(i*10,3);
    } else {
      v.append(i*10,1);
      for (int j=1;j<count-1;j++) {
	v.append(i*10,0);
      }
      v.append(i*10,2);
    }
  }

  FILE_POSITION_TYPE k;
  
  for (int i=0;i<10000;i++) {
    k = (int)floor(1000*DRAND48);
    try {
      it1 = v.locate_first_at_least(k);
      if (10*((k+9)/10) != it1->key() || it1->value() == 2 || it1->value() == 0) {
	cerr << "First at least " << k << " is (" << it1->key() << "," << it1->value() << ")" << endl;
      }
    } 
    catch (sortedvector<FILE_POSITION_TYPE,FILE_POSITION_TYPE>::KeyOutOfRange &) {
      if (k >= 0 && k <= 990) {
	cerr << "First at least " << k << " is out of range." << endl;      
      }
      it1 = v.end();
    }
    try {
      it2 = v.locate_last_at_most(k);
      if (10*(k/10) != it2->key() || it2->value() == 1 || it2->value()==0) {
	cerr << "Last at most " << k << " is (" << it2->key() << "," << it2->value() << ")" << endl;
      }
    } 
    catch (sortedvector<FILE_POSITION_TYPE,FILE_POSITION_TYPE>::KeyOutOfRange &) {
      cerr << "Last at most " << k << " is out of range." << endl;      
      it2 = v.end();
    }
    for (int j=0; j<1000; j++) {
      // int f = (int)floor(102*drand48()-1);
      int f = (int)floor(100*DRAND48);
      sortedvector<FILE_POSITION_TYPE,FILE_POSITION_TYPE>::const_iterator fit,it11,it21;
      fit = v.iter(f);
      try {
	if ((it11=v.finger_locate_first_at_least(fit,k)) != it1) {
	  checkpoint;
	  cerr << "Finger " << f << " first at least " << k 
	       << " == " << it11->key() << " != " << it1->key() << endl;
	}
      } 
      catch (sortedvector<FILE_POSITION_TYPE,FILE_POSITION_TYPE>::KeyOutOfRange &) {
	if (it1 != v.end()) {
	  checkpoint;
	  cerr << "Finger " << f << " first at least " << k << " is out of range, but non-finger is " << it1->key() << endl;      
	}
      }
      catch (sortedvector<FILE_POSITION_TYPE,FILE_POSITION_TYPE>::InvalidFinger &) {
	checkpoint;
	cerr << "Finger " << f << " is invalid" << endl;
      }
      fit = v.iter(f);
      try {
	if ((it21=v.finger_locate_last_at_most(fit,k)) != it2) {
	  checkpoint;
	  cerr << "Finger " << f << " last at most " << k << " == " << it21->key() << " != " << it2->key() << endl;
	}
      } 
      catch (sortedvector<FILE_POSITION_TYPE,FILE_POSITION_TYPE>::KeyOutOfRange &) {
	if (it2 != v.end()) {
	  checkpoint;
	  cerr << "Finger " << f << " last at most " << k << " is out of range, but non-finger is " << it2->key() << endl;      
	}
      }
      catch (sortedvector<FILE_POSITION_TYPE,FILE_POSITION_TYPE>::InvalidFinger &) {
	checkpoint;
	cerr << "Finger " << f << " is invalid" << endl;
      }
    }
  }

  return 0;
}
