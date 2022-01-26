#include <ctype.h>
#include "xspacefsm.h"
#include <string.h>

#define WRAPZ(x) ((x>255)?255:x)

void XSpaceFSM::initialize() {
  Z = new unsigned char[Tree->num_suffix()];
  for(unsigned i=0; i < Tree->num_suffix(); ++i) Z[i] = 0;
  
  nZ = new unsigned char[Tree->num_nodes()];
  for(unsigned i=0; i < Tree->num_nodes(); ++i) nZ[i] = 0;
}

bool XSpaceFSM::interesting(st_index n,char cl,char cr) const {
  if (Tree->str(n)[-1]!=cl || Tree->str(n)[-1]==Tree->str()[0]
      ||
      Tree->str(n)[mersize]!=cr || Tree->str(n)[mersize]==Tree->str()[0])
    return true;
  
  if (!n.is_leaf()) {
    n = Tree->first_child(n);
    for(n.set_good(); n.is_good(); n = Tree->next_child(n)) {
      if (interesting(n,cl,cr)) return true;
    }
  }
  
  return false;
}

void XSpaceFSM::selfprocess(st_index n) const {
  if (!n.is_leaf()) {
    if (Tree->node(n).length() >= mersize) {
      if (interesting(n,Tree->str(n)[-1],Tree->str(n)[mersize])) {
	nZ[n.index()] = 1;
      }
    }
    else {
      n = Tree->first_child(n);
      for(n.set_good(); n.is_good(); n = Tree->next_child(n)) selfprocess(n);
    }
  }
}

void XSpaceFSM::selfstream() {
  selfprocess(Tree->root());
}

void XSpaceFSM::output(FILE* o, FILE_POSITION_TYPE offset) const {
  output_nodes(o,offset,Tree->root());
}

void XSpaceFSM::output_nodes(FILE* o, FILE_POSITION_TYPE offset, st_index n) const {
  if (!n.is_leaf()) {
    if (nZ[n.index()]) {
      print_nodes(o,offset,n);
      fputc('\n',o);
    }
    else {
      n = Tree->first_child(n);
      for(n.set_good(); n.is_good(); n = Tree->next_child(n)) 
	output_nodes(o,offset,n);
    }
  }
  else {
    if (Z[n.index()]) {
      print_nodes(o,offset,n);
      fputc('\n',o);
    }
  }
}

void XSpaceFSM::print_nodes(FILE *o, FILE_POSITION_TYPE offset, st_index n) const {
  if (!n.is_leaf()) {
    n = Tree->first_child(n);
    for(n.set_good(); n.is_good(); n = Tree->next_child(n)) print_nodes(o,offset,n);
  }
  else {
    offset += Tree->str(n)+mersize-Tree->str();
    fprintf(o," %llu.%c",
	    offset,
	    Tree->str(n)[mersize]);
  }
}

void XSpaceFSM::printnode(st_index n, bool parent=true, unsigned seen=0) const {
  if (parent) {
    cerr << "Current: ";
  } else {
    cerr << "  Child: ";
  }
  n.print(cerr);
  cerr << '\t';
  if (!parent) {
    for (int i=0;i<seen;i++) {
      cerr << " ";
    }
  }
  if (n.is_leaf()) {
    for (int i=0;i<20;i++) {
      if (!Tree->str(n)[seen+i]) break;
      cerr << Tree->str(n)[seen+i];
    }
    cerr << endl;
  } else { /* is node */
    for (int i=seen;i<Tree->node(n).length();i++) {
      cerr << Tree->str(n)[i];
    }
    cerr << endl;
    if (parent) {
      st_index ch = Tree->first_child(n);
      for(ch.set_good(); ch.is_good(); ch = Tree->next_child(ch)) {
	printnode(ch,false,Tree->node(n).length());
      }
    }
  }
}

void XSpaceFSM::stream(const char *map,FILE *f,long unsigned len) {
  // creating some variable shortcuts
  const char *head=Tree->str();
  const st_index Root = Tree->root();
  const SufTree &T = *Tree;

  char c; // holds the current character
  char c0; 
  st_index n = Root; // holds current node
  unsigned dep = 0; // current depth
  unsigned char curZ=0;

  unsigned flen = 0; // fast match length for internal quick traversal

  unsigned int bufpos = 0;
  unsigned int bufpos_ahead = 0;
  unsigned int bufsize = mersize+2;
  char *buffer = new char[bufsize];

  for (unsigned int i=0;i<bufsize;i++) {
    buffer[i] = Tree->str()[0];
  }
  
  //  char hist[20];
  //  for(int i=0; i < 20; ++i) hist[i] = 0;

  // Tree->dump();

  for(long unsigned pos = 0; pos < len; ++pos) {

//     cerr << "Valid state!" << endl;
//     cerr << "Current buffer: ";
//     for (int i=0;i<bufsize;i++) {
//       cerr << buffer[(bufpos_ahead+i)%bufsize];
//     }
//     cerr << endl;
//     printnode(n);
//     cerr << "dep: " << dep << " flen: " << flen << endl;

    c = map[(unsigned char) fgetc(f)];
//     cerr << "Read character: " << c << endl;
    bufpos = bufpos_ahead;
    ++bufpos_ahead;
    bufpos_ahead %= bufsize;
    buffer[bufpos] = c;
    c0 = buffer[bufpos_ahead];

    for(;;) { // update state

//       cerr << "Current buffer: ";
//       for (int i=0;i<bufsize;i++) {
//  	cerr << buffer[(bufpos_ahead+i)%bufsize];
//       }
//       cerr << endl;
//       printnode(n);
//       cerr << "dep: " << dep << " flen: " << flen << endl;

      if (T.is_suffix(n)) {
	if (dep >= flen) {
	  if (dep < mersize && T.str(n)[dep] == c) {
	    ++dep;
	    break; // we accept char so we leave the for()
	  }
	  if (dep >= mersize) {
	    // If we are at the correct depth, then we must decide
	    // whether or not this suffix is interesting...
	    if (!Z[n.index()] && interesting(n,c0,c)) {
	      Z[n.index()] = 1;
	      // printnode(n);
	      // cerr << "Interesting buffer: ";
	      // for (int i=1;i<bufsize-1;i++) {
	      // cerr << buffer[(bufpos_ahead+i)%bufsize];
	      // }
	      // cerr << endl;
	    }
	  }
	  // Take suffix, whether or not we match...
	  head = T.str(n)+1; 
	  // go to suffix link of parent
	  for(n.set_good();n.is_good();n=T.next_child(n)) ;
	  n.set_good();
	  flen = dep? dep-1 : 0;
	  dep = T.node(n).length();
#ifdef DEBUG2
	  fprintf(stdout,"suffixing to ");
	  print_node(n,dep);
	  fprintf(stdout," with fastlen %u\n",flen);
#endif
	  // if (flen == mersize) break;
	}
	else { // (dep < flen) // in leaf so dep can be moved to flen
	  dep = flen;
	}
      } 
      else { // non-leaf
	if (dep >= flen) { // no recovery needed
	  if (dep < T.node(n).length()) {
	    if (dep < mersize && T.str(n)[dep] == c) {
	      ++dep;
	      break; // we accept char to we leave the for()
	    }
	    if (dep >= mersize) {
	      // If we are at the correct depth, then we must decide
	      // whether or not this node is interesting...
	      if (!nZ[n.index()] && interesting(n,c0,c)) {
		nZ[n.index()] = 1;
		// printnode(n);
		// cerr << "Interesting buffer: ";
		// for (int i=1;i<bufsize-1;i++) {
		// cerr << buffer[(bufpos_ahead+i)%bufsize];
		// }
		// cerr << endl;
	      }
	    }
	    if (dep <= 1) { // fail to root
	      n = Root;
	      flen = 0;
	      dep = 0;
	    }
	    else { // fail to suffix
	      head = T.str(n)+1; // remember the current string
	      // go to suffix link of parent
	      for(n.set_good(); n.is_good(); n = T.next_child(n)) ;
	      n.set_good();
	      flen = dep? dep-1 : 0;
	      dep = T.node(n).length();
	    }
	  }
	  else { // dep == length(), must branch based on c
#ifdef DEBUG2
	    fprintf(stdout,"branching on child of ");
	    print_node(n,dep); fprintf(stdout,"\n");
#endif
	    if (dep < mersize+1) { // must branch based on c
	      
	      st_index ch = T.first_child(n);
	      if (ch.is_good()) {
		for(ch.set_good();
		    ch.is_good() && T.str(ch)[dep] != c; 
		    ch = T.next_child(ch)) ;
	      }
	      else {
		ch = n;
		for(ch.set_good();ch.is_good(); ch = T.next_child(ch)) ;
	      }
	      if (ch.is_nogood()) {
		head = T.str(n)+1; 
		flen = dep? dep-1 : 0;
		dep = T.node(ch).length();
		if (T.is_root(n)) { n = ch; break;}
	      }
	      n = ch;
	    } else {
	      // don't care what children are here, check if this node is interesting...
	      if (interesting(n,c0,c)) {
		nZ[n.index()] = 1;
	      }
	      head = T.str(n)+1; // remember the current string
	      // go to suffix link of parent
	      for(n.set_good(); n.is_good(); n = T.next_child(n)) ;
	      n.set_good();
	      flen = dep? dep-1 : 0;
	      dep = T.node(n).length();
	    }
	  }
	}
	else { // (dep < flen) { // some recovery needed
#ifdef DEBUG2
	  print_node(n,dep);
	  fprintf(stdout," needs recovery\n");
#endif
	  if (T.node(n).length() < flen) {
	    dep = T.node(n).length();
	    
	    st_index ch = T.first_child(n);
	    if (ch.is_good()) {
	      n = ch;
	      for(n.set_good();
		  n.is_good() && T.str(n)[dep] != head[dep]; 
		  n = T.next_child(n)) ;
	    }
	    else {
	      for(n.set_good();n.is_good(); n = T.next_child(n)) ;
	    }
	    if (n.is_nogood()) {
	      ++head;
	      for(n.set_good(); n.is_good(); n = T.next_child(n)) ;
	      n.set_good();
	      flen = dep? dep-1 : 0;
	      dep = T.node(n).length();
#ifdef DEBUG2
	      fprintf(stdout,"followed suffix to ");
	      print_node(n,dep);
	      fprintf(stdout," with fastlen %u\n",flen);
#endif
	    }
#ifdef DEBUG2
	    fprintf(stdout,"found child at ");
	    print_node(n,dep);
	    fprintf(stdout," with fastlen %u\n",flen);
#endif
	  }
	  else {
	    // loop around
	    dep = flen;
#ifdef DEBUG2
	    fprintf(stdout,"extend depth at ");
	    print_node(n,dep);
	    fprintf(stdout," with fastlen %u\n",flen);
#endif
	  }
	}
      }
    } // for(;;)

    // exactly 1 new char was processed by the loop above
    // and dep and n are now up to date

  } // for (pos < len)

  delete [] buffer;
}
