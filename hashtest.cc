#include <iostream>
#include <string>
#include "hash.h"
#include "char_io.t"

int
main(int argc, char **argv) {
  
#if defined(PROFILE)
  set_profile_signal_handler();
#endif

  string file = argv[1];
  Normalized<MapFileChars> cp(file+".seq",'$');
  string tmpl = argv[2];
  // hashset h(cp,4,tmpl.c_str());
  // taghashset h(cp,6,4,tmpl.c_str());
  hash *h = spacedselect(cp,4,tmpl.c_str());

  tic;
  while (h->next()) {
    hash_t v = h->value();
    // cerr << h.pos() << '\t' << h.str(v) << '\t' << v << '\t' << h.rcvalue() << '\t' << h.cvalue() << endl;
  }
  toc;

  exit(0);

  cp.reset();
  spaced h1(cp,4,tmpl.c_str());
  
  tic;
  while (h1.next()) {
    hash_t v = h1.value();
    // cerr << h1.pos() << '\t' << h1.str(v) << '\t' << v << '\t' << h1.rcvalue() << '\t' << h1.cvalue() << endl;
  }
  toc;

  exit(0);

//   contigshift h2(cp,4,h.weight());
//   tic;
//   while (h2.next()) {
//     hash_t v = h2.value();
//     // cerr << h.str(v) << endl;
//   }
//   timestamp("Contiguous hash using shift");
//   toc;

//   h.reset();
//   tic;
//   while (h.next()) {
//     hash_t v = h.value();
//     // cerr << h.str(v) << endl;
//   }
//   timestamp("Tagged Spaced seed set");
//   toc;

//   contigshift h1(cp,4,h.weight());
//   tic;
//   while (h1.next()) {
//     hash_t v = h1.value();
//     // cerr << h.str(v) << endl;
//   }
//   timestamp("Contiguous seed using shift");
//   toc;

}
