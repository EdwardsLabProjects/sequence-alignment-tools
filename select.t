
#ifndef _SELECT_T
#define _SELECT_T

#include "util.h"
#include "fasta_io.h"
#include "fasta_io.t"

template<class H> 
void charmapping(FastaFile<H> *ff, int mapindex) {
  if (mapindex < 2) {
    return;
  }
  for (char f=0;f<127;f++) {
    char t;
    if ((t=charmap(mapindex,f)) != f) {
      ff->mapto(f,t);
    }
  }
}

template<class H>
FastaFile<H> * pick_fasta_file(std::string database,
			       int ffselect,
			       bool memmap,
			       bool alignments,
			       fasta_file_seq_params const & ffp,
			       bool verbose) {
  FastaFile<H> * ff=0;
  if ((ffselect == 0 && exist(database+".sqn")) || ffselect == 3) {
    if (!ffp.translate) {
      if (verbose) timestamp("Normalized sequence database...");
      if (memmap) {
	if (verbose) timestamp("Using mmap for sequence I/O...");
	if (!ffp.mapindex) {
	  ff = new IndexedFastaFile<Normalized<MapFileChars>,H>(database,alignments,ffp);
	} else {
	  ff = new IndexedFastaFile<Mapped<Normalized<MapFileChars> >,H>(database,alignments,ffp);
	  charmapping(ff,ffp.mapindex);
	  if (verbose) timestamp("Character mapping applied...");
        }
      } else {
	if (verbose) timestamp("Not using mmap for sequence I/O...");
	if (!ffp.mapindex) {
	  ff = new IndexedFastaFile<Normalized<BufferedFileChars>,H>(database,alignments,ffp);
	} else {
	  ff = new IndexedFastaFile<Mapped<Normalized<BufferedFileChars> >,H>(database,alignments,ffp);
	  charmapping(ff,ffp.mapindex);
	  if (verbose) timestamp("Character mapping applied...");
	}
      }
    } else {
      if (verbose) timestamp("Translated normalized sequence database...");
      if (memmap) {
	if (verbose) timestamp("Using mmap for sequence I/O...");
	if (!ffp.mapindex) {
	  ff = new IndexedFastaFile<Translated<Normalized<MapFileChars> >,H>(database,alignments,ffp);
	} else {
	  ff = new IndexedFastaFile<Mapped<Translated<Normalized<MapFileChars> > >,H>(database,alignments,ffp);
	  charmapping(ff,ffp.mapindex);
	  if (verbose) timestamp("Character mapping applied...");
	}
      } else {
	if (verbose) timestamp("Not using mmap for sequence I/O...");
	if (!ffp.mapindex) {
	  ff = new IndexedFastaFile<Translated<Normalized<BufferedFileChars> >,H>(database,alignments,ffp);
	} else {
	  ff = new IndexedFastaFile<Mapped<Translated<Normalized<BufferedFileChars> > >,H>(database,alignments,ffp);
	  charmapping(ff,ffp.mapindex);
	  if (verbose) timestamp("Character mapping applied...");
	}
      }
    } 
  } else if ((ffselect == 0 && exist(database+".sqz")) || ffselect == 4) {
    if (!ffp.translate) {
      if (verbose) timestamp("Compressed sequence database...");
      if (memmap) {
	if (verbose) timestamp("Using mmap for sequence I/O...");
	if (!ffp.mapindex) {
	  ff = new IndexedFastaFile<Compressed<MapFileChars>,H>(database,alignments,ffp);
	} else {
	  ff = new IndexedFastaFile<Mapped<Compressed<MapFileChars> >,H>(database,alignments,ffp);
	  charmapping(ff,ffp.mapindex);
	  if (verbose) timestamp("Character mapping applied...");
	}
      } else {
	if (verbose) timestamp("Not using mmap for sequence I/O...");
	if (!ffp.mapindex) {
	  ff = new IndexedFastaFile<Compressed<BufferedFileChars>,H>(database,alignments,ffp);
	} else {
	  ff = new IndexedFastaFile<Mapped<Compressed<BufferedFileChars> >,H>(database,alignments,ffp);
	  charmapping(ff,ffp.mapindex);
	  if (verbose) timestamp("Character mapping applied...");
	}
      }
    } else {
      if (verbose) timestamp("Translated compressed sequence database...");
      if (memmap) {
	if (verbose) timestamp("Using mmap for sequence I/O...");
	if (!ffp.mapindex) {
	  ff = new IndexedFastaFile<Translated<Compressed<MapFileChars> >,H>(database,alignments,ffp);
	} else {
	  ff = new IndexedFastaFile<Translated<Compressed<MapFileChars> >,H>(database,alignments,ffp);
	  charmapping(ff,ffp.mapindex);
	  if (verbose) timestamp("Character mapping applied...");	  
	}
      } else {
	if (verbose) timestamp("Not using mmap for sequence I/O...");
	if (!ffp.mapindex) {
	  ff = new IndexedFastaFile<Translated<Compressed<BufferedFileChars> >,H>(database,alignments,ffp);
	} else {
	  ff = new IndexedFastaFile<Mapped<Translated<Compressed<BufferedFileChars> > >,H>(database,alignments,ffp);
	  charmapping(ff,ffp.mapindex);
	  if (verbose) timestamp("Character mapping applied...");
	}
      }
    }
  } else if ((ffselect == 0 && exist(database+".seq")) || ffselect == 2) {
    if (!ffp.translate) {
      if (verbose) timestamp("Indexed sequence database...");
      if (memmap) {
	if (verbose) timestamp("Using mmap for sequence I/O...");
	if (!ffp.mapindex) {
	  ff = new IndexedFastaFile<MapFileChars,H>(database,alignments,ffp);
	} else {
	  ff = new IndexedFastaFile<Mapped<MapFileChars>,H>(database,alignments,ffp);
	  charmapping(ff,ffp.mapindex);
	  if (verbose) timestamp("Character mapping applied...");
	}
      } else {
	if (verbose) timestamp("Not using mmap for sequence I/O...");
	if (!ffp.mapindex) {
	  ff = new IndexedFastaFile<BufferedFileChars,H>(database,alignments,ffp);
	} else {
	  ff = new IndexedFastaFile<Mapped<BufferedFileChars>,H>(database,alignments,ffp);	  
	  charmapping(ff,ffp.mapindex);
	  if (verbose) timestamp("Character mapping applied...");
	}
      }
    } else {
      if (verbose) timestamp("Translated indexed sequence database...");
      if (memmap) {
	if (verbose) timestamp("Using mmap for sequence I/O...");
	if (!ffp.mapindex) {
	  ff = new IndexedFastaFile<Translated<MapFileChars>,H>(database,alignments,ffp);
	} else {
	  ff = new IndexedFastaFile<Mapped<Translated<MapFileChars> >,H>(database,alignments,ffp);
	  charmapping(ff,ffp.mapindex);
	  if (verbose) timestamp("Character mapping applied...");
	}
      } else {
	if (verbose) timestamp("Not using mmap for sequence I/O...");
	if (!ffp.mapindex) {
	  ff = new IndexedFastaFile<Translated<BufferedFileChars>,H>(database,alignments,ffp);
	} else {
	  ff = new IndexedFastaFile<Mapped<Translated<BufferedFileChars> >,H>(database,alignments,ffp);
	  charmapping(ff,ffp.mapindex);
	  if (verbose) timestamp("Character mapping applied...");
	}
      }      
    }
  } else { /* either all files didn't exist or ffselect == 1 */
    if (ffp.translate) {
      timestamp("Can\'t translate from raw sequence database");
    }
    if (verbose) timestamp("Raw sequence database...");
    if (memmap) {
      if (verbose) timestamp("Using mmap for sequence I/O...");
      if (!ffp.mapindex) {
	ff = new StreamedFastaFile<MapFileChars,H>(database,alignments,ffp);
      } else {
	ff = new StreamedFastaFile<Mapped<MapFileChars>,H>(database,alignments,ffp);
	charmapping(ff,ffp.mapindex);
	if (verbose) timestamp("Character mapping applied...");
      }
    } else {
      if (verbose) timestamp("Not using mmap for sequence I/O...");
      if (!ffp.mapindex) {
	ff = new StreamedFastaFile<BufferedFileChars,H>(database,alignments,ffp);
      } else {
	ff = new StreamedFastaFile<Mapped<BufferedFileChars>,H>(database,alignments,ffp);
	charmapping(ff,ffp.mapindex);
	if (verbose) timestamp("Character mapping applied...");
      }
    }
  }
  return ff;
}
  
#endif

