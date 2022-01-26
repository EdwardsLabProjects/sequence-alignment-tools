
# External depedencies

ATACINC = .../include
ATACLIB = .../lib

CS2 = ../CS2
PGEN = ../primegen

#
# Automaticly chose the platform based on the value of uname...
#

OS = $(shell uname)

ifeq ($(strip $(OS)),OSF1) 
  PLATFORM=alpha-cxx
endif

ifeq ($(strip $(OS)),AIX) 
  PLATFORM=aix-xlc
endif

ifeq ($(strip $(OS)),Linux) 
  PLATFORM=linux-gpp
endif

ifeq ($(strip $(OS)),CYGWIN_NT-5.1) 
  PLATFORM=cygwin-gpp
endif

ifeq ($(strip $(OS)),SunOS) 
  PLATFORM=sunos-CC
endif

ifeq ($(strip $(OS)),Darwin) 
  PLATFORM=darwin-gpp
endif

#
# Set the platform/machine that this'll be compiled on
#
# PLATFORM={alpha-cxx,alpha-gpp,linux-gpp,aix-xlc,aix-gpp,cygwin-gpp}
#

# PLATFORM=aix-gpp
# PLATFORM=aix-xlc
PLATFORM=linux-gpp
# PLATFORM=linux-icc
# PLATFORM=alpha-cxx
# PLATFORM=alpha-gpp
# PLATFORM=cygwin-gpp
# PLATFORM=sunos-CC
# PLATFORM=darwin-gpp

ifeq ($(strip $(PLATORM)),)
  CXX=cannot determine platform
endif

ifeq ($(strip $(PLATFORM)),aix-xlc) 
  CXX      = xlC -+ -q64 -D_LARGE_FILES -bhalt:8 -qnolm 
  CC       = xlc -+ -q64 -D_LARGE_FILES -bhalt:8 -qnolm 
  # PROFILE  = -pg -g -DPROFILE -O -qarch=auto -qtune=auto -qcache=auto -qmaxmem=-1
  OPTIMIZE = -O -qinline -qarch=auto -qtune=auto -qcache=auto -qmaxmem=-1 
  # DEBUG = -g # -DDEBUG_MEMORY
endif

ifeq ($(strip $(PLATFORM)),aix-gpp) 
  CXX         = /usr/bin/g++ -maix64 -D_LARGE_FILES
  CC          = /usr/bin/gcc -maix64 -D_LARGE_FILES
  OPTIMIZE    = -O
endif

ifeq ($(strip $(PLATFORM)),alpha-cxx) 
  CXX           = cxx # -ieee
  CC            = cc # -ieee
  # PROFILE       = -pg -g3 -DPROFILE
  # PROFILE       = -p -g3 -DPROFILE
  OPTIMIZE      = -O
  EXTRACXXFLAGS = -DNO_STD_NAMESPACE
  # DEBUG = -DDEBUG_MEMORY
  # DEBUG = -g
endif

ifeq ($(strip $(PLATFORM)),alpha-gpp) 
  CXX           = g++
  CC            = gcc
  OPTIMIZE      = -O
endif

ifeq ($(strip $(PLATFORM)),linux-gpp) 
  CXX          = g++ -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE
  CC           = gcc -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE
  OPTIMIZE     = -O3 
  DEBUG        = -g # -DDEBUG_MEMORY
  # PROFILE      = -g -pg -DPROFILE
endif

ifeq ($(strip $(PLATFORM)),linux-icc) 
  CXX          = icpc -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE
  CC           = icc -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE
  OPTIMIZE     = -fast -v
  DEBUG        = # -g # -DDEBUG_MEMORY
  # PROFILE = -g -pg -DPROFILE
endif

ifeq ($(strip $(PLATFORM)),cygwin-gpp) 
  CXX          = g++ -mno-cygwin -DNO_ZLIB
  CC           = gcc -mno-cygwin -DNO_ZLIB
  # CXX          = /cygdrive/c/MinGW/bin/g++.exe
  OPTIMIZE     = -O 
  # DEBUG        = -g 
  E            = .exe
endif

ifeq ($(strip $(PLATFORM)),sunos-CC) 
  CXX          = /opt/SUNWspro/bin/CC -DSIXTYFOURBITS -xarch=v9b -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
  CC           = /opt/SUNWspro/bin/cc -DSIXTYFOURBITS -xarch=v9b -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
  OPTIMIZE     = -O  
  # PROFILE      = -g -xpg -DPROFILE
endif

ifeq ($(strip $(PLATFORM)),darwin-gpp) 
  CXX          = g++
  OPTIMIZE     = -O 
  E            = 
endif                                                                                                                        

#
# The rest of these options should be left as is...
#

CXXFLAGS = $(EXTRACXXFLAGS) $(WARN) $(OPTIMIZE) $(DEBUG) $(PROFILE) $(INCLUDEPATH) $(PARAM)
LDFLAGS =  $(EXTRALDFLAGS) $(WARN) $(OPTIMIZE) $(DEBUG) $(PROFILE) $(LIBRARYPATH) $(PARAM)
LDLIBS = -lz -lm 
CFLAGS = $(CXXFLAGS)

SOURCES = keyword_tree.cc keyword_tree.h \
	  shift_and.cc shift_and.h \
	  shift_and_inexact.cc shift_and_inexact.h \
          mapFile.cc mapFile.h \
	  bufferedFile.cc bufferedFile.h \
	  fileStar.cc fileStar.h \
	  char_io.cc char_io.h \
	  fasta_io.cc fasta_io.h \
	  util.h util.cc \
	  primer_alignment.h primer_alignment.cc \
	  pattern_alignment.h pattern_alignment.cc \
	  pattern_match.h pattern_match.cc \
	  compress_seq.cc \
          exact_match.cc \
	  primer_match.cc \
	  inexact_match.cc \
	  types.cc \
	  sortedvector_test.cc \
	  atac_seq.cc extract_seq.cc select.cc \
	  pcr_match.cc sts_io.cc \
	  kmer_count.cc kmer_annotate.cc polyrun.cc \
	  exact_bases.cc exact_halves.cc filter_bitvec.cc \
	  aacomp.cc aacomplookup.cc hash_table.cc \
	  gs_hash_table.cc rand_hash_table.cc polyrun.cc \
	  walk_graph.cc netflo.cc cannon_csbh_graph.cc csbh_annotate.cc \
	  nrdb.cc peptide_scan.cc protein_mw.cc \
	  peptide_mult.cc \
	  Indexer.cc IndexerAA.cc WordGraph.cc Xspace.cc XspaceLo.cc \
	  xspacefsm.cc rl_index.cc word_graph.cc solid_simulation.cc \
	  solid_assembly.cc allvall.cc hash.cc perfposht.cc allvall_dump.cc \
	  allvall_tobm.cc genome_simulation.cc pairscan.cc suftree.cc rlst.cc hashtest.cc

PROGS = exact_match$(E) compress_seq$(E) inexact_match$(E) \
	 primer_match$(E) \
	 chario$(E) sortedvector_test$(E) kmer_count$(E) \
	 kmer_annotate$(E) \
	 polyrun$(E) pcr_match$(E) \
	 aacomp$(E) aacomplookup$(E) \
	 cannon_csbh_graph$(E) csbh_annotate$(E) nrdb$(E)\
	 peptide_scan$(E) protein_mw$(E) peptide_mult$(E)\
	 Indexer$(E) IndexerAA$(E) WordGraph$(E) Xspace$(E) XspaceLo$(E) \
	 xmers$(E) solid_simulation$(E) solid_assembly$(E) allvall$(E) pairscan$(E) allvall_dump$(E) \
	 allvall_tobm$(E) genome_simulation$(E) test$(E)
# 	 walk_graph$(E) \


PRIMER_MATCH_OBJ = primer_match.o pattern_alignment.o keyword_tree.o \
	 shift_and.o mapFile.o bufferedFile.o char_io.o fasta_io.o util.o \
	 pattern_match.o primer_alignment.o fileStar.o types.o select.o \
	 shift_and_inexact.o exact_bases.o exact_halves.o filter_bitvec.o \
	 sts_io.o hash_table.o gs_hash_table.o rand_hash_table.o rlst.o  oligotm.o

NRDB_OBJ = nrdb.o pattern_alignment.o keyword_tree.o \
	   shift_and.o mapFile.o bufferedFile.o char_io.o fasta_io.o util.o \
	   pattern_match.o primer_alignment.o fileStar.o types.o select.o \
	   shift_and_inexact.o exact_bases.o exact_halves.o filter_bitvec.o \
	   sts_io.o hash_table.o gs_hash_table.o rand_hash_table.o

EXACT_MATCH_OBJ = exact_match.o pattern_alignment.o keyword_tree.o \
	 shift_and.o mapFile.o bufferedFile.o char_io.o fasta_io.o util.o \
	 pattern_match.o fileStar.o types.o select.o \
	 shift_and_inexact.o exact_bases.o primer_alignment.o exact_halves.o \
	 filter_bitvec.o hash_table.o gs_hash_table.o rand_hash_table.o	rlst.o

PEPSCAN_OBJ = peptide_scan.o pattern_alignment.o keyword_tree.o \
	 shift_and.o mapFile.o bufferedFile.o char_io.o fasta_io.o util.o \
	 pattern_match.o fileStar.o types.o select.o \
	 shift_and_inexact.o exact_bases.o primer_alignment.o exact_halves.o \
	 filter_bitvec.o hash_table.o gs_hash_table.o rand_hash_table.o rlst.o

INEXACT_MATCH_OBJ = inexact_match.o keyword_tree.o shift_and.o \
		     mapFile.o bufferedFile.o char_io.o fasta_io.o util.o \
		     shift_and_inexact.o pattern_alignment.o pattern_match.o \
		     fileStar.o types.o select.o exact_bases.o \
		     primer_alignment.o exact_halves.o filter_bitvec.o \
		     hash_table.o gs_hash_table.o rand_hash_table.o rlst.o

COMPRESS_SEQ_OBJ = compress_seq.o fasta_io.o char_io.o util.o mapFile.o \
			 bufferedFile.o fileStar.o types.o

CHARIO_OBJ = chario.o fasta_io.o char_io.o util.o mapFile.o bufferedFile.o fileStar.o types.o

ATASEQ_OBJ = fasta_io.o char_io.o util.o mapFile.o bufferedFile.o fileStar.o types.o atac_seq.o

EXSEQ_OBJ = fasta_io.o char_io.o util.o mapFile.o bufferedFile.o fileStar.o types.o extract_seq.o

KMER_OBJ = fasta_io.o char_io.o util.o mapFile.o bufferedFile.o fileStar.o types.o kmer_count.o

KMERAN_OBJ = fasta_io.o char_io.o util.o mapFile.o bufferedFile.o fileStar.o types.o kmer_annotate.o

AACOMP_OBJ = fasta_io.o char_io.o util.o mapFile.o bufferedFile.o fileStar.o types.o aacomp.o

PROMW_OBJ = fasta_io.o char_io.o util.o mapFile.o bufferedFile.o fileStar.o types.o protein_mw.o

PEPMULT_OBJ = fasta_io.o char_io.o util.o mapFile.o bufferedFile.o fileStar.o types.o peptide_mult.o

AACOMPL_OBJ = fasta_io.o char_io.o util.o mapFile.o bufferedFile.o fileStar.o types.o aacomplookup.o

SVT_OBJ = sortedvector_test.o

TST_OBJ = hashtest.o hash.o bits.o char_io.o mapFile.o util.o

TANDEM_MATCH_OBJ = tandem_match.o pattern_alignment.o keyword_tree.o \
	 shift_and.o mapFile.o bufferedFile.o char_io.o fasta_io.o util.o \
	 pattern_match.o primer_alignment.o fileStar.o types.o select.o \
	 shift_and_inexact.o exact_bases.o  exact_halves.o filter_bitvec.o

XMER_OBJ = xmers.o shift_and.o mapFile.o bufferedFile.o char_io.o fasta_io.o util.o \
           pattern_match.o fileStar.o types.o select.o keyword_tree.o \
	   exact_bases.o exact_halves.o filter_bitvec.o shift_and_inexact.o \
	   rand_hash_table.o hash_table.o gs_hash_table.o primer_alignment.o \
	   pattern_alignment.o rlst.o

AVA_OBJ =  shift_and.o mapFile.o bufferedFile.o char_io.o fasta_io.o util.o \
           pattern_match.o fileStar.o types.o primer_alignment.o \
	   pattern_alignment.o allvall.o hash.o bits.o perfposht.o oligotm.o

PS_OBJ =   shift_and.o mapFile.o bufferedFile.o char_io.o fasta_io.o util.o \
           pattern_match.o fileStar.o types.o primer_alignment.o \
	   pattern_alignment.o pairscan.o hash.o bits.o perfposht.o

AVAD_OBJ =  shift_and.o mapFile.o bufferedFile.o char_io.o fasta_io.o util.o \
            pattern_match.o fileStar.o types.o primer_alignment.o \
	    pattern_alignment.o allvall_dump.o hash.o bits.o perfposht.o

AVATBM_OBJ = shift_and.o mapFile.o bufferedFile.o char_io.o fasta_io.o util.o \
             pattern_match.o fileStar.o types.o primer_alignment.o \
	     pattern_alignment.o allvall_tobm.o hash.o bits.o perfposht.o

AVAMER_OBJ = util.o types.o allvall_merge.o 

GS_OBJ = genome_simulation.o mapFile.o bufferedFile.o char_io.o fasta_io.o util.o \
           fileStar.o types.o 

SS_OBJ = solid_simulation.o mapFile.o bufferedFile.o char_io.o fasta_io.o util.o \
           fileStar.o types.o 

PCR_MATCH_OBJ = pcr_match.o pattern_alignment.o keyword_tree.o		 \
	 shift_and.o mapFile.o bufferedFile.o char_io.o fasta_io.o util.o \
	 pattern_match.o primer_alignment.o fileStar.o types.o select.o \
	 shift_and_inexact.o sts_io.o exact_bases.o exact_halves.o \
	 filter_bitvec.o hash_table.o gs_hash_table.o rand_hash_table.o rlst.o

PRUN_OBJ = fasta_io.o char_io.o util.o mapFile.o bufferedFile.o fileStar.o \
	 types.o polyrun.o

WG_OBJ = char_io.o util.o mapFile.o bufferedFile.o fileStar.o \
	 types.o walk_graph.o netflo.o keyword_tree.o pattern_match.o 

CCSBHG_OBJ = char_io.o util.o mapFile.o bufferedFile.o fileStar.o \
	 types.o cannon_csbh_graph.o keyword_tree.o pattern_match.o \
	 word_graph.o netflo.o

SA_OBJ = char_io.o util.o mapFile.o bufferedFile.o fileStar.o \
	 types.o solid_assembly.o keyword_tree.o pattern_match.o \
	 word_graph.o netflo.o

CSBHGAN_OBJ = char_io.o util.o mapFile.o bufferedFile.o fileStar.o \
	 types.o csbh_annotate.o keyword_tree.o pattern_match.o \
	 word_graph.o netflo.o

IDX_OBJ = Indexer.o rl_index.o

IDXAA_OBJ = IndexerAA.o rl_index.o

RLWG_OBJ = WordGraph.o rl_index.o

ST_OBJ = suftree.o util.o mapFile.o

XS_OBJ = Xspace.o rl_index.o

XSL_OBJ = XspaceLo.o xspacefsm.o rl_index.o

all: $(PROGS)
	 if [ "$(PLATFORM)" != "" ]; then \
	   if [ ! -d "$(PLATFORM)" ]; then \
	     mkdir $(PLATFORM); \
	   fi; \
	   cp -f $(PROGS) $(PLATFORM); \
	 fi

keyword_tree$(E): $(KEYWORD_TREE_OBJ)
	 $(CXX) $(LDFLAGS) -o $@ $(KEYWORD_TREE_OBJ) $(LDLIBS)

exact_match$(E): $(EXACT_MATCH_OBJ) $(PGEN)/primegen.a
	 $(CXX) $(LDFLAGS) -o $@ $(EXACT_MATCH_OBJ) $(PGEN)/primegen.a $(LDLIBS)


rand_hash_table.o:   INCLUDEPATH += -I$(PGEN)
gs_hash_table.o:     INCLUDEPATH += -I$(PGEN)

primer_match.o: EXTRACXXFLAGS += -DPRIMER3TM
primer_match$(E): $(PRIMER_MATCH_OBJ) $(PGEN)/primegen.a
	 $(CXX) $(LDFLAGS) -o $@ $(PRIMER_MATCH_OBJ) $(PGEN)/primegen.a $(LDLIBS)

xmers.o:   INCLUDEPATH += -I$(PGEN)
xmers$(E): $(XMER_OBJ) $(PGEN)/primegen.a
	 $(CXX) $(LDFLAGS) -o $@ $(XMER_OBJ) $(PGEN)/primegen.a $(LDLIBS)

allvall.o: EXTRACXXFLAGS += -DPRIMER3TM
allvall$(E): $(AVA_OBJ)
	 $(CXX) $(LDFLAGS) -o $@ $(AVA_OBJ) $(LDLIBS)

pairscan$(E): $(PS_OBJ)
	 $(CXX) $(LDFLAGS) -o $@ $(PS_OBJ) $(LDLIBS)

allvall_dump$(E): $(AVAD_OBJ)
	 $(CXX) $(LDFLAGS) -o $@ $(AVAD_OBJ) $(LDLIBS)

allvall_tobm$(E): $(AVATBM_OBJ)
	 $(CXX) $(LDFLAGS) -o $@ $(AVATBM_OBJ) $(LDLIBS)

allvall_merge$(E): $(AVAMER_OBJ)
	 $(CXX) $(LDFLAGS) -o $@ $(AVAMER_OBJ) $(LDLIBS)

genome_simulation$(E): $(GS_OBJ)
	 $(CXX) $(LDFLAGS) -o $@ $(GS_OBJ) $(LDLIBS)

solid_simulation$(E): $(SS_OBJ)
	 $(CXX) $(LDFLAGS) -o $@ $(SS_OBJ) $(LDLIBS)

inexact_match$(E): $(INEXACT_MATCH_OBJ) $(PGEN)/primegen.a
	 $(CXX) $(LDFLAGS) -o $@ $(INEXACT_MATCH_OBJ) $(PGEN)/primegen.a $(LDLIBS)

compress_seq$(E): $(COMPRESS_SEQ_OBJ)
	 $(CXX) $(LDFLAGS) -o $@ $(COMPRESS_SEQ_OBJ) $(LDLIBS)

chario$(E): $(CHARIO_OBJ)
	 $(CXX) $(LDFLAGS) -o $@ $(CHARIO_OBJ) $(LDLIBS)

sortedvector_test$(E): $(SVT_OBJ)
	 $(CXX) $(LDFLAGS) -o $@ $(SVT_OBJ) $(LDLIBS)

hashtest$(E): $(TST_OBJ)
	 $(CXX) $(LDFLAGS) -o $@ $(TST_OBJ) $(LDLIBS)

# atac_seq.o: INCLUDEPATH += -I$(ATACINC)
# atac_seq$(E): LIBRARYPATH += -L$(ATACLIB)
# atac_seq$(E): LDLIBS += -lATAC

atac_seq$(E): $(ATASEQ_OBJ)
	 $(CXX) $(LDFLAGS) -o $@ $(ATASEQ_OBJ) $(LDLIBS)

extract_seq$(E): $(EXSEQ_OBJ)
	 $(CXX) $(LDFLAGS) -o $@ $(EXSEQ_OBJ) $(LDLIBS)

kmer_count$(E): $(KMER_OBJ)
	 $(CXX) $(LDFLAGS) -o $@ $(KMER_OBJ) $(LDLIBS)

aacomp$(E): $(AACOMP_OBJ)
	 $(CXX) $(LDFLAGS) -o $@ $(AACOMP_OBJ) $(LDLIBS)

aacomplookup$(E): $(AACOMPL_OBJ)
	 $(CXX) $(LDFLAGS) -o $@ $(AACOMPL_OBJ) $(LDLIBS)

# walk_graph.o:   INCLUDEPATH += -I$(CS2)/$(PLATFORM)
# walk_graph$(E): LIBRARYPATH += -L$(CS2)/$(PLATFORM)
walk_graph.o:   INCLUDEPATH += -I$(CS2)
walk_graph$(E): LIBRARYPATH += -L$(CS2)
walk_graph$(E): LDLIBS += -lCS2

walk_graph$(E): $(WG_OBJ)
	 $(CXX) $(LDFLAGS) -o $@ $(WG_OBJ) $(LDLIBS)

# word_graph.o:          INCLUDEPATH += -I$(CS2)/$(PLATFORM)
# cannon_csbh_graph$(E): LIBRARYPATH += -L$(CS2)/$(PLATFORM)
word_graph.o:          INCLUDEPATH += -I$(CS2)
cannon_csbh_graph$(E): LIBRARYPATH += -L$(CS2)
cannon_csbh_graph$(E): LDLIBS += -lCS2

cannon_csbh_graph$(E): $(CCSBHG_OBJ) 
	 $(CXX) $(LDFLAGS) -o $@ $(CCSBHG_OBJ) $(LDLIBS)

solid_assembly$(E): LIBRARYPATH += -L$(CS2)
solid_assembly$(E): LDLIBS += -lCS2

solid_assembly$(E): $(SA_OBJ) 
	 $(CXX) $(LDFLAGS) -o $@ $(SA_OBJ) $(LDLIBS)

csbh_annotate$(E): LIBRARYPATH += -L$(CS2)
csbh_annotate$(E): LDLIBS += -lCS2

csbh_annotate$(E): $(CSBHGAN_OBJ) 
	 $(CXX) $(LDFLAGS) -o $@ $(CSBHGAN_OBJ) $(LDLIBS)

nrdb$(E): $(NRDB_OBJ)
	 $(CXX) $(LDFLAGS) -o $@ $(NRDB_OBJ) $(LDLIBS)

kmer_annotate$(E): $(KMERAN_OBJ)
	 $(CXX) $(LDFLAGS) -o $@ $(KMERAN_OBJ) $(LDLIBS)

tandem_match$(E): $(TANDEM_MATCH_OBJ)
	 $(CXX) $(LDFLAGS) -o_OBJ) $(LDLIBS)

pcr_match$(E): $(PCR_MATCH_OBJ) $(PGEN)/primegen.a
	 $(CXX) $(LDFLAGS) -o $@ $(PCR_MATCH_OBJ) $(PGEN)/primegen.a $(LDLIBS)

polyrun$(E): $(PRUN_OBJ)
	 $(CXX) $(LDFLAGS) -o $@ $(PRUN_OBJ) $(LDLIBS)

peptide_scan$(E): $(PEPSCAN_OBJ) $(PGEN)/primegen.a
	$(CXX) $(LDFLAGS) -o $@ $(PEPSCAN_OBJ) $(PGEN)/primegen.a $(LDLIBS)

protein_mw$(E): $(PROMW_OBJ)
	$(CXX) $(LDFLAGS) -o $@ $(PROMW_OBJ) $(LDLIBS)

peptide_mult$(E): $(PEPMULT_OBJ)
	$(CXX) $(LDFLAGS) -o $@ $(PEPMULT_OBJ) $(LDLIBS)

Indexer$(E): $(IDX_OBJ)
	$(CXX) $(LDFLAGS) -o $@ $(IDX_OBJ) $(LDLIBS)

IndexerAA$(E): $(IDXAA_OBJ)
	$(CXX) $(LDFLAGS) -o $@ $(IDXAA_OBJ) $(LDLIBS)

WordGraph$(E): $(RLWG_OBJ)
	$(CXX) $(LDFLAGS) -o $@ $(RLWG_OBJ) $(LDLIBS)

suftree$(E): $(ST_OBJ)
	$(CXX) $(LDFLAGS) -o $@ $(ST_OBJ) $(LDLIBS)

Xspace$(E): $(XS_OBJ)
	$(CXX) $(LDFLAGS) -o $@ $(XS_OBJ) $(LDLIBS)

XspaceLo$(E): $(XSL_OBJ)
	$(CXX) $(LDFLAGS) -o $@ $(XSL_OBJ) $(LDLIBS)

.PHONY: depend clean realclean dbclean test verify

depend: Makefile.d

clean:
	 - rm -f *% *~ *.o *.bak core gmon.out *% *.d *.lst *.org0 *.rpl0 primer_match_z.cc
	 - rm -rf cxx_repository
	 - rm -rf tempinc
	 - rm -rf test

realclean: clean
	 - rm -f $(PROGS) a.out *.exe

dbclean: 
	 - rm -f db/*.{seq,idx,hdr,tbl,sqn,tbz,sqz}

Makefile.d: $(SOURCES)
	 - touch Makefile.d; makedepend -- -Y -- -fMakefile.d $(SOURCES) >/dev/null 2>&1 

test: primer_match compress_seq
	 - rm -rf test
	 ./testscript.sh

verify: $(PROGS)
	 ./testscript.sh -v

%.dist: %.files %.Makefile
	- rm -rf $*.tar.gz 
	if [ ! -d tmp ] ; then \
		 mkdir tmp; \
	fi
	- rm -rf tmp/$*
	mkdir tmp/$* tmp/$*/src tmp/$*/bin
	for f in `cat $*.files` ; do \
		echo "Adding to $*: $$f" ; \
		case "$$f" in \
			*.c)  cat $*.LICENSE.code $$f > tmp/$*/src/$$f;; \
			*.cc) cat $*.LICENSE.code $$f > tmp/$*/src/$$f;; \
			*.h)  cat $*.LICENSE.code $$f > tmp/$*/src/$$f;; \
			*.t)  cat $*.LICENSE.code $$f > tmp/$*/src/$$f;; \
			*)    cat $$f > tmp/$*/src/$$f;; \
		esac; \
	done
	@ echo "Adding to $*: $*.Makefile"
	cp -f $*.Makefile tmp/$*/src/Makefile
	@ echo "Adding to $*: $*.LICENSE.txt"
	cp -f $*.LICENSE.txt tmp/$*/src/LICENSE.txt
	( cd tmp/$*/src; make ; mv -f $(PLATFORM)/* ../bin; make realclean; rm -rf $(PLATFORM) )
	tar -cvf - -C tmp $* | gzip -9 -c > $*.tar.gz

include Makefile.d
# DO NOT DELETE
