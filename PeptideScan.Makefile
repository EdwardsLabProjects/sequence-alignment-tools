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

#
# Set the platform/machine that this'll be compiled on
#
# PLATFORM={alpha-cxx,alpha-gpp,linux-gpp,aix-xlc,aix-gpp,cygwin-gpp}
#

# PLATFORM=aix-gpp
# PLATFORM=aix-xlc
# PLATFORM=linux-gpp
# PLATFORM=alpha-cxx
# PLATFORM=alpha-gpp
# PLATFORM=cygwin-gpp

ifeq ($(strip $(PLATORM)),)
  CXX=cannot determine platform
endif

ifeq ($(strip $(PLATFORM)),aix-xlc) 
  CXX      = xlC -+ -q64 -D_LARGE_FILES -bhalt:8 -qnolm 
  OPTIMIZE = -O -qinline -qarch=auto -qtune=auto -qcache=auto -qmaxmem=-1 
endif

ifeq ($(strip $(PLATFORM)),aix-gpp) 
  CXX         = /usr/bin/g++ -maix64 -D_LARGE_FILES
  OPTIMIZE    = -O
endif

ifeq ($(strip $(PLATFORM)),alpha-cxx) 
  CXX           = cxx 
  OPTIMIZE      = -O
  EXTRACXXFLAGS = -DNO_STD_NAMESPACE
endif

ifeq ($(strip $(PLATFORM)),alpha-gpp) 
  CXX           = g++
  OPTIMIZE      = -O
endif

ifeq ($(strip $(PLATFORM)),linux-gpp) 
  CXX          = g++ -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE
  OPTIMIZE     = -O
endif

ifeq ($(strip $(PLATFORM)),cygwin-gpp) 
  CXX          = g++
  OPTIMIZE     = -O 
  E            = .exe
endif

#
# The rest of these options should be left as is...
#

CXXFLAGS = -DNOPRIMEGEN $(EXTRACXXFLAGS) $(WARN) $(OPTIMIZE) $(DEBUG) $(PROFILE) $(INCLUDEPATH) $(PARAM)
LDFLAGS =  -DNOPRIMEGEN $(EXTRALDFLAGS) $(WARN) $(OPTIMIZE) $(DEBUG) $(PROFILE) $(LIBRARYPATH) $(PARAM)
LDLIBS = -lm -lz

PROGS = compress_seq$(E) peptide_scan$(E)

COMPRESS_SEQ_OBJ = compress_seq.o fasta_io.o char_io.o util.o	\
mapFile.o bufferedFile.o fileStar.o types.o

PEPSCAN_OBJ = peptide_scan.o pattern_alignment.o keyword_tree.o		\
shift_and.o mapFile.o bufferedFile.o char_io.o fasta_io.o util.o	\
pattern_match.o fileStar.o types.o select.o shift_and_inexact.o		\
exact_bases.o primer_alignment.o exact_halves.o filter_bitvec.o		\
hash_table.o rlst.o

all: $(PROGS)
	 if [ "$(PLATFORM)" != "" ]; then \
	   if [ ! -d "$(PLATFORM)" ]; then \
	     mkdir $(PLATFORM); \
	   fi; \
	   cp -f $(PROGS) $(PLATFORM); \
	 fi

compress_seq$(E): $(COMPRESS_SEQ_OBJ)
	 $(CXX) $(LDFLAGS) -o $@ $(COMPRESS_SEQ_OBJ) $(LDLIBS)

peptide_scan$(E): $(PEPSCAN_OBJ)
	$(CXX) $(LDFLAGS) -o $@ $(PEPSCAN_OBJ) $(LDLIBS)

.PHONY: depend clean realclean

depend: Makefile.d

clean:
	 - rm -f *% *~ *.o *.bak core gmon.out *% *.d *.lst *.org0 *.rpl0
	 - rm -rf cxx_repository
	 - rm -rf tempinc
	 - rm -rf test

realclean: clean
	 - rm -f $(PROGS) a.out *.exe

Makefile.d:
	 - touch Makefile.d; makedepend -- -Y -- -fMakefile.d *.cc >/dev/null 2>&1 

include Makefile.d
