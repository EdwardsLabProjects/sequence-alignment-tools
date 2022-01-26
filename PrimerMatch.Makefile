
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

PARAM=-DNOPRIMEGEN

CXXFLAGS = $(EXTRACXXFLAGS) $(WARN) $(OPTIMIZE) $(DEBUG) $(PROFILE) $(INCLUDEPATH) $(PARAM)
LDFLAGS =  $(EXTRALDFLAGS) $(WARN) $(OPTIMIZE) $(DEBUG) $(PROFILE) $(LIBRARYPATH) $(PARAM)
LDLIBS = -lz -lm 
CFLAGS = $(CXXFLAGS)

PROGS = compress_seq$(E) primer_match$(E) pcr_match$(E) 

COMPRESS_SEQ_OBJ = compress_seq.o fasta_io.o char_io.o util.o	\
mapFile.o bufferedFile.o fileStar.o types.o

PRIMER_MATCH_OBJ = primer_match.o pattern_alignment.o keyword_tree.o	\
shift_and.o mapFile.o bufferedFile.o char_io.o fasta_io.o util.o	\
pattern_match.o primer_alignment.o fileStar.o types.o select.o		\
shift_and_inexact.o exact_bases.o exact_halves.o filter_bitvec.o	\
sts_io.o hash_table.o rlst.o

PCR_MATCH_OBJ = pcr_match.o pattern_alignment.o keyword_tree.o		\
shift_and.o mapFile.o bufferedFile.o char_io.o fasta_io.o util.o	\
pattern_match.o primer_alignment.o fileStar.o types.o select.o		\
shift_and_inexact.o sts_io.o exact_bases.o exact_halves.o		\
filter_bitvec.o hash_table.o rlst.o

all: $(PROGS)
	 if [ "$(PLATFORM)" != "" ]; then \
	   if [ ! -d "$(PLATFORM)" ]; then \
	     mkdir $(PLATFORM); \
	   fi; \
	   cp -f $(PROGS) $(PLATFORM); \
	 fi

compress_seq$(E): $(COMPRESS_SEQ_OBJ)
	 $(CXX) $(LDFLAGS) -o $@ $(COMPRESS_SEQ_OBJ) $(LDLIBS)

primer_match$(E): $(PRIMER_MATCH_OBJ)
	 $(CXX) $(LDFLAGS) -o $@ $(PRIMER_MATCH_OBJ) $(LDLIBS)

pcr_match$(E): $(PCR_MATCH_OBJ)
	 $(CXX) $(LDFLAGS) -o $@ $(PCR_MATCH_OBJ) $(LDLIBS)

.PHONY: depend clean realclean

depend: Makefile.d

clean:
	 - rm -f *% *~ *.o *.bak core gmon.out *% *.d *.lst *.org0 *.rpl0
	 - rm -rf cxx_repository
	 - rm -rf tempinc
	 - rm -rf test

realclean: clean
	 - rm -f $(PROGS) a.out *.exe

Makefile.d: $(SOURCES)
	 - touch Makefile.d; makedepend -- -Y -- -fMakefile.d $(SOURCES)>/dev/null 2>&1

include Makefile.d
