#!/bin/sh

if [ "$1" = "-v" ]; then
  VERIFY=1;
else
  VERIFY=0;
fi

if [ ! -d test ]; then 
    if [ $VERIFY -ne 1 ]; then
	mkdir test
    else 
	echo "Can't cd to test! Abort!";
	exit 1;
    fi
fi

cd test

if [ $VERIFY -ne 1 ]; then
    if [ ! -e pri_non_hum_nt.fasta ]; then
	ln -s ../db/pri_non_hum_nt.fasta
    fi
    if [ ! -e db ]; then
	ln -s ../db
    fi
    rm -f primer_match*
    rm -f compress_seq*
    cp ../primer_match ../primer_match.exe .
    cp ../compress_seq ../compress_seq.exe .
    PM=./primer_match
    CS=./compress_seq
else
    PM="echo ./primer_match"
    CS="echo ./compress_seq"
fi

CKSUMPROG=/usr/bin/cksum

checkident() {
    JIDSET=1;
    for a in $IDSET; do
	KIDSET=1;
	for b in $IDSET; do
	    if [ $JIDSET -lt $KIDSET ]; then 
		if cmp $a $b; then 
		    true;
		else
		    true # exit 1;
		fi
	    fi
	    KIDSET=`expr $KIDSET + 1`
	done
	JIDSET=`expr $JIDSET + 1`
    done
};

checkcksum() {
    JCKSUM=1;
    FILECKSUMS="";
    for a in $CKSUMDATA; do
	case $JCKSUM in 
	    1) FILENAME=$a;;
	    *) FILECKSUMS="$FILECKSUMS $a";;
	esac
	JCKSUM=`expr $JCKSUM + 1`
    done
    CKSUM=`$CKSUMPROG $FILENAME | awk '{print $1}'`
    VALIDCKSUM=0;
    for a in $FILECKSUMS; do 
	if [ $a = $CKSUM ]; then
	    VALIDCKSUM=1;
	fi
    done
    if [ $VALIDCKSUM -ne 1 ]; then
	echo "Cksum of $FILENAME is $CKSUM, not in $FILECKSUMS";
	# exit 1;
    fi
}

checklines() {
    JCKSUM=1;
    FILELINES="";
    for a in $LINESDATA; do
	case $JCKSUM in 
	    1) FILENAME=$a;;
	    *) FILELINES="$FILELINES $a";;
	esac
	JCKSUM=`expr $JCKSUM + 1`
    done
    LINES=`wc -l $FILENAME | awk '{print $1}'`
    VALIDLINES=0;
    for a in $FILELINES; do 
	if [ $a = $LINES ]; then
	    VALIDLINES=1;
	fi
    done
    if [ $VALIDLINES -ne 1 ]; then
	echo "File line count of $FILENAME is $LINES, not in $FILELINES";
	# exit 1;
    fi
}

checksize() {
    JCKSUM=1;
    FILESIZES="";
    for a in $SIZEDATA; do
	case $JCKSUM in 
	    1) FILENAME=$a;;
	    *) FILESIZES="$FILESIZES $a";;
	esac
	JCKSUM=`expr $JCKSUM + 1`
    done
    SIZE=`wc -c $FILENAME | awk '{print $1}'`
    VALIDSIZE=0;
    for a in $FILESIZES; do 
	if [ $a = $SIZE ]; then
	    VALIDSIZE=1;
	fi
    done
    if [ $VALIDSIZE -ne 1 ]; then
	echo "File size of $FILENAME is $SIZE, not in $FILESIZES";
	# exit 1;
    fi
}

function cmd() {
  echo "$@" 1>&2;
  "$@"
}

ALLALIGN='-A %h\t%H\t%f\t%s\t%e\t%5\t%3\t%S\t%E\t%i\t%d\t%p\t%q\t%Q\t%t\t%T\t%A\t%r\t%R\t%%\n'
ALLCOUNT='-C %i\t%p\t%q\t%r\t%R\t%c\t%C\t%+\t%%\n'
MINOUT='-A %h\t%i\t%r\t%d\n -C %i\t%r\t%c\t%C\t\n'
K0="";
K1="-k 1"
KK1="-K 1"
K2="-k 2"
O38="-3 8"
R="-r"
DB="-i pri_non_hum_nt.fasta"
PAT1="-P db/pat.txt"
PAT2="-F db/pat.fasta"


echo "Checking various forms of pattern/primer input..."

# 1.1
cmd $PM $DB $PAT1 $K0 $ALLALIGN $ALLCOUNT -o out.1.1
# 1.2
cmd $PM $DB $PAT2 $K0 $ALLALIGN $ALLCOUNT -o out.1.2
# 1.3
cmd $PM $DB -p "AGAAGCGAGTTCT CGCCAGCAGAGTT TTTTCTGAGAATCAAG CTATTGATAAGGGAGTGC ATGGCGGTTTTGTCGAA AAGAAAAGGGGGAAA TCATGAAGTAAAC TTGGCTGCTGCCCCCAG AGAAAAGGGGGAAA CTATTGATAAGGGAGTG" $K0 $ALLALIGN $ALLCOUNT -o out.1.3

IDSET="out.1.1 out.1.2 out.1.3";
checkident;

# 1.4
cmd $PM $DB $PAT1 $K1 $ALLALIGN $ALLCOUNT -o out.1.4
# 1.5
cmd $PM $DB $PAT2 $K1 $ALLALIGN $ALLCOUNT -o out.1.5
# 1.6
cmd $PM $DB -p "AGAAGCGAGTTCT CGCCAGCAGAGTT TTTTCTGAGAATCAAG CTATTGATAAGGGAGTGC ATGGCGGTTTTGTCGAA AAGAAAAGGGGGAAA TCATGAAGTAAAC TTGGCTGCTGCCCCCAG AGAAAAGGGGGAAA CTATTGATAAGGGAGTG" $K1 $ALLALIGN $ALLCOUNT -o out.1.6

IDSET="out.1.4 out.1.5 out.1.6";
checkident

# 1.4a
cmd $PM $DB $PAT1 $KK1 $ALLALIGN $ALLCOUNT -o out.1.4a
# 1.5a
cmd $PM $DB $PAT2 $KK1 $ALLALIGN $ALLCOUNT -o out.1.5a
# 1.6a
cmd $PM $DB -p "AGAAGCGAGTTCT CGCCAGCAGAGTT TTTTCTGAGAATCAAG CTATTGATAAGGGAGTGC ATGGCGGTTTTGTCGAA AAGAAAAGGGGGAAA TCATGAAGTAAAC TTGGCTGCTGCCCCCAG AGAAAAGGGGGAAA CTATTGATAAGGGAGTG" $KK1 $ALLALIGN $ALLCOUNT -o out.1.6a

IDSET="out.1.4a out.1.5a out.1.6a";
checkident

# 1.7
cmd $PM $DB $PAT1 $K2 $ALLALIGN $ALLCOUNT -o out.1.7
# 1.8
cmd $PM $DB $PAT2 $K2 $ALLALIGN $ALLCOUNT -o out.1.8
# 1.9
cmd $PM $DB -p "AGAAGCGAGTTCT CGCCAGCAGAGTT TTTTCTGAGAATCAAG CTATTGATAAGGGAGTGC ATGGCGGTTTTGTCGAA AAGAAAAGGGGGAAA TCATGAAGTAAAC TTGGCTGCTGCCCCCAG AGAAAAGGGGGAAA CTATTGATAAGGGAGTG" $K2 $ALLALIGN $ALLCOUNT -o out.1.9

IDSET="out.1.7 out.1.8 out.1.9";
checkident

echo "...done."

echo "Checking I/O possibilities..."

# 2.0
cmd $CS $DB -n true -z true -C false

CKSUMDATA="pri_non_hum_nt.fasta.hdr 234836286 797228482"
checkcksum;
SIZEDATA="pri_non_hum_nt.fasta.hdr 1294825"
checksize;

# CKSUMDATA="pri_non_hum_nt.fasta.idx 3266324250 2226310952 3721755287"
# checkcksum;
# SIZEDATA="pri_non_hum_nt.fasta.idx 350183 350433"
# checksize;

CKSUMDATA="pri_non_hum_nt.fasta.seq 2763003605 2832676483"
checkcksum;
SIZEDATA="pri_non_hum_nt.fasta.seq 15948018"
checksize;

CKSUMDATA="pri_non_hum_nt.fasta.sqn 2779833989 309891982"
checkcksum;
SIZEDATA="pri_non_hum_nt.fasta.sqn 15948018"
checksize;

CKSUMDATA="pri_non_hum_nt.fasta.sqz 312782519 3216638554"
checkcksum;
SIZEDATA="pri_non_hum_nt.fasta.sqz 7974016"
checksize;

CKSUMDATA="pri_non_hum_nt.fasta.tbl 4093313041 2419063880"
checkcksum;
SIZEDATA="pri_non_hum_nt.fasta.tbl 15"
checksize;

CKSUMDATA="pri_non_hum_nt.fasta.tbz 4093313041 2419063880"
checkcksum;
SIZEDATA="pri_non_hum_nt.fasta.tbz 15"
checksize;

# 2.1
cmd $PM $DB $PAT2 $K0 $ALLALIGN $ALLCOUNT -o out.2.1
# 2.2
cmd $PM $DB $PAT2 $K0 $ALLALIGN $ALLCOUNT -B -o out.2.2
# 2.3
cmd $PM $DB $PAT2 $K0 $ALLALIGN $ALLCOUNT -D 1 -o out.2.3
# 2.4
cmd $PM $DB $PAT2 $K0 $ALLALIGN $ALLCOUNT -B -D 1 -o out.2.4
# 2.5
cmd $PM $DB $PAT2 $K0 $ALLALIGN $ALLCOUNT -D 2 -o out.2.5
# 2.6
cmd $PM $DB $PAT2 $K0 $ALLALIGN $ALLCOUNT -B -D 2 -o out.2.6
# 2.7
cmd $PM $DB $PAT2 $K0 $ALLALIGN $ALLCOUNT -D 3 -o out.2.7
# 2.8
cmd $PM $DB $PAT2 $K0 $ALLALIGN $ALLCOUNT -B -D 3 -o out.2.8
# 2.9
cmd $PM $DB $PAT2 $K0 $ALLALIGN $ALLCOUNT -D 4 -o out.2.9
# 2.10
cmd $PM $DB $PAT2 $K0 $ALLALIGN $ALLCOUNT -B -D 4 -o out.2.10

IDSET="out.1.1 out.2.1 out.2.2 out.2.3 out.2.4 out.2.5 out.2.6 out.2.7 out.2.8 out.2.9 out.2.10";
checkident;

# 3.1
cmd $PM $DB $PAT2 $K1 $ALLALIGN $ALLCOUNT -o out.3.1
# 3.2
cmd $PM $DB $PAT2 $K1 $ALLALIGN $ALLCOUNT -B -o out.3.2
# 3.3
cmd $PM $DB $PAT2 $K1 $ALLALIGN $ALLCOUNT -D 1 -o out.3.3
# 3.4
cmd $PM $DB $PAT2 $K1 $ALLALIGN $ALLCOUNT -B -D 1 -o out.3.4
# 3.5
cmd $PM $DB $PAT2 $K1 $ALLALIGN $ALLCOUNT -D 2 -o out.3.5
# 3.6
cmd $PM $DB $PAT2 $K1 $ALLALIGN $ALLCOUNT -B -D 2 -o out.3.6
# 3.7
cmd $PM $DB $PAT2 $K1 $ALLALIGN $ALLCOUNT -D 3 -o out.3.7
# 3.8
cmd $PM $DB $PAT2 $K1 $ALLALIGN $ALLCOUNT -B -D 3 -o out.3.8
# 3.9
cmd $PM $DB $PAT2 $K1 $ALLALIGN $ALLCOUNT -D 4 -o out.3.9
# 3.10
cmd $PM $DB $PAT2 $K1 $ALLALIGN $ALLCOUNT -B -D 4 -o out.3.10

IDSET="out.1.4 out.3.1 out.3.2 out.3.3 out.3.4 out.3.5 out.3.6 out.3.7 out.3.8 out.3.9 out.3.10";
checkident;

# 3a.1
cmd $PM $DB $PAT2 $K2 $ALLALIGN $ALLCOUNT -o out.3a.1
# 3a.2
cmd $PM $DB $PAT2 $K2 $ALLALIGN $ALLCOUNT -B -o out.3a.2
# 3a.3
cmd $PM $DB $PAT2 $K2 $ALLALIGN $ALLCOUNT -D 1 -o out.3a.3
# 3a.4
cmd $PM $DB $PAT2 $K2 $ALLALIGN $ALLCOUNT -B -D 1 -o out.3a.4
# 3a.5
cmd $PM $DB $PAT2 $K2 $ALLALIGN $ALLCOUNT -D 2 -o out.3a.5
# 3a.6
cmd $PM $DB $PAT2 $K2 $ALLALIGN $ALLCOUNT -B -D 2 -o out.3a.6
# 3a.7
cmd $PM $DB $PAT2 $K2 $ALLALIGN $ALLCOUNT -D 3 -o out.3a.7
# 3a.8
cmd $PM $DB $PAT2 $K2 $ALLALIGN $ALLCOUNT -B -D 3 -o out.3a.8
# 3a.9
cmd $PM $DB $PAT2 $K2 $ALLALIGN $ALLCOUNT -D 4 -o out.3a.9
# 3a.10
cmd $PM $DB $PAT2 $K2 $ALLALIGN $ALLCOUNT -B -D 4 -o out.3a.10

IDSET="out.1.7 out.3a.1 out.3a.2 out.3a.3 out.3a.4 out.3a.5 out.3a.6 out.3a.7 out.3a.8 out.3a.9 out.3a.10";
checkident;

echo "...done."

# rm -f out.1.* out.2.* out.3.*

echo "Checking primer indexing strategies..."

# 4.1
cmd $PM $DB $PAT2 $K0 $ALLALIGN $ALLCOUNT -o out.4.1
# 4.2
cmd $PM $DB $PAT2 $K0 $ALLALIGN $ALLCOUNT -D 1 -N 1 -o out.4.2
# 4.3
cmd $PM $DB $PAT2 $K0 $ALLALIGN $ALLCOUNT -D 1 -N 2 -o out.4.3
# 4.4
cmd $PM $DB $PAT2 $K0 $ALLALIGN $ALLCOUNT -D 1 -N 3 -o out.4.4
# 4.5
cmd $PM $DB $PAT2 $K0 $ALLALIGN $ALLCOUNT -D 1 -N 4 -o out.4.5
# 4.6
cmd $PM $DB $PAT2 $K0 $ALLALIGN $ALLCOUNT -D 2 -o out.4.6
# 4.7
cmd $PM $DB $PAT2 $K0 $ALLALIGN $ALLCOUNT -D 2 -N 1 -o out.4.7
# 4.8
cmd $PM $DB $PAT2 $K0 $ALLALIGN $ALLCOUNT -D 2 -N 2 -o out.4.8
# 4.9
cmd $PM $DB $PAT2 $K0 $ALLALIGN $ALLCOUNT -D 2 -N 3 -o out.4.9
# 4.10
cmd $PM $DB $PAT2 $K0 $ALLALIGN $ALLCOUNT -D 2 -N 4 -o out.4.10
# 4.11
cmd $PM $DB $PAT2 $K0 $ALLALIGN $ALLCOUNT -D 3 -N 1 -o out.4.11
# 4.12
cmd $PM $DB $PAT2 $K0 $ALLALIGN $ALLCOUNT -D 3 -N 2 -o out.4.12
# 4.13
cmd $PM $DB $PAT2 $K0 $ALLALIGN $ALLCOUNT -D 3 -N 3 -o out.4.13
# 4.14
cmd $PM $DB $PAT2 $K0 $ALLALIGN $ALLCOUNT -D 3 -N 4 -o out.4.14
# 4.15
cmd $PM $DB $PAT2 $K0 $ALLALIGN $ALLCOUNT -D 4 -o out.4.15
# 4.16
cmd $PM $DB $PAT2 $K0 $ALLALIGN $ALLCOUNT -D 4 -N 1 -o out.4.16
# 4.17
cmd $PM $DB $PAT2 $K0 $ALLALIGN $ALLCOUNT -D 4 -N 2 -o out.4.17
# 4.18
cmd $PM $DB $PAT2 $K0 $ALLALIGN $ALLCOUNT -D 4 -N 3 -o out.4.18
# 4.19
cmd $PM $DB $PAT2 $K0 $ALLALIGN $ALLCOUNT -D 4 -N 4 -o out.4.19

IDSET="out.4.1 out.4.2 out.4.3 out.4.4 out.4.5 out.4.6 out.4.7 out.4.8 out.4.9 out.4.10 out.4.11 out.4.12 out.4.13 out.4.14 out.4.15 out.4.16 out.4.17 out.4.18 out.4.19"
checkident;

# rm -f out.4.*

# 5.1
cmd $PM $DB $PAT2 $K1 $MINOUT | sort > out.5.1
# 5.2
cmd $PM $DB $PAT2 $K1 $MINOUT -N 5 | sort > out.5.2
# 5.3
cmd $PM $DB $PAT2 $K1 $MINOUT -N 10 | sort > out.5.3
# 5.4
cmd $PM $DB $PAT2 $K1 $MINOUT -N 11 | sort > out.5.4
# 5.5
cmd $PM $DB $PAT2 $K1 $MINOUT -N 12 | sort > out.5.5
# 5.6
cmd $PM $DB $PAT2 $K1 $MINOUT -N 13 | sort > out.5.6

# 5.1
cmd $PM $DB $PAT2 $KK1 $MINOUT  | sort > out.5.1a
# 5.2
cmd $PM $DB $PAT2 $KK1 $MINOUT -N 5 | sort > out.5.2a
# 5.3
cmd $PM $DB $PAT2 $KK1 $MINOUT -N 10 | sort > out.5.3a
# 5.4
cmd $PM $DB $PAT2 $KK1 $MINOUT -N 11 | sort > out.5.4a
# 5.5
cmd $PM $DB $PAT2 $KK1 $MINOUT -N 12 | sort > out.5.5a
# 5.6
cmd $PM $DB $PAT2 $KK1 $MINOUT -N 13 | sort > out.5.6a

# 5.1
cmd $PM $DB $PAT2 $K1 $O38 $R $MINOUT | sort > out.5.1b
# 5.2
cmd $PM $DB $PAT2 $K1 $O38 $R $MINOUT -N 5 | sort > out.5.2b
# 5.3
cmd $PM $DB $PAT2 $K1 $O38 $R $MINOUT -N 6 | sort > out.5.3b
# 5.4
cmd $PM $DB $PAT2 $K1 $O38 $R $MINOUT -N 7 | sort > out.5.4b
# 5.5
cmd $PM $DB $PAT2 $K1 $O38 $R $MINOUT -N 8 | sort > out.5.5b
# 5.6
cmd $PM $DB $PAT2 $K1 $O38 $R $MINOUT -N 9 | sort > out.5.6b
# 5.3
cmd $PM $DB $PAT2 $K1 $O38 $R $MINOUT -N 10 | sort > out.5.7b
# 5.4
cmd $PM $DB $PAT2 $K1 $O38 $R $MINOUT -N 11 | sort > out.5.8b
# 5.5
cmd $PM $DB $PAT2 $K1 $O38 $R $MINOUT -N 12 | sort > out.5.9b
# 5.6
cmd $PM $DB $PAT2 $K1 $O38 $R $MINOUT -N 13 | sort > out.5.10b

# 5.1
cmd $PM $DB $PAT2 $K2 $O38 $R $MINOUT | sort > out.5.1c
# 5.2
cmd $PM $DB $PAT2 $K2 $O38 $R $MINOUT -N 5 | sort > out.5.2c
# 5.3
cmd $PM $DB $PAT2 $K2 $O38 $R $MINOUT -N 6 | sort > out.5.3c
# 5.4
cmd $PM $DB $PAT2 $K2 $O38 $R $MINOUT -N 7 | sort > out.5.4c
# 5.5
cmd $PM $DB $PAT2 $K2 $O38 $R $MINOUT -N 8 | sort > out.5.5c
# 5.6
cmd $PM $DB $PAT2 $K2 $O38 $R $MINOUT -N 9 | sort > out.5.6c

IDSET="out.5.1 out.5.2 out.5.4 out.5.5 out.5.6"
checkident;

IDSET="out.5.1a out.5.2a out.5.4a out.5.5a out.5.6a"
checkident;

IDSET="out.5.1b out.5.2b out.5.3b out.5.4b out.5.5b out.5.6b out.5.7b out.5.8b out.5.9b out.5.10b"
checkident;

IDSET="out.5.1c out.5.2c out.5.3c out.5.4c out.5.5c out.5.6c"
checkident;

echo "...done."

# rm -f out.5.*

echo "Checking for previous known bugs..."

#
# This is the bug for "large initial/final exact position" constraints
# 6.1
cmd $PM $DB -p ATCCTTTTCAGCACTTTTTCT $K1 -s 15 $ALLALIGN $ALLCOUNT -o out.6.1

# Checksums: (aix-xlc);(alpha-cxx);(linux-gcc);(cygwin-gcc)
CKSUMDATA="out.6.1 1750922134 215977123 1750922134 123444061";
checkcksum;
# Sizedata: (aix-xlc,alpha-cxx,linux-gcc);(cygwin-gcc)
SIZEDATA="out.6.1 1743 1750";
checksize;

echo "...done."

# rm -f out.6.*

#
# Check wildcard behavior...

echo "Checking wildcard behaviour..."

# 7.1
cmd $PM $DB $PAT2 $K1 $ALLALIGN $ALLCOUNT -w -o out.7.1
# 7.2
cmd $PM $DB $PAT2 $K1 $ALLALIGN $ALLCOUNT -W -o out.7.2

# Checksums: (aix-xlc);(alpha-cxx);(linux-gcc);(cygwin-gcc)
CKSUMDATA="out.7.1 3957987070 2823332668 3957987070 2311392354";
checkcksum;
# Sizedata: (aix-xlc,alpha-cxx,linux-gcc);(cygwin-gcc)
SIZEDATA="out.7.1 85463 85813";
checksize;
# Checksums: (aix-xlc);(alpha-cxx);(linux-gcc);(cygwin-gcc)
CKSUMDATA="out.7.2 665225986 997893636 11 3103467215";
checkcksum;
# Sizedata: (aix-xlc,alpha-cxx,linux-gcc);(cygwin-gcc)
SIZEDATA="out.7.2 92174 92547";
checksize;

# rm -f out.7.*

echo "...done"
