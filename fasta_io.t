// -*- c++ -*-

#ifndef _FASTA_IO_T_
#define _FASTA_IO_T_

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

#include <string>
#include <assert.h>
#include <ctype.h>
#include "types.h"
#include "sortedvector.t"
#include "char_io.h"
#include "char_io.t"

struct fasta_file_seq_params {
  fasta_file_seq_params() : 
    check_params(true), upper_case(true), eos_start(true), eos_char('\n'), translate(false), mapindex(0), frame(0), offset(0) {}
  bool check_params;
  bool upper_case;
  bool eos_start;
  char eos_char;
  bool translate;
  int  mapindex;
  int  frame;
  FILE_POSITION_TYPE offset;
};

class Header {
public:
  Header(long unsigned int i=0,
	 CharacterProducer * cp=0,
	 FILE_POSITION_TYPE fpos=0,
	 FILE_POSITION_TYPE len=0) : _header("") {
    if (cp) {
      _header.reserve(len);
      cp->pos(fpos);
      for (FILE_POSITION_TYPE i=0;i<len;i++) {
      _header += cp->getch();
      }
    }
  }
  std::string const & header() const {
    return _header;
  }
private:
  std::string _header;
};

class Header_SI {
public:
  Header_SI(long unsigned int i=0,
	    CharacterProducer * cp=0,
	    FILE_POSITION_TYPE fpos=0,
	    FILE_POSITION_TYPE len=0) 
    : _index(i), _header(""), _sheader("") {
    if (cp) {
      _header.reserve(len);
      cp->pos(fpos);
      for (FILE_POSITION_TYPE i=0;i<len;i++) {
      _header += cp->getch();
      }
      const std::string ws(" \t");
      int p = anypos(_header,ws);
      if (p > 0) {
	_sheader = _header.substr(0,p);
      } else if (p == 0) {
	_sheader = "";
      } else {
	_sheader = _header;
      }
    }
  }
  std::string const & header() const {
    return _header;
  }
  std::string const & short_header() const {
    return _sheader;
  }
  long unsigned int const & index() const {
    return _index;
  }
private:
  long unsigned int _index;
  std::string _header;
  std::string _sheader;
};

class Lazy_Header_SI {
public:
  Lazy_Header_SI(long unsigned int i=0,
		 CharacterProducer * cp=0,
		 FILE_POSITION_TYPE fpos=0,
		 FILE_POSITION_TYPE len=0) 
    : _cp(cp), _fpos(fpos), _len(len), _index(i), 
      _header(""), _read(false), _sheader("") {}
  std::string const & header() const {
    if (!_read) read_header();
    return _header;
  }
  std::string const & short_header() const {
    if (!_read) read_header();
    return _sheader;
  }
  long unsigned int const & index() const {
    return _index;
  }
private:
  void read_header() const {
    if (_cp) {
      ((Lazy_Header_SI*)this)->_header.reserve(_len);
      _cp->pos(_fpos);
      for (FILE_POSITION_TYPE i=0;i<_len;i++) {
	((Lazy_Header_SI*)this)->_header += _cp->getch();
      }
      const std::string ws(" \t");
      int p = anypos(_header,ws);
      if (p > 0) {
	((Lazy_Header_SI*)this)->_sheader = _header.substr(0,p);
      } else if (p == 0) {
	((Lazy_Header_SI*)this)->_sheader = "";
      } else {
	((Lazy_Header_SI*)this)->_sheader = _header;
      }
      ((Lazy_Header_SI*)this)->_read = true;
    } 
  }
  CharacterProducer * _cp;
  FILE_POSITION_TYPE _fpos;
  FILE_POSITION_TYPE _len;
  long unsigned int _index;
  bool _read;
  std::string _header;
  std::string _sheader;
};

template<class H>
class FastaFile : public CharacterProducer {
public:
  FastaFile() {}
  virtual ~FastaFile() {};
  virtual H const & get_header_data(FILE_POSITION_TYPE)=0;
  virtual FILE_POSITION_TYPE get_seq_pos(FILE_POSITION_TYPE)=0;
  virtual bool is_subseq(FILE_POSITION_TYPE,FILE_POSITION_TYPE)=0;
  virtual void insert(FILE_POSITION_TYPE,H const &)=0;
  virtual bool fasta_pos(long unsigned int const &, FILE_POSITION_TYPE const &)=0;
  virtual void getbasepos(FILE_POSITION_TYPE p0, FILE_POSITION_TYPE& p1, int& f) const=0;
  virtual char getbasech()=0;
  virtual void mapto(char,char)=0;
};

template<class T, class H>
class SortSeqFastaFile : public FastaFile<H> {
private:
  T chars_;
  FILE_POSITION_TYPE offset_;
  sortedvector<FILE_POSITION_TYPE,H> headers_;
  long unsigned int headers_index_;
  typename sortedvector<FILE_POSITION_TYPE,H>::iterator headers_it_;
  bool first_lookup;
  void set_header_item(FILE_POSITION_TYPE pos) {
    try {
      if (first_lookup) {
	headers_it_ = headers_.locate_last_at_most(pos-1);
	headers_index_ = headers_.index(headers_it_);
	first_lookup = false;
      } else {
	headers_it_ = headers_.iter(headers_index_);
	headers_it_ = headers_.finger_locate_last_at_most(headers_it_,pos-1);
	headers_index_ = headers_.index(headers_it_);
      }
    }
    catch (typename sortedvector<FILE_POSITION_TYPE,H>::KeyOutOfRange &) {
      throw NoHeaderData();
    }
  }
  H null;
public:
  SortSeqFastaFile(std::string const & filename, 
		   fasta_file_seq_params const & ffp) 
    : chars_(filename.c_str(),ffp.eos_char,ffp.frame), offset_(ffp.offset), first_lookup(true) {};
  class NoHeaderData{};
  H const & get_header_data(FILE_POSITION_TYPE pos) {
    try {
      set_header_item(pos);
    }
    catch (NoHeaderData &) {
      return null;
    }
    return headers_it_->value();
  }
  FILE_POSITION_TYPE get_seq_pos(FILE_POSITION_TYPE pos) {
    try {
      set_header_item(pos);
    } 
    catch (NoHeaderData &) {
      return 0;
    }
    return (FILE_POSITION_TYPE)(pos-headers_it_->key());
  }
  bool is_subseq(FILE_POSITION_TYPE start, FILE_POSITION_TYPE end) {
    try {
      set_header_item(start+1);
      typename sortedvector<FILE_POSITION_TYPE,H>::iterator it1 = headers_it_;
      set_header_item(end);
      return (it1 == headers_it_);
    } 
    catch (NoHeaderData &) {
      return false;
    }
  }
  void insert(FILE_POSITION_TYPE p, H const & header) {
    headers_.append(p,header);
  }
  FILE_POSITION_TYPE lastkey() const {
    if (headers_.size() > 0) {
      return headers_.rbegin()->key();
    } else {
      return 0;
    }
  }
  void reserve(long unsigned int s) {
    headers_.reserve(s);
  }
  inline char getch() { return chars_.getch(); }
  inline char getcodonid() { return chars_.getcodonid(); }
  inline unsigned char getnch() { return chars_.getnch(); }
  inline char ch(unsigned char c) { return chars_.ch(c); }
  inline int nch(char c) { return chars_.nch(c); }
  inline unsigned int size() const { return chars_.size(); }
  inline bool eof() const { return chars_.eof(); }
  inline FILE_POSITION_TYPE pos() const { return chars_.pos()+offset_; }
  inline void pos(FILE_POSITION_TYPE p) { chars_.pos(p-offset_); }
  inline void reset() { chars_.reset(); }
  inline float progress() const { return chars_.progress(); }
  inline FILE_POSITION_TYPE length() const { return chars_.length(); }
  inline char const * c_str() const { return chars_.c_str(); }
  inline char const * filename() const { return chars_.filename(); }
  inline bool has_filename() const { return chars_.has_filename(); }
  bool fasta_pos(long unsigned int const & fe, 
		 FILE_POSITION_TYPE const & rp) {
    if (/* fe < 0 || */ fe >= headers_.size()) {
      return false;
    }
    if (fe+1 < headers_.size() && 
	headers_.iter(fe+1)->key() <= rp) {
      return false;
    }
    pos(headers_.iter(fe)->key() + rp);
    return true;
  }
  void getbasepos(FILE_POSITION_TYPE p0, FILE_POSITION_TYPE& p1, int& f) const {
    chars_.getbasepos(p0,p1,f);
  }
  char getbasech() {
    return chars_.getbasech();
  }
  void mapto(char f, char t) {
    chars_.mapto(f,t);
  }
};

template <class T, class H>
class IndexedFastaFile : public SortSeqFastaFile<T,H> {
private:
  FileStarChars _hdrcp;
  void check_fasta_file_params(std::string const & filename_, 
		    fasta_file_seq_params const & ffp) {
    std::string idbfn = change_extn(filename_,"idb");
    if (exist(idbfn)) {
#if ! defined(__CYGWIN__) and ! defined(__MINGW32__)
      ifstream idb(idbfn.c_str());
#else 
      ifstream idb(idbfn.c_str(),ios::binary);
#endif
      sortedvector<FILE_POSITION_TYPE,FILE_POSITION_TYPE> svind;
      sortedvector<FILE_POSITION_TYPE,FILE_POSITION_TYPE>::iterator svit;
      svind.bread(idb);
      svit=svind.begin();
      FILE_POSITION_TYPE last_value=svit->value();
      FILE_POSITION_TYPE last_key=svit->key();
      ++svit;

      if (ffp.eos_start && last_key == 0) {
	timestamp("Bad format for indexed sequence database.");
	timestamp("Parameter indicates EOS as first character, but first sequence starts at 0.");
	exit(1);
      }

      if ((!ffp.eos_start) && last_key == (ffp.translate?3:1)) {
	timestamp("Bad format for indexed sequence database.");
	timestamp("Parameter indicates no EOS as first character, but first sequence starts at 1.");
	exit(1);
      }

      if (last_key > (ffp.translate?3:1)) {
        timestamp("Bad format for indexed sequence database.");
	timestamp("First sequence starts at position > 1.")
	exit(1);
      }

      if (!ffp.eos_start) {
	this->pos(svit->key()-1);
      }
      char ch = this->getch();
      if (ch != ffp.eos_char) {
	timestamp("Bad format for indexed sequence database.");
	timestamp("EOS character mismatch.");
	fprintf(stderr,"From indexed sequence database: %c\nFrom primer_match config: %c\n",ch,ffp.eos_char);
	exit(1);
      }
      this->pos(0);

      if (ffp.upper_case && this->nch('a') >= 0) {
	timestamp("Bad format for indexed sequence database.");
	timestamp("Parameter indicates uppercase, but lowercase characters permitted.");
	exit(1);
      }

    } else {
      std::string idxfn = change_extn(filename_,"idx");
      ifstream idx(idxfn.c_str());
//       if (!idx.good()) {
// 	checkpoint;
// 	std::cerr << idxfn.c_str() << endl;
//       }

      long unsigned int counter=0;
      FILE_POSITION_TYPE start_position=0,end_position=0;
      FILE_POSITION_TYPE start_position1=0,end_position1=0;
      FILE_POSITION_TYPE start_position2=0,end_position2=0;

      idx >> counter >> start_position >> start_position1 >> start_position2;

//       std::cerr << counter << " "
// 		<< start_position << " "
// 		<< start_position1 << " "
// 		<< start_position2 << std::endl;

      idx >> counter >> end_position >> end_position1 >> end_position2;

//       std::cerr << counter << " "
// 		<< end_position << " "
// 		<< end_position1 << " "
// 		<< end_position2 << std::endl;

      if (ffp.eos_start && start_position1 == 0) {
	timestamp("Bad format for indexed sequence database.");
	timestamp("Parameter indicates EOS as first character, but first sequence starts at 0.");
	exit(1);
      }

      if (! ffp.eos_start && start_position1 == 1) {
	timestamp("Bad format for indexed sequence database.");
	timestamp("Parameter indicates no EOS as first character, but first sequence starts at 1.");
	exit(1);
      }

      if (start_position1 > (ffp.translate?3:1)) {
        timestamp("Bad format for indexed sequence database.");
	timestamp("First sequence starts at position > 1.")
	exit(1);
      }

      if (!ffp.eos_start) {
	this->pos(end_position1-1);
      }
      char ch = this->getch();
      if (ch != ffp.eos_char) {
	timestamp("Bad format for indexed sequence database.");
	timestamp("EOS character mismatch.");
	fprintf(stderr,"From indexed sequence database: %c\nFrom primer_match config: %c\n",ch,ffp.eos_char);
	exit(1);
      }
      this->pos(0);

      if (ffp.upper_case && this->nch('a') >= 0) {
	timestamp("Bad format for indexed sequence database.");
	timestamp("Parameter indicates uppercase, but lowercase characters permitted.");
	exit(1);
      }

    }
  }
  void index_headers(std::string const & filename_) {
    std::string idbfn = change_extn(filename_,"idb");
    if (exist(idbfn)) {
#if ! defined(__CYGWIN__) and ! defined(__MINGW32__)
      ifstream idb(idbfn.c_str());
#else 
      ifstream idb(idbfn.c_str(),ios::binary);
#endif
      sortedvector<FILE_POSITION_TYPE,FILE_POSITION_TYPE> svind;
      sortedvector<FILE_POSITION_TYPE,FILE_POSITION_TYPE>::iterator svit;
      svind.bread(idb);
      long unsigned int counter=1;
      this->reserve(svind.size());
      svit=svind.begin();
      FILE_POSITION_TYPE last_value=svit->value();
      FILE_POSITION_TYPE last_key=svit->key();
      ++svit;
      while (svit!=svind.end()) {
	this->insert(last_key,
	       H(counter,&_hdrcp,last_value,svit->value()-last_value-1));
	last_key = svit->key();
	last_value = svit->value();
	++counter;
	++svit;
      }
    } else {
      std::string idxfn = change_extn(filename_,"idx");
      ifstream idx(idxfn.c_str());
      long unsigned int counter=0;
      FILE_POSITION_TYPE start_position=0,end_position=0;
      FILE_POSITION_TYPE start_position1=0,end_position1=0;
      FILE_POSITION_TYPE start_position2=0,end_position2=0;
      idx >> counter >> start_position >> start_position1 >> start_position2;
      while (!idx.eof()) {
	idx >> counter >> end_position >> end_position1 >> end_position2;
	if (idx.eof()) break;
	this->insert(start_position1,
	       H(counter,&_hdrcp,
		 start_position,
		 end_position-start_position-1));
	start_position = end_position;
	start_position1 = end_position1;
	start_position2 = end_position2;
      }
    }
  }
public:
  IndexedFastaFile(std::string const & filename, bool headers, 
		   fasta_file_seq_params const & ffp) : 
    SortSeqFastaFile<T,H>(filename+".seq",ffp), 
    _hdrcp(change_extn(filename,"hdr").c_str(),ffp.eos_char) {
    if (ffp.check_params) {
      check_fasta_file_params(filename,ffp);
    }
    if (headers) {
      index_headers(filename);
    } 
  }
  ~IndexedFastaFile() {};
};

template <class T, class H>
class StreamedFastaFile : public SortSeqFastaFile<T,H> {
private:
  bool _store_headers;
  bool _eof;
  bool _uc;
  bool _eos_start;
  char _eos;
  unsigned int _stride;
  unsigned int _stride_inc;
  FILE_POSITION_TYPE _prev_line_begin;
  FILE_POSITION_TYPE _pos;
  long unsigned int _hcount;
  sortedvector<FILE_POSITION_TYPE,FILE_POSITION_TYPE> _filepos;
  long unsigned int _filepos_index;
  bool first_lookup;
  FILE_POSITION_TYPE _max_pos_seen;
  T _hdrcp;
public:
  StreamedFastaFile(std::string const & filename, bool sthead, 
		    fasta_file_seq_params const & ffp) : 
    SortSeqFastaFile<T,H>(filename,ffp), 
    _store_headers(sthead), _eof(false), _uc(ffp.upper_case), _eos(ffp.eos_char), 
    _eos_start(ffp.eos_start),
    _stride(0), _prev_line_begin(0), _pos(0), _hcount(1), first_lookup(true),
    _max_pos_seen(0), _hdrcp(filename.c_str(),ffp.eos_char) {}
  ~StreamedFastaFile() {};
  unsigned char getnch() {
    return (unsigned char) (this->getch());
  }
  char getch() {
    if (SortSeqFastaFile<T,H>::eof() && !_eof) {
      _eof = true;
      _pos++;
      // cerr << "1: Output char: " << _eos << " Position: " << _pos << " T::pos(): " << T::pos() << endl;
      return _eos;
    } else if (_eos_start && _pos == 0) {
      _pos++;
      return _eos;
    } else {
      char ch = SortSeqFastaFile<T,H>::getch();
      if (ch != '\n' && ch != '\r' && ch != ' ' && ch != '>') {
	_pos++;
	// cerr << "2: Output char: " << ch << " Position: " << _pos << " T::pos(): " << T::pos() << endl;
	if (_uc) {
	  return toupper(ch);
	} else {
	  return ch;
	}
      } else {
	bool seen_newline=(false||_pos==0);
	bool looped = false;
	while (ch == '\n' || ch == '\r' || ch == ' ') {
	  looped = true;
	  if (SortSeqFastaFile<T,H>::eof()) {
	    _eof = true;
	    _pos++;
	    // cerr << "3: Output char: " << _eos << " Position: " << _pos << " T::pos(): " << T::pos() << endl;
	    return _eos;
	  }
	  if (ch == '\n') seen_newline = true;
	  ch = SortSeqFastaFile<T,H>::getch();
	}
	if (!seen_newline && looped) {
	  checkpoint;
	  ostrstream ss0;
	  ss0 << "WARNING: Irregular FASTA file format at character " 
	      << SortSeqFastaFile<T,H>::pos()-2
	      << "." << ends;
	  std::string v0(ss0.str());
	  timestamp(v0.c_str());
	  ostrstream ss1;
	  ss1 << "WARNING: Results may be incorrect! Please use compress_seq." 
	      << ends;
	  std::string v1(ss1.str());
	  timestamp(v1.c_str());
	}
	if (ch == '>') {
	  // checkpoint;
	  FILE_POSITION_TYPE headerstart = SortSeqFastaFile<T,H>::pos();
	  while (!SortSeqFastaFile<T,H>::eof() && 
		 (ch=SortSeqFastaFile<T,H>::getch())!='\n' && ch!='\r');
	  FILE_POSITION_TYPE headerend = SortSeqFastaFile<T,H>::pos();
	  if (ch == '\r') {
	    ch = SortSeqFastaFile<T,H>::getch();
	    assert(ch == '\n');
	  }
	  if (_pos != (_eos_start?1:0)) {
	    _pos++;
	  }
	  if (_store_headers) {
	    try {
 	      // cerr << ":" << _hcount << " " << endl;
// 	      << headerstart << " "
// 	      << headerend-headerstart-1 << " "
// 	      << endl;	
	      // checkpoint;
	      // cerr << this->lastkey() << " " << _pos << endl;
	      if (this->lastkey() < _pos || _pos == (_eos_start?1:0)) {
		// checkpoint;
		// cerr << this->lastkey() << " " << _pos << endl;
		this->insert(_pos, H(_hcount,&_hdrcp,headerstart,
			       headerend-headerstart-1));
	      }
	    } 
	    catch (typename sortedvector<FILE_POSITION_TYPE,H>::BadAppend &) {
	      // checkpoint;
	      /* No error, silently tolerated */ ;
	    }
	  } 
	  try {
	    // checkpoint;
	    // cerr << _filepos.rbegin()->key() << " " << _pos << endl;
	    if (_pos == (_eos_start?1:0) || _filepos.rbegin()->key() < _pos) { 
	      _filepos.append(_pos,SortSeqFastaFile<T,H>::pos());
	      // cerr << _pos << " " << SortSeqFastaFile<T,H>::pos() << endl;
	      _hcount++;
	    }
	  }
	  catch (sortedvector<FILE_POSITION_TYPE,FILE_POSITION_TYPE>::BadAppend &) {
	    // checkpoint;
	    /* No error, silently tolerated */ ;	      
	  }
	  // checkpoint;
	  // cerr << _pos << " " << SortSeqFastaFile<T,H>::pos() << endl;
	  // checkpoint;
	  // cerr << _prev_line_begin << " " << _pos << endl;
	  if (_pos == (_eos_start?1:0)) {
	    ch = SortSeqFastaFile<T,H>::getch();
	    _prev_line_begin = _pos;
	    _pos++;
	    if (_uc) {
	      return toupper(ch);
	    } else {
	      return ch;
	    }	  
	  } else {
	    _prev_line_begin = _pos;
	    return _eos;
	  }
	} else {
	  if (_stride==0) {
	    sortedvector<FILE_POSITION_TYPE,FILE_POSITION_TYPE>::const_iterator filepos_it = _filepos.end()-1;
	    _stride = (_pos - filepos_it->key());
	    _stride_inc = (SortSeqFastaFile<T,H>::pos() - filepos_it->value()) - _stride-1;
	    // checkpoint;
	    // cerr << "_pos == " << _pos << endl;
	    // cerr << "_filepos.key(_filepos_it) == " << _filepos.key(_filepos_it) << endl;
	    // cerr << "_stride == " << _stride << endl;
	    // cerr << "_stride_inc == " << _stride_inc << endl;
	  } else {
	    if (_pos - _prev_line_begin != _stride && _prev_line_begin != 0) {
	      // checkpoint;
	      // cerr << _pos << " " << _prev_line_begin << " " << _stride << " " << _stride_inc << endl;
	      // checkpoint;
	      ostrstream ss0;
	      ss0 << "WARNING: Irregular FASTA file format at character " 
		  << SortSeqFastaFile<T,H>::pos()-2
		  << "." << ends;
	      std::string v0(ss0.str());
	      timestamp(v0.c_str());
	      ostrstream ss1;
	      ss1 << "WARNING: Results may be incorrect! Please use compress_seq." 
		  << ends;
	      std::string v1(ss1.str());
	      timestamp(v1.c_str());
	    }
	  }
	  // checkpoint;
	  // cerr << _prev_line_begin << " " << _pos << endl;
	  _prev_line_begin = _pos;
	  _pos++;
	  // cerr << "6: Output char: " << ch << " Position: " << _pos << " T::pos(): " << T::pos() << endl;
	  if (_uc) {
	    return toupper(ch);
	  } else {
	    return ch;
	  }
	}
      }
    }
  }
  bool eof() const {
    return _eof;
  }
  bool fasta_pos(long unsigned int const & fe, 
		 FILE_POSITION_TYPE const & rp) {
    // checkpoint;
    /* if (fe < 0) {
      return false;
    }
    */
    if (fe+1 < _filepos.size() && _filepos.iter(fe+1)->key() <= rp) {
      return false;
    }
    if (fe >= _filepos.size()) {
      // We haven't seen this fasta entry yet!
      // Read forward 100000 bases...
      // checkpoint;
      // cerr << pos() << endl;
      while (fe >= _filepos.size() && !eof()) {
	pos(pos() + 100000);
	// checkpoint;
	// cerr << pos() << endl;
      }
    }
    if (fe >= _filepos.size()) {
      // We tried! Return false...
      // checkpoint;
      return false;
    }
    pos(_filepos.iter(fe)->key() + rp);
    return true;
  }
  FILE_POSITION_TYPE pos() const {
    return _pos;
  }
  void pos(FILE_POSITION_TYPE p) {
    // checkpoint;
    // cerr << p << "!!" << _pos << " " << SortSeqFastaFile<T,H>::pos() << endl;
    if (_pos > _max_pos_seen) {
      _max_pos_seen = _pos;
    }
    if (p > _max_pos_seen) {
      // checkpoint;
      // set position to the beginning of the last entry...
      if (_filepos.size() != 0) {
	fasta_pos(_filepos.size()-1,0);
      }
      // and get characters until our position hits p...
      while (!eof() && _pos < p) {
	getnch();
      }
      if (eof()) return;
      // then proceed, since we'll have seen all the headers before
      // p... this could be made more efficient, but why? If the user
      // really wants this to be efficient, they should be using
      // compress_seq...
    }
    sortedvector<FILE_POSITION_TYPE,FILE_POSITION_TYPE>::const_iterator 
      filepos_it;
    if (first_lookup) {
      try {
	if (p > 0) {
	  filepos_it = _filepos.locate_last_at_most(p);
	} else {
	  filepos_it = _filepos.iter(0);
	}
	_filepos_index = _filepos.index(filepos_it);
      }
      catch (sortedvector<FILE_POSITION_TYPE,FILE_POSITION_TYPE>::KeyOutOfRange &) {
	// checkpoint;
	assert(0);
      }
      first_lookup = false;
    } else {
      try { 
	filepos_it = _filepos.iter(_filepos_index);
	filepos_it = _filepos.finger_locate_last_at_most(filepos_it,p);
	_filepos_index = _filepos.index(filepos_it);
      } 
      catch (sortedvector<FILE_POSITION_TYPE,FILE_POSITION_TYPE>::InvalidFinger &) {
	assert(0);
      }	
      catch (sortedvector<FILE_POSITION_TYPE,FILE_POSITION_TYPE>::KeyOutOfRange &) {
	assert(0);
      }
    }
    FILE_POSITION_TYPE seq_offset=p-filepos_it->key();
    // cerr << p << " " << seq_offset << endl;
  
    // cerr << "!!" << _pos << " " << SortSeqFastaFile<T,H>::pos() << endl;
    // if (_filepos.value(_filepos_it) < _filepos.key(_filepos_it)) {
    // cerr << _filepos << endl;
    // }
    // cerr << "!!" << _pos << " " << SortSeqFastaFile<T,H>::pos() << endl;
    if (_stride > 0) {
      seq_offset += (seq_offset/_stride)*_stride_inc;
    }
    // checkpoint;
    // cerr << "!!" << _pos << " " << SortSeqFastaFile<T,H>::pos() << endl;
    // cerr << _filepos.key(_filepos_it) << " " << _filepos.value(_filepos_it) << " " << seq_offset << endl;
    // checkpoint;
    // cerr << "!!" << _pos << " " << SortSeqFastaFile<T,H>::pos() << endl;
    SortSeqFastaFile<T,H>::pos(filepos_it->value()+seq_offset);
    _pos = p;
    _eof = false;
    _prev_line_begin=0;
    // checkpoint;
    // cerr << "!!" << _pos << " " << SortSeqFastaFile<T,H>::pos() << endl;
    // int i=0;
    // while (i<10) {
    // cerr << getch();
    // i++;
    // }
    // checkpoint;
    // cerr << "!!" << _pos << " " << SortSeqFastaFile<T,H>::pos() << endl;
    // SortSeqFastaFile<T,H>::pos(_filepos.value(_filepos_it)+seq_offset);
    // _pos = p;
    // checkpoint;
    // cerr << "!!" << _pos << " " << SortSeqFastaFile<T,H>::pos() << endl;
    // checkpoint;
  }
};


#endif
