#include "rl_index.h"
#include <ctype.h>

void index_elt::fprint(FILE *f) const {
  fprintf(f,">%s\n",defline.c_str());
  fprintf(f,"%llu %llu %llu %llu\n",
	  cstart,cstop,start,stop);
}

int index_list::iload_fasta(FILE *f,char fl) {
  if (!f) return -1;

  index_elt ind;
  ind.flags = fl;
  
  FILE_POSITION_TYPE cpos=0;
  FILE_POSITION_TYPE pos=0;
  int count = 0;

  for(int c=0; c != EOF; ) {
    c = fgetc(f);
    if (isspace(c))
      continue;
    else if (c == '>' || c == EOF) {

      if (count) {
	cpos = FTELL(f);
	if (cpos == -1) {
	  fprintf(stderr,"index_list::iload_fasta: unable to ftell\n");
	  return -1;
	}
	ind.cstop = cpos;
	ind.stop  = pos;
	push_back(ind);
      }
      
      ++count;

      ind.defline.resize(0);
      while((c=fgetc(f))!=EOF && c != '\n') {
	ind.defline += char(c);
      }
      if (c == EOF) break;

      cpos = FTELL(f);
      if (cpos == -1) {
	fprintf(stderr,"index_list::iload_fasta: unable to ftell\n");
	return -1;
      }
      ind.start = pos++;
      ind.cstart = cpos;
    }
    else {
      pos++;
    }
  }
  return 0;
}

int index_list::isave(FILE *f) const {
  if (!f) return -1;

  for(index_list::const_iterator i=this->begin(); i != this->end(); ++i) {
    i->fprint(f);
  }
  return 0;
}

int index_list::iload(FILE *f,char fl) {
  if (!f) return -1;

  index_elt ind;
  ind.flags = fl;

  int scanc;

  for(int c=0; c != EOF; ) {
    c = fgetc(f);
    if (isspace(c))
      continue;
    else if (c == '>' || c == EOF) {

      if (c == EOF) return 0;

      ind.defline.resize(0);
      while((c=fgetc(f))!=EOF && c != '\n') {
	ind.defline += char(c);
      }
      if (c == EOF) return -1;

      scanc = fscanf(f,"%llu %llu %llu %llu\n",
		     &ind.cstart,&ind.cstop,&ind.start,&ind.stop);

      if (scanc < 4) return -1;

      push_back(ind);
    }
    else {
      return -1;
    }
  }
  return 0;
}

int min_index_list::iload(FILE *f,char fl) {
  if (!f) return -1;

  min_index_elt ind;

  int scanc;

  for(int c=0; c != EOF; ) {
    c = fgetc(f);
    if (isspace(c))
      continue;
    else if (c == '>' || c == EOF) {

      if (c == EOF) return 0;

      while((c=fgetc(f))!=EOF && c != '\n');
      if (c == EOF) return -1;

      FILE_POSITION_TYPE d1,d2;

      scanc = fscanf(f,"%llu %llu %llu %llu\n",
		     &d1,&d2,&ind.start,&ind.stop);

      if (scanc < 4) return -1;

      push_back(ind);
    }
    else {
      return -1;
    }
  }
  return 0;
}

int sequence::sload_fasta(FILE *f) {
  const FILE_POSITION_TYPE len = stop - start;
  FILE_POSITION_TYPE pos=1;

  if (seek_fasta(f) < 0) {
    fprintf(stderr,"sequence::sload_fasta: unable to fseek\n");
    return -1;
  }
  for (int c=0; (c = fgetc(f)) != EOF && pos < len; ) {
    if (isspace(c)) {
      continue;
    }
    else if (c == '>') {
      fprintf(stderr,"sequence::sload_fasta: premature end of sequence\n");
      fprintf(stderr,"sequence::sload_fasta: was %ld long and should have been %ld\n",pos-1,stop-start);
      return -1;
    }
    else {
      chars[pos++] = c;
    }
  }
  if (pos < len) {
    fprintf(stderr,"sequence::sload_fasta: premature end of fasta file\n");
    return -1;
  }
  return 0;
}

int sequence::sload(FILE *f,bool ends) {
  const FILE_POSITION_TYPE len = stop - start+(ends?1:-1);

  if (FSEEK(f,start,SEEK_SET) < 0) {
    fprintf(stderr,"sequence::sload: unable to fseek\n");
    return -1;
  }
  size_t fread_res;
  if ( (fread_res=fread(chars,sizeof(char),len,f)) < size_t(len) ) {
    fprintf(stderr,"sequence::sload: premature end of sequence\n");
    fprintf(stderr,"sequence::sload: %lu < %lu\n",fread_res,size_t(len));
    return -1;
  }
  return 0;
}

int min_sequence::sload(FILE *f,bool ends) {
  const FILE_POSITION_TYPE len = stop - start+(ends?1:-1);

  if (FSEEK(f,start,SEEK_SET) < 0) {
    fprintf(stderr,"sequence::sload: unable to fseek\n");
    return -1;
  }
  size_t fread_res;
  if ( (fread_res=fread(chars,sizeof(char),len,f)) < size_t(len) ) {
    fprintf(stderr,"sequence::sload: premature end of sequence\n");
    fprintf(stderr,"sequence::sload: %lu < %lu\n",fread_res,size_t(len));
    return -1;
  }
  return 0;
}

int sequence::ssave(FILE *f,bool ends) const {
  const FILE_POSITION_TYPE len = stop - start+(ends?1:-1);

  if (FSEEK(f,start,SEEK_SET) < 0) {
    fprintf(stderr,"sequence::ssave: unable to fseek\n");
    return -1;
  }
  unsigned sz=0;
  if ((sz=fwrite(chars,sizeof(char),len,f)) < size_t(len)) {
    fprintf(stderr,"sequence::ssave: weird fwrite failure\n");
    fprintf(stderr,"sequence::ssave: length was %u and should have been %ld\n",sz,len);
    return -1;
  }
  return 0;
}
