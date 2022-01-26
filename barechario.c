
#include <stdio.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/mman.h>


typedef off_t FILE_POSITION_TYPE;
typedef off_t FILE_SIZE_TYPE;

#define timestamp(m) { time_t t;\
                       t = time(NULL); \
                       sbuf = asctime_r(localtime(&t),sbuf); \
                       fprintf(stderr,"[%24.*s] %s\n",(int)strlen(sbuf)-1,sbuf,m); \
                     }

#define tic {tictoctime_ = time(NULL);}

#define toc { time_t t; time(&t); \
              t-=tictoctime_; \
              fprintf(stderr,"Elapsed time: %d:%02d:%02d\n",t/3600,(t%3600)/60,t%60);\
            } 
#define tocassgn(tsec) { time_t t; \
                         time(&t); \
                         t-=tictoctime_; \
                         fprintf(stderr,"Elapsed time: %d:%02d:%02d\n",t/3600,(t%3600)/60,t%60); \
                         tsec = t;\
                       } 

int main(int argc,char *argv[]) {

  time_t tictoctime_=0;
  char *sbuf;
  char *filename;
  char ch;
  char *fileString;
  FILE_SIZE_TYPE mapsize=0;
  FILE_POSITION_TYPE pos=0;
  int f;
  struct stat  sb;
  time_t tsec;
  double mapsize1;
  double cps;
  unsigned int pagesize;

  sbuf=(char*)malloc(1024*sizeof(char));

  if (argc < 2) {
    fprintf(stderr,"Need filename as first command-line arguement\n");
    exit(1);
  }
  filename = strdup(argv[1]);
  
  f = open(filename, O_RDONLY);
  if ((f < 0) || (errno)) {
    fprintf(stderr, "Couldn't open '%s'\n", filename);
    perror("open");
    exit(1);
  }
  
  fstat(f, &sb);
  if (errno) {
    fprintf(stderr, "Couldn't stat '%s'\n", filename);
    perror("fstat\n");
    exit(1);
  }
  
  mapsize = sb.st_size;

  fileString = (char *)mmap(0, mapsize, PROT_READ, 
				  (MAP_FILE|MAP_SHARED),f, 0);
  
  if (errno) {
    fprintf(stderr, "Couldn't map '%s'\n", filename);
    perror("mmap");
    exit(1);
  }
  
  close(f);
  
  ch=0;

  pagesize=getpagesize();
  // fprintf(stderr,"Page size is %d\n",pagesize);

  timestamp("Start scan");
  tic;
  while (pos < mapsize) {
    if (fileString[pos] > ch) {
      ch = fileString[pos];
    }
    pos+=pagesize;
  }
  fprintf(stderr,"Max. character: %c\n",ch);
  timestamp("End scan");

  tocassgn(tsec);

  if (mapsize > 1024*1024*1024) {
    mapsize1 = ((double)mapsize)/(1024*1024*1024);
    fprintf(stderr,"File size: %.2f GB\n",mapsize1);
  } else if (mapsize > 1024*1024) {
    mapsize1 = ((double)mapsize)/(1024*1024);
    fprintf(stderr,"File size: %.2f MB\n",mapsize1);
  } else if (mapsize > 1024) {
    mapsize1 = ((double)mapsize)/(1024);
    fprintf(stderr,"File size: %.2f kB\n",mapsize1);
  } else {
    mapsize1 = ((double)mapsize);
    fprintf(stderr,"File size: %.2f bytes\n",mapsize1);
  }
  cps = ((double)mapsize)/tsec;
  if (cps > 1024*1024*1024) {
    cps /= 1024*1024*1024;
    fprintf(stderr,"Scan rate: %.2f GB/s, %.2f Gb/s\n",cps,cps*8);
  } else if (cps > 1024*1024) {
    cps /= 1024*1024;
    fprintf(stderr,"Scan rate: %.2f MB/s, %.2f Mb/s\n",cps,cps*8);
  } else if (cps > 1024) {
    cps /= 1024;
    fprintf(stderr,"Scan rate: %.2f kB/s, %.2f kb/s\n",cps,cps*8);
  } else {
    fprintf(stderr,"Scan rate: %.2f B/s, %.2f b/s\n",cps,cps*8);
  }

}

