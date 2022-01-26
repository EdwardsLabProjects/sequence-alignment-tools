#ifndef NATHANS_TIME_HACKS_H
#define NATHANS_TIME_HACKS_H

// Nathan's time counter
#include <time.h>
static time_t tictoctime_;
#define TIC {tictoctime_ = time(NULL);}
#define TOC(F) { time_t t; time(&t); t-=tictoctime_; fprintf(F,"Elapsed time: %d:%02d:%02d",t/3600,(t%3600)/60,t%60);} 

#endif
