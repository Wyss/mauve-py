#ifndef _timing_h_
#define _timing_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// an opaque timer type.
typedef struct dmtimer_s dmtimer_t;

// starts the timer
dmtimer_t * StartTimer();

// reads the timer (msec)
unsigned int ReadTimer( dmtimer_t * t );

// stops the timer.
void StopTimer( dmtimer_t * t );

// initialize the timing code.
void InitTime();


#endif /* _timing_h_ */
