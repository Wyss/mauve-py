#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef WIN32
#include <sys/time.h>
#include <unistd.h>
//#include <malloc.h>
#endif

#include <stdlib.h>
#include <stdio.h>


#include "libMems/dmSML/util.h"
#include "libMems/dmSML/timing.h"


struct dmtimer_s {
#ifdef WIN32    
    unsigned int last;
#else
    struct timeval tv;
#endif
};



typedef int Int;
typedef unsigned int UInt;
typedef double Float64;

#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <mmsystem.h>
// keep this many significant bits from the PerformanceCounter values.
#define NUM_FREQ_BITS   (14)
static Int ShiftAmt;
static Int TicksPerSecond;
static Int LastReadValue;
static Int BaseTime;
#endif /* WIN32 */


dmtimer_t * StartTimer() {
#ifdef WIN32
    dmtimer_t * t = malloc( sizeof( *t ) );
    t->last = timeGetTime();
    return( t );
#else
    dmtimer_t * t = malloc( sizeof( *t ) );
    gettimeofday( &(t->tv), NULL );
    return( t );
#endif /* WIN32 */
}



unsigned int ReadTimer( dmtimer_t * t ) {
#ifdef WIN32
    /*
    Int ticks;
    LARGE_INTEGER pcnow;
    Float64 seconds;
    QueryPerformanceCounter( &pcnow );
    Shift64( ShiftAmt, (int*)&pcnow.HighPart, (int*)&pcnow.LowPart );
    ticks = pcnow.LowPart;
    LastReadValue = ticks;
    if( ticks < BaseTime ) {
        // handle wraparound.
        ticks += ((1 << NUM_FREQ_BITS)) - BaseTime;
    } else {
        ticks -= BaseTime;
    }
    seconds = (Float64)ticks / (Float64)TicksPerSecond;
    return( (int)(seconds * 10000 + 0.5) );
    */
    unsigned int cur = timeGetTime();
    return( cur - t->last );
#else
    struct timeval current;
    struct timezone dummy;
    unsigned int begintime, endtime;
    gettimeofday( &current, &dummy );
    begintime = 1000 * t->tv.tv_sec + (t->tv.tv_usec/1000);
    endtime = 1000 * current.tv_sec + (current.tv_usec/1000);
    return( endtime - begintime );
#endif
}



void StopTimer( dmtimer_t * t ) {
    free( t );
}



#ifdef WIN32
static void InitTimeWIN32() {

    timeBeginPeriod( 1 );

    /*
    LARGE_INTEGER pcfreq;
    UInt pchi, pclow, hihibit, lowhibit, highbit;
    UInt i;
    ShiftAmt = 0;
    QueryPerformanceFrequency( &pcfreq );
    pchi = pcfreq.HighPart;
    pclow = pcfreq.LowPart;
    // we want to look at the most significant 14 bits of the counter,
    // so we get about 1/10000th second accuracy
    // (between 8192ths - 16383ths second accuracy to be exact).
    // find the highest bit set in the high part.
    for( i = sizeof( pchi ) * 8; i ; i-- ) {
        if( pchi & 0x80000000 ) {
            break;
        }
        pchi = pchi << 1;
    }
    hihibit = i;
    // find the highest bit set in the low part.
    for( i = sizeof( pclow ) * 8; i ; i-- ) {
        if( pclow & 0x80000000 ) {
            break;
        }
        pclow = pclow << 1;
    }
    lowhibit = i;
    if( hihibit ) {
        highbit = hihibit + 32;
    } else {
        highbit = lowhibit;
    }
    pchi = pcfreq.HighPart;
    pclow = pcfreq.LowPart;
    if( highbit <= NUM_FREQ_BITS ) {
        ShiftAmt = 0;
    } else {
        ShiftAmt = highbit - NUM_FREQ_BITS;
    }
    Shift64( ShiftAmt, (int*)&pchi, (int*)&pclow );
    // now we have the most significant 14 bits of frequency.
    TicksPerSecond = pclow;
    // now actually read the counter, compute the ticks and store it away
    // so we have a base for the first call.
    QueryPerformanceCounter( &pcfreq );
    // this demonstrates the procedure for converting a LARGE_INTEGER
    // to ticks.
    Shift64( ShiftAmt, (int*)&pcfreq.HighPart, (int*)&pcfreq.LowPart );
    LastReadValue = pcfreq.LowPart;
    BaseTime = LastReadValue;
    */
}
#endif /* WIN32 */


void InitTime() {
#ifdef WIN32    
    InitTimeWIN32();
#endif
}
