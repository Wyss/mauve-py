#ifndef __TWISTER_H__
#define __TWISTER_H__

#ifdef __cplusplus
extern "C" {
#endif

void SetTwisterSeed (unsigned long seed);
unsigned long CreateTwisterSeed(void);
double RandTwisterDouble (void);
unsigned RandTwisterUnsigned(void);

#ifdef __cplusplus
}
#endif

#endif //  __TWISTER_H__

