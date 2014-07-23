#ifndef __libMems_configuration_h__
#define __libMems_configuration_h__

#if defined(WIN32)||defined(WIN64)

// set the mems library name to include based on the configuration...

#if defined(WIN64)&&defined(NDEBUG)&&!defined(FASTDEBUG)&&defined(_OPENMP)
#pragma comment(lib, "mems64omp.lib")
#endif
#if defined(WIN64)&&defined(FASTDEBUG)&&defined(_OPENMP)
#pragma comment(lib, "mems64fdomp.lib")
#endif
#if defined(WIN32)&&!defined(WIN64)&&defined(NDEBUG)&&!defined(FASTDEBUG)&&defined(_OPENMP)
#pragma comment(lib, "memsomp.lib")
#endif
#if defined(WIN32)&&!defined(WIN64)&&defined(FASTDEBUG)&&defined(_OPENMP)
#pragma comment(lib, "memsfdomp.lib")
#endif
#if defined(WIN64)&&defined(NDEBUG)&&!defined(FASTDEBUG)&&!defined(_OPENMP)
#pragma comment(lib, "mems64.lib")
#endif
#if defined(WIN64)&&defined(FASTDEBUG)&&!defined(_OPENMP)
#pragma comment(lib, "mems64fd.lib")
#endif
#if defined(WIN32)&&!defined(WIN64)&&defined(NDEBUG)&&!defined(FASTDEBUG)&&!defined(_OPENMP)
#pragma comment(lib, "mems.lib")
#endif
#if defined(WIN32)&&!defined(WIN64)&&defined(FASTDEBUG)&&!defined(_OPENMP)
#pragma comment(lib, "memsfd.lib")
#endif


#endif

#endif // __libMems_configuration_h__

