
#ifndef _IBPEP_memory_debug_h_
#define _IBPEP_memory_debug_h_

#include <stdio.h>

#define QUOTE(a) #a

#define outputsizeof(os,type) os << "sizeof(" << QUOTE(type) << ") = " << sizeof(type) << endl

#ifdef DEBUG_MEMORY

#define NEW_DIAGNOSTIC(t,p,b,f,l) \
  fprintf(stderr,"0x%p = %s::new(%ld)[%s:%d]\n",p,t,b,f,l)

#define DELETE_DIAGNOSTIC(t,p,b,f,l) \
     fprintf(stderr,"%s::delete(0x%p,%ld)[%s:%d]\n",t,p,b,f,l)

#define ARRAYNEW_DIAGNOSTIC(t,p,b,f,l) \
  fprintf(stderr,"0x%p = %s::new[](%ld)[%s:%d]\n",p,t,b,f,l)

#define ARRAYDELETE_DIAGNOSTIC(t,p,b,f,l) \
     fprintf(stderr,"%s::delete[](0x%p,%ld)[%s:%d]\n",t,p,b,f,l)

#define MEMORY_DEBUG(type) \
\
  void* operator new(size_t bytes)\
{ \
  void *p = ::operator new(bytes); \
  NEW_DIAGNOSTIC(QUOTE(type),p,bytes,__FILE__,__LINE__);\
  return p;\
 }\
  void* operator new [] (size_t bytes)\
{ \
  void *p = ::operator new[](bytes); \
  ARRAYNEW_DIAGNOSTIC(QUOTE(type),p,bytes,__FILE__,__LINE__);\
  return p;\
 }\
\
  void* operator new (unsigned int bytes, void * const & p0)\
{ \
  void *p = ::operator new (bytes,p0); \
  NEW_DIAGNOSTIC(QUOTE(type),p,bytes,__FILE__,__LINE__);\
  return p;\
 }\
\
  void  operator delete(void* p, size_t bytes)\
{ ::operator delete(p);\
  DELETE_DIAGNOSTIC(QUOTE(type),p,bytes,__FILE__,__LINE__);\
 }\
\
  void  operator delete[](void* p, size_t bytes)\
{ ::operator delete[](p);\
  ARRAYDELETE_DIAGNOSTIC(QUOTE(type),p,bytes,__FILE__,__LINE__);\
 }\


#else 

#define MEMORY_DEBUG(type)

#endif

#endif
