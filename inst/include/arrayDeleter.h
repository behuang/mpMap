#ifndef _ARRAY_DELETER_H
#define _ARRAY_DELETER_H

//helper struct for a reference-counter array
template<typename T> struct arrayDeleter
{
   void operator()(T* p)
   {
      delete [] p;
   }
};

#endif