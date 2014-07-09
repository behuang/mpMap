#ifndef _SHARED_ARRAY_H
#define _SHARED_ARRAY_H

#include <memory>
//helper struct for a reference-counter array
template<typename T> struct arrayDeleter
{
   void operator()(T* p)
   {
      delete [] p;
   }
};
template<typename T> class sharedArray : public std::shared_ptr<T>
{
public:
	sharedArray(T* t)
		: std::shared_ptr<T>(t, arrayDeleter<T>())
	{}
	sharedArray()
		: std::shared_ptr<T>()
	{}
	T& operator[](int i)
	{
		return this->get()[i];
	}
};

#endif
