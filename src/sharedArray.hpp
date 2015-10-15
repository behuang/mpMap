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
#ifdef __linux__
#include <sys/mman.h>
template<typename T> struct mmapDeleter
{
public:
	mmapDeleter(std::size_t length)
		:length(length)
	{}
	void operator()(T* p)
	{
		munmap(p, sizeof(T)*length);
	}
private:
	std::size_t length;
};
template<typename T> class sharedMmapArray : protected std::shared_ptr<T>
{
public:
	sharedMmapArray(std::size_t length)
		: std::shared_ptr<T>(), length(length)
	{
		void* result = mmap(NULL, sizeof(T)*length, PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);
		if(result == MAP_FAILED)
		{
			throw std::runtime_error("Unable to allocate sufficient memory");
		}
		*static_cast<std::shared_ptr<T>* >(this) = std::shared_ptr<T>((T*)result, mmapDeleter<T>(length));
	}
	sharedMmapArray()
		: std::shared_ptr<T>()
	{}
	T* get()
	{
		return std::shared_ptr<T>::get();
	}
private:
	std::size_t length;
};
#endif
#endif
