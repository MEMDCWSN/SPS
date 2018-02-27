#ifndef VECTOR_H
#define VECTOR_H

#include <string>
#include <vector>
#include <sstream>
#include <stdexcept>


namespace sps
{

template < typename T, typename = std::allocator<T> >
  struct vector : std::vector<T>
  {
      typedef typename std::vector<T>::size_type size_type;

      std::string sfilename;

      explicit vector() {}

      explicit vector(size_t ns, const T & cv = T()) : std::vector<T>(ns, cv) {}

      vector(const vector & cv) : std::vector<T>(cv) {}

      template<typename I>
        vector(I ifirst, I ilast) : std::vector<T>(ifirst, ilast) {}

      T & operator [] (size_t ns)
      {
        range_check(ns);

        return std::vector<T>::operator [] (ns);
      }

      const T & operator [] (size_t ns) const
      {
        range_check(ns);

        return std::vector<T>::operator [] (ns);
      }

  private:
      void range_check(size_type ns) const
      {
        if (std::vector<T>::size() < ns)
        {
          std::ostringstream smessage;

          smessage << sfilename << ": " << std::vector<T>::size() << " < " << ns;

          throw std::out_of_range(smessage.str().c_str());
        }
      }
  };

template <>
  struct vector<bool> : vector<char>
  {
      using vector<char>::size_type;

      explicit vector() {}

      explicit vector(size_t ns, const char & cv = char()) : vector<char>(ns, cv) {}

      vector(const vector & cv) : vector<char>(cv) {}

      template<typename I>
        vector(I ifirst, I ilast) : vector<char>(ifirst, ilast) {}
  };

} // namespace sps

#endif
