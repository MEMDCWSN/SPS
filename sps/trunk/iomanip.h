#ifndef IOMANIP_H
#define IOMANIP_H

#include <iostream>


namespace sps
{

inline std::ostream & clear(std::ostream & sout)
{
  sout.clear(sout.rdstate() & ~std::ios::failbit & ~std::ios::badbit);

  return sout;
}

inline std::istream & rewind(std::istream & sin)
{
  sin.clear();
  sin.seekg(0, std::ios::beg);

  return sin;
}

} // namespace sps


#endif
