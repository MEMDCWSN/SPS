/**
    @file uncertain.h

    @note
    Copyright 2006, The Regents of the University of California
    All Rights Reserved
    
    Permission to use, copy, modify and distribute any part of this 
    program for educational, research and non-profit purposes, without fee, 
    and without a written agreement is hereby granted, provided that the 
    above copyright notice, this paragraph and the following three paragraphs 
    appear in all copies.
    
    Those desiring to incorporate this work into commercial 
    products or use for commercial purposes should contact the Technology 
    Transfer & Intellectual Property Services, University of California, 
    San Diego, 9500 Gilman Drive, Mail Code 0910, La Jolla, CA 92093-0910, 
    Ph: (858) 534-5815, FAX: (858) 534-7345, E-MAIL:invent@ucsd.edu.
    
    IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY 
    FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, 
    INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN 
    IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY 
    OF SUCH DAMAGE.
    
    THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE UNIVERSITY 
    OF CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, 
    ENHANCEMENTS, OR MODIFICATIONS.  THE UNIVERSITY OF CALIFORNIA MAKES NO 
    REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND, EITHER IMPLIED OR 
    EXPRESS, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
    MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF 
    THE SOFTWARE WILL NOT INFRINGE ANY PATENT, TRADEMARK OR OTHER RIGHTS.
*/


#ifndef UNCERTAIN_H
#define UNCERTAIN_H

#include <cmath>
#include <map>
//TR1 includes. GCC 4.0 and above only!
#ifdef __GLIBCXX__
#  include <tr1/unordered_map>
#else
#  ifdef __IBMCPP__
#    define __IBMCPP_TR1__
#  endif
#  include <unordered_map>
#endif

#define hash_map tr1::unordered_map
#define data_type mapped_type

/**
    Uncertainty number.

    Measurement including relative error.
*/

template <typename T>
  class uncertain : private std::pair<T, T>
  {
    typedef std::pair<T, T> base;

  public:
    using base::first;
    using base::second;

    /**
        Constructor.

        @param  t       Value.
        @param  e       Relative error.

        @note Half of the relative error is stored because it is later doubled on comparisons.
    */

    explicit uncertain(const T &t, const T &e) : base(t - e / 2, t + e / 2) {}

    T value() const
    {
      return (first + second) / 2;
    }

    T e() const
    {
      return std::abs(second - first) / 2;
    }

    friend bool operator < (const uncertain &nu, const uncertain &nv)
    {
      return nu.second < nv.first;
    }

    friend bool operator == (const uncertain &nu, const uncertain &nv)
    {
      return nu.first <= nv.second && nv.first <= nu.second;
    }

    friend bool operator != (const uncertain &nu, const uncertain &nv)
    {
      return ! operator == (nu, nv);
    }
  };


namespace std { namespace tr1
{
template <typename T>
  struct hash< uncertain<T> >
  {
    size_t operator() (const uncertain<T> & t) const
    {
      // size_t(std::floor(t.value() / (t.e() * 2) + 1) / 2))
      return size_t(floor((t.value() / abs(t.second - t.first) + 1) / 2));
    }
  };
}}


namespace std
{

template <typename T, typename U>
  struct map< uncertain<T>, U>;

template <typename T, typename U>
  struct multimap< uncertain<T>, U>;

}


#endif // UNCERTAIN_H
