/**
    @file dbg_print.h

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


#ifndef DBG_PRINT_H
#define DBG_PRINT_H

#include <vector>
#include <iostream>


template <typename T>
  inline std::ostream & DBG_PRINT(std::ostream & out, const T &a, int s = 0);

template <typename T,typename U>
  inline std::ostream & DBG_PRINT(std::ostream & out, const std::pair<T, U> &a, int s = 0);

template <template <typename> class C, typename T>
  inline std::ostream & DBG_PRINT(std::ostream & out, const C<T> &a, int s = 0);


template <typename T>
  inline std::ostream & DBG_PRINT(std::ostream & out, const T &a, int s)
  {
    out << a;

    return out;
  }

template <typename T,typename U>
  inline std::ostream & DBG_PRINT(std::ostream & out, const std::pair<T, U> &a, int s)
  {
    out << std::endl;
    for (int i = 0; i < s; i ++)
      out << ' ';

    out << "{";
    DBG_PRINT(out, a.first, s + 1);

    out << " ";

    DBG_PRINT(out, a.second, s + 1);

    out << "}";

    return out;
  }

template <template <typename> class C, typename T>
  inline std::ostream & DBG_PRINT(std::ostream & out, const C<T> &a, int s)
  {
    out << std::endl;
    for (int i = 0; i < s; i ++)
      out << ' ';

    out << "{";
    for (typename C<T>::const_iterator i = a.begin(); i != a.end(); i ++)
      if (i == a.begin())
        DBG_PRINT(out, * i, s + 1);
      else
        DBG_PRINT(out << " ", * i, s + 1);

    out << "}";

    return out;
  }


/**
    Log utils.
*/

template <typename T>
  inline T & note(T & sout)
  {
    return static_cast<T &>(sout << std::hex << std::uppercase << getpid() << std::dec << ": ");
  }

template <typename T>
  inline T & note(T & sout, const char * pf, int nl)
  {
    return static_cast<T &>(note(sout) << pf << ": " << nl << ": ");
  }


#endif // DBG_PRINT_H
