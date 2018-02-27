/**
    @file ioscopy.h

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


#ifndef IOS_H
#define IOS_H

#include <sstream>
#include <iostream>

#include "cstream.h"


namespace sps
{

/**
    Iostream with copy constructor.
*/

template <typename T>
  struct ios;


/**
    Iostream with copy constructor.

    Handles stringstream.
*/

template <typename T>
  struct ios : public std::stringstream
  {
    using std::stringstream::copyfmt;
    using std::stringstream::clear;

    /**
        Constructor.

        @param  in    Iostream

        @note   The buffer is copied.
    */

    ios(const std::stringstream & sin) : std::stringstream(sin.str())
    {
      copyfmt(sin);
      clear(sin.rdstate());
    }

    ios(const ios & sin) : std::stringstream(sin.str())
    {
      copyfmt(sin);
      clear(sin.rdstate());
    }

    ios()
    {
    }

    ios & operator = (const ios & sin)
    {
      copyfmt(sin);
      clear(sin.rdstate());
      str(sin.str());

      return * this;
    }
  };


/**
    Iostream with copy constructor.

    Handles istream.
*/

template <>
  struct ios<std::istream> : public std::istream
  {
    using std::istream::copyfmt;
    using std::istream::clear;
    using std::ios::rdbuf;

    /**
        Constructor.

        @param  in    istream

        @note   Automatic rewind involved.
    */

    ios(const std::istream & sin) : std::istream(sin.rdbuf())
    {
      copyfmt(sin);
      clear(sin.rdstate());

      if (rdbuf())
        rdbuf()->pubseekpos(0);
    }

    ios(const ios & sin) : std::istream(sin.rdbuf())
    {
      copyfmt(sin);
      clear(sin.rdstate());

      if (rdbuf())
        rdbuf()->pubseekpos(0);
    }

    ios() : std::istream(0)
    {
    }

    ios & operator = (const ios & sin)
    {
      copyfmt(sin);
      clear(sin.rdstate());
      rdbuf(sin.rdbuf());

      if (rdbuf())
        rdbuf()->pubseekpos(0);

      return * this;
    }
  };


/**
    Iostream with copy constructor.

    Handles ostream.
*/

template <>
  struct ios<std::ostream> : public std::ostream
  {
    using std::ostream::copyfmt;
    using std::ostream::clear;
    using std::ios::rdbuf;

    /**
        Constructor.

        @param  in    ostream

        @note   Automatic rewind involved.
    */

    ios(const std::ostream & sout) : std::ostream(sout.rdbuf())
    {
      copyfmt(sout);
      clear(sout.rdstate());

      if (rdbuf())
        rdbuf()->pubseekpos(0);
    }

    ios(const ios & sin) : std::ostream(sin.rdbuf())
    {
      copyfmt(sin);
      clear(sin.rdstate());

      if (rdbuf())
        rdbuf()->pubseekpos(0);
    }

    ios() : std::ostream(0)
    {
    }

    ios & operator = (const ios & sin)
    {
      copyfmt(sin);
      clear(sin.rdstate());
      rdbuf(sin.rdbuf());

      if (rdbuf())
        rdbuf()->pubseekpos(0);

      return * this;
    }
  };

} // namespace sps


#endif
