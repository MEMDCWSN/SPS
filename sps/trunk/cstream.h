/**
    @file cstream.h

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


#ifndef CSTREAM_H
#define CSTREAM_H

#include <errno.h>
#include <stdio.h>
#include <unistd.h>

#include <iostream>
#include <ext/stdio_filebuf.h>


/**
    C++ stdio wrapper

    Glues streams to system IO.

    @note
    Uses GCC extension.
*/

class cstream_base
{
public:
  cstream_base(std::ios_base::openmode nmode): np(0), pf(::tmpfile()), buf(pf, nmode) {}
  cstream_base(FILE * pf, std::ios_base::openmode nmode): np(0), pf(pf), buf(pf, nmode) {}

  ~cstream_base()
  {
    if (np)
      close(np);
  }

  int dup2(int ni)
  {
    np = ni;

    ::rewind(pf);

    if (::dup2(buf.fd(), np) == -1)
      return errno;

    if (::fclose(pf))
      return errno;

    return 0;
  }

  FILE * file_pointer() const
  {
    return pf;
  }

  int file_descriptor() const
  {
    return fileno(pf);
  }

protected:
  int np;
  FILE * pf;
  __gnu_cxx::stdio_filebuf<char> buf;
};

class cstream : public cstream_base, public std::iostream
{
  cstream(int);

public:
  cstream() : cstream_base(std::ios::in | std::ios::out), std::iostream(& buf) {}
  cstream(FILE * pf) : cstream_base(pf, std::ios::in | std::ios::out), std::iostream(& buf) {}
  cstream(const cstream & cs) : cstream_base(::fdopen(::dup(buf.fd()), "w+b"), std::ios::in | std::ios::out), std::iostream(& buf) {}

  int dup2(int ni)
  {
    flush();

    return cstream_base::dup2(ni);
  }

  int dup2(FILE * np)
  {
    return dup2(fileno(np));
  }
};


#endif // CSTREAM_H
