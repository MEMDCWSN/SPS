/**
    @file ion.h

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


#ifndef ION_H
#define ION_H

#include <map>
#include <string>
#include <iostream>


/**
    Ion representation.
*/

class Ion
{
public:
  enum Type {a, b, c, x, y, z, phos};

  /**
      Constructor.

      @param  t   Type of ion
      @param  m   Mass
      @param  a   Amino acid
      @param  p   Likelihood
  */

  Ion(bool c, Type t, double m, const std::string & a, double p) : charge(c), type(t), mass(m), acid(a), probability(p) {}

  /**
      Comparison operator.

      @param  t   Operand
      @return     True if lesser

      @note For use with sorted containers.
  */

  bool operator < (const Ion &t) const
  {
    return mass < t.mass;
  }

  bool charge;
  Type type;
  double mass;
  std::string acid;
  double probability;

  static const Ion ion[15];
  static std::map<Type, std::string> cname;
  static std::map<std::string, Type> ctype;
};


#endif
