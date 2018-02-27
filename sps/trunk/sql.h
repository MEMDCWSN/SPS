/**
    @file sql.h

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


#ifndef SQL_H
#define SQL_H

#include <map>
#include <set>
#include <list>


struct SQLDerefTagType {};


/**
    Distinct rows.

    Selects distinct rows out of a sequence.
*/

template <typename T, typename U, U T::* P>
  class SQLDistinct1
  {
    typedef std::map<U, T *> map1;
    typedef typename map1::iterator iterator1;

  public:
    typedef std::list<T *> list;


    /**
        Constructor.

        @param  begin   Begining input iterator
        @param  end     Ending input iterator

        @note Complexity O(n)
    */
    //@{
    SQLDistinct1(T * begin, T * end)
    {
      for (T * i = begin; i != end; i ++)
        index1_.insert(std::pair<U, T *>(i->*P, &*i));
    }

    template <typename InputIterator>
      SQLDistinct1(InputIterator begin, InputIterator end)
      {
        for (InputIterator i = begin; i != end; i ++)
          index1_.insert(std::pair<U, T *>((*i).*P, &*i));
      }
    //@}

    /**
        Constructor.

        @param  begin   Begining input iterator
        @param  end     Ending input iterator
        @param  tag     Dereference iterator tag

        @note Complexity O(n)
    */

    template <typename InputIterator>
      SQLDistinct1(InputIterator begin, InputIterator end, SQLDerefTagType tag)
      {
        for (InputIterator i = begin; i != end; i ++)
          index1_.insert(std::pair<U, T *>((*i)->*P, &**i));
      }

    /**
        SQL select.

        Same as: @code SELECT DISTINCT *;

        @return         Selected records

        @note Complexity O(log(n))
    */

    list operator () ()
    {
      list res;

      for (iterator1 i = index1_.begin(); i != index1_.end(); i ++)
        res.push_back(i->second);

      return res;
    }

  protected:
    map1 index1_;
  };


/**
    Single index.

    Indexes a table using 1 identifier.
*/

template <typename T, typename U, U T::* P>
  class SQLSelect1
  {
    typedef std::multimap<U, T *> map1;
    typedef typename map1::iterator iterator1;

  public:
    typedef std::list<T *> list;


    /**
        Constructor.

        @param  begin   Begining input iterator
        @param  end     Ending input iterator

        @note Complexity O(n)
    */
    //@{
    SQLSelect1(T * begin, T * end)
    {
      for (T * i = begin; i != end; i ++)
        index1_.insert(std::pair<U, T *>(i->*P, &*i));
    }

    template <typename InputIterator>
      SQLSelect1(InputIterator begin, InputIterator end)
      {
        for (InputIterator i = begin; i != end; i ++)
          index1_.insert(std::pair<U, T *>((*i).*P, &*i));
      }
    //@}

    /**
        Constructor.

        @param  begin   Begining input iterator
        @param  end     Ending input iterator
        @param  tag     Dereference iterator tag

        @note Complexity O(n)
    */

    template <typename InputIterator>
      SQLSelect1(InputIterator begin, InputIterator end, SQLDerefTagType tag)
      {
        for (InputIterator i = begin; i != end; i ++)
          index1_.insert(std::pair<U, T *>((*i)->*P, &**i));
      }

    /**
        SQL select.

        Same as: @code SELECT * WHERE P = u;

        @param    u     First operand
        @return         Selected records

        @note Complexity O(log(n))
    */

    list operator () (const U & u)
    {
      list res;

      std::pair<iterator1, iterator1> r1 = index1_.equal_range(u);

      for (iterator1 i = r1.first; i != r1.second; i ++)
        res.push_back(i->second);

      return res;
    }

  protected:
    map1 index1_;
  };


/**
    Double indexes.

    Indexes a table with 2 identifiers.
*/

template <typename T, typename U, U T::* P, typename V, V T::* Q>
  class SQLSelect2 : public SQLSelect1<T, U, P>
  {
    typedef SQLSelect1<T, U, P> base;

    typedef std::multimap<U, T *> map1;
    typedef std::multimap<V, T *> map2;
    typedef typename map1::iterator iterator1;
    typedef typename map2::iterator iterator2;

  public:
    typedef std::list<T *> list;


    /**
        Constructor.

        @param  begin   Begining input iterator
        @param  end     Ending input iterator

        @note Complexity O(n)
    */
    //@{
    SQLSelect2(T * begin, T * end) : base(begin, end)
    {
      for (T * i = begin; i != end; i ++)
        index2_.insert(std::pair<V, T *>(i->*Q, &*i));
    }

    template <typename InputIterator>
      SQLSelect2(InputIterator begin, InputIterator end) : base(begin, end)
      {
        for (InputIterator i = begin; i != end; i ++)
          index2_.insert(std::pair<V, T *>((*i).*Q, &*i));
      }
    //@}

    /**
        Constructor.

        @param  begin   Begining input iterator
        @param  end     Ending input iterator
        @param  tag     Dereference iterator tag

        @note Complexity O(n)
    */

    template <typename InputIterator>
      SQLSelect2(InputIterator begin, InputIterator end, SQLDerefTagType tag) : base(begin, end, tag)
      {
        for (InputIterator i = begin; i != end; i ++)
          index2_.insert(std::pair<V, T *>((*i)->*Q, &**i));
      }

    /**
        SQL select.

        Same as: @code SELECT * WHERE P = u AND Q = v;

        @param    u     First operand
        @param    v     Second operand
        @return         Selected records

        @note Complexity O(log(n))
    */
    //@{
    using base::operator ();

    list operator () (const U & u, const V & v)
    {
      list res;

      std::pair<iterator1, iterator1> r1 = index1_.equal_range(u);
      std::pair<iterator2, iterator2> r2 = index2_.equal_range(v);

      std::set<T *> set1;
      std::set<T *> set2;

      for (iterator1 i = r1.first; i != r1.second; i ++)
        set1.insert(i->second);

      for (iterator2 i = r2.first; i != r2.second; i ++)
        set2.insert(i->second);

      std::set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(), std::back_inserter(res));

      return res;
    }
    //@}

  protected:
    using base::index1_;

    map2 index2_;
  };


/**
    Triple indexes.

    Indexes a table with 3 identifiers.
*/

template <typename T, typename U, U T::* P, typename V, V T::* Q, typename W, W T::* R>
  class SQLSelect3 : public SQLSelect2<T, U, P, V, Q>
  {
    typedef SQLSelect2<T, U, P, V, Q> base;

    typedef std::multimap<U, T *> map1;
    typedef std::multimap<V, T *> map2;
    typedef std::multimap<W, T *> map3;
    typedef typename map1::iterator iterator1;
    typedef typename map2::iterator iterator2;
    typedef typename map3::iterator iterator3;

  public:
    typedef std::list<T *> list;


    /**
        Constructor.

        @param  begin   Begining input iterator
        @param  end     Ending input iterator

        @note Complexity O(n)
    */
    //@{
    SQLSelect3(T * begin, T * end) : base(begin, end)
    {
      for (T * i = begin; i != end; i ++)
        index3_.insert(std::pair<W, T *>(i->*R, &*i));
    }

    template <typename InputIterator>
      SQLSelect3(InputIterator begin, InputIterator end) : base(begin, end)
      {
        for (InputIterator i = begin; i != end; i ++)
          index3_.insert(std::pair<W, T *>((*i).*R, &*i));
      }
    //@}

    /**
        Constructor.

        @param  begin   Begining input iterator
        @param  end     Ending input iterator
        @param  tag     Dereference iterator tag

        @note Complexity O(n)
    */

    template <typename InputIterator>
      SQLSelect3(InputIterator begin, InputIterator end, SQLDerefTagType tag) : base(begin, end, tag)
      {
        for (InputIterator i = begin; i != end; i ++)
          index3_.insert(std::pair<W, T *>((*i)->*R, &**i));
      }

    /**
        SQL select.

        Same as: @code SELECT * WHERE P = u AND Q = v AND R = w;

        @param    u     First operand
        @param    v     Second operand
        @param    w     Third operand
        @return         Selected records

        @note Complexity O(log(n))
    */
    //@{
    using base::operator ();

    list operator () (const U & u, const V & v, const W & w)
    {
      list res;

      std::pair<iterator1, iterator1> r1 = index1_.equal_range(u);
      std::pair<iterator2, iterator2> r2 = index2_.equal_range(v);
      std::pair<iterator3, iterator3> r3 = index3_.equal_range(w);

      std::set<T *> set1;
      std::set<T *> set2;
      std::set<T *> set3;
      std::set<T *> set4;

      for (iterator1 i = r1.first; i != r1.second; i ++)
        set1.insert(i->second);

      for (iterator2 i = r2.first; i != r2.second; i ++)
        set2.insert(i->second);

      for (iterator3 i = r3.first; i != r3.second; i ++)
        set3.insert(i->second);

      std::set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(), std::inserter(set4, set4.end()));
      std::set_intersection(set4.begin(), set4.end(), set3.begin(), set3.end(), std::back_inserter(res));

      return res;
    }
    //@}

  protected:
    using base::index1_;
    using base::index2_;

    map3 index3_;
  };


#endif
