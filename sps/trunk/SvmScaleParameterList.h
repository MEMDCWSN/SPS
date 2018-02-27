/*
 * SvmScaleParameterList.h
 *
 *  Created on: Jan 1, 2011
 *      Author: jsnedecor
 */

#ifndef SVMSCALEPARAMETERLIST_H_
#define SVMSCALEPARAMETERLIST_H_

//system includes
#include <string>
#include <vector>

//external includes
#include "DelimitedTextReader.h"
#include "Logger.h"

namespace specnets
{
  class SvmScaleParameterList
  {
  public:
    /*!
     * \brief Load in Svm parameters by name
     * @param filename - filename
     */

    //! \name CONSTRUCTORS
     //@{
    SvmScaleParameterList();

     //! \name DESTRUCTOR
     //@{
     virtual ~SvmScaleParameterList();

    bool loadSvmKeyRange(const char * keyFile);

    /*!
     * \brief returns name of svm key
     * @param index - index within file (0 to numLines -1)
     */
    string getKeyName(unsigned int lineIndex) const;

    /*!
     * \brief returns index of svm key
     * @param index - index within file (0 to numLines -1)
     */
    unsigned int getKeyIndex(unsigned int lineIndex) const;

    /*!
     * \brief returns start range of svm key
     * @param index - index within file (0 to numLines -1)
     */
    double getStartRange(unsigned int lineIndex) const;

    /*!
     * \brief returns start range of svm key
     * @param index - index within file (0 to numLines -1)
     */
    double getEndRange(unsigned int lineIndex) const;

    /*! \brief Returns size of m_lines vector

     */
    unsigned int size() const;

    /*! \brief Returns key index based on name
     *
     */

    bool getKeyIndexByName(string &name, unsigned int &outputIndex);


  protected:
    //indices of header
    vector<int> m_requiredHeaderIndex;

    //expected header
    vector<string> m_requiredHeader;

    //header to read in from file
    map<string, unsigned int> m_header;

    //values from delimited text file
    vector<vector<string> > m_lines;

    map<string, unsigned int> m_nameToKey; //mapping of name to key index.
  };
}
#endif /* SVMSCALEPARAMETERLIST_H_ */
