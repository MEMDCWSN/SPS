/*
 * FalseDiscoveryRatePeptide.h
 *
 *  Created on: Feb 14, 2011
 *      Author: jsnedecor
 *
 *      THIS CLASS IS DEPRECATED! PLEASE USE PeptideSpectrumMatchSet class INSTEAD!
 */

#ifndef FALSEDISCOVERYRATEPEPTIDES_H_
#define FALSEDISCOVERYRATEPEPTIDES_H_

#include "inspect_parse.h"

#include <string>
#include <vector>
#include <algorithm>

#include "Logger.h"

using namespace std;

namespace specnets
{
  struct FalseDiscoveryRatePeptide
  {
    InspectResultsLine * m_peptideResult;
    unsigned int m_scanNum;
    double m_score;
    bool m_isDecoy;
  };

  class FalseDiscoveryRatePeptides
  {
  public:
    //! \name CONSTRUCTORS
    //@{
    FalseDiscoveryRatePeptides();

    //! \name DESTRUCTOR
    //@{
    virtual ~FalseDiscoveryRatePeptides();

    void addFdrPeptide(FalseDiscoveryRatePeptide  &input);

    bool calculatePValues();

    bool getPValueByScan(unsigned int scanNum, double &score);

    bool getPValueByIndex(unsigned int indexNum, double &score);

    bool getMQScoreByIndex(unsigned int indexNum, double &score);

    bool getScanNumByIndex(unsigned int indexNum, unsigned int &scanNum);

    InspectResultsLine * getResultByIndex(unsigned int indexNum);

    unsigned int size();

  protected:
    vector<FalseDiscoveryRatePeptide> m_inputPeptides;
    vector<double> m_pValues;
    bool m_calculatedPValue;
  };
}
#endif /* FALSEDISCOVERYRATEPEPTIDE_H_ */
