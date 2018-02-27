/*
 * SvmScaleParameterList.h
 *
 *  Created on: Jan 1, 2011
 *      Author: jsnedecor
 */

#ifndef SVMMODEL_H_
#define SVMMODEL_H_

//system includes
#include <string>
#include <vector>

//external includes
#include "utils.h"
#include "DelimitedTextReader.h"
#include "Logger.h"
#include "SvmScaleParameterList.h"

namespace specnets
{
  class SvmModel
  {
  public:
    //! \name CONSTRUCTORS
    //@{
    SvmModel();

    SvmModel(SvmScaleParameterList * modelKeyRanges);

    //! \name DESTRUCTOR
    //@{
    virtual ~SvmModel();

    /*!
     * \brief sets SvmScaleParameter list
     *
     */
    void setKeyRanges(SvmScaleParameterList * modelKeyRanges);

    bool loadSvmModel(const char * modelFile);

    /*!
     * \brief returns value of rho
     */
    double getRho() const;

    /*!
     * \brief return value of gamma
     *
     */
    double getGamma() const;

    /*!
     * \brief returns number of model vectors
     */
    unsigned int size() const;

    /*!
     * \brief return scaling factor for vector at index i
     */
    double getScalingFactor(unsigned int i) const;

    /*!
     * \brief return value of key keyIndex for vector i
     */
    bool getKeyScalingFactor(unsigned int i,
                             unsigned int keyIndex,
                             double &returnDouble);
    /*!
     * \brief rescale vector to match expected range
     */
    void scaleVector(vector<float> &input, vector<float> &output);

    /*!
     * \brief get Svm score for model
     */
    double getSvmScore(vector<float> &inputVector);

    /*!
     * \brief get key index by name
     */
    bool getKeyIndexByName(string &name, unsigned int &outputIndex);

    /*!
     * \brief get key index size
     */
    unsigned int getKeyIndexSize(void);
  protected:

    //rho value
    double m_rho;

    //gamma value
    double m_gamma;

    //Scaling factor for each vector
    vector<double> m_scalingFactors;

    //values from delimited text file
    vector<map<unsigned int, double> > m_model;

    //Scaling ranges for model
    SvmScaleParameterList * m_modelKeyRanges;

  };
}
#endif /* SVMMODEL_H_ */
