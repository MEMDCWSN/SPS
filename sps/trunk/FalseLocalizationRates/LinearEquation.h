/*
 * LinearEquation.h
 *
 *  Created on: Oct 12, 2010
 *      Author: jsnedecor
 */

#ifndef LINEAREQUATION_H_
#define LINEAREQUATION_H_

//module includes

//external includes
#include "Logger.h"

//system includes
#include <iostream>
#include <iomanip>

#include <fstream>
#include <string>
#include <vector>

using namespace std;

namespace specnets
{
  class LinearEquation
  {
  public:
    /**
     * Loads in information from CPLEX LP formatted file.
     * @params filename - filename to input LP. Only one LP allowed per file
       @return True if data was loaded successfully, false otherwise.
     */
    //bool loadCPLEXLP(const char * filename);
    static bool readGlpsolOutput(const char * filename,
                          vector<float> &peakActivity,
                          vector<float> &peakError,
                          vector<float> &varientQuantities,
                          vector<string> &varientNames,
                          double &objectiveError);

    /**
     * Writes out LP into a CPLEX LP formatted file
     * @params filename - filename to output LP.
     */
    bool saveCPLEXLP(const char * filename);


    vector<vector<float> > m_lhsConstraints; //! Each row represents the lhs of a single equation. Each column is a single variable constraint.
    vector<float> m_rhsConstraints; //! Each value represents the right hand side constant for the equation.
  };
}


#endif /* LINEAREQUATION_H_ */
