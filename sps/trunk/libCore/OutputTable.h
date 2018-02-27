/*
 * Table.h
 *
 *  Created on: Mar 3, 2011
 *      Author: aguthals
 */

#ifndef OUTPUTTABLE_H_
#define OUTPUTTABLE_H_

#include <cstring>
#include <string>
#include <vector>

#include "aminoacid.h"
#include "Logger.h"
#include "utils.h"

namespace specnets
{

  /*
   * 2D Table that is printed to file
   */
  class OutputTable
  {
  public:

    // contents of table. cell.first = text, cell.second = whether cell is a header
    vector<vector<pair<string, bool> > > values;

    OutputTable(void);

    void setValue(int row, int col, string value, bool header = false);

    void loadHistogram(vector<TwoValues<float> > &bins);

    /**
     * Override this method to prepare table
     */
    void prepareTable(void);

    /**
     * Prints contents of values to file
     * @param filename
     * @return true of file was written successfully, false if not
     */
    bool printToCSV(const char* filename, string delim = ",");

    /**
     * Loads values from a file
     * @param filename
     * @return true of file was loaded successfully, false if not
     */
    bool loadFromCSV(const char* filename, string delim = ",");

    /**
     * TODO: printToHTML
     */

  };
}

#endif /* OUTPUTTABLE_H_ */
