/*
 * Table.cpp
 *
 *  Created on: Mar 3, 2011
 *      Author: aguthals
 */

#include "OutputTable.h"

namespace specnets
{
  OutputTable::OutputTable(void)
  {

  }

  void OutputTable::setValue(int row, int col, string value, bool header)
  {

    if (row < 0 || row > values.size())
    {
      DEBUG_MSG("Setting (" << row << ", " << col << ") = " << value);
      ERROR_MSG("Invalid assignment of row " << row << " to table with size " << values.size());
      return;
    }
    if (col < 0 || col > values[row].size())
    {
      DEBUG_MSG("Setting (" << row << ", " << col << ") = " << value);
      ERROR_MSG("Invalid assignment of col " << col << " to table row " << row << " with size " << values[row].size());
      return;
    }

    values[row][col].first = value;
    values[row][col].second = header;
  }

  void OutputTable::loadHistogram(vector<TwoValues<float> > &bins)
  {
    values.resize(bins.size() + 1);
    values[0].resize(2);
    setValue(0, 0, "Bin", true);
    setValue(0, 1, "Frequency", true);

    for (unsigned int i = 0; i < bins.size(); i++)
    {
      values[i + 1].resize(2);
      setValue(i + 1, 0, parseFloat(bins[i][0], 3), false);
      setValue(i + 1, 1, parseFloat(bins[i][1], 2), false);
    }
  }

  /**
   * Prints contents of values to file
   * @param filename
   * @return true of file was written successfully, false if not
   */
  bool OutputTable::printToCSV(const char* filename, string delim)
  {
    FILE* output = fopen(filename, "wb");

    if (output == NULL)
    {
      ERROR_MSG("Could not write to file " << filename);
      return false;
    }
    for (int i = 0; i < values.size(); i++)
    {
      if (values[i].size() == 0)
      {
        fprintf(output, "\n");
        continue;
      }
      else
      {
        fprintf(output, "%s", values[i][0].first.c_str());
      }
      for (int j = 1; j < values[i].size(); j++)
      {
        fprintf(output, "%s%s", delim.c_str(), values[i][j].first.c_str());
      }
      fprintf(output, "\n");
    }
    fclose(output);
    return true;
  }

  /**
   * Loads values from a file
   * @param filename
   * @return true of file was loaded successfully, false if not
   */
  bool OutputTable::loadFromCSV(const char* filename, string delim)
  {
    BufferedLineReader blr;
    if (blr.Load(filename) <= 0)
    {
      ERROR_MSG("Could not read file \'" << filename << "\'");
      return false;
    }
    values.resize(blr.size());

    int rowIdx = 0, colIdx = 0;
    for (int l = 0; l < blr.size(); l++)
    {
      if (strlen(blr.getline(l)) == 0)
      {
        continue;
      }
      vector<string> rowValues;
      splitText(blr.getline(l), rowValues, delim.c_str());

      values[rowIdx].resize(rowValues.size());
      for (int col = 0; col < rowValues.size(); col++)
      {
        values[rowIdx][col].first = rowValues[col];
        values[rowIdx][col].second = false;
      }
      rowIdx++;
    }
    values.resize(rowIdx);
    return true;
  }
}
