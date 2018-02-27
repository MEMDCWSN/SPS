/**
 * Helper functions related to Inspect result parsing
 *
 *  Created on: Aug 18, 2010
 *      Author: jsnedecor
 *
 *      NOTE: THIS CLASS IS DEPRECATED! Please use PeptideSpectrumMatch class for any future work.
 */

#ifndef INSPECT_PARSE_H_
#define INSPECT_PARSE_H_

#include "utils.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>
#include <cstring>

namespace specnets
{
  class InspectResultsLine
  {
  public:
    bool parseFromFields(const vector<string>& fields,
                         const vector<string>& field_names);

    string SpectrumFile;
    int scan;
    string Annotation;
    string Protein;
    int Charge;
    float MQScore;
    float Score;
    int Length;
    float TotalPRMScore;
    float MedianPRMScore;
    float FractionY;
    float FractionB;
    float Intensity;
    int NTT;
    float p_value;
    float F_Score;
    float DeltaScore;
    float DeltaScoreOther;
    int RecordNumber;
    int DBFilePos;
    int SpecFilePos;
  };

  /**Inspect parsing and translation between SpecNets and Inspect formats
   */
  class InspectResultsSet
  {
  public:
    vector<InspectResultsLine> results; //scan number, InspectResultsLine

    /** Inspect annotation parsing
     * Changes Inspect annotations to SpecNets (should add converse function)
     * @param inspect = inspect string
     * @param specnets = output string for specnets
     */
    void inspectToSpecNets(string &inspect, string &specnets);

    /** Inspect result file parsing
     *
     */
    bool loadInspectResultsFile(const char *inspect_results_file);

    bool getResultByScan(int scan, vector<InspectResultsLine>& output_results);
  };
}
#endif /* INSPECT_PARSE_H_ */
