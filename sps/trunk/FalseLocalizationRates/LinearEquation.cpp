#include "LinearEquation.h"

namespace specnets
{
  // -------------------------------------------------------------------------
  bool LinearEquation::readGlpsolOutput(const char * filename,
                                        vector<float> &peakActivity,
                                        vector<float> &peakError,
                                        vector<float> &varientQuantities,
                                        vector<string> &varientNames,
                                        double &objectiveError)
  {
    ifstream ipFilehandle(filename, ios::binary);

    if (!ipFilehandle.is_open() || !ipFilehandle.good())
    {
      ERROR_MSG("Unable to read file! " << filename);
      return false;
    }
    //initialize variables
    stringstream ss;
    string currLine;

    int rowNum;
    int columnNum;
    //CPLEX variables

    //read from file
    getline(ipFilehandle, currLine); //ignore first line
    ipFilehandle.ignore(12); //ignore "row" indicator;
    getline(ipFilehandle, currLine);
    ss << currLine;
    ss >> rowNum;

    //Get columns
    ipFilehandle.ignore(12); //ignore "column" indicator;
    ss.clear();
    ss.str(""); //clear stringbuffer
    currLine = "";
    getline(ipFilehandle, currLine);
    ss << currLine;
    ss >> columnNum;

    //ignore non-zeros
    ipFilehandle.ignore(256, '\n');
    //ignore status
    ipFilehandle.ignore(256, '\n');
    //ignore objective error
    ipFilehandle.ignore(20);
    ss.clear();
    ss.str(""); //clear stringbuffer
    currLine = "";
    getline(ipFilehandle, currLine);
    size_t objPos = currLine.find_first_of(" ");
    ss << currLine.substr(0, objPos);
    ss >> objectiveError;

    //ignore next three lines.
    ipFilehandle.ignore(256, '\n');
    ipFilehandle.ignore(256, '\n');
    ipFilehandle.ignore(256, '\n');

    //peak activity
    peakActivity.resize(rowNum);
    vector<int> tempPeakActivity(rowNum);

    float total = 0;
    for (int i = 0; i < rowNum; i++)
    {
      int tempPeak;
      //ignore names
      ipFilehandle.ignore(23);
      char * str = new char[256];
      //index
      ipFilehandle.read(str, 13);
      ss.clear();
      ss.str(""); //clear stringbuffer
      ss << str;
      ss >> tempPeak;
      tempPeakActivity[i] = tempPeak;
      total += tempPeak;
      ipFilehandle.ignore(256, '\n');
      delete [] str;
      str = NULL;
    }

    for (int i = 0; i < rowNum; i++)
    {
      peakActivity[i] = tempPeakActivity[i] / total;
    }
    //ignore columns
    ipFilehandle.ignore(256, '\n');
    ipFilehandle.ignore(256, '\n');
    ipFilehandle.ignore(256, '\n');

    vector<double> tempPeakError;
    vector<double> tempVarientQuantities;

    bool ticker = false;
    double tempError = 0;
    double tempQuantile = 0;
    double errorTotal = 0;
    double varientTotal = 0;

    for (int i = 0; i < columnNum; i++)
    {
      //read in peak error
      ipFilehandle.ignore(7);
      char * colname = new char[256];

      ipFilehandle.read(colname, 12); //Column name
      string col;
      col.assign(colname);

      char * str = new char[256];
      ipFilehandle.ignore(4); //next col

      ipFilehandle.read(str, 13);

      //clear string buffer
      ss.clear();
      ss.str("");

      if ('e' == colname[0]) //error column
      {
        ticker = (ticker) ? false : true;
        if (ticker)
        {
          tempPeakError.push_back(tempError);
          errorTotal += tempError;
          tempError = 0;
          ss << str;
          ss >> tempError;
        }
        else
        {
          int add;
          ss << str;
          ss >> add;
          tempError += add;

        }
      }
      else
      {
        ss << str;
        ss >> tempQuantile;
        tempVarientQuantities.push_back(tempQuantile);
        varientNames.push_back(col);
        varientTotal += tempQuantile;
      }
      ipFilehandle.ignore(256, '\n');
      delete [] str;
      str = NULL;
      delete [] colname;
      colname = NULL;
    }

    for (int i = 0; i < tempPeakError.size(); i++)
    {
      peakError.push_back(tempPeakError[i]/errorTotal);
    }

    for (int i = 0; i < tempVarientQuantities.size(); i++)
    {
      if (varientTotal > 0)
      {
        float newVarient = (float) tempVarientQuantities[i]/(float) varientTotal;
        varientQuantities.push_back(newVarient);
      }
    }

    return true;
  }
  // -------------------------------------------------------------------------
  bool LinearEquation::saveCPLEXLP(const char * filename)
  {
    if (m_lhsConstraints.size() != m_rhsConstraints.size())
    {
      ERROR_MSG("Left hand constraints and right hand constraints are unequal sizes!");
      return false;
    }
    //open delimited file
    ofstream lpFilehandle(filename, ios::binary);

    if (!lpFilehandle.is_open() || !lpFilehandle.good())
    {
      ERROR_MSG("Unable to read file! " << filename);
      return false;
    }

    lpFilehandle << "Minimize" << endl;
    lpFilehandle << "error: " << endl;

    //generate error variables
    lpFilehandle << "err" << 1 << "_pos + err" << 1 << "_neg";

    for (int i = 2; i <= m_rhsConstraints.size(); i++)
    {
      lpFilehandle << " + " << endl;
      lpFilehandle << "err" << i << "_pos + err" << i << "_neg";
    }
    lpFilehandle << endl;
    lpFilehandle << "Subject To" << endl;

    for (int i = 0; i < m_rhsConstraints.size(); i++)
    {
      int constraintIndex = i + 1;
      lpFilehandle << "c" << constraintIndex << ": ";

      lpFilehandle << "err" << constraintIndex << "_pos - err"
          << constraintIndex << "_neg";

      for (int j = 0; j < m_lhsConstraints[i].size(); j++)
      {
        if (m_lhsConstraints[i][j] != 0)
        {
          int variableIndex = j + 1;
          lpFilehandle << " + ";
          lpFilehandle << std::fixed << std::setprecision(3) << m_lhsConstraints[i][j];
          lpFilehandle << "Qv";
          lpFilehandle << setw(3) << setfill('0') << variableIndex;
        }
      }
      lpFilehandle << " = ";
      lpFilehandle << std::fixed << std::setprecision(3) << m_rhsConstraints[i];
      lpFilehandle << endl;
    }
    lpFilehandle << "end";
    lpFilehandle.close();
    return true;
  }
}
