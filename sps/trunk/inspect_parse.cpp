/** Parsing and reformatting for inspect results
 */
#include "inspect_parse.h"

namespace specnets
{
  //-----------------------------------------------------------------------------

  bool InspectResultsSet::getResultByScan(int scan,
                                     vector<InspectResultsLine>& output_results)
  {
    vector<InspectResultsLine>::iterator it;

    bool found_results = false;
    for (it = results.begin(); it != results.end(); it++)
    {
      if (it->scan == scan)
      {
        output_results.push_back(*it);
        found_results = true;
      }
    }
    return found_results;
  }

  //-----------------------------------------------------------------------------
  bool InspectResultsSet::loadInspectResultsFile(const char *inspect_results_file)
  {
    //open Inspect file
    ifstream inspect_results_handle(inspect_results_file, ios::binary);

    if (!inspect_results_handle.is_open() || !inspect_results_handle.good())
    {
      cout << "Error: couldn't open original inspect results file for reading:"
          << inspect_results_file << endl;
      return false;
    }

    //read fields for Inspect file
    char line_buff[1024];
    inspect_results_handle.getline(line_buff, 1024);

    int i = 0;
    while(line_buff[i] != '\0')
    {
      if (line_buff[i] == '\r')
      {
        line_buff[i] = '\0';
      }
      else
      {
        i++;
      }
    }

    bool read_line = true;
    vector<string> field_names;

    if (line_buff[0] != '#')
    {
      read_line = false;
    }
    else
    {
      splitText(line_buff, field_names, (const char*)"\t");

      int i;

      //for (i=0; i<field_names.size(); i++)
      //std::cout << i << "\t" << field_names[i] << endl;
    }

    int curr_line = 0;
    //load in InSpect file.
    while (!inspect_results_handle.eof() && inspect_results_handle.good())
    {
      vector<string> fields;

      if (read_line)
      {
        inspect_results_handle.getline(line_buff, 1024);
        if (inspect_results_handle.gcount() < 5)
          continue;
      }
      else
      {
        read_line = true;
      }


      splitText(line_buff, fields, (const char*)"\t");


      InspectResultsLine res;

      res.parseFromFields(fields, field_names);

      //add line to results vector
      results.push_back(res);
    }
    return true;
  }

  //-----------------------------------------------------------------------------
  void InspectResultsSet::inspectToSpecNets(string &inspect, string &specnets)
  {
    int inspect_length = inspect.length();
    bool modification = false;
    specnets = "*.";

    const char * ptr_inspect = inspect.c_str();

    int i;
    //note: inspect_length - 3 means we ignore the trailing characters around annotation
    for (i = 2; i <= inspect_length - 3; i++)
    {
      char curr_char = ptr_inspect[i];
      if ('p' == curr_char)
      {
        //phosphorylation
        specnets += "[80]";
      }
      else if ('+' == curr_char || '-' == curr_char)
      {
        modification = true;
        specnets += '[';
        if ('-' == curr_char)
          specnets += curr_char; //only add - mod. + isn't included in SpecNets
      }
      else if (modification)
      {
        if (0 == curr_char) // check for null
        {
          // do nothing
        }
        else if ('9' >= curr_char || '.' == curr_char)
        { //this is a number (or decimal)
          specnets += curr_char;
        }
        else
        {
          specnets += ']';
          specnets += curr_char;
          modification = false;
        }
      }
      else if ('a' > curr_char) //this is a capital letter
      {
        specnets += curr_char;
      }
    }

    if (modification)
    {
      specnets += ']';
    }

    specnets += ".*";

  }

  //-----------------------------------------------------------------------------
  bool InspectResultsLine::parseFromFields(const vector<string>& fields,
                                           const vector<string>& field_names)
  {
    /* if (fields.size() != 23)
     {
     cout << "Error: inspect results line has " << fields.size()
     << ", expecting 20" << endl;
     return false;
     }
     */

    for (int i = 0; i < field_names.size(); i++)
    {
      if (field_names[i].compare("#SpectrumFile") == 0)
      {
        SpectrumFile = fields[i];
      }
      else if (field_names[i].compare("Scan#") == 0)
      {
        if (sscanf(fields[i].c_str(), "%d", &scan) != 1 || scan < 0 || scan
            > 100000000)
        {
          std::cerr << "ERROR! scan" << endl;
          return false;
        }
      }
      else if (field_names[i].compare("Annotation") == 0)
      {
        Annotation = fields[i];
      }
      else if (field_names[i].compare("Protein") == 0)
      {
        Protein = fields[i];
      }
      else if (field_names[i].compare("Charge") == 0)
      {
        if (sscanf(fields[i].c_str(), "%d", &Charge) != 1 || Charge < 0
            || Charge > 20)
        {
          std::cerr << "ERROR! Charge" << endl;
          return false;
        }
      }
      else if (field_names[i].compare("MQScore") == 0)
      {
        if (sscanf(fields[i].c_str(), "%f", &MQScore) != 1 || MQScore
            < -999999999 || MQScore > 999999999)
        {
          std::cerr << "ERROR! MQScore" << endl;
          return false;
        }
      }
      else if (field_names[i].compare("Length") == 0)
      {
        if (sscanf(fields[i].c_str(), "%d", &Length) != 1 || Length < 1
            || Length > 999999999)
        {
          std::cerr << "ERROR! Length" << endl;
          return false;
        }
      }
      else if (field_names[i].compare("TotalPRMScore") == 0)
      {
        if (sscanf(fields[i].c_str(), "%f", &TotalPRMScore) != 1
            || TotalPRMScore < -999999999 || TotalPRMScore > 999999999)
        {
          std::cerr << "ERROR! TotalPRMScore" << endl;
          return false;
        }
      }
      else if (field_names[i].compare("MedianPRMScore") == 0)
      {
        if (sscanf(fields[i].c_str(), "%f", &MedianPRMScore) != 1
            || MedianPRMScore < -999999999 || MedianPRMScore > 999999999)
        {
          std::cerr << "ERROR! MedianPRMScore" << endl;
          return false;
        }
      }
      else if (field_names[i].compare("FractionY") == 0)
      {
        if (sscanf(fields[i].c_str(), "%f", &FractionY) != 1 || FractionY < 0
            || FractionY > 1000)
        {
          std::cerr << "ERROR! FractionY" << endl;
          return false;
        }
      }
      else if (field_names[i].compare("FractionB") == 0)
      {
        if (sscanf(fields[i].c_str(), "%f", &FractionB) != 1 || FractionB < 0
            || FractionB > 1000)
        {
          std::cerr << "ERROR! FractionB" << endl;
          return false;
        }
      }
      else if (field_names[i].compare("Intensity") == 0)
      {
        if (sscanf(fields[i].c_str(), "%f", &Intensity) != 1 || Intensity < 0)
        {
          std::cerr << "ERROR! Intensity" << endl;
          return false;
        }
      }
      else if (field_names[i].compare("NTT") == 0)
      {
        if (sscanf(fields[i].c_str(), "%d", &NTT) != 1 || NTT < 0 || NTT > 3)
        {
          std::cerr << "ERROR! NTT" << endl;
          return false;
        }
      }
      else if (field_names[i].compare("F_Score") == 0)
      {
        if (sscanf(fields[i].c_str(), "%f", &F_Score) != 1)
        {
          std::cerr << "ERROR! F_Score" << endl;
          return false;
        }
      }
      else if (field_names[i].compare("p-value") == 0)
      {
        if (sscanf(fields[i].c_str(), "%f", &p_value) != 1)
        {
          std::cerr << "ERROR! p-value" << endl;
          return false;
        }
      }
      else if (field_names[i].compare("DeltaScore") == 0)
      {
        if (sscanf(fields[i].c_str(), "%f", &DeltaScore) != 1)
        {
          std::cerr << "ERROR! DeltaScore" << endl;
          return false;
        }
      }
      else if (field_names[i].compare("DeltaScoreOther") == 0)
      {
        if (sscanf(fields[i].c_str(), "%f", &DeltaScoreOther) != 1)
        {
          std::cerr << "ERROR! DeltaScoreOther" << endl;
          return false;
        }
      }
      else if (field_names[i].compare("RecordNumber") == 0)
      {
        if (sscanf(fields[i].c_str(), "%d", &RecordNumber) != 1)
        {
          std::cerr << "ERROR! RecordNumber" << endl;
          return false;
        }
      }
      else if (field_names[i].compare("DBFilePos") == 0)
      {
        if (sscanf(fields[i].c_str(), "%d", &DBFilePos) != 1)
        {
          std::cerr << "ERROR! DBFilePos" << endl;
          return false;
        }
      }
      else if (field_names[i].compare("SpecFilePos") == 0)
      {
        if (sscanf(fields[i].c_str(), "%d", &SpecFilePos) != 1)
        {
          std::cerr << "ERROR! SpecFilePos" << endl;
          return false;
        }
      }
    }

    Score = MQScore;

    const int ann_length = Annotation.length();

    if ((Annotation[1] != '.') || (Annotation[ann_length - 2] != '.'))
    {
      cerr << "Error: bad annotation format: " << Annotation << endl;
      cerr << "Expecting X.XXXXXXXXX.X" << endl;
      return false;
    }

    return true;
  }
}
