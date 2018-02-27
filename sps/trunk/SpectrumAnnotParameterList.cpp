/**
 * Parsing for annotation statistics input files
 */
#include <iostream>
#include <fstream>
#include <vector>
#include "SpectrumAnnotParameterList.h"
#include "utils.h"

/** Spectrum annotation parameter  parsing
 *	@param header - string vector of column names for input
 *	@param fields - string vector of current line split by "\t"
 */
namespace specnets
{
  // -------------------------------------------------------------------------
  bool SpectrumAnnotParameter::parseFromFields(const vector<string>& header,
                                               const vector<string>& fields)
  {
    int i;
    for (i = 0; i < fields.size(); i++)
    {
      if (header[i].compare("ion types") == 0)
      {
        ionNames = fields[i];
      }
      else if (header[i].compare("statistic") == 0)
      {
        statistic = fields[i];
      }
      else if (header[i].compare("name") == 0)
      {
        statisticName = fields[i];
      }
      else
      {
        std::cerr << "Warning: Unknown field name " << header[i] << endl;
        return 0;
      }
    }
    return true;
  }

  string SpectrumAnnotParameter::getFragSpecificIonNames(const string& fragType) const
  {
    string fragTypeUse(fragType);
    if (fragTypeUse.length() == 0)
    {
      fragTypeUse = "CID";
    }

    fragTypeUse += "=";

    unsigned int fragPos = ionNames.find(fragTypeUse);

    if (fragPos == string::npos)
    {
      return ionNames;
    }

    unsigned int startPos = ionNames.find("=", fragPos);

    if (startPos == string::npos)
    {
      ERROR_MSG("Invalid ion string \'" << ionNames << "\', coud not locate \'" << fragTypeUse << "\'");
      return ionNames;
    }

    unsigned int len = ionNames.find(";", startPos);

    if (len != string::npos)
    {
      len = len - startPos - 1;
    }

    return ionNames.substr(startPos + 1, len);
  }
  // -------------------------------------------------------------------------
  bool SpectrumAnnotParameterList::loadSpectrumAnnotFile(const char *spectrum_annot_file)
  {
    //open stats parameters
    ifstream spectra_stats_handle(spectrum_annot_file, ios::binary);

    if (!spectra_stats_handle.is_open() || !spectra_stats_handle.good())
    {
      cout << "Error: couldn't open original stats results file for reading:"
          << spectrum_annot_file << endl;
      return false;
    }

    //read fields for file
    string line_buff;
    getline(spectra_stats_handle, line_buff);

    if ('\r' == line_buff[line_buff.size() - 1])
    { //strip off carriage return for dos
      line_buff.resize(line_buff.size() - 1);
    }

    bool read_header = false;
    vector<string> field_names;
    vector<string> header;

    if (line_buff.compare("Spectrum Statistics V1") == 0)
    { //check file format
      while (!spectra_stats_handle.eof())
      {
        getline(spectra_stats_handle, line_buff);

        if ('\r' == line_buff[line_buff.size() - 1])
        { //strip off carriage return for dos
          line_buff.resize(line_buff.size() - 1);
        }

        if (!read_header)
        { //haven't yet gotten header
          if (line_buff[0] == '#' || line_buff.length() == 0)
          { //skip initial comments
            //do nothing
          }
          else
          {
            splitText(line_buff.c_str(), header, "\t");
            read_header = true;
            int i;
#ifdef DEBUG
            for (i=0; i<header.size(); i++)
            std::cout << "header " << i << "\t" << header[i] << endl;
#endif
          }
        }
        else if (line_buff.size() == 0)
        {
          continue;
        }
        else
        {
          //load in  file, header already defined
          vector<string> fields;
          SpectrumAnnotParameter stat_param;
          splitText(line_buff.c_str(), fields, "\t");
          stat_param.parseFromFields(header, fields);
          m_params.push_back(stat_param);
        }
      }
    }
    return true;
  }

  // -------------------------------------------------------------------------
  SpectrumAnnotParameter & SpectrumAnnotParameterList::operator[](unsigned int i)
  {
    return m_params[i];
  }

  // -------------------------------------------------------------------------
  const SpectrumAnnotParameter & SpectrumAnnotParameterList::operator[](unsigned int i) const
  {
    return m_params[i];
  }

  // -------------------------------------------------------------------------
  SpectrumAnnotParameterList & SpectrumAnnotParameterList::operator=(SpectrumAnnotParameterList &other)
  {
    m_params.resize(other.m_params.size());
    m_params = other.m_params;
  }

  // -------------------------------------------------------------------------
  unsigned int SpectrumAnnotParameterList::size()
  {
    return (unsigned int)m_params.size();
  }

  // -------------------------------------------------------------------------
  unsigned int SpectrumAnnotParameterList::resize(unsigned int newSize)
  {
    m_params.resize(newSize);
    return m_params.size();
  }
}

