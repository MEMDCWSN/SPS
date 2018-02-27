// Header Include
#include "CommandLineParser.h"
#include "Logger.h"

// System Includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DEBUG_OPTIONS 0

using namespace std;

namespace specnets
{

  //----------------------------------------------------------------------------
  CommandLineParser::CommandLineParser(int argc,
                                       char ** argv,
                                       int numNonOptions,
                                       const std::vector<Option> & listOptions)
  {
    for (size_t i = numNonOptions + 1; i < argc; i++)
    {
      if (DEBUG_OPTIONS) DEBUG_VAR(i);
      if (DEBUG_OPTIONS) DEBUG_VAR(argv[i]);

      bool unknownOption = true;
      if (DEBUG_OPTIONS) DEBUG_VAR(listOptions.size());
      for (size_t j = 0; j < listOptions.size(); j++)
      {
        if (DEBUG_OPTIONS) DEBUG_VAR(j);

        string option("-");
        option += listOptions[j].m_flag;
        if (DEBUG_OPTIONS) DEBUG_VAR(option);

        if (::strcmp(argv[i], option.c_str()) == 0)
        {
          if (DEBUG_OPTIONS) DEBUG_TRACE;
          unknownOption = false;
          if (DEBUG_OPTIONS) DEBUG_VAR(listOptions[j].m_numValues);
          if (listOptions[j].m_numValues != 0)
          {
            if (DEBUG_OPTIONS) DEBUG_TRACE;
            for (unsigned int v = 0; v < listOptions[j].m_numValues; v++)
            {
              i++;
              if (i == argc)
              {
                if (DEBUG_OPTIONS) DEBUG_TRACE;
                char buf[256];
                sprintf(buf, "%u", listOptions[j].m_numValues);
                m_error = "Option [" + option + "] requires [" + buf + "] argument(s)";
                return;
              }
              else
              {
                if (DEBUG_OPTIONS) DEBUG_VAR(listOptions[j].m_name);
                if (DEBUG_OPTIONS) DEBUG_VAR(argv[i]);
                m_listOptions.push_back(make_pair<string, string> (listOptions[j].m_name,
                                                                   argv[i]));
              }
            }
          }
          else
          {
            if (DEBUG_OPTIONS) DEBUG_VAR(listOptions[j].m_name);
            m_listOptions.push_back(make_pair<string, string> (listOptions[j].m_name,
                                                               ""));
          }
        }

      } // for (size_t j = 0; j < listOptions.size(); j++)

      if (unknownOption)
      {
        // Check for "automatic" parameters
        string theArg = argv[i];
        if (DEBUG_OPTIONS) DEBUG_VAR(theArg);
        if (theArg.find("-auto") == 0)
        {
          unknownOption = false;
          string theOption = theArg.substr(5);
          if (DEBUG_OPTIONS) DEBUG_VAR(theOption);
          if (theOption.length() == 0) {
            if (DEBUG_OPTIONS) DEBUG_TRACE;
            m_error = "Automatic option requires additional characters beyond [auto]";
            return;
          }
          i++;
          if (i == argc)
          {
            m_error = "Automatic option requires [1] argument";
            return;
          }
          else
          {
            string theValue = argv[i];
            if (DEBUG_OPTIONS) DEBUG_VAR(theValue);
            m_listOptions.push_back(make_pair<string, string> (theOption,
                                                               theValue));
                                                               
          }
        }
      }
      
      if (unknownOption)
      {
        m_error = "Unknown option [";
        m_error += argv[i];
        m_error += "]";
        return;
      }

    } // for (size_t i = numNonOptions + 1; i < argc; i++)

    return;
  }

  CommandLineParser::~CommandLineParser()
  {
    // EMPTY
  }

  //----------------------------------------------------------------------------
  void CommandLineParser::getOptionsAsParameterList(ParameterList & params) const
  {
    for (size_t i = 0; i < m_listOptions.size(); i++)
    {
      //DEBUG_VAR(m_listOptions[i].first);
      //DEBUG_VAR(m_listOptions[i].second);
      params.setValue(m_listOptions[i].first, m_listOptions[i].second);
    }
    return;
  }

  //----------------------------------------------------------------------------
  bool CommandLineParser::validate(string & error) const
  {
    if (m_error.empty())
    {
      return true;
    }
    error = m_error;
    return false;
  }

} //namespace specnets
