// Header Include
#include "ReportModuleContplot.h"

// Module Includes
#include "Logger.h"
//#include "base64.h"


// System Includes
#include <stdio.h>
#include <string.h>

// Defines used by the module


using namespace spsReports;
using namespace std;

// -------------------------------------------------------------------------
ReportModuleContplot::ReportModuleContplot(void)
{
    m_name            = "contplot";
    m_type            = "contplot";
}
// -------------------------------------------------------------------------
ReportModuleContplot::ReportModuleContplot(const ParameterList & params)
{
    m_params          = params;
    m_name            = "contplot";
    m_type            = "contplot";
}
// -------------------------------------------------------------------------
ReportModuleContplot::~ReportModuleContplot(void)
{
}
// -------------------------------------------------------------------------
ReportModuleBase * ReportModuleContplot::clone(const ParameterList & input_params) const
{
  return new ReportModuleContplot(input_params);
}
// -------------------------------------------------------------------------
ReportModuleBase * ReportModuleContplot::clone(void) const
{
  return new ReportModuleContplot();
}
// -------------------------------------------------------------------------
bool ReportModuleContplot::invoke(void)
{
  return true;
}

// -------------------------------------------------------------------------
bool ReportModuleContplot::invoke(int argc, char ** argv)
{
  //DEBUG_MSG("Entering  ReportModuleContplot::invoke()");

  // process the options for the image
  contplot.processOptions(argc, argv);
  // execute (generate the image)
  //contplot.plot();

/*
  string str = argv[0];
  //convert path
  size_t found;
  found = str.find_last_of("/\\");
  string exePath;
  exePath = '.';
  if(found != string::npos)
    exePath = str.substr(0, found);

  if(exePath[exePath.length()-1] != '/')
    exePath += '/';

  string fn = exePath + "contig";
  for(int i = 1 ; i < argc ; i++) {
    str = argv[i];
    found=str.find("zoom");
    if(found!=string::npos) {
      fn += 's';
      break;
    }
  }

  fn += ".png";

  // load the image, and store it in internal image
  int length;
  char * buffer;

  ifstream is;
  is.open (fn.c_str(), ios::binary );

  // check for file open
  if(!is.is_open())
    return false;

  // get length of file:
  is.seekg (0, ios::end);
  length = is.tellg();
  is.seekg (0, ios::beg);

  // allocate memory:
  buffer = new char [length+1];

  // read data as a block:
  is.read(buffer, length);
  is.close();
  // terminator
  buffer[length] = 0;

  // if uuencoded, encode
  m_image = base64_encode(reinterpret_cast<const unsigned char*>(buffer), length);
*/

  //DEBUG_MSG("Exiting  ReportModuleContplot::invoke()");

  return true;
}
// -------------------------------------------------------------------------
bool ReportModuleContplot::loadInputData(void)
{
  return true;
}
// -------------------------------------------------------------------------
bool ReportModuleContplot::saveOutputData(void)
{
  return true;
}
// -------------------------------------------------------------------------
bool ReportModuleContplot::validateParams(std::string & error)
{
  m_isValid = true;
  return true;
}
////////////////////////////////////////////////////////////////////////////////
// -------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////


