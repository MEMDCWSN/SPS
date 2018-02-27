////////////////////////////////////////////////////////////////////////////////
// Header Include
#include "ReportModuleSpecplot.h"

// Module Includes
#include "Logger.h"



// System Includes
#include <stdio.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////
// Defines used by the module

using namespace spsReports;
using namespace std;

////////////////////////////////////////////////////////////////////////////////
ReportModuleSpecplot::ReportModuleSpecplot(void)
{
    m_name            = "specplot";
    m_type            = "specplot";
}
////////////////////////////////////////////////////////////////////////////////
ReportModuleSpecplot::ReportModuleSpecplot(const ParameterList & params)
{
    m_params          = params;
    m_name            = "specplot";
    m_type            = "specplot";
}
////////////////////////////////////////////////////////////////////////////////
ReportModuleSpecplot::~ReportModuleSpecplot(void)
{
}
////////////////////////////////////////////////////////////////////////////////
ReportModuleBase * ReportModuleSpecplot::clone(const ParameterList & input_params) const
{
  return new ReportModuleSpecplot(input_params);
}
// -------------------------------------------------------------------------
ReportModuleBase * ReportModuleSpecplot::clone(void) const
{
  return new ReportModuleSpecplot();
}
// -------------------------------------------------------------------------
bool ReportModuleSpecplot::invoke(void)
{
  return true;
}

// -------------------------------------------------------------------------
bool ReportModuleSpecplot::invoke(int argc, char ** argv)
{
  //DEBUG_MSG("Entering  ReportModuleSpecplot::invoke()");


  specplot.processOptions(argc, argv);

  specplot.plot();


  //DEBUG_MSG("Exiting  ReportModuleSpecplot::invoke()");

  return true;
}
// -------------------------------------------------------------------------
bool ReportModuleSpecplot::loadInputData(void)
{
  return true;
}
// -------------------------------------------------------------------------
bool ReportModuleSpecplot::saveOutputData(void)
{
  return true;
}
// -------------------------------------------------------------------------
bool ReportModuleSpecplot::validateParams(std::string & error)
{
  m_isValid = true;
  return true;
}
////////////////////////////////////////////////////////////////////////////////
// -------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////


