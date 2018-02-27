// Module Includes
#include "ReportModuleBase.h"

// System Includes
#include <stdio.h>

using namespace std;

namespace spsReports {

// -------------------------------------------------------------------------
ReportModuleBase::ReportModuleBase(void) :
  m_isValid(false), m_name("ReportModuleBase")
{
}

// -------------------------------------------------------------------------
ReportModuleBase::ReportModuleBase(const ParameterList & params) :
  m_params(params), m_name("ReportModuleBase"), m_isValid(false)
{
}

// -------------------------------------------------------------------------
ReportModuleBase::~ReportModuleBase(void)
{
}

// -------------------------------------------------------------------------
bool ReportModuleBase::isValid(void) const
{
  return m_isValid;
}

// -------------------------------------------------------------------------
string ReportModuleBase::getName(void) const
{
  return m_name;
}

// -------------------------------------------------------------------------
string ReportModuleBase::getType(void) const
{
  return m_type;
}

// -------------------------------------------------------------------------
void ReportModuleBase::setName(std::string name)
{
  m_name = name;
}

// -------------------------------------------------------------------------
std::string ReportModuleBase::makeName(std::string baseFilename,
                               int numSplit) const
{
  std::string filename(baseFilename);
  char buf[256];
  sprintf(buf, "%d", numSplit + 1);

  filename += "_";
  filename += buf;

  return filename;
}

}; // namespace
