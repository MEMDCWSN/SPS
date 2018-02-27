#include "SpecSet.h"

#include <string>

namespace specnets
{
  class ParameterList;

  bool generateRelaunchScript(ParameterList & ip);
  bool performGenoMS(ParameterList & ip);
  bool performMergeOfCSPSAndGenoMS(ParameterList & ip);
  void performGenoMSRename(std::string & statusFileName);
}

