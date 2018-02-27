#include <string>


namespace specnets
{
  class ParameterList;
  
  std::string getCurrentDirectory(void);
  std::string getProjPath(const ParameterList & pl, 
                          const std::string & addPath,
                          bool forceAbsolutePath = false);
  std::string getCurrentTimeString(void);
}

