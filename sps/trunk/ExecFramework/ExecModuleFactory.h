#ifndef __ExecModuleFactory_H__
#define __ExecModuleFactory_H__

//----------------------------------------------------------------------------
//  This header contains both the factory methods and the method to register 
//    all SPS module types with the SpsModuleFactory. There is no good way 
//    to do this outside of dynamic loading which seems overkill for what we are 
//    trying to do. So creators of new modules must manually add them here and 
//    the recompile 'spsmodulexec'.
//  NOTE: All SpsModules MUST be added to the RegisterAllModules() implementation
//        in order to be run using 'spsmodulexec' (for remote execution)
//----------------------------------------------------------------------------

#include "ExecBase.h"

// System Includes
#include <map>
#include <string>

namespace specnets
{
  class ExecModuleFactory
  {
  public:
    static void cleanup(void);

    //! Registers a single module by name
    static void Register(std::string moduleName, ExecBase * module);

    //! Registers all the known modules with the factory
    static void RegisterAllModules(void);

    //! Retrieves a module from the factory by its name
    static ExecBase * getModule(std::string moduleName);

  private:
    static std::map<std::string, ExecBase *> & getMap(void);

  }; // class ExecModuleFactory

} // namespace specnets

#endif // __ExecModuleFactory_H__
