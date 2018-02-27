////////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_MODULE_FACTORY_H__
#define __REPORT_MODULE_FACTORY_H__
////////////////////////////////////////////////////////////////////////////////

//----------------------------------------------------------------------------
//  This header contains both the factory methods and the method to register
//    all SPS module types with the SpsModuleFactory. There is no good way
//    to do this outside of dynamic loading which seems overkill for what we are
//    trying to do. So creators of new modules must manually add them here and
//    the recompile 'spsmodulexec'.
//  NOTE: All SpsModules MUST be added to the RegisterAllModules() implementation
//        in order to be run using 'spsmodulexec' (for remote execution)
//----------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
#include "ReportModuleBase.h"

// System Includes
#include <map>
#include <string>

////////////////////////////////////////////////////////////////////////////////

namespace spsReports {


  class ReportModuleFactory {

   public:

    static void cleanup(void);

    //! Registers a single module by name
    static void RegisterObject(std::string moduleName, ReportModuleBase * module);

    //! Registers all the known modules with the factory
    static void RegisterAllModules(void);

    //! Retrieves a module from the factory by its name
    static ReportModuleBase * getModule(std::string moduleName);

   private:

    static std::map<std::string, ReportModuleBase *> & getObjectMap(void);

  };
////////////////////////////////////////////////////////////////////////////////
} // namespace spsReports
////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////
