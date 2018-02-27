////////////////////////////////////////////////////////////////////////////////
// Module Includes
#include "Logger.h"
#include "ReportModuleFactory.h"

//----------------------------------------
// Add new module include files here
#include "ReportModuleSpecplot.h"
#include "ReportModuleContplot.h"

////////////////////////////////////////////////////////////////////////////////

using namespace std;

////////////////////////////////////////////////////////////////////////////////
namespace spsReports {

  void ReportModuleFactory::cleanup(void) {
    map<string, ReportModuleBase *> & moduleMap = getObjectMap();
    //cout << "cur_size = " << moduleMap.size() << "\n";
    for (map<string, ReportModuleBase *>::iterator modIt = moduleMap.begin(); modIt != moduleMap.end(); modIt++) {
      //cout << "deleting " << modIt->first << "(" << modIt->second << ")\n";
      delete modIt->second;
    }
  }
  // Registers a singe exec module with the factory
  void ReportModuleFactory::RegisterObject( std::string moduleName, ReportModuleBase * module )
  {
    map<string, ReportModuleBase *> & moduleMap = getObjectMap();
    moduleMap[moduleName] = module;
    //cout << "registering " << moduleName << "(" << module << ")\n";
  };

  // Registers all the known modules with the factory
  void ReportModuleFactory::RegisterAllModules( void )
  {
    //--------------------------------------------------------
    // Add new module registrations here
    //--------------------------------------------------------
    RegisterObject( ReportModuleSpecplot().getName(), new ReportModuleSpecplot() );
    RegisterObject( ReportModuleContplot().getName(), new ReportModuleContplot() );
  };


  // Retrieves a module from the factory by its name
  ReportModuleBase * ReportModuleFactory::getModule( std::string moduleName )
  {
    map<string, ReportModuleBase *> & moduleMap = getObjectMap();
    if ( moduleMap.find( moduleName ) != moduleMap.end() )  {
      return moduleMap[moduleName];
    }

    return 0;
  };

  //----------------------------------------------------------------------------
  // NOTE: One might be tempted to remove this method and define the map as
  //       a static variable. Don't! C++ does not define when static variables
  //       will be instatiated. By making it a method static it is guaranteed
  //       to be instantiated when the method is first called.
  //----------------------------------------------------------------------------
  static map<string, ReportModuleBase *> objectMap;

  map<string, ReportModuleBase *> & ReportModuleFactory::getObjectMap( void )
  {
//    static map<string, ReportModuleBase *> objectMap;
    return objectMap;
  }

////////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
////////////////////////////////////////////////////////////////////////////////


