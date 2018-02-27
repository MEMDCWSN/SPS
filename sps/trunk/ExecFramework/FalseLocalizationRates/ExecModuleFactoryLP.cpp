// Module Includes
#include "Logger.h"
#include "ExecModuleFactoryLP.h"

//----------------------------------------
// Add new module include files here
#include "ExecPairSpectra.h"
#include "ExecFlr.h"
#include "ExecLpSolver.h"
//----------------------------------------

using namespace std;
using namespace specnets;


//----------------------------------------------------------------------------
namespace specnets
{
  void ExecModuleFactoryLP::cleanup(void) {
    map<string, ExecBase *> & moduleMap = getMap();
    //cout << "cur_size = " << moduleMap.size() << "\n";
    for (map<string, ExecBase *>::iterator modIt = moduleMap.begin(); modIt != moduleMap.end(); modIt++) {
      //cout << "deleting " << modIt->first << "(" << modIt->second << ")\n";
      delete modIt->second;
    }
  }
  // Registers a singe exec module with the factory
  void ExecModuleFactoryLP::Register( std::string moduleName, ExecBase * module )
  {
    map<string, ExecBase *> & moduleMap = getMap();
    moduleMap[moduleName] = module;
    //cout << "registering " << moduleName << "(" << module << ")\n";
  };

  // Registers all the known modules with the factory
  void ExecModuleFactoryLP::RegisterAllModules( void )
  {
    //--------------------------------------------------------
    // Add new module registrations here
    //--------------------------------------------------------
    Register( ExecPairSpectra().getName(), new ExecPairSpectra() );
    Register( ExecFlr().getName(), new ExecFlr() );
    Register( ExecLpSolver().getName(), new ExecLpSolver());
  };

  // Retrieves a module from the factory by its name
  ExecBase * ExecModuleFactoryLP::getModule( std::string moduleName )
  {
    map<string, ExecBase *> & moduleMap = getMap();
    if ( moduleMap.find( moduleName ) != moduleMap.end() )
    {
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
  map<string, ExecBase *> & ExecModuleFactoryLP::getMap( void )
  {
    static map<string, ExecBase *> theMap;
    return theMap;
  }

}; // namespace specnets


