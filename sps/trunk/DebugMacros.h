#ifndef DEBUG_MACROS_H_
#define DEBUG_MACROS_H_

#include <iostream>

#define DEBUG_TRACE std::cout << "DEBUG: " << __FILE__ << " : " << __LINE__ << std::endl; std::cout.flush();
#define DEBUG_MSG( msg ) std::cout << "DEBUG: " << __FILE__ << " : " << __LINE__ <<  " : " << msg << std::endl; std::cout.flush();
#define DEBUG_VAR( var ) std::cout << "DEBUG: " << __FILE__ << " : " << __LINE__ <<  " : " << #var << " = " << var << endl; std::cout.flush();

#define ERROR_TRACE std::cout << "ERROR: " << __FILE__ << " : " << __LINE__ << endl; cout.flush();
#define ERROR_MSG( msg ) std::cout << "ERROR: " << __FILE__ << " : " << __LINE__ <<  " : " << msg << endl; cout.flush();
#define ERROR_VAR( var ) std::cout << "ERROR: " << __FILE__ << " : " << __LINE__ <<  " : " << #var << " = " << var << endl; std::cout.flush();

#define WARN_MSG( msg ) std::cout << "WARNING: " << __FILE__ << " : " << __LINE__ <<  " : " << msg << std::endl; std::cout.flush();

#endif /* DEBUG_MACROS_H_ */
