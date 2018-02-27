// Header Include
#include "Logger.h"
#include "ParallelProcessExecution.h"

// System Includes
#include <vector>

using namespace specnets;
using namespace std;

// -------------------------------------------------------------------------
ParallelProcessExecution::ParallelProcessExecution(ExecBase * moduleExec) :
  ParallelExecution(moduleExec)
{
  // EMPTY
}

// -------------------------------------------------------------------------
ParallelProcessExecution::~ParallelProcessExecution(void)
{
  // EMPTY
}

// -------------------------------------------------------------------------
bool ParallelProcessExecution::invoke(int numSplits)
{
  DEBUG_TRACE;

  vector<ExecBase *> const & subModules = m_moduleExec->split(numSplits);

  DEBUG_TRACE;

  // DEBUG: THIS IS A DUMMY IMPLEMENTATION TO MAKE SURE EXECUTING IN A SEPARATE PROCESS WORKS
  std::string commandLine("/home/lars/sps/trunk/ExecFramework/main_execmodule ");
  commandLine += subModules[i]->getType();
  commandLine += " ";
  commandLine += filenames[0];
  int returnValue = system(commandLine.c_str());
  DEBUG_VAR(commandLine);
  if (returnValue)
  {
    LOG_ERROR_MSG(log, "Return code from system() is " << returnValue );
    return false;
  }

  DEBUG_TRACE;

  return true;
}

