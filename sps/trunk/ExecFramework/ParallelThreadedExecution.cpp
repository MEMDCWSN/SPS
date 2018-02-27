// Header Include
#include "ExecBase.h"
#include "Logger.h"
#include "ParallelThreadedExecution.h"
#include "omp.h"

// System Includes
#include <vector>

using namespace specnets;
using namespace std;

// -------------------------------------------------------------------------
ParallelThreadedExecution::ParallelThreadedExecution(ExecBase * moduleExec) :
    ParallelExecution(moduleExec)
{
  // EMPTY
}

// -------------------------------------------------------------------------
ParallelThreadedExecution::~ParallelThreadedExecution(void)
{
  // EMPTY
}

// -------------------------------------------------------------------------
bool ParallelThreadedExecution::invoke(int numSplits)
{
  DEBUG_VAR(numSplits);

  vector<ExecBase *> const & subModules = m_moduleExec->split(numSplits);

  DEBUG_VAR(subModules.size());

  omp_set_num_threads(numSplits);

#pragma omp parallel for
  for (int i = 0; i < subModules.size(); i++)
  {
    subModules[i]->invoke();
  }

  DEBUG_MSG("All threads completed");

  // Merge back the results
  m_moduleExec->merge();

  DEBUG_TRACE;

  for (int i = 0; i < subModules.size(); i++)
  {
    delete subModules[i];
  }

  return true;
}

