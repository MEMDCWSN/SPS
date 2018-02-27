// Header Include
#include "Logger.h"
#include "FileUtils.h"
#include "ParallelPbsExecution.h"
#include "SgeGridMonitor.h"
#include "Specific.h"

// System Includes
#include <stdlib.h>
#include <fstream>
#include <vector>

const unsigned int GRID_WAIT_SECONDS = 10;

using namespace specnets;
using namespace std;

// -------------------------------------------------------------------------
ParallelPbsExecution::ParallelPbsExecution(ExecBase * moduleExec,
                                           bool gridExecution,
                                           bool manual,
                                           bool resume) :
  ParallelExecution(moduleExec), m_gridExecution(gridExecution),
      m_manual(manual), m_resume(resume)
{
  // EMPTY
}

// -------------------------------------------------------------------------
ParallelPbsExecution::~ParallelPbsExecution(void)
{
  // EMPTY
}

// -------------------------------------------------------------------------
bool ParallelPbsExecution::invoke(int numSplits)
{
  DEBUG_TRACE;

  vector<ExecBase *> const & subModules = m_moduleExec->split(numSplits);
  vector<string> vecJobIds(subModules.size());

  string sge_path = m_moduleExec->m_params.getValue("GRID_SGE_EXE_DIR");
  DEBUG_VAR(sge_path);

  string gridType = m_moduleExec->m_params.getValue("GRID_TYPE");
  string gridParams = m_moduleExec->m_params.getValue("GRID_PARAMS");
  DEBUG_VAR(gridParams);
  int numCpus = m_moduleExec->m_params.getValueInt("GRID_NUMCPUS");
  DEBUG_VAR(numCpus);
  string dataDir = m_moduleExec->m_params.getValue("GRID_DATA_DIR");
  DEBUG_VAR(dataDir);
  if (dataDir.empty())
  {
    dataDir = ".";
  }

  DEBUG_VAR(m_gridExecution);

  if (!m_resume)
  {
    string execFilename(dataDir);
    execFilename += "/run_jobs.sh";
    DEBUG_VAR(execFilename);
    ofstream ofs_exec(execFilename.c_str(), ios::out | ios::binary);
    if (!ofs_exec)
    {
      ERROR_MSG("Can not open  " << execFilename);
      return false;
    }

    SgeGridMonitor sgm;

    for (int i = 0; i < subModules.size(); i++)
    {
      // Save the input data will return the files we need to copy to target machine
      vector<std::string> filenamesIn;
      subModules[i]->saveInputData(filenamesIn);
      if (filenamesIn.size() == 0)
      {
        ERROR_MSG("No files returned from saveInputData");
        return false;
      }

      string jobFilename(dataDir);
      jobFilename += "/";
      jobFilename += subModules[i]->getName();
      jobFilename += "_job.sh";

      DEBUG_VAR(jobFilename);
      ofstream ofs(jobFilename.c_str(), ios::out | ios::binary);
      if (!ofs)
      {
        ERROR_MSG("Can not open  " << jobFilename);
        return false;
      }
      ofs << "#!/bin/bash" << endl;
      ofs << "#" << endl;
      ofs << "#$ -cwd" << endl;
      ofs << "#$ -j y" << endl;
      ofs << "#$ -S /bin/bash" << endl;
      ofs << "#" << endl;

      ofs << m_moduleExec->m_params.getValue("GRID_EXE_DIR");
      ofs << "/main_execmodule " << subModules[i]->getType() << " "
          << filenamesIn[0] << " ";

      if (numCpus > 1)
      {
        ofs << "-t ";
        char buf[128];
        sprintf(buf, "%d", numCpus);
        ofs << buf;
      }
      ofs << "&" << endl;
      ofs << "wait" << endl;
      ofs.close();

      string vmem = m_moduleExec->m_params.getValue("GRID_VMEM");

      if (m_gridExecution)
      {
        string newJobId = sgm.submitJob(sge_path, jobFilename, gridParams, gridType);

        if (!newJobId.empty())
        {
          DEBUG_VAR(newJobId);

          SgeGridMonitor::JobStatus js = sgm.getJobStatus(newJobId);
          if (js == SgeGridMonitor::JS_ERROR)
          {
            ERROR_MSG("Submitting job");
            return false;
          }
        }
        else
        {
          ERROR_MSG("Problem creating job");
          return false;
        }
        vecJobIds[i] = newJobId;
      }

      //if (!m_gridExecution)
      {
        ofs_exec << "qsub " << gridParams << " " << jobFilename << endl;
      }
    }

    //DEBUG_TRACE;
    ofs_exec.close();

    string command("chmod a+x " + execFilename);
    //DEBUG_TRACE;
    //system(command.c_str());

    if (!m_gridExecution)
    {
      return true;
    }

    // Wait for the processing to be done (this could take a while)
    //   The method we use is to wait for all the output files to exist and stop changing
    bool done = false;
    while (!done)
    {
      // Give the files a chance to change
      sleep(GRID_WAIT_SECONDS);

      sgm.refreshInfo(sge_path);

      int jobsRemaining = 0;
      for (int i = 0; i < subModules.size(); i++)
      {
        if (sgm.getJobStatus(vecJobIds[i]) != SgeGridMonitor::JS_UNKNOWN)
        {
          jobsRemaining++;
        }
      }

      if (jobsRemaining != 0)
      {
        DEBUG_MSG(jobsRemaining << " jobs still in queue");
        continue;
      }

      done = true;
    }
  }

  if (m_resume || m_gridExecution)
  {
    //DEBUG_TRACE;
    for (int i = 0; i < subModules.size(); i++)
    {
      DEBUG_MSG("Loading results for module " << i);
      // Load back in the results
      subModules[i]->loadOutputData();
    }

    //DEBUG_TRACE;

    // Merge back the results
    int returnValue = m_moduleExec->merge();
    if (!returnValue)
    {
      ERROR_MSG("Problem merging results from sub-modules" );
      return false;
    }
    //DEBUG_TRACE;
    return true;
  } // if (!m_resume)

  //DEBUG_TRACE;
  return true;
}

