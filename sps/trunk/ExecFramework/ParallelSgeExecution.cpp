// Header Include
#include "Logger.h"
#include "FileUtils.h"
#include "ParallelSgeExecution.h"

// System Includes
#include <stdlib.h>
#include <fstream>
#include <vector>

const unsigned int GRID_WAIT_SECONDS = 10;

using namespace specnets;
using namespace std;

// -------------------------------------------------------------------------
ParallelSgeExecution::ParallelSgeExecution(ExecBase * moduleExec,
                                           bool gridExecution,
                                           bool manual,
                                           bool resume,
                                           bool debug) :
  ParallelExecution(moduleExec), m_gridExecution(gridExecution),
      m_manual(manual), m_resume(resume), m_debug(debug)
{
  // EMPTY
}

// -------------------------------------------------------------------------
ParallelSgeExecution::~ParallelSgeExecution(void)
{
  // EMPTY
}

// -------------------------------------------------------------------------
bool ParallelSgeExecution::invoke(int numSplits)
{
  DEBUG_TRACE;

  vector<ExecBase *> const & subModules = m_moduleExec->split(numSplits);

  string sge_path = m_moduleExec->m_params.getValue("GRID_SGE_EXE_DIR");
  if (sge_path.empty())
  {
    sge_path = ".";
  } 

  string execFilename("run_jobs.sh");
  ofstream ofs_exec(execFilename.c_str(), ios::out);
  if (!ofs_exec)
  {
    ERROR_MSG("Can not open:  " << execFilename);
    return false;
  }

  string gridParams = m_moduleExec->m_params.getValue("GRID_PARAMS");
  int numCpus = m_moduleExec->m_params.getValueInt("GRID_NUMCPUS");    
  string dataDir = m_moduleExec->m_params.getValue("GRID_DATA_DIR_OUT");
  if (dataDir.empty())
  {
    dataDir = ".";
  }

  if (!m_resume)
  {
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

      if (!m_debug) {
        string jobFilename(dataDir);
        jobFilename += "/";
        jobFilename += subModules[i]->getName();
        jobFilename += "_job.sh";
        
        DEBUG_VAR(jobFilename);

        ofstream ofs(jobFilename.c_str(), ios::out);
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
            
        //string command("chmod a+x " + jobFilename);
        string command("chmod a+x "); command += jobFilename;
        //string command("chmod 755 "); command += jobFilename;
        int returnValue = system(command.c_str());

        if (!m_gridExecution)
        {
          ofs_exec << "qsub " << gridParams << " " << jobFilename << endl;
        }
      } else {
        string command = m_moduleExec->m_params.getValue("GRID_EXE_DIR");
        char buf[128];
        sprintf(buf, "%d", (int)i+1);
      
        command += "/main_execmodule ";
        command += subModules[i]->getType();
        command += " ";
        command += filenamesIn[0];
        command += " -lf ";
        command += dataDir;
        command += "/";
        command += subModules[i]->getType();
        command += "_log_";
        command += buf;
        command += ".txt ";
        if (numCpus > 1)
        {
          command += "-t ";
          char buf[128];
          sprintf(buf, "%d", numCpus);
          command += buf;
        }           
        
        DEBUG_VAR(command);
        if (system(command.c_str())) {
          ERROR_MSG("Executing debug jobs!");
          return false;
        }
      }
      
    } // for (int i = 0; i < numSplit; i++)

    if (m_gridExecution && !m_debug)
    {
      ofs_exec << "#!/bin/bash" << endl;
      ofs_exec << "#" << endl;
      ofs_exec << "#$ -cwd" << endl;
      ofs_exec << "#$ -j y" << endl;
      ofs_exec << "#$ -S /bin/bash" << endl;
      ofs_exec << "#" << endl;
      ofs_exec << ". $(echo ";
      ofs_exec << dataDir;
      ofs_exec << "/$1\"_\"$SGE_TASK_ID\"_job.sh\")";
    }
    
    if (!m_debug) {
      DEBUG_TRACE;
      ofs_exec.close();

      string command("chmod a+x ") ; command += execFilename;
      DEBUG_TRACE;
      int returnValue = system(command.c_str());
    }
    
    DEBUG_MSG("Submodules: " << subModules.size());

    DEBUG_VAR(m_gridExecution);
    if (!m_gridExecution || subModules.size() == 0)
    {
      return true;
    }

    if (!m_debug) {
      char buf[128];
      sprintf(buf, "%d", (int)subModules.size());

      string command2(sge_path);
      command2 += "/qsub -sync yes ";
      command2 += gridParams;
      command2 += " -t 1-";
      command2 += buf;
      command2 += " ";
      command2 += execFilename;
      command2 += " ";
      command2 += subModules[0]->getType();  // all the same type.. use first
      DEBUG_VAR(command2);
      
      int returnValue = system(command2.c_str());
      DEBUG_VAR(returnValue);
    
      if (returnValue != 0) {
        ERROR_MSG("Executing batch jobs!");
        return false;
      }
    } // if (!m_debug) {
    
  } // if (!m_resume)

  if (m_resume || m_gridExecution)
  {
    DEBUG_TRACE;
    for (int i = 0; i < subModules.size(); i++)
    {
      DEBUG_MSG("Loading results for module " << i);
      // Load back in the results
      if (!subModules[i]->loadOutputData())
      {
        DEBUG_MSG("Could not load results for module " << i);
        return false;
      }
    }

    DEBUG_TRACE;

    // Merge back the results
    int returnValue = m_moduleExec->merge();
    if (!returnValue)
    {
      ERROR_MSG("Problem merging results from sub-modules" );
      return false;
    }
    DEBUG_TRACE;
    return true;
  } // if (!m_resume)

  DEBUG_TRACE;
  return true;
}

