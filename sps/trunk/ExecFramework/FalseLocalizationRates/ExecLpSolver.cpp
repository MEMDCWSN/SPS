// Module Includes
#include "ExecLpSolver.h"

// System Includes
#include <stdio.h>

using namespace specnets;
using namespace std;

namespace specnets
{
// -------------------------------------------------------------------------
  ExecLpSolver::ExecLpSolver(void) :
      ownInput(true), ownOutput(true), m_lpFilenames(0x0)
  {
    m_name = "ExecLpSolver";
    m_type = "ExecLpSolver";
  }

// -------------------------------------------------------------------------
  ExecLpSolver::ExecLpSolver(const ParameterList & params) :
      ExecBase(params), ownInput(true), ownOutput(true), m_lpFilenames(0x0)
  {
    m_name = "ExecLpSolver";
    m_type = "ExecLpSolver";
  }
// -------------------------------------------------------------------------
  ExecLpSolver::ExecLpSolver(vector<string> * lpFilenames,
                             const ParameterList & params) :
      ExecBase(params), ownInput(true), ownOutput(true), m_lpFilenames(lpFilenames)
  {
    m_name = "ExecLpSolver";
    m_type = "ExecLpSolver";
  }
// -------------------------------------------------------------------------
  ExecLpSolver::~ExecLpSolver(void)
  {
    if (ownInput)
    {
      delete m_lpFilenames;
      m_lpFilenames = 0x0;
    }
  }
// -------------------------------------------------------------------------
  ExecBase * ExecLpSolver::clone(const ParameterList & inputParams) const
  {
    return new ExecLpSolver(inputParams);
  }
// -------------------------------------------------------------------------
  bool ExecLpSolver::invoke(void)
  {

    if (m_lpFilenames == 0x0)
    {
      ERROR_MSG("Input pointers not defined!");
      return false;
    }

    int startBaseIdx;
    if (m_params.exists("IDX_START"))
    {
      startBaseIdx = max(0, m_params.getValueInt("IDX_START"));
    }
    else
    {
      startBaseIdx = 0;
    }
    DEBUG_TRACE;

    int endBaseIdx;
    if (m_params.exists("IDX_END"))
    {
      endBaseIdx = max(0, m_params.getValueInt("IDX_END"));
    }
    else
    {
      endBaseIdx = -1;
    }

    DEBUG_VAR(m_params.getValueInt("IDX_START"));
    DEBUG_VAR(m_params.getValueInt("IDX_END"));

    if (endBaseIdx == -1 || endBaseIdx > m_lpFilenames->size() - 1)
    {
      endBaseIdx = m_lpFilenames->size() - 1;
    }

    for (int i = startBaseIdx; i <= endBaseIdx; i++)
    {
      string inputFile = (*m_lpFilenames)[i];
      FilenameManager lpOutputFm(inputFile.c_str());
      lpOutputFm.extension = "txt";
      lpOutputFm.joinFilename();
      string outputFile = (lpOutputFm.filenameFull);

      string glpSolExe = m_params.getValue("GLPSOL_EXE");

      DEBUG_VAR(inputFile);
      DEBUG_VAR(outputFile);

      string lpRunCommand = glpSolExe + " --lp \"" + inputFile + "\" -o \""
          + outputFile + "\"";

      DEBUG_VAR(lpRunCommand);
      int returnValue = spsSystem(lpRunCommand.c_str());
      if (returnValue != 0)
      {
        ERROR_MSG("Unable to run command! " << lpRunCommand << " code " << returnValue);
        return false;
      }
    }
    return true;
  }
// -------------------------------------------------------------------------

  bool ExecLpSolver::saveOutputData(void)
  {
    return true;
  }
// -------------------------------------------------------------------------
  bool ExecLpSolver::saveInputData(std::vector<std::string> & filenames)
  {
    return true;
  }
// -------------------------------------------------------------------------
  bool ExecLpSolver::loadInputData(void)
  {
    if (m_lpFilenames == NULL)
    {
      m_lpFilenames = new vector<string>;
      ownInput = true;
    }

    if (m_params.exists("OUTPUT_LP_DIRECTORY"))
    {
      if (!directoryContents(m_params.getValue("OUTPUT_LP_DIRECTORY"),
                             "",
                             "lp",
                             *m_lpFilenames,
                             false))
      {
        ERROR_MSG("Unable to read from directory! " << m_params.getValue("OUTPUT_LP_DIRECTORY"));
        return false;
      }
      //append directory names to output
      for (int i = 0; i < m_lpFilenames->size(); i++)
      {
        (*m_lpFilenames)[i] = m_params.getValue("OUTPUT_LP_DIRECTORY") + '/' + (*m_lpFilenames)[i];
      }
    }

    DEBUG_VAR(m_lpFilenames->size());
    return true;
  }
// -------------------------------------------------------------------------
  bool ExecLpSolver::loadOutputData(void)
  {
    return true;
  }
// -------------------------------------------------------------------------
  vector<ExecBase*> const & ExecLpSolver::split(int numSplit)
  {
    DEBUG_VAR(numSplit);

    if (numSplit < 2)
    {
      DEBUG_MSG("Number split [" << numSplit << "] must be at least 2");
      return m_subModules;
    }
    int numLpLines = m_lpFilenames->size();
    DEBUG_VAR(numLpLines);
    if (numLpLines == 0)
    {
      DEBUG_MSG("Must have at least one lp");
      return m_subModules;
    }

    int startBaseIdx;
    if (m_params.exists("IDX_START"))
    {
      startBaseIdx = max(0, m_params.getValueInt("IDX_START"));
    }
    else
    {
      startBaseIdx = 0;
    }
    DEBUG_VAR(startBaseIdx);

    int endBaseIdx;
    if (m_params.exists("IDX_END"))
    {
      endBaseIdx = max(0, m_params.getValueInt("IDX_END"));
      endBaseIdx = min(endBaseIdx, (int)m_lpFilenames->size() - 1);
    }
    else
    {
      endBaseIdx = m_lpFilenames->size() - 1;
    }
    DEBUG_VAR(endBaseIdx);

    if (startBaseIdx == endBaseIdx)
    {
      DEBUG_MSG("Must have more than one row");
      return m_subModules;
    }

    DEBUG_TRACE;

    //find remainder;
    int numPerSplit = (int)ceil(((float)endBaseIdx - (float)startBaseIdx + 1)/ (float)numSplit);

    DEBUG_VAR(numPerSplit);

    int startIndex = startBaseIdx;
    int endIndex = startBaseIdx + numPerSplit - 1;

    bool splitAll = false;

    for (int i = 0; i < numSplit; i++)
    {
      if (splitAll)
      {
        continue;
      }
      // We have enough comparisons for the CPU
      // Copy the parameters
      ParameterList childParams(m_params);
      // But set the start and end indices
      char buf[128];
      sprintf(buf, "%d", startIndex);
      childParams.setValue("IDX_START", buf);
      sprintf(buf, "%d", endIndex);
      childParams.setValue("IDX_END", buf);

      //DEBUG_MSG("Start [" << startIndex << "] End [" << i << "] Split [" << numSplit << "]");
      ExecBase * theClone = new ExecLpSolver(m_lpFilenames,
                                        childParams);

      theClone->setName(makeName(m_name, i));

      std::string suffix("");
      char bufSplit[128];
      sprintf(bufSplit, "%d", i + 1);
      theClone->m_params.setValue("NUM_SPLIT", bufSplit);

      m_subModules.push_back(theClone);
      startIndex = endIndex + 1;
      endIndex = startIndex + numPerSplit - 1;

      if (endIndex > (int)m_lpFilenames->size() - 1)
      {
        endIndex = (int)m_lpFilenames->size() - 1;
      }
      if (endIndex < startIndex)
      {
        splitAll = true;
      }
    }
    DEBUG_VAR(m_subModules.size());
    return m_subModules;
  }
// -------------------------------------------------------------------------
  bool ExecLpSolver::merge(void)
  {
    //no need to merge, all results are produced by split
    return true;
  }
// -------------------------------------------------------------------------
  bool ExecLpSolver::validateParams(std::string & error)
  {
    m_isValid = false;

    VALIDATE_PARAM_EXIST("GLPSOL_EXE");

    m_isValid = true;
    return true;
  }
}
