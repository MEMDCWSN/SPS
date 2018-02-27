//Header Include
#include "ExecGenoMS.h"

// Module Includes
#include "ExecBase.h"
#include "utils.h"
#include "Logger.h"

// System Includes
#include <stdio.h>
#include <string.h>

// Other Includes
#include "Specific.h"


using namespace specnets;
using namespace std;

// -------------------------------------------------------------------------
ExecGenoMS::ExecGenoMS(void)
{
  m_name = "ExecGenoMS";
  m_type = "ExecGenoMS";
  m_dataOwned = true;
  m_isValid = false;
}

// -------------------------------------------------------------------------
ExecGenoMS::ExecGenoMS(const ParameterList & params)
{
  m_name = "ExecGenoMS";
  m_type = "ExecGenoMS";
  m_dataOwned = true;
  m_isValid = false;
  m_params = params;
}

// -------------------------------------------------------------------------
ExecGenoMS::~ExecGenoMS(void)
{
  if(m_dataOwned)
    {
      //delete shit
    }

}

// -------------------------------------------------------------------------
bool ExecGenoMS::isValid(void) const
{
  return m_isValid;
}

// -------------------------------------------------------------------------
string ExecGenoMS::getName(void) const
{
  return m_name;
}

// -------------------------------------------------------------------------
string ExecGenoMS::getType(void) const
{
  return m_type;
}

// -------------------------------------------------------------------------
void ExecGenoMS::setName(std::string name)
{
  m_name = name;
}

// -------------------------------------------------------------------------

//Full validation of the parameters happens within GenoMS code
bool ExecGenoMS::validateParams(std::string & error)
{

  m_isValid = false;

  VALIDATE_PARAM_EXIST("OUTPUT_FILE");
  VALIDATE_PARAM_EXIST("EXE_DIR");

  m_isValid = true;


  return true;
}



ExecBase * ExecGenoMS::clone(const ParameterList & input_params) const
{
  return new ExecGenoMS(input_params);
}

bool ExecGenoMS::loadInputData(void)
{
  return true;
}


bool ExecGenoMS::saveOutputData(void)
{
  return true;
}

bool ExecGenoMS::saveInputData(std::vector<std::string> & filenames)
{
  return true;
}

bool ExecGenoMS::loadOutputData(void)
{
  return true;
}

vector<ExecBase *> const & ExecGenoMS::split(int numSplit)
{
  m_subModules.resize(0);
  return m_subModules;
}

bool ExecGenoMS::merge(void)
{
}



bool ExecGenoMS::invoke(void)
{

  DEBUG_MSG("Entering ExecGenoMS::invoke()");

  std::string exeDir = m_params.getValue("EXE_DIR");
  rtrim(exeDir);

  int ret = callGenoMS(exeDir);
  if(ret != 0)
    {
      ERROR_MSG("Error executing GenoMS!");
      return false;
    }


  return true;

}



int ExecGenoMS::callGenoMS(string & exeDir)
{
#if INCLUDE_GENOMS == 0
  ERROR_MSG("GenoMS is not available. Returning in error.");
  return false;
#endif

  //GenoMS command
  std::string GenoMSCommand;

  //Because I don't know how to get the name of the original config file,
  //We have to write a new one for GenoMS
  //m_params.setValue("MSGFDBPATH",exeDir);
  m_params.writeToFile("genoMS.params");


#if defined(__linux__) || defined(__CYGWIN__)
  GenoMSCommand = "unset LD_LIBRARY_PATH && ";
#endif
  GenoMSCommand += "java -jar -Xmx2G ";

  GenoMSCommand += exeDir;
  GenoMSCommand += "/GenoMS.jar";

  GenoMSCommand += " -i genoMS.params";

  GenoMSCommand += " -o ";
  GenoMSCommand += m_params.getValue("OUTPUT_FILE");

  GenoMSCommand += " -r ";
  GenoMSCommand += exeDir;
  GenoMSCommand += "/DBs_GenoMS";

  if(m_params.exists("PROJECT_DIR"))
    {
      GenoMSCommand += " -k ";
      GenoMSCommand += m_params.getValue("PROJECT_DIR");
    }

  GenoMSCommand += " -e ";
  GenoMSCommand += exeDir;


  if(m_params.exists("RUN_DBSEARCH"))
    GenoMSCommand += " -x";

  if(m_params.exists("PEAK_PENALTY"))
    GenoMSCommand += " -p";

  if(m_params.exists("ALL2ALL_SIMILARITY"))
    GenoMSCommand += " -s";

  if(m_params.exists("HMM_LATE_ADD"))
    GenoMSCommand += " -a";

  if(m_params.exists("FDR_CUTOFF"))
    {
      GenoMSCommand += " -w ";
      GenoMSCommand += m_params.getValue("FDR_CUTOFF");
    }

  if(m_params.exists("MUTATION_MODE"))
    GenoMSCommand += " -f";

  if(m_params.exists("GENERATE_REPORTS"))
    {
      GenoMSCommand += " -q";
    }
  if(m_params.exists("LOG_FILE"))
    {
      GenoMSCommand += " -l ";
      GenoMSCommand += m_params.getValue("LOG_FILE");
    }

  DEBUG_MSG("call genoMS: " << GenoMSCommand);

  int ret = spsSystem(GenoMSCommand.c_str());
  DEBUG_MSG("GenoMS return value: " << ret);
  return ret;

}
