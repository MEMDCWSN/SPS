// Module Includes
#include "Logger.h"
#include "SgeGridMonitor.h"

// External Module Includes
#include "utils.h"

// System Includes
#include <fstream>
#include <stdio.h>

using namespace specnets;
using namespace std;

static string getUsername(void)
{
  int retVal = system("whoami > username.tmp");
  ifstream ifs("username.tmp", ios::binary);
  char buf[256];
  ifs.getline(buf, 256);
  return buf;
}

namespace specnets
{
  // -------------------------------------------------------------------------
  SgeGridMonitor::SgeGridMonitor(void)
  {
    /* EMPTY */
  }

  // -------------------------------------------------------------------------
  SgeGridMonitor::~SgeGridMonitor(void)
  {
    jobs.clear();
  }

  // -------------------------------------------------------------------------
  SgeGridMonitor::JobStatus SgeGridMonitor::getJobStatus(const string & jobId)
  {
#if 0
    // This is the original code that used to work
    //   However now it seems the job ID gets truncated in the info list
    //   so they wont match exactly
    if (jobs.find(jobId) != jobs.end()) {
      return jobs[jobId];
    }
#else
    // So now we are using code that provides for a sub-match
    std::map<std::string, JobStatus>::iterator itrJob = jobs.begin();
    std::map<std::string, JobStatus>::iterator itrJobEnd = jobs.end();

    for (; itrJob != itrJobEnd; itrJob++) {
      if (jobId.find(itrJob->first) != string::npos) {
        return itrJob->second;
      }
    }
#endif
    return SgeGridMonitor::JS_UNKNOWN;
  }

  // -------------------------------------------------------------------------
  SgeGridMonitor::JobStatus SgeGridMonitor::getJobStatus(int jobId)
  {
    return SgeGridMonitor::JS_UNKNOWN;
  }

  // -------------------------------------------------------------------------
  string SgeGridMonitor::submitJob(const string & sgePath,
                                   const string & jobFilename,
                                   const string & params,
                                   const string & gridType)
  {
    string command(sgePath);
    if (!command.empty())
    {
      command += "/";
    }
    command += "qsub ";
    command += params;
    command += " ";
    command += jobFilename;
    command += " > ./qstat.tmp";
    DEBUG_VAR(command);
    int retVal = system(command.c_str());

    ifstream ifs("./qstat.tmp", ios::binary);
    char buf[256];
    int linecount = 0;
    if (ifs.good() && !ifs.eof())
    {
      ifs.getline(buf, 256);
      DEBUG_VAR(buf);

      string jobId;
      if (gridType == "pbs") {
        jobId = buf;
        if (jobId.empty()) {
          ERROR_MSG("Unable to parse qsub output.");
          ERROR_MSG("qstat.tmp is empty?.");
          return "";
        }
        DEBUG_VAR(jobId);
        return jobId;
      } else if (gridType == "sge") {
        list<string> listString;
        splitText(buf, listString, " \t");
        list<string>::iterator itr = listString.begin();
        list<string>::iterator itr_end = listString.end();
        //DEBUG_VAR(listString.size());
        if (listString.size() < 5) {
          ERROR_MSG("Unable to parse qsub output.");
          ERROR_MSG("Not enough [" << listString.size() << "] entries in file.");
          return "";
        }
        for (int count = 0; itr != itr_end; itr++,count++) {
          if (count == 2) {
            jobId = *itr;
            return jobId;
          }
        }
      } else {
        ERROR_MSG("Unknown grid type");
        return "";
      }

    } // if (ifs.good() && !ifs.eof())

    return "";
  }


  // -------------------------------------------------------------------------
  bool SgeGridMonitor::refreshInfo(const string & sgePath)
  {
    //DEBUG_TRACE;
    jobs.clear();

    string command(sgePath);
    if (!command.empty())
    {
      command += "/";
    }
    command += "qstat > ./qstat.tmp";
    //DEBUG_VAR(command);
    int retVal = system(command.c_str());

    ifstream ifs("./qstat.tmp", ios::binary);
    char buf[256];
    int linecount = 0;
    while (ifs.good() && !ifs.eof())
    {
      ifs.getline(buf, 256);
      // Skip the first two lines.. they are headed info
      linecount++;
      if (linecount < 3)
      {
        continue;
      }
      if (strlen(buf) == 0)
      {
        continue;
      }
      //DEBUG_VAR(buf);
      list<string> listString;
      splitText(buf, listString, " \t");
      //DEBUG_VAR(listString.size());

      if (listString.size() < 5)
      {
        ERROR_MSG("Unable to parse qstat output");
        return false;
      }
      list<string>::iterator itr = listString.begin();
      list<string>::iterator itr_end = listString.end();
      string jobId;
      string status;
      for (int count = 0; itr != itr_end; itr++,count++)
      {
        if (count == 0)
        {
          jobId = *itr;
          //DEBUG_VAR(jobId);
        }
        if (count == 4)
        {
          status = *itr;
          //DEBUG_VAR(status);
        }
      }

      SgeGridMonitor::JobStatus js = JS_UNKNOWN;
      if (status == "d")
      {
        js = JS_DELETION;
      }
      else if (status == "E")
      {
        js = JS_ERROR;
      }
      else if (status == "h")
      {
        js = JS_HOLD;
      }
      else if (status == "r" || status == "R" || status == "t" )
      {
        js = JS_RUNNING;
      }
      else if (status == "R")
      {
        js = JS_RESTARTED;
      }
      else if (status == "s" || status == "S"|| status == "T")
      {
        js = JS_SUSPENDED;
      }
      else if (status == "w" || status == "qw")
      {
        js = JS_WAITING;
      }

      jobs[jobId] = js;
    }

    return false;
  }

} // namespace specnets

