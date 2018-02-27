#include "main_specnets_helpers.h"
#include "FileUtils.h"
#include "ParameterList.h"

// System Includes
#include <fstream>
#include <iostream>
#include <string>
#include <unistd.h>
#include <stdlib.h>

using namespace std;

namespace specnets
{
  //-----------------------------------------------------------------------------
  string getCurrentDirectory(void)
  {
    char buf[1024];
    getcwd(buf, 1024);
    string currentWorkindDir(buf);
    return currentWorkindDir;

    //int dummy;
    //dummy = spsSystem("pwd > pwd.tmp");
    //ifstream ifs("pwd.tmp", ios::binary);
    //char buf[1024];
    //ifs.getline(buf, 1024);
    //return buf;
  }

  //-----------------------------------------------------------------------------
  string getProjPath(const ParameterList & pl, 
                     const string & addPath, 
                     bool forceAbsolutePath)
  {
    string projDir = pl.getValue("PROJECT_DIR", "");
    bool relDir = pl.getValueBool("RELATIVE_DIR", true);
    if (forceAbsolutePath) {
      relDir = false;
    }

    return getPath(projDir, addPath, relDir);
  }

  //-----------------------------------------------------------------------------
  string getCurrentTimeString(void)
  {
    time_t rawtime;
    struct tm * timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);

    return asctime(timeinfo);
  }

}

