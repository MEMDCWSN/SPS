#include <stdio.h>
#include <unistd.h>
#include <signal.h>
#include <dirent.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
//#include <sys/mman.h>

#if defined(__linux__)
#include <cxxabi.h>
#include <execinfo.h>
#include <sys/wait.h>
#include <sys/sysinfo.h>
#elif defined(__MINGW32__) || defined(__CYGWIN__)
#include <windows.h>
#endif

#include <ctime>
#include <cctype>
#include <cstdlib>
#include <set>
#include <map>
#include <vector>
#include <limits>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <iterator>
#include <typeinfo>
#include <algorithm>
#include "array.h"
#include "range.h"

#include "mzxml.h"
#include "spsplot.h"
#include "cstream.h"
#include "iomanip.h"
#include "dbg_print.h"
#include "spectrum.h"
#include "inputParams.h"
//#include "utils2.h"

#include "ExecFramework/copyright.h"

#if defined(__MINGW32__)
#define abort() return -9

#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif // defined(__MINGW32__) || defined(__CYGWIN__)

#if !defined(ACCESSPERMS)
#define ACCESSPERMS (S_IRWXU|S_IRWXG|S_IRWXO)
#endif


using namespace std;
using namespace sps;
using namespace specnets;


namespace specnets
{
const char * shtmlhead =
  "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\"\n"
  "  \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n"

  "<HTML xmlns=\"http://www.w3.org/1999/xhtml\">\n"

  "<head>\n"

  "  <meta http-equiv=\"Content-Type\" content=\"text/shtml; charset=ISO-8859-1\" />\n"
//   "  <meta http-equiv=\"PRAGMA\" content=\"NO-CACHE\" />\n"
//   "  <meta http-equiv=\"CACHE-CONTROL\" content=\"NO-STORE, NO-CACHE, MUST-REVALIDATE, POST-CHECK=0, PRE-CHECK=0\" />\n"
//   "  <meta http-equiv=\"EXPIRES\" content=\"01 Jan 1970 00:00:00 GMT\" />\n"
//   "  <meta http-equiv=\"Last-Modified\" content=\"01 Jan 1970 00:00:00 GMT\" />\n"
//   "  <meta http-equiv=\"If-Modified-Since\" content=\"01 Jan 1970 00:00:00 GMT\" />\n"

  "  <title>UCSD Computational Mass Spectrometry Website</title>\n";
  
const char * shtmlhead2 =

  "  <script type=\"JavaScript\">\n"
  "  <!--hide\n"
  "  function modalWin(filename, name)\n"
  "  {\n"
  "    if (! name)\n"
  "      name = \"name\";\n"

  "    if (window.showModalDialog)\n"
  "    {\n"
  "      window.showModalDialog(filename, name, \"dialogWidth:640px; dialogHeight:480px\");\n"
  "    }\n"
  "    else\n"
  "    {\n"
  "      window.open(filename, name, 'width=640, height=480, toolbar=no, directories=no, status=no, menubar=no, scrollbars=no, resizable=no, modal=yes');\n"
  "    }\n"
  "  }\n"

  "  function modalForm(form, name)\n"
  "  {\n"
  "    if (! name)\n"
  "      name = \"name\";\n"

  "    if (window.focus)\n"
  "    {\n"
  "      modalWin(\"\", name);\n"
  "      form.target = name;\n"
  "    }\n"
  "  }\n"
  "  //-->\n"
  "  </script>\n"

  "  <script src=\"../../js/pager.js\" language=\"javascript\" type=\"text/javascript\"></script>\n"
  "  <noscript>\n"
  "    <meta http-equiv=\"refresh\" content=\"10\">\n"
  "  </noscript>\n"

  "  <script language=\"JavaScript\">\n"
  "  <!--\n"
  "  var sURL = unescape(window.location.pathname);\n"
  "  function doLoad()\n"
  "  {\n"
  "      setTimeout( \"refresh()\", 10*1000 );\n"
  "  }\n"

  "  function refresh()\n"
  "  {\n"
  "      window.location.href = sURL;\n"
  "  }\n"
  "  //-->\n"
  "  </script>\n"

  "  <script language=\"JavaScript1.1\">\n"
  "  <!--\n"
  "  function refresh()\n"
  "  {\n"
  "      window.location.replace( sURL );\n"
  "  }\n"
  "  //-->\n"
  "  </script>\n"

  "  <script language=\"JavaScript1.2\">\n"
  "  <!--\n"
  "  function refresh()\n"
  "  {\n"
  "      window.location.reload(true);\n"
  "  }\n"
  "  //-->\n"
  "  </script>\n"

//   "  <script language=\"JavaScript\">\n"
//   "  <!--\n"
//   "  function needReload()\n"
//   "  {\n"
//   "      var loc=window.location.href.toString()\n"
//   "      if(loc.indexOf(\"?\") == -1)\n"
//   "      {\n"
//   "          location.replace(loc+\"?1\");\n"
//   "      }\n"
//   "      else\n"
//   "      {\n"
//   "          return false;\n"
//   "      }\n"
//   "  }\n"
//   "  //-->\n"
//   "  </script>\n"

//   "  <script language=\"JavaScript1.2\">\n"
//   "  <!--\n"
//   "  function refreshScroll()\n"
//   "  {\n"
//   "    sp = document.body.scrollTop;\n"
//   "    refresh();\n"
//   "    window.scrollTo(0, \"+sp+\");\n"
//   "  }\n"
//
//   "  function doLoadScroll()\n"
//   "  {\n"
//   "      setTimeout( \"refreshScroll()\", 10*1000 );\n"
//   "  }\n"
//
//   // function saves scroll position
//   "  function fScroll(val)\n"
//   "  {\n"
//   "      var hidScroll = document.getElementById('hidScroll');\n"
//   "      hidScroll.value = val.scrollTop;\n"
//   "  }\n"
//
//   // function moves scroll position to saved value
//   "  function fScrollMove(what)\n"
//   "  {\n"
//   "      var hidScroll = document.getElementById('hidScroll');\n"
//   "      document.getElementById(what).scrollTop = hidScroll.value;\n"
//   "  }\n"
//   "  //-->\n"
//   "  </script>\n"

  "</head>\n";

#if defined(__MINGW32__)
inline void symlink (const char * sold, const char * snew)
{
  CopyFile(sold, snew, false);
}

inline TCHAR * cuserid(char *)
{
  static TCHAR sbuffer[32767];
  DWORD ncount = sizeof(sbuffer);

  GetUserName(sbuffer, & ncount);

  return sbuffer;
}
#endif



void SpsPlot::genHtmlHeader(ostream & shtml)
{
  shtml << shtmlhead << '\n';
  
  string path = ".";
  if(sopt[htmlDefs].length() > 0)
    path = sopt[htmlDefs];
  
  shtml << "  <link href=\"";
  shtml << path;
  shtml << "/styles/main.css\" rel=\"stylesheet\" type=\"text/css\" />\n";
  shtml << "  <link rel=\"shortcut icon\" href=\"";
  shtml << path;
  shtml << "/images/favicon.ico\" type=\"image/icon\" />\n";
  shtml << "  <script src=\"";
  shtml << path;
  shtml << "/scripts/util.js\" language=\"javascript\" type=\"text/javascript\"></script>\n";
  shtml << "  <script src=\"";
  shtml << path;
  shtml << "/scripts/download.js\" language=\"javascript\" type=\"text/javascript\"></script>\n";
  shtml << "  <script src=\"";
  shtml << path;
  shtml << "/scripts/render.js\" language=\"javascript\" type=\"text/javascript\"></script>\n";
  shtml << "  <script src=\"";
  shtml << path;
  shtml << "/scripts/inspect.js\" language=\"javascript\" type=\"text/javascript\"></script>\n";
  shtml << "  <script src=\"";
  shtml << path;
  shtml << "/js/mootools.js\" type=\"text/javascript\"></script>\n";
  shtml << "  <script src=\"";
  shtml << path;
  shtml << "/js/slimbox.js\" type=\"text/javascript\"></script>\n";
  shtml << "  <link href=\"";
  shtml << path;
  shtml << "/css/slimbox.css\" rel=\"stylesheet\" type=\"text/css\" media=\"screen\" />\n";
  shtml << "  <script src=\"";
  shtml << path;
  shtml << "/js/sorttable.js\"></script>\n";

  shtml << shtmlhead2 << '\n';
}


void SpsPlot::genHtmlHeader2(ostream & shtml)
{
/*  string path = ".";
  if(sopt[htmlDefs].length() > 0)
    path = sopt[htmlDefs];

  shtml << "  <div id=\"rb_logos\">\n";
  shtml << "    <h3 align=center><img src=\"";
  shtml << path;
  shtml << "/images/pagelogo_ccms_left.jpg\" border=\"0\"/></h3>\n";
  shtml << "  </div>\n";*/
  
  string title = "";
  if(sopt[ePageTitle].length() > 0)
    title = sopt[ePageTitle];

  shtml << "  <div id=\"rb_logos\">\n";
  shtml << "    <h3 align=center>";
  shtml << title;
  shtml << "</h3>\n";
  shtml << "  </div>\n";
  
  
}



void SpsPlot::dumpAbruijn(char *fname)
{
  ofstream outFile(fname);
  Data::ccluster_t::iterator it;
  for(it = cdata.ccluster.begin() ; it != cdata.ccluster.end() ; it++) {
    outFile << " ------- " << it->first << " ------- " << endl;
    list< pair<unsigned, unsigned> >::iterator it2;
    for(it2 = it->second.begin() ; it2 != it->second.end() ; it2++) {
      outFile << "\t" << it2->first << " : " << it2->second << endl;
    }
  }
  outFile.close();
}






bool spsplot_sort(const string &a, const string &b)
{
  if(a.length() == b.length())
    return a.compare(b) < 0;
  return (a.length() < b.length());
}

sps::vector<string> entries(const string & str, bool sort)
{
  // get directory and file prefix
  int found = str.find_last_of("/\\");
  string directory = str.substr(0,found);
  string file = str.substr(found+1);

   // return values
  sps::vector<string> centry;

  // open directory
  struct dirent * pe;
  DIR * pd = opendir(directory.c_str());
  // if there is something wrong, return empty results
  if (! pd)
    return sps::vector<string>(centry.begin(), centry.end());

  // cycle thru found files
  for (struct stat cs; pe = readdir(pd); ) {

    // get the file name ( + directory)
    string aux = pe->d_name;
    // find filename start
    int found = aux.find_last_of("/\\");
    // get the file name
    string str = aux.substr(found+1);
  
    // test for directory items    
    if ( ( str.compare(".") == 0 ) || ( str.compare("..") == 0 ) )
      continue;

    // test for prefix and sufix
    if(str.compare(0, file.length(), file) == 0 ) {

      // compose string with path
      string aux2 = directory;
      aux2 += '/';
      aux2 += str;
      // add to list if conditions matched
      centry.push_back(aux2);
    }

  }
  closedir(pd);

  if(sort)
    std::sort(centry.begin(), centry.end(), spsplot_sort);

  return sps::vector<string>(centry.begin(), centry.end());}

/*
sps::vector<string> entries(const string & spattern, bool sort)
{
  std::vector<string> centry;
  const size_t nslash = spattern.rfind('/');
  size_t nstar[2];
  nstar[0] = spattern.find_first_of("*", nslash);
  nstar[1] = spattern.find_first_of("*", nstar[0]);

  struct dirent * pe;
  DIR * pd = opendir(spattern.substr(0, nslash).c_str());

  if (! pd)
    return sps::vector<string>(centry.begin(), centry.end());

  for (struct stat cs; pe = readdir(pd); )
  {
    if ( ( strcmp(".", pe->d_name) == 0 ) || ( strcmp("..", pe->d_name) == 0 ) )
      continue;

//      if (spattern.substr(nslash + 1, nstar[0] - nslash - 1) + ssni + spattern.substr(nstar[1] + 1) == pe->d_name)

    string fileName = wildcardToRegex(spattern.substr(nslash + 1, - 1) );

    if(regularExpressionMatch(fileName, pe->d_name) == 0)
      centry.push_back(spattern.substr(0, nslash) + "/" + pe->d_name);
  }
  closedir(pd);

  if(sort)
    std::sort(centry.begin(), centry.end(), spsplot_sort);
  
  return sps::vector<string>(centry.begin(), centry.end());
} */


/*
sps::vector<string> entries(const string & spattern)
{
  list<string> centry;
  const size_t nslash = spattern.rfind('/');
  size_t nstar[2];
  nstar[0] = spattern.find_first_of("*", nslash);
  nstar[1] = spattern.find_first_of("*", nstar[0]);

  for (int ni = 1; true; ++ ni)
  {
    struct dirent * pe;
    DIR * pd = opendir(spattern.substr(0, nslash).c_str());

    if (! pd)
      return sps::vector<string>(centry.begin(), centry.end());

    ostringstream sni;
    sni << ni;
    const string ssni = sni.str();

    for (struct stat cs; pe = readdir(pd); )
    {
      if ( ( strcmp(".", pe->d_name) == 0 ) || ( strcmp("..", pe->d_name) == 0 ) )
        continue;

      if (spattern.substr(nslash + 1, nstar[0] - nslash - 1) + ssni + spattern.substr(nstar[1] + 1) == pe->d_name)
      {
        centry.push_back(spattern.substr(0, nslash) + "/" + pe->d_name);
        break;
      }
    }
    closedir(pd);

    if (! pe)
      break;
  }

  return sps::vector<string>(centry.begin(), centry.end());
} */
}

#if defined(__MINGW32__) || defined(__CYGWIN__)
int (* gnuplot_main)(int, char **);
#else
extern "C" int gnuplot_main(int argc, char **argv);
#endif

int main(int argc, char **argv)
{
#if defined(__MINGW32__) || defined(__CYGWIN__)
  HINSTANCE pdll = LoadLibraryA("gnuplot.dll");

  if (! pdll)
  {
    cerr << "error: cannot load 'gnuplot.dll': " << GetLastError() << endl;
    return -3;
  }

  gnuplot_main = (int (*)(int, char **)) (GetProcAddress(pdll , "gnuplot_main"));

  if (! gnuplot_main)
  {
    cerr << "error: cannot find 'gnuplot_main'" << endl;
    return -3;
  }
#endif

  SpsPlot m;

  // path
  if (argv[0][0] == '/')
    m.spath = argv[0];
  else
  {
    char sCwd[FILENAME_MAX];

    if (! GetCurrentDir(sCwd, sizeof(sCwd)))
    {
      cerr << "error: cannot get current directory: " << strerror(errno) << endl;
      return -2;
    }

    m.spath = string(sCwd) + "/" + argv[0];
  }

  // process options
  for (int i = 1; i < argc; ++ i)
    if (! strcmp(argv[i], "--help"))
      return m.help();
    else if (! strcmp(argv[i], "--version"))
      return m.version();
    else if (! strcmp(argv[i], m.sbopt[m.eannot]))
      m.bopt[m.eannot] = true;
    else if (! strcmp(argv[i], m.sbopt[m.everbose]))
      m.bopt[m.everbose] = true;
    else if (! strcmp(argv[i], m.sbopt[m.eqtok]))
      m.bopt[m.eqtok] = true;
    else if (! strcmp(argv[i], m.sbopt[m.ektoq]))
      m.bopt[m.ektoq] = true;
    else if (! strcmp(argv[i], m.sbopt[m.eshift]))
      m.bopt[m.eshift] = true;
    else if (! strcmp(argv[i], m.sbopt[m.estats]))
      m.bopt[m.estats] = true;
    else if (! strcmp(argv[i], m.sbopt[m.enotitle]))
      m.bopt[m.enotitle] = true;
    // double
    else if (++ i < argc)
      if (! strcmp(argv[i - 1], m.ssopt[m.eproject]))
        m.sopt[m.eproject] = argv[i];
      else if (! strcmp(argv[i - 1], m.ssopt[m.eaa]))
        m.sopt[m.eaa] = argv[i];
      else if (! strcmp(argv[i - 1], m.ssopt[m.epklbin]))
        m.sopt[m.epklbin] = argv[i];
      else if (! strcmp(argv[i - 1], m.ssopt[m.emgf]))
        m.sopt[m.emgf] = argv[i];
      else if (! strcmp(argv[i - 1], m.ssopt[m.epkl]))
        m.sopt[m.epkl] = argv[i];
      else if (! strcmp(argv[i - 1], m.ssopt[m.emzxml]))
        m.sopt[m.emzxml] = argv[i];
#if ! defined(SPECPLOT)
      else if (! strcmp(argv[i - 1], m.ssopt[m.estars]))
        m.sopt[m.estars] = argv[i];
      else if (! strcmp(argv[i - 1], m.ssopt[m.ecomp]))
        m.sopt[m.ecomp] = argv[i];
      else if (! strcmp(argv[i - 1], m.ssopt[m.eseqs]))
        m.sopt[m.eseqs] = argv[i];
      else if (! strcmp(argv[i - 1], m.ssopt[m.eindex]))
        m.sopt[m.eindex] = argv[i];
      else if (! strcmp(argv[i - 1], m.scopt[m.ecluster]))
      {
        list<string> s;

        for (; i < argc; ++ i)
          if (argv[i][0] == '-')
            { i --; break; }
          else if (strlen(argv[i]))
            s.push_back(argv[i]);

        m.copt[m.ecluster].assign(s.begin(), s.end());
      }
      else if (! strcmp(argv[i - 1], m.scopt[m.eclusterms]))
      {
        list<string> s;

        for (; i < argc; ++ i)
          if (argv[i][0] == '-')
            { i --; break; }
          else if (strlen(argv[i]))
            s.push_back(argv[i]);

        m.copt[m.eclusterms].assign(s.begin(), s.end());
      }
      else if (! strcmp(argv[i - 1], m.scopt[m.eclusterscan]))
      {
        list<string> s;

        for (; i < argc; ++ i)
          if (argv[i][0] == '-')
            { i --; break; }
          else if (strlen(argv[i]))
            s.push_back(argv[i]);

        m.copt[m.eclusterscan].assign(s.begin(), s.end());
      }
      else if (! strcmp(argv[i - 1], m.ssopt[m.emp]))
        m.sopt[m.emp] = argv[i];
      else if (! strcmp(argv[i - 1], m.ssopt[m.emidx]))
        m.sopt[m.emidx] = argv[i];
      else if (! strcmp(argv[i - 1], m.ssopt[m.efasta]))
        m.sopt[m.efasta] = argv[i];
      else if (! strcmp(argv[i - 1], m.ssopt[m.erefindex]))
        m.sopt[m.erefindex] = argv[i];
      else if (! strcmp(argv[i - 1], m.ssopt[m.erefmp]))
        m.sopt[m.erefmp] = argv[i];
      else if (! strcmp(argv[i - 1], m.ssopt[m.erefmidx]))
        m.sopt[m.erefmidx] = argv[i];
#else
      else if (! strcmp(argv[i - 1], m.ssopt[m.elabel]))
        m.sopt[m.elabel] = argv[i];
#endif // SPECPLOT
      else if (! strcmp(argv[i - 1], m.ssopt[m.efont]))
        m.sopt[m.efont] = argv[i];
      else if (! strcmp(argv[i - 1], m.ssopt[m.eoutdir]))
        m.sopt[m.eoutdir] = argv[i];
      else if (! strcmp(argv[i - 1], m.ssopt[m.eoutfile]))
        m.sopt[m.eoutfile] = argv[i];
      else if (! strcmp(argv[i - 1], m.ssopt[m.econtig]))
        m.sopt[m.econtig] = argv[i];
      else if (! strcmp(argv[i - 1], m.ssopt[m.enodes]))
        m.sopt[m.enodes] = argv[i];
      else if (! strcmp(argv[i - 1], m.ssopt[m.ep]))
        m.sopt[m.ep] = argv[i];
      else if (! strcmp(argv[i - 1], m.ssopt[m.espectrum]))
        m.sopt[m.espectrum] = argv[i];
      else if (! strcmp(argv[i - 1], m.ssopt[m.espectrumscan]))
        m.sopt[m.espectrumscan] = argv[i];
      else if (! strcmp(argv[i - 1], m.ssopt[m.espectruminfo]))
        m.sopt[m.espectruminfo] = argv[i];
      else if (! strcmp(argv[i - 1], m.scopt[m.epeptide]))
      {
        list<string> s;

        for (; i < argc; ++ i)
          if (argv[i][0] == '-')
            { i --; break; }
          else if (strlen(argv[i]))
          {
            s.push_back(argv[i]);

            // constant
            for (string::size_type np; np = s.back().find(Ion::cname[Ion::phos]), np != string::npos; )
              s.back().replace(np, Ion::cname[Ion::phos].size(), "+80");

            transform(s.back().begin(), s.back().end(), s.back().begin(), ::toupper);
          }

        m.copt[m.epeptide].assign(s.begin(), s.end());
      }
      else if (! strcmp(argv[i - 1], m.ssopt[m.eformat]))
      {
        m.sopt[m.eformat] = argv[i];
        transform(m.sopt[m.eformat].begin(), m.sopt[m.eformat].end(), m.sopt[m.eformat].begin(), ::tolower);
      }
      else if (! strcmp(argv[i - 1], m.ssopt[m.ezoom]))
        m.sopt[m.ezoom] = argv[i];
      else if (! strcmp(argv[i - 1], m.ssopt[m.eprefix]))
        m.sopt[m.eprefix] = argv[i];
      else if (! strcmp(argv[i - 1], m.ssopt[m.epeakmasstol]))
        m.sopt[m.epeakmasstol] = argv[i];
      else if (! strcmp(argv[i - 1], m.ssopt[m.eparentmasstol]))
        m.sopt[m.eparentmasstol] = argv[i];
      else if (! strcmp(argv[i - 1], m.ssopt[m.einspect]))
        m.sopt[m.einspect] = argv[i];
      else if (! strcmp(argv[i - 1], m.ssopt[m.erange]))
        m.sopt[m.erange] = argv[i];
      else if (! strcmp(argv[i - 1], m.ssopt[m.htmlDefs]))
        m.sopt[m.htmlDefs] = argv[i];
      else if (! strcmp(argv[i - 1], m.ssopt[m.ePageTitle]))
        m.sopt[m.ePageTitle] = argv[i];
      else
        return m.error(argv[i - 1]);
    else
      return m.error(argv[i - 1]);

  // traditional parameter file format
  if (! m.sopt[m.ep].empty())
    for (unsigned i = 0; i < 2; ++ i)
    {
      InputParams cparams;

      switch (i)
      {
      case 0:
        continue;

      case 1:
        cparams.readParams(const_cast<char *>(m.sopt[m.ep].c_str()));
        break;
      }

      list<string> clargv;

      if (cparams.paramPresent("PROJECT"))
        m.sopt[m.eproject] = cparams.getValue("PROJECT");
      if (cparams.paramPresent("FILE_AA"))
        m.sopt[m.eaa] = cparams.getValue("FILE_AA");
      if (cparams.paramPresent("FILE_STARS"))
        m.sopt[m.estars] = cparams.getValue("FILE_STARS");
      if (cparams.paramPresent("FILE_COMP"))
        m.sopt[m.ecomp] = cparams.getValue("FILE_COMP");
      if (cparams.paramPresent("FILE_SEQS"))
        m.sopt[m.eseqs] = cparams.getValue("FILE_SEQS");
      if (cparams.paramPresent("FILE_MS"))
        m.sopt[m.epklbin] = cparams.getValue("FILE_MS");
      if (cparams.paramPresent("FILE_INDEX"))
        m.sopt[m.eindex] = cparams.getValue("FILE_INDEX");
      if (cparams.paramPresent("FILE_CLUSTER"))
//        readFilesFromFile(cparams.getValue("FILE_CLUSTER"), m.copt[m.ecluster], true);
        m.copt[m.ecluster] = entries(cparams.getValue("FILE_CLUSTER"), true);
      if (cparams.paramPresent("FILE_CLUSTERMS"))
        readFilesFromFile(cparams.getValue("FILE_CLUSTERMS"), m.copt[m.eclusterms], true);
//        m.copt[m.eclusterms] = entries(cparams.getValue("FILE_CLUSTERMS"), true);
      if (cparams.paramPresent("FILE_CLUSTERSCAN"))
        readFilesFromFile(cparams.getValue("FILE_CLUSTERSCAN"), m.copt[m.eclusterscan], true);
//        m.copt[m.eclusterscan] = entries(cparams.getValue("FILE_CLUSTERSCAN"), true);
      if (cparams.paramPresent("FILE_MP"))
        m.sopt[m.emp] = cparams.getValue("FILE_MP");
      if (cparams.paramPresent("FILE_MIDX"))
        m.sopt[m.emidx] = cparams.getValue("FILE_MIDX");
      if (cparams.paramPresent("FILE_FASTA"))
        m.sopt[m.efasta] = cparams.getValue("FILE_FASTA");
      if (cparams.paramPresent("FILE_REFINDEX"))
        m.sopt[m.erefindex] = cparams.getValue("FILE_REFINDEX");
      if (cparams.paramPresent("FILE_REFMP"))
        m.sopt[m.erefmp] = cparams.getValue("FILE_REFMP");
      if (cparams.paramPresent("FILE_REFMIDX"))
        m.sopt[m.erefmidx] = cparams.getValue("FILE_REFMIDX");
      if (cparams.paramPresent("OUTDIR"))
        m.sopt[m.eoutdir] = cparams.getValue("OUTDIR");
      if (cparams.paramPresent("GRID_NUMNODES"))
        m.sopt[m.enodes] = cparams.getValue("GRID_NUMNODES");
      if (cparams.paramPresent("TOLERANCE_PEAK"))
        m.sopt[m.epeakmasstol] = cparams.getValue("TOLERANCE_PEAK");
      if (cparams.paramPresent("TOLERANCE_PM"))
        m.sopt[m.eparentmasstol] = cparams.getValue("TOLERANCE_PM");
      if (cparams.paramPresent("PARTIAL_OVERLAPS"))
        m.sopt[m.epartialoverlaps] = cparams.getValue("PARTIAL_OVERLAPS");
      if (cparams.paramPresent("HTML_DEFS"))
        m.sopt[m.htmlDefs] = cparams.getValue("HTML_DEFS");
      if (cparams.paramPresent("REPORT_TITLE"))
        m.sopt[m.ePageTitle] = cparams.getValue("REPORT_TITLE");
      if (cparams.paramPresent("FONT_PATH"))
        m.sopt[m.efont] = cparams.getValue("FONT_PATH");

      if (cparams.paramPresent("FILE_INSPECT"))
        m.sopt[m.einspect] = cparams.getValue("FILE_INSPECT");
      if (cparams.paramPresent("FILE_PRM"))
        m.sopt[m.eprm] = cparams.getValue("FILE_PRM");
      if (cparams.paramPresent("FILE_STARSINDEX"))
        m.sopt[m.estarsindex] = cparams.getValue("FILE_STARSINDEX");
      if (cparams.paramPresent("FILE_PAIRS"))
        m.sopt[m.epairs] = cparams.getValue("FILE_PAIRS");
      if (cparams.paramPresent("FILE_MATCHPAIRS"))
        m.sopt[m.ematchpairs] = cparams.getValue("FILE_MATCHPAIRS");
      if (cparams.paramPresent("FILE_STARSPAIRS"))
        m.sopt[m.estarspairs] = cparams.getValue("FILE_STARSPAIRS");
    }

  if (int nz = m())
  {
    const string sfilename = m.sopt[m.eoutdir] + "/spsplot.log";

    cout << "error: see '" << sfilename << "'" << endl;
    note(cerr) << nz << endl;
    return -2;
  }

  return 0;
}


SpsPlot::ccron_t SpsPlot::ccron;


namespace
{

pair<SpsPlot::bopt_t, const char *> sbopt_value[] =
{
  pair<SpsPlot::bopt_t, const char *>(SpsPlot::eannot,        "--annot"),
  pair<SpsPlot::bopt_t, const char *>(SpsPlot::everbose,      "--verbose"),
  pair<SpsPlot::bopt_t, const char *>(SpsPlot::eqtok,         "--qtok"),
  pair<SpsPlot::bopt_t, const char *>(SpsPlot::ektoq,         "--ktoq"),
  pair<SpsPlot::bopt_t, const char *>(SpsPlot::eshift,        "--shift"),
  pair<SpsPlot::bopt_t, const char *>(SpsPlot::estats,        "--stats"),
  pair<SpsPlot::bopt_t, const char *>(SpsPlot::enotitle,      "--notitle"),
};

pair<SpsPlot::sopt_t, const char *> ssopt_value[] =
{
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::eproject,      "--project"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::eaa,           "--aa"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::estars,        "--stars"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::ecomp,         "--comp"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::eseqs,         "--seqs"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::epklbin,       "--pklbin"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::eindex,        "--index"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::emp,           "--mp"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::emidx,         "--midx"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::efasta,        "--fasta"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::erefindex,     "--refindex"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::erefmp,        "--refmp"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::erefmidx,      "--refmidx"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::efont,         "--font"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::eoutdir,       "--outdir"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::eoutfile,      "--outfile"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::econtig,       "--contig"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::enodes,        "--nodes"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::ep,            "--p"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::espectrum,     "--spectrum"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::espectrumscan, "--spectrumscan"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::eformat,       "--format"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::ezoom,         "--zoom"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::emgf,          "--mgf"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::epkl,          "--pkl"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::emzxml,        "--mzxml"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::espectruminfo, "--spectruminfo"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::eprefix,       "--prefix"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::epeakmasstol,  "--peakmasstol"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::eparentmasstol,"--parentmasstol"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::elabel,        "--label"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::einspect,      "--inspect"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::erange,        "--range"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::htmlDefs,      "--htmlDefs"),
  pair<SpsPlot::sopt_t, const char *>(SpsPlot::ePageTitle,    "--pageTitle")  
};

pair<SpsPlot::copt_t, const char *> scopt_value[] =
{
  pair<SpsPlot::copt_t, const char *>(SpsPlot::ecluster,      "--cluster"),
  pair<SpsPlot::copt_t, const char *>(SpsPlot::eclusterms,    "--clusterms"),
  pair<SpsPlot::copt_t, const char *>(SpsPlot::eclusterscan,  "--clusterscan"),
  pair<SpsPlot::copt_t, const char *>(SpsPlot::epeptide,      "--peptide"),
};

} // namespace


std::map<SpsPlot::bopt_t, const char *> SpsPlot::sbopt(sbopt_value, sbopt_value + sizeof(sbopt_value) / sizeof(* sbopt_value));
std::map<SpsPlot::sopt_t, const char *> SpsPlot::ssopt(ssopt_value, ssopt_value + sizeof(ssopt_value) / sizeof(* ssopt_value));
std::map<SpsPlot::copt_t, const char *> SpsPlot::scopt(scopt_value, scopt_value + sizeof(scopt_value) / sizeof(* scopt_value));

extern void handle_signal_abort(int sig);

int SpsPlot::operator () ()
{
  // integrity check
  if (sopt[epklbin].empty() && sopt[emgf].empty() && sopt[epkl].empty() && sopt[emzxml].empty())
    return error(ssopt[epklbin]);

  if (sopt[eoutfile].find('/') != string::npos)
    return error(ssopt[eoutfile]);

#if ! defined(SPECPLOT)
  if (! bopt[estats])
  {
    if (sopt[ecomp].empty())
      return error(ssopt[ecomp]);
    if (sopt[eseqs].empty())
      return error(ssopt[eseqs]);
  }
#endif

#if defined(__linux__)
  // less restrictive file permissions
  umask(S_IWOTH);
#endif

  if (! sopt[eoutdir].empty())
    for (string::size_type i = 0, j = min(sopt[eoutdir].find('/', i), sopt[eoutdir].size()); i < sopt[eoutdir].size(); i = j + 1, j = min(sopt[eoutdir].find('/', i), sopt[eoutdir].size()))
      if (! sopt[eoutdir].substr(0, j).empty())
#if defined(__linux__) || defined(__CYGWIN__)
        if (mkdir(sopt[eoutdir].substr(0, j).c_str(), ACCESSPERMS) == -1)
#else
        if (mkdir(sopt[eoutdir].substr(0, j).c_str()) == -1)
#endif
          if (errno != EEXIST)
          {
            cerr << "error: cannot create directory '" << sopt[eoutdir].substr(0, j) << "': " << strerror(errno) << endl;
            return -2;
          }

  // stderr pipe
  const string sfilename = sopt[eoutdir] + "/spsplot.log";

  cstream sout(::fopen(sfilename.c_str(), "w"));

  if (! sout.file_pointer())
  {
    cerr << "error: cannot write to '" << sfilename << "': " << strerror(errno) << endl;
    return -2;
  }

  if (int nz = sout.dup2(stderr))
  {
    cout << "error: see '" << sfilename << "'" << endl;
    note(cerr) << nz << endl;
    return -2;
  }

  // init
  ios_base::sync_with_stdio(false);

  signal(SIGFPE, sigimp);
  signal(SIGILL, sigimp);
  signal(SIGABRT, sigimp);
  signal(SIGSEGV, sigimp);
  signal(SIGINT, sigext);
  signal(SIGTERM, sigext);


#if defined(__linux__)
  signal(SIGALRM, & sigexp);
#endif

  ntime[0] = time(0);
  string stime = ctime(& ntime[0]);

  cerr << endl;
  cerr << "*** " << stime.erase(stime.size() - 1) << " ***" << endl;

  // process options
  double npeakmasstol;
  istringstream(sopt[epeakmasstol]) >> npeakmasstol;

  double nparentmasstol;
  istringstream(sopt[eparentmasstol]) >> nparentmasstol;

  // load files
  if (! sopt[ecomp].empty())
    switch (Load_abinfo(sopt[ecomp].c_str(), cdata.cabinfo))
    {
    case -1:
      cerr << "error: cannot open '" << sopt[ecomp] << "': " << strerror(errno) << endl;
      return -2;

    case -2:
      cerr << "error: cannot read '" << sopt[ecomp] << "'" << endl;
      return -2;
    }

  if (! sopt[eseqs].empty())
    if (! cdata.cconsset.LoadSpecSet_pklbin(sopt[eseqs].c_str()))
      return -2;

  if (sopt[enodes].empty())
  {
    if (! sopt[estars].empty())
      if (! cdata.cspecset.LoadSpecSet_pklbin(sopt[estars].c_str()))
        return -2;

    if (! sopt[epklbin].empty())
      if (! cdata.cmsset.LoadSpecSet_pklbin(sopt[epklbin].c_str()))
        return -2;

    if (! sopt[emgf].empty())
      if (! cdata.cmsset.LoadSpecSet_mgf(sopt[emgf].c_str()))
        return -2;

    if (! sopt[epkl].empty())
      if (! cdata.cmsset.LoadSpecSet_pkl(sopt[epkl].c_str()))
        return -2;

    if (! sopt[emzxml].empty())
    {
      sps::vector<short> msLevel;

      if (sopt[espectrumscan].empty())
      {
        if (! LoadMzxml(sopt[emzxml].c_str(), cdata.cmsset, & msLevel, 2))
          return -2;
      }
      else
      {
        if (! LoadMzxml(sopt[emzxml].c_str(), cdata.cmsset, sopt[espectrumscan], & msLevel, 2))
          return -2;
      }
    }

    cdata.cclusterscan.resize(copt[eclusterscan].size());

    for (int i = 0; i < copt[eclusterscan].size(); ++ i)
      if (! Load_binArray(copt[eclusterscan][i].c_str(), cdata.cclusterscan[i]))
        return -2;

    cdata.cspecmap.resize(copt[eclusterms].size());

    for (int i = 0; i < copt[eclusterms].size(); ++ i) {
      if (! cdata.cspecmap[i].LoadSpecSet_pklbin(copt[eclusterms][i].c_str()))
        return -2;
    }

    if (! sopt[emp].empty())
      if (! Load_binArray(sopt[emp].c_str(), cdata.cproteinidx))
        return -2;

    if (! sopt[emidx].empty())
      if (! cdata.cprotset.LoadSpecSet_pklbin(sopt[emidx].c_str()))
        return -2;

    if (! sopt[efasta].empty())
      if (! cdata.cfasta.Load(sopt[efasta].c_str()))
      {
        cerr << "error: cannot open '" << sopt[efasta] << "': " << strerror(errno) << endl;
        return -2;
      }

    ifstream sclusteridx;

    if (! sopt[eindex].empty())
      if (sclusteridx.open(sopt[eindex].c_str()), ! sclusteridx)
      {
        cerr << "error: cannot open '" << sopt[eindex] << "': " << strerror(errno) << endl;
        return -2;
      }

    if (! sopt[eindex].empty())
      for (unsigned i = 0; sclusteridx.peek() != EOF; ++ i)
      {
        string s;

        getline(sclusteridx, s);

        const string::size_type ni = s.rfind('/');

        if (ni != string::npos)
          s = s.substr(ni + 1);

        cdata.cclusteridx[i] = s;
      }

    ifstream scluster[copt[ecluster].size()];

    for (unsigned i = 0; i < copt[ecluster].size(); ++ i)
      if (scluster[i].open(copt[ecluster][i].c_str()), ! scluster[i])
      {
        cerr << "error: cannot open '" << copt[ecluster][i] << "': " << strerror(errno) << endl;
        return -2;
      }

    unsigned nspectrum = 0;

    for (unsigned i = 0; i < copt[ecluster].size(); ++ i) {
      for (unsigned j = 0; scluster[i].peek() != EOF; ++ j)
      {
        string sline, sname;
        char cplot;
        unsigned nclusteridx, nspectrum3;

        getline(scluster[i], sname, '.');
        scluster[i] >> nclusteridx >> cplot >> nspectrum3;
        getline(scluster[i], sname);

        Data::ccluster_t::iterator ic = cdata.ccluster.insert(cdata.ccluster.end(), make_pair(nspectrum, list< pair<unsigned, unsigned> >()));

        unsigned k = 0;

        for (; scluster[i] && scluster[i].peek() != '\n'; ++ k)
        {
          unsigned nclusteridx2, nspectrum2, aux;

          scluster[i] >> aux >> nclusteridx2 >> nspectrum2;
          getline(scluster[i], sline);

          ic->second.push_back(make_pair(nclusteridx2, nspectrum2));
        }
        if ( k == 0 ) {
          cdata.ccluster.erase(ic);
        } //else
        nspectrum++;
      }
      nspectrum--;
    }

//dumpAbruijn("Abruijn.txt");


    if (! sopt[erefmp].empty())
      if (! Load_binArray(sopt[erefmp].c_str(), cdata.creferenceidx))
        return -2;

    if (! sopt[erefmidx].empty())
      if (! cdata.crefset.LoadSpecSet_pklbin(sopt[erefmidx].c_str()))
        return -2;

    ifstream srefindex;

    if (! sopt[erefindex].empty())
      if (srefindex.open(sopt[erefindex].c_str()), ! srefindex)
      {
        cerr << "error: cannot open '" << sopt[erefindex] << "': " << strerror(errno) << endl;
        return -2;
      }

    if (! sopt[erefindex].empty())
      for (unsigned i = 0; srefindex.peek() != EOF; ++ i)
      {
        string s;
        unsigned nc;

        getline(srefindex, s, ':');
        srefindex >> nc;

        if (s == sopt[eproject])
          cdata.crefindex[nc - 1] = i;

        getline(srefindex, s);
      }

    AAJumps caa(1, 0.01, -1, AAJumps::NO_MODS, false);

    if (! sopt[eaa].empty())
      if (! caa.loadJumps(sopt[eaa].c_str()))
      {
        cerr << "error: cannot open '" << sopt[eaa] << "': " << strerror(errno) << endl;
        return -2;
      }

    for (unsigned i = 0; i < caa.masses.size(); ++ i)
      cdata.caa[uf(caa.masses[i], npeakmasstol)] = caa.aaLetters[i];

    for (unsigned i = 0; i < caa.masses.size(); ++ i)
      cdata.ciaa[caa.aaLetters[i]] = caa.masses[i];

    ifstream slabel;

    if (! sopt[elabel].empty())
    {
      if (slabel.open(sopt[elabel].c_str()), ! slabel)
      {
        cerr << "error: cannot open '" << sopt[elabel] << "': " << strerror(errno) << endl;
        return -2;
      }

      map< double, map<double, string> > ccolor;

      for (unsigned i = 0; slabel.peek() != EOF; ++ i)
      {
        float nf[2];
        unsigned nu[2];

        string sline;
        getline(slabel, sline);
        istringstream ssline(sline);

        string sfield[7];

        for (unsigned j = 0; j < sizeof(sfield) / sizeof(* sfield); ++ j)
          getline(ssline, sfield[j], '\t');

        if (sfield[6].empty())
          sfield[6] = "black";

        if (sfield[0].find("BEGIN IONS") != string::npos)
          continue;
        else if (sfield[0].find("END IONS") != string::npos)
          break;

        istringstream(sfield[0]) >> nf[0];
        istringstream(sfield[1]) >> nf[1];
        istringstream(sfield[4]) >> nu[0];
        istringstream(sfield[5]) >> nu[1];

        cdata.clabel[uf(nf[0], npeakmasstol)][uf(nf[1], npeakmasstol)] = sps::make_tuple(sfield[2], sfield[3], nu[0], nu[1], sfield[6]);

        ccolor[nf[0]][nf[1]] = sfield[6];
      }

      // depth
      for (map< double, map<double, string> >::iterator i = ccolor.begin(); i != ccolor.end(); ++ i)
        for (map<double, string>::reverse_iterator j = i->second.rbegin(); j != i->second.rend() && &* j != &* i->second.begin(); ++ j)
          if (find(cdata.ccolor.begin(), cdata.ccolor.end(), j->second) == cdata.ccolor.end())
            cdata.ccolor.push_back(j->second);

      // breadth
      for (map< double, map<double, string> >::iterator i = ccolor.begin(); i != ccolor.end(); ++ i)
        if (! i->second.empty())
          if (find(cdata.ccolor.begin(), cdata.ccolor.end(), i->second.begin()->second) == cdata.ccolor.end())
            cdata.ccolor.push_back(i->second.begin()->second);
    }

    if (! sopt[eprm].empty())
      if (! cdata.cprm.LoadSpecSet_pklbin(sopt[eprm].c_str()))
        return -2;

    if (! sopt[estarsindex].empty())
      if (! Load_binArray(sopt[estarsindex].c_str(), cdata.cstarsindex))
        return -2;

    if (! sopt[ematchpairs].empty())
      if (! cdata.cmatchpairs.LoadSpecSet_pklbin(sopt[ematchpairs].c_str()))
        return -2;

    if (sopt[epartialoverlaps] == "0")
    {
      sps::vector<Results_ASP> cpairs, cstarspairs;

      if (! sopt[epairs].empty())
        if (! Load_results_bin(sopt[epairs].c_str(), cpairs))
          return -2;

      if (! sopt[estarspairs].empty())
        if (! Load_results_bin(sopt[estarspairs].c_str(), cstarspairs))
          return -2;

      cdata.cpairs.resize(cpairs.size());
      cdata.cstarspairs.resize(cstarspairs.size());

      copy(cpairs.begin(), cpairs.end(), cdata.cpairs.begin());
      copy(cstarspairs.begin(), cstarspairs.end(), cdata.cstarspairs.begin());
    }
    else
    {
      if (! sopt[epairs].empty())
        if (! Load_results_bin(sopt[epairs].c_str(), cdata.cpairs))
          return -2;

      if (! sopt[estarspairs].empty())
        if (! Load_results_bin(sopt[estarspairs].c_str(), cdata.cstarspairs))
          return -2;
    }

    ifstream sinspect;

    if (! sopt[einspect].empty())
    {
      if (sinspect.open(sopt[einspect].c_str()), ! sinspect)
      {
        cerr << "error: cannot open '" << sopt[einspect] << "': " << strerror(errno) << endl;
        return -2;
      }

      // skip header
      string s;
      getline(sinspect, s);

      // iterate lines
      for (unsigned i = 0; sinspect.peek() != EOF; ++ i)
      {
        unsigned nu[2];

        string sline;
        getline(sinspect, sline);
        istringstream ssline(sline);

        string sfield[5];

        for (unsigned j = 0; j < sizeof(sfield) / sizeof(* sfield); ++ j)
          getline(ssline, sfield[j], '\t');

        istringstream(sfield[1]) >> nu[0];
        istringstream(sfield[4]) >> nu[1];

        string sacid[2];

        if (sfield[2].size() >= 4)
          sacid[0] = sfield[2].substr(2, sfield[2].size() - 4);

        // alter acids
        replace(sacid[0].begin(), sacid[0].end(), 'I', 'L');

        if (bopt[eqtok])
          replace(sacid[0].begin(), sacid[0].end(), 'Q', 'K');
        else if (bopt[ektoq])
          replace(sacid[0].begin(), sacid[0].end(), 'K', 'Q');

        // strip out ptms
        sacid[1] = sacid[0];
        sacid[1].erase(remove_if(sacid[1].begin(), sacid[1].end(), ::ispunct), sacid[1].end());
        sacid[1].erase(remove_if(sacid[1].begin(), sacid[1].end(), ::isdigit), sacid[1].end());

        cdata.cinspect.insert(make_pair(nu[0], sps::make_tuple(sacid[0], sacid[1], nu[1])));
      }

      // iterate proteins
      for (unsigned i = 0; i < cdata.cfasta.size(); ++ i)
      {
        // alter acids
        replace(cdata.cfasta[i], cdata.cfasta[i] + strlen(cdata.cfasta[i]), 'I', 'L');

        if (bopt[eqtok])
          replace(cdata.cfasta[i], cdata.cfasta[i] + strlen(cdata.cfasta[i]), 'Q', 'K');
        else if (bopt[ektoq])
          replace(cdata.cfasta[i], cdata.cfasta[i] + strlen(cdata.cfasta[i]), 'K', 'Q');
      }
    }
  }

  // filename debug info
  cdata.cspecset.specs.sfilename = sopt[estars];
  cdata.cconsset.specs.sfilename = sopt[eseqs];
  cdata.cprotset.specs.sfilename = sopt[emidx];
  cdata.crefset.specs.sfilename = sopt[erefmidx];
  cdata.cmsset.specs.sfilename = sopt[epklbin];

  for (int i = 0; i < copt[eclusterscan].size(); ++ i)
    cdata.cclusterscan[i].sfilename = copt[eclusterscan][i];

  for (int i = 0; i < copt[eclusterms].size(); ++ i)
    cdata.cspecmap[i].specs.sfilename = copt[eclusterms][i];

  cdata.cproteinidx.sfilename = sopt[emp];
  cdata.creferenceidx.sfilename = sopt[erefmp];

  cdata.cprm.specs.sfilename = sopt[eprm];
  cdata.cmatchpairs.specs.sfilename = sopt[ematchpairs];
  cdata.cstarsindex.sfilename = sopt[estarsindex];
  cdata.cpairs.sfilename = sopt[epairs];
  cdata.cstarspairs.sfilename = sopt[estarspairs];

#if defined(SPECPLOT)
  // explicit spectrum generation
  for (unsigned ns = 0; ns < cdata.cmsset.specs.size(); ++ ns)
  {
    Data::cpeptide_t cpeptide;
    Data::cstat_t cstat;

    cstream splot;
    ofstream shtml;
    SpecPlot cplot(cdata, cdata.cmsset, splot, shtml, cpeptide, cstat);

    cplot.sopt[SpecPlot::ezoom] = sopt[ezoom];
    cplot.sopt[SpecPlot::efont] = sopt[efont];
    cplot.sopt[SpecPlot::eoutdir] = sopt[eoutdir];
    cplot.sopt[SpecPlot::eoutfile] = sopt[eoutfile];
    cplot.sopt[SpecPlot::eprefix] = sopt[eprefix];
    cplot.sopt[SpecPlot::eformat] = sopt[eformat];
    cplot.sopt[SpecPlot::erange] = sopt[erange];
    cplot.bopt[SpecPlot::eannot] = bopt[eannot];
    cplot.bopt[SpecPlot::esuperpose] = true;
    cplot.sopt[SpecPlot::epeakmasstol] = sopt[epeakmasstol];
    cplot.bopt[SpecPlot::everbose] = bopt[everbose];
    cplot.bopt[SpecPlot::eshift] = bopt[eshift];
    cplot.bopt[SpecPlot::enotitle] = bopt[enotitle];

    // no header hack
    if (copt[epeptide].empty())
    {
      cplot.bopt[SpecPlot::enoheader] = true;
      copt[epeptide].push_back("");
    }

    unsigned ni = 0;

    // index
    if (! sopt[espectrum].empty())
    {
      unsigned nsi = 0;
      istringstream(sopt[espectrum]) >> nsi;
      -- nsi;

      if (nsi < cdata.cmsset.specs.size())
        ni = nsi;
      else
      {
        cerr << "error: invalid spectrum '" << sopt[espectrum] << "' " << endl;
        return -3;
      }
    }
    // scan
    else if (! sopt[espectrumscan].empty())
    {
      unsigned nss = 0;
      istringstream(sopt[espectrumscan]) >> nss;

      unsigned j;
      for (j = 0; j < cdata.cmsset.specs.size(); ++ j)
        if (nss == cdata.cmsset.specs[j].scan)
          break;

      if (j < cdata.cmsset.specs.size())
        ni = j;
      else
      {
        cerr << "error: invalid spectrum scan '" << sopt[espectrumscan] << "'" << endl;
        return -3;
      }

      cplot.sopt[SpecPlot::etitle] = "Scan " + sopt[espectrumscan];
    }
    // info
    else if (! sopt[espectruminfo].empty())
    {
      unsigned j;
      for (j = 0; j < cdata.cmsset.specs.size(); ++ j)
        if (cdata.cmsset.specs[j].info)
          if (sopt[espectruminfo] == cdata.cmsset.specs[j].info)
            break;

      if (j < cdata.cmsset.specs.size())
        ni = j;
      else
      {
        cerr << "error: invalid spectrum info '" << sopt[espectruminfo] << "'" << endl;
        return -3;
      }
    }
    // loop
    else
    {
      ni = ns;
    }

    for (unsigned i = 0; i < copt[epeptide].size(); ++ i)
      cpeptide[i][ni] = make_pair(string(), copt[epeptide][i]);

    if (cplot())
    {
      cout << "error: see '" << sopt[eoutdir] + "/spsplot.log'" << endl;
      return -4;
    }

    if (splot.dup2(stdin))
      abort();

    if (gnuplot_main(0, 0))
    {
      cout << "error: see '" << sopt[eoutdir] + "/spsplot.log'" << endl;
      return -4;
    }

    if (bopt[everbose])
    {
      cout << "m/z: " << cstat[ni].second[0] << endl;
      cout << "Annotation: " << cstat[ni].first << endl;
      cout << "Parent Mass Experimental: " << cstat[ni].second[1] << endl;

      if (std::isfinite(cstat[ni].second[2]))
        cout << "Parent Mass Theoretical: " << cstat[ni].second[2] << endl;
      else
        cout << "Parent Mass Theoretical: N/A" << endl;

      if (std::isfinite(cstat[ni].second[3]))
        cout << "Precursor m/z error: " << cstat[ni].second[0] - cstat[ni].second[3] << " Da, " << fabs(cstat[ni].second[0] - cstat[ni].second[3]) * 1e6 / cstat[ni].second[3] << " ppm" << endl;
      else
        cout << "Precursor m/z error: N/A" << endl;

      cout << "B (%): " << cstat[ni].second[4] << endl;
      cout << "Y (%): " << cstat[ni].second[5] << endl;
      cout << "BY Intensity (%): " << cstat[ni].second[6] << endl;
    }

    if (! sopt[espectrum].empty() || ! sopt[espectrumscan].empty() || ! sopt[espectruminfo].empty())
      break;
  }

  return 0;

#else // SPECPLOT
  /// statistics
  if (bopt[estats])
    return stats();

  // erase null consensus
  for (abinfo_t::iterator i = cdata.cabinfo.begin(), j = i; i != cdata.cabinfo.end(); i = j)
  {
    j ++;

    if (! cdata.cconsset.specs[i->first].peakList.size())
      cdata.cabinfo.erase(i->first);
  }

  // save requested peptide
  if (! copt[epeptide].empty())
    if (! sopt[espectrum].empty())
    {
      ofstream speptide_user((sopt[eoutdir] + "/spectrum." + sopt[espectrum] + ".txt").c_str());

      speptide_user << copt[epeptide].back() << "\n";
    }
    else if (! sopt[econtig].empty())
    {
      ofstream speptide_user((sopt[eoutdir] + "/contig." + sopt[econtig].substr(0, sopt[econtig].find('-')) + ".txt").c_str());

      speptide_user << copt[epeptide].back() << "\n";
    }

#if defined(__linux__)
  struct
  {
    int operator () ()
    {
      pid_t nw;
      int nstatus;

      do
      {
        nw = waitpid(-1, & nstatus, WUNTRACED | WCONTINUED);

        if (nw == -1)
          return nstatus;
        if (! WIFEXITED(nstatus))
          return nstatus;
      } while (! WIFEXITED(nstatus) && ! WIFSIGNALED(nstatus));

      return 0;
    }
  } cwait;
#endif

  // multiprocessor node and grid node
  if (sopt[enodes].empty())
  {
    unsigned nprocess = 0, nc = 0, ndone = 0;

    string s[] = {"running", "0"};
    ccron_t::iterator ic = ccron.insert(make_pair(& typeid(& SpsPlot::genindex), make_pair(this, make_pair(& SpsPlot::genindex, sps::vector<string>(s, s + sizeof(s) / sizeof(* s)))))).first;

    if (sopt[enodes].empty() || sopt[econtig] == "0")
    {
      nprocess ++;

      switch (fork())
      {
      case -1:
        abort();

      case 0: 
        if (genhtml())
          abort();

        _exit(0);
      }

      // skeleton generated
      if (sopt[econtig] == "0")
        return 0;

      // console
      cout << sopt[eoutdir] << "/index.html" << flush << endl;

#if defined(__linux__)
      ualarm(1, 10 * 1000000);
#endif
    } 

    abinfo_t::iterator ir[] =
    {
      cdata.cabinfo.begin(),
      -- cdata.cabinfo.end()
    };

    if (! sopt[econtig].empty())
    {
      unsigned nr[2] = {0, 0};
      const size_t np = sopt[econtig].find(':');

      if (np == string::npos)
      {
        istringstream(sopt[econtig]) >> nr[0];

        ir[0] = cdata.cabinfo.lower_bound(nr[0] - 1);
        ir[1] = ir[0];
      }
      else
      {
        istringstream(sopt[econtig].substr(0, np)) >> nr[0];
        istringstream(sopt[econtig].substr(np + 1)) >> nr[1];

        ir[0] = cdata.cabinfo.lower_bound(nr[0] - 1);
        ir[1] = cdata.cabinfo.lower_bound(max(nr[0], nr[1]) - 1);
      }
    }

    const unsigned nsize = distance(ir[0], ir[1]) + 1;

#if defined(__linux__)
    const unsigned ncpu = 2; //get_nprocs() * 2;
#endif

    // iterate contigs
    for (abinfo_t::iterator i = ir[0]; i != cdata.cabinfo.end() && i->first <= ir[1]->first; ++ i)
    {
      // index html
      {
        stringstream ss[2];
        ss[0] << ndone * 100 / nsize;

        if (ndone)
        {
          unsigned nremaining = (nsize - ndone) * (time(0) - ntime[0]) / ndone;
          ss[1] << nremaining / 60 << ":" << setfill('0') << setw(2) << nremaining % 60;
        }

        string s[] = {"running", ss[0].str(), ss[1].str()};
        ic->second.second.second = sps::vector<string>(s, s + sizeof(s) / sizeof(* s));
      }

#if defined(__linux__)
      nprocess ++;

      switch (fork())
      {
      case -1:
        abort();

      case 0: 
        if (gencontig(i->first))
          abort();

        _exit(0);
      }

      if (nprocess > ncpu) // parent
      {
        if (cwait())
          cout << "error: see '" << sopt[eoutdir] + "/spsplot.log'" << endl;

        nprocess --; 
        ndone ++;
      } 
#else 
      if (gencontig(i->first))
        abort();

      ndone ++;
#endif
    }

#if defined(__linux__)
    for (; nprocess; -- nprocess, ++ ndone)
    {
      if (cwait())
        cout << "error: see '" << sopt[eoutdir] + "/spsplot.log'" << endl;

      // index html
      {
        stringstream ss[2];
        ss[0] << ndone * 100 / nsize;

        if (ndone)
        {
          unsigned nremaining = (nsize - ndone) * (time(0) - ntime[0]) / ndone;
          ss[1] << nremaining / 60 << ":" << setfill('0') << setw(2) << nremaining % 60;
        }

        string s[] = {"running", ss[0].str(), ss[1].str()};
        ic->second.second.second = sps::vector<string>(s, s + sizeof(s) / sizeof(* s));
      }
    }
#endif

    // index html
    {
      string s[] = {"done"};
      ic->second.second.second = sps::vector<string>(s, s + sizeof(s) / sizeof(* s));
    }

    sigexp();
  }

#if defined(__linux__)
  // head node
  else
  {
    unsigned nnodes = 0;
    istringstream(sopt[enodes]) >> nnodes;

    /**
        Grid node spawning

        Runs spsplot on a grid node.
    */

    struct qsub_t
    {
      string & spath;
      string (& sopt)[esopts];
      sps::vector<string> (& copt)[ecopts];

      qsub_t(string & spath, string (& sopt)[esopts], sps::vector<string> (& copt)[ecopts]) : spath(spath), sopt(sopt), copt(copt)
      {
      }

      int operator () (const string & scontig)
      {
        switch (fork())
        {
        case -1:
          abort();

        case 0:
          list<string> cs;
          cs.push_back("qsub");
          cs.push_back("-l");
          cs.push_back("arch=lx26-x86");
          cs.push_back("-cwd");
          cs.push_back("-j");
          cs.push_back("y");
          cs.push_back("-S");
          cs.push_back("/bin/bash");
          cs.push_back("-b");
          cs.push_back("y");
          cs.push_back(spath);
          //cs.push_back("--");
          cs.push_back("--contig");

          cs.push_back(scontig);

          for (unsigned j = 0; j < sizeof(sopt) / sizeof(* sopt); ++ j)
            if (! sopt[j].empty())
            {
              cs.push_back(ssopt[(SpsPlot::sopt_t) (j)]);
              cs.push_back(sopt[j]);
            }

          for (unsigned j = 0; j < sizeof(copt) / sizeof(* copt); ++ j)
            if (! copt[j].empty())
            {
              cs.push_back(scopt[(SpsPlot::copt_t) (j)]);

              for (unsigned k = 0; k < copt[j].size(); ++ k)
                cs.push_back(copt[j][k]);
            }

          sps::vector<const char *> cp;
          cp.reserve(cs.size() + 1);

          for (list<string>::iterator j = cs.begin(); j != cs.end(); ++ j)
            cp.push_back(j->c_str());

          cp.push_back(0);

          //execvp(spath.c_str(), (char * const *) (& cp[0]));
          execvp("qsub", (char * const *) (& cp[0]));

          _exit(__LINE__);
        }

        return 0;
      }
    } cqsub(spath, sopt, copt);

    // spawns genhtml
    if (cqsub("0"))
      abort();

    if (cwait())
      abort();

    // iterate contigs
    double nd = max<double>(cdata.cabinfo.size() / double(nnodes), 1.);
    unsigned ni = 0;
    for (abinfo_t::iterator i = cdata.cabinfo.begin(), j = i, k; i != cdata.cabinfo.end(); i = j, ni = unsigned(nd), nd += cdata.cabinfo.size() / double(nnodes))
    {
      advance(j, min<unsigned>(distance(j, cdata.cabinfo.end()), unsigned(nd) - ni));
      k = j;
      k --;

      // spawns gencontig
      if (i->first <= k->first)
      {
        if (cqsub(static_cast<ostringstream &>(ostringstream() << i->first + 1).str() + ":" + static_cast<ostringstream &>(ostringstream() << k->first + 1).str()))
          abort();

        if (cwait())
          abort();
      }
    }
  }
#endif // __linux__

  return 0;
#endif // SPECPLOT
}


int SpsPlot::genheader(ostream & shtml)
{
  genHtmlHeader(shtml);

  shtml << "<head>\n";
  shtml << "  <script type=\"text/javascript\">\n";
  shtml << "  function common()\n";
  shtml << "  {\n";
  shtml << "    document.write('<input type=\"hidden\" name=\"" << spath << "\" value=\"\">');\n";

  for (unsigned j = 0; j < sizeof(sopt) / sizeof(* sopt); ++ j)
    if (! sopt[j].empty())
      switch ((SpsPlot::sopt_t) (j))
      {
      case eproject:
      case estars:
      case ecomp:
      case eseqs:
      case epklbin:
      case eindex:
      case emp:
      case emidx:
      case efasta:
      case efont:
      case eoutdir:
      case erefindex:
      case erefmp:
      case erefmidx:
      case ePageTitle:
      case htmlDefs:
        shtml << "    document.write('<input type=\"hidden\" name=\"" << ssopt[(SpsPlot::sopt_t) (j)] << "\" value=\"" << sopt[j] << "\">');\n";
        break;
      }

  for (unsigned j = 0; j < sizeof(copt) / sizeof(* copt); ++ j)
    if (! copt[j].empty())
      switch ((SpsPlot::copt_t) (j))
      {
      case ecluster:
      case eclusterms:
      case eclusterscan:
        shtml << "    document.write('<input type=\"hidden\" name=\"" << scopt[(SpsPlot::copt_t) (j)] << "\" value=\"";

        for (unsigned k = 0; k < copt[j].size(); ++ k)
          shtml << copt[j][k] << " ";

        shtml << "\">');\n";
        break;
      }

  shtml << "  }\n";
  shtml << "  </script>\n";
  shtml << "</head>\n";

  genHtmlHeader2(shtml);
  
  return 0;
}


void SpsPlot::genFooter(ostream & page)
{
//  page << "<div class=\"footer\">";
//  page << "<p>";
//  page << FOOTER_TEXT;
//  page << "</p>";
//  page << "</div>";
}

void SpsPlot::genTableFooter(ostream & page)
{
  page << "<table><tr><td>";
  page << "<p>";
  page << FOOTER_TEXT;
  page << "</p>";
  page << "</td></tr></table>";
}


#if ! defined(SPECPLOT)
/**
    Generates wrapper HTML pages

    Generates contigs.*.html, protein.*.html, cluster.*.html and cluster.txt.
*/

int SpsPlot::genhtml()
{
  // contigs
  for (abinfo_t::iterator i = cdata.cabinfo.begin(), j = i; i != cdata.cabinfo.end(); i = j)
  {
    advance(j, min<unsigned>(distance(j, cdata.cabinfo.end()), 20));

    ostringstream snc[2];
    snc[0] << i->first + 1;
    snc[1] << (-- abinfo_t::iterator(j))->first + 1;

    const string sfilename = sopt[eoutdir] + "/contigs." + snc[0].str() + ".html";
    ofstream shtml(sfilename.c_str());

    if (! shtml)
      return -2;


    genheader(shtml);
    shtml << "<body>\n";
    shtml << "  <div id=\"bodyWrapper\">\n";
    shtml << "  <div id=\"textWrapper\">\n";
    shtml << "  <br>\n";
    shtml << "    <div id=\"cyclelinks\" align=\"center\">";
    for (abinfo_t::iterator i = cdata.cabinfo.begin(), j = i; i != cdata.cabinfo.end(); i = j)
    {
      advance(j, min<unsigned>(distance(j, cdata.cabinfo.end()), 20));

      ostringstream snc[2];
      snc[0] << i->first + 1;
      snc[1] << (-- abinfo_t::iterator(j))->first + 1;

      shtml << "<a href=\"contigs." << snc[0].str() << ".html\">[" << snc[0].str() << "-" << snc[1].str() << "]</a> ";
    }
    shtml << "</div>\n";
    shtml << "  <br>\n";
    //shtml << "  <table class=\"sortable\" border=1 bordercolor=\"#003399\" cellpadding=\"4\" align=center style=\"width: 80%; border-collapse: collapse;border-top: 0.5pt solid black; border-bottom: 0.5pt solid black; border-left: 0.5pt solid black; border-right: 0.5pt solid black;\">\n";
    shtml << "  <table class=\"result sortable\" width=\"100%\" cellspacing=\"3\" cellpadding=\"1\">\n";
    shtml << "    <tr bgcolor=\"#003399\">\n";
    shtml << "      <th><span style=\"color:white\">Index</span></th>\n";
    shtml << "      <th><span style=\"color:white\">Spectra</span></th>\n";
    shtml << "      <th><span style=\"color:white\">Contig</span></th>\n";
    shtml << "      <th><span style=\"color:white\">Contig Sequence</span></th>\n";
    if (! sopt[efasta].empty())
      shtml << "      <th><span style=\"color:white\">Protein</span></th>\n";
    shtml << "    </tr>\n";

    ofstream splot;
    ContPlot cplot(cdata, splot, shtml);

    cplot.sopt[ContPlot::econtig] = snc[0].str() + ":" + snc[1].str();
    cplot.sopt[ContPlot::ezoom] = ".5";
    cplot.sopt[ContPlot::efont] = sopt[efont];
    cplot.sopt[ContPlot::eoutdir] = sopt[eoutdir];
    cplot.sopt[ContPlot::eprefix] = "contigs";
    cplot.sopt[ContPlot::eformat] = sopt[eformat];
    cplot.sopt[ContPlot::emp] = sopt[emp];
    cplot.sopt[ContPlot::efasta] = sopt[efasta];
    cplot.sopt[ContPlot::erefmp] = sopt[erefmp];
    cplot.sopt[ContPlot::efollow] = "contig";
    cplot.sopt[ContPlot::epeakmasstol] = sopt[epeakmasstol];

    try
    {
      if (cplot())
        throw;
    }
    catch (exception & e)
    {
      note(cerr, __PRETTY_FUNCTION__, __LINE__) << "exception caught on " << cplot << endl;
      cerr << e.what() << endl;
    }

    shtml << "  </table>\n";
    shtml << "  <br>\n";
    shtml << "  <div id=\"cyclelinks\" align=\"center\">";
    for (abinfo_t::iterator i = cdata.cabinfo.begin(), j = i; i != cdata.cabinfo.end(); i = j)
    {
      advance(j, min<unsigned>(distance(j, cdata.cabinfo.end()), 20));

      ostringstream snc[2];
      snc[0] << i->first + 1;
      snc[1] << (-- abinfo_t::iterator(j))->first + 1;

      shtml << "<a href=\"contigs." << snc[0].str() << ".html\">[" << snc[0].str() << "-" << snc[1].str() << "]</a> ";
    }
    genTableFooter(shtml);
    shtml << "</div>\n";
    shtml << "  <br>\n";
    shtml << "  </div>\n";
    shtml << "  </div>\n";
    genFooter(shtml);
    shtml << "</body>\n";

    shtml.flush();

    if (i == cdata.cabinfo.begin())
    {
      unlink((sopt[eoutdir] + "/contigs.html").c_str());
      string aux;
      aux = "cp ";
      aux += sfilename;
      aux += " ";
      aux += sopt[eoutdir];
      aux += "/contigs.html";
      system(aux.c_str());
      //symlink(sfilename.c_str(), (sopt[eoutdir] + "/contigs.html").c_str());
    }
  }

  // protein html
  if (! sopt[efasta].empty())
  {
    map< int, set<unsigned> > cfastamap;
    map< int, array<unsigned, 4> > nstat;

    // fasta / contig association
    for (abinfo_t::iterator i = cdata.cabinfo.begin(); i != cdata.cabinfo.end(); i ++)
    {
      int np;

      if (sopt[erefmp].empty() || ! cdata.reference(i->first))
        np = cdata.cproteinidx[i->first][0];
      else
        np = cdata.creferenceidx[cdata.crefindex[i->first]][0];

      cfastamap[np].insert(i->first);
    }

    // iterate proteins
    for (map< int, set<unsigned> >::iterator ip = cfastamap.begin(); ip != cfastamap.end(); ++ ip)
    {
      // skip missing fasta
      if (cdata.cfasta.size() <= ip->first)
      {
        cerr << "warning: missing fasta for contig " << ip->first << endl;
        continue;
      }

      ostringstream sfilename;
      sfilename << ip->first;

      const string ssfilename = sfilename.str();

      ofstream shtml((sopt[eoutdir] + "/protein." + ssfilename + ".html").c_str());

      if (! shtml)
        return -2;

      genheader(shtml);
      shtml << "<body>\n";
      shtml << "  <div id=\"bodyWrapper\">\n";
      shtml << "  <div id=\"textWrapper\">\n";
      shtml << "  <br>\n";
      shtml << "    <div id=\"cyclelinks\" align=\"center\"></div>\n";
      shtml << "  <br>\n";
      shtml << "    <div align=\"center\">\n";

      string sprotein = cdata.cfasta[ip->first];
      unsigned nb[] = {0, 0, 0, 0};
      map< range1d<unsigned>, set<unsigned> > cb;

      // iterate contigs
      for (set<unsigned>::iterator ic = ip->second.begin(); ic != ip->second.end(); ++ ic)
      {
        unsigned ni[2];

        if (sopt[erefmp].empty() || ! cdata.reference(* ic))
        {
          if (cdata.cprotset.specs[* ic].peakList.empty())
            continue;

          ni[0] = unsigned(cdata.cprotset.specs[* ic][0][1]);
          ni[1] = unsigned((* cdata.cprotset.specs[* ic].peakList.rbegin())[1]);
        }
        else
        {
          if (cdata.crefset.specs[cdata.crefindex[* ic]].peakList.empty())
            continue;

          ni[0] = unsigned(cdata.crefset.specs[cdata.crefindex[* ic]][0][1]);
          ni[1] = unsigned((* cdata.crefset.specs[cdata.crefindex[* ic]].peakList.rbegin())[1]);
        }

        nb[0] += cdata.cabinfo[* ic].first.first.size();
        cb[range1d<unsigned>(ni[0], ni[1] - 1)].insert(* ic);
      }

      stringstream scell[2];

      // iterate acids
      for (unsigned na = 0, ncol[2] = {0xAAAAAA}; na < sprotein.size(); ++ na)
      {
        if (na % 50 == 0)
        {
          scell[0] << "<br>" << na + 1;
          scell[1] << "<br>";
        }
        else if (na % 10 == 0)
          scell[1] << " ";

        ncol[1] = 0xAAAAAA;
        for (map< range1d<unsigned>, set<unsigned> >::iterator ib = cb.begin(); ib != cb.end(); ++ ib)
          if (ib->first == range1d<unsigned>(na, na))
          {
            ncol[1] = 0x000000;
            break;
          }

        if (ncol[0] != ncol[1])
        {
          scell[1] << "</font><font color=#" << hex << ncol[1] << ">";
          ncol[0] = ncol[1];
        }

        if (ncol[0] == 0x000000)
          nb[1] ++;

        scell[1] << sprotein[na];
      }

      string sname = cdata.cfasta.getID(ip->first);

      replace(sname.begin(), sname.end(), '|', ' ');
      transform(sname.begin(), sname.end(), sname.begin(), ::toupper);

      nstat[ip->first][0] = ip->second.size();
      nstat[ip->first][1] = nb[0];
      nstat[ip->first][2] = nb[1];
      nstat[ip->first][3] = nb[1] * 1000 / sprotein.size();

      shtml << "  <table class=\"result\" width=\"100%\" style=\"border-spacing: 0px;\">\n";
      shtml << "    <tr>\n";
      shtml << "      <td colspan=\"0\"><h2><i>" << sname << "</i></h2>\n";
      shtml << "      <hr><b>" << nstat[ip->first][0] << " contigs, " << nstat[ip->first][1] << " spectra, " << nstat[ip->first][2] << " amino acids, " << nstat[ip->first][3] / 10. << "% coverage</b></td>\n";
      shtml << "    <td></td>\n";
      shtml << "    </tr>\n";
      shtml << "    <tr>\n";
      shtml << "      <td align=\"right\"><tt>" << scell[0].rdbuf() << sps::clear << "</tt></td>\n";
      shtml << "      <td><tt><font color=#AAAAAA>" << scell[1].rdbuf() << sps::clear << "</font></tt></td>\n";
      shtml << "    </tr>\n";
      shtml << "  </table>\n";
      shtml << "  <a href=\"protein_details." << ssfilename << ".html\">Protein coverage</a>\n";

      shtml << "    </div>\n";
      shtml << "  <br>\n";
      //shtml << "  <table class=\"sortable\" border=1 bordercolor=\"#003399\" cellpadding=\"4\" align=center style=\"width: 80%; border-collapse: collapse;border-top: 0.5pt solid black; border-bottom: 0.5pt solid black; border-left: 0.5pt solid black; border-right: 0.5pt solid black;\">\n";
      shtml << "  <table class=\"result sortable\" width=\"100%\" cellspacing=\"3\" cellpadding=\"1\">\n";
      shtml << "    <tr bgcolor=\"#003399\">\n";
      shtml << "      <th><span style=\"color:white\">Index</span></th>\n";
      shtml << "      <th><span style=\"color:white\">Spectra</span></th>\n";
      shtml << "      <th><span style=\"color:white\">Contig</span></th>\n";
      shtml << "      <th><span style=\"color:white\">Contig Sequence</span></th>\n";
      shtml << "      <th><span style=\"color:white\">Protein</span></th>\n";
      shtml << "    </tr>\n";

      SpsPlot::Data::csequence_t csequence;

      // iterate contigs
      for (set<unsigned>::iterator ic = ip->second.begin(); ic != ip->second.end(); ++ ic)
      {
        stringstream sc;
        sc << * ic + 1;

        ofstream splot;
        ContPlot cplot(cdata, splot, shtml, csequence);

        cplot.sopt[ContPlot::ezoom] = ".5";
        cplot.sopt[ContPlot::efont] = sopt[efont];
        cplot.sopt[ContPlot::eoutdir] = sopt[eoutdir];
        cplot.sopt[ContPlot::eprefix] = "contigs";
        cplot.sopt[ContPlot::eformat] = sopt[eformat];
        cplot.sopt[ContPlot::emp] = sopt[emp];
        cplot.sopt[ContPlot::efasta] = sopt[efasta];
        cplot.sopt[ContPlot::erefmp] = sopt[erefmp];
        cplot.sopt[ContPlot::efollow] = "contig." + sc.str();
        cplot.sopt[ContPlot::epeakmasstol] = sopt[epeakmasstol];

        cplot.sopt[ContPlot::econtig] = sc.str();

        try
        {
          if (cplot())
            throw;
        }
        catch (exception & e)
        {
          note(cerr, __PRETTY_FUNCTION__, __LINE__) << "exception caught on " << cplot << endl;
          cerr << e.what() << endl;
        }
      }

      shtml << "  </table>\n";
      genTableFooter(shtml);
      shtml << "  <br>\n";
      shtml << "  </div>\n";
      shtml << "  </div>\n";
      genFooter(shtml);
      shtml << "</body>\n";

      // details
      {
        ofstream shtml((sopt[eoutdir] + "/protein_details." + ssfilename + ".html.old").c_str());

        genheader(shtml);
        shtml << "<body>\n";
        shtml << "  <div id=\"bodyWrapper\">\n";
        shtml << "  <div id=\"textWrapper\">\n";
        shtml << "  <table class=\"result\" width=\"100%\" style=\"border-spacing: 0px\">\n";
        shtml << "    <tr>\n";
        shtml << "      <td colspan=\"0\"><h2><i>" << sname << "</i></h2><hr></td>\n";
        shtml << "    </tr>\n";
        shtml << "  </table>\n";

        const unsigned nwidth = 20;

        // table parsing
        map< unsigned, map<unsigned, pair<unsigned, string> > > ctable;

        // table bounds
        unsigned nbound[] = {numeric_limits<unsigned>::max(), numeric_limits<unsigned>::min()};

        // iterate ranges
        for (map< range1d<unsigned>, set<unsigned> >::iterator ir = cb.begin(); ir != cb.end(); ++ ir)
        {
          if (nbound[0] > ir->first.value[0])
            nbound[0] = ir->first.value[0];
          if (nbound[1] < ir->first.value[1])
            nbound[1] = ir->first.value[1];
        }

        // iterate ranges
        for (map< range1d<unsigned>, set<unsigned> >::iterator ir = cb.begin(); ir != cb.end(); ++ ir)
          // iterate contigs
          for (set<unsigned>::iterator ic = ir->second.begin(); ic != ir->second.end(); ++ ic)
          {
            unsigned na = 0;

            // padding
            for (; na < ir->first.value[0]; ++ na)
              ctable[* ic][na] = make_pair(1U, string("<td />"));

            // segment
            bool bfirst = false;
            for (map< unsigned, pair<unsigned, string> >::iterator is = csequence[* ic].find(na - ir->first.value[0]); is != csequence[* ic].end(); na += is->second.first, ++ is)
            {
              string sstyle = "background-color: transparent;";

              if (! bfirst && is->second.first >= 1)
              {
                bfirst = true;
                sstyle += " border-left: solid 1px;";
              }

              if (is->second.first == 1)
                ctable[* ic][na] = make_pair(1U, string("<td style=\"" + sstyle + "\" align=\"center\">" + is->second.second + "</td>"));
              else if (is->second.first > 1)
              {
                ostringstream sspan;
                sspan << is->second.first;

                ctable[* ic][na] = make_pair(is->second.first, string("<td style=\"" + sstyle + "\" align=\"center\" colspan=\"" + sspan.str() + "\">" + is->second.second + "</td>"));
              }
            }

            // padding
            for (; na <= nbound[1]; ++ na)
              ctable[* ic][na] = make_pair(1U, string("<td />"));
          }

        // iterate tables
        for (unsigned nt = nbound[0]; nt <= nbound[1]; nt += nwidth)
        {
          shtml << "  <table class=\"result\" width=\"100%\" style=\"border-spacing: 1px; background-color: #CCCFFF\">\n";

          shtml << "    <tr>";
          shtml << "<td width=\"5%\"/>";

          // iterate acids
          for (unsigned na = nt; na < nt + nwidth && na <= nbound[1]; ++ na)
            shtml << "<th align=\"center\">" << sprotein[na] << "</th>";

          shtml << "</tr>\n";

          // iterate ranges
          for (map< range1d<unsigned>, set<unsigned> >::iterator ir = cb.begin(); ir != cb.end(); ++ ir)
            // iterate contigs
            for (set<unsigned>::iterator ic = ir->second.begin(); ic != ir->second.end(); ++ ic)
            {
              // skip row
              if (nt + nwidth <= ir->first.value[0])
                continue;

              // skip row
              if (ir->first.value[1] + 1 <= nt)
                continue;

              shtml << "    <tr><th align=\"right\"><a href=\"contig." << *ic + 1 << ".html\" style=\"color: white\">" << * ic + 1 << "</a></th>";

              unsigned na = nt;

              for (map< unsigned, pair<unsigned, string> >::iterator ia = ctable[* ic].find(na); na < nt + nwidth && na <= nbound[1]; ia = ctable[* ic].find(na))
                if (ia != ctable[* ic].end())
                {
                  shtml << ia->second.second;
                  na += ia->second.first;
                }
                else
                {
                  shtml << "<td style=\"background-color: transparent\"/>\n";
                  ++ na;
                }

              shtml << "</tr>\n";
            }

          shtml << "  </table>\n";
          shtml << "  <br>\n";
        }
        genTableFooter(shtml);

        shtml << "  </div>\n";
        shtml << "  </div>\n";
        genFooter(shtml);
        shtml << "</body>\n";
      }
    }

    // proteins
    ofstream shtml((sopt[eoutdir] + "/proteins.html").c_str());

    if (! shtml)
      return -2;

    genheader(shtml);
    shtml << "<body>\n";
    shtml << "  <div id=\"bodyWrapper\">\n";
    shtml << "  <div id=\"textWrapper\">\n";
    shtml << "  <br>\n";
    shtml << "    <div id=\"cyclelinks\" align=\"center\"></div>\n";
    shtml << "  <br>\n";
    //shtml << "  <table class=\"sortable\" border=1 bordercolor=\"#003399\" cellpadding=\"4\" align=center style=\"width: 80%; border-collapse: collapse;border-top: 0.5pt solid black; border-bottom: 0.5pt solid black; border-left: 0.5pt solid black; border-right: 0.5pt solid black;\">\n";
    shtml << "  <table class=\"result sortable\" width=\"100%\" cellspacing=\"3\" cellpadding=\"1\">\n";
    shtml << "    <tr bgcolor=\"#003399\">\n";
    shtml << "      <th><span style=\"color:white\">Protein</span></th>\n";
    shtml << "      <th><span style=\"color:white\">Description</span></th>\n";
    shtml << "      <th><span style=\"color:white\">Contigs</span></th>\n";
    shtml << "      <th><span style=\"color:white\">Spectra</span></th>\n";
    shtml << "      <th><span style=\"color:white\">Amino Acids</span></th>\n";
    shtml << "      <th><span style=\"color:white\">Coverage (%)</span></th>\n";
    shtml << "    </tr>\n";

    // iterate proteins
    for (map< int, set<unsigned> >::iterator ip = cfastamap.begin(); ip != cfastamap.end(); ++ ip)
    {
      // skip missing fasta
      if (cdata.cfasta.size() <= ip->first)
      {
        cerr << "warning: missing fasta for " << ip->first << endl;
        continue;
      }

      ostringstream sfilename;
      sfilename << ip->first;

      const string ssfilename = sfilename.str();

      shtml << "  <tr>\n";
      shtml << "    <td><a href=\"protein." << ssfilename << ".html\"><i>" << cdata.cfasta.getID(ip->first) << "</i></a></td>\n";
      shtml << "    <td>" << cdata.cfasta.getDesc(ip->first) << "</td>\n";
      shtml << "    <td>" << nstat[ip->first][0] << "</td>\n";
      shtml << "    <td>" << nstat[ip->first][1] << "</td>\n";
      shtml << "    <td>" << nstat[ip->first][2] << "</td>\n";
      shtml << "    <td>" << nstat[ip->first][3] / 10. << "</td>\n";
      shtml << "  </tr>\n";
    }

    shtml << "  </table>\n";
    genTableFooter(shtml);
    shtml << "  <br>\n";
    shtml << "  </div>\n";
    shtml << "  </div>\n";
    genFooter(shtml);
    shtml << "</body>\n";
  }

  // input filenames
  if (! copt[ecluster].empty())
  {
    ofstream stxt((sopt[eoutdir] + "/cluster.txt").c_str());
    stxt << "Scan\tContig\tFilename\n";

//cout << "cdata.cclusteridx.size(): " << cdata.cclusteridx.size() << endl;

    map<string, string> chtml[cdata.cclusteridx.size()];

    // iterate contigs
    for (abinfo_t::iterator i = cdata.cabinfo.begin(); i != cdata.cabinfo.end(); ++ i)
    {
      unsigned nc = i->first;

      const unsigned nsize = cdata.cabinfo[nc].first.first.size();

      for (unsigned nb = 0; nb <= nsize / 100; ++ nb)
      {
        Data::cpeptide_t cpeptide;

        stringstream snc;
        snc << nc + 1;

        if (nsize > 100)
          snc << "-" << nb;

        const string ssnc = snc.str();

        {
          ofstream splot;
          ofstream shtml;
          ContPlot cplot(cdata, splot, shtml, cpeptide);

          cplot.sopt[ContPlot::econtig] = ssnc;
          cplot.sopt[ContPlot::efont] = sopt[efont];
          cplot.sopt[ContPlot::eoutdir] = sopt[eoutdir];
          cplot.sopt[ContPlot::eprefix] = "contig";
          cplot.sopt[ContPlot::eformat] = sopt[eformat];
          cplot.sopt[ContPlot::emp] = sopt[emp];
          cplot.sopt[ContPlot::emidx] = sopt[emidx];
          cplot.sopt[ContPlot::efasta] = sopt[efasta];
          cplot.sopt[ContPlot::erefmp] = sopt[erefmp];
          cplot.sopt[ContPlot::erefmidx] = sopt[erefmidx];
          cplot.sopt[ContPlot::epeakmasstol] = sopt[epeakmasstol];

          try
          {
            if (cplot())
              throw;
          }
          catch (exception & e)
          {
            note(cerr, __PRETTY_FUNCTION__, __LINE__) << "exception caught on " << cplot << endl;
            cerr << e.what() << endl;
          }
        }

        // iterate spectra
        map<unsigned, pair<string, string> >::iterator is[] = {cpeptide[0].begin(), cpeptide[0].begin()};

        if (is[0] != cpeptide[0].end())
          advance(is[0], min(nsize, nb * 100) - nb * 100);

        if (is[1] != cpeptide[0].end())
          advance(is[1], min(nsize, (nb + 1) * 100) - nb * 100);

        for (; is[0] != is[1]; ++ is[0])
        {
          const unsigned ns = is[0]->first;

          stringstream sns;
          sns << ns + 1 << "." << cpeptide.rbegin()->first;

          const string ssns = sns.str();
          const unsigned nsize = cdata.cabinfo[nc].first.first.size();

          stringstream snc;
          snc << nc + 1;

          if (nsize > 100)
            snc << "-" << nb;

          // iterate cluster
          for (list< pair<unsigned, unsigned> >::iterator ie = cdata.ccluster.find(ns)->second.begin(); ie != cdata.ccluster.find(ns)->second.end(); ++ ie)
          {
            // spectrum
            Spectrum & spec = cdata.cspecmap[ie->first].specs[ie->second];

            if (cdata.cclusterscan.empty())
              stxt << ie->first << "-" << ie->second << '\t' << nc + 1 << '\t' << cdata.cclusteridx[ie->first] << '\n';
            else
              stxt << cdata.cclusterscan[ie->first][ie->second][0] << '\t' << nc + 1 << '\t' << cdata.cclusteridx[ie->first] << '\n';

            // html
            ostringstream srow;

            ostringstream sns2;
            sns2 << "." << ie->first << "." << ie->second;

            const string ssns2 = sns2.str();

            ostringstream sscan;
            sscan << ie->first << "-" << ie->second;

            {
              ofstream splot;
              SpecPlot cplot(cdata, cdata.cspecset, splot, srow, cpeptide);

              cplot.sopt[SpecPlot::eprefix] = "c_spectrum";
              cplot.sopt[SpecPlot::espectrum] = ssns;
              cplot.sopt[SpecPlot::eindex] = sopt[eindex];
              cplot.sopt[SpecPlot::escan] = sscan.str();
              cplot.sopt[SpecPlot::efasta] = sopt[efasta];
              cplot.sopt[SpecPlot::eformat] = sopt[eformat];
              cplot.sopt[SpecPlot::epeakmasstol] = sopt[epeakmasstol];

              cplot.copt[SpecPlot::ecluster] = copt[ecluster];
              cplot.copt[SpecPlot::eclusterms] = copt[eclusterms];
              cplot.copt[SpecPlot::eclusterscan] = copt[eclusterscan];

              try
              {
                if (cplot())
                  throw;
              }
              catch (exception & e)
              {
                note(cerr, __PRETTY_FUNCTION__, __LINE__) << "exception caught on " << cplot << endl;
                cerr << e.what() << endl;
              }
            }
  
//            cout << "ie->first: " << ie->first << endl;
//            cout << "sscan: " << sscan.str() << endl;
//            cout << "srow: " << srow.str() << endl;
            

            chtml[ie->first].insert(make_pair(sscan.str(), srow.str()));

//            cout << "#2" << endl;
          }
        }
      }
    }


    ofstream shtml[cdata.cclusteridx.size()];

    // iterate cluster
    for (unsigned ic = 0; ic < cdata.cclusteridx.size(); ++ ic)
    {
      ostringstream snf;
      snf << ic;

      // iterate filename
      for (map<string, string>::iterator i = chtml[ic].begin(), j = i; i != chtml[ic].end(); i = j)
      {
        advance(j, min<unsigned>(distance(j, chtml[ic].end()), 20));

        ostringstream snc[2];
        snc[0] << i->first;
        snc[1] << (-- map<string, string>::iterator(j))->first;

        const string sfilename = sopt[eoutdir] + "/r_cluster." + snf.str() + "." + snc[0].str() + ".html";

        ofstream shtml(sfilename.c_str());

        if (! shtml)
          return -2;


        genheader(shtml);
        shtml << "<body>\n";
        shtml << "  <div id=\"bodyWrapper\">\n";
        shtml << "  <div id=\"textWrapper\">\n";
        shtml << "  <br>\n";
        shtml << "    <div align=\"center\">\n";
        shtml << "    <table class=\"result\" width=\"100%\" style=\"border-spacing: 0px;\">\n";
        shtml << "      <tr>\n";
        shtml << "        <td colspan=\"0\"><h2><i>" << cdata.cclusteridx[ic] << "</i></h2>\n";
        shtml << "        <hr></td>\n";
        shtml << "      </tr>\n";
        shtml << "    </table>\n";

        shtml << "    </div>\n";
        shtml << "  <br>\n";
        shtml << "    <div id=\"cyclelinks\" align=\"center\">";
        for (map<string, string>::iterator m = chtml[ic].begin(), n = m; m != chtml[ic].end(); m = n)
        {
          advance(n, min<unsigned>(distance(n, chtml[ic].end()), 20));

          ostringstream snc[2];
          snc[0] << m->first;
          snc[1] << (-- map<string, string>::iterator(n))->first;

          shtml << "<a href=\"r_cluster." << snf.str() << "." << snc[0].str() << ".html\">[" << snc[0].str() << "-" << snc[1].str() << "]</a> ";
        }
        shtml << "</div>\n";
        shtml << "  <br>\n";
        //shtml << "  <table class=\"sortable\" border=1 bordercolor=\"#003399\" cellpadding=\"4\" align=center style=\"width: 80%; border-collapse: collapse;border-top: 0.5pt solid black; border-bottom: 0.5pt solid black; border-left: 0.5pt solid black; border-right: 0.5pt solid black;\">\n";
        shtml << "  <table class=\"result sortable\" width=\"100%\" cellspacing=\"3\" cellpadding=\"1\">\n";
        shtml << "    <tr bgcolor=\"#003399\">\n";
        shtml << "      <th><span style=\"color:white\">Scan</span></th>\n";
        shtml << "      <th><span style=\"color:white\">Spectrum</span></th>\n";
        shtml << "      <th><span style=\"color:white\">Peptide</span></th>\n";
        shtml << "      <th><span style=\"color:white\">Mass (m)</span></th>\n";
        shtml << "      <th><span style=\"color:white\">Charge (z)</span></th>\n";
        shtml << "      <th><span style=\"color:white\">B (%)</span></th>\n";
        shtml << "      <th><span style=\"color:white\">Y (%)</span></th>\n";
        shtml << "      <th><span style=\"color:white\">BY Intensity (%)</span></th>\n";
        shtml << "    </tr>\n";

        for (map<string, string>::iterator k = i; k != j; ++ k)
          shtml << k->second;

        shtml << "  </table>\n";
        shtml << "  <br>\n";
        shtml << "  <div id=\"cyclelinks\" align=\"center\">";
        for (map<string, string>::iterator m = chtml[ic].begin(), n = m; m != chtml[ic].end(); m = n)
        {
          advance(n, min<unsigned>(distance(n, chtml[ic].end()), 20));

          ostringstream snc[2];
          snc[0] << m->first;
          snc[1] << (-- map<string, string>::iterator(n))->first;

          shtml << "<a href=\"r_cluster." << snf.str() << "." << snc[0].str() << ".html\">[" << snc[0].str() << "-" << snc[1].str() << "]</a> ";
        }
        genTableFooter(shtml);
        shtml << "</div>\n";
        shtml << "  <br>\n";
        shtml << "  </div>\n";
        shtml << "  </div>\n";
        genFooter(shtml);
        shtml << "</body>\n";

        shtml.flush();

        if (i == chtml[ic].begin())
        {
          unlink((sopt[eoutdir] + "/r_cluster." + snf.str() + ".html").c_str());
          string aux;
          aux = "cp ";
          aux += sfilename;
          aux += " ";
          aux += sopt[eoutdir];
          aux += "/r_cluster.";
          aux += snf.str();
          aux += ".html";
          system(aux.c_str());
//          symlink(sfilename.c_str(), (sopt[eoutdir] + "/r_cluster." + snf.str() + ".html").c_str());
        }
      }
    }
  }

  return 0;
}
#endif // SPECPLOT


#if ! defined(SPECPLOT)
/**
    Generates index.html

    Generates index.html periodically.

    @param  argv  Dynamic parameters (state, status, remaining time)

    @note   Async-safe function.
*/

int SpsPlot::genindex(const sps::vector<std::string> & argv)
{
  ntime[1] = time(0) - ntime[0];

  char sfilename[sopt[eoutdir].size() + sizeof("/index.html") + 1];

  strncpy(sfilename, sopt[eoutdir].data(), sopt[eoutdir].size());
  sfilename[sopt[eoutdir].size()] = '\0';
  strcat(sfilename, "/index.html");

  ofstream shtml(sfilename);

  if (! shtml)
    return -2;

  genHtmlHeader(shtml);
  genHtmlHeader2(shtml);

  shtml << (argv[0] == "running" ? "<body onload=\"doLoad()\">" : "<body>") << '\n';
  shtml << "  <br>\n";
  shtml << "  <div id=\"bodyWrapper\">\n";
  shtml << "  <div id=\"textWrapper\">\n";
  shtml << "  <hr>\n";
  //shtml << "  <table class=\"mainform\" border=\"0\" cellspacing=\"0\" cellpadding=\"2\" align=\"center\" width=\"100%\">\n";
  //shtml << "  <table class=\"mainform\" border=0 cellpadding=\"4\" align=center style=\"width: 80%; border-collapse: collapse;\">\n";
  shtml << "  <table class=\"mainform\">\n";
  shtml << "    <tr>\n";
  shtml << "      <th colspan=\"0\">Job Status</th>\n";
  shtml << "    </tr>\n";
  shtml << "    <tr><td>\n";
  //shtml << "  <table class=\"sched\" border=\"1\" cellspacing=\"1\" cellpadding=\"4\" width=\"100%\">\n";
  //shtml << "  <table class=\"sched\" border=1 bordercolor=\"#003399\" cellpadding=\"4\" align=center style=\"width: 100%; border-collapse: collapse;border-top: 0.5pt solid black; border-bottom: 0.5pt solid black; border-left: 0.5pt solid black; border-right: 0.5pt solid black;\">\n";
  shtml << "  <table class=\"sched\" width=\"100%\">\n";
  shtml << "    <tr>\n";
  shtml << "      <th width=\"25%\" bgcolor=\"#003399\"><span style=\"color:white\">Job</span></th>\n";
  shtml << "      <td>" << sopt[eproject] << "</td>\n";
  shtml << "    </tr>\n";
  shtml << "    <tr>\n";
  shtml << "      <th width=\"25%\" bgcolor=\"#003399\"><span style=\"color:white\">User</span></th>\n";
  shtml << "      <td>" << cuserid(0) << "</td>\n";
  shtml << "    </tr>\n";
  shtml << "    <tr>\n";
  shtml << "      <th width=\"25%\" bgcolor=\"#003399\"><span style=\"color:white\">Status</span></th>\n";
  if (argv[0] == "running")
    shtml << "      <td style=\"background-color:" << "lightyellow" << ";\">" << argv[1] << " %</td>\n";
  else if (argv[0] == "done")
    shtml << "      <td style=\"background-color:" << "lightgreen" << ";\">" << 100 << " %</td>\n";
  else
    shtml << "      <td style=\"background-color:" << "red" << ";\">" << argv[1] << " %</td>\n";
  shtml << "    </tr>\n";
  shtml << "    <tr>\n";
  shtml << "      <th width=\"25%\" bgcolor=\"#003399\"><span style=\"color:white\">Elapsed</span></th>\n";
  shtml << "      <td>" << ntime[1] / 60 << ":" << setfill('0') << setw(2) << ntime[1] % 60;
  if (argv.size() > 2 && ! argv[2].empty())
    shtml << " (" << argv[2] << " remaining)";
  shtml << "</td>\n";
  shtml << "    </tr>\n";
  shtml << "    <tr>\n";
  shtml << "      <th width=\"25%\" bgcolor=\"#003399\"><span style=\"color:white\">Log</span></th>\n";
  if (argv[0] == "running")
    shtml << "      <td><a href=\"spsplot.log\">Warnings</a></td>\n";
  else if (argv[0] == "done")
    shtml << "      <td><a href=\"spsplot.log\">Warnings</a></td>\n";
  else
    shtml << "      <td><a href=\"spsplot.log\">Errors</a></td>\n";
  shtml << "    </tr>\n";
  shtml << "    <tr>\n";
  if (sopt[efasta].empty())
  {
    shtml << "      <th width=\"25%\" bgcolor=\"#003399\"><span style=\"color:white\">Data</span></th>\n";
    shtml << "      <td><a href=\"contigs.html\">Group by Contig</a></td>\n";
  }
  else
  {
    shtml << "      <th width=\"25%\" bgcolor=\"#003399\" rowspan=\"2\"><span style=\"color:white\">Data</span></th>\n";
    shtml << "      <td><a href=\"contigs.html\">Group by Contig</a></td>\n";
    shtml << "    </tr>\n";
    shtml << "    <tr>\n";
    shtml << "      <td><a href=\"proteins.html\">Group by Protein</a></td>\n";
  }
  shtml << "    </tr>\n";
  if (! copt[ecluster].empty())
  {
    shtml << "    <tr>\n";
    shtml << "      <th width=\"25%\" bgcolor=\"#003399\"><span style=\"color:white\">Cluster Data</span></th>\n";
    shtml << "      <td>\n";
    shtml << "        <a href=\"cluster.txt\">All Clusters (txt)</a><br>\n";

    for (unsigned i = 0; i < cdata.cclusteridx.size(); ++ i)
      shtml << "      <a href=\"r_cluster." << i << ".html\">Group by <i>" << cdata.cclusteridx[i] << "</i></a><br>\n";

    shtml << "      </td>\n";
    shtml << "    </tr>\n";
  }
  genTableFooter(shtml);
  shtml << "  </table>\n";
  shtml << "    </td></tr>\n";
  shtml << "    <tr>\n";
  shtml << "      <td colspan=\"0\" class=\"bottomline\">&nbsp;</td>\n";
  shtml << "    </tr>\n";
  shtml << "  </table>\n";
  if (argv[0] != "running")
  {
    shtml << "  <br>\n";
    shtml << "  <div align=\"center\"><input type=button value=\"Refresh\" onClick=\"refresh()\"></div>\n";
  }
  shtml << "  <br>\n";
  shtml << "  </div>\n";
  shtml << "  </div>\n";
  genFooter(shtml);
  shtml << "</body>\n";

  return 0;
}
#endif // SPECPLOT


#if ! defined(SPECPLOT)
int SpsPlot::gencontig(unsigned nc)
{
  unsigned nspectrum = 0;
  istringstream(sopt[espectrum]) >> nspectrum;
  -- nspectrum;

  const unsigned nsize = cdata.cabinfo[nc].first.first.size();

  for (unsigned nb = 0; nb <= nsize / 100; ++ nb)
  {
    Data::cpeptide_t cpeptide;

    stringstream snc;
    snc << nc + 1;

    if (nsize >= 100)
      snc << "-" << nb;

    const string ssnc = snc.str();

    ofstream shtml((sopt[eoutdir] + "/contig." + ssnc + ".html").c_str());

    if (! shtml)
      return -2;

    genheader(shtml);
    shtml << "<body>\n";
    shtml << "  <div id=\"bodyWrapper\">\n";
    shtml << "  <div id=\"textWrapper\">\n";

    {
      cstream splot;
      ofstream shtml;
      ContPlot cplot(cdata, splot, shtml);

      cplot.sopt[ContPlot::econtig] = ssnc;
      cplot.sopt[ContPlot::ezoom] = ".5";
      cplot.sopt[ContPlot::efont] = sopt[efont];
      cplot.sopt[ContPlot::eoutdir] = sopt[eoutdir];
      cplot.sopt[ContPlot::eprefix] = "contigs";
      cplot.sopt[ContPlot::eformat] = sopt[eformat];
      cplot.sopt[ContPlot::emp] = sopt[emp];
      cplot.sopt[ContPlot::efasta] = sopt[efasta];
      cplot.sopt[ContPlot::erefmp] = sopt[erefmp];
      cplot.sopt[ContPlot::epeakmasstol] = sopt[epeakmasstol];

      if (cplot() || splot.dup2(stdin) || gnuplot_main(0, 0))
      {
        note(cerr) << cplot << endl;
        abort();
      }
    }

    {
      cstream splot;
      ofstream shtml;
      ContPlot cplot(cdata, splot, shtml, cpeptide);

      cplot.sopt[ContPlot::econtig] = ssnc;
      cplot.sopt[ContPlot::efont] = sopt[efont];
      cplot.sopt[ContPlot::eoutdir] = sopt[eoutdir];
      cplot.sopt[ContPlot::eprefix] = "contig";
      cplot.sopt[ContPlot::eformat] = sopt[eformat];
      cplot.sopt[ContPlot::emp] = sopt[emp];
      cplot.sopt[ContPlot::emidx] = sopt[emidx];
      cplot.sopt[ContPlot::efasta] = sopt[efasta];
      cplot.sopt[ContPlot::erefmp] = sopt[erefmp];
      cplot.sopt[ContPlot::erefmidx] = sopt[erefmidx];
      cplot.sopt[ContPlot::epeakmasstol] = sopt[epeakmasstol];

      if (cplot() || splot.dup2(stdin) || gnuplot_main(0, 0))
      {
        note(cerr) << cplot << endl;
        abort();
      }
    }

    shtml << "  <br>\n";
    shtml << "  <table align=\"center\">\n";
    shtml << "    <tr align=\"center\">\n";
    shtml << "      <td><img src=\"contig." << ssnc << ".png\"/></td>\n";
    shtml << "    </tr>\n";
    shtml << "    <tr align=\"center\">\n";
    shtml << "      <td>\n";
    shtml << "        <br>\n";
//    shtml << "        <form method=\"POST\" action=\"/cgi-bin/spsplot.fcgi\">";
//    shtml << "          <script type=\"text/javascript\">common();</script>\n";
//    shtml << "          <input type=\"hidden\" name=\"--contig\" value=\"" << ssnc << "\" />\n";
//    shtml << "          <input type=\"text\" name=\"--peptide\" style=\"text-transform: uppercase; width:50%\" /><br>\n";
//    shtml << "          <input type=submit value=\"Update\"/>\n";
//    shtml << "        </form>\n";
    shtml << "      </td>\n";
    shtml << "    </tr>\n";
    shtml << "  </table>\n";
    shtml << "  <br>\n";
    //shtml << "  <table border=1 bordercolor=\"#003399\" cellpadding=\"4\" align=center class=\"sortable\" style=\"width: 80%; border-collapse: collapse;border-top: 0.5pt solid black; border-bottom: 0.5pt solid black; border-left: 0.5pt solid black; border-right: 0.5pt solid black;\">\n";
    shtml << "  <table class=\"result sortable\" width=\"100%\" cellspacing=\"3\" cellpadding=\"1\">\n";
    shtml << "    <tr bgcolor=\"#003399\">\n";
    shtml << "      <th><span style=\"color:white\">Index</span></th>\n";
    shtml << "      <th><span style=\"color:white\">Spectrum</span></th>\n";
    shtml << "      <th><span style=\"color:white\">Peptide</span></th>\n";
    shtml << "      <th><span style=\"color:white\">Mass (m)</span></th>\n";
    shtml << "      <th><span style=\"color:white\">Charge (z)</span></th>\n";
    shtml << "      <th><span style=\"color:white\">B (%)</span></th>\n";
    shtml << "      <th><span style=\"color:white\">Y (%)</span></th>\n";
    shtml << "      <th><span style=\"color:white\">BY Intensity (%)</span></th>\n";
    shtml << "    </tr>\n";

    {
      cstream splot;
      SpecPlot cplot(cdata, cdata.cmsset, splot, shtml, cpeptide);

      cplot.sopt[SpecPlot::econtig] = ssnc;
      cplot.sopt[SpecPlot::ezoom] = ".4";
      cplot.sopt[SpecPlot::efont] = sopt[efont];
      cplot.sopt[SpecPlot::eoutdir] = sopt[eoutdir];
      cplot.sopt[SpecPlot::eprefix] = "spectrum";
      cplot.sopt[SpecPlot::eformat] = sopt[eformat];
      cplot.sopt[SpecPlot::epeakmasstol] = sopt[epeakmasstol];

      if (! copt[ecluster].empty())
        cplot.sopt[SpecPlot::efollow] = "cluster";

      if (cplot() || splot.dup2(stdin) || gnuplot_main(0, 0))
      {
        note(cerr) << cplot << endl;
        abort();
      }
    }

    shtml << "  </table>\n";
    genTableFooter(shtml);
    shtml << "  <br>\n";
    shtml << "  </div>\n";
    shtml << "  </div>\n";
    genFooter(shtml);
    shtml << "</body>\n";

    {
      cstream splot;
      ofstream shtml;
      SpecPlot cplot(cdata, cdata.cmsset, splot, shtml, cpeptide);

      cplot.sopt[SpecPlot::econtig] = ssnc;
      cplot.sopt[SpecPlot::ezoom] = "1";
      cplot.sopt[SpecPlot::efont] = sopt[efont];
      cplot.sopt[SpecPlot::eoutdir] = sopt[eoutdir];
      cplot.sopt[SpecPlot::eprefix] = "l_spectrum";
      cplot.sopt[SpecPlot::eformat] = sopt[eformat];
      cplot.sopt[SpecPlot::epeakmasstol] = sopt[epeakmasstol];

      if (cplot() || splot.dup2(stdin) || gnuplot_main(0, 0))
      {
        note(cerr) << cplot << endl;
        abort();
      }
    }

    // iterate spectra
    map<unsigned, pair<string, string> >::iterator is[] = {cpeptide[0].begin(), cpeptide[0].begin()};

    if (is[0] != cpeptide[0].end())
      advance(is[0], min(nsize, nb * 100) - nb * 100);

    if (is[1] != cpeptide[0].end())
      advance(is[1], min(nsize, (nb + 1) * 100) - nb * 100);

    for (; is[0] != is[1]; ++ is[0])
    {
      const unsigned ns = is[0]->first;

      // speedup
      if (! sopt[espectrum].empty() && nspectrum != ns)
        continue;

      stringstream sns;
      sns << ns + 1 << "." << cpeptide.rbegin()->first;

      const string ssns = sns.str();
      const unsigned nsize = cdata.cabinfo[nc].first.first.size();

      stringstream snc;
      snc << nc + 1;

      if (nsize > 100)
        snc << "-" << nb;

      const string ssnc = snc.str();

      ofstream shtml((sopt[eoutdir] + "/cluster." + ssnc + "." + ssns + ".html").c_str());

      if (! shtml)
        return -2;

      genheader(shtml);
      shtml << "<body>\n";
      shtml << "  <div id=\"bodyWrapper\">\n";
      shtml << "  <div id=\"textWrapper\">\n";
      shtml << "  <br>\n";
      shtml << "      <h3 align=\"center\"><img src=\"l_spectrum." << ssnc << "." << ssns << ".png\"/></h3>\n";
      shtml << "  <br>\n";
      //shtml << "  <table border=1 bordercolor=\"#003399\" cellpadding=\"4\" align=center class=\"sortable\" style=\"width: 80%; border-collapse: collapse;border-top: 0.5pt solid black; border-bottom: 0.5pt solid black; border-left: 0.5pt solid black; border-right: 0.5pt solid black;\">\n";
      shtml << "  <table class=\"result sortable\" width=\"100%\" cellspacing=\"3\" cellpadding=\"1\">\n";
      shtml << "    <tr bgcolor=\"#003399\">\n";
      shtml << "      <th><span style=\"color:white\">Scan</span></th>\n";
      shtml << "      <th><span style=\"color:white\">Spectrum</span></th>\n";
      shtml << "      <th><span style=\"color:white\">Peptide</span></th>\n";
      shtml << "      <th><span style=\"color:white\">Mass (m)</span></th>\n";
      shtml << "      <th><span style=\"color:white\">Charge (z)</span></th>\n";
      shtml << "      <th><span style=\"color:white\">B (%)</span></th>\n";
      shtml << "      <th><span style=\"color:white\">Y (%)</span></th>\n";
      shtml << "      <th><span style=\"color:white\">BY Intensity (%)</span></th>\n";
      shtml << "    </tr>\n";

      {
        cstream splot;
        SpecPlot cplot(cdata, cdata.cspecset, splot, shtml, cpeptide);

        cplot.sopt[SpecPlot::econtig] = ssnc;
        cplot.sopt[SpecPlot::ezoom] = ".4";
        cplot.sopt[SpecPlot::efont] = sopt[efont];
        cplot.sopt[SpecPlot::eoutdir] = sopt[eoutdir];
        cplot.sopt[SpecPlot::eprefix] = "c_spectrum";
        cplot.sopt[SpecPlot::espectrum] = ssns;
        cplot.sopt[SpecPlot::eindex] = sopt[eindex];
        cplot.sopt[SpecPlot::eformat] = sopt[eformat];
        cplot.sopt[SpecPlot::epeakmasstol] = sopt[epeakmasstol];

        cplot.copt[SpecPlot::ecluster] = copt[ecluster];
        cplot.copt[SpecPlot::eclusterms] = copt[eclusterms];
        cplot.copt[SpecPlot::eclusterscan] = copt[eclusterscan];

        if (cplot() || splot.dup2(stdin) || gnuplot_main(0, 0))
        {
          note(cerr) << cplot << endl;
          abort();
        }
      }

      shtml << "  </table>\n";
      genTableFooter(shtml);
      shtml << "  <br>\n";
      shtml << "  </div>\n";
      shtml << "  </div>\n";
      genFooter(shtml);
      shtml << "</body>\n";

      {
        cstream splot;
        ofstream shtml;
        SpecPlot cplot(cdata, cdata.cspecset, splot, shtml, cpeptide);

        cplot.sopt[SpecPlot::econtig] = ssnc;
        cplot.sopt[SpecPlot::ezoom] = "1";
        cplot.sopt[SpecPlot::efont] = sopt[efont];
        cplot.sopt[SpecPlot::eoutdir] = sopt[eoutdir];
        cplot.sopt[SpecPlot::eprefix] = "l_c_spectrum";
        cplot.sopt[SpecPlot::espectrum] = ssns;
        cplot.sopt[SpecPlot::eformat] = sopt[eformat];
        cplot.sopt[SpecPlot::epeakmasstol] = sopt[epeakmasstol];

        cplot.copt[SpecPlot::ecluster] = copt[ecluster];
        cplot.copt[SpecPlot::eclusterms] = copt[eclusterms];

        if (cplot() || splot.dup2(stdin) || gnuplot_main(0, 0))
        {
          note(cerr) << cplot << endl;
          abort();
        }
      }
    }
  }

  return 0;
}
#endif // SPECPLOT


void SpsPlot::sigext(int nsig)
{
#if 0
//#if ! defined(SPECPLOT)
  // index html
  ccron_t::iterator i = ccron.find(& typeid(& SpsPlot::genindex));
  if (i != ccron.end())
  {
    i->second.second.second[0] = "canceled";
    i->second.second.second.resize(2);

    sigexp();
  }
#endif // SPECPLOT

  if (! nsig)
    return;

  // chain
  signal(nsig, SIG_DFL);

#if defined(__linux__)
  cout << strsignal(nsig) << endl;
  note(cerr) << strsignal(nsig) << endl;
  raise(nsig);
  //pause();
  kill(0, SIGKILL);
#else
  raise(nsig);
#endif

  return;
}


void SpsPlot::sigimp(int nsig)
{
#if ! defined(SPECPLOT)
#if defined(__linux__)
  if (dbg_handler)
    dbg_handler(nsig);
#endif

  // index html
  ccron_t::iterator i = ccron.find(& typeid(& SpsPlot::genindex));
  if (i != ccron.end())
  {
    i->second.second.second[0] = "canceled";
    i->second.second.second.resize(2);

    sigexp();
  }
#endif // SPECPLOT

  if (! nsig)
    return;

  // chain
  signal(nsig, SIG_DFL);
  raise(nsig);
}


void SpsPlot::sigexp(int nsig)
{
  for (ccron_t::reverse_iterator i = SpsPlot::ccron.rbegin(), j; j = i, i != SpsPlot::ccron.rend(); ++ i)
    (i->second.first->*i->second.second.first)(i->second.second.second);
}


int SpsPlot::help(ostream & sout)
{
  sout << "Usage: spsplot [OPTION]\n";
  sout << "Options:\n";
  sout << "  --pklbin FILE               Spectrum data (pklbin format)\n";
  sout << "  --mgf FILE                  Spectrum data (mgf format)\n";
  sout << "  --pkl FILE                  Spectrum data (pkl format)\n";
  sout << "  --mzxml FILE                Spectrum data (mzxml format)\n";
  sout << '\n';
  sout << "  --aa FILE                   Amino acids file (txt format)\n";
#if ! defined(SPECPLOT)
  sout << "  --stars FILE                Star file (pklbin format)\n";
  sout << "  --comp FILE                 Component file (bin format)\n";
  sout << "  --seqs FILE                 Consensus spectra file (plkbin format)\n";
  sout << "  --index FILE                Cluster input index file (txt format)\n";
  sout << "  --cluster FILES             Cluster file (txt format)\n";
  sout << "  --clusterms FILES           Cluster MS/MS files (pklbin format)\n";
  sout << "  --clusterscan FILES         Cluster scan index file (bin format)\n";
  sout << "  --mp FILE                   Matched proteinidx file (bin format)\n";
  sout << "  --midx FILE                 Matched contig proteinidx file (pklbin format)\n";
  sout << "  --fasta FILE                Proteins database file (txt format)\n";
  sout << "  --refindex FILE             Reference input index file (txt format)\n";
  sout << "  --refmp FILE                Matched reference proteinidx file (bin format)\n";
  sout << "  --refmidx FILE              Matched reference contig proteinidx file (pklbin format)\n";
  sout << '\n';
  sout << "  --font PATH                 Font path\n";
  sout << "  --contig INDEX              Specific contig (X) or range (X:Y) 1-based\n";
  sout << "  --nodes NUMBER              Number of grid nodes\n";
  sout << '\n';
  sout << "  --spectrum INDEX            Spectrum (X) 1-based\n";
#else
  sout << "  --label FILE                Explicit labels [mass, height, type, order, charge]\n";
  sout << '\n';
  sout << "  --spectrum INDEX            Spectrum (X) 1-based\n";
  sout << "  --spectrumscan INDEX        Spectrum scan (X)\n";
  sout << "  --spectruminfo STRING       Spectrum info\n";
#endif // SPECPLOT
  sout << "  --peptide STRING            Contiguous amino acid letters\n";
  sout << "  --annot                     Write annotations [mass, height, type, order, charge]\n";
  sout << "  --shift                     Apply mass shift\n";
  sout << "  --notitle                   Disables title generation\n";
  sout << "  --zoom NUMBER               Zoom factor\n";
  sout << '\n';
  sout << "  --p FILE                    Parameters to read from file\n";
  sout << "  --format STRING             Image output format (png or eps)\n";
  sout << "  --outdir PATH               Output directory\n";
  sout << "  --prefix STRING             Output filename prefix\n";
  sout << "  --outfile FILE              Output filename (overrides prefix)\n";
  sout << "  --peakmasstol NUMBER        Amino acid peak mass tolerance\n";
  sout << "  --parentmasstol NUMBER      Parent mass tolerance\n";
  sout << "  --htmlDefs STRING           CSS and Javascript parent file path\n";
  sout << "  --pageTitle STRING          Report header\n";

#if ! defined(SPECPLOT)
  sout << '\n';
  sout << "  --stats                     Compile statistics\n";
  sout << "  --inspect FILE              Inspect file used for statistics\n";
  sout << "  --qtok                      Replace Q amino acids to K\n";
  sout << "  --ktoq                      Replace K amino acids to Q\n";
/*  sout << "  --prm FILE                  PRM scored file (pklbin format)\n";
  sout << "  --starsindex FILE           Matched pairs file (bin format)\n";
  sout << "  --pairs FILE                Filtered pairs file (bin format)\n";
  sout << "  --matchpairs FILE           Filtered pairs match file (pklbin format)\n";
  sout << "  --starspairs FILE           Matched peaks file (bin format)\n";*/
#endif // SPECPLOT

  sout << '\n';
  sout << "  --help                      Display this help and exit\n";
  sout << "  --version                   Output version information and exit\n";

  return 0;
}


int SpsPlot::version(ostream & sout)
{
  sout << PROGRAM_NAME << endl;
#if ! defined(SPECPLOT)
  sout << "spsplot 1.0." << XSTR(SPS_VERSION) << endl << endl;
#else
  sout << "specplot 1.0." << XSTR(SPS_VERSION) << endl << endl;
#endif

  sout << COPYRIGHT1 << endl;
  sout << COPYRIGHT2 << endl;
  sout << endl;

  return 0;
}


int SpsPlot::error(const string & a)
{
  cout << PROGRAM_NAME << endl;
#if ! defined(SPECPLOT)
  cout << "spsplot 1.0." << XSTR(SPS_VERSION) << endl << endl;
#else
  cout << "specplot 1.0." << XSTR(SPS_VERSION) << endl << endl;
#endif

  cout << COPYRIGHT1 << endl;
  cout << COPYRIGHT2 << endl;
  cout << endl;

  cerr << "Invalid or missing option " << a << endl;

  cerr << "Type 'spsplot --help' for more information." << endl;

  return -1;
}
