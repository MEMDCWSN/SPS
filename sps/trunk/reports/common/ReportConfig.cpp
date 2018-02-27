///////////////////////////////////////////////////////////////////////////////
#include "ReportConfig.h"

#include <fstream>

#include "ReportBase64.h"

///////////////////////////////////////////////////////////////////////////////

namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
// File name composition
string ReportConfig::composeFileName(const string &projectDir, const string &fileName)
{
  // Compose output path
  string aux = projectDir;
  if(aux[aux.length()-1] != '/')
    aux += '/';
  // add filename
  aux += fileName;
  // return composed filename
  return aux;
}
///////////////////////////////////////////////////////////////////////////////
void ReportConfig::clearPwd(void)
{
  m_passwd = "";
  m_pwd = false;
}
///////////////////////////////////////////////////////////////////////////////
void ReportConfig::setPwd(const string &p)
{
  // checks if password size is 0. If so, the password is cleared
  if(p.length() == 0) {
    clearPwd();
    return;
  }

  // uu64Encodes the password. This is a simple encoding schema.
  string encoded = base64_encode(reinterpret_cast<const unsigned char*>(p.c_str()), p.length());
  m_passwd = encoded;
  // validate password usage
  m_pwd = true;
}
///////////////////////////////////////////////////////////////////////////////
int ReportConfig::readFile(void)
{
  // string for holding a single line of text read
  string line;

  // open the file
  ifstream f(m_fname.c_str());

  // check for file open errors.
  if(!f)
    return ERROR;

  // check of file was open
  if (f.is_open()) {

    // Check if the file has readable contents at this point
    if( f.good() ) {
      // read a single line
      getline (f, line);
      // as we only have on field, it is the password flag
      m_pwd = (line == "1" ? true : false);
    }

    // if the file is still valid, read the password (the next line)
    if( f.good() ) {
      getline (f, line);
      m_passwd = line;
    }

    // nothing more to read ; close the file
    f.close();
  }

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
int ReportConfig::writeFile(void)
{
  // open the file for writing
  ofstream f(m_fname.c_str());

  // if the file is not open...
  if(!f)
    return ERROR;

  // dump the password flag and the password to the file
  f << m_pwd << endl;
  f << m_passwd << endl;

  // close the file
  f.close();

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
};
///////////////////////////////////////////////////////////////////////////////
