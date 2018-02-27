///////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_CONFIG_H__
#define __REPORT_CONFIG_H__
///////////////////////////////////////////////////////////////////////////////
#include <string>

#include "ReportDefines.h"

///////////////////////////////////////////////////////////////////////////////
using namespace std;
///////////////////////////////////////////////////////////////////////////////
namespace spsReports {
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
 /*! \brief ReportConfig table class

   Defines a report configuration table class.
   Currently, holds password for dynamic reports.

   */
class ReportConfig {

 protected:


   /*! \brief Compose the full file name
   */
  string composeFileName(const string &projectDir, const string &fileName);

 public:

  // Password section

   /*! \brief Password usage flag
   */
  bool    m_pwd;

   /*! \brief Password contents
   */
  string  m_passwd;

   /*! \brief Table file name
   */
  string  m_fname;

  // Constructors and destructor
  ReportConfig(string &d, string &f) : m_pwd(false) {m_fname = composeFileName(d,f);};
  ~ReportConfig() {};

   /*! \brief Clears the password

     Clears the password

     */
  void clearPwd(void);

   /*! \brief Sets the password

     Sets the report password

    @param pwd The new password
     */
  void setPwd(const string &pwd);

   /*! \brief Sets the password

     Sets the report password

    @param pwd The new password
     */
  int readFile(void);

   /*! \brief Read the file

     Reeads the configuration file
     */
  int writeFile(void);

};
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
