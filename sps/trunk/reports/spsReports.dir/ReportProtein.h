///////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_PROTEIN_H__
#define __REPORT_PROTEIN_H__
///////////////////////////////////////////////////////////////////////////////

#include "ReportBase.h"

////////////////////////////////////////////////////////////////////////////////
namespace spsReports {

///////////////////////////////////////////////////////////////////////////////
  /*! \brief Builds protein pages
   */
class ReportProtein : public ReportBase {


 public:

  // Constructors and destructor
  /*! \brief Constructor
  This constructor, which does not recieve a contig table filename, generates the protein list page.
   */
  ReportProtein(const string &projectPath, const string &proteinTableFilename);

  /*! \brief Constructor
  This constructor, which recieves a contig table filename, generates the single protein page. The page
  contains protein specific information, plus the list of contigs that map to it below.
   */
  ReportProtein(const string &projectPath, const string &proteinTableFilename, const string &contigTableFilename);

  // destructor
  /*! \brief Destructor
   */
  virtual ~ReportProtein() {};

};

};
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
