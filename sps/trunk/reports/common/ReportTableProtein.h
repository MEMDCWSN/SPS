////////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_TABLE_PROTEIN_H__
#define __REPORT_TABLE_PROTEIN_H__
////////////////////////////////////////////////////////////////////////////////
#include <string>
#include <vector>

#include "ReportTableBase.h"

////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
#define TABLE_PROTEIN_FILTER_COL_PROTEIN  0
////////////////////////////////////////////////////////////////////////////////
// Defines a view for proteins (default view)
////////////////////////////////////////
//
// Table colTypes has the following structure:
//
// colTypes[0] -> (CTstring) Protein name in fasta file
//   --> text = Protein name = "<1>"
//   --> columnLabel = "Protein"
//   --> link = "/cgi-bin/spsplot --"
// colTypes[1] -> (CTstring) Protein description
//   --> text = protein description in fasta file = "<2>"
//   --> columnLabel = "Description"
// colTypes[2] -> (CTstring) contigs
//   --> text = number of assotiated contigs = "<3>"
//   --> columnLabel = "Contigs"
// colTypes[3] -> (CTstring) spectra
//   --> text = number of associated spectra = "<4>"
//   --> columnLabel = "Spectra"
// colTypes[4] -> (CTstring) Amino acids
//   --> text = "<5>"
//   --> columnLabel = "Amino acids"
// colTypes[5] -> (CTstring) coverage %
//   --> text = "<6>"
//   --> columnLabel = "Coverage (%)"
//
////////////////////////////////////////////////////////////////////////////////
// Generated table, per row:
//
// cells[row][0] -> text   --> index of protein
// cells[row][1] -> text   --> protein name
// cells[row][2] -> text   --> protein description
// cells[row][3] -> text   --> number of contigs
// cells[row][4] -> text   --> number of spectra
// cells[row][5] -> text   --> amino acids
// cells[row][6] -> text   --> coverage %
//
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
 /*! \brief Report table class for Proteins

   Defines Proteins views in reports, based on proteins tables

   */
class ReportTableProtein : public ReportTableBase {

 protected:


 public:

  //////////////////////////////////////////////////////////////////////////////
  // Constructors and destructor

  //! \name CONSTRUCTORS
  //@{

  /*! \brief Contructor. Takes a project path and the table file name as parameters
   */
  ReportTableProtein(const string &projectPath, const string &tableFilename, int columnFilter = TABLE_FILTER_COL_NONE);

  /*! \brief Default destructor
   */
  virtual ~ReportTableProtein() {};


  //////////////////////////////////////////////////////////////////////////////
  // Methods to build views

  //default view
  /*! \brief
   */
  virtual void defineView(void);
  // use a different view
  /*! \brief
   */
  virtual void defineView2(void);

};
////////////////////////////////////////////////////////////////////////////////
}; // namespace
////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////
