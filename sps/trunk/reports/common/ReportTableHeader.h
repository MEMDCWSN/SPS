////////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_TABLE_HEADER_H__
#define __REPORT_TABLE_HEADER_H__
////////////////////////////////////////////////////////////////////////////////
// includes
#include <string>
#include <vector>

#include "spectrum.h"
#include "ReportTableBase.h"

////////////////////////////////////////////////////////////////////////////////
// defines
#define TABLE_INPUT_HEADER_FILENAME "tableHeader.txt"

#define FILTER_COLUMN 0

namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
// View for header page
////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////
// Generated table, per column:
////////////////////////////////////////
//
// cells[0][0] -> text   --> Table index
// cells[1][0] -> text   --> index in specset
// cells[2][0] -> text   --> scan #
// cells[3][0] -> text   --> cluster index (used when filtered by cluster consensus)
// cells[4][0] -> text   --> Protein index
// cells[5][0] -> text   --> Protein name
// cells[row][6] -> text   --> spectrm file name
// cells[row][7] -> text   --> Reference sequence --> generateSequence()
// cells[row][8] -> text   --> Homolog sequence --> generateSequence()
// cells[row][9] -> text   --> DeNovo sequence --> generateSequence()
// cells[row][10] -> text  --> User sequence
// cells[row][11] -> text  --> mass value, from specs[i].parentMass
// cells[row][12] -> text  --> charge value, from specs[i].parentCharge
// cells[row][13] -> text  --> B% 
// cells[row][14] -> text  --> Y&
// cells[row][15] -> text  --> BY Int %
// cells[row][16] -> text  --> input file extension: mgf or mzxml
//
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
 /*! \brief Report Header class for Proteins

   Defines Header class, which holds data to build the initial report page
   
   */
class ReportTableHeader : public ReportTableBase {

 protected:


 public:
  
  // Constructors and destructor

  //! \name CONSTRUCTORS
  //@{

  /*! \brief Contructor. Takes a project path and the table file name as parameters
   */
  ReportTableHeader(const string &projectPath, const string &tableFilename);

    //@}

    //! \name DESTRUCTOR
    //@{

  /*! \brief Default destructor
   */
  ~ReportTableHeader() {};
  
  //default view
  /*! \brief Default view, which consists of the initial report table
   */
  virtual void defineView(void);
  
};
////////////////////////////////////////////////////////////////////////////////
}; //namespace
////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////
