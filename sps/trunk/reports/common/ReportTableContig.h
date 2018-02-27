///////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_TABLE_CONTIG_H__
#define __REPORT_TABLE_CONTIG_H__
///////////////////////////////////////////////////////////////////////////////
#include <string>
#include <vector>

#include "ReportTableBase.h"

///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
// Specific defines
#define TABLE_CONTIG_FILTER_COL_CONTIG    0
#define TABLE_CONTIG_FILTER_COL_PROTEIN   1

///////////////////////////////////////////////////////////////////////////////
// Contig views
//
///////////////////////////////////////////////////////////////////////////////
// View for contig list
////////////////////////////////////////
//
// Table colTypes has the following structure:
//
// colTypes[0] -> (CTstring) Contig index
//   --> text = index in specset = "<0>"
//   --> columnLabel = "Index"
// colTypes[1] -> (CTstring) Number of spectra
//   --> text = number of spectra for this contig = "<1>"
//   --> columnLabel = "Spectra"
// colTypes[2] -> (CTIOD) Contig image
//   --> icon  --> "/cgi/contplot --projectdir <projectdir> --p <paramsfile> --contig <0> --reference <2> --homolog <3> --consensus <5> --user <5> --zoom .4 --target cout --spectra<8>"
//   --> columnLabel = "Contig"
//   --> url = "/cgi/spsplot --projectdir <projectdir> --p <paramsfile> --table contig --contig <0> --output cout"
//   --> label = ""
// colTypes[2] -> (CTseqsBox)
//   --> columnLabel = "Sequences"
//     sequences[0]->(CTIOD) : Reference
//       --> icon  --> ""
//       --> label --> Reference sequence = "<2>"
//       --> url   --> ""
//     sequences[1]->(CTIOD) : Homolog
//       --> icon  --> ""
//       --> label --> Homolog sequence = "<3>"
//       --> url   --> ""
//     sequences[2]->(CTIOD) : consensus
//       --> icon  --> ""
//       --> label --> consensus sequence = "<4>"
//       --> url   --> ""
//     sequences[3]->(CTIOD) : User
//       --> icon  --> ""
//       --> label --> User sequence = "<5>"
//       --> url   --> ""
//       --> id    --> 'user_<row>_<col>'
//     sequences[4]-> isInput  = true
//       --> text = ""
//       --> id = 'user_<row>_<col>'
//     sequences[5]-> isButton = true
//       --> text = "Update"
//       --> onClick="javascript:DoOnCick('input_<row>_<col>', <0>, 'user_<row>_<col>')";
//      #--> url --> --update <0> --pklbin <1> --peptide <6>
// colTypes[3] -> (CTstring) protein
//       --> text = protein name & description = "<6><nl><7>"  # <nl> means 'newline'
//       --> columnLabel = "Protein"
//
///////////////////////////////////////////////////////////////////////////////
// View for contig image with edit box
////////////////////////////////////////
//
// Table colTypes has the following structure:
// colTypes[0] -> (CTseqsBox)
//   --> columnLabel = ""
//     sequences[0]->(CTIOD) : Reference
//       --> icon  --> "/cgi/contplot --projectdir <projectdir> --contig <0> --reference <1> --homolog <2> --consensus <3> --user <4> --zoom 1 --output-format uu64 --output cout"
//       --> columnLabel = ""
//       --> url = ""
//       --> label = ""
//     sequences[1]-> isInput  = true
//       --> text = ""
//       --> id    --> 'user_<row>_<col>'
//     sequences[2]-> isButton = true
//       --> text = "Update"
//       --> onClick='javascript:DoOnCick('input_<row>_<col>', <0>, user_<row>_<col>);';
//       --> url --> --update --contig <0> --peptide <4>
//
////////////////////////////////////////////////////////////////////////////////
// Generated table, per row:
//
// cells[row][0] -> text   --> Contig index (1-based)
// cells[row][1] -> text   --> Protein index (1-based)
// cells[row][2] -> text   --> Number of spectra
// cells[row][3] -> text   --> Reference sequence --> generateSequence(i, return)
// cells[row][4] -> text   --> Homolog sequence   --> generateSequence(i, return)
// cells[row][5] -> text   --> DeNov sequence
// cells[row][6] -> text   --> User sequence
// cells[row][7] -> text   --> Protein name
// cells[row][8] -> text   --> Protein description
// cells[row][9] -> text   --> spectra (coma separated indexes)
// cells[row][10] -> text  --> Reference intervals
// cells[row][11] -> text  --> Homolog intervals
// cells[row][12] -> text  --> Reference offset
// cells[row][13] -> text  --> Homolog offset
// cells[row][14] -> text  --> reverse flag
// cells[row][15] -> text  --> file names needed: 15 abruijn
// cells[row][16] -> text  --> file names needed: 16 stars
// cells[row][17] -> text  --> file names needed: 17 sps_seqs
// cells[row][18] -> text  --> tool used
// cells[row][19] -> text  --> Grouped homolog sequence
// cells[row][20] -> text  --> Grouped reference sequence
//
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
 /*! \brief Report table class for Contigs

   Defines Contigs views in reports, based on contigs tables

   */
class ReportTableContig : public ReportTableBase {

 protected:

 public:

  //////////////////////////////////////////////////////////////////////////////
  // Constructors and destructor

  //! \name CONSTRUCTORS
  //@{

  /*! \brief Contructor. Takes a project path and the table file name as parameters
   */
  ReportTableContig(const string &projectPath, const string &tableFilename, int columnFilter = TABLE_FILTER_COL_NONE);

    //@}

    //! \name DESTRUCTOR
    //@{

  /*! \brief Default destructor
   */
  virtual ~ReportTableContig() {};

  /*! \brief Element find method, given a string as input
   */
  virtual int find(string &data);
  /*! \brief Element find method, given an integer as input
   */
  virtual int find(int data);

  // IDs[0] = proteinID, IDs[1] = contigID
  /*! \brief Returns the composite ID string vector for the table
   */
  virtual void getId(vector<string> &IDs) {getFilteredId(IDs, 1);ReportTableBase::getId(IDs);};

  //////////////////////////////////////////////////////////////////////////////
  // Methods to build views

  //default view
  /*! \brief Defines the default view for the clusters table, which consists of the clusters table
   */
  virtual void defineView(void);
  // use a different view
  /*! \brief Defines the alternative view, which consist of a single large contig image
   */
  virtual void defineView2(void);
  // view to generate images
  /*! \brief Defines the view used to generate images for the single file reports
   */
  virtual void defineViewImages(void);

};
////////////////////////////////////////////////////////////////////////////////
}; // namespace
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
