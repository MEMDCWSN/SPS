////////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_RENDERER_HTML_DYNAMIC_H__
#define __REPORT_RENDERER_HTML_DYNAMIC_H__
////////////////////////////////////////////////////////////////////////////////
#include "ReportRendererHtml.h"

////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
// ReportRendererHtml
//
// According to each columnType, a specific HTML sequence is output to output stream.
//
// While building column headers
//
//
// <td>
//   <columnLabel>
// </td>
//
//
// -----------------------------------------------------------------------------
// While building table
//
// ReportColumnTypeString:
// <td class="<cssClass>">
//   <text> or
//   <input type="button" value="<text>" onClick="<onClick>" /> or
//   <input ID="<ID>" type="text" style="text-transform: uppercase; width:100%" />
// </td>
//
//
// ReportColumnTypeImageOnDemand:
// <td class="<cssClass>">
//   <a href="<>" onclick="<onClick>">
//       <label>  or <img src='<icon>'>
//   </a>
// </td>
//
//
// ReportColumnTypeSequencesBox:
// <td class="<cssClass>">
//   <table>
//                --- begining of repeat block
//     <tr>
//       call to ReportColumnType
//     </tr>
//                --- repeat per ReportColumnType entry on vector
//   </table>
// </td>
//
////////////////////////////////////////////////////////////////////////////////
  /*! \brief Class to render an HTML dynamic report
   */
class ReportRendererHtmlDynamic : public ReportRendererHtml {

 protected:


  //////////////////////////////////////////////////////////////////////////////
  // Report page prolog and epilogue

  // prolog
  /*! \brief Renders page prolog
   */
  virtual int renderProlog(ReportBase *table, ostream &outstream);
  // epilog
  /*! \brief Renders page epilog
   */
  virtual int renderEpilog(ReportBase *table, ostream &outstream);

  //////////////////////////////////////////////////////////////////////////////
  // Table rendering exceptions

  // main page
  /*! \brief Table expections rendering method
   */
  virtual int renderTableExceptionMainPage(ReportTableBase *table, ostream &outstream);


 public:


  // Constructors and destructor
  /*! \brief Constructors
   */
  ReportRendererHtmlDynamic()  {};
  /*! \brief destructor
   */
  ~ReportRendererHtmlDynamic() {};


  /*! \brief Report generation entry point
   */
  virtual int generateReport(ReportGeneratorData &data);


};
///////////////////////////////////////////////////////////////////////////////
}; //namespace
////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////
