////////////////////////////////////////////////////////////////////////////////
#include "ReportRendererHtmlDynamic.h"
#include "ReportDefines.h"
#include "copyright.h"

////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
int ReportRendererHtmlDynamic::generateReport(ReportGeneratorData &data)
{
  //////////////////////////////////////////////////////////////////////////////
  // data needed to generate tables

  string              empty;
  vector<string>      aux;

  setTablesDir(data.tablesDir);
  // Specify project directory (needed for input datat)
  setProjectDir(data.projectDir);
  // Specify server cgi-bin location
  setServerLocation(data.server);
  // Specify cells per line
  setCellsPerLine(data.cellsPerLine);
  // specify fasta filename
  setFastaFilename(data.filenameProteins);
  // specify if clusters layer is present
  setClusterLayerState(data.noClusters);
  // set the MS/MS filename
  setMsFilename(data.filenameConsensusSpectra);
  // specify the realign flag
  setAllowRealign(data.allowRealign);
  // specify the dynamic flag
  setDynamic(data.dynamic);
  // Specify AAs file
  setAasFile(data.aasFile);
  setAasFileDir(data.aasFileDir);
  setAasFilename();

  //set password
  setPwd(data.pwd.length());

  // Specify project directory (needed to get data to generate reports dynamically)
  if(data.targetProjectDir.length())
    setProjectDirRel(data.targetProjectDir);
  else
    setProjectDirRel(data.projectDir);


  //////////////////////////////////////////////////////////////////////////////
  // header report generation

  // define report cluster object, and load tables
  ReportHeader rh(data.tablesDir, m_tableNameHeader);
  // load tables
  rh.load();
  // define file name for HTML report page
  string fn = "index.html";
  // add specified target directory
  fn = composeFileName(data.outDir, fn);
  // open file to write to
  ofstream ofh(fn.c_str(), ios::out | ios::binary);
  // build the report page
  renderReport(&rh, ofh, empty);
  // close output file
  ofh.close();

  //////////////////////////////////////////////////////////////////////////////
  // End

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Page prolog and epilogue
////////////////////////////////////////////////////////////////////////////////
// prolog
int ReportRendererHtmlDynamic::renderProlog(ReportBase *table, ostream &outstream)
{
  // document begin
  outstream << "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Strict//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd\">";
  outstream << "<html xmlns=\"http://www.w3.org/1999/xhtml\" xml:lang=\"en\" lang=\"en\">";
  // head section
  outstream << "<head>";
  outstream << "<meta http-equiv=\"content-type\" content=\"text/html; charset=iso-8859-1\" />";
  outstream << "<meta http-equiv=\"Content-Language\" content=\"en-us\" />";
  // title and icon
  outstream << "  <title>" << REPORT_TITLE << "</title>\n";
  outstream << "<link rel=\"shortcut icon\" href=\"images/favicon.ico\" type=\"image/icon\" />";

  // spsReports code
  outstream << "<script type=\"text/javascript\" language=\"javascript\" src=\"js/bst.js\"></script>";
  outstream << "<script type=\"text/javascript\" language=\"javascript\" src=\"js/spsReport.js\"></script>";
  outstream << "<script type=\"text/javascript\" language=\"javascript\" src=\"js/XHConn.js\"></script>";

  // main styles
  outstream << "<link href=\"styles/main.css\" rel=\"stylesheet\" type=\"text/css\" />";

  // lightbox code
  outstream << "<script type='text/javascript' language=\"javascript\" src='js/prototype.js'></script>";
  outstream << "<script type='text/javascript' language=\"javascript\" src='js/scriptaculous.js?load=effects,builder'></script>";
  outstream << "<script type='text/javascript' language=\"javascript\" src='js/lightbox.js'></script>";
  outstream << "<link rel='stylesheet' href='css/lightbox.css' type='text/css' media='screen' />";

  // jQuery
	outstream << "<script src='js/jquery-1.8.2.min.js' type='text/javascript'></script>";
	outstream << "<script src='js/jquery.cookie.js' type='text/javascript'></script>";
	outstream << "<script src='js/jquery.contextMenu.js' type='text/javascript'></script>";
  outstream << "<script src='js/jquery.base64.js' type='text/javascript'></script>";
	outstream << "<link href='js/jquery.contextMenu.css' rel='stylesheet' type='text/css' />";
  outstream << "<script>jQuery.noConflict();</script>";

  // pwd form
	outstream << "<script src='js/pwd.js' type='text/javascript'></script>";
  outstream << "<link rel='stylesheet' href='css/pwd.css' type='text/css' media='screen' />";

  // end of head section
  outstream << "</head>";

  //body section
  outstream << "<body onload='javascript:init();'> ";

  // containers
  outstream << "<div class='wrapper' align='center'>";
  // logo
  //outstream << "<div><h3 align=center><img src='data:image/jpg;base64," << ICON_LOGO << "' /></h3></div><br /><br />";
  outstream << "<div><h3 align=center><img src='images/logo.jpg' /></h3></div><br /><br />";
}
////////////////////////////////////////////////////////////////////////////////
// epilog
int ReportRendererHtmlDynamic::renderEpilog(ReportBase *table, ostream &outstream)
{
  // contanier section end
  outstream << "</div></div>";

  // page footers
  outstream << "<div class='push'></div>";
  outstream << "</div>";

  outstream << "<div class='footer'><table width='100%'><tr><td class='VHSep'></td></tr><tr><td class='HSep'></td><td class='ln'></td><td class='HSep'></td></tr><tr><td class='HSep'></td><td class='Footer'>";
  outstream << REPORT_FOOTER_TEXT;
  outstream << "</td></tr></table></div>";

  // body end
  outstream << "</body>\n";

  // the box for the navigation icons
  addNavBox(outstream);

  // Add context menus
  addContextMenu(outstream);

  // add password form
  addPwdForm(outstream);

  // html end
  outstream << "</html>\n";

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
// Table rendering exceptions
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
int ReportRendererHtmlDynamic::renderTableExceptionMainPage(ReportTableBase *table, ostream &outstream)
{
  // declare and initialize table row interator, to allow direct access to the table
  TableIterator ti = table->begin();

  if(ti == table->end())
    return ERROR;

  // get table data
  string job          = (*(*ti))[0];
  string user         = (*(*ti))[1];
  string status       = (*(*ti))[2];
  string elapsed      = (*(*ti))[3];
  string log          = (*(*ti))[4];
  string contig       = (*(*ti))[5];
  string protein      = (*(*ti))[6];
  string cluster      = (*(*ti))[7];
  string fileStr      = (*(*ti))[8];
  string indicesStr   = (*(*ti))[9];

  vector<string> files, indices;
  // split indices and filenames into strings
  stringSplit2(fileStr, files, "|");
  stringSplit2(indicesStr, indices, "|");

  // output tables

  outstream << "<div id='bodyWrapper'>";
  outstream << "<div id='textWrapper'>";

  outstream << "<input type='hidden' id='projectDir' value='" << m_projectDirRel << "' />";
  outstream << "<input type='hidden' id='tablesDir' value='" << m_tablesDir << "' />";
  outstream << "<input type='hidden' id='cellsPerLine' value='" << m_cellPerLine << "' />";
  outstream << "<input type='hidden' id='serverLocation' value='" << m_serverLocation << "' />";
  outstream << "<input type='hidden' id='fastaFilename' value='" << m_fastaFilename << "' />";
  outstream << "<input type='hidden' id='noClusters' value='" << m_noClusters << "' />";
  outstream << "<input type='hidden' id='contigRefs' value='" << m_pwd << "' />";
  outstream << "<input type='hidden' id='msFilename' value='" << m_msFilename << "' />";
  outstream << "<input type='hidden' id='allowRealign' value='" << m_allowRealign << "' />";
  outstream << "<input type='hidden' id='dynamic' value='" << m_dynamic << "' />";
  outstream << "<input type='hidden' id='aas_file' value='" << m_aasFilename << "' />";

  outstream << "<div id='mainDiv2'>";


  outstream << "<table align='center'><tr><td></td><td width='990px'>";
  outstream << "<table class='mainform'>";
  outstream << "<tr>";
  outstream << "<th colspan='0'>Job Status</th>";
  outstream << "</tr>";
  outstream << "<tr><td>";

  outstream << "<table class='sched ' width='100%'>";

  outstream << "  <tr>";
  outstream << "    <th width='25%' bgcolor='#003399'><span style='color:white'>Job</span></th>";
  outstream << "    <td>" << job <<"</td>";
  outstream << "  </tr>";

  outstream << "  <tr>";
  outstream << "    <th width='25%' bgcolor='#003399'><span style='color:white'>User</span></th>";
  outstream << "    <td>" << user << "</td>";
  outstream << "  </tr>";

  outstream << "  <tr>";
  outstream << "    <th width='25%' bgcolor='#003399'><span style='color:white'>Status</span></th>";
  outstream << "    <td style='background-color:#ffffff;' id='status'></td>";
  outstream << "  </tr>";

  //outstream << "  <tr>";
  //outstream << "    <th width='25%' bgcolor='#003399'><span style='color:white'>Elapsed</span></th>";
  //outstream << "    <td>" << elapsed << "</td>";
  //outstream << "  </tr>";

//  outstream << "  <tr>";
//  outstream << "    <th width='25%' bgcolor='#003399'><span style='color:white'>Log</span></th>";
//  outstream << "    <td><a href='spsplot.log'>" << log << "</a></td>";
//  outstream << "  </tr>";

  outstream << "  <tr>";
  outstream << "    <th width='25%' bgcolor='#003399' rowspan='2'><span style='color:white'>Data</span></th>";
  outstream << "    <td><a href='#' onclick='javascript:TablesAll.loadPage(" << PAGE_CONTIGS << ",0);'>Group by Contig</a></td>";
  outstream << "  </tr>";

  outstream << "  <tr>";
  outstream << "    <td><a href='#' onclick='javascript:TablesAll.loadPage(" << PAGE_PROTEINS << ",0);'>Group by Protein</a></td>";
  outstream << "  </tr>";

  outstream << "  <tr>";
  outstream << "    <th width='25%' bgcolor='#003399'><span style='color:white'>Cluster Data</span></th>";
  outstream << "    <td>";
  //outstream << "      <a href='" << cluster << ".txt'>All Clusters (txt)</a><br>";

  for(int i = 0 ; i < files.size() && i < indices.size() ; i++) {
    int idx = getInt(indices[i].c_str());
    outstream << "      <a  href='#' onclick='javascript:TablesAll.loadPage(" << PAGE_SPECTRA << "," << idx << ");'>Group by <i>" << files[i] << "</i></a><br />";
  }
  outstream << "    </td>";
  outstream << "  </tr>";

  outstream << "  </table></td></tr>";

  outstream << "  <tr>";
  outstream << "    <td colspan='0' class='bottomline'>&nbsp;</td>";
  outstream << "  </tr>";

  outstream << "</table>";

  outstream << "</td><td></td></tr></table>";

  outstream << "</div>";

  outstream << "<div id='mainDiv' style='display: none'></div>";

  //outstream << "</div></div>";
}
///////////////////////////////////////////////////////////////////////////////
}; //namespace
////////////////////////////////////////////////////////////////////////////////
