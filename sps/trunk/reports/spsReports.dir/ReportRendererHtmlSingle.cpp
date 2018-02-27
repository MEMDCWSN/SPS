////////////////////////////////////////////////////////////////////////////////

#include "ReportRendererHtmlSingle.h"
#include "Tokenizer.h"
#include "copyright.h"
#include "Timer.h"
#include "ReportDefines.h"

////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
int ReportRendererHtmlSingle::generateReport(ReportGeneratorData &data)
{
  //////////////////////////////////////////////////////////////////////////////
  // data needed to generate tables

  // Specify executables directory (needed for specplot module)
  setExeDir(data.exeDir);
  // specify display level
  setDisplayLevel(data.displayLevel);
  // tables location
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
  // Specfiy AAs file
  setAasFile(data.aasFile);
  setAasFileDir(data.aasFileDir);
  setAasFilename();
  // Specify project directory (needed to get data to generate reports dynamically)
  if(data.targetProjectDir.length())
    setProjectDirRel(data.targetProjectDir);
  else
    setProjectDirRel(data.projectDir);

  //set password
  setPwd(data.pwd.length());

  //////////////////////////////////////////////////////////////////////////////
  // file generation

  // define file name for HTML report page
  string fn = "index.html";
  // add specified target directory
  fn = composeFileName(data.outDir, fn);
  // open file to write to
  ofstream of(fn.c_str(), ios::out | ios::binary);


  //////////////////////////////////////////////////////////////////////////////
  // report generation

  // render prolog
  renderProlog(of);

  // render the code to hold the images
  of << "<script language='javascript'>";

  // Proteins table
  ReportTableProtein t4(data.tablesDir, m_tableNameProtein);
  // load the table
  t4.loadTable();
  // write the table
  writeTable(&t4, of, 0);

  ReportTableProteinCoverage t5(data.tablesDir, m_tableNameProteinCoverage);
  // load the table
  t5.loadTable();
  // write the table
  writeTable(&t5, of, 1);

  // Define contig list view, not filtered, for images
  ReportTableContig t1(data.tablesDir, m_tableNameContig);
  // Define view for generating images
  t1.defineViewImages();
  // no borders and no header
  t1.noBorders(); t1.noHeaders();
  // load the table
  t1.loadTable();
  // build the table
  buildTable(&t1, -1, -1);
  // render the table
  renderTable(&t1, of);
  // write the table
  writeTable(&t1, of, 2);

  // Define cluster list view, not filtered, for images
  ReportTableClusterConsensus t2(data.tablesDir, m_tableNameCluster);
  // Define view for generating images
  t2.defineViewImages();
  // no borders and no header
  t2.noBorders(); t2.noHeaders();
  // load the table
  t2.loadTable();
  // build the table
  buildTable(&t2, -1, -1);
  // render the table
  renderTable(&t2, of);
  // write the table
  writeTable(&t2, of, 3);

  // Define spectra list view, not filtered, for images
  ReportTableInputSpectra t3(data.tablesDir, m_tableNameSpectra);
  // Define view for generating images
  t3.defineViewImages();
  // no borders and no header
  t3.noBorders(); t3.noHeaders();
  // load the table
  t3.loadTable();
  // build the table
  buildTable(&t3, -1, -1);
  // render the table
  renderTable(&t3, of);
  // write the table
  writeTable(&t3, of, 4);


  // end of image load section
  of << "</script>";

  // render prolog
  renderProlog2(of);

  // The index page table
  ReportTableBase *t0 = new ReportTableHeader(data.tablesDir, m_tableNameHeader);
  // load tables
  t0->loadTable();
  // render main page
  renderTableExceptionMainPage(t0, of);

  // render epilog
  renderEpilog(of);


  // close output file
  of.close();


  //////////////////////////////////////////////////////////////////////////////
  // End

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Page prolog and epilogue
////////////////////////////////////////////////////////////////////////////////
// prolog
int ReportRendererHtmlSingle::renderProlog(ostream &outstream)
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


  // language initializer section
  outstream << "<script language='javascript'>";
  outstream << "var I=[];";
  outstream << "var R=0;";
  outstream << "var T=[];";
  outstream << "</script>";

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
int ReportRendererHtmlSingle::renderProlog2(ostream &outstream)
{
  // end of head section
  outstream << "</head>";

  //body section
  outstream << "<body onload='javascript:init();'> ";

  // containers
  outstream << "<div class='wrapper' align='center'>";
  // logo
  //outstream << "<div><h3 align=center><img src='data:image/jpg;base64," << ICON_LOGO << "' /></h3></div><br /><br />";
  outstream << "<div><h3 align=center><img src='images/logo.jpg' /></h3></div><br /><br />";

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// epilog
int ReportRendererHtmlSingle::renderEpilog(ostream &outstream)
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
// main page
int ReportRendererHtmlSingle::renderTableExceptionMainPage(ReportTableBase *table, ostream &outstream)
{
  // declare and initialize table row interator, to allow direct access to the table
  TableIterator ti = table->begin();

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
  outstream << "<input type='hidden' id='aas_file' value='" << m_aasFilename << "' />";
  outstream << "<input type='hidden' id='dynamic' value='false' />";
  outstream << "<input type='hidden' id='tableLoad' value='false' />";
  outstream << "<input type='hidden' id='tableReload' value='false' />";
  outstream << "<input type='hidden' id='serverUpdate' value='false' />";
  outstream << "<input type='hidden' id='serverImages' value='false' />";

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
// Table rendering methods
///////////////////////////////////////////////////////////////////////////////
int ReportRendererHtmlSingle::renderTable(ReportTableBase *table, ostream &outstream)
{
  // render contents
  return ReportRendererBase::renderTable(table, outstream);
}
////////////////////////////////////////////////////////////////////////////////
// build a row comon method. Cycles thu all columnType cells and buildCell for each ColumnType item
int ReportRendererHtmlSingle::renderTableRow(vector<string> &row, ostream &outstream)
{
  return ReportRendererBase::renderTableRow(row, outstream);
}
////////////////////////////////////////////////////////////////////////////////
// Table Building methods
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// render cell comon method. Renders a specific cell based on ColumnType specifications and a row of data
int ReportRendererHtmlSingle::buildTableCell(ReportColumnTypeBase *base, vector<string> *row, stringstream &ss)
{
  return ReportRendererBase::buildTableCell(base, row, ss);
}
///////////////////////////////////////////////////////////////////////////////
// render image
int ReportRendererHtmlSingle::buildCellImageOnDemand(ReportColumnTypeImageOnDemand *ct, vector<string> *row, stringstream &ss)
{
  // auxiliary variables
  stringstream icon, label, url;

  //  cout << ct->iconDisplayLevel << " < " << m_displayLevel << " -- " << ct->iconRenderer << endl;
  string id = parseTemplates(ct->id, row);
  string tag  = IMAGE_ICON_ID_PREFIX; tag += id;
  string tag3 = IMAGE_LARGE_ID_PREFIX; tag3 += id;

  // Icon path/image
  if(ct->iconParams.size()) {
    // parse the needed files
    parseFilesVector(row, ct->files);
    // get parsed parameters / URL
    vector<string> pars; // = parseTemplatesAll(ct->iconParams, row);
    parseParamsVector(pars, ct->iconParams, row);
    // if there is a renderer, use it to render the image
    if(ct->iconRenderer.size()  && ct->iconDisplayLevel < m_displayLevel) {
      // render the image
      renderImage(ct->iconRenderer, pars, ct->files);

      stringstream aux;
      for(int i = 0 ; i < pars.size() ; i++)
        aux << pars[i] << ' ';;

      // output the image
      ss << "I[R]=[];";
      ss << "I[R][0]='';"; //='"   << ct->iconRenderer << "';";
      ss << "I[R][1]='';"; //='"   << aux.str() << "';";
      ss << "I[R][2]='"   << tag <<"';";
      ss << "I[R][3]='';";
      ss << "I[R++][4]='" << m_image << "';";
    }
  }

  // URL template to be used to get the image
  if(ct->renderer.size() && ct->linkDisplayLevel < m_displayLevel) {
    // parse the needed files
    parseFilesVector(row, ct->files);
    // parse the params vector
    vector<string> pars;
    parseParamsVector(pars, ct->params, row);
    // render the image
    renderImage(ct->renderer, pars, ct->files);

    stringstream aux;
    for(int i = 0 ; i < pars.size() ; i++)
      aux << pars[i] << ' ';;

    // output the image
    ss << "I[R]=[];";
    ss << "I[R][0]='';"; //='"   << ct->renderer << "';";
    ss << "I[R][1]='';"; //='"   << aux.str() << "';";
    ss << "I[R][2]='"   << tag3 <<"';";
    ss << "I[R][3]='';";
    ss << "I[R++][4]='" << m_image << "';";
  }

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
int ReportRendererHtmlSingle::writeTable(ReportTableBase *table, ostream &outstream, int idx)
{
  outstream << "T[" << idx << "]=\"";
  table->writeTable(outstream, FILE_SEPARATOR, '%');
  outstream << "\";";
}
///////////////////////////////////////////////////////////////////////////////
}; //namespace
////////////////////////////////////////////////////////////////////////////////
