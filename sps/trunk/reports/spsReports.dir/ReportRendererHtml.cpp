////////////////////////////////////////////////////////////////////////////////

#include "ReportRendererHtml.h"
#include "Tokenizer.h"
#include "copyright.h"
#include "Timer.h"
#include "ReportDefines.h"
#include "Specific.h"

////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////
int ReportRendererHtml::generateReport(ReportGeneratorData &data)
{
  // performance timers
  Timer_c timer, timer2;

  // set current number of child processes
  m_nprocess = 0;
  // set number of available CPUs
  m_ncpu = data.cpu;
  // reset the file numbering
  m_fileNumber = 0;
  //////////////////////////////////////////////////////////////////////////////
  // data needed to generate tables

  // Specify executables directory (needed for specplot module)
  setExeDir(data.exeDir);
  // Specify project directory (needed for input data and report output)
  setProjectDir(data.projectDir);
  // Specify AAs file
  setAasFile(data.aasFile);
  setAasFileDir(data.aasFileDir);
  setAasFilename();

  // specify display level
  setDisplayLevel(data.displayLevel);

  //////////////////////////////////////////////////////////////////////////////
  // header report generation

  DEBUG_MSG("---- Generating header page ----");
  generateReportPageHeader(data);
  if(m_ncpu == 1) DEBUG_MSG("Took " << timer.restart());

  //////////////////////////////////////////////////////////////////////////////
  // proteins report generation

  DEBUG_MSG("---- Generating proteins page ----");
  generateReportPageProteins(data);
  if(m_ncpu == 1) DEBUG_MSG("Took " << timer.restart());


  //////////////////////////////////////////////////////////////////////////////
  // Proteins main page

  DEBUG_MSG("---- Generating protein pages ----");
  generateReportPageProtein(data);
  if(m_ncpu == 1) DEBUG_MSG("Took " << timer.restart());

  //////////////////////////////////////////////////////////////////////////////
  // Protein coverage main page

  DEBUG_MSG("---- Generating protein coverage pages ----");
  generateReportPageProteinCoverage(data);
  if(m_ncpu == 1) DEBUG_MSG("Took " << timer.restart());


  //////////////////////////////////////////////////////////////////////////////
  // Protein coverage main page

  DEBUG_MSG("---- Generating protein coverage CSV text files ----");
  generateReportPageProteinCoverageCSV(data);
  if(m_ncpu == 1) DEBUG_MSG("Took " << timer.restart());

  //////////////////////////////////////////////////////////////////////////////
  // contigs report generation

  DEBUG_MSG("---- Generating contigs page ----");
  generateReportPageContigs(data);
  if(m_ncpu == 1) DEBUG_MSG("Took " << timer.restart());


  //////////////////////////////////////////////////////////////////////////////
  // contig pages report generation

  DEBUG_MSG("---- Generating contig pages ----");
  generateReportPageContig(data);
  if(m_ncpu == 1) DEBUG_MSG("Took " << timer.restart());


  //////////////////////////////////////////////////////////////////////////////
  // cluster report generation

  DEBUG_MSG("---- Generating cluster pages ----");
  generateReportPageCluster(data);
  if(m_ncpu == 1) DEBUG_MSG("Took " << timer.restart());


  //////////////////////////////////////////////////////////////////////////////
  // spectra report generation

  DEBUG_MSG("---- Generating spectra pages ----");
  generateReportPageSpectra(data);
  if(m_ncpu == 1) DEBUG_MSG("Took " << timer.restart());


  //////////////////////////////////////////////////////////////////////////////
  // wait for processes to finish

  if(m_ncpu > 1 && forkable())
    for (; m_nprocess; -- m_nprocess)
      waitProcess();

  //////////////////////////////////////////////////////////////////////////////
  // Timer

  //DEBUG_MSG("---- Done! ----");
  //DEBUG_MSG(timer2.stop());


  //////////////////////////////////////////////////////////////////////////////
  // End

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
int ReportRendererHtml::generateReportPageHeader(ReportGeneratorData &data)
{
  string empty;

  clearNavigationBar();
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
  renderReport(&rh, ofh, empty, true);
  // close output file
  ofh.close();
}
////////////////////////////////////////////////////////////////////////////////
int ReportRendererHtml::generateReportPageProteins(ReportGeneratorData &data)
{
  // define report protein list object, and load tables
  ReportProtein pl(data.tablesDir, m_tableNameProtein);
  // load table
  pl.load();
  // generate pages with pagination
  paginate(data, pl, -1, "proteins.", "proteins", 10);
}
////////////////////////////////////////////////////////////////////////////////
int ReportRendererHtml::generateReportPageProteinCoverage(ReportGeneratorData &data)
{
  // clear the page split navigation bar
  clearNavigationBar();
  string empty;
  // define report protein list object, and load tables
  ReportProteinCoverage pl(data.tablesDir, m_tableNameProteinCoverage);
  // load table
  pl.load();
  // get the proteins ID column from the proteins coverage table
  vector<string> aux = pl.getTableColumn(0, TABLE_PROTEIN_COVERAGE_FILTER_COL_PROTEIN);
  // cycle tru all proteins
  for(int i = 0 ; i < aux.size() ; i++) {
    // set the filter for a specific contig in cluster table
    pl.applyFilter(aux[i]);
    // set filename
    string fnamePrefix = PROTEIN_DETAILS_NAME_PREFIX;  fnamePrefix +=  aux[i]; fnamePrefix += ".html";
    // generate pages with pagination
    string fn = composeFileName(data.outDir, fnamePrefix);
    // open file to write to
    ofstream of(fn.c_str(), ios::out | ios::binary);
    // build THE navigation bar
    stringstream navBarTop;
    buildNavigationBar(pl, navBarTop, 12);
    string nbt = navBarTop.str();
    // build a page
    renderReport(&pl, of, nbt, true);
    //renderTableExceptionProteinCoverage(&pl, of);
    // close output file
    of.close();
  }
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
int ReportRendererHtml::generateReportPageProteinCoverageCSV(ReportGeneratorData &data)
{
  string empty;
  // define report protein list object, and load tables
  ReportProteinCoverageCsv pl(data.tablesDir, m_tableNameProteinCoverage);
  // load table
  pl.load();
  // get the proteins ID column from the proteins coverage table
  vector<string> aux = pl.getTableColumn(0, TABLE_PROTEIN_COVERAGE_FILTER_COL_PROTEIN);
  // cycle tru all proteins
  for(int i = 0 ; i < aux.size() ; i++) {
    // set the filter for a specific contig in cluster table
    pl.applyFilter(aux[i]);
    // set filename
    string fnamePrefix = PROTEIN_DETAILS_NAME_PREFIX;  fnamePrefix +=  aux[i]; fnamePrefix += ".txt";
    // generate pages with pagination
    string fn = composeFileName(data.outDir, fnamePrefix);
    // open file to write to
    ofstream of(fn.c_str(), ios::out | ios::binary);
    // build a page
    renderReport(&pl, of, empty, false);
    // close output file
    of.close();
  }
  return OK;
}////////////////////////////////////////////////////////////////////////////////
int ReportRendererHtml::generateReportPageProtein(ReportGeneratorData &data)
{
  // define specific protein report object
  ReportProtein p2(data.tablesDir, m_tableNameProtein, m_tableNameContig);
  // load table
  p2.load();
  // get the proteins column from the proteins table
  vector<string> aux = p2.getTableColumn(0, TABLE_PROTEIN_FILTER_COL_PROTEIN);
  // cycle tru all proteins
  for(int i = 0 ; i < aux.size() ; i++) {
    // set the filter for a specific contig in cluster table
    p2.applyFilter(aux[i]);
    // set filename
    string fnamePrefix = "protein.";  fnamePrefix +=  aux[i]; fnamePrefix += ".";
    // generate pages with pagination
    paginate(data, p2, DEFAULT_ROWS_PER_TABLE, fnamePrefix, "protein", 11);
  }
}
////////////////////////////////////////////////////////////////////////////////
int ReportRendererHtml::generateReportPageContigs(ReportGeneratorData &data)
{
  // define report object for contig list
  ReportContig cl(data.tablesDir, m_tableNameContig);
  // load table
  cl.load();
  // generate pages with pagination
  paginate(data, cl, DEFAULT_ROWS_PER_TABLE, "contigs.", "contigs", 20);
}
////////////////////////////////////////////////////////////////////////////////
int ReportRendererHtml::generateReportPageContig(ReportGeneratorData &data)
{
  // Child (specific) contig report pages
  clearNavigationBar();
  // define specific contig report object
  ReportContig rc(data.tablesDir, m_tableNameContig, m_tableNameCluster, m_tableNameSpectra, data.noClusters);
  // load table
  rc.load();
  // get the contig id column from the contigs table
  vector<string> aux = rc.getTableColumn(0, TABLE_CONTIG_FILTER_COL_CONTIG);
  // cycle tru all contigs
  for(int i = 0 ; i < aux.size() ; i++) {
    // set the filter for a specific contig in cluster/spectra table
    rc.applyFilter(aux[i]);
    // set filename
    string fnamePrefix = "contig.";  fnamePrefix +=  aux[i]; fnamePrefix += ".";
    // generate pages with pagination
    paginate(data, rc, DEFAULT_ROWS_PER_TABLE, fnamePrefix, "contig", 21);
  }
}
////////////////////////////////////////////////////////////////////////////////
int ReportRendererHtml::generateReportPageCluster(ReportGeneratorData &data)
{
  if(data.noClusters)
    return OK;
  // cluster pages are generated depending on flag
  if(!data.noClusters) {
    clearNavigationBar();
    // define report cluster object, and load tables
    ReportCluster clust(data.tablesDir, m_tableNameCluster, m_tableNameSpectra);
    // load tables
    clust.load();
    // get the cluster id column from the cluster table
    vector<string> aux = clust.getTableColumn(0, TABLE_CLUSTER_FILTER_COL_CLUSTER);
    // cycle tru all clusters
    for(int i = 0 ; i < aux.size() ; i++) {
      // set the filter for a specific contig in cluster table
      clust.applyFilter(aux[i]);
      // set filename
      string fnamePrefix = "cluster.";  fnamePrefix +=  aux[i]; fnamePrefix += ".";
      // generate pages with pagination
      paginate(data, clust, DEFAULT_ROWS_PER_TABLE, fnamePrefix, "cluster", 31);
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
int ReportRendererHtml::generateReportPageSpectra(ReportGeneratorData &data)
{
  clearNavigationBar();
  // define report cluster object, and load tables
  ReportInputSpectra is(data.tablesDir, m_tableNameSpectra, data.noClusters);
  // load tables
  is.load();
  // get the file index column from the input spectra table
  vector<string> aux = is.getTableColumn(0, TABLE_SPECTRA_FILTER_COL_FILE);
  // sort elements for duplicate removal
  sort(aux.begin(), aux.end(), stringSortCompare);
  // duplicate find
  vector<string>::iterator it;
  // using default comparison
  it = unique (aux.begin(), aux.end(), stringUniqueCompare);
  // remove extra items
  aux.resize( it - aux.begin() );
  // cycle tru all input files
  for(int i = 0 ; i < aux.size() ; i++) {
    // set the filter for a specific contig in cluster table
    is.applyFilter(aux[i]);
    // set filename
    string fnamePrefix = "spectra.";  fnamePrefix +=  aux[i]; fnamePrefix += ".";
    // generate pages with pagination
    paginate(data, is, DEFAULT_ROWS_PER_TABLE, fnamePrefix, "spectra", 40);
  }
}
////////////////////////////////////////////////////////////////////////////////
int ReportRendererHtml::addBarLink(stringstream &out, string &img, string &text, string &link)
{
  string sep = "&nbsp;&nbsp;&nbsp;";

  out << "<a href='" << link << "'>";
  out << "<img height='24' width='24' src='data:image/png;base64," << img << "' />";
  out << "<span style='font-family:Calibri;font-size:140%;color:blue'><u>" << text << "</u></span>";
  out << "</a>";
  out << sep;
}
///////////////////////////////////////////////////////////////////////////////
int ReportRendererHtml::buildNavigationBar(ReportBase &rep, stringstream &out, int barType)
{
  string img, text, link;

  string sep = "&nbsp;&nbsp;&nbsp;";
  // initial page
  img = ICON_HOME;
  text = "Initial page";
  link = "index.html";
  addBarLink(out, img, text, link);

  // get the IDs
  vector<string> IDs;
  rep.getId(IDs);

  // protein --> add protein list
  //if(IDs.size() == 1) {
    img = ICON_PROTEIN_LIST;
    text = "Protein list";
    link = "proteins.0.html";
    addBarLink(out, img, text, link);
  //}

  img = ICON_CONTIG_LIST;
  text = "Contig list";
  link = "contigs.0.html";
  addBarLink(out, img, text, link);

  // protein coverage --> add protein list
  if(barType == 12) {
    stringstream aux;
    // link to protein
    aux <<  "protein." << IDs[0] << ".0.html";
    img = ICON_PROTEIN;
    text = "Protein";
    link = aux.str();
    addBarLink(out, img, text, link);
  }


  // contig --> add jump to protein
  if(IDs.size() > 1) {

    // In case the contig maps to a protein
    if(IDs[0].compare("-1"))  {

      stringstream aux, aux2;

      // link to protein list
      //img = ICON_PROTEIN_LIST;
      //text = "Protein list";
      //link = "proteins.0.html";
      //addBarLink(out, img, text, link);

      // link to protein
      aux <<  "protein." << IDs[0] << ".0.html";
      img = ICON_PROTEIN;
      text = "Protein";
      link = aux.str();
      addBarLink(out, img, text, link);

      // link to protein coverage
      aux2 << PROTEIN_DETAILS_NAME_PREFIX << IDs[0] << ".html";
      img = ICON_PROTEIN_COVERAGE;
      text = "Protein coverage";
      link = aux2.str();
      addBarLink(out, img, text, link);

    // in case the contig doesn't map to a protein
    } //else {

      //img = ICON_CONTIG_LIST;
      //text = "Contig list";
      //link = "contigs.0.html";
     // addBarLink(out, img, text, link);
    //}
  }

  if(IDs.size() > 2) {
    stringstream aux;
    aux << "contig." << IDs[1] << ".0.html";
    img = ICON_CONTIG;
    text = "Contig";
    link = aux.str();
    addBarLink(out, img, text, link);
  }

  if(IDs.size() > 3) {
  }
}
///////////////////////////////////////////////////////////////////////////////
int ReportRendererHtml::paginate(ReportGeneratorData &data, ReportBase &rep, int ipp, string fnamePrefix, string pageName, int barType)
{
  //Timer_c timer;

  // build THE navigation bar
  stringstream navBarTop;
  buildNavigationBar(rep, navBarTop, barType);
  string nbt = navBarTop.str();

  if(ipp == -1) {
    string fn = fnamePrefix;
    fn += "0.html";
    fn = composeFileName(data.outDir, fn);
    buildPage(rep, fn, 0, -1, nbt);
    return OK;
  }

  // get the total number of elements, plus the vector of IDs
  vector<string> IDs, IDsEnd, fNames;
  // navigation bar name and suffix
  string fnameSuffix = ".html";
  int nPages = rep.getPageSplitStats(ipp, IDs, IDsEnd, fNames, fnamePrefix, fnameSuffix);

  // Declare filename holder and navigation bar holder variables
  string navBar;
  // build the navigation bar (pagination)
  buildNavigationBar(navBar, IDs, IDsEnd, fNames);
  // set the navigation bar
  setNavigationBar(navBar);
  // cycle thru all pages
  for(int j = 0 ; j < nPages ; j++) {
    // verbose output
    //stringstream aux; aux << "Generating " << pageName << " page: ";
    //verboseOutput(cout, aux.str().c_str(), fNames[j].c_str(), "...", false);
    // add specified target directory
    string fn = composeFileName(data.outDir, fNames[j]);
    // build a page
    buildPage(rep, fn, j * ipp, ipp, nbt);
    // verbose output (time)
    //verboseOutput(cout, timer.restart());
  }
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
//int ReportRendererHtml::buildPage(ReportBase &rep, string &fn, int start, int len)
int ReportRendererHtml::buildPage(ReportBase &rep, string &fn, int start, int len, string &navBar)
{
  // if number of CPUs greater then 2, run in multiprocess mode
  if(m_ncpu > 1 && forkable()) {

    // add a new process
    m_nprocess ++;

    //before a fork, add a file numbering displacement to avoid collisions
    m_fileNumber += 1000;

    switch (myfork()) {

    case -1:
      // something went wrong
      abort();
    case 0:
      // open file to write to
      ofstream of(fn.c_str(), ios::out | ios::binary);
      // build the report page
      renderReport(&rep, of, navBar, true, start, len);
      // close output file
      of.close();
      // exit process
      _exit(0);
    }

    // parent process

    // test if process # limit has been reached
    if (m_nprocess >= m_ncpu) {
      // if so, wait for one to finish
      waitProcess();
      // and decrease the number of active processes
      m_nprocess --;
    }

    // if we are not on linux, let's keep it single threaded
    return OK;
  }

  // open file to write to
  ofstream of(fn.c_str(), ios::out | ios::binary);
  // render the page
  int ret = renderReport(&rep, of, navBar, true, start, len);
  // close output file
  of.close();
  // exit
  return ret;
}
////////////////////////////////////////////////////////////////////////////////
// Page component generation methods
////////////////////////////////////////////////////////////////////////////////
void ReportRendererHtml::buildNavigationBar(string &navBar, vector<string> &IDs, vector<string> &IDsEnd, vector<string> &fNames)
{
  navBar = "<div align='center'><p>";

  for(int i = 0 ; i < IDs.size() ; i++, navBar += "   ") {

    navBar += "<a href='";
    navBar += fNames[i];
    navBar += "'>[";
    navBar += IDs[i];
    navBar += '-';
    navBar += IDsEnd[i];
    navBar += "]</a>";
  }
  navBar += "</p></div>";
}
////////////////////////////////////////////////////////////////////////////////
// Page prolog and epilogue
////////////////////////////////////////////////////////////////////////////////
// prolog
int ReportRendererHtml::renderProlog(ReportBase *table, ostream &outstream)
{
  outstream << "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\"\n";
  outstream << "  \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n";
  outstream << "<HTML xmlns=\"http://www.w3.org/1999/xhtml\">\n";

  outstream << "<head>\n";
  outstream << "  <meta http-equiv=\"Content-Type\" content=\"text/shtml; charset=ISO-8859-1\" />\n";
  // page title
  outstream << "  <title>" << REPORT_TITLE << "</title>\n";
  // main styles and icon
  outstream << "  <link href=\"styles/main.css\" rel=\"stylesheet\" type=\"text/css\" />\n";
  outstream << "  <link rel=\"shortcut icon\" href=\"images/favicon.ico\" type=\"image/icon\" />\n";

  // light box
  outstream << "  <script type='text/javascript' src='js/prototype.js'></script>\n";
  outstream << "  <script type='text/javascript' src='js/scriptaculous.js?load=effects,builder'></script>\n";
  outstream << "  <script type='text/javascript' src='js/lightbox.js'></script>\n";
  outstream << "  <link rel='stylesheet' href='css/lightbox.css' type='text/css' media='screen' />\n";

  // sort table
  outstream << "  <script src=\"js/sorttable.js\"></script>\n";

  // end of head section
  outstream << "</head>\n";

  // body section
  outstream << "<body>\n";

  outstream << "<div class='wrapper' align='center'>";

  // Logo image
  //outstream << "<div><h3 align=center><img src='data:image/jpg;base64," << ICON_LOGO << "' /></h3></div><br /><br />";
  outstream << "<div><h3 align=center><img src='images/logo.jpg' /></h3></div><br /><br />";

  // contanier section
  outstream << "<div id='bodyWrapper'>";
  outstream << "<div id='textWrapper'>";

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// epilog
int ReportRendererHtml::addFooterStyles(ReportBase *table, ostream &outstream)
{
  // footer styles
  outstream << "<style>\n";
  outstream << "* {margin-top: 0; margin-bottom: 0;}\n";
  outstream << "html, body {height: 100%;}\n";
  outstream << ".wrapper {min-height: 100%; height: auto !important; height: 100%;margin: 0 auto -4.5em;}\n";
  outstream << ".footer, .push {height: 4.5em;}\n";
  outstream << "TD.ln {	height: 1px; background-color: #0055ff; PADDING: 0pt; MARGIN:	0pt; line-height: 1px;}\n";
  outstream << "TD.HSep 	{width: 	20pt;}\n";
  outstream << "TD.VHSep 	{height: 	20pt;}\n";
  outstream << "</style>\n";

  return 0;
}
////////////////////////////////////////////////////////////////////////////////
// epilog
int ReportRendererHtml::renderEpilog(ReportBase *table, ostream &outstream)
{
  // contanier section end
  outstream << "</div></div>";

  // footer
  outstream << "<div class='push'></div>";
  outstream << "</div>";

  outstream << "<div class='footer'><table width='100%'><tr><td class='VHSep'></td></tr><tr><td class='HSep'></td><td class='ln'></td><td class='HSep'></td></tr><tr><td class='HSep'></td><td class='Footer'>";
  outstream << REPORT_FOOTER_TEXT;
  outstream << "</td></tr></table></div>";

  outstream << "</body>\n";
  outstream << "</html>\n";

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
int ReportRendererHtml::addNavBox(ostream &outstream)
{
  // the box for the navigation icons
  outstream << "<div id='navButtons' style='position: fixed; top: 0px; left: 0px;'></div>";

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
int ReportRendererHtml::addContextMenu(ostream &outstream)
{
  // context menu
  outstream << "<ul id='myMenu' class='contextMenu'>";
	outstream << "<li class='newtab'><a href='#newtab'>Open in new tab</a></li>";
	outstream << "<li class='newwin'><a href='#newwin'>Open in new window</a></li>";
	outstream << "</ul>";

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
int ReportRendererHtml::addPwdForm(ostream &outstream)
{
  outstream << "<div><a id='login-link' href='#login-box' class='login-window'></a></div>";
  outstream << "<div id='login-box' class='login-popup'>";
  outstream << "<form id='pwdForm' method='post' class='signin' action='#'>";
  outstream << "<fieldset class='textbox' id='pwdField'>";
  outstream << "<label class='password'>";
  outstream << "<span>Password</span>";
  outstream << "<input id='password' name='password' value='' type='password' placeholder='Password'>";
  outstream << "</label>";
  outstream << "<button class='button' type='button' onclick='checkIds2();return false;'>Unlock</button>";
  outstream << "</fieldset>";
  outstream << "</form>";
  outstream << "</div>";

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// between tables
int ReportRendererHtml::renderInterTable(ReportBase *table, ostream &outstream)
{
  // output a line break in a table row
  outstream << "<tr><td><br /></td></tr>";

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
// Table rendering exceptions
///////////////////////////////////////////////////////////////////////////////
void ReportRendererHtml::breakProteinIntoChunks(vector<string> &in, int &count)
{
  count = 0;
  for(int i = 0 ; i < in.size() ; i++) {
    // string container
    string aux;
    for(int j = 0 ; j < in[i].length() ; j++) {
      // line break;
      if(count % 50 == 0) {
        aux += "<br />";
      // spacer
      } else if(count % 10 == 0) {
        aux += "&nbsp;";
      }
      // copy element
      aux += in[i][j];
      count++;
    }
    in[i] = aux;
  }
}
///////////////////////////////////////////////////////////////////////////////
void ReportRendererHtml::colorProteinString(vector<string> &in, string &out)
{
  for(int i = 0 ; i < in.size() ; i++) {
    if(in[i].size()) {
      out += "<font color='";
      if(i % 2) {
        out += "0";
      } else {
        out += "#AAAAAA";
      }
      out += "'>";
      out += in[i];
      out += "</font>";
    }
  }
}
///////////////////////////////////////////////////////////////////////////////
int ReportRendererHtml::renderTableExceptionMainPage(ReportTableBase *table, ostream &outstream)
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

  //outstream << "<div id='bodyWrapper'>";
  //outstream << "<div id='textWrapper'>";
  outstream << "<table align='center'><tr><td></td><td width='990px'>";
  outstream << "<table class='mainform'>";
  outstream << "<tr> ";
  outstream << "    <th colspan='0'>Job Status</th>";
  outstream << "  </tr>";
  outstream << "  <tr><td>";
  outstream << "<table class='sched' width='100%'>";
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
  outstream << "    <td style='background-color:lightgreen;'>" << status << "</td>";
  outstream << "  </tr>";
  //outstream << "  <tr>";
  //outstream << "    <th width='25%' bgcolor='#003399'><span style='color:white'>Elapsed</span></th>";
  //outstream << "    <td>" << elapsed << "</td>";
  //outstream << "  </tr>";
  //outstream << "  <tr>";
  //outstream << "    <th width='25%' bgcolor='#003399'><span style='color:white'>Log</span></th>";
  //outstream << "    <td><a href='spsplot.log'>" << log << "</a></td>";
  //outstream << "  </tr>";
  outstream << "  <tr>";
  outstream << "    <th width='25%' bgcolor='#003399' rowspan='2'><span style='color:white'>Data</span></th>";
  outstream << "    <td><a href='" << contig << ".0.html'>Group by Contig</a></td>";
  outstream << "  </tr>";
  outstream << "  <tr>";
  outstream << "    <td><a href='" << protein << ".0.html'>Group by Protein</a></td>";
  outstream << "  </tr>";
  outstream << "  <tr>";
  outstream << "    <th width='25%' bgcolor='#003399'><span style='color:white'>Cluster Data</span></th>";
  outstream << "    <td>";
  //outstream << "      <a href='" << cluster << ".txt'>All Clusters (txt)</a><br />";

  for(int i = 0 ; i < files.size() && i < indices.size() ; i++)
    outstream << "      <a href='spectra." << indices[i] << ".0.html'>Group by <i>" << files[i] << "</i></a><br />";

  outstream << "    </td>";
  outstream << "  </tr>";
  outstream << "  </td></tr></table></td></tr>";
  outstream << "  <tr>";
  outstream << "    <td colspan='0' class='bottomline'>&nbsp;</td>";
  outstream << "  </tr>";
  outstream << "</table>";
  outstream << "</td><td></td></tr></table>";
  //outstream << "</div></div>";
}
///////////////////////////////////////////////////////////////////////////////
int ReportRendererHtml::renderTableExceptionProteinHeader(ReportTableBase *table, ostream &outstream)
{
  // declare and initialize table row interator, to allow direct access to the table
  TableIterator ti = table->begin();

  if(ti == table->end())
    return ERROR;

  // get the protein name, at indice 1
  string protienId    = (*(*ti))[0];
  string proteinName  = (*(*ti))[1];
  string proteinDesc  = (*(*ti))[2];
  string contigs      = (*(*ti))[3];
  string spectra      = (*(*ti))[4];
  string aas          = (*(*ti))[5];
  string coverage     = (*(*ti))[6];
  string sequence     = (*(*ti))[7];

  // remove | from protein name
  for(int i = 0 ; i < proteinName.size() ; i++)
    if(proteinName[i] == '|')
      proteinName[i] = ' ';

  // format protein sequence
  int count;
  vector<string> breaked;
  string coloredProtein;
  stringSplit2(sequence, breaked, "|");
  breakProteinIntoChunks(breaked, count);
  colorProteinString(breaked, coloredProtein);

  // build legend for protein sequence
  string legend;
  int current = 1;
  while(current <= count) {
    legend += "<br />";
    legend += parseInt(current);
    current += 50;
  }

  // protein header HTML code
  outstream << "<table class='result' width='90%' style='border-spacing: 0px;' align='center'>";
  outstream << "<tr>";
  outstream << "<td colspan='0' width='50%'><h2><i>" << proteinName << "</i></h2>";
  outstream << "<hr><b>" << contigs << " contigs, " << spectra << " spectra, " << aas << " amino acids, " << coverage << " coverage" << "</b></td>";
  outstream << "<td></td><td></td>";
  outstream << "</tr>";
  outstream << "<tr> ";
  outstream << "<td>" << proteinDesc << "</td>";
  outstream << "<td align='right'><tt>" << legend << "</tt></td>";
  outstream << "<td><tt>" << coloredProtein << "</tt></td>";
  outstream << "</tr>";
  outstream << "</table>";

  // link to protein detains page
  outstream << "<div align='center'>";
  outstream << "<br /><a href='" << PROTEIN_DETAILS_NAME_PREFIX << protienId << ".html'>Protein coverage</a>";
  outstream << "</div>";
  outstream << "<br />";

    // return status OK
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Table Protein Coverage
////////////////////////////////////////////////////////////////////////////////
// Generated table, per row:
//
// cells[row][0] -> text          --> Protein ID
// cells[row][1] -> text          --> Protein name
// cells[row][2] -> text          --> Protein length (AAs)
// cells[row][3] -> text list     --> Protein sequence, separated by |
// cells[row][4] -> Contig data   --> CSPS Contigs, separated by |
// cells[row][5] -> Contig data   --> SPS Contigs, separated by |
//      Contig data: : items separated by &
//        0 -> Contig ID
//        1 -> Contig name
//        2 -> Contig start
//        3 -> Contig end
//        4 -> Contig items
//           Contig Item: items separated by @. Contents separated by !
//              0 -> Beginning
//              1 -> Span
//              0 -> Content
////////////////////////////////////////////////////////////////////////////////
int ReportRendererHtml::renderTableExceptionProteinCoverage(ReportTableBase *table, ostream &outstream)
{
  // declare and initialize table row interator, to allow direct access to the table
  TableIterator ti = table->begin();

  // get the protein coverage data
  string proteinId       = (*(*ti))[TABLE_COVERAGE_FIELD_ID];
  string proteinName     = (*(*ti))[TABLE_COVERAGE_FIELD_NAME];
  string proteinLength   = (*(*ti))[TABLE_COVERAGE_FIELD_SEQ_REFERENCE];
  string proteinSequence = (*(*ti))[TABLE_COVERAGE_FIELD_PROT_SEQUENCE];
  string contigCsps      = (*(*ti))[TABLE_COVERAGE_CSPS_DATA];
  string contigSps       = (*(*ti))[TABLE_COVERAGE_SPS_DATA];

  // split the protein
  vector<string> proteinData;
  stringSplit(proteinSequence, proteinData, TABLE_SEP_L1);

  // variables to hold contig data
  vector<ContigData> contigDataCsps;
  vector<ContigData> contigDataSps;

  // get CSPS contig data
  buildContigDataStructure(contigDataCsps, contigCsps);
  // get SPS contig data
  buildContigDataStructure(contigDataSps, contigSps);

  // protein header HTML code
  outstream << "<table class='result' width='100%' style='border-spacing: 0px;' align='center'>";
  outstream << "<tr>";
  outstream << "<td colspan='0'><h2><i>" + proteinName + "</i></h2><hr></td>";
  outstream << "</tr>";
  outstream << "<tr><td>&nbsp;</td></tr>";
  outstream << "</table>";

  //outstream << "<div id='bodyWrapper'>";
  //outstream << "<div id='textWrapper'>";

  // CSV protein coverage info
  outstream << "<table><tr><td><a href='" << PROTEIN_DETAILS_NAME_PREFIX << proteinId << ".txt" << "'>Protein coverage as Excel-ready format (TXT file)</a></td></tr><tr><td>&nbsp;</td></tr></table>";

  // general position indexer
  int i = 0;
  // protein length as an integer
  int proteinLengthi = getInt(proteinLength.c_str());

  // Keep it under protein size
  while(i < proteinLengthi) {

    // Build a map key index. This is used to maintain the contig index order when outputing them under the protein sequence
    vector<int> spsID;
    vector<int> cspsID;
    // get the csps contig indexes
    getOrder(contigDataCsps, i, CELLS_PER_LINE, cspsID);
    // get the sps contig indexes
    getOrder(contigDataSps, i, CELLS_PER_LINE, spsID);

    // if we are starting a new table, add table beggining
    outstream << "<table class=\"result2\" width=\"100%\" style=\"background-color: #CCCFFF\">\n";

    // output protein
    generateProteinSequence(outstream, i, proteinData, proteinLengthi);

    // Add CSPS contig information (if exists)
    generateOutputContig(i, proteinLengthi, cspsID, outstream, CELLS_PER_LINE, contigDataCsps, false);

    // Add SPS contig information (if exists)
    generateOutputContig(i, proteinLengthi, spsID, outstream, CELLS_PER_LINE, contigDataSps, true);

    // HTML table terminator
    outstream << "</table><br />\n";

    i += CELLS_PER_LINE;
  }

  //page  = "</div></div>";
  //page += "</body>";

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
int ReportRendererHtml::renderTableExceptionProteinCoverageCSV(ReportTableBase *table, ostream &outstream)
{
  // declare and initialize table row interator, to allow direct access to the table
  TableIterator ti = table->begin();

  // get the protein coverage data
  string proteinId       = (*(*ti))[TABLE_COVERAGE_FIELD_ID];
  string proteinName     = (*(*ti))[TABLE_COVERAGE_FIELD_NAME];
  string proteinLength   = (*(*ti))[TABLE_COVERAGE_FIELD_SEQ_REFERENCE];
  string proteinSequence = (*(*ti))[TABLE_COVERAGE_FIELD_PROT_SEQUENCE];
  string contigCsps      = (*(*ti))[TABLE_COVERAGE_CSPS_DATA];
  string contigSps       = (*(*ti))[TABLE_COVERAGE_SPS_DATA];


  vector<string> proteinData;
  stringSplit(proteinSequence, proteinData, TABLE_SEP_L1);

  // variables to hold contig data
  vector<ContigData> contigDataCsps;
  vector<ContigData> contigDataSps;

  // get CSPS contig data
  buildContigDataStructure(contigDataCsps, contigCsps);
  // get SPS contig data
  buildContigDataStructure(contigDataSps, contigSps);

  // general position indexer
  int i = 0;
  // protein length as an integer
  int proteinLengthi = getInt(proteinLength.c_str());

  // Keep it under protein size
  while(i < proteinLengthi) {

    // Build a map key index. This is used to maintain the contig index order when outputing them under the protein sequence
    vector<int> spsID;
    vector<int> cspsID;
    // get the csps contig indexes
    getOrder(contigDataCsps, i, CELLS_PER_LINE, cspsID);
    // get the sps contig indexes
    getOrder(contigDataSps, i, CELLS_PER_LINE, spsID);

    // output protein
    generateProteinSequenceCSV(outstream, i, proteinData, proteinLengthi);

    // Add CSPS contig information (if exists)
    generateOutputContigCSV(i, proteinLengthi, cspsID, outstream, CELLS_PER_LINE, contigDataCsps, false);

    // Add SPS contig information (if exists)
    generateOutputContigCSV(i, proteinLengthi, spsID, outstream, CELLS_PER_LINE, contigDataSps, true);

    i += CELLS_PER_LINE;
  }
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
void ReportRendererHtml::generateProteinSequence(ostream &outstream, int i, vector<string> &proteinData, int proteinLength)
{
  outstream << "<tr>";
  outstream << "<td class='rc3' align='center'>";
  outstream << parseInt(i+1);
  outstream << "&nbsp;</td>";

  for(int j = i ; ( j < i + CELLS_PER_LINE ) && ( j < proteinLength ) ; j++) {

    outstream << "<td align='center'";
    // if an empty cell, it's a separator column. It should be 1 pixel wide
    if( proteinData[j].length() == 0 )
      outstream << "class='rh2' ";
    else {
      outstream << "class='rh1' ";
    }
    outstream << ">";
    // The AA from the protein sequence
    outstream <<  proteinData[j];
    // header cell terminator
    outstream << "</td>";
  }
  // Header row terminator
  outstream << "</tr>\n";
}
////////////////////////////////////////////////////////////////////////////////
void ReportRendererHtml::generateProteinSequenceCSV(ostream &outstream, int i, vector<string> &proteinData, int proteinLength)
{
  // if we are starting a new table, add table beggining
  outstream << '\n';
  outstream << parseInt(i+1);
  outstream << CSV_SEP;

  for(int j = i ; ( j < i + CELLS_PER_LINE ) && ( j < proteinLength ) ; j++)
    outstream << CSV_SEP << proteinData[j] << CSV_SEP;

  // Header row terminator
  outstream << "\n";
}
////////////////////////////////////////////////////////////////////////////////
void ReportRendererHtml::generateOutputContig(int i, int proteinSize, vector<int> &vectorID, ostream &outstream, int cellPerLine, vector<ContigData> &contig, bool link)
{
  // Add contig information (if exists)
  for(int j = 0 ; j < vectorID.size() ; j++)  {
    // get the contig sequence info
    const ContigData &contigSequence = contig[vectorID[j]];
    // Check if sequence in the range we are outputing now
    if( (contigSequence.start < i + cellPerLine) &&
        (contigSequence.end > i              )    ) {

      // Write the contig id and link
      if(link) {
        outstream << "<tr><th align='right'><a href='contig.";
        outstream << getIntFromSeqName(contigSequence.name);
        outstream << ".0.html' style='color: white'>";
        outstream << contigSequence.name;
        outstream << "</a></th>";

        //outstream << "<tr><th align=\"right\">";
        //outstream << contigSequence.name;
        //outstream << "</th>";
      } else {
        outstream << "<tr>";
        outstream << "<td class='rc4'>";
        outstream << "CSPS ";
        //outstream << contigSequence.name;
        outstream << parseInt(contigSequence.id + 1);
        outstream << "&nbsp;";
        outstream << "</TD>";
      }

      // find first cell to output
      int l = 0;
      while( (l < contigSequence.contigMatch.size()) &&
             (i > contigSequence.contigMatch[l].start + contigSequence.contigMatch[l].colSpan) )
        l++;

      // cycle thru
      for(int k = i ; (k < i + cellPerLine) && (k < proteinSize) ; k++) {

        // if start position is lower than current position, output an empty cell
        if(k < contigSequence.start)
          outstream << "<td class='rc2' />\n";
        else if( (l >= contigSequence.contigMatch.size()) )
          outstream << "<td class='rc2' />\n";
        // otherwise, the content
        else {

          outstream << "<td ";
          // outstream << " class=\"rc1\" style=\"background-color: transparent; border: solid 0 #060; border-left-width:1px;border-right-width:1px;\" ";

          int border = 0;

          string outputString = contigSequence.contigMatch[l].data;

          // Calc colspan nr of occupied cells
          int colspan = contigSequence.contigMatch[l].colSpan + 1;
          // careful with split cells at the beggining or end -- end
          if(contigSequence.contigMatch[l].start + contigSequence.contigMatch[l].colSpan >= i + cellPerLine) {
            colspan -= contigSequence.contigMatch[l].start + contigSequence.contigMatch[l].colSpan - i - cellPerLine + 1;
            border += 2;
            if(colspan < (contigSequence.contigMatch[l].colSpan + 1) / 2 )
              outputString = "";
          }
          // beggining
          if(contigSequence.contigMatch[l].start < i) {
            colspan -= i - contigSequence.contigMatch[l].start;
            border++;
            if(colspan <= (contigSequence.contigMatch[l].colSpan + 1) / 2 )
              outputString = "";
          }

          outstream << " class='rca";
          outstream << parseInt(border);
          outstream << "' ";

          if(colspan > 1) {
            outstream << " colspan='";
            outstream << parseInt(colspan);
            outstream << "'";
            k += colspan-1;
          }
          outstream << '>';
          outstream <<  outputString;
          outstream << "</td>";
          l++;
        }
      }
      outstream << "</tr>\n";
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
void ReportRendererHtml::generateOutputContigCSV(int i, int proteinSize, vector<int> &vectorID, ostream &outstream, int cellPerLine, vector<ContigData> &contig, bool link)
{
  // Add contig information (if exists)
  for(int j = 0 ; j < vectorID.size() ; j++)  {
    // get the contig sequence info
    const ContigData &contigSequence = contig[vectorID[j]];
    // Check if sequence in the range we are outputing now
    if( (contigSequence.start < i + cellPerLine) &&
        (contigSequence.end > i              )    ) {

      // Write the contig id and link
      if(link) {
        outstream << contigSequence.name;
      } else {
        outstream << "CSPS " << parseInt(contigSequence.id + 1);
      }
      // separator
      outstream << CSV_SEP;

      // find first cell to output
      int l = 0;
      while( (l < contigSequence.contigMatch.size()) &&
             (i > contigSequence.contigMatch[l].start + contigSequence.contigMatch[l].colSpan) )
        l++;

      // cycle thru
      for(int k = i ; (k < i + cellPerLine) && (k < proteinSize) ; k++) {

        // if start position is lower than current position, output an empty cell
        if(k < contigSequence.start) {
          outstream << CSV_SEP;
          outstream << CSV_SEP;
        } else if( (l >= contigSequence.contigMatch.size()) ) {
          outstream << CSV_SEP;
          outstream << CSV_SEP;
        // otherwise, the content
        } else {

          int border = 0;

          string outputString = contigSequence.contigMatch[l].data;

          // Calc colspan nr of occupied cells
          int colspan = contigSequence.contigMatch[l].colSpan + 1;
          // careful with split cells at the beggining or end -- end
          if(contigSequence.contigMatch[l].start + contigSequence.contigMatch[l].colSpan >= i + cellPerLine) {
            colspan -= contigSequence.contigMatch[l].start + contigSequence.contigMatch[l].colSpan - i - cellPerLine + 1;
            border += 2;
            if(colspan < (contigSequence.contigMatch[l].colSpan + 1) / 2 )
              outputString = "";
          }
          // beggining
          if(contigSequence.contigMatch[l].start < i) {
            colspan -= i - contigSequence.contigMatch[l].start;
            border++;
            if(colspan <= (contigSequence.contigMatch[l].colSpan + 1) / 2 )
              outputString = "";
          }

          if(!(border & 0x01))
            outstream << '|';

          if(colspan == 1) {
            // cell content
            outstream << CSV_SEP;
            outstream << outputString;
          } else {

            if((outputString.length() > 0) && (outputString[0] == '(')) {
              // outputing (xx, yy), and empty cells following
              outstream << CSV_SEP;
              outstream << outputString;
              // until the end of cell
              while(--colspan && (k < i+cellPerLine)) {
                outstream << CSV_SEP;
                outstream << CSV_SEP;
                k++;
              }

            } else {
              // outputing A.B.C.D
              int cells = 0;
              while( (cells < colspan-1) && (k < i+cellPerLine)) {
                outstream << CSV_SEP;
                outstream << outputString[cells];
                outstream << CSV_SEP;
                outstream << '.';
                k++;
                cells++;
              }
              outstream << CSV_SEP;
              outstream << outputString[cells];

            }

          }
          // cell end
          outstream << CSV_SEP;
          // right border content (separator)
          if( (!(border & 0x02)) && (true) && (true) ) ;
//            outstream << '|';

          l++;
        }
      }
      outstream << '\n';
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
void ReportRendererHtml::getOrder(vector<ContigData> &contig, int i, int size, vector<int> &order)
{
  vector<ContigOrdering> aux;
  // Make sure ordes is empty
  order.clear();
  // cycle thru all contigs
  for(int j = 0 ; j < contig.size() ; j++ ) {
    // Check for empty contigs
    if(contig[j].contigMatch.size() == 0) continue;
    // Check if sequence in the range we are outputing now
    if( (contig[j].start < i + size) &&
        (contig[j].end > i              )    ) {
      //order.push_back(j);
      // create ordering cell with contig info
      ContigOrdering orderingCell;
      orderingCell.contigIndex  = j;
      orderingCell.startIndex   = contig[j].start;
      orderingCell.endIndex     = contig[j].end;
      aux.push_back(orderingCell);
    }
  }

  // Order
  //sort(order.begin(), order.end());

  // Order
  sort(aux.begin(), aux.end());
  // set order in order vector
  for(int k = 0 ; k < aux.size() ; k++)
    order.push_back(aux[k].contigIndex);
}

////////////////////////////////////////////////////////////////////////////////
string ReportRendererHtml::getIntFromSeqName(const string &seq)
{
  vector<string> aux;
  stringSplit(seq, aux, ":");
  if(aux.size() > 1)
    return aux[1];
  return "";
}
////////////////////////////////////////////////////////////////////////////////
void ReportRendererHtml::buildContigDataStructure(vector<ContigData> &contigDataArray, string &contig)
{
  // split by contig
  vector<string> aux;
  stringSplit(contig, aux, TABLE_SEP_L1);
  // cycle thru all contigs
  for(int i = 0 ; i < aux.size() ; i++) {
    // array to hold data for one contig
    ContigData contigData;
    // get contig header info + elements
    vector<string> aux2;
    stringSplit(aux[i], aux2, TABLE_SEP_L2);
    // add header data
    // ID
    contigData.id = getInt(aux2[0].c_str());
    // Name
    contigData.name = aux2[1];
    // start
    contigData.start = getInt(aux2[2].c_str());
    // end
    contigData.end = getInt(aux2[3].c_str());
    // array for contig elements
    //ContigMatch contigMatch;
    // get elements
    if(aux2.size() > 4) {
      vector<string> elems;
      stringSplit(aux2[4], elems, TABLE_SEP_L3);
      // cycle thru elems
      for(int j = 0 ; j < elems.size() ; j++) {
        // elemet storage space
        contigMatchElem matchElem;
        // split elements
        vector<string> elem;
        stringSplit(elems[j], elem, TABLE_SEP_L4);
        // copy start position
        matchElem.start   = getInt(elem[0].c_str());
        // copy colspan
        matchElem.colSpan = getInt(elem[1].c_str());
        // copy element data
        matchElem.data    = elem[2];
        // store the element
        contigData.contigMatch.push_back(matchElem);
      }
    }
    contigDataArray.push_back(contigData);
  }
}
////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////
// Table rendering methods
///////////////////////////////////////////////////////////////////////////////
  // render table comon method. Cycles thru all the table rows and invokes buildTableRow() method to render a row
int ReportRendererHtml::renderTable(ReportTableBase *table, ostream &outstream)
{
  switch(table->getRenderingException()) {

  case RENDERING_EXCEPTION_NONE:
    break;

  case RENDERING_EXCEPTION_PROTEIN_HEADER:
    return renderTableExceptionProteinHeader(table, outstream);
    break;

  case RENDERING_EXCEPTION_MAIN_PAGE:
    return renderTableExceptionMainPage(table, outstream);

  case RENDERING_EXCEPTION_PROTEIN_COVERAGE_PAGE:
    return renderTableExceptionProteinCoverage(table, outstream);

  case RENDERING_EXCEPTION_PROTEIN_COVERAGE_CSV_PAGE:
    return renderTableExceptionProteinCoverageCSV(table, outstream);

  default:
    break;
  }


  // HTML class definition
  string aux = "class='result sortable'";
  // boder definition
  if(!table->doBorders())
    aux = "border='0'";
  // HTML table start
  outstream << "<table " << aux << " align='center'>";
  // render contents
  ReportRendererBase::renderTable(table, outstream);
  // HTML table terminator
  outstream << "</table>";

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// render a header row method. Renders a table header row.
int ReportRendererHtml::renderTableHeaderRow(vector<string> &row, ostream &outstream)
{
  // render the header
  outstream << "<tr>";
  ReportRendererBase::renderTableHeaderRow(row, outstream);
  outstream << "</tr>";

  // return status OK
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// build a row comon method. Cycles thu all columnType cells and buildCell for each ColumnType item
int ReportRendererHtml::renderTableRow(vector<string> &row, ostream &outstream)
{
  outstream << "<tr align='center' valign='middle'>";
  ReportRendererBase::renderTableRow(row, outstream);
  outstream << "</tr>";

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Table Building methods
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// build header cell comon method. Renders all ColumnTypes
int ReportRendererHtml::buildTableHeaderCell(ReportColumnTypeBase *base, stringstream &ss)
{
  // cell begin HTML tag with class
  ss << "<th><span style='color:white'>";
  // cell content
  ss << base->columnLabel;
  // cell end HTML tag
  ss << "</span></th>";

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
// render cell comon method. Renders a specific cell based on ColumnType specifications and a row of data
int ReportRendererHtml::buildTableCell(ReportColumnTypeBase *base, vector<string> *row, stringstream &ss)
{
  // don't render cells with dynamic content
  if(base->dynamic)
    return OK;

  // check if this cell is a StringMultiple. If it is, we shouldn't issue <a> tags here.
  ReportColumnTypeStringMultiple *auxM = dynamic_cast<ReportColumnTypeStringMultiple*>(base);

  // auxiliary variables
  stringstream cls, link;

  // gather needed attributes
  if(base->link.size() && (auxM == NULL))
    link << parseTemplates(base->link, row); // parseTemplatesAll

  if(base->cssClass.size())
    cls << " class='" << base->cssClass << "'";


  // cell begin HTML tag with class
  ss << "<td align='center' valign='middle'" << cls.str() << ">";
  // Link section
  if(base->link.size()  && (auxM == NULL))
    ss << "<a href='" << link.str() << "'>";
  // process base class cell renderer
  ReportRendererBase::buildTableCell(base, row, ss);
  // Link section
  if(base->link.size()  && (auxM == NULL))
    ss << "</a>";
  // cell end HTML tag
  ss << "</td>";

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
// Renderers for the diferent cell types
int ReportRendererHtml::buildCellImageOnDemand(ReportColumnTypeImageOnDemand *ct, vector<string> *row, stringstream &ss)
{
  // auxiliary variables
  stringstream icon, label, url;

  //  cout << ct->iconDisplayLevel << " < " << m_displayLevel << " -- " << ct->iconRenderer << endl;

  // Icon path/image
  if(ct->iconParams.size()) {

    // result holder
    string image64;
    // parse the needed files
    parseFilesVector(row, ct->files);
    // get parsed parameters / URL
    vector<string> pars; // = parseTemplatesAll(ct->iconParams, row);
    parseParamsVector(pars, ct->iconParams, row);
    // if there is a renderer, use it to render the image
    if(ct->iconRenderer.size()  && ct->iconDisplayLevel < m_displayLevel) {
      //renderImage(ct->iconRenderer, parseTemplatesAll(pars, row));
      renderImage(ct->iconRenderer, pars, ct->files);
      image64 = "data:image/png;base64,";
      image64 += m_image;
      ss << "<img src='" << image64 << "' />";
    } else
      ss << "<p>" << ct->alt <<"</p>";
  }


  // label to show for the link
  if(ct->label.size())
    label << parseTemplates(ct->label, row); // parseTemplatesAll
  // URL template to be used to get the image
  if(ct->renderer.size() && ct->linkDisplayLevel < m_displayLevel) {
    //url << ct->renderer << ' ' <<  parseTemplatesAll(ct->params, row);
    //renderImage(ct->renderer, parseTemplatesAll(ct->params, row));
    //outstream << "<img src='data:image/png;base64," << m_image << "' />";
    // parse the needed files
    parseFilesVector(row, ct->files);
    // parse the params vector
    vector<string> pars;
    parseParamsVector(pars, ct->params, row);
    // render the image
    renderImage(ct->renderer, pars, ct->files);

    //outstream << "<a href='" << url.str() << "'>";
    string image64 = "data:image/png;base64,";
    image64 += m_image;
    //outstream << "<a href='#' onClick='JavaScript: var a=\"" << m_image << "\";showImage(a);'>";
    ss << "<a href='" << image64 << "' rel='lightbox'>";
  }

  // button
  if(ct->label.size()) {
    string aux = label.str();
    if(ct->splitLabel)
      aux = splitSequence(aux, "<wbr>");
    ss << aux;
  }
    //outstream << "<div onClick='javascript:alert(\"#1\");var a='" << m_image << "';showImage(a);'>" << label.str() << "</div>";

  // url for IOD - terminator
  if(ct->renderer.size())
    ss << "<a>";

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
int ReportRendererHtml::buildCellString(ReportColumnTypeString *ct, vector<string> *row, stringstream &ss)
{
  // auxiliary variables
  string textProcessed;
  stringstream onclick, id, text;

  // gather needed attributes
  if(ct->text.size()) {
    text << parseTemplates(ct->text, row); // parseTemplatesAll
    textProcessed = text.str();
    if(ct->splitText)
      textProcessed = splitText(textProcessed, "<wbr>");
  }

  if(ct->onClick.size())
    onclick << " onclick='" << parseTemplates(ct->onClick, row) << "'"; // parseTemplatesAll

  if(ct->id.size())
    id << " ID='" << parseTemplates(ct->id, row) << "'"; // parseTemplatesAll


  // edit box
  if(ct->isInput) {
    ss << "<input type='text' style='text-transform: uppercase; width:100%'" << id.str() << onclick.str() << " />";
  }

  // button
  else if(ct->isButton) {
    ss << "<input type='button' value='" << textProcessed <<"'" << onclick.str() << id.str() << "' />";

  // regular text
  } else {
    ss << "<p" << id.str() << onclick.str() << ">" << textProcessed << "</p>";
  }

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
// Method to process 'multiple' table cells, to have several lines of data in the same table cell
int ReportRendererHtml::buildCellStringMultiple(ReportColumnTypeStringMultiple *ct, vector<string> *row, stringstream &ss)
{
  // auxiliary variables
  string text, link;
  vector<string> links, texts;

  // gather needed attributes
  if(ct->link.size())
    link = parseTemplates(ct->link, row); // parseTemplatesAll

  if(ct->text.size())
    text = parseTemplates(ct->text, row); // parseTemplatesAll

  StringExplode(link, links, true, "|");
  StringExplode(text, texts, true, "|");

  for(int i = 0 ; i < texts.size() ; i++) {

    // the ith link
    if(links.size() > i)
      ss << "<a href='" << ct->linkPrefix << links[i] << ct->linkSuffix << "'>";

    // regular text
    ss << "<p>" << texts[i] << "</p>";

    // the ith link
    if(links.size() > i)
      ss << "</a>";

    // line break
    if(i < texts.size()-1)
      ss << "<br />";
  }

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
int ReportRendererHtml::buildCellBox(ReportColumnTypeBox *ct, vector<string> *row, stringstream &ss)
{
  // sequences box begin sequence
  ss << "<table align='center' border='0'><tr align='center'><td>";

  // sequences box cycle thru all cells
  vector<ReportColumnTypeBase *>::iterator it = ct->sequences.begin();
  for(bool begin = true; it != ct->sequences.end() ; it++) {

    // skip dynamic
    if( (*it)->dynamic )
      continue;

    // validate cell
    if( (*it)->validator.length()) {
      string aux = parseTemplates((*it)->validator, row); // parseTemplatesAll
      if(!aux.length()) {
        continue;
      }
    }

    // if there is a column label
    if((*it)->columnLabel.size()) {
      //new line, if it is not the first one
      if(!begin)
        ss << "<br /><br />";
      // the column label, in bold
      ss << "<b>" << (*it)->columnLabel << "</b>";
      //new line
      ss << "<br />";
      // subsequent new lines between items will be inserted
      begin = false;
    }

    stringstream link;

    // gather needed attributes
    if((*it)->link.size())
      link << parseTemplates((*it)->link, row); // parseTemplatesAll

    // Link section
    if((*it)->link.size())
      ss << "<a href='" << link.str() << "'>";

    // call base class to render cell
    ReportRendererBase::buildTableCell(*it, row, ss);

    // Link section
    if((*it)->link.size())
      ss << "</a>";
  }

  // sequences box end sequence
  ss << "</td></tr></table>";

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
void ReportRendererHtml::genRandomString(string &s, const int len)
{
  static const char alphanum[] =
      "0123456789"
      "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
      "abcdefghijklmnopqrstuvwxyz";

  for (int i = 0; i < len; ++i) {
      s += alphanum[rand() % (sizeof(alphanum) - 1)];
  }
}
///////////////////////////////////////////////////////////////////////////////
void ReportRendererHtml::parseFilesVector(vector<string> *row, vector<ReportParamsFiles> &files)
{
  // check if files were loaded
  if(!spsFiles)
    return;

  // Set the files
  for(int i = 0 ; i < files.size() ; i++) {

    if( files[i].validator.length()) {
      string aux = parseTemplates(files[i].validator, row);
      if(!aux.length())
        continue;
    }

    string param = parseTemplates(files[i].param, row);

    files[i].data = spsFiles->getData(files[i].file, getInt(param.c_str()) );
    //cout << "FILE is " << files[i].file << endl;
    if(files[i].data)
      files[i].valid = true;
  }
}
///////////////////////////////////////////////////////////////////////////////
int ReportRendererHtml::renderImage(const string &object, vector<string> &params, vector<ReportParamsFiles> &files)
{
  //cout << object << " : " << paramList << endl;

  // clear image container
  m_image.clear();

  // Get the exemplar from the factory
  ReportModuleBase * module = ReportModuleFactory::getModule(object);
  if (module == 0)
    return ERROR;

  // get a clean module
  ReportModuleBase *moduleReal = module->clone();
  if(!moduleReal)
    return ERROR;

  // Set the files
  for(int i = 0 ; i < files.size() ; i++)
    if(files[i].valid)
      moduleReal->setData(files[i].type, files[i].data);


  // get params as a vector
  //vector<string> params;
  //StringExplode(paramList, params, true);
  //stringSplit(paramList, params);

  // build parameter list for module
  int count = params.size();
  char *aux[count+4];

  // Fill parameters structure
  for(int i = 0 ; i < params.size() ; i++)
    aux[i+1] = (char *)params[i].c_str();

  // add executable location/name
  string exeAux = m_exeDir; exeAux += '/' ; exeAux += object;
  aux[0] = (char *)exeAux.c_str();

  // set the output file specification command
  aux[count + 1] = OUTFILE_COMMAND;
  // set the output file name
  string fileName = "file";
  //genRandomString(fileName, 20);
  fileName += parseInt(m_fileNumber++);
  fileName += ".png";
  aux[count + 2] = (char *)fileName.c_str();

  //cout << fileName << endl;

  // set terminator
  aux[count+3] = 0;

  //string a_aa;
  //for(int i = 0 ; i < count + 3 ; i++) {
  //  a_aa += aux[i];
  //  a_aa += ' ';
  //}
  //DEBUG_MSG("Image Command : " << a_aa);

  // invoke
  try {
    moduleReal->invoke(count+3, aux);
  }
  catch(...) {
    return ERROR;
  }

  // get the generated image
  moduleReal->getData(m_image);

  // delete the object
  delete moduleReal;

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
}; //namespace
////////////////////////////////////////////////////////////////////////////////
