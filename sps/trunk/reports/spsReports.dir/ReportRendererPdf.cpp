////////////////////////////////////////////////////////////////////////////////
#include "ReportRendererPdf.h"
#include "Tokenizer.h"
#include "ReportBase64.h"
#include "copyright.h"

////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
// Page component generation methods
////////////////////////////////////////////////////////////////////////////////
PdfPoint PdfPoint::operator+(const PdfPoint &o)
{
  PdfPoint ret;
  ret.m_x = m_x + o.m_x;
  ret.m_y = m_y + o.m_y;
  return ret;
}
////////////////////////////////////////////////////////////////////////////////
//PdfPoint &PdfPoint::operator=(const PdfPoint &o)
//{
//  m_x = o.m_x;
//  m_y = o.m_y;
//  return *this;
//}
////////////////////////////////////////////////////////////////////////////////
void PdfText::calcSize(HPDF_Doc &pdf)
{
  // get font
  HPDF_Font pdfFont = HPDF_GetFont(pdf, m_font.c_str(), 0);
  // get font metrics
  m_font_descent = HPDF_Font_GetDescent(pdfFont);
  m_font_ascent = HPDF_Font_GetAscent(pdfFont);
  // get page
  HPDF_Page page = HPDF_GetCurrentPage(pdf);
  //get font size
  HPDF_Page_SetFontAndSize(page, pdfFont, m_fontSize);
  // get text width
  HPDF_REAL real_width = HPDF_Page_TextWidth (page, m_text.c_str());
  // set text height
  m_size.m_y = (float)m_fontSize;
  // set text width
  m_size.m_x = real_width;
}
////////////////////////////////////////////////////////////////////////////////
void PdfImage::calcSize(HPDF_Doc &pdf)
{
  m_size.m_x = m_size.m_y = 0.0;

  if(m_image) {
    m_size.m_x = HPDF_Image_GetWidth(m_image) * m_factor;
    m_size.m_y = HPDF_Image_GetHeight(m_image) * m_factor;
  }
}
////////////////////////////////////////////////////////////////////////////////
// Calculates the size of a table cell
void PdfTableCell::calcSize(HPDF_Doc &pdf, float cellPadding)
{
  // table cell size must include cell padding
  m_size.m_x = m_size.m_y = 0.0;

  for(int i = 0 ; i < m_text.size() ; i++)
    m_text[i].calcSize(pdf);

  for(int i = 0 ; i < m_text.size() ; i++) {
    PdfPoint aux = m_text[i].getSize();
    if(m_size.m_x < aux.m_x)
      m_size.m_x = aux.m_x;
    if(m_size.m_y < aux.m_y)
      m_size.m_y = aux.m_y;
  }

  for(int i = 0 ; i < m_image.size() ; i++)
    m_image[i].calcSize(pdf);

  for(int i = 0 ; i < m_image.size() ; i++) {
    PdfPoint aux = m_image[i].getSize();
    if(m_size.m_x < aux.m_x)
      m_size.m_x = aux.m_x;
    if(m_size.m_y < aux.m_y)
      m_size.m_y = aux.m_y;
  }

  // table cell size must include cell padding
  m_size.m_x += 2 * cellPadding;
  m_size.m_y += 2 * cellPadding;
}
////////////////////////////////////////////////////////////////////////////////
// Calculates the size of a row
void PdfTableRow::calcSize(HPDF_Doc &pdf, float cellSpacing, float cellPadding)
{
  // auxiliary object for dimensions keeping
  PdfPoint sz;
  // cycle thru all cells
  for(int i = 0 ; i < m_cell.size() ; i++) {
    // calc cell size
    m_cell[i].calcSize(pdf, cellPadding);
    // get cell size
    PdfPoint aux = m_cell[i].getSize();
    // add it to row size
    sz.m_x += aux.m_x;
    // expand row height, if cell bigger then current row size
    if(sz.m_y < aux.m_y)
      sz.m_y = aux.m_y;
  }
  // add cell spacing between cells in the same row
  sz.m_x += (m_cell.size() - 1) * cellSpacing;
  // store calculated size
  m_size = sz;
}
////////////////////////////////////////////////////////////////////////////////
// calculates the size of the table
void PdfTable::calcSizes(HPDF_Doc &pdf)
{
  HPDF_Page page = HPDF_GetCurrentPage(pdf);
  float w = HPDF_Page_GetWidth(page);
  float h = HPDF_Page_GetHeight(page);

  for(int i = 0 ; i < m_row.size() ; i++) {
    // calc row sizes
    m_row[i].calcSize(pdf, m_cellSpacing, m_cellPadding);

    // check all row cells sizes - find the larger cells to estabelish row width
    for(int j = 0 ; j < m_row[i].length() ; j++) {
      // get cell size
      PdfPoint aux = m_row[i].getSizeOf(j);

      if(j >= m_Cols.size())
        m_Cols.push_back(aux.m_x);
      else
        if(m_Cols[j] < aux.m_x)
          m_Cols[j] = aux.m_x;
    }

    // get row size
    PdfPoint aux = m_row[i].getSize();

    if(i >= m_Rows.size())
      m_Rows.push_back(aux.m_y);
    else
      if(m_Rows[i] < aux.m_y)
        m_Rows[i] = aux.m_y;
  }

  // determine table width
  m_tableWidth = m_cellSpacing * (m_Cols.size() - 1);
  for(int i = 0 ; i < m_Cols.size() ; i++)
    m_tableWidth += m_Cols[i];

  // determine table height
  m_tableHeight = m_cellSpacing * (m_Rows.size() - 1);
  for(int i = 0 ; i < m_Rows.size() ; i++)
    m_tableHeight += m_Rows[i];
}


////////////////////////////////////////////////////////////////////////////////
// Page component generation methods
////////////////////////////////////////////////////////////////////////////////
int ReportRendererPdf::paginate(ReportGeneratorData &data, ReportRendererBase &rr, ReportBase &rep)
{
/*  Timer_c timer;

  // get the total number of elements, plus the vector of IDs
  vector<string> IDs, IDsEnd, fNames;
  // navigation bar name and suffix
  string fnameSuffix = ".html";
  int nPages = rep.getPageSplitStats(DEFAULT_ROWS_PER_TABLE, IDs, IDsEnd, fNames, fnamePrefix, fnameSuffix);
  // Declare filename holder and navigation bar holder variables
  string navBar;
  // build the navigation bar
  rr.buildNavigationBar(navBar, IDs, IDsEnd, fNames);
  // set the navigation bar
  rr.setNavigationBar(navBar);
  // cycle thru all pages
  for(int j = 0 ; j < nPages ; j++) {
    // verbose output
    stringstream aux; aux << "Generating " << pageName << " page: ";
    verboseOutput(cout, aux.str().c_str(), fNames[j].c_str(), "...", false);

    // Create page
    HPDF_Page page_1;
    page_1 = HPDF_AddPage (pdf);
    HPDF_Page_SetSize (page_1, HPDF_PAGE_SIZE_LETTER, HPDF_PAGE_PORTRAIT);

    // build the report page
    rr.renderReport(&rep, page_1, j * DEFAULT_ROWS_PER_TABLE, DEFAULT_ROWS_PER_TABLE);
    // close output file
    of.close();
    // verbose output (time)
    verboseOutput(cout, timer.restart());
  } */
}
///////////////////////////////////////////////////////////////////////////////
int ReportRendererPdf::generateReport(ReportGeneratorData &data)
{
/*
  // performance timers
  Timer_c timer, timer2;


  //////////////////////////////////////////////////////////////////////////////
  // data needed to generate tables

  //data.projectDir     = m_projectDir;

  vector<string>      aux;

  // Create PDF Object
  HPDF_Doc pdf;
  pdf = HPDF_New (error_handler, NULL);
  // check for errors
  if (!pdf) {
    printf ("ERROR: cannot create pdf object.\n");
    return 1;
  }

  // static report renderer for contigs page
  ReportRendererPdf rrPDF;
  // init report
  rrPDF.initReport(pdf);
  // Specify executables directory (needed for specplot module)
  rrPDF.setExeDir(data.exeDir);
  // Specify project directory (needed for input data and report output)
  rrPDF.setProjectDir(data.projectDir);

  // specify fasta filename
  rrPDF.setDisplayLevel(data.displayLevel);


  //////////////////////////////////////////////////////////////////////////////
  // header report generation

  // define report cluster object, and load tables
  ReportHeader rh(data.tablesDir, m_tableNameHeader);
  // load tables
  rh.load();
  // build the report page
  rrPDF.renderReport(&rh, pdf);

  //////////////////////////////////////////////////////////////////////////////
  // proteins report generation

  // Proteins main page

  // verbose output
  verboseOutput(cout, "---- Generating protien pages ----");
  verboseOutput(cout, "Generating protiens page...", false);
  // define report protein list object, and load tables
  ReportProtein pl(data.tablesDir, m_tableNameProtein);
  // load table
  pl.load();

  // build the report page
  rrPDF.renderReport(&pl, pdf);


  // verbose output (time)
  verboseOutput(cout, timer.restart());


  // Child (specific) protein report pages

  // define specific protein report object
  ReportProtein p2(data.tablesDir, m_tableNameProtein, m_tableNameContig);
  // load table
  p2.load();
  // get the proteins column from the proteins table
  aux = p2.getTableColumn(0, TABLE_PROTEIN_FILTER_COL_PROTEIN);
  // cycle tru all proteins
  for(int i = 0 ; i < aux.size() ; i++) {
    // set the filter for a specific contig in cluster table
    p2.applyFilter(aux[i]);
    // build the report page
    rrPDF.renderReport(&p2, pdf);
    // generate pages with pagination
    //paginatePdf(data, rrPDF, p2);
  }


  //////////////////////////////////////////////////////////////////////////////
  // contigs report generation

  // Contigs main page

  // verbose output
  verboseOutput(cout, "---- Generating contig pages ----");
  // define report object for contig list
  ReportContig cl(data.tablesDir, m_tableNameContig);
  // load table
  cl.load();
  // generate pages with pagination
  rrPDF.renderReport(&cl, pdf);
  //paginatePdf(data, rrPDF, cl);


  // Child (specific) contig report pages

  // define specific contig report object
  ReportContig rc(data.tablesDir, m_tableNameContig, m_tableNameCluster);
  // load table
  rc.load();
  // get the contig id column from the contigs table
  aux = rc.getTableColumn(0, TABLE_CONTIG_FILTER_COL_CONTIG);
  // cycle tru all contigs
  for(int i = 0 ; i < aux.size() ; i++) {
    // set the filter for a specific contig in cluster table
    rc.applyFilter(aux[i]);
    // generate pages with pagination
    rrPDF.renderReport(&rc, pdf);
    //paginatePdf(data, rrPDF, rc);
  }
*/

/*
  //////////////////////////////////////////////////////////////////////////////
  // cluster report generation

  // Consensus cluster clusters page

  verboseOutput(cout, "---- Generating cluster pages ----");

  // define report cluster object, and load tables
  ReportCluster clust(data.tablesDir, m_tableNameCluster, m_tableNameSpectra);
  // load tables
  clust.load();
  // get the cluster id column from the cluster table
  aux = clust.getTableColumn(0, TABLE_CLUSTER_FILTER_COL_CLUSTER);
  // cycle tru all clusters
  for(int i = 0 ; i < aux.size() ; i++) {
    // set the filter for a specific contig in cluster table
    clust.applyFilter(aux[i]);
    // generate pages with pagination
    paginatePdf(data, rrPDF, clust);
  }


  //////////////////////////////////////////////////////////////////////////////
  // spectra report generation

  // verbose output
  verboseOutput(cout, "---- Generating input spectra pages ----");
  // define report cluster object, and load tables
  ReportInputSpectra is(data.tablesDir, m_tableNameSpectra);
  // load tables
  is.load();
  // get the file index column from the input spectra table
  aux = is.getTableColumn(0, TABLE_SPECTRA_FILTER_COL_FILE);
  // sort elements for duplicate removal
  sort (aux.begin(), aux.end(), stringSortCompare);
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
    // generate pages with pagination
    paginatePdf(data, rrPDF, is);
  }
*/

/*
  //////////////////////////////////////////////////////////////////////////////

  verboseOutput(cout, "---- Done! ----");
  // verbose output (time)
  verboseOutput(cout, timer2.stop());


  //////////////////////////////////////////////////////////////////////////////
  // End

  // save PDF to file
  string fn = composeFileName(data.outDir, "report.pdf");
  verboseOutput(cout, "Saving PDF report to " , fn.c_str());
  HPDF_SaveToFile (pdf, fn.c_str());
  // Clean up
  HPDF_Free (pdf);
*/
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Page component generation methods
////////////////////////////////////////////////////////////////////////////////
int ReportRendererPdf::initReport(HPDF_Doc &pdf)
{
  m_pdf = &pdf;

  // report initializer
  renderReportProlog(pdf);

  // set the page margin
  m_pageMargin = PDF_PAGE_MARGIN;

  // create first page
  createPdfPage(pdf);
}
////////////////////////////////////////////////////////////////////////////////
// Page component generation methods
////////////////////////////////////////////////////////////////////////////////
int ReportRendererPdf::renderReport(ReportBase *report, HPDF_Doc &pdf, bool renderDecoration, int start, int count)
{
  // render the page prolog
  renderTableProlog(report, pdf);

  // render the header section

  // render the tables
  ReportIterator ri = report->begin();
  for(bool first = true ; ri != report->end() ; ri++ , first = false) {
    // if not the first table being rendered, render the section between tables
    if(!first)
      renderInterTable(report, pdf);

    // clear table data container
    reportTableCells.clear();

    PdfTable table;

    // build the table
    buildTable(*ri, table);

    // set column sizes
    table.calcSizes(pdf);

    // render the table
    renderTable(*ri, table, pdf, m_pageMargin);
  }

  // render the page epilog
  renderTableEpilog(report, pdf);
}
////////////////////////////////////////////////////////////////////////////////
// Page prolog and epilogue
////////////////////////////////////////////////////////////////////////////////
// prolog
int ReportRendererPdf::renderReportProlog(HPDF_Doc &pdf)
{
  /* set compression mode */
  HPDF_SetCompressionMode (pdf, HPDF_COMP_ALL);
  /* set page mode to use outlines. */
  HPDF_SetPageMode (pdf, HPDF_PAGE_MODE_USE_OUTLINE);
  // show 2 pages
  HPDF_SetPageLayout (pdf, HPDF_PAGE_LAYOUT_TWO_COLUMN_RIGHT );
  /* set password */
  //HPDF_SetPassword (pdf, "owner", "user");

  string logoRight = "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";

  string logoLeft = "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";

  string logoMiddle = "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";

  // convert UU64 encoded image to binary
  string img2 = base64_decode(logoRight);

  // Logo image
  imageHeaderRight = HPDF_LoadJpegImageFromMem(pdf, (const HPDF_BYTE *)img2.c_str(), img2.length());

  // convert UU64 encoded image to binary
  img2 = base64_decode(logoLeft);

  // Logo image
  imageHeaderLeft = HPDF_LoadJpegImageFromMem(pdf, (const HPDF_BYTE *)img2.c_str(), img2.length());

  // convert UU64 encoded image to binary
  img2 = base64_decode(logoMiddle);

  // Logo image
  imageHeaderCenter = HPDF_LoadJpegImageFromMem(pdf, (const HPDF_BYTE *)img2.c_str(), img2.length());

  // set page number
  m_pageIndex = 0;

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Create new page
////////////////////////////////////////////////////////////////////////////////
int ReportRendererPdf::createPdfPage(HPDF_Doc &pdf)
{
  // Create page
  HPDF_Page page;
  page = HPDF_AddPage(pdf);
  // use 300 dpi
  //HPDF_Page_Concat (page, 72.0f / 300.0f, 0, 0, 72.0f / 300.0f, 0, 0);
  // define page size
  HPDF_Page_SetSize(page, HPDF_PAGE_SIZE_LETTER, HPDF_PAGE_PORTRAIT);
   // set total page size on counter
  positionH =  HPDF_Page_GetHeight(page);

  // render page header
  PdfRect r;

  if(!m_pageIndex) {
    r.top = 100;
    r.left = 10;
    r.right = 10;
    renderPageProlog(pdf, page, imageHeaderCenter, r);
  } else if(m_pageIndex % 2) {
    r.top = 0;
    r.left = 20;
    r.right = 0;
    renderPageProlog(pdf, page, imageHeaderRight, r);
  } else {
    r.top = 0;
    r.left = 0;
    r.right = 20;
    renderPageProlog(pdf, page, imageHeaderLeft, r);
  }

  // render page footer
  renderPageEpilog(pdf, page);

  // update page counter
  m_pageIndex++;

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Page prolog and epilogue
////////////////////////////////////////////////////////////////////////////////
// prolog
int ReportRendererPdf::renderPageProlog(HPDF_Doc &pdf, HPDF_Page &page, HPDF_Image image, PdfRect &r)
{
  // image dimensions
  double iw = HPDF_Image_GetWidth (image);
  double ih = HPDF_Image_GetHeight (image);

  HPDF_Page_SetLineWidth (page, 0.5);
  // page dimensions
  double pw = HPDF_Page_GetWidth (page);
  double ph = HPDF_Page_GetHeight (page);
  // multiply factor to adjust image to end of page
  double factor = (pw - r.left - r.right) / iw;
  // y start position, based on factor
  double y = ph - ih * factor;

  /* Draw image to the canvas. (normal-mode with actual size.)*/
  HPDF_Page_DrawImage(page, image, r.left, y - r.top, iw * factor, ih * factor);

  // starting position
  positionH -= ih * factor + r.top + 20;

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// epilog
int ReportRendererPdf::renderPageEpilog(HPDF_Doc &pdf, HPDF_Page &page)
{
  // page dimensions
  double pw = HPDF_Page_GetWidth (page);
  double ph = HPDF_Page_GetHeight (page);

  // line structure
  PdfLine line;

  line.m_start.m_x  = 20;
  line.m_start.m_y  = 20;
  line.m_end.m_x    = pw - 20;
  line.m_end.m_y    = 20;
  line.m_color_red  = 0.0;
  line.m_color_gree = 0.0;
  line.m_color_blue = 1.0;
  line.m_width      = 1.0;

  drawLine(page, line);


  HPDF_Font font = HPDF_GetFont(pdf, "Helvetica", 0);
  HPDF_Page_SetFontAndSize(page, font, 8);

  HPDF_Page_BeginText(page);
  HPDF_Page_MoveTextPos(page, line.m_start.m_x, line.m_start.m_y - 8);
  HPDF_Page_ShowText(page, REPORT_FOOTER_TEXT_PDF);
  HPDF_Page_EndText(page);

  // bottom limit
  PdfLimitBottom = 40;

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Table prolog and epilogue
////////////////////////////////////////////////////////////////////////////////
// prolog
int ReportRendererPdf::renderTableProlog(ReportBase *table, HPDF_Doc &pdf)
{
}
////////////////////////////////////////////////////////////////////////////////
// epilog
int ReportRendererPdf::renderTableEpilog(ReportBase *table, HPDF_Doc &pdf)
{
}
////////////////////////////////////////////////////////////////////////////////
// between tables
int ReportRendererPdf::renderInterTable(ReportBase *table, HPDF_Doc &pdf)
{

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
// PDF rendering
///////////////////////////////////////////////////////////////////////////////
int ReportRendererPdf::drawLine(HPDF_Page &page, PdfLine &line)
{
  HPDF_Page_SetLineWidth(page, line.m_width);
  HPDF_Page_SetRGBStroke(page, line.m_color_red, line.m_color_gree, line.m_color_blue);

  HPDF_Page_MoveTo(page, line.m_start.m_x, line.m_start.m_y);
  HPDF_Page_LineTo(page, line.m_end.m_x, line.m_end.m_y);
  HPDF_Page_Stroke(page);

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
int ReportRendererPdf::drawText(HPDF_Doc &pdf, HPDF_Page &page, PdfTextData &data)
{
  HPDF_Font font = HPDF_GetFont(pdf, data.font.c_str(), 0);
  HPDF_Page_SetFontAndSize(page, font, data.fontSize);
  //HPDF_Page_SetRGBStroke(page, data.color.m_red, data.color.m_green, data.color.m_blue);
  HPDF_Page_SetRGBFill (page, data.color.m_red, data.color.m_green, data.color.m_blue);

  HPDF_Page_BeginText(page);
  // show all text items
  HPDF_Page_MoveTextPos(page, data.rect.left, data.rect.top );
  HPDF_Page_ShowText(page, data.text.c_str() );
  // close text output
  HPDF_Page_EndText(page);  // render the header
}
///////////////////////////////////////////////////////////////////////////////
int ReportRendererPdf::drawTextRect(HPDF_Doc &pdf, HPDF_Page &page, PdfTextData &data)
{
  //HPDF_Font font = HPDF_GetFont(pdf, data.font.c_str(), 0);
  //HPDF_Page_SetFontAndSize(page, font, data.fontSize);
  HPDF_Page_SetRGBFill (page, 0.75, 0.0, 0.0);

  HPDF_Page_Rectangle(page, data.area.left, data.area.bottom, data.area.right - data.area.left, data.area.top - data.area.bottom);
  //HPDF_Page_Stroke (page);
  HPDF_Page_Fill (page);

  //HPDF_Page_SetRGBFill(page, data.color.m_red, data.color.m_green, data.color.m_blue);

  // show all text items
  //HPDF_Page_BeginText(page);

  //HPDF_Page_MoveTextPos(page, data.rect.left, data.rect.bottom );
  //HPDF_Page_TextRect (page, data.rect.left, data.rect.top, data.rect.right, data.rect.bottom, data.text.c_str(), (HPDF_TextAlignment)data.align, NULL);
  // close text output
  //HPDF_Page_EndText(page);  // render the header
}
///////////////////////////////////////////////////////////////////////////////
int ReportRendererPdf::drawImage(HPDF_Doc &pdf, HPDF_Page &page, PdfImageData &data)
{
  /* Draw image to the canvas. (normal-mode with actual size.)*/
  HPDF_Page_DrawImage(page, data.image, data.rect.left, data.rect.bottom, data.rect.right - data.rect.left, data.rect.top - data.rect.bottom);
}
///////////////////////////////////////////////////////////////////////////////
// Table rendering exceptions
///////////////////////////////////////////////////////////////////////////////
void ReportRendererPdf::breakProteinIntoChunks(vector<string> &in, int &count)
{
  count = 0;
  for(int i = 0 ; i < in.size() ; i++) {
    // string container
    string aux;
    for(int j = 0 ; j < in[i].length() ; j++) {
      // line break;
      if(count % 50 == 0) {
        aux += "<br>";
      // spacer
      } else if(count % 10 == 0) {
        aux += " ";
      }
      // copy element
      aux += in[i][j];
      count++;
    }
    in[i] = aux;
  }
}
///////////////////////////////////////////////////////////////////////////////
void ReportRendererPdf::colorProteinString(vector<string> &in, string &out)
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
int ReportRendererPdf::renderTableExceptionMainPage(ReportTableBase *table, HPDF_Doc &pdf)
{
/*  // declare and initialize table row interator, to allow direct access to the table
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
  outstream << "    <td style='background-color:lightgreen;'>" << status << " %</td>";
  outstream << "  </tr>";
  outstream << "  <tr>";
  outstream << "    <th width='25%' bgcolor='#003399'><span style='color:white'>Elapsed</span></th>";
  outstream << "    <td>" << elapsed << "</td>";
  outstream << "  </tr>";
  outstream << "  <tr>";
  outstream << "    <th width='25%' bgcolor='#003399'><span style='color:white'>Log</span></th>";
  outstream << "    <td><a href='spsplot.log'>" << log << "</a></td>";
  outstream << "  </tr>";
  outstream << "  <tr>";
  outstream << "    <th width='25%' bgcolor='#003399' rowspan='2'><span style='color:white'>Data</span></th>";
  outstream << "    <td><a href='" << contig << ".0.html'>Group by Contig</a></td>";
  outstream << "  </tr>";
  outstream << "  <tr>";
  outstream << "    <td><a href='" << protein << ".html'>Group by Protein</a></td>";
  outstream << "  </tr>";
  outstream << "  <tr>";
  outstream << "    <th width='25%' bgcolor='#003399'><span style='color:white'>Cluster Data</span></th>";
  outstream << "    <td>";
  outstream << "      <a href='" << cluster << ".txt'>All Clusters (txt)</a><br>";

  for(int i = 0 ; i < files.size() && i < indices.size() ; i++)
    outstream << "      <a href='spectra." << indices[i] << ".0.html'>Group by <i>" << files[i] << "</i></a><br>";

  outstream << "    </td>";
  outstream << "  </tr>";
  outstream << "  </td></tr></table></td></tr>";
  outstream << "  <tr>";
  outstream << "    <td colspan='0' class='bottomline'>&nbsp;</td>";
  outstream << "  </tr>";
  outstream << "</table>";
  outstream << "</div></div>";
*/
}
///////////////////////////////////////////////////////////////////////////////
int ReportRendererPdf::renderTableExceptionProteinHeader(ReportTableBase *table, HPDF_Doc &pdf)
{

  // declare and initialize table row interator, to allow direct access to the table
  TableIterator ti = table->begin();

  // get the protein name, at indice 1
  string protienId    = (*(*ti))[0];
  string proteinName  = (*(*ti))[1];
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
    legend += "<br>";
    legend += parseInt(current);
    current += 50;
  }

  // protein header HTML code
  //outstream << "<table class='result' width='100%' style='border-spacing: 0px;'>";
  //outstream << "<tr>";
  //outstream << "<td colspan='0'><h2><i>" << proteinName << "</i></h2>";
  //outstream << "<hr><b>" << contigs << " contigs, " << spectra << " spectra, " << aas << " amino acids, " << coverage << " //coverage" << "</b></td>";
  //outstream << "<td></td>";
  //outstream << "</tr>";
  //outstream << "<tr> ";
  //outstream << "<td align='right'><tt>" << legend << "</tt></td>";
  //outstream << "<td><tt>" << coloredProtein << "</tt></td>";
  //outstream << "</tr>";
  //outstream << "</table>";

  // link to protein detains page
  //outstream << "<div>";
  //outstream << "<a href='protein_details." << protienId << ".html'>Protein coverage</a>";
  //outstream << "</div>";
  //outstream << "<br>";

  // return status OK
  return OK;
}
///////////////////////////////////////////////////////////////////////////////
// Table rendering methods
///////////////////////////////////////////////////////////////////////////////
  // render table comon method. Cycles thru all the table rows and invokes buildTableRow() method to render a row
int ReportRendererPdf::renderTable(ReportTableBase *table, PdfTable &pdfTable, HPDF_Doc &pdf, float pageMargin)
{

  switch(table->getRenderingException()) {

  case RENDERING_EXCEPTION_NONE:
    break;

  case RENDERING_EXCEPTION_PROTEIN_HEADER:
    return renderTableExceptionProteinHeader(table, pdf);
    break;

  case RENDERING_EXCEPTION_MAIN_PAGE:
    return renderTableExceptionMainPage(table, pdf);

  default:
    break;
  }

  // test for empty tables
  if(!pdfTable.length())
    return OK;

  if(pdfTable[0].length()) {
    renderTableHeaderRow(pdfTable[0], pdf, pdfTable.getRowAt(0), pdfTable.getCols(), pageMargin, pdfTable.getTableWidth(), pdfTable.getCellSpacing(), pdfTable.getCellPadding());
    positionH -= pdfTable.getRowAt(0) + pdfTable.getCellSpacing();
  }

  // render content rows
  for(int i = 1 ; i < pdfTable.length() ; i++)
    if(pdfTable[i].length()) {

      if(pdfTable.getRowAt(i) > positionH - PdfLimitBottom) {
        createPdfPage(pdf);
      }

      renderTableRow(pdfTable[i], pdf, pdfTable.getRowAt(i), pdfTable.getCols(), pageMargin, pdfTable.getTableWidth(), pdfTable.getCellSpacing(), pdfTable.getCellPadding());
      positionH -= pdfTable.getRowAt(i) + pdfTable.getCellSpacing();
    }

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// render a header row method. Renders a table header row.
int ReportRendererPdf::renderTableHeaderRow(PdfTableRow &row, HPDF_Doc &pdf, float rowHeight, vector<float> &cols, float pageMargin, float tableWidth, float cellSpacing, float cellPadding)
{
  // get current page
  HPDF_Page page = HPDF_GetCurrentPage(pdf);

  // render all table cells
  float aux = pageMargin;
  for(int i = 0 ; i < row.length() ; i++) {
    renderTableHeaderCell(row[i], pdf, page, aux, cellPadding, rowHeight, cols[i]);
    aux += (i < cols.size() ? cols[i] : 40) + cellSpacing;
  }

  // line structure
  PdfLine line;

  line.m_start.m_x  = pageMargin;
  line.m_start.m_y  = positionH;
  line.m_end.m_x    = tableWidth + pageMargin;
  line.m_end.m_y    = positionH;
  line.m_color_red  = 0.0;
  line.m_color_gree = 0.0;
  line.m_color_blue = 0.0;
  line.m_width      = 1.0;

  drawLine(page, line);

  line.m_start.m_y  = positionH - rowHeight;
  line.m_end.m_y    = positionH - rowHeight;

  drawLine(page, line);

  aux = pageMargin;
  for(int i = 0 ; i <= row.length() ; i++) {

    line.m_start.m_x  = aux;
    line.m_start.m_y  = positionH;
    line.m_end.m_x    = aux;
    line.m_end.m_y    = positionH - rowHeight;

    drawLine(page, line);

    aux += (i < cols.size() ? cols[i] : 40) + cellSpacing;
  }

  // return status OK
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// render a header cell method. Renders a table header cell.
int ReportRendererPdf::renderTableHeaderCell(PdfTableCell &cell, HPDF_Doc &pdf, HPDF_Page &page, float posX, float cellPadding, float height, float width)
{
  // show all text items
  for(int i = 0 ; i < cell.getTextLength() ; i++) {

    PdfTextData data;
    data.rect.left    = posX + cellPadding;
    data.rect.top     = positionH - cellPadding - cell.getTextAt(i).getFontSize() - cell.getTextAt(i).getFontSize() * (float)(cell.getTextAt(i).getDescent() / (float)cell.getTextAt(i).getAscent());
    data.rect.right   = posX + width - cellPadding;
    data.rect.bottom  = positionH - height + cellPadding;
    data.area.left    = posX;
    data.area.top     = positionH;
    data.area.right   = posX + width;
    data.area.bottom  = positionH - height;
    data.font         = cell.getTextAt(i).getFont();
    data.fontSize     = cell.getTextAt(i).getFontSize();
    data.text         = cell.getTextAt(i).getText();

    drawTextRect(pdf, page, data);

    drawText(pdf, page, data);
  }

  // return status OK
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// build a row comon method. Cycles thu all columnType cells and buildCell for each ColumnType item
int ReportRendererPdf::renderTableRow(PdfTableRow &row, HPDF_Doc &pdf, float rowHeight, vector<float> &cols, float pageMargin, float tableWidth, float cellSpacing, float cellPadding)
{
  // get current page
  HPDF_Page page = HPDF_GetCurrentPage(pdf);

  // render all table rows
  float aux = pageMargin;
  for(int i = 0 ; i < row.length() ; i++) {
    renderTableCell(row[i], pdf, page, aux, cellPadding);
    aux += (i < cols.size() ? cols[i] : 40) + cellSpacing;
  }

  // line structure
  PdfLine line;

  line.m_start.m_x  = pageMargin;
  line.m_start.m_y  = positionH - rowHeight;
  line.m_end.m_x    = tableWidth + pageMargin;
  line.m_end.m_y    = positionH - rowHeight;
  line.m_color_red  = 0.0;
  line.m_color_gree = 0.0;
  line.m_color_blue = 0.0;
  line.m_width      = 1.0;

  drawLine(page, line);

  aux = pageMargin;
  for(int i = 0 ; i <= row.length() ; i++) {

    line.m_start.m_x  = aux;
    line.m_start.m_y  = positionH;
    line.m_end.m_x    = aux;
    line.m_end.m_y    = positionH - rowHeight;

    drawLine(page, line);

    aux += (i < cols.size() ? cols[i] : 40) + cellSpacing;
  }

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// render a cell method. Renders a table cell.
int ReportRendererPdf::renderTableCell(PdfTableCell &cell, HPDF_Doc &pdf, HPDF_Page &page, float posX, float cellPadding)
{
  for(int i = 0 ; i < cell.getTextLength() ; i++) {

    PdfTextData data;
    data.rect.left    = posX + cellPadding;
    data.rect.top     = positionH - cellPadding - cell.getTextAt(i).getFontSize() - cell.getTextAt(i).getFontSize() * (float)(cell.getTextAt(i).getDescent() / (float)cell.getTextAt(i).getAscent());
    data.font         = cell.getTextAt(i).getFont();
    data.fontSize     = cell.getTextAt(i).getFontSize();
    data.text         = cell.getTextAt(i).getText();

    drawText(pdf, page, data);
  }

  for(int i = 0 ; i < cell.getImageLength() ; i++) {

    PdfImageData data;
    data.rect.left    = posX + cellPadding;
    data.rect.top     = positionH - cellPadding;
    data.rect.right   = posX + cellPadding + cell.getImageAt(i).getSize().m_x;
    data.rect.bottom  = positionH - cellPadding - cell.getImageAt(i).getSize().m_y;
    data.image        = cell.getImageAt(i).getImage();
    drawImage(pdf, page, data);
  }

  // return status OK
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Table Building methods
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// build table comon method. Renders all ColumnTypes
int ReportRendererPdf::buildTable(ReportTableBase *table, PdfTable &pdfTable)
{
  // set current row
  currentRow = -1;

  // declare target row
  PdfTableRow pdfTableRow;
  // Render the table headers
  if(table->doHeaders()) {
    // render it
    buildTableHeaderRow(table, pdfTableRow);
  }
  // store row in rendered table
  pdfTable.addRow(pdfTableRow);


  // set current row
  currentRow = 0;

  // declare and initialize table row interator
  TableIterator ti = table->begin();

  // cycle thru all rows
  for(; ti != table->end() ; ti++, currentRow++) {
    // get the row
    vector<string> *row = *ti;
    // declare target row
    PdfTableRow  pdfTableRow;
    // render it
    if(buildTableRow(table, row, pdfTableRow) == ERROR)
      return ERROR;
    // store row in rendered table
    pdfTable.addRow(pdfTableRow);
  }

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// build header row comon method. Renders all ColumnTypes
int ReportRendererPdf::buildTableHeaderRow(ReportTableBase *table, PdfTableRow  &pdfTableRow)
{
  // declare and initialize table column interator
  TableCtIterator tci = table->ctBegin();
  // cycle thru all column types
  for( ; tci != table->ctEnd() ; tci++) {
    // get the column type item
    ReportColumnTypeBase *base = *tci;
    // declare string to hold cell data
    PdfTableCell pdfTableCell;
    // render it
    if(buildTableHeaderCell(base, pdfTableCell) == ERROR)
      return ERROR;
    //store cell
    pdfTableRow.addCell(pdfTableCell);
  }

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// build header cell comon method. Renders all ColumnTypes
int ReportRendererPdf::buildTableHeaderCell(ReportColumnTypeBase *base, PdfTableCell &pdfTableCell)
{
  PdfText t(base->columnLabel);

  pdfTableCell.addText(t);

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// build header row comon method. Renders all ColumnTypes
int ReportRendererPdf::buildTableRow(ReportTableBase *table, vector<string> *row, PdfTableRow &pdfTableRow)
{
  // set current column
  currentCol = 0;

  // declare and initialize table column interator
  TableCtIterator tci = table->ctBegin();
  // cycle thru all column types
  for( ; tci != table->ctEnd() ; tci++, currentCol++) {
    // get the column type item
    ReportColumnTypeBase *base = *tci;
    // declare string to hold cell data
    PdfTableCell pdfTableCell;
    // render it
    if(buildTableCell(base, row, pdfTableCell) == ERROR)
      return ERROR;
    //store cell
    pdfTableRow.addCell(pdfTableCell);
  }

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
// render cell comon method. Renders a specific cell based on ColumnType specifications and a row of data
int ReportRendererPdf::buildTableCell(ReportColumnTypeBase *base, vector<string> *row, PdfTableCell &pdfTableCell)
{
 // check for valid cells
  if(base->validator.length()) {
    string aux = parseTemplates(base->validator, row);
    if(!aux.length())
      return OK;
  }

  // test for Image On Demand column type
  ReportColumnTypeImageOnDemand *auxI = dynamic_cast<ReportColumnTypeImageOnDemand *>(base);
  if(auxI != NULL)
    return buildCellImageOnDemand(auxI, row, pdfTableCell);

  // test for String column type
  ReportColumnTypeString *auxS = dynamic_cast<ReportColumnTypeString*>(base);
  if(auxS != NULL)
    return buildCellString(auxS, row, pdfTableCell);

  // test for String column type
  ReportColumnTypeStringMultiple *auxM = dynamic_cast<ReportColumnTypeStringMultiple*>(base);
  if(auxM != NULL)
    return buildCellStringMultiple(auxM, row, pdfTableCell);

  // test for Box column type
  ReportColumnTypeBox *auxB = dynamic_cast<ReportColumnTypeBox*>(base);
  if(auxB != NULL)
    return buildCellBox(auxB, row, pdfTableCell);

  // something's wrong...
  return ERROR;
}
///////////////////////////////////////////////////////////////////////////////
// Renderers for the diferent cell types
int ReportRendererPdf::buildCellImageOnDemand(ReportColumnTypeImageOnDemand *ct, vector<string> *row, PdfTableCell &pdfTableCell)
{
  // auxiliary variables
  stringstream icon, label, url;

  // Icon path/image
  if(ct->iconParams.size()) {

    // get parsed parameters / URL
    vector<string> pars; // = parseTemplatesAll(ct->iconParams, row);
    parseParamsVector(pars, ct->iconParams, row);
    // if there is a renderer, use it to render the image
    if(ct->iconRenderer.size()  && ct->iconDisplayLevel < m_displayLevel) {
      //renderImage(ct->iconRenderer, parseTemplatesAll(pars, row));
      renderImage(ct->iconRenderer, pars);
      // convert UU64 encoded image to binary
      string img2 = base64_decode(m_image);
      // addImage
      HPDF_Image image = HPDF_LoadPngImageFromMem(*m_pdf, (const HPDF_BYTE *)img2.c_str(), img2.length());

      PdfImage img(image);
      img.setFactor(0.5);
      pdfTableCell.addImage(img);
    }
  }


  // label to show for the link
  if(ct->label.size())
    label << parseTemplates(ct->label, row); // parseTemplatesAll
  // URL template to be used to get the image
  if(ct->renderer.size() && ct->linkDisplayLevel < m_displayLevel) {
    //url << ct->renderer << ' ' <<  parseTemplatesAll(ct->params, row);
    //renderImage(ct->renderer, parseTemplatesAll(ct->params, row));
    //outstream << "<img src='data:image/png;base64," << m_image << "' />";
    vector<string> pars;
    parseParamsVector(pars, ct->params, row);
    //renderImage(ct->renderer, pars);

    //outstream << "<a href='" << url.str() << "'>";
    string image64 = "data:image/png;base64,";
    image64 += m_image;
    //outstream << "<a href='#' onClick='JavaScript: var a=\"" << m_image << "\";showImage(a);'>";
    //ss << "<a href=\"" << image64 << "\" rel=\"lightbox\">";
  }

  // button
  if(ct->label.size())
    //ss << label.str();
    //outstream << "<div onClick='javascript:alert(\"#1\");var a='" << m_image << "';showImage(a);'>" << label.str() << "</div>";

  // url for IOD - terminator
  if(ct->renderer.size())
    //ss << "<a>";

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
int ReportRendererPdf::buildCellString(ReportColumnTypeString *ct, vector<string> *row, PdfTableCell &pdfTableCell)
{
  // auxiliary variables
  string text;

  // gather needed attributes
  if(ct->text.size())
    text = parseTemplates(ct->text, row); // parseTemplatesAll

  // generate text object
  PdfText t(text);
  // add to cell object
  pdfTableCell.addText(t);

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
// Method to process 'multiple' table cells, to have several lines of data in the same table cell
int ReportRendererPdf::buildCellStringMultiple(ReportColumnTypeStringMultiple *ct, vector<string> *row, PdfTableCell &pdfTableCell)
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
    //if(links.size() > i)
    //  ss << "<a href='" << ct->linkPrefix << links[i] << ct->linkSuffix << "'>";

    // regular text
    //ss << "<p>" << texts[i] << "</p>";

    // the ith link
    //if(links.size() > i)
    //  ss << "</a>";

    // line break
    //if(i < texts.size()-1)
    //  ss << "<br>";
  }

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
int ReportRendererPdf::buildCellBox(ReportColumnTypeBox *ct, vector<string> *row, PdfTableCell &pdfTableCell)
{
  // sequences box begin sequence
  //ss << "<table align='center' border='0'><tr align='center'><td>";

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
      //if(!begin)
      //  ss << "<br><br>";
      // the column label, in bold
      //ss << "<b>" << (*it)->columnLabel << "</b>";
      //new line
      //ss << "<br>";
      // subsequent new lines between items will be inserted
      begin = false;
    }

    stringstream link;

    // gather needed attributes
    if((*it)->link.size())
      link << parseTemplates((*it)->link, row); // parseTemplatesAll

    // Link section
    //if((*it)->link.size())
    //  ss << "<a href='" << link.str() << "'>";

    // call base class to render cell
    //ReportRendererBase::buildTableCell(*it, row, ss);

    // Link section
    //if((*it)->link.size())
    //  ss << "</a>";
  }

  // sequences box end sequence
  //ss << "</td></tr></table>";

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
int ReportRendererPdf::renderImage(const string &object, vector<string> & params)
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

  // get params as a vector
  //vector<string> params;
  //StringExplode(paramList, params, true);
  //stringSplit(paramList, params);

  // build parameter list for module
  int count = params.size();
  char *aux[count+2];

  // Fill parameters structure
  for(int i = 0 ; i < params.size() ; i++)
    aux[i+1] = (char *)params[i].c_str();

  // add executable location/name
  string exeAux = m_exeDir; exeAux += '/' ; exeAux += object;
  aux[0] = (char *)exeAux.c_str();
  // set terminator
  aux[count+1] = 0;

  // invoke
  try {
    moduleReal->invoke(count+1, aux);
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
