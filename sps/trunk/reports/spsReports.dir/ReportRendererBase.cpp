////////////////////////////////////////////////////////////////////////////////
#include "ReportRendererBase.h"
#include "utils.h"
#include "Tokenizer.h"
#include "ReportBase64.h"

////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
// render report. Cycles thru all the report tables and renders them, by invoing renderTable()
ReportRendererBase::ReportRendererBase() :
   m_tableNameHeader("tableHeader.txt"),
   m_tableNameProtein("tableProtein.txt"),
   m_tableNameProteinCoverage("tableProteinCoverage.txt"),
   m_tableNameContig("tableContig.txt"),
   m_tableNameCluster("tableCluster.txt"),
   m_tableNameSpectra("tableSpectra.txt"),
   m_displayLevel(100),
   m_noClusters(false),
   m_allowRealign(true),
   m_dynamic(true),
   m_ncpu(1),
   m_verbose(false),
   m_pwd(false)
{
  // Fill the factory with all known module types
  ReportModuleFactory::RegisterAllModules();
}
///////////////////////////////////////////////////////////////////////////////
void ReportRendererBase::verboseOutput(ostream &os, const char *str, bool nl)
{
  if(!m_verbose) return;
  os << str;
  if(nl)
    os << endl;
  else
    os.flush();
}
///////////////////////////////////////////////////////////////////////////////
void ReportRendererBase::verboseOutput(ostream &os, const char *str, const char *val, bool nl)
{
  if(!m_verbose) return;
  os << str << val;
  if(nl)
    os << endl;
  else
    os.flush();
}
///////////////////////////////////////////////////////////////////////////////
void ReportRendererBase::verboseOutput(ostream &os, const char *str, const char *val, const char *term, bool nl)
{
  if(!m_verbose) return;
  os << str << val << term;
  if(nl)
    os << endl;
  else
    os.flush();
}
void ReportRendererBase::verboseOutput(ostream &os, const double &v, bool nl)
{
  if(!m_verbose) return;
  os << v;
  if(nl)
    os << endl;
  else
    os.flush();
}
///////////////////////////////////////////////////////////////////////////////
string ReportRendererBase::composeFileName(const string &projectDir, const string &fileName)
{
  // Compose output path
  string aux = projectDir;
  if(aux[aux.length()-1] != '/')
    aux += '/';
  // add filename
  aux += fileName;
  // return composed filename
  return aux;
}
////////////////////////////////////////////////////////////////////////////////
int ReportRendererBase::renderReport(ReportBase *report, ostream &outstream, string &navBar, bool renderDecoration, int start, int count)
{
  // render the page prolog
  if(renderDecoration)
    renderProlog(report, outstream);

  // render the header section

  // render the inter-page navigation bar
  if(navBar.length() && renderDecoration) {
    outstream << navBar;
    outstream << "<div><table width='100%'><tr></td><td class='ln'></td></tr><tr><td class='VHSep'></td></tr></table></div>";
  }

  // render the navigation bar (for split pages)
  if(m_navBar.length() && renderDecoration)
    outstream << m_navBar;

  // render the tables
  ReportIterator ri = report->begin();
  for(bool first = true ; ri != report->end() ; ri++ , first = false) {
    // if not the first table being rendered, render the section between tables
    if(!first && renderDecoration)
      renderInterTable(report, outstream);

    // clear table data container
    reportTableCells.clear();

    // build the table
    buildTable(*ri, start, count);

    // render the table
    renderTable(*ri, outstream);
  }

  // render the navigation bar (for split pages) - bottom of page
  if(m_navBar.length() && renderDecoration)
    outstream << m_navBar;

  // render the page epilog
  if(renderDecoration)
    renderEpilog(report, outstream);
}
////////////////////////////////////////////////////////////////////////////////
ReportRendererBase::~ReportRendererBase()
{
  ReportModuleFactory::cleanup();
}
////////////////////////////////////////////////////////////////////////////////
// Header section rendering (top of page)
////////////////////////////////////////////////////////////////////////////////
int ReportRendererBase::buildElement(ReportElementsBase *elem, ostream &outstream)
{
  // test for Element Div type
  ReportElementsDiv *eleDiv = dynamic_cast<ReportElementsDiv *>(elem);
  if(eleDiv != NULL)
    return buildElementDiv(eleDiv, outstream);

  // test for Element Table type
  ReportElementsTable *eleTable = dynamic_cast<ReportElementsTable *>(elem);
  if(eleTable != NULL)
    return buildElementTable(eleTable, outstream);

  // test for Element Row type
  ReportElementsRow *eleRow = dynamic_cast<ReportElementsRow*>(elem);
  if(eleRow != NULL)
    return buildElementRow(eleRow, outstream);

  // test for Element Cell type
  ReportElementsCell *eleCell = dynamic_cast<ReportElementsCell*>(elem);
  if(eleCell != NULL)
    return buildElementCell(eleCell, outstream);

  // test for Element String type
  ReportElementsString *eleStr = dynamic_cast<ReportElementsString*>(elem);
  if(eleStr != NULL)
    return buildElementString(eleStr, outstream);

  // something's wrong...
  return ERROR;
}
///////////////////////////////////////////////////////////////////////////////
// Table rendering methods
///////////////////////////////////////////////////////////////////////////////
void ReportRendererBase::transposeRow(int col, vector<string> &row)
{
  row.clear();
  for(int i = 0 ; i < reportTableCells.size() ; i++)
    row.push_back(reportTableCells[i][col]);
}
///////////////////////////////////////////////////////////////////////////////
// render table comon method. Cycles thru all the table rows and invokes buildTableRow() method to render a row
int ReportRendererBase::renderTable(ReportTableBase *table, ostream &outstream)
{
  // test for empty tables
  if(!reportTableCells.size())
    return OK;

  if(table->getDirection() == 1)
    return renderTableVertically(outstream);

  if(reportTableCells[0].size())
    renderTableHeaderRow(reportTableCells[0], outstream);

  // render content rows
  for(int i = 1 ; i < reportTableCells.size() ; i++)
    if(reportTableCells[i].size())
      renderTableRow(reportTableCells[i], outstream);

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
// render table comon method. Cycles thru all the table rows and invokes buildTableRow() method to render a row
int ReportRendererBase::renderTableVertically(ostream &outstream)
{
  // test for empty tables
  if(!reportTableCells.size())
    return OK;

  // transposed row container
  vector<string> row;

  for(int i = 0 ; i < reportTableCells[0].size() ; i++) {
    transposeRow(i, row);
    if(row.size())
      renderTableHeaderRow(row, outstream);
  }

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// render a header row method. Renders a table header row.
int ReportRendererBase::renderTableHeaderRow(vector<string> &row, ostream &outstream)
{
  // render all table rows
  for(int i = 0 ; i < row.size() ; i++)
    renderTableHeaderCell(row[i], outstream);

  // return status OK
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// render header cell method. Renders a header cell
int ReportRendererBase::renderTableHeaderCell(string &cell, ostream &outstream)
{
  // output the string
  outstream << cell;

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// build a row comon method. Cycles thu all columnType cells and buildCell for each ColumnType item
int ReportRendererBase::renderTableRow(vector<string> &row, ostream &outstream)
{
  // render all table rows
  for(int i = 0 ; i < row.size() ; i++)
    renderTableCell(row[i], outstream);

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// build cell comon method. builds a specific cell based on ColumnType specifications and a row of data
int ReportRendererBase::renderTableCell(string &cell, ostream &outstream)
{
  // output the string
  outstream << cell;
}
////////////////////////////////////////////////////////////////////////////////
// Table Building methods
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// build table. Cycles thru all the table elements and invokes renderCell() method to render a cell
int ReportRendererBase::buildTable(ReportTableBase *table, int start, int count)
{
  // set current row
  currentRow = -1;

  // declare target row
  vector<string>  rowHeader;
  // Render the table headers
  if(table->doHeaders()) {
    // render it
     buildTableHeaderRow(table, rowHeader);
  }
  // store row in rendered table
  reportTableCells.push_back(rowHeader);



  // set current row
  currentRow = 0;

  // declare and initialize table row interator
  TableIterator ti = table->begin();
  if(start > 0)
    ti.moveTo(start);
  // cycle thru all rows
  for(; ti != table->end() && ( (currentRow < count) || (count == -1) ) ; ti++, currentRow++) {
    // get the row
    vector<string> *row = *ti;
    // declare target row
    vector<string>  tableRow;
    // render it
    if(buildTableRow(table, row, tableRow) == ERROR)
      return ERROR;
    // store row in rendered table
    reportTableCells.push_back(tableRow);
  }

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// render table. Cycles thru all the table elements and invokes renderCell() method to render a cell
int ReportRendererBase::buildTableHeaderRow(ReportTableBase *table, vector<string> &renderedRow)
{
  // declare and initialize table column interator
  TableCtIterator tci = table->ctBegin();
  // cycle thru all column types
  for( ; tci != table->ctEnd() ; tci++) {
    // get the column type item
    ReportColumnTypeBase *base = *tci;
    // declare string to hold cell data
    stringstream ss;
    // render it
    if(buildTableHeaderCell(base, ss) == ERROR)
      return ERROR;
    //store cell
    renderedRow.push_back(ss.str());
  }

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// render table. Cycles thru all the table elements and invokes renderCell() method to render a cell
int ReportRendererBase::buildTableRow(ReportTableBase *table, vector<string> *row, vector<string> &renderedRow)
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
    stringstream ss;
    // render it
    if(buildTableCell(base, row, ss) == ERROR)
      return ERROR;
    //store cell
    renderedRow.push_back(ss.str());
  }

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// render a cell. This method is renderer specific, and is the one who specifies how each cell should be rendered
int ReportRendererBase::buildTableCell(ReportColumnTypeBase *base, vector<string> *row, stringstream &str)
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
    return buildCellImageOnDemand(auxI, row, str);

  // test for String column type
  ReportColumnTypeString *auxS = dynamic_cast<ReportColumnTypeString*>(base);
  if(auxS != NULL)
    return buildCellString(auxS, row, str);

  // test for String column type
  ReportColumnTypeStringMultiple *auxM = dynamic_cast<ReportColumnTypeStringMultiple*>(base);
  if(auxM != NULL)
    return buildCellStringMultiple(auxM, row, str);

  // test for Box column type
  ReportColumnTypeBox *auxB = dynamic_cast<ReportColumnTypeBox*>(base);
  if(auxB != NULL)
    return buildCellBox(auxB, row, str);

  // something's wrong...
  return ERROR;
}
///////////////////////////////////////////////////////////////////////////////
void ReportRendererBase::parseParamsVector(vector<string> &params, vector<ReportParamsOption> &options, vector<string> *row)
{
  for(int i = 0 ; i < options.size() ; i++) {

    if(options[i].validator.length()) {
      string aux = parseTemplates(options[i].validator, row); // parseTemplatesAll
      if(!aux.length()) {
        continue;
      }
    }

    params.push_back(options[i].param);
    if(options[i].option.length())
      params.push_back(parseTemplates(options[i].option, row)); // parseTemplatesAll
  }
}
///////////////////////////////////////////////////////////////////////////////
string ReportRendererBase::parseTemplatesAll(const string &str, vector<string> *row,
   char primaryOpen, char primaryClose,
   char secundaryOpen, char secundaryClose,
   int add)
{
  // string ret = parseTemplates(str, row);
  // return parseTemplates(ret, row, '«', '»', 1);

  // string to hold open tag options
  string normal; normal = primaryOpen;
  string both(normal); both += secundaryOpen;

  // return string and auxiliary string to store tag contents
  string ret, tag;
  // store the position after each tag
  size_t lastPosition = 0;
  // get the first tag
  size_t first = str.find_first_of(both, lastPosition);
  // declare variable to hold second delimiter
  size_t second;
  // get next tag start
  bool bNormal = normal.find_first_of(str[first]) != string::npos;
  // which was it?
  if(bNormal)
    // if normal tag
    second = str.find_first_of(primaryClose, first+1);
  else
    // if add tag
    second = str.find_first_of(secundaryClose, first+1);

  // repeat until all the tags are processed
  while(first != string::npos && second != string::npos) {
    // copy in beetwen tags
    ret += str.substr(lastPosition, first - lastPosition);
    // copy the tag contents to an auxiliary string
    tag = str.substr(first+1, second - first - 1);

    // translate the tag contents
    if(bNormal)
      ret += getTag(tag, row, primaryOpen, primaryClose);
    else
      ret += getTag(tag, row, secundaryOpen, secundaryClose, add);

    // update last position to past the tag
    lastPosition = second + 1;

    // get the first tag
    first = str.find_first_of(both, lastPosition);
    // get next tag
    bNormal = normal.find_first_of(str[first]) != string::npos;
    // which was it?
    if(bNormal)
      // if normal tag
      second = str.find_first_of(primaryClose, first+1);
    else
      // if add tag
      second = str.find_first_of(secundaryClose, first+1);
  }

  // the remaining of the string
  ret += str.substr(lastPosition);
  // return the string
  return ret;
}
///////////////////////////////////////////////////////////////////////////////
string ReportRendererBase::parseTemplates(const string &str, vector<string> *row, char tagOpen, char tagClose, int add)
{
  // return string and auxiliary string to store tag contents
  string ret, tag;
  // store the position after each tag
  size_t lastPosition = 0;
  // get the first tag
  size_t first = str.find_first_of(tagOpen, lastPosition);
  size_t second = str.find_first_of(tagClose, first+1);
  // repeat until all the tags are processed
  while(first != string::npos && second != string::npos) {
    // copy in beetwen tags
    ret += str.substr(lastPosition, first - lastPosition);
    // copy the tag contents to an auxiliary string
    tag = str.substr(first+1, second - first - 1);
    // translate the tag contents
    ret += getTag(tag, row, tagOpen, tagClose, add);
    // update last position to past the tag
    lastPosition = second + 1;
    // get next tag
    first = str.find_first_of(tagOpen, lastPosition);
    second = str.find_first_of(tagClose, first+1);
  }

  // the remaining of the string
  ret += str.substr(lastPosition);
  // return the string
  return ret;
}
///////////////////////////////////////////////////////////////////////////////
string ReportRendererBase::getTag(const string &tag, vector<string> *row, char tagOpen, char tagClose, int add)
{
  // set default return value
  string ret = "";
  //ret = tagOpen;
  //ret += tag;
  //ret += tagClose;

  // find conditional tags
  vector<string> auxM;
  StringExplode(tag, auxM, false, "?");
  // check if we have conditional tags. If so, use the form a = b ? c : d
  if(auxM.size() == 2) {
    // split the condition by the equal sign
    vector<string> ab, cd;
    StringExplode(auxM[0], ab, false, "=");
    // split the option by the : sign
    StringExplode(auxM[1], cd, false, ":");
    // check if we have what we need
    if(ab.size() == 2 && cd.size() == 2) {
      string a = translateTag(ab[0], row, add, true) ;
      string b = translateTag(ab[1], row, add, true);
      string c = translateTag(cd[0], row, add, true);
      string d = translateTag(cd[1], row, add, true);
      // check if the condition is true
      if(a == b)
        return c;
      return d;
    }
  }

  // find multi-tag contents
  vector<string> aux;
  StringExplode(tag, aux, false, "|");

  for(int i = 0 ; i < aux.size() ; i++) {
    // translate the tag contents
    string content = translateTag(aux[i], row, add);
    // check for non-empty contents, and stop earching if not empty
    if(content.length())
      return content;
  }

  // return the translated attribute
  return ret;
}

///////////////////////////////////////////////////////////////////////////////
string ReportRendererBase::translateTag(const string &tag, vector<string> *row, int add, bool preserveEmpty)
{
  // set default return value
  string ret = (preserveEmpty ? tag : "");
  //ret = tagOpen;
  //ret += tag;
  //ret += tagClose;

  // search for a number
  string seq("0123456789");
  if(tag.find_first_not_of(seq) == string::npos) {
    int aux = getInt(tag.c_str());
    if(aux < row->size())
      ret = (*row)[aux];

    if(add) {
      int aux2 = getInt(ret.c_str());
      aux2 += add;
      ret = parseInt(aux2);
    }
  }

  // seach for row
  if(tag.compare("row") == 0)
    ret = parseInt(currentRow);
  // search for col
  else if(tag.compare("col") == 0)
    ret = parseInt(currentCol);
  // projectdir tag
  else if(tag.compare("projectdir") == 0)
    ret = m_projectDir;
  // projectdir tag
  else if(tag.compare("aas_file") == 0)
    ret = m_aasFilename;

  // return the translated attribute
  return ret;
}
///////////////////////////////////////////////////////////////////////////////
string ReportRendererBase::splitSequence(const string &sequence, const char *sep)
{
  string ret = "";
  int count = 0;

  for(int i = 0 ; i < sequence.length() ; i++) {

    if(sequence[i] == '(')
      while(sequence[i] != ')' && i < sequence.length()) {
        ret += sequence[i];
        count++;
        i++;
      }

    if(sequence[i] == '[')
      while(sequence[i] != ']' && i < sequence.length()) {
        ret += sequence[i];
        count++;
        i++;
      }

    if(i < sequence.length())
      ret += sequence[i];
    count++;


    if(count >= 10) {
      ret += sep;
      count = 0;
    }
  }
  return ret;
}
///////////////////////////////////////////////////////////////////////////////
string ReportRendererBase::splitText(const string &sequence, const char *sep)
{
  string ret = "";
  int count = 0;

  for(int i = 0 ; i < sequence.length() ; i++) {

    ret += sequence[i];
    count++;


    if(count >= 10) {
      ret += sep;
      count = 0;
    }
  }
  return ret;
}
///////////////////////////////////////////////////////////////////////////////
}; //namespace
////////////////////////////////////////////////////////////////////////////////
