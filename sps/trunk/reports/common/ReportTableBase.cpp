////////////////////////////////////////////////////////////////////////////////
#include <fstream>
#include <algorithm>

#include "ReportTableBase.h"

////////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace specnets;

namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
// ReportIterator methods
////////////////////////////////////////////////////////////////////////////////
TableIterator::TableIterator(ReportTableBase &r)
{
  // keep iterated class pointer
  m_reportTableBase = &r;
  // add iterator pointer to iterated object iterators list
  m_reportTableBase->m_allIterators.push_back(this);
  // define iterator to use, and initialize it
  if(!m_reportTableBase->m_filtered) {
    index = 0;
    m_tableIterator == m_reportTableBase->m_cells.begin();
  } else {
    index = -1;
    m_rowIterator = m_reportTableBase->m_filteredTable.begin();
  }
}
////////////////////////////////////////////////////////////////////////////////
TableIterator::~TableIterator()
{
  // remove iterator from the iterated object iterators list
  m_reportTableBase->m_allIterators.remove(this);
}
////////////////////////////////////////////////////////////////////////////////
TableIterator &TableIterator::begin(void)
{
  // check if the table is being filtered
  if(!m_reportTableBase->m_filtered)
    // if it is not, just return an iterator to the table begin
    m_tableIterator = m_reportTableBase->m_cells.begin();
  else
    // if it is, return the iterator to 1st item after applying the filter
    m_rowIterator = m_reportTableBase->m_filteredTable.begin();
  return *this;
}
////////////////////////////////////////////////////////////////////////////////
TableIterator &TableIterator::end(void)
{
  // check if the table is being filtered
  if(!m_reportTableBase->m_filtered)
    // if it is not, just return an iterator to the table end
    m_tableIterator = m_reportTableBase->m_cells.end();
  else
    // if it is, return the iterator to last item after applying the filter
    m_rowIterator = m_reportTableBase->m_filteredTable.end();
  return *this;
}
////////////////////////////////////////////////////////////////////////////////
void TableIterator::moveTo(int count)
{
  // move to the 1st item
  begin();

  // cycle to the ith item. This needs to be done in case the table is being filtered, and navigation needs to be done through the list
  for(int i = 0 ; i < count ; i++)
    (*this)++;
}
////////////////////////////////////////////////////////////////////////////////
TableIterator TableIterator::operator++(int i)
{
  // The ++ operator over the iterator will depend on if the table is being filtered.
  if(!m_reportTableBase->m_filtered)
    m_tableIterator++;
  else
    m_rowIterator++;
  return *this;
}
////////////////////////////////////////////////////////////////////////////////
TableIterator TableIterator::operator--(int i)
{
  // The -- operator over the iterator will depend on if the table is being filtered.
  if(!m_reportTableBase->m_filtered)
    m_tableIterator--;
  else
    m_rowIterator--;
  return *this;
}
////////////////////////////////////////////////////////////////////////////////
TableIterator TableIterator::operator+=(int i)
{
  // The += operator over the iterator will depend on if the table is being filtered.
  if(!m_reportTableBase->m_filtered) {
    if(m_tableIterator + i >= m_reportTableBase->m_cells.end())
      m_tableIterator = m_reportTableBase->m_cells.end();
    else if(m_tableIterator + i < m_reportTableBase->m_cells.begin())
      m_tableIterator = m_reportTableBase->m_cells.begin();
    else
      m_tableIterator += i;
  } else {
    if(i > 0)
      for(int j = 0 ; j < i && m_rowIterator != m_reportTableBase->m_filteredTable.end() ; j++)
        m_rowIterator++;
    else if(i < 0)
      for(int j = 0 ; j > i && m_rowIterator != m_reportTableBase->m_filteredTable.begin() ; j--)
        m_rowIterator--;
  }
  return *this;
}
////////////////////////////////////////////////////////////////////////////////
bool TableIterator::operator==(const TableIterator &o) const
{
  // different fltered list sizes means the iterators are different
  if(o.m_reportTableBase->m_filteredTable.size() != m_reportTableBase->m_filteredTable.size())
    return false;

  // if list has size greater than 0, then compare row iterators
  if(m_reportTableBase->m_filtered)
    return m_rowIterator == o.m_rowIterator;

  // else compare table iterators
  return m_tableIterator == o.m_tableIterator;
}
////////////////////////////////////////////////////////////////////////////////
bool TableIterator::operator!=(const TableIterator &o) const
{
  // different fltered list sizes means the iterators are different
  if(o.m_reportTableBase->m_filteredTable.size() != m_reportTableBase->m_filteredTable.size())
    return true;

  // if list has size greater than 0, then compare row iterators
  if(m_reportTableBase->m_filtered)
    return m_rowIterator != o.m_rowIterator;

  // else compare table iterators
  return m_tableIterator != o.m_tableIterator;
}
////////////////////////////////////////////////////////////////////////////////
vector<string> *TableIterator::operator*()
{
  if(!m_reportTableBase->m_filtered)
    return &(*m_tableIterator);
  return *m_rowIterator;
}
////////////////////////////////////////////////////////////////////////////////
vector<string> *TableIterator::operator->()
{
  if(!m_reportTableBase->m_filtered)
    return &(*m_tableIterator);
  return *m_rowIterator;
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
// table sorter methods
////////////////////////////////////////////////////////////////////////////////
bool ReportTableSortHandler::operator()(const vector<string> &a, const vector<string> &b)
{
  //if(a.size() <= m_sortColumn || b.size() <= m_sortColumn) return false;

  //string sa = a[m_sortColumn];
  //string sb = b[m_sortColumn];

  int ai = getInt(a[m_sortColumn].c_str());
  int bi = getInt(b[m_sortColumn].c_str());

  if(m_sortDirection)
    return bi < ai;

  return ai < bi;
  //if(sa.length() < sb.length())
  //  return true;
  //return (sa.compare(sb) < 0);
}



////////////////////////////////////////////////////////////////////////////////
// Report table base methods
////////////////////////////////////////////////////////////////////////////////
ReportTableBase::ReportTableBase()
: m_renderHeader(true),
  m_renderBorders(true),
  m_direction(0),
  m_idColumn(-1),
  m_sortColumn(-1),
  m_renderingException(RENDERING_EXCEPTION_NONE),
  m_filtered(false)
{
  m_colTypes.resize(0);
}
////////////////////////////////////////////////////////////////////////////////
ReportTableBase::ReportTableBase(const string &projectPath, const string &tableFilename)
: m_renderHeader(true),
  m_renderBorders(true),
  m_tableFilename(tableFilename),
  m_projectDir(projectPath),
  m_direction(0),
  m_idColumn(-1),
  m_sortColumn(-1),
  m_renderingException(RENDERING_EXCEPTION_NONE),
  m_filtered(false)
{
  m_colTypes.resize(0);
}
////////////////////////////////////////////////////////////////////////////////
ReportTableBase::~ReportTableBase()
{
  clearView();
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableBase::clearView(void)
{
  // Delete all colType objects
  for(int i = 0; i < m_colTypes.size() ; i++)
    if(m_colTypes[i])
      delete m_colTypes[i];
  // set size to zero
  m_colTypes.clear();
}
////////////////////////////////////////////////////////////////////////////////
ReportColumnTypeBase *ReportTableBase::getColType(unsigned x) const
{
  // Check the column. If out-of-bounds, return NULL
  if(x >= m_colTypes.size())
    return NULL;
  // return a pointer to the column type object
  return m_colTypes[x];
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableBase::setValue(unsigned row, unsigned col, string &value)
{
  // check bonds
  if(row >= m_cells.size())
    return ERROR;
  if(col >= m_cells[row].size())
    return ERROR;

  m_cells[row][col] = value;
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
string *ReportTableBase::getCell(unsigned x, unsigned y)
{
  // check row (Y) location. If out-of-bounds, return NULL
  if(y >= m_cells.size())
    return NULL;
  // get the row int an auxiliary variable
  vector<string> &aux = m_cells[y];
  // check the column within the row. If out-of-bounds, return NULL
  if(x >= aux.size())
    return NULL;
  // return a pointer to the string
  return &(aux[x]);
}
////////////////////////////////////////////////////////////////////////////////
string ReportTableBase::getColumnId(int index)
{
  if(index < m_cells.size()) {
    if((m_cells[index].size() > m_idColumn) && (m_idColumn >= 0))
      return m_cells[index][m_idColumn];
  }
  string ret;
  return ret;
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableBase::getFilteredId(vector<string> &IDs, int column)
{
  TableIterator it = begin();
  // get the row
  vector<string> &row = **it;
  IDs.push_back(row[column]);
}
////////////////////////////////////////////////////////////////////////////////
vector<string> ReportTableBase::getColumn(unsigned column)
{
  // define return structure
  vector<string> ret;

  // cycle thru all rows
  for(int i = 0 ; i < m_cells.size() ; i++) {
    // check row (column) location. If out-of-bounds, ignore
    if(column < m_cells[i].size())
      ret.push_back(m_cells[i][column]);
  }
  // return the structure
  return ret;
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableBase::find(int column, string &data)
{
  for(int i = 0 ; i < m_cells.size() ; i++) {
    // test if 'column' is out of bonds for this row
    if(column >= m_cells[i].size())
      return ERROR;
    // check if row satisfyes the condition. If so, return it's index
    if(m_cells[i][column].compare(data) == 0)
      return i;
  }
  return ERROR;
}
////////////////////////////////////////////////////////////////////////////////
// Sort the table
////////////////////////////////////////////////////////////////////////////////
void ReportTableBase::sortTable(void)
{
  ReportTableSortHandler sorter(m_sortColumn, m_direction);
  std::sort(m_cells.begin(), m_cells.end(), sorter);
}
////////////////////////////////////////////////////////////////////////////////
// operations that return an iterator position
TableIterator ReportTableBase::begin(void)
{
  TableIterator tableIterator(*this);
  return tableIterator.begin();
}
////////////////////////////////////////////////////////////////////////////////
TableIterator ReportTableBase::end(void)
{
  TableIterator tableIterator(*this);
  return tableIterator.end();
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableBase::applyFilter(const unsigned column, const string &s)
{
  // clear the current filtered rows list
  m_filteredTable.clear();
  // set filtered flag to false
  m_filtered = false;
  // check for undefined filter column, which means no filtering
  if(column == TABLE_FILTER_COL_NONE)
    return OK;
  // cycle thru all rows. If strings are equal, add to filtered rows list
  for(int i = 0 ; i < m_cells.size() ; i++) {
    // test if 'column' is out of bonds for this row
    if(column >= m_cells[i].size())
      return ERROR;
    // check if row satisfyes the condition. If so, add it to the list
    if(m_cells[i][column].compare(s) == 0)
      m_filteredTable.push_back(&(m_cells[i]));
  }

  // all existing iterators should point to "end" if filter changed
  list<TableIterator *>::iterator it;
  for(it = m_allIterators.begin() ; it != m_allIterators.end() ; it++)
    *(*it) = end();

  // set filtered flag to 'filtered'
  m_filtered = true;

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Table element count (when filtered)
int ReportTableBase::getElementCount(void)
{
  if(m_filtered)
    return m_filteredTable.size();
  else
    return m_cells.size();
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableBase::getPageSplitStats(int itemsPerPage, vector<string> &IDs, vector<string> &IDsEnd, vector<string> &fNames, string &fnamePrefix, string &fnameSuffix)
{
  // calculate # of pages
  int elems = getElementCount();
  int nPages = elems / itemsPerPage;
  int elemsLeft = elems % itemsPerPage;
  if(elemsLeft > 0)
    nPages++;

  int i = 0;
  for(TableIterator it = begin() ; it != end() ; i++) {

    string fname = fnamePrefix;
    fname += parseInt(i);
    fname += fnameSuffix;
    fNames.push_back(fname);

    string start;
    if(m_idColumn < 0)
      start = parseInt(i);
    else {
      const vector<string> &aux = *(*it);
      start = parseInt(getInt(aux[m_idColumn].c_str()));
    }

    IDs.push_back(start);

    for(int j = 0 ; it != end() && j < itemsPerPage ; j++ , it++);

    it--;
    string end;
    if(m_idColumn < 0)
      end = parseInt(i);
    else {
      const vector<string> &aux = *(*it);
      end = parseInt(getInt(aux[m_idColumn].c_str()));
    }
    IDsEnd.push_back(end);
    it++;

  }

  return nPages;
}
////////////////////////////////////////////////////////////////////////////////
// Table loading and saving methods
////////////////////////////////////////////////////////////////////////////////
int ReportTableBase::loadTable(char separator)
{
  // make sure the table is empty;
  m_cells.empty();

  // get filename with full path
  string fn = composeFileName(m_projectDir, m_tableFilename);
  // open file to read
  ifstream file(fn.c_str(), ios::in | ios::binary);
    // if error, say so
  if(!file.is_open()) {
    //cerr << "ERROR: Report table: could not open file " << filename << endl;
    return ERROR;
  }

  // delimiter used to parse the file
  string delim;
  delim = separator;
  // parse the file
  string aux;
  // mark the first line
  bool firstLine = true;
  // get one line from file
  while(getline(file, aux)) {

    // remove comments
    //int found = aux.find_last_of("#");
    //if(found != string::npos)
    //  aux = str.substr(0, found);
    if(firstLine) {
      // vector to store split cells
      vector<string> headers;
      // split the line
      stringSplit2(aux, headers, ";");
      // add the strings to the table
      for(int i = 0 ; i < headers.size() ; i++)
        m_colHeadings.push_back(headers[i]);
      // set first line to false
      firstLine = false;
      // resume to following lines
      continue;
    }

    // vector to store split cells
    vector<string> parts;
    // split the line
    stringSplit2(aux, parts, delim);
    // add the strings to the table
    m_cells.push_back(parts);
  }

  // close the file
  file.close();

  // return OK
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Saves the table in text format
int ReportTableBase::saveTable(char separator, bool noEmptyCells)
{
  // get filename with full path
  string fn = composeFileName(m_projectDir, m_tableFilename);
  // open file to write to
  ofstream file(fn.c_str(), ios::out | ios::binary);
    // if error, say so
  if(!file.is_open()) {
    //cerr << "ERROR: Report table: could not open file for writing: " << filename << endl;
    return ERROR;
  }

  // output file

  //output the table headings, if they exist
  if(m_colHeadings.size()) {
    string line;
    for(int j = 0 ; j < m_colHeadings.size() ; j++) {
      if(j) line += separator;

      string cell = m_colHeadings[j];
      if(noEmptyCells && (cell.length() == 0) )
        cell = " ";

      line += cell; //m_colHeadings[j];
    }
    file << line << endl;
  }

  writeTable(file, noEmptyCells, separator);

  // close the file
  file.close();
  // return OK status
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableBase::writeTable(ostream &outstream, bool noEmptyCells, char separator, char endline)
{
  // output the table, seprated by FILE_SEPARATOR, one line per row
  for(int i = 0 ; i < m_cells.size() ; i++) {
    string line;
    if(i)
      line += endline;
    for(int j = 0 ; j < m_cells[i].size() ; j++) {
      if(j) line += separator;

      string cell = m_cells[i][j];
      if(noEmptyCells && (cell.length() == 0) )
        cell = " ";

      line += cell; //m_cells[i][j];
    }
    outstream << line;
  }
  // return OK status
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableBase::update(const ReportTableData &reportTableData)
{
  // load table
  if(loadTable() == ERROR)
    return ERROR;

  // check boundaries
//  if(reportTableData.updateRow < 0 || reportTableData.updateRow >= m_cells.size())
//    return ERROR;
//  if(reportTableData.updateCol < 0 || reportTableData.updateCol >= m_cells[reportTableData.updateRow].size())
//    return ERROR;

  // change data
//  m_cells[reportTableData.updateRow][reportTableData.updateCol] = reportTableData.updateText;

  // filter table, for a subset of rows
  applyFilter(reportTableData.filterCol, reportTableData.filterText);

  // iterate data
  TableIterator ti = begin();
  // cycle thru all rows
  for(; ti != end() ; ti++) {
    // get the row
    vector<string> &row = **ti;

    //cout << "Old value: " <<  row[reportTableData.updateCol] << endl;

    // process the row
    if(row.size() > reportTableData.updateCol)
      row[reportTableData.updateCol] = reportTableData.updateText;

    //cout << "New value: " <<  row[reportTableData.updateCol] << endl;
  }

  // save the table for future use
  return saveTable();
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableBase::getData(ostream &ss, const ReportTableData &reportTableData, char fieldSep, char rowSep)
{
  // load table
  if(loadTable() == ERROR)
    return ERROR;

  m_sortColumn = reportTableData.sortColumn;
  m_direction  = reportTableData.direction;

  // sort the table, if needed
  if(m_sortColumn >= 0 &&  m_cells.size() && m_sortColumn < m_cells[0].size() )
    sortTable();

  // filter table, for a subset of rows
  applyFilter(reportTableData.filterCol, reportTableData.filterText);

  // iterate data
  bool first = true;
  TableIterator ti = begin();
  // move to desired item
  ti += reportTableData.startRow;
  int rows = reportTableData.rows;
  // cycle thru all rows
  for(; ti != end() && (rows > 0 || rows == -1); ti++, (rows != -1 ? rows-- : rows) ) {
    // put row separator if not the first row
    if(!first)
      ss << rowSep;
    // get the row
    vector<string> &row = **ti;

    // process the row, for the required items
    for(int i = 0 ; i < row.size() ; i++) {
      // output col separator if not the first item
      if(i) ss << fieldSep;
      // output the cell contents
      ss << row[i];
    }
    // tell it's not the first row
    first = false;
  }
}
////////////////////////////////////////////////////////////////////////////////
// File name composition
string ReportTableBase::composeFileName(const string &projectDir, const string &fileName)
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
// Compare two tables contents
int ReportTableBase::compare(ReportTableBase &o, vector<pair<int, int> > &res)
{
  int size1 = m_cells.size();
  int size2 = o.m_cells.size();

  // check for different sizes
  if(size1 != size2) {
    res.push_back(make_pair<int,int>(size1, size2));
    return -1;
  }

  // check for empty tables
  if(!size1)
    return -2;

  int width1 = m_cells[0].size();
  int width2 = o.m_cells[0].size();

  // check for different rows sizes
  if(width1 != width2) {
    res.push_back(make_pair<int,int>(width1, width2));
    return -3;
  }

  // check for empty rows
  if(!width1)
    return -4;

  // check table contents
  for(int i = 0 ; i < size1 ; i++)
    // check rows
    for(int j = 0 ; j < width1 ; j++)
      // check contents
      if(m_cells[i][j] != o.m_cells[i][j])
        // store differences
        res.push_back(make_pair<int,int>(i, j));

  // return number of found differences
  return res.size();
}
////////////////////////////////////////////////////////////////////////////////
}; // namespace
////////////////////////////////////////////////////////////////////////////////
