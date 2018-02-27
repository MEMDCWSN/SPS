////////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_TABLE_BASE_H__
#define __REPORT_TABLE_BASE_H__
////////////////////////////////////////////////////////////////////////////////
#include <iterator>
#include <vector>

#include "ReportColumnTypes.h"
#include "spectrum.h"

////////////////////////////////////////////////////////////////////////////////
// Constants used by this class

#define OK      1
#define ERROR  -1

#define FILE_SEPARATOR      ';'
#define FILE_SEPARATOR_TSV  '\t'


#define DEFAULT_ROOT_DIRECTORY  "."

#define TABLE_FILTER_COL_NONE     -1


#define RENDERING_EXCEPTION_NONE                        0
#define RENDERING_EXCEPTION_PROTEIN_HEADER              1
#define RENDERING_EXCEPTION_MAIN_PAGE                   2
#define RENDERING_EXCEPTION_PROTEIN_COVERAGE_PAGE       3
#define RENDERING_EXCEPTION_PROTEIN_COVERAGE_CSV_PAGE   4

////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
 /*! \brief ReportTableData class

   Helper class used as container to transport data used in update and filtering operations.

   */
class ReportTableData {

 public:

  // Update data -- used for updating a table
  /*! \brief Column to apply the filter on
   */
  int     filterCol;
  /*! \brief Filter string
   */
  string  filterText;

  /*! \brief Column to be updated
   */
  int     updateCol;
  /*! \brief Text to update column with
   */
  string  updateText;

  /*! \brief Column index to sort by
   */
  int     sortColumn;
  /*! \brief Sort direction
   */
  int     direction;

  /*! \brief Row to start on when applying an operation
   */
  int     startRow;
  /*! \brief Number of rows to apply an operation on
   */
  int     rows;

  // constructor
  /*! \brief Default constructor. Sets default values.
   */
  ReportTableData() :
      filterCol(-1),
      updateCol(-1),
      sortColumn(-1),
      direction(0),
      startRow(0),
      rows(-1)
  {};

};

////////////////////////////////////////////////////////////////////////////////
// We need to declare class ReportTableBase before declaring TableIterator
class ReportTableBase;
////////////////////////////////////////////////////////////////////////////////
typedef vector<ReportColumnTypeBase *>::iterator   TableCtIterator;
////////////////////////////////////////////////////////////////////////////////
/*
 *  The iterator class knows the internals of the ReportTableBase, so that it
 *  may move from one element to the next.
 */
 /*! \brief Report Table Iterator class

   Defines a Report Table Class Interator, used to iterate through a set of selected table rows.

   */
class TableIterator {

  // if we are iterating the list, we have the pointer to the current element on the list
  /*! \brief Pointer to the current element on the list
  */
  list<vector<string> *>::iterator m_rowIterator;
  // If not, we use the table iterator
  /*! \brief Table iterator
  */
  vector<vector<string> >::iterator m_tableIterator;
  // if the list is empty, we are iterating the table. we keep the index of the current position
  /*! \brief Current table index, when iteratinh the table
  */
  int index;

  /*! \brief  pointer to itererated class
  */
  ReportTableBase *m_reportTableBase;

 public:

  // constructor and destructor
  TableIterator(ReportTableBase &x);        // --> adds iterator to iterator list in base class
  ~TableIterator();                         // --> removes iterator from iterator list in ReportTable class

  // Move the iterator to specific locations
  /*! \brief Move the iterator to the begin
  */
  TableIterator &begin(void);
  /*! \brief Move the iterator to the end
  */
  TableIterator &end(void);

  /*! \brief Move the iterator to a specified position
  */
  void moveTo(int);

  // increment and decrement operators
  /*! \brief Increment operator
  */
  TableIterator operator++(int);
  /*! \brief Decrement operator
  */
  TableIterator operator--(int);

  /*! \brief += operator
  */
  TableIterator operator+=(int);
  /*! \brief -= operator
  */
  TableIterator operator-=(int);

  // Comparison operators
  /*! \brief Comparison operator: ==
  */
  bool operator==(const TableIterator&) const;
  /*! \brief Comparison operator: !=
  */
  bool operator!=(const TableIterator&) const;

  // Data access operators
  /*! \brief Data access operator
  */
  vector<string> *operator*();
  /*! \brief Data access operator
  */
  vector<string> *operator->();

};
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
 /*! \brief Report Table Sort Handles

   Helper class used to sort tables.

   */
class ReportTableSortHandler {

  /*! \brief Column to sort
   */
  int m_sortColumn;

  /*! \brief Direction to sort
   */
  int m_sortDirection;


 public:

  /*! \brief Constructor
   */
  ReportTableSortHandler(int col, int dir) : m_sortColumn(col), m_sortDirection(dir) {};
  /*! \brief Destructor
   */
  ~ReportTableSortHandler() {};

  /*! \brief The sort method
   */
  bool operator()(const vector<string> &a, const vector<string> &b);

};
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
 /*! \brief Report Table Base class

   Defines base structure for report tables classes.

   */
class ReportTableBase {

 protected:

  //////////////////////////////////////////////////////////////////////////////
  // Table and column type holders, and iterator control

  // colHeadings contains the descriptive text that appears in the table files.
  /*! \brief ontains the descriptive text that appears in the table files.
  */
  vector<string>                  m_colHeadings;

  // Table and column type holders
  /*! \brief Table and column type holders
  */
  vector<ReportColumnTypeBase *>  m_colTypes;
  /*! \brief The table cells
  */
  vector<vector<string> >         m_cells;

  // Filtered list, containing a pointer to the filtered table
  /*! \brief Filtered list, containing a pointer to the filtered table
  */
  list<vector<string> *>          m_filteredTable;

  // List of iterators. Used to keep track of interators, in case of a table filter or order change
  /*! \brief List of iterators. Used to keep track of interators, in case of a table filter or order change
  */
  list<TableIterator *>           m_allIterators;

  // table filename
  /*! \brief table filename
  */
  string  m_tableFilename;

  // Directories and data file names
  /*! \brief Directories and data file names
  */
  string  m_projectDir;

  // specifies if header should be rendered. Default is true.
  /*! \brief specifies if header should be rendered. Default is true.
  */
  bool m_renderHeader;

  // specifies if borders should be rendered. Default is true.
  /*! \brief specifies if borders should be rendered. Default is true.
  */
  bool m_renderBorders;

  // filter column
  /*! \brief filter by column index
  */
  unsigned m_filterColumn;

  // ID column - used to to get IDs
  /*! \brief ID column - used to to get IDs
  */
  int m_idColumn;

  // Sort column - used to to sort the table
  /*! \brief Sort column - used to to sort the table
  */
  int m_sortColumn;

  // table direction
  /*! \brief table direction
  */
  int m_direction;

  // table rendering exception. This is used to tell the renderer to use a specific method to render this table
  /*! \brief table rendering exception. This is used to tell the renderer to use a specific method to render this table
  */
  int m_renderingException;

  // specifies if the table is filtered
  /*! \brief  specifies if the table is filtered
  */
  bool m_filtered;


  //////////////////////////////////////////////////////////////////////////////
  // Method to change iterators, based on the iterator list, when a change occurs
  /*! \brief Method to change iterators, based on the iterator list, when a change occurs
  */
  int invalidateIterators(void) {};

  //////////////////////////////////////////////////////////////////////////////
  // General filtering method. Filter a column by a string
  /*! \brief General filtering method. Filter a column by a string
  */
  virtual int applyFilter(const unsigned column, const string &);


  //////////////////////////////////////////////////////////////////////////////
  // Methods to build and edit the table

  // filename composition helper method
  /*! \brief filename composition helper method
  */
  virtual string composeFileName(const string &projectDir, const string &fileName);
  // edit the table
  /*! \brief edit the table
  */
  virtual int setValue(unsigned row, unsigned col, string &value);


 public:


  // Constructors and destructor
  /*! \brief Constructor
  */
  ReportTableBase();
  /*! \brief Constructor
  */
  ReportTableBase(const string &projectPath, const string &tableFilename);

  // Destructor deletes colTypes vector
  /*! \brief Destructor deletes colTypes vector
  */
  virtual ~ReportTableBase();

  // add a data row to the data table
  /*! \brief add a data row to the data table
  */
  virtual void addDataRow(vector<string> &row) {m_cells.push_back(row);};

  // find a data row in the table. Use 'column' as match column and data as comparison term. Returns row index or -1
  /*! \brief find a data row in the table. Use 'column' as match column and data as comparison term. Returns row index or -1
  */
  virtual int find(int column, string &data);

  //////////////////////////////////////////////////////////////////////////////
  // Methods to load & save data files

  // Load a table from file, given a filename
  /*! \brief Load a table from file, given a filename
  */
  virtual int loadTable(char separator = FILE_SEPARATOR);
  // Save a table to file, given a filename
  /*! \brief Save a table to file, given a filename
  */
  virtual int saveTable(char separator = FILE_SEPARATOR, bool noEmptyCells = false);
  // write the table to a stream
  /*! \brief write the table to a stream
  */
  virtual int writeTable(ostream &outstream, bool noEmptyCells = false, char separator = FILE_SEPARATOR, char endline = '\n');

  /*! \brief Sets the table filename
  */
  virtual void setTableFilename(const char *fn) {m_tableFilename = fn;};
  /*! \brief Sets the table filename
  */
  virtual void setTableFilename(string &fn)     {m_tableFilename = fn;};

  /*! \brief Sets the project directory
  */
  virtual void setProjectDir(const char *fn) {m_projectDir = fn;};
  /*! \brief Sets the project directory
  */
  virtual void setProjectDir(string &fn)     {m_projectDir = fn;};


  //////////////////////////////////////////////////////////////////////////////
  // Input/output methods

  // update table
  /*! \brief updates a table cell
  */
  virtual int update(const ReportTableData &reportTableData);

  // get data
  /*! \brief gets a table cell data
  */
  virtual int getData(ostream &ss, const ReportTableData &reportTableData, char fieldSep, char rowSep);


  //////////////////////////////////////////////////////////////////////////////
  // Data comparison methods

  // compare the contents of two tables
  /*! \brief compare the contents of two tables
  */
  virtual int compare(ReportTableBase &o, vector<pair<int, int> > &res);


  //////////////////////////////////////////////////////////////////////////////
  // Methods to build views

  // clear view
  /*! \brief clear view
  */
  virtual void clearView(void);
  // default view
  /*! \brief default view
  */
  virtual void defineView(void)  {};
  // alternative view
  /*! \brief alternative view
  */
  virtual void defineView2(void) {};
  // for generating images list
  /*! \brief View for generating images list
  */
  virtual void defineViewImages(void) {};
  // specific view for no contigs mapping
  /*! \brief specific view for no contigs mapping
  */
  virtual void defineViewNoClusters(void) {};

  // add a heading
  /*! \brief adds a heading
  */
  virtual void addHeading(char *h)    {string a(h);m_colHeadings.push_back(a);};
  /*! \brief adds a heading
  */
  virtual void addHeading(string &h)  {m_colHeadings.push_back(h);};

  // method to switch off table header rendering
  /*! \brief method to switch off table header rendering
  */
  virtual void noHeaders(void)        {m_renderHeader = false;};
  // query about headers rendering
  /*! \brief query about headers rendering
  */
  virtual bool doHeaders(void)        {return m_renderHeader;};
  // method to switch off table border rendering
  /*! \brief method to switch off table border rendering
  */
  virtual void noBorders(void)        {m_renderBorders = false;};
  // query about border rendering
  /*! \brief query about border rendering
  */
  virtual bool doBorders(void)        {return m_renderBorders;};

  /*! \brief Sets rendeting direction
  */
  virtual void setDirection(int st)   {m_direction = st;};
  // query about rendeting direction
  /*! \brief query about rendeting direction
  */
  virtual int  getDirection(void)     {return m_direction;};

  // returns the # of elements (filtered)
  /*! \brief  returns the # of elements (filtered)
  */
  virtual int  getElementCount(void);
  virtual int  getPageSplitStats(int itemsPerPage, vector<string> &IDs, vector<string> &IDsEnd, vector<string> &fNames, string &fnamePrefix, string &fnameSuffix);


  // filter related methods. (Re)builds filtered table structure
  /*! \brief filter related methods. (Re)builds filtered table structure
  */
  virtual int  applyFilter(const string &value)   {return applyFilter(m_filterColumn, value);};
  // define column to apply filter to, by index
  /*! \brief define column to apply filter to, by index
  */
  virtual void setFilterColumn(unsigned c)  {m_filterColumn = c;};

  // set the column ID (used for table sorting and to get column id, given it's index)
  /*! \brief set the column ID (used for table sorting and to get column id, given it's index)
  */
  virtual void setIdColumn(int id)          {m_idColumn = id;   };
  /*! \brief Gets the column ID
  */
  virtual string getColumnId(int index);
  // gets ID
  /*! \brief gets table ID
  */
  virtual void getId(vector<string> &IDs) {getFilteredId(IDs);};

  /*! \brief Get IDs from filtered table. If no column ID is supplyed, the column ID index is used.
  */
  virtual void getFilteredId(vector<string> &IDs) {getFilteredId(IDs, m_idColumn);};
  /*! \brief Get IDs from filtered table, given a column
  */
  virtual void getFilteredId(vector<string> &IDs, int column);


  // set/get the table rendering exception flag
  /*! \brief set/get the table rendering exception flag
  */
  virtual void setRenderingException(int v) {m_renderingException = v;    };
  /*! \brief set/get the table rendering exception flag
  */
  virtual int  getRenderingException(void)  {return m_renderingException; };

  // sort the table
  /*! \brief sorts the table
  */
  virtual void sortTable(void);
  // sort the table
  /*! \brief sorts the table
  */
  virtual void sortTable(int col)           {m_sortColumn=col; sortTable();};

  // Retrieve a columType from vector
  /*! \brief Retrieves a columType from vector
  */
  virtual ReportColumnTypeBase *getColType(unsigned x) const;
  // Retrieve a cell from table
  /*! \brief Retrieves a cell from table
  */
  virtual string *getCell(unsigned x, unsigned y);
  // get a list of elements of a column for the entire table. This is usefull for iteration when the table is also being filtered for other purposes
  /*! \brief get a list of elements of a column for the entire table.
  */
  virtual vector<string> getColumn(unsigned column);

  // operations that return an iterator position for the table contents
  /*! \brief operations that return an iterator position for the table contents
  */
  virtual TableIterator begin(void);
  virtual TableIterator end(void);
  /*! \brief
  */
  // operations that return an iterator position for the column types vector
  /*! \brief operations that return an iterator position for the column types vector
  */
  virtual TableCtIterator ctBegin() {return m_colTypes.begin();};
  /*! \brief operations that return an iterator position for the column types vector
  */
  virtual TableCtIterator ctEnd()   {return m_colTypes.end();  };

  // must make iterator class friend
  friend class TableIterator;

};
////////////////////////////////////////////////////////////////////////////////
}; // namespace
////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////
