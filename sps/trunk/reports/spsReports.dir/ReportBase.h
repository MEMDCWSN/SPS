///////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_BASE_H__
#define __REPORT_BASE_H__
///////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>

#include "ReportTableBase.h"

////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
// Report Elements
////////////////////////////////////////////////////////////////////////////////
  /*! \brief Report Elements: base class
   */
class ReportElementsBase {
 public:

  ReportElementsBase() {};
  virtual ~ ReportElementsBase() {};
};
/*----------------------------------------------------------------------------*/
  /*! \brief Report Elements: DIV
   */
class ReportElementsDiv : public ReportElementsBase {
 public:

  vector<ReportElementsBase *>  items;

  virtual ~ReportElementsDiv() {for(int i=0;i<items.size();i++) delete items[i];};
};
/*----------------------------------------------------------------------------*/
  /*! \brief Report Elements: string
   */
class ReportElementsString : public ReportElementsBase {
 public:

  string  str;
};
/*----------------------------------------------------------------------------*/
  /*! \brief Report Elements: table cell
   */
class ReportElementsCell : public ReportElementsBase {
 public:

  vector<ReportElementsBase *>  data;

  virtual ~ ReportElementsCell() {for(int i=0;i<data.size();i++) delete data[i];};
};
/*----------------------------------------------------------------------------*/
  /*! \brief Report Elements: table row
   */
class ReportElementsRow : public ReportElementsBase {
 public:

  vector<ReportElementsCell *>  cells;

  virtual ~ ReportElementsRow() {for(int i=0;i<cells.size();i++) delete cells[i];};
};
/*----------------------------------------------------------------------------*/
  /*! \brief Report Elements: table
   */
class ReportElementsTable : public ReportElementsBase {
 public:

  vector<ReportElementsRow *>  rows;

  virtual ~ ReportElementsTable() {for(int i=0;i<rows.size();i++) delete rows[i];};
};
////////////////////////////////////////////////////////////////////////////////
/*
 *  The iterator class knows the internals of the ReportTableBase, so that it
 *  may move from one element to the next.
 */
typedef  vector<ReportTableBase *>::iterator ReportIterator;

///////////////////////////////////////////////////////////////////////////////
  /*! \brief Report base class
   */
class ReportBase {

 protected:

  // where to store report headers.
  /*! \brief where to store report headers.
   */
  vector<ReportElementsBase *> elements;

  // where report tables are stored
  /*! \brief where report tables are stored
   */
  vector<ReportTableBase *> table;



  // adds a table to the report
  /*! \brief adds a table to the report
   */
  virtual int addTable(ReportTableBase *t)  {table.push_back(t);};


 public:

  // Constructors and destructor
  /*! \brief Constructors
   */
  ReportBase();
  /*! \brief destructor
   */
  virtual ~ReportBase();

  // load data
  /*! \brief load data
   */
  virtual int load(void);
  // update table
  /*! \brief update table
   */
  virtual int update(const ReportTableData &reportTableData);

  virtual void sortTables(void);

  // Apply a filter
  /*! \brief Apply a filter
   */
  virtual int applyFilter(string &value);
  /*! \brief Apply a filter
   */
  virtual int applyFilter(unsigned tableIdx, string &value);

  // get a specific table by index
  /*! \brief get a specific table by index
   */
  virtual ReportTableBase *getTable(unsigned tableIdx);

  // get a specific column of a specific table
  /*! \brief get a specific column of a specific table
   */
  virtual vector<string> getTableColumn(unsigned tableIdx, unsigned column);

  // operations that return an iterator position
  /*! \brief Returns iterator position: begin
   */
  virtual ReportIterator begin(void) {return table.begin();};
  /*! \brief Returns iterator position: end
   */
  virtual ReportIterator end(void)   {return table.end();  };

  // get the total number of elements for this report
  /*! \brief get the total number of elements for this report
   */
  virtual int getElementCount(void);

  // gets ID
  /*! \brief gets table ID
   */
  virtual void getId(vector<string> &IDs) {if(table.size()) table[0]->getId(IDs);};

  /*! \brief Get page split statistics
   */
  virtual int getPageSplitStats(int elemsPerPage, vector<string> &IDs, vector<string> &IDsEnd, vector<string> &fNames, string &fnamePrefix, string &fnameSuffix);

};
///////////////////////////////////////////////////////////////////////////////
}; // namespace
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
