///////////////////////////////////////////////////////////////////////////////
#include "ReportBase.h"
///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
ReportBase::ReportBase()
{
  table.resize(0);
  elements.resize(0);
}
////////////////////////////////////////////////////////////////////////////////
ReportBase::~ReportBase()
{
  for(int i=0;i<table.size();i++)
    if(table[i])
      delete table[i];

  for(int i=0;i<elements.size();i++)
    if(elements[i])
      delete elements[i];
}
////////////////////////////////////////////////////////////////////////////////
int ReportBase::load(void)
{
  for(int i = 0 ; i < table.size() ; i++)
    table[i]->loadTable();
}
///////////////////////////////////////////////////////////////////////////////
int ReportBase::update(const ReportTableData &reportTableData)
{
}
///////////////////////////////////////////////////////////////////////////////
void ReportBase::sortTables(void)
{
  for(int i = 0 ; i < table.size() ; i++)
    table[i]->sortTable();
}
///////////////////////////////////////////////////////////////////////////////
int ReportBase::applyFilter(string &value)
{
  for(int i = 0 ; i < table.size() ; i++)
    table[i]->applyFilter(value);
  // return code ok
  return OK;
}
///////////////////////////////////////////////////////////////////////////////
int ReportBase::applyFilter(unsigned tableIdx, string &value)
{
  if(tableIdx >= table.size())
    return ERROR;
  return table[tableIdx]->applyFilter(value);
}
///////////////////////////////////////////////////////////////////////////////
ReportTableBase *ReportBase::getTable(unsigned tableIdx)
{
  if(tableIdx >= table.size())
     return NULL;
  return table[tableIdx];
}
///////////////////////////////////////////////////////////////////////////////
vector<string> ReportBase::getTableColumn(unsigned tableIdx, unsigned column)
{
  // get table
  ReportTableBase *aux = getTable(tableIdx);
  // get the row
  if(aux)
    return aux->getColumn(column);

  // return empty results, if table not found
  vector<string> ret;
  return ret;
}
///////////////////////////////////////////////////////////////////////////////
// get the total number of elements for this report
int ReportBase::getElementCount(void)
{
  // declare return varable
  int ret = 0;
  // cycle thru all tables
  for(int i = 0 ; i < table.size() ; i++)
  // for each table, get element count
    ret += table[i]->getElementCount();
  // return value
  return ret;
}
///////////////////////////////////////////////////////////////////////////////
int ReportBase::getPageSplitStats(int elemsPerPage, vector<string> &IDs, vector<string> &IDsEnd, vector<string> &fNames, string &fnamePrefix, string &fnameSuffix)
{
  // only the last page might need pagination
  int i = table.size();
  if(i)
    return table[i-1]->getPageSplitStats(elemsPerPage, IDs, IDsEnd, fNames, fnamePrefix, fnameSuffix);
  return 0;
}
///////////////////////////////////////////////////////////////////////////////
}; //namespace
///////////////////////////////////////////////////////////////////////////////
