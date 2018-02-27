/*
 * MappedContigStatTable.h
 *
 *  Created on: Mar 4, 2011
 *      Author: aguthals
 */

#ifndef MAPPEDCONTIGSTATTABLE_H_
#define MAPPEDCONTIGSTATTABLE_H_

#include <cstring>
#include <string>
#include <list>
#include <map>
#include <set>
#include <vector>

#include "MappedSpecnets.h"
#include "OutputTable.h"
#include "utils.h"

using namespace std;

namespace specnets
{
  class MappedContigStatTable: public OutputTable
  {
  public:
    MappedSpecnets* mapped_sps_proj;

    MappedContigStatTable();

    MappedContigStatTable(MappedSpecnets* _mapped_sps_proj);

    /**
     * Prepares output table with all necessary statistics
     * @return
     */
    void prepareTable(int index);

  };
}

#endif /* MAPPEDCONTIGSTATTABLE_H_ */
