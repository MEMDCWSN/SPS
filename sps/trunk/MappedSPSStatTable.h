/*
 * MappedSpecnetsStatTable.h
 *
 *  Created on: Mar 4, 2011
 *      Author: aguthals
 */

#ifndef MAPPEDSPSSTATTABLE_H_
#define MAPPEDSPSSTATTABLE_H_

#include <cstring>
#include <string>
#include <list>
#include <map>
#include <set>
#include <vector>

#include "MappedSpecnets.h"
#include "OutputTable.h"
#include "utils.h"

namespace specnets
{

  class MappedSPSStatTable: public OutputTable
  {
  public:

    MappedSpecnets* mapped_sps_proj;

    MappedSPSStatTable();

    MappedSPSStatTable(MappedSpecnets* _mapped_sps_proj);

    /**
     * Prepares output table with all necessary statistics
     * @return
     */
    virtual void prepareTable(void);
  };
}

#endif /* MAPPEDSPSSTATTABLE_H_ */
