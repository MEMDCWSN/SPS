/*
 * MappedContigSetTable.h
 *
 *  Created on: Dec 7, 2012
 *      Author: aguthals
 */

#ifndef MAPPEDCONTIGSETTABLE_H_
#define MAPPEDCONTIGSETTABLE_H_

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
  class MappedContigSetTable: public OutputTable
  {
  public:
    MappedSpecnets* mapped_sps_proj;

    MappedContigSetTable();

    MappedContigSetTable(MappedSpecnets* _mapped_sps_proj);

    /**
     * Prepares output table with all necessary statistics
     * @return
     */
    void prepareTable(void);

  };
}


#endif /* MAPPEDCONTIGSETTABLE_H_ */
