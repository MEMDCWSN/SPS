#ifndef __REPORT_DEFINES
#define __REPORT_DEFINES

#include <vector>

#define DEFAULT_ANNOTATION_MODEL  "model_cid.txt"
#define ANNOTATION_MODEL_CID      "model_cid.txt"
#define ANNOTATION_MODEL_ETD      "model_etd.txt"
#define ANNOTATION_MODEL_PRM      "model_prm.txt"

////////////////////////////////////////////////////////////////////////////////
typedef std::pair<std::pair<std::vector<int>, std::vector<int> >, std::vector< std::pair<std::vector<int>, std::vector<double> > > > abContig_t;
typedef std::vector< std::pair<std::vector<int>, std::vector<double> > > abContigData_t;
typedef std::pair<std::vector<int>, std::vector<int> > abContigState_t;



#endif

