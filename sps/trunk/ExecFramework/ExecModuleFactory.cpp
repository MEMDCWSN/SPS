// Module Includes
#include "Logger.h"
#include "ExecModuleFactory.h"

//----------------------------------------
// Add new module include files here
#include "ExecAlignment.h"
#include "ExecAssembly.h"
#include "ExecCreateSpectralLibrary.h"
#include "ExecDeconvoluteMS2.h"
#include "ExecFdrPeptide.h"
#include "ExecFilterAligns.h"
#include "ExecFilterContigPairs.h"
#include "ExecFilterPairs.h"
#include "ExecFilterStarPairs.h"
#include "ExecGenoMS.h"
//#include "ExecGFNetwork.h"
#include "ExecHomologyAssembly.h"
//#include "ExecMainSpecnets.h"
#include "ExecMergeConvert.h"
#include "ExecMetaAssembly.h"
#include "ExecCCMSMetaSPSSpecnetsParams.h"
#include "ExecMsCluster.h"
#include "ExecParallelAssembly.h"
#include "ExecPrmClustering.h"
#include "ExecPrmScoring.h"
#include "ExecProjection.h"
#include "ExecProjectionStatistics.h"
#include "ExecProteosafeToSpecnetsParams.h"
#include "ExecProtProtAlign.h"
#include "ExecQualityControl.h"
#include "ExecQCSpectrum.h"
//#include "ExecReportProteinCoverage.h"
#include "ExecReportSpsplot.h"
#include "ExecReportSPSStats.h"
//#include "ExecSpecNetsPropagation.h"
#include "ExecSpecProtAlign.h"
#include "ExecSpecTagGen.h"
#include "ExecPenaltyGen.h"
#include "ExecSpectralLibrarySearch.h"
#include "ExecSpectralLibrarySearchMolecular.h"
#include "ExecSpectralLibraryDecoyCreation.h"
#include "ExecSpectralLibrarySLGFTraining.h"
#include "FirstCorrectTags.h"
#include "ExecSpectralLibrarySLGFCreation.h"
#include "ExecSpectraExtraction.h"
#include "ExecSpectralProbability.h"
#include "ExecStatistics.h"
#include "ExecQCSpectralPairs.h"
#include "ExecSvmStatistics.h"
#include "ExecTagSearch.h"
//#include "ExecSpecNetworkEval.h"
#include "ExecSpectralLibrarySearchSLGF.h"
#include "ExecCCMSMetabolomicsSpecnetsParams.h"
#include "ExecLibraryView.h"
#include "ExecSpectraExtractionSingle.h"
#include "ExecSpectraUpdateSingle.h"
#include "ExecCreateSpectralLibraryMolecular.h"
#include "ExecSpectraExtractionTable.h"
#include "ExecFilterSpectra.h"
#include "ExecSpectralLibraryBatchRemove.h"
#include "ExecSpectraSyntheticSpectraScanExtraction.h"
#include "ExecSpectralLibraryPeakStatistics.h"
#include "ExecIntersectingSpecProb.h"
#include "ExecStarGF.h"
#include "ExecTrainMassError.h"
#include "ExecNetworksPropagation.h"
#include "ExecTagFilterPairs.h"
#include "ExecMolecularParallelPairs.h"
#include "ExecKLFilter.h"
#include "ExecGenerateTags.h"
#include "ExecClusterBins.h"

//----------------------------------------

using namespace std;
using namespace specnets;

//----------------------------------------------------------------------------
namespace specnets {
static map<string, ExecBase *>* theMap = 0;

void ExecModuleFactory::cleanup(void) {
	map<string, ExecBase *> & moduleMap = getMap();
	//cout << "cur_size = " << moduleMap.size() << "\n";
	for (map<string, ExecBase *>::iterator modIt = moduleMap.begin();
			modIt != moduleMap.end(); modIt++) {
		//cout << "deleting " << modIt->first << "(" << modIt->second << ")\n";
		delete modIt->second;
	}
	delete theMap;
	theMap = 0;
}
// Registers a singe exec module with the factory
void ExecModuleFactory::Register(std::string moduleName, ExecBase * module) {
	map<string, ExecBase *> & moduleMap = getMap();
	moduleMap[moduleName] = module;
	//cout << "registering " << moduleName << "(" << module << ")\n";
}
;

// Registers all the known modules with the factory
void ExecModuleFactory::RegisterAllModules(void) {
	//--------------------------------------------------------
	// Add new module registrations here
	//--------------------------------------------------------
	Register(ExecAlignment().getName(), new ExecAlignment());
	Register(ExecAssembly().getName(), new ExecAssembly());
	Register(ExecCreateSpectralLibrary().getName(),
			new ExecCreateSpectralLibrary());
	Register(ExecDeconvoluteMS2().getName(), new ExecDeconvoluteMS2());
	Register(ExecFdrPeptide().getName(), new ExecFdrPeptide());
	Register(ExecFilterAligns().getName(), new ExecFilterAligns());
	Register(ExecFilterContigPairs().getName(), new ExecFilterContigPairs());
	Register(ExecFilterPairs().getName(), new ExecFilterPairs());
	Register(ExecFilterStarPairs().getName(), new ExecFilterStarPairs());
	Register(ExecGenoMS().getName(), new ExecGenoMS());
	//Register(ExecGFNetwork().getName(), new ExecGFNetwork());
	Register(ExecHomologyAssembly().getName(), new ExecHomologyAssembly());
	//    Register( ExecMainSpecnets().getName(), new ExecMainSpecnets() );
	Register(ExecMergeConvert().getName(), new ExecMergeConvert());
	Register(ExecMetaAssembly().getName(), new ExecMetaAssembly());
	Register(ExecCCMSMetaSPSSpecnetsParams().getName(),
			new ExecCCMSMetaSPSSpecnetsParams());
	Register(ExecProteosafeToSpecnetsParams().getName(),
			new ExecProteosafeToSpecnetsParams());
	Register(ExecMsCluster().getName(), new ExecMsCluster());
	Register(ExecParallelAssembly().getName(), new ExecParallelAssembly());
	Register(ExecPrmClustering().getName(), new ExecPrmClustering());
	Register(ExecPrmScoring().getName(), new ExecPrmScoring());
	Register(ExecProjection().getName(), new ExecProjection());
	Register(ExecQualityControl().getName(), new ExecQualityControl());
	Register(ExecQCSpectrum().getName(), new ExecQCSpectrum());
	Register(ExecSpectraExtraction().getName(), new ExecSpectraExtraction());
	Register(ExecSpectralLibrarySLGFTraining().getName(),
			new ExecSpectralLibrarySLGFTraining());
    Register(FirstCorrectTags().getName(),
			new FirstCorrectTags());
	Register(ExecSpectralLibrarySLGFCreation().getName(),
			new ExecSpectralLibrarySLGFCreation());
	Register(ExecProjectionStatistics().getName(),
			new ExecProjectionStatistics());
	Register(ExecProtProtAlign().getName(), new ExecProtProtAlign());
	//Register(ExecReportProteinCoverage().getName(), new ExecReportProteinCoverage());
	Register(ExecReportSpsplot().getName(), new ExecReportSpsplot());
	Register(ExecReportSPSStats().getName(), new ExecReportSPSStats());
	//    Register( ExecSpecNetsPropagation().getName(), new ExecSpecNetsPropagation() );
	Register(ExecSpecProtAlign().getName(), new ExecSpecProtAlign());
	Register(ExecSpecTagGen().getName(), new ExecSpecTagGen());
	Register(ExecPenaltyGen().getName(), new ExecPenaltyGen());

	Register(ExecSpectralLibrarySearch().getName(),
			new ExecSpectralLibrarySearch());
	Register(ExecSpectralLibrarySearchMolecular().getName(),
			new ExecSpectralLibrarySearchMolecular());
	Register(ExecStatistics().getName(), new ExecStatistics());
	Register(ExecQCSpectralPairs().getName(), new ExecQCSpectralPairs());
	Register(ExecNetworksPropagation().getName(), new ExecNetworksPropagation());
	Register(ExecTagFilterPairs().getName(), new ExecTagFilterPairs());
	Register(ExecSvmStatistics().getName(), new ExecSvmStatistics());
	Register(ExecTagSearch().getName(), new ExecTagSearch());
	//Register(ExecSpecNetworkEval().getName(), new ExecSpecNetworkEval());
	Register(ExecSpectralLibraryDecoyCreation().getName(),
			new ExecSpectralLibraryDecoyCreation());
	Register(ExecSpectralLibrarySearchSLGF().getName(),
			new ExecSpectralLibrarySearchSLGF());
	Register(ExecCCMSMetabolomicsSpecnetsParams().getName(),
			new ExecCCMSMetabolomicsSpecnetsParams());
	Register(ExecLibraryView().getName(), new ExecLibraryView());
	Register(ExecSpectraExtractionSingle().getName(),
			new ExecSpectraExtractionSingle());
	Register(ExecSpectraUpdateSingle().getName(),
			new ExecSpectraUpdateSingle());
	Register(ExecCreateSpectralLibraryMolecular().getName(),
			new ExecCreateSpectralLibraryMolecular());
	Register(ExecSpectraExtractionTable().getName(),
			new ExecSpectraExtractionTable());
	Register(ExecFilterSpectra().getName(), new ExecFilterSpectra());
	Register(ExecSpectralLibraryBatchRemove().getName(),
			new ExecSpectralLibraryBatchRemove());
	Register(ExecSpectraSyntheticSpectraScanExtraction().getName(),
			new ExecSpectraSyntheticSpectraScanExtraction());
	Register(ExecSpectralLibraryPeakStatistics().getName(),
			new ExecSpectralLibraryPeakStatistics());
	Register(ExecTrainMassError().getName(), new ExecTrainMassError());
	Register(ExecIntersectingSpecProb().getName(),
			new ExecIntersectingSpecProb());
	Register(ExecStarGF().getName(), new ExecStarGF());
	Register(ExecSpectralProbability().getName(),
			new ExecSpectralProbability());
		Register(ExecKLFilter().getName(),
				 new ExecKLFilter());
               Register(ExecGenerateTags().getName(),
				 new ExecGenerateTags());
    Register(ExecMolecularParallelPairs().getName(),
			new ExecMolecularParallelPairs());
			Register(ExecClusterBins().getName(),
					 new ExecClusterBins());
}

// Retrieves a module from the factory by its name
ExecBase * ExecModuleFactory::getModule(std::string moduleName) {
	map<string, ExecBase *> & moduleMap = getMap();
	if (moduleMap.find(moduleName) != moduleMap.end()) {
		return moduleMap[moduleName];
	}

	return 0;
}
;

//----------------------------------------------------------------------------
// NOTE: One might be tempted to remove this method and define the map as
//       a static variable. Don't! C++ does not define when static variables
//       will be instatiated. By making it a method static it is guaranteed
//       to be instantiated when the method is first called.
//----------------------------------------------------------------------------
map<string, ExecBase *> & ExecModuleFactory::getMap(void) {
	if (theMap == 0) {
		theMap = new map<string, ExecBase *>;
	}
	return *theMap;
}

}
;
// namespace specnets
