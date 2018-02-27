#include "aminoacid.h"
#include "batch.h"
#include "db_fasta.h"
#include "inputParams.h"
#include "spectrum.h"

using namespace specnets;

int main(int argc, char **argv) {
	InputParams params; 	vector<const char *> paramStrs;   paramStrs.resize(3);
	paramStrs[0] = "INPUT_FASTA";
	paramStrs[1] = "OUTPUT_FASTA";
	paramStrs[2] = "INPUT_MATCHES_LIST";
//	paramStrs[2] = "INPUT_MATCHES2";
//	paramStrs[3] = "INPUT_MATCHES2";

	if(argc==1)	params.readParams("protid.params"); else params.readParams(argv[1]);
	if(!params.confirmParams(paramStrs)) {
		cerr << "ERROR: Parameters file ";
		if(argc==1) cerr<<"protid.params"; else cerr<<argv[1];
		cerr << " is incomplete. One of the following is missing: INPUT_FASTA, OUTPUT_FASTA, INPUT_MATCHES_LIST\n";
		return -1;
	}
	if(params.paramPresent("AMINO_ACID_MASSES")) { AAJumps jumps(-1); jumps.loadJumps(params.getValue("AMINO_ACID_MASSES"),true); }

	DB_fasta db; if(db.Load(params.getValue("INPUT_FASTA"))<=0) return -1;
    vector<bool> matchedProts(db.size());   for(unsigned int pIdx=0; pIdx<matchedProts.size(); pIdx++) matchedProts[pIdx]=false;

    // Load set of tagsearch results
	vector<list<pair<int,float> > > simpleProteinHits, curProteinHits;
	unsigned int pivot1, pivot2, lineIdx;
/*	if(Load_binListArray<int,list<int>,list<int>::iterator>(params.getValue("INPUT_MATCHES1"), simpleProteinHits)<0) { cerr<<"Error reading "<<params.getValue("INPUT_MATCHES1")<<"!\n"; exit(-1); }
	if(Load_binListArray<int,list<int>,list<int>::iterator>(params.getValue("INPUT_MATCHES2"), curProteinHits)<0) { cerr<<"Error reading "<<params.getValue("INPUT_MATCHES2")<<"!\n"; exit(-1); }
	unsigned int pivot1=simpleProteinHits.size(), pivot2;
	simpleProteinHits.resize(simpleProteinHits.size()+curProteinHits.size());
	for(pivot2=0; pivot2<curProteinHits.size(); pivot2++, pivot1++)
		simpleProteinHits[pivot1].swap(curProteinHits[pivot2]);
*/
	BufferedLineReader blr;
	if(blr.Load(params.getValue("INPUT_MATCHES_LIST"))<=0) {
	  cerr<<"ERROR reading "<<params.getValue("INPUT_MATCHES_LIST")<<"!\n";
	  return -1;
	}
	if(Load_binListArray<int,float,list<pair<int,float> >,list<pair<int,float> >::iterator>(blr.getline(0), simpleProteinHits)<0) {
	  cerr<<"Error reading "<<blr.getline(0)<<"!\n";
	  return -1;
	}
	for(lineIdx=1; lineIdx<blr.size(); lineIdx++) {
		if (strlen(blr.getline(lineIdx))>0) {
			if(Load_binListArray<int,float,list<pair<int,float> >,list<pair<int,float> >::iterator>(blr.getline(lineIdx), curProteinHits)<0) {
			  cerr<<"Error reading "<<blr.getline(lineIdx)<<"!\n";
			  return -1;
			}
			pivot1=simpleProteinHits.size();
			simpleProteinHits.resize(simpleProteinHits.size()+curProteinHits.size());
			for(pivot2=0; pivot2<curProteinHits.size(); pivot2++, pivot1++)
				simpleProteinHits[pivot1].swap(curProteinHits[pivot2]);
		}
  }

	MaximumParsimony(simpleProteinHits);

	for(pivot1=0; pivot1<simpleProteinHits.size(); pivot1++) {
		if(not simpleProteinHits[pivot1].empty()) {
		  matchedProts[simpleProteinHits[pivot1].front().first]=true;
		}
  }

	FILE *protsFID = (FILE *)fopen(params.getValue("OUTPUT_FASTA"), "wb");
	for (unsigned int pIdx=0; pIdx<matchedProts.size(); pIdx++)
		if (matchedProts[pIdx])
			fprintf(protsFID, ">%s %s\n%s\n", db.getID(pIdx),db.getDesc(pIdx), db[pIdx]);
	fclose(protsFID);

	return 0;

}
