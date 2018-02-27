#include <cstdio>
#include <cstring>
#include "spectrum.h"
#include "inputParams.h"
#include "batch.h"

int main(int argc, char** argv) {
    InputParams params; bool paramsOk;
	if(argc<=1) paramsOk=params.readParams("intersectPairs.params");	else paramsOk=params.readParams(argv[1]);
	if(!paramsOk) {
		cerr << "Error opening parameters file ";
		if(argc<=1) cerr << "intersectPairs.params\n"; else cerr << argv[1] << endl;
		return -1;
	}
	vector<char *> paramStrs;   paramStrs.resize(3);
	paramStrs[0] = "INPUT_IDX_TO_KEEP";
	paramStrs[1] = "INPUT_ALIGNS";
	paramStrs[2] = "OUTPUT_ALIGNS";
	if(!params.confirmParams(paramStrs)) {
		cerr << "ERROR: Parameters file ";
		if(argc==1) cerr<<"findconnected.params"; else cerr<<argv[1];
		cerr << " is incomplete. One of the following is missing: INPUT_IDX_TO_KEEP, INPUT_ALIGNS, OUTPUT_ALIGNS\n";
		return -1;
	}

	char *idxFN = params.getValue("INPUT_IDX_TO_KEEP");
	char *alignsFN  = params.getValue("INPUT_ALIGNS");
	char *alignsOutFN   = params.getValue("OUTPUT_ALIGNS");

	// Load subset of indices to keep
	ifstream inIdx(idxFN, ios::binary);
	int idxCount, idxMax, idxValue;
	vector<bool> idxOk;

    cout << "Processing "<<idxFN<<"... "; cout.flush();
	if (inIdx.eof() or inIdx.fail()) { cerr<<"Error opening "<<idxFN<<endl; return(-1); }
	inIdx>>idxCount>>idxMax;  idxOk.resize(idxMax+1);
	for(int i=0; i<idxCount; i++) { inIdx>>idxValue; idxOk[idxValue-1]=true; }
	inIdx.close();
	cout << "done. Number of acceptable indices is " << idxCount << " out of "<<idxMax<<endl;

	// Load pairs
	vector<Results_ASP> aligns;
    cout << "Processing "<<alignsFN<<"... "; cout.flush();
	if (Load_resultsASP(alignsFN, aligns)==0) { cerr << "Error reading "<<alignsFN<<"!\n"; return -1; }
    else  cout << "Loading aligns complete. Num pairs: " << aligns.size() << "\n";

	// Count number of pairs between spectra with indices to keep
	int okCount = 0;
	for(unsigned int i=0; i<aligns.size(); i++)
		if(idxOk[aligns[i].spec1] and idxOk[aligns[i].spec2]) okCount++;
    cout<<okCount<<" pairs between spectra from the given subset. Writing output file..."; cout.flush();

	// Output all subset pairs
/*	ofstream output(alignsOutFN, ios::binary);  if (output.fail()) { cerr<<"Error opening "<<alignsOutFN<<endl; return(-1); }
//	output<<okCount<<endl;
	for(unsigned int i=0; i<aligns.size(); i++)
		if(idxOk[aligns[i].spec1] and idxOk[aligns[i].spec2])
			output<<i+1<<"\n";
//			output<<aligns[i].spec1+1<<";"<<aligns[i].spec2+1<<";"<<aligns[i].shift1<<";"<<aligns[i].score1<<";"<<aligns[i].score2<<"\n";
	output.close();
*/

	// Output the indices of all subset pairs - binary format
	FILE *fp = fopen(alignsOutFN,"wb");  if ((int)fp==0) { cerr<<"Error opening "<<alignsOutFN<<endl; return(-1); }
	int *index = (int *)malloc(sizeof(int)*okCount);
	fwrite(&okCount,sizeof(int),1,fp);  okCount=0;
	for(unsigned int i=0; i<aligns.size(); i++)
		if(idxOk[aligns[i].spec1] and idxOk[aligns[i].spec2]) index[okCount++]=i+1;
	fwrite(index,sizeof(int),okCount,fp);
	fclose(fp);
	cout<<"done.\n";

	return 0;
}
