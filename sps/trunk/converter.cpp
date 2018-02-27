#include <cstdio>
#include <cstring>
#include "spectrum.h"
#include "inputParams.h"
#include "batch.h"

int main(int argc, char** argv) {
    InputParams params; bool paramsOk;
	if(argc<=1) paramsOk=params.readParams("converter.params");	else paramsOk=params.readParams(argv[1]);
	if(!paramsOk) {
		cerr << "Error opening parameters file ";
		if(argc<=1) cerr << "converter.params\n"; else cerr << argv[1] << endl;
		return -1;
	}

	char *specSetFN = params.getValue("INPUT_SPECS");
	char *alignsFN  = params.getValue("INPUT_ALIGNS");
	char *twocol1FN = params.getValue("INPUT_TWOCOL1");
	char *twocol2FN = params.getValue("INPUT_TWOCOL2");
	int twocolCount1 = params.getValueInt("COUNT_TWOCOL1");
	int twocolCount2 = params.getValueInt("COUNT_TWOCOL2");
	char *outputSpecsFN = params.getValue("OUTPUT_PKLBIN");
	char *alignsOutFN   = params.getValue("OUTPUT_ALIGNS");
	char *twocolOut1FN  = params.getValue("OUTPUT_TWOCOL1");
	char *twocolOut2FN  = params.getValue("OUTPUT_TWOCOL2");

	if(strlen(specSetFN)>0 and strlen(outputSpecsFN)>0) {
	    SpecSet specSet;
	    cout << "Processing "<<specSetFN<<"... "; cout.flush();
	    if (specSet.LoadSpecSet_pkl(specSetFN)<0) { cerr << "Error reading "<<specSetFN<<"!\n"; return -1; }
	    else  cout << "Loading specs complete. Num specs: " << specSet.size() << "\n";

	    if (specSet.SaveSpecSet_pklbin(outputSpecsFN)<0) { cerr << "Error writing "<<outputSpecsFN<<"!\n"; return -1; }
	}

	if(strlen(alignsFN)>0 and strlen(alignsOutFN)>0) {
		vector<Results_ASP> aligns;
	    cout << "Processing "<<alignsFN<<"... "; cout.flush();
		if (Load_resultsASP(alignsFN, aligns)==0) { cerr << "Error reading "<<alignsFN<<"!\n"; return -1; }
	    else  cout << "Loading aligns complete. Num pairs: " << aligns.size() << "\n";

		// Compensate for Load_resultsASP subtracting 1 from every spectrum index. Conversion doesn't change the contents of the files.
	    for(unsigned int i=0; i<aligns.size(); i++) { aligns[i].spec1++; aligns[i].spec2++; }
	    if (Save_resultsASPbin(alignsOutFN, aligns)<0) { cerr << "Error writing "<<alignsOutFN<<"!\n"; return -1; }
	}

	if(strlen(twocol1FN)>0 and strlen(twocolOut1FN)>0 and twocolCount1>0) {
	    cout << "Processing "<<twocol1FN<<"... "; cout.flush();
		vector<TwoValues<float> > values(twocolCount1);
		unsigned int curEntry=0;
		ifstream ifs(twocol1FN, ios::binary);
		while (!ifs.eof()) { ifs>>values[curEntry][0]>>values[curEntry][1]; curEntry++; }
		ifs.close();
		Utils::save_tcb(twocolOut1FN,values);
		cout<<"done.\n";
	}

	if(strlen(twocol2FN)>0 and strlen(twocolOut2FN)>0 and twocolCount2>0) {
	    cout << "Processing "<<twocol2FN<<"... "; cout.flush();
		vector<TwoValues<float> > values(twocolCount2);
		unsigned int curEntry=0;
		ifstream ifs(twocol2FN, ios::binary);
		while (!ifs.eof()) { ifs>>values[curEntry][0]>>values[curEntry][1]; curEntry++; }
		ifs.close();
		Utils::save_tcb(twocolOut2FN,values);
		cout<<"done.\n";
	}

	return 0;
}
