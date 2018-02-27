#include "inputParams.h"
#include "batch.h"
#include "ms1.h"
#include "mzxml.h"

#include <iostream>
#include <fstream>
#include <cstdio>

using namespace std;

int main(int argc, char **argv){
	InputParams params;

	if(argc==1)	params.readParams("procms1.params"); else params.readParams(argv[1]);
	if(not(params.paramPresent("INPUT_SPECS1") and params.paramPresent("INPUT_SPECS2")) and not(params.paramPresent("INPUT_MZXML1")) ) {
		cerr << "ERROR: Parameters file ";
		if(argc==1) cerr<<"procms1.params"; else cerr<<argv[1];
		cerr << " is incomplete. One of the following is missing: INPUT_SPECS1/INPUT_SPECS2 or INPUT_MZXML1\n";
		return -1;
	}

	float peakTol = params.paramPresent("TOLERANCE_PEAK")?(float) params.getValueDouble("TOLERANCE_PEAK"):0.5;
/*	char *specSetFN = params.getValue(paramStrs[0]);
	char *ms3setsFN = params.getValue(paramStrs[1]);
	char *outputSpecsFN = params.getValue(paramStrs[2]);
	float pmTol = params.paramPresent("TOLERANCE_PM")?(float) params.getValueDouble("TOLERANCE_PM"):1;
	float resolution = params.paramPresent("RESOLUTION")?(float) params.getValueDouble("RESOLUTION"):0.1;
    int specType = params.paramPresent("SPEC_TYPE_MSMS")?((int) params.getValueInt("SPEC_TYPE_MSMS")?1:0):0;
	float ionOffset = specType?AAJumps::massHion:0;
  */
    SpecSet specs1, specs2, *mzxmlSpecs;
	vector<unsigned int> scanNums1, scanNums2;   vector<short> msLevel;

	if(params.paramPresent("INPUT_MZXML1")) {
		char filename[4096];   mzxmlSpecs = new SpecSet;
		short minMsLevel=2; if(params.paramPresent("MIN_MSLEVEL")) minMsLevel=(short) params.getValueInt("MIN_MSLEVEL");
		if(params.paramPresent("ISOTOPES_MODEL")) minMsLevel=1;  // Needed for guessing parent charges

		if( LoadMzxml(params.getValue("INPUT_MZXML1"), *mzxmlSpecs, &msLevel, minMsLevel) == 0 )
			{ cout<<"No spectra found in "<<params.getValue("INPUT_MZXML1")<<" (min MS level set to "<<minMsLevel<<"). Exiting...\n"; return -1; }
		else cout<<"Found "<<mzxmlSpecs->size()<<" specs (min MS level set to "<<minMsLevel<<") in "<<params.getValue("INPUT_MZXML1")<<"\n";

		if(not params.paramPresent("INPUT_MZXML2")) {  // Just convert mzxml to pklbin
			if(not params.paramPresent("OUTPUT_SPECS1")) { cerr<<"ERROR: Missing output filename (parameter OUTPUT_SPECS1 missing in .params file)!\n"; return -1; }

			if(params.paramPresent("ISOTOPES_MODEL")) {
				IsoEnvelope isoEnvs;
				if(!isoEnvs.LoadModel(params.getValue("ISOTOPES_MODEL"))) return -1;

				// Guess charges for all MS2 spectra based on the MS1 isotopic envelopes
				int idxMS1 = -1;
				for(unsigned int specIdx=0; specIdx<mzxmlSpecs->size(); specIdx++)
					if(msLevel[specIdx]==1) idxMS1=(int)specIdx;
					else if(idxMS1>=0) (*mzxmlSpecs)[specIdx].guessPrecursorZPM((*mzxmlSpecs)[idxMS1],.1,3,isoEnvs,true);

				// Restore minMsLevel to what was indicated in the .params file
				if(params.paramPresent("MIN_MSLEVEL") and params.getValueInt("MIN_MSLEVEL")>1)
					mzxmlSpecs->extract(msLevel,(short)params.getValueInt("MIN_MSLEVEL"),*mzxmlSpecs);
			}

			sprintf(filename,"%s.pklbin",params.getValue("OUTPUT_SPECS1"));         mzxmlSpecs->SaveSpecSet_pklbin(filename);
			scanNums1.resize(mzxmlSpecs->size()); for(unsigned int specIdx=0; specIdx<mzxmlSpecs->size(); specIdx++) scanNums1[specIdx]=(*mzxmlSpecs)[specIdx].scan;
			sprintf(filename,"%s_scannums.bin",params.getValue("OUTPUT_SPECS1"));   Save_binArray(filename, scanNums1);
			sprintf(filename,"%s_mslevels.bin",params.getValue("OUTPUT_SPECS1"));   Save_binArray(filename, msLevel);
			return 0;
		} else {  // Processing the mzxml files, no need to output anything
			mzxmlSpecs->extract(msLevel, (short)1, specs1);
			if( LoadMzxml(params.getValue("INPUT_MZXML2"), *mzxmlSpecs, &msLevel, minMsLevel) == 0 )
				{ cout<<"No spectra found in "<<params.getValue("INPUT_MZXML2")<<" (min MS level set to "<<minMsLevel<<"). Exiting...\n"; return -1; }
			else cout<<"Found "<<mzxmlSpecs->size()<<" specs (min MS level set to "<<minMsLevel<<") in "<<params.getValue("INPUT_MZXML2")<<"\n";
			mzxmlSpecs->extract(msLevel, (short)1, specs2);
		}
		delete mzxmlSpecs;
	} else {  // Attempt to load spectra from INPUT_SPECS1/INPUT_SPECS2
		if(not(params.paramPresent("INPUT_SPECS1"))) { cerr<<"ERROR: missing parameter INPUT_SPECS1 in .params file!\n"; return -1; }
		if(not(params.paramPresent("INPUT_SPECS2"))) { cerr<<"ERROR: missing parameter INPUT_SPECS2 in .params file!\n"; return -1; }

	    if(specs1.LoadSpecSet_pklbin(params.getValue("INPUT_SPECS1"))<0) { cerr<<"ERROR reading "<<params.getValue("INPUT_SPECS1")<<"!\n"; return -1; }
	    cout << "Loading "<<params.getValue("INPUT_SPECS1")<<" complete. Num specs: " << specs1.size() <<"\n";
		scanNums1.resize(specs1.size());  for(unsigned int pivot=0; pivot<scanNums1.size(); pivot++) scanNums1[pivot]=pivot;

	    if(specs2.LoadSpecSet_pklbin(params.getValue("INPUT_SPECS2"))<0) { cerr<<"ERROR reading "<<params.getValue("INPUT_SPECS2")<<"!\n"; return -1; }
	    cout << "Loading "<<params.getValue("INPUT_SPECS2")<<" complete. Num specs: " << specs2.size() <<"\n";
		scanNums2.resize(specs2.size());  for(unsigned int pivot=0; pivot<scanNums2.size(); pivot++) scanNums2[pivot]=pivot;
	}

	vector<TwoValues<unsigned int> > matchedScans;

	if(params.paramPresent("OUTPUT_ALIGNS")) {
		AlignChromatography(specs1, specs2, peakTol, matchedScans, true);

	    FILE *fp = fopen(params.getValue("OUTPUT_ALIGNS"),"wb");
	    if ((int)fp==0) { cerr << "ERROR: cannot open " << params.getValue("OUTPUT_ALIGNS") << "\n";   return -1; }
	    if(params.paramPresent("INPUT_MZXML1")) {
		    fprintf(fp,"%s\n",params.getValue("INPUT_MZXML1"));
		    fprintf(fp,"%s\n",params.getValue("INPUT_MZXML2"));
	    } else {
		    fprintf(fp,"%s\n",params.getValue("INPUT_MZXML1"));
		    fprintf(fp,"%s\n",params.getValue("INPUT_MZXML2"));
	    }
		for(int matchIdx=(int)matchedScans.size(); matchIdx>=0; matchIdx--)
			fprintf(fp,"%d %d\n",specs1[matchedScans[matchIdx][0]].scan,specs2[matchedScans[matchIdx][1]].scan);
	    fclose(fp);
	}

	return 0;
}
