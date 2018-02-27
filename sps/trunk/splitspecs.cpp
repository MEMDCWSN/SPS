#include "spectrum.h"
#include "utils.h"
#include <string>

using namespace std;

int main(int argc, char **argv){
	unsigned int pivot, specIdx, fIdx, keepIdx;
	SpecSet specs, partial;
	vector<Spectrum *> chunk;
//	specs.LoadSpecSet_pklbin("pairs.pklbin");
//	specs.LoadSpecSet_pklbin("pairs_pv005_mp6.pklbin");
	specs.LoadSpecSet_pklbin("pairs_mp6d.pklbin");
	cout<<"Done loading "<<specs.size()<<" spectra\n";
	
	// Computes/filters-by number of matched peaks
	vector<vector<unsigned int> > pairsToKeep;
//	Load_binArray("idx_pairs_cast.bin",pairsToKeep);
//	Load_binArray("idx_pairs_annot.bin",pairsToKeep);

	vector<float> modpos; Load_binArray("pairs_modpos_mp6d.bin",modpos);
	
	vector<unsigned int> numPeaks(specs.size());
	if(specs.size()%2!=0) cerr<<"ERROR: expected even number of spectra!\n";
	chunk.resize(specs.size());   pivot=0;
	for(specIdx=0; specIdx<specs.size(); specIdx+=2) {
		numPeaks[specIdx]=specs[specIdx].size();
		numPeaks[specIdx+1]=specs[specIdx+1].size();
		if(min(numPeaks[specIdx],numPeaks[specIdx+1])>=6) {
			modpos[pivot/2] = modpos[specIdx/2];
			chunk[pivot++] = &specs[specIdx];
			chunk[pivot++] = &specs[specIdx+1];
		}
	}
	Save_binArray("numPeaks_mp6d.bin",numPeaks);
	modpos.resize(pivot/2); Save_binArray("pairs_modpos_mp6dF.bin",modpos);
	chunk.resize(pivot); SaveSpecSet_pklbin("pairs_mp6dF.pklbin",chunk);

	return 0;
	
	chunk.resize(2*pairsToKeep.size());   pivot=0;
	for(keepIdx=0; keepIdx<pairsToKeep.size(); keepIdx++) {
		specIdx = 2*(pairsToKeep[keepIdx][0]-1);
		chunk[pivot++] = &specs[specIdx++];
		chunk[pivot++] = &specs[specIdx];
	}
	cout<<"Outputting "<<pivot<<" spectra... "; cout.flush();
//	SaveSpecSet_pklbin("pairs_castW.pklbin",chunk);
	SaveSpecSet_pklbin("pairs_pv005_mp6_annot.pklbin",chunk);
	cout<<"done\n"; cout.flush();

	return 0;
	
	specs.LoadSpecSet_pklbin("../../../clusters/shewMEms.pklbin");
	cout<<"Done loading "<<specs.size()<<" spectra\n";
	Load_binArray("specidx_pairs_cast.bin",pairsToKeep);
	chunk.resize(pairsToKeep.size());   pivot==0;
	for(keepIdx=0; keepIdx<pairsToKeep.size(); keepIdx++) chunk[keepIdx] = &specs[pairsToKeep[keepIdx][0]-1];
	cout<<"Outputting "<<pivot<<" spectra... "; cout.flush();
	SaveSpecSet_pklbin("specs_pairs_castW.pklbin",chunk);
	cout<<"done\n"; cout.flush();

	return 0;

	// Computes/filters-by number of matched peaks
	for(specIdx=0; specIdx<specs.size(); specIdx++) numPeaks[specIdx]=specs[specIdx].size();
	Save_binArray("numPeaks.bin",numPeaks);
	
	chunk.resize(specs.size());   pivot=0;
	for(specIdx=0; specIdx<specs.size(); specIdx++)
		if(specs[specIdx].size()>=6) chunk[pivot++]=&specs[specIdx];
	chunk.resize(pivot);
	cout<<"Outputting "<<pivot<<" spectra... "; cout.flush();
	SaveSpecSet_pklbin("pairs_pv005_mp6.pklbin",chunk);
	cout<<"done\n"; cout.flush();
	return 0;
	
	// Splits pairs into subsets of pairs connecting at least one spectrum with an index in the corresponding range 
	vector<vector<unsigned int> > index;
	Load_binListArray<unsigned int,vector<unsigned int>,vector<unsigned int>::iterator>("pairs_split.bla",index);
	cout<<"Loaded "<<index.size()<<" lists of spectrum indices\n";
	
	char filename[256];
	vector<unsigned int>::iterator iter;
	for(pivot=0; pivot<index.size(); pivot++) {
		cout<<"Chunk "<<pivot<<": "<<index[pivot].size()<<" spectra... "; cout.flush();
		chunk.resize(index[pivot].size());
		for(specIdx=0, iter=index[pivot].begin(); specIdx<index[pivot].size(); specIdx++, iter++)
			chunk[specIdx]=&specs[(*iter)-1];
		sprintf(filename,"pairs_%2d.pklbin",pivot+1); if(filename[6]==' ') filename[6]='0';
		cout<<"saving "<<filename<<"... "; cout.flush();
		SaveSpecSet_pklbin(filename,chunk);
		cout<<"done\n";
	}
	
	return 0;
}
