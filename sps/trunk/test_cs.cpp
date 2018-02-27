#include "spectrum.h"
#include "alignment_scoring.h"

int main(int argc, char *argv) {
	SpecSet contigs, specs;

	contigs.LoadSpecSet_pklbin("/home/nbandeira/data/karl/6mix/orbcal_p1/assembly/sps_seqs.pklbin");
	cout<<"Got "<<contigs.size()<<" contigs\n";
	specs.LoadSpecSet_pklbin("/home/nbandeira/data/karl/6mix/spectra/orbcal_contigs.pklbin");
	cout<<"Got "<<specs.size()<<" specs\n";

	TwoValues<int> res, bestCandidateMP;
	TwoValues<float> bestCandidateScores;
	AAJumps validJumps(0);
	list<float> shiftScores;
	list<TwoValues<unsigned int> > shiftPairs;
	vector<list<TwoValues<int> > > shiftMatchedPeaks;

	res=computeShifts(specs[556], contigs[203], 0.4, 1.0, 0, 0, 0, validJumps,
			shiftScores, shiftPairs, shiftMatchedPeaks, bestCandidateScores, bestCandidateMP,
			0, false);

	cout << "Got "<< shiftScores.size() << " shifts\n";
	for(list<TwoValues<unsigned int> >::iterator iter=shiftPairs.begin(); iter!=shiftPairs.end(); iter++) {
		if(shiftMatchedPeaks[(*iter)[0]].size()>=6) {
			cout<<" ("<<(float((*iter)[0])-float(res[0]))*0.1<<","<<(float((*iter)[1])-float(res[0]))*0.1<<"), ";
			cout<<shiftMatchedPeaks[(*iter)[0]].size()<<" matched peaks: ";
			for(list<TwoValues<int> >::iterator i2=shiftMatchedPeaks[(*iter)[0]].begin(); i2!=shiftMatchedPeaks[(*iter)[0]].end(); i2++)
				cout<<"["<<(*i2)[0]<<","<<(*i2)[1]<<"]";
			cout<<endl;
		}
	}

	int nzCount=0;
	for(list<float>::iterator iter=shiftScores.begin(); iter!=shiftScores.end(); iter++)
		nzCount += (*iter)>0;
	cout <<"Got "<<nzCount<<" shifts with a score higher than zero\n";

	return 0;
}
