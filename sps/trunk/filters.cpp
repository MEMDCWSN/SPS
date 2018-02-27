#include "filters.h"

namespace specnets
{

	void instantiate_FilterTriangles() {
		vector<unsigned int> selected;
		vector<Results_SP> aligns_SP;    FilterTriangles(aligns_SP, 0, 0, 0, selected);
		vector<Results_PA> aligns_PA;    FilterTriangles(aligns_PA, 0, 0, 0, selected);
		vector<Results_ASP> aligns_ASP;  FilterTriangles(aligns_ASP, 0, 0, 0, selected);
	}

	void FilterAligns_instantiate() {
		vector<unsigned int> vui;   vector<TwoValues<float> > vtvf;   vector<vector<float> > vvf;
		vector<Results_ASP> asp;   vector<Results_PA> pa;
		FilterAligns<Results_ASP>(asp,vui,vtvf,vtvf,vvf,0,0,0,false);
		FilterAligns<Results_PA>(pa,vui,vtvf,vtvf,vvf,0,0,0,false);
	}

	void FilterMatchRatioASP(SpecSet &specSet, vector<Results_ASP> &results, short ratioType, vector<TwoValues<float> > &ratios) {
			float szOverlap, shift, ratio1, ratio2, totalScore1, totalScore2;
			int spec1, spec2, i, j;

			ratios.resize(results.size());
			for(i=0; i<results.size(); i++) {
					spec1 = results[i].spec1;     spec2 = results[i].spec2;     shift = results[i].shift1;
					szOverlap = min(specSet[spec1].parentMass, specSet[spec2].parentMass);

					totalScore1 = 0;
					for(j=0; j<=specSet[spec1].size() && specSet[spec1][j][0]<=szOverlap-56; j++) if (specSet[spec1][j][0]>56) totalScore1+=specSet[spec1][j][1];
					while (j<=specSet[spec1].size() && specSet[spec1][j][0]<=specSet[spec1].parentMass-szOverlap+56) j++;
					for(; j<=specSet[spec1].size() && specSet[spec1][j][0]<=specSet[spec1].parentMass-56; j++) totalScore1+=specSet[spec1][j][1];

					totalScore2 = 0;
					for(j=0; j<=specSet[spec2].size() && specSet[spec2][j][0]<=szOverlap-56; j++) if (specSet[spec2][j][0]>56) totalScore2+=specSet[spec2][j][1];
					while (j<=specSet[spec2].size() && specSet[spec2][j][0]<=specSet[spec2].parentMass-szOverlap+56) j++;
					for(; j<=specSet[spec2].size() && specSet[spec2][j][0]<=specSet[spec2].parentMass-56; j++) totalScore2+=specSet[spec2][j][1];

	//cerr << "score1 = " << results[i].score1 << ", score2 = " << results[i].score2 << endl;
	//cerr << "totalScore1 = " << totalScore1 << ", totalScore2 = " << totalScore2 << endl;

					switch (ratioType) {
							case 6: if (totalScore1-results[i].score1>1) ratio1 = results[i].score1/(totalScore1-results[i].score1); else ratio1=1;
								 if (totalScore2-results[i].score2>1) ratio2 = results[i].score2/(totalScore2-results[i].score2); else ratio2=1;
								 break;
							case 7: ratio1 = results[i].score1/totalScore1;  ratio2 = results[i].score2/totalScore2; break;
							default: ratio1 = 0; ratio2 = 0;
					}

					if (ratio1<0) ratios[i][0]=0; else ratios[i][0]=ratio1;
					if (ratio2<0) ratios[i][1]=0; else ratios[i][1]=ratio2;
	//cerr << " --- ratio 1 = " << ratios[i][0] << ", ratio 2 = " << ratios[i][1] << endl;
			}
	}

	void ComputeScoreMeans(SpecSet &specSet, vector<Results_ASP> &results, vector<TwoValues<float> > &scoreMeans) {
			int i,specIdx,n,numSpecs = specSet.size();

			scoreMeans.resize(numSpecs);
			for(i=0; i<numSpecs; i++) { scoreMeans[i][0]=0; scoreMeans[i][1]=0; }

			for(i=0; i<results.size(); i++) {
					specIdx = results[i].spec1;   n = (int)scoreMeans[specIdx][1];
					scoreMeans[specIdx][0] = scoreMeans[specIdx][0]*n/(n+1) + results[i].score1/(n+1);
					scoreMeans[specIdx][1] = n+1;

					specIdx = results[i].spec2;   n = (int)scoreMeans[specIdx][1];
					scoreMeans[specIdx][0] = scoreMeans[specIdx][0]*n/(n+1) + results[i].score2/(n+1);
					scoreMeans[specIdx][1] = n+1;
			}
	}
}
