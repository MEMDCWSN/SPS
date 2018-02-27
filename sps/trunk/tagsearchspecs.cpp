#include "spectrum.h"
#include "inputParams.h"
#include "batch.h"
#include "alignment_scoring.h"
#include "db_fasta.h"
#include "tags.h"
#include <ctime>

struct ProteinMatch {
  unsigned int protIdx;
  float score, percScore, tagPos;
  string peptide;

  ProteinMatch() { protIdx=0; score=0; percScore=0; peptide=""; tagPos=0; }
  void set(int in_protIdx, float in_score, float in_percScore, string &in_peptide, float in_tagPos=0)
    { protIdx=in_protIdx; score=in_score; percScore=in_percScore; peptide=in_peptide; tagPos=in_tagPos; }
  ProteinMatch &operator=(const ProteinMatch &other)
    { protIdx=other.protIdx; score=other.score; percScore=other.percScore; peptide=other.peptide; tagPos=other.tagPos; return *this; }
};

int main(int argc, char **argv) {
	InputParams params; 	vector<char *> paramStrs;   paramStrs.resize(2);
	paramStrs[0] = "INPUT_FASTA";
	paramStrs[1] = "OUTPUT_FILES_PREFIX";

	if(argc==1)	params.readParams("tagsearchspecs.params"); else params.readParams(argv[1]);
	if(!params.confirmParams(paramStrs)) {
		cerr << "ERROR: Parameters file ";
		if(argc==1) cerr<<"tagsearchspecs.params"; else cerr<<argv[1];
		cerr << " is incomplete. One of the following is missing: INPUT_FASTA, OUTPUT_FILES_PREFIX\n";
		return -1;
	}

	char *dbFN = params.getValue(paramStrs[0]);
	char *outputPrefix = params.getValue(paramStrs[1]);

	unsigned int tagLen = params.paramPresent("TAG_LEN")?(unsigned int) params.getValueInt("TAG_LEN"):6;
	unsigned int maxNumMissedPeaks = params.paramPresent("MAX_MISSED_PEAKS")?(unsigned int) params.getValueInt("MAX_MISSED_PEAKS"):1;
//	int numJumps = params.paramPresent("DOUBLE_AA_JUMPS")?(int) params.getValueInt("DOUBLE_AA_JUMPS"):1;
//	short minMatchFlanking = params.paramPresent("MATCH_TAG_FLANKING_MASSES")?(short) params.getValueInt("MATCH_TAG_FLANKING_MASSES"):1;
//	float minRatio = params.paramPresent("MIN_RATIO")?(float) params.getValueDouble("MIN_RATIO"):0;
	float peakTol = params.paramPresent("TOLERANCE_PEAK")?(float) params.getValueDouble("TOLERANCE_PEAK"):0.5;
//	float pmTol = params.paramPresent("TOLERANCE_PM")?(float) params.getValueDouble("TOLERANCE_PM"):1;
    int specType = params.paramPresent("SPEC_TYPE_MSMS")?((int) params.getValueInt("SPEC_TYPE_MSMS")?1:0):0;
	float ionOffset = specType?AAJumps::massHion:0;
    bool maxParsimony = (bool)params.paramPresent("MAX_PARSIMONY")?((int) params.getValueInt("MAX_PARSIMONY")?1:0):0;
	float maxCharge = params.paramPresent("MAX_PEAK_CHARGE")?(float) params.getValueDouble("MAX_PEAK_CHARGE"):1.0;
	float missingPeakPenalty = params.paramPresent("PENALTY_MISSING_PEAK")?(float) params.getValueDouble("PENALTY_MISSING_PEAK"):-1.0;
	float noisePeakPenalty = params.paramPresent("PENALTY_NOISE_PEAK")?(float) params.getValueDouble("PENALTY_NOISE_PEAK"):-1.0;

    SpecSet specSet;   short loadOk=0;
    if(params.paramPresent("INPUT_SPECS")) loadOk=specSet.LoadSpecSet_pkl(params.getValue("INPUT_SPECS"));
    else if(params.paramPresent("INPUT_SPECS_PKLBIN")) loadOk=specSet.LoadSpecSet_pklbin(params.getValue("INPUT_SPECS_PKLBIN"));
    if (loadOk<=0 or specSet.size()==0) { cerr<<"Error loading spectra!\n"; return -1; }
	vector<float> specScores(specSet.size());
	for(unsigned int specIdx=0; specIdx<specSet.size(); specIdx++) {
		specSet[specIdx].addZPMpeaks(peakTol, ionOffset, true);
		specScores[specIdx]=0; for(unsigned int peakIdx=0; peakIdx<specSet[specIdx].size();peakIdx++) specScores[specIdx]+=specSet[specIdx][peakIdx][1];
	}

	DB_fasta db; if(db.Load(dbFN)<=0) return -1;
	db.replaceAA('L','I');   db.replaceAA('K','Q');

	FILE *pepsFID = (FILE *)fopen("tmp_matches.txt","wb");
	fprintf(pepsFID,"SpecIdx\tPeptide\tTag pos\tTag score\tPerc peptide score\tProtIdx\tProtein\n");

	vector<Tag> tags;  AAJumps jumps(1);
	list<pair<int,string> > matches;
    vector<list<ProteinMatch> > proteinHits(specSet.size());   ProteinMatch curMatch;   string tmpString;

	vector<list<unsigned int> > matchedTagsIdx;  list<unsigned int>::iterator iterIdx;
	vector<list<float> > matchedTagsScores;      list<float>::iterator iterScores;
	vector<list<unsigned int> > matchedTagsPos;  list<unsigned int>::iterator iterPos;
    for(unsigned int protIdx=0; protIdx<db.size(); protIdx++) {
		cout<<"Processing protein "<<protIdx+1<<" out of "<<db.size()<<"..."; cout.flush();
    	if(ExtractTags(db[protIdx], tags, tagLen)) {
//cerr<<"Tags: "; for(unsigned int i=0;i<tags.size();i++) cerr<<"("<<i<<","<<tags[i].asString(tmpString)<<")"; cerr<<"\n"; cerr.flush();
    		MatchTagsToSpectra(specSet,tags,peakTol,maxCharge,missingPeakPenalty,noisePeakPenalty,maxNumMissedPeaks,matchedTagsIdx,matchedTagsScores,matchedTagsPos);

			for(unsigned int specIdx=0; specIdx<specSet.size(); specIdx++)
				if(matchedTagsIdx[specIdx].size()) {
					iterScores  = matchedTagsScores[specIdx].begin();
					iterPos     = matchedTagsPos[specIdx].begin();
					for(iterIdx=matchedTagsIdx[specIdx].begin(); iterIdx!=matchedTagsIdx[specIdx].end(); iterIdx++,iterScores++,iterPos++) {
						curMatch.set(protIdx,*iterScores,(*iterScores)/specScores[specIdx],tags[*iterIdx].asString(tmpString),specSet[specIdx][*iterPos][0]);
						proteinHits[specIdx].push_back(curMatch);
					}
				}
    	}
    	cout<<"done\n"; cout.flush();
    }

	vector<list<int> > simpleProteinHits(specSet.size());
	if(maxParsimony) {
	 	for(unsigned int specIdx=0; specIdx<specSet.size(); specIdx++) { simpleProteinHits[specIdx].clear();
	    	for(list<ProteinMatch>::iterator iter=proteinHits[specIdx].begin(); iter!=proteinHits[specIdx].end(); iter++)
	      		simpleProteinHits[specIdx].push_back(iter->protIdx);
	    }

		MaximumParsimony(simpleProteinHits);
	}

	for(unsigned int specIdx=0; specIdx<specSet.size(); specIdx++) { if(proteinHits[specIdx].empty()) continue;
		for(list<ProteinMatch>::iterator iter=proteinHits[specIdx].begin(); iter!=proteinHits[specIdx].end(); iter++) {
			if(maxParsimony and simpleProteinHits[specIdx].front()!=iter->protIdx) continue;
			curMatch = *iter;
	    	fprintf(pepsFID,"%d\t%s\t%.1f\t%.1f\t%.1f\t%d",specIdx+1,curMatch.peptide.c_str(),curMatch.tagPos,curMatch.score,curMatch.percScore,curMatch.protIdx);
	    	fprintf(pepsFID,"\t%s %s\n",db.getID(curMatch.protIdx),db.getDesc(curMatch.protIdx));
		}
	}

	fclose(pepsFID);

	return 0;
}
