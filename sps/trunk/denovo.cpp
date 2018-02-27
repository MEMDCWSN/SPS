#include "denovo.h"
#include "utils.h"
#include <cmath>
#include <limits.h>

namespace specnets
{

	static const bool DENOVO_DEBUG=false;
	static const unsigned int MAX_NUM_SEQS=100000;      // Stops generation of additional de novo sequences
	static const unsigned int SOFT_MAX_NUM_SEQS=2;  // Triggers removal of repeated de novo sequences

	static unsigned int module_num_rejected_seqs=0;  // Used by denovo() to count number of sequences rejected by denovo_buildSeqs_enforcePM()

	class DenovoSeq {
		public:
		list<TwoValues<float> > peaks;
		float score;

		DenovoSeq() { peaks.clear(); score=0; }
		unsigned int size() { return peaks.size(); }
		bool empty() { return peaks.empty(); }
		bool operator<(DenovoSeq &other) {
			if(peaks.size()<other.size()) return true;    if(peaks.size()>other.size()) return false;
			list<TwoValues<float> >::iterator iter1 = peaks.begin(), iter2 = other.peaks.begin();
			while(iter1!=peaks.end()) {
				double diff = ((double)(*iter1)[0])-((double)(*iter2)[0]);
				if(diff<-0.0001) return true;
				if(diff>0.0001) return false;
				iter1++; iter2++;
			}
			return true;
		}

		DenovoSeq &operator=(DenovoSeq &other) { peaks=other.peaks; score=other.score; return *this; }

		bool same_peaks(DenovoSeq &other, float massTol=0.001) {
			if(peaks.size()!=other.size()) return false;
			list<TwoValues<float> >::iterator iter1 = peaks.begin(), iter2 = other.peaks.begin();
			while(iter1!=peaks.end()) {
				if(fabs(((double)(*iter1)[0])-((double)(*iter2)[0]))>massTol) return false;
				iter1++; iter2++;
			}
			return true;
		}
	};


	bool cmp_peak_lists_ui(list<unsigned int> &L1, list<unsigned int> &L2) {
		if(L1.size()<L2.size()) return true;
		if(L1.size()>L2.size()) return false;
		list<unsigned int>::iterator iter1 = L1.begin(), iter2 = L2.begin();
		while(iter1!=L1.end()) {
			if((*iter1)<(*iter2)) return true;
			if((*iter1)>(*iter2)) return false;
			iter1++; iter2++;
		}
		return true;
	}

	bool same_peak_lists_ui(list<unsigned int> &L1, list<unsigned int> &L2) {
		if(L1.size()!=L2.size()) return false;
		list<unsigned int>::iterator iter1 = L1.begin(), iter2 = L2.begin();
		while(iter1!=L1.end()) {
			if((*iter1)!=(*iter2)) return false;
			iter1++; iter2++;
		}
		return true;
	}

	/*
	bool cmp_peak_lists_tvf(list<TwoValues<float> > &L1, list<TwoValues<float> > &L2) {
		if(L1.size()<L2.size()) return true;
		if(L1.size()>L2.size()) return false;
		list<TwoValues<float> >::iterator iter1 = L1.begin(), iter2 = L2.begin();
		while(iter1!=L1.end()) {
			double diff = ((double)(*iter1)[0])-((double)(*iter2)[0]);
			if(diff<-0.0001) return true;
			if(diff>0.0001) return false;
			iter1++; iter2++;
		}
		return true;
	}

	bool same_peak_lists_tvf(list<TwoValues<float> > &L1, list<TwoValues<float> > &L2, float massTol=0.001) {
		if(L1.size()!=L2.size()) return false;
		list<TwoValues<float> >::iterator iter1 = L1.begin(), iter2 = L2.begin();
		while(iter1!=L1.end()) {
			if(fabs(((double)(*iter1)[0])-((double)(*iter2)[0]))>massTol) return false;
			iter1++; iter2++;
		}
		return true;
	}
	*/

	unsigned int unique_peak_lists(list<DenovoSeq> &seqs, float massTol=0.001) {
		list<DenovoSeq>::iterator seqsIter, seqNext;

		seqs.sort();   seqsIter = seqs.begin();
		while(seqsIter!=seqs.end()) {
			seqNext=seqsIter; seqNext++;
			if(seqNext!=seqs.end() and seqsIter->same_peaks(*seqNext,massTol))
				seqsIter=seqs.erase(seqsIter); else seqsIter++;
		}

		return seqs.size();
	}

	void denovo_buildSeqs_idx(Spectrum &spec, vector<vector<list<TwoValues<unsigned int> > > > &predPairs,
												vector<vector<float> > &maxScores,
												TwoValues<unsigned int> &endsPair, int newPeak,
												unsigned int firstSuff, float minPathScore, float peakTol, float pmOffset,
												list<list<unsigned int> > &seqsPeaks, list<float> &seqsScores) {
		list<TwoValues<unsigned int> >::iterator predIter;
		list<list<unsigned int> >::iterator seqsIter;
		list<float>::iterator scoresIter;
		float maxSeqScore=0, endsPairScore=0;

	/*cerr<<"Entering denovo_buildSeqs with endsPair ("<<endsPair[0]<<","<<endsPair[1]<<"): ";
	if(seqsPeaks.size()==0) cerr<<"(empty)"; else
		for(list<unsigned int>::iterator iter=seqsPeaks.front().begin(); iter!=seqsPeaks.front().end(); iter++)
			cerr<<*iter<<" ";
	cerr<<", seqScore = "<<seqsScores.front()<<"\n";
	*/
		if(fabs(spec.parentMass+(pmOffset-AAJumps::massHion)*spec.idDist-spec[endsPair[0]][0]-spec[endsPair[1]][0])<=2*peakTol)
			endsPairScore = max(spec[endsPair[0]][1],spec[endsPair[1]][1]); else endsPairScore = spec[endsPair[0]][1]+spec[endsPair[1]][1];

		list<unsigned int> curSeq; // curSeq is the intermediate list of peaks where new peaks are added to
		if(newPeak<0) {
			curSeq.push_back(endsPair[0]);   curSeq.push_back(endsPair[1]);
			seqsPeaks.push_back(curSeq);     seqsScores.push_back(endsPairScore);
		}

		list<list<unsigned int> > expandedSeqs;   list<float> expandedScores;  // Collection of all expanded de-novo seqs
		list<list<unsigned int> > curNewSeqs;     list<float> curNewScores;    // Intermediate recursion current input seq -> output expanded seqs
		float curScore = seqsScores.front();
		int oldNewPeak=-1;   bool suffixExtension=false;
		for(predIter=predPairs[endsPair[0]][endsPair[1]-firstSuff].begin(); predIter!=predPairs[endsPair[0]][endsPair[1]-firstSuff].end(); predIter++) {
			curSeq = seqsPeaks.front();  unsigned int otherPeak;   // on input, seqsPeaks contains at most one list of peaks
			if(endsPair[0]==(*predIter)[0]) {
				curSeq.push_back((*predIter)[1]); newPeak=(int)(*predIter)[1]; otherPeak=endsPair[0];
			} else {
				curSeq.push_front((*predIter)[0]); newPeak=(int)(*predIter)[0]; otherPeak=endsPair[1];
			}
			float scoreAdd=0,       // Score that needs to be added to curScore to get the score of the path from (*predIter)[0] to (*predIter)[1]
						predsPairScore=0; // Combined score of the predecessor peaks
			if(fabs(spec.parentMass+(pmOffset-AAJumps::massHion)*spec.idDist-spec[newPeak][0]-spec[otherPeak][0])<=2*peakTol) {
				predsPairScore=max(spec[newPeak][1],spec[otherPeak][1]); if(spec[newPeak][1]>spec[otherPeak][1]) scoreAdd = spec[newPeak][1]-spec[otherPeak][1]; else scoreAdd=0;
			} else { predsPairScore=spec[newPeak][1]+spec[otherPeak][1]; scoreAdd = spec[newPeak][1]; }
			maxSeqScore = curScore + scoreAdd + maxScores[(*predIter)[0]][(*predIter)[1]-firstSuff] - predsPairScore;

	/*if(endsPair[0]==6 and endsPair[1]==30) {
		cerr<<"--- endsPair (6,24) score predecessor ("<<(*predIter)[0]<<","<<(*predIter)[1]<<") at maxSeqScore = "<<maxSeqScore<<"\n";
		cerr<<"--- curScore = "<<curScore<<"\n";
		cerr<<"--- scoreAdd = "<<scoreAdd<<"\n";
		cerr<<"--- maxScores[(*predIter)[0]][(*predIter)[1]-firstSuff] = "<<maxScores[(*predIter)[0]][(*predIter)[1]-firstSuff]<<"\n";
		cerr<<"--- predsPairScore = "<<predsPairScore<<"\n";
	}*/
			if(maxSeqScore<minPathScore) continue;

			curNewSeqs.clear();   curNewSeqs.push_back(curSeq);   curNewScores.clear();   curNewScores.push_back(curScore+scoreAdd);

			denovo_buildSeqs_idx(spec,predPairs,maxScores,*predIter,newPeak,firstSuff,minPathScore,peakTol,pmOffset,curNewSeqs,curNewScores);
			expandedSeqs.splice(expandedSeqs.end(),curNewSeqs);
			expandedScores.splice(expandedScores.end(),curNewScores);
		}

		// Quick-not-so-elegant fix to remove repeated sequences generated from symmetric peaks (A,B) via
		//   predecessors (I,B), (A,J), both reachable from (I,J)

		expandedSeqs.sort(cmp_peak_lists_ui);
		list<list<unsigned int> >::iterator repIter = expandedSeqs.begin(), repNext;
		while(repIter!=expandedSeqs.end()) {
			repNext=repIter; repNext++;
			if(repNext!=expandedSeqs.end() and same_peak_lists_ui(*repIter,*repNext))
				repIter=expandedSeqs.erase(repIter); else repIter++;
		}

		if(predPairs[endsPair[0]][endsPair[1]-firstSuff].size()>0) {
			seqsPeaks.clear();   seqsPeaks.splice(seqsPeaks.end(),expandedSeqs);
			seqsScores.clear();  seqsScores.splice(seqsScores.end(),expandedScores);
		}

	//cerr<<"Finished with endsPair ("<<endsPair[0]<<","<<endsPair[1]<<") with "<<seqsPeaks.size()<<" sequences\n";
	}

	void denovo_buildSeqs(Spectrum &spec, vector<vector<list<TwoValues<unsigned int> > > > &predPairs,
												vector<vector<float> > &maxScores,
												TwoValues<unsigned int> &endsPair, int newPeak,
												unsigned int firstSuff, float minPathScore, float peakTol, float pmOffset,
												list<list<TwoValues<float> > > &seqsPeaks, list<float> &seqsScores) {
		list<TwoValues<unsigned int> >::iterator predIter;
		list<list<unsigned int> >::iterator seqsIter;
		list<float>::iterator scoresIter;
		float maxSeqScore=0, endsPairScore=0;

		if(fabs(spec.parentMass+(pmOffset-AAJumps::massHion)*spec.idDist-spec[endsPair[0]][0]-spec[endsPair[1]][0])<=2*peakTol)
			endsPairScore = max(spec[endsPair[0]][1],spec[endsPair[1]][1]); else endsPairScore = spec[endsPair[0]][1]+spec[endsPair[1]][1];

		list<TwoValues<float> > curSeq; // curSeq is the intermediate list of peaks where new peaks are added to
		if(newPeak<0) {
			curSeq.push_back(spec[endsPair[0]]);   curSeq.push_back(spec[endsPair[1]]);
			seqsPeaks.clear();               seqsScores.clear();
			seqsPeaks.push_back(curSeq);     seqsScores.push_back(endsPairScore);
		}

		list<list<TwoValues<float> > > expandedSeqs;   list<float> expandedScores;  // Collection of all expanded de-novo seqs
		list<list<TwoValues<float> > > curNewSeqs;     list<float> curNewScores;    // Intermediate recursion current input seq -> output expanded seqs
		float curScore = seqsScores.front();
		for(predIter=predPairs[endsPair[0]][endsPair[1]-firstSuff].begin(); predIter!=predPairs[endsPair[0]][endsPair[1]-firstSuff].end(); predIter++) {
			curSeq = seqsPeaks.front();  unsigned int otherPeak;   // on input, seqsPeaks contains at most one list of peaks
			if(endsPair[0]==(*predIter)[0]) {
				curSeq.push_back(spec[(*predIter)[1]]); newPeak=(int)(*predIter)[1]; otherPeak=endsPair[0];
			} else {
				curSeq.push_front(spec[(*predIter)[0]]); newPeak=(int)(*predIter)[0]; otherPeak=endsPair[1];
			}
			float scoreAdd=0,       // Score that needs to be added to curScore to get the score of the path from (*predIter)[0] to (*predIter)[1]
						predsPairScore=0; // Combined score of the predecessor peaks
			if(fabs(spec.parentMass+(pmOffset-AAJumps::massHion)*spec.idDist-spec[newPeak][0]-spec[otherPeak][0])<=2*peakTol) {
				predsPairScore=max(spec[newPeak][1],spec[otherPeak][1]); if(spec[newPeak][1]>spec[otherPeak][1]) scoreAdd = spec[newPeak][1]-spec[otherPeak][1]; else scoreAdd=0;
			} else { predsPairScore=spec[newPeak][1]+spec[otherPeak][1]; scoreAdd = spec[newPeak][1]; }
			maxSeqScore = curScore + scoreAdd + maxScores[(*predIter)[0]][(*predIter)[1]-firstSuff] - predsPairScore;

			if(maxSeqScore<minPathScore) continue;
	//for(list<TwoValues<float> >::iterator pivotIter=seqsPeaks.front().begin(); pivotIter!=seqsPeaks.front().end(); pivotIter++) cerr<<"..";
	//cerr<<"("<<endsPair[0]<<","<<endsPair[1]<<")<-("<<(*predIter)[0]<<","<<(*predIter)[1]<<"), masses ("<<spec[endsPair[0]][0]<<","<<spec[endsPair[1]][0]<<")<-("<<spec[(*predIter)[0]][0]<<","<<spec[(*predIter)[1]][0]<<")\n"; cerr.flush();

			curNewSeqs.clear();   curNewSeqs.push_back(curSeq);   curNewScores.clear();   curNewScores.push_back(curScore+scoreAdd);

			denovo_buildSeqs(spec,predPairs,maxScores,*predIter,newPeak,firstSuff,minPathScore,peakTol,pmOffset,curNewSeqs,curNewScores);
			expandedSeqs.splice(expandedSeqs.end(),curNewSeqs);
			expandedScores.splice(expandedScores.end(),curNewScores);
		}

		// Quick-not-so-elegant fix to remove repeated sequences generated from symmetric peaks (A,B) via
		//   predecessors (I,B), (A,J), both reachable from (I,J)
	/*
		expandedSeqs.sort(cmp_peak_lists);
		list<list<unsigned int> >::iterator repIter = expandedSeqs.begin(), repNext;
		while(repIter!=expandedSeqs.end()) {
			repNext=repIter; repNext++;
			if(repNext!=expandedSeqs.end() and same_peak_lists(*repIter,*repNext))
				repIter=expandedSeqs.erase(repIter); else repIter++;
		}
	*/

		if(predPairs[endsPair[0]][endsPair[1]-firstSuff].size()>0) {
			seqsPeaks.clear();   seqsPeaks.splice(seqsPeaks.end(),expandedSeqs);
			seqsScores.clear();  seqsScores.splice(seqsScores.end(),expandedScores);
		}

	//cerr<<"Finished with endsPair ("<<endsPair[0]<<","<<endsPair[1]<<") with "<<seqsPeaks.size()<<" sequences\n";
	}

	// denovo_buildSeqs_enforcePM - reconstructs all de novo sequences with a score
	//   higher than minPathScore or only the highest scoring sequence (if topSequenceOnly
	//   is set to true). Also, only returns de novo sequences where the reconstructed
	//   peptide's parent mass is within pmTol of the experimental parentMass
	//
	// spec - spectrum used to build the spectrum graph
	// predPairs  - all predecessors for a pair of [prefix,suffix] peaks
	// maxScores  - maximum score possible for each [prefix,suffix] pair (over all possible predecessors in predPairs)
	// endsPair   - [prefix,suffix] peak indices at the edge of the current sequence in seqsPeaks
	// newPeak    - index of the peak just changed in endsPair (set to -1 if just initialized)
	// firstSuff  - index of the first peak in the spectrum whose mass is larger than half the precursor
	// minPathScore - minimum score for any reported de novo sequence
	// peakTol    - peak mass tolerance (in Daltons)
	// pmOffset   - symmetric peak's masses should add up to spec.parentMass-1+pmOffset
	// pmTol      - parent mass tolerance (in Daltons)
	// cumSequenceOffset - cumulative mass offset between initial peak masses and the theoretical jump mass (for the jump crossing (M+H)/2)
	// jumps      - valid amino acid mass jumps between peaks (exact masses used to compute parent mass offset)
	// topSequenceOnly - if set to true, return only the maximal-scoring de novo sequence
	// seqsPeaks  - peak masses in the current de novo sequences (set to empty if just initialized)
	// seqsScores - scores for the current de novo sequences (set to empty if just initialized)
	float denovo_buildSeqs_enforcePM(Spectrum &spec, vector<vector<list<TwoValues<unsigned int> > > > &predPairs,
												vector<vector<float> > &maxScores,
												TwoValues<unsigned int> &endsPair, int newPeak,
												unsigned int firstSuff, float minPathScore, float peakTol, float pmOffset,
												float pmTol, float cumSequenceOffset, AAJumps &jumps, bool topSequenceOnly,
												list<list<TwoValues<float> > > &seqsPeaks, list<float> &seqsScores) {
		list<TwoValues<unsigned int> >::iterator predIter;
		list<float>::iterator scoresIter;
		float maxSeqScore=0, endsPairScore=0;

	//cerr<<"  >>> starting with pair ("<<endsPair[0]<<","<<endsPair[1]<<"), # predecessors = "<<predPairs[endsPair[0]][endsPair[1]-firstSuff].size()<<", cumSequenceOffset = "<<cumSequenceOffset<<"\n"; cerr.flush();

		if(predPairs[endsPair[0]][endsPair[1]-firstSuff].size()==0) {  // If there are no more predecessors then the recursion ends here
	//cerr<<"  >>> No predecessors for peaks ("<<endsPair[0]<<","<<endsPair[1]<<"), cumSequenceOffset = "<<cumSequenceOffset<<", score/minPathScore = "<<seqsScores.front()<<"/"<<minPathScore<<"\n"; cerr.flush();
			TwoValues<float> first, last;
			first = seqsPeaks.front().front();    last = seqsPeaks.front().back();
	//		if(fabs(cumSequenceOffset)>pmTol+0.00001 or seqsScores.front()<=minPathScore+0.00001) {
	//*		float pmMassDifference = fabs(last[0]-first[0]-cumSequenceOffset-(spec.parentMass-AAJumps::massMH));
	//*		if(pmMassDifference>pmTol+0.00001 or fabs(cumSequenceOffset)>=peakTol+0.00001 or seqsScores.front()<=minPathScore+0.00001) {
			float pmMassDifference = fabs(last[0]-first[0]-(spec.parentMass-AAJumps::massMH));
			if(pmMassDifference>pmTol+0.00001 or seqsScores.front()<=minPathScore-0.00001) {
	if(pmMassDifference>pmTol+0.00001) module_num_rejected_seqs++;
	//	cerr<<"  >>> Rejected sequence with pmMassDifference = "<<pmMassDifference<<", score = "<<seqsScores.front()<<"\n"; cerr.flush();
	//for(list<TwoValues<float> >::iterator iter=seqsPeaks.front().begin(); iter!=seqsPeaks.front().end(); iter++)
	//	cerr<<"  --- ("<<(*iter)[0]<<","<<(*iter)[1]<<")"; cerr<<"\n"; cerr.flush();
				seqsPeaks.clear(); seqsScores.clear(); return minPathScore;
			}
			float massOffset = seqsPeaks.front().front()[0];  // Correct eventual mass offsets caused by basing de novo sequence on spec[endsPair[0]][0] (see curSeq initialization below)
			for(list<TwoValues<float> >::iterator massIter=seqsPeaks.front().begin(); massIter!=seqsPeaks.front().end(); massIter++)
				(*massIter)[0]-=massOffset;
	//cerr<<"  >>-> Got final sequence with pmMassDifference = "<<pmMassDifference<<", cumSequenceOffset = "<<cumSequenceOffset<<", score = "<<seqsScores.front()<<"\n"; cerr.flush();
	//for(list<TwoValues<float> >::iterator iter=seqsPeaks.front().begin(); iter!=seqsPeaks.front().end(); iter++)
	//	cerr<<"("<<(*iter)[0]<<","<<(*iter)[1]<<")"; cerr<<"\n"; cerr.flush();
	//		return max(topSequenceOnly?seqsScores.front():0,minPathScore);
			return seqsScores.front();
		}

		if(fabs(spec.parentMass+pmOffset-AAJumps::massHion*spec.idDist-spec[endsPair[0]][0]-spec[endsPair[1]][0])<=2*peakTol)
			endsPairScore = max(spec[endsPair[0]][1],spec[endsPair[1]][1]); else endsPairScore = spec[endsPair[0]][1]+spec[endsPair[1]][1];

		list<TwoValues<float> > curSeq; // curSeq is the intermediate list of peaks where new peaks are added to
		if(newPeak<0) {
			curSeq.push_back(spec[endsPair[0]]);   curSeq.push_back(TwoValues<float>(spec[endsPair[1]][0]+cumSequenceOffset,spec[endsPair[1]][1]));
	//cerr<<" -- denovo_buildSeqs_enforcePM: starting with ("<<spec[endsPair[0]][0]<<","<<spec[endsPair[0]][1]<<"), ("<<spec[endsPair[1]][0]+cumSequenceOffset<<","<<spec[endsPair[1]][1]<<"), endsPairScore="<<endsPairScore<<"\n"; cerr.flush();
			seqsPeaks.clear();               seqsScores.clear();
			seqsPeaks.push_back(curSeq);     seqsScores.push_back(endsPairScore);
		}

		list<list<TwoValues<float> > > expandedSeqs;   list<float> expandedScores;  // Collection of all expanded de-novo seqs
		list<list<TwoValues<float> > > curNewSeqs;     list<float> curNewScores;    // Intermediate recursion current input seq -> output expanded seqs
		float curJump;                            // Current jump between peaks
		TwoValues<unsigned int> jumpsBounds;      // Possible exact masses for curJump
		float curScore = seqsScores.front();
		for(predIter=predPairs[endsPair[0]][endsPair[1]-firstSuff].begin(); predIter!=predPairs[endsPair[0]][endsPair[1]-firstSuff].end(); predIter++) {
	//		curSeq = seqsPeaks.front();   unsigned int otherPeak;   // on input, seqsPeaks contains at most one list of peaks
			unsigned int otherPeak;   // on input, seqsPeaks contains at most one list of peaks
			if(endsPair[0]==(*predIter)[0]) {
	//			curSeq.push_back((*predIter)[1]);
				newPeak=(int)(*predIter)[1]; otherPeak=endsPair[0]; curJump=spec[newPeak][0]-spec[endsPair[1]][0];
			} else {
	//			curSeq.push_front((*predIter)[0]);
				newPeak=(int)(*predIter)[0]; otherPeak=endsPair[1]; curJump=spec[endsPair[0]][0]-spec[newPeak][0];
			}
			float scoreAdd=0, predsPairScore=0;
			if(fabs(spec.parentMass-AAJumps::massHion*spec.idDist-spec[newPeak][0]-spec[otherPeak][0])<=2*peakTol) {
				predsPairScore=max(spec[newPeak][1],spec[otherPeak][1]); if(spec[newPeak][1]>spec[otherPeak][1]) scoreAdd = spec[newPeak][1]-spec[otherPeak][1]; else scoreAdd=0;
			} else { predsPairScore=spec[newPeak][1]+spec[otherPeak][1]; scoreAdd = spec[newPeak][1]; }
			maxSeqScore = curScore + scoreAdd + maxScores[(*predIter)[0]][(*predIter)[1]-firstSuff] - predsPairScore;

	//cerr<<"  ...> expanding predecessor with mass "<<curJump<<", scores = "<<maxSeqScore<<"/"<<minPathScore<<"\n";
			if(maxSeqScore<minPathScore-0.00001) continue;

			jumps.find(curJump,peakTol,jumpsBounds);
	//cerr<<"  ...> jumps mass "<<curJump<<" hits to ["<<jumpsBounds[0]<<":"<<jumpsBounds[1]<<"]\n";
			for(unsigned int jumpIdx=jumpsBounds[0]; jumpIdx<jumpsBounds[1]; jumpIdx++) {
				curSeq = seqsPeaks.front();
				if(newPeak==(int)(*predIter)[1]) curSeq.push_back(TwoValues<float>(curSeq.back()[0]+jumps[jumpIdx],spec[newPeak][1]));
				else curSeq.push_front(TwoValues<float>(curSeq.front()[0]-jumps[jumpIdx],spec[newPeak][1]));
	//cerr<<"  ...> expanding predecessor of ("<<endsPair[0]<<","<<endsPair[1]<<") with jump mass "<<jumps[jumpIdx]<<", incremental score = "<<curScore+scoreAdd<<"\n";
				curNewSeqs.clear();   curNewSeqs.push_back(curSeq);   curNewScores.clear();   curNewScores.push_back(curScore+scoreAdd);

	//*		    minPathScore=denovo_buildSeqs_enforcePM(spec,predPairs,maxScores,*predIter,newPeak,firstSuff,minPathScore,peakTol,pmOffset,pmTol,cumSequenceOffset+(jumps[jumpIdx]-curJump),jumps,topSequenceOnly,curNewSeqs,curNewScores);
					float tmpScore=denovo_buildSeqs_enforcePM(spec,predPairs,maxScores,*predIter,newPeak,firstSuff,minPathScore,peakTol,pmOffset,pmTol,0,jumps,topSequenceOnly,curNewSeqs,curNewScores);

				if(curNewScores.size()>0) {  // If new sequences were reported
	// Do this somewhere else:				if(topSequenceOnly) { expandedSeqs.clear(); expandedScores.clear(); }
					minPathScore = tmpScore;
					expandedSeqs.splice(expandedSeqs.end(),curNewSeqs);
					expandedScores.splice(expandedScores.end(),curNewScores);
	//cerr<<"  --- Intermediate count = "<<expandedSeqs.size()<<" seqs\n";
				}
			}
		}

		seqsPeaks.clear();   seqsScores.clear();
		seqsPeaks.splice(seqsPeaks.end(),expandedSeqs);   seqsScores.splice(seqsScores.end(),expandedScores);
		return minPathScore;
	}

	// denovo_buildSeqs_enforcePM - reconstructs all de novo sequences with a score
	//   higher than minPathScore or only the highest scoring sequence (if topSequenceOnly
	//   is set to true). Also, only returns de novo sequences where the reconstructed
	//   peptide's parent mass is within pmTol of the experimental parentMass
	//
	// spec - spectrum used to build the spectrum graph
	// predPairs  - all predecessors for a pair of [prefix,suffix] peaks
	// maxScores  - maximum score possible for each [prefix,suffix,pmError] tuple (over all possible predecessors in predPairs)
	// endsPair   - [prefix,suffix] peak indices at the edge of the current sequence in seqsPeaks
	// newPeak    - index of the peak just changed in endsPair (set to -1 if just initialized)
	// firstSuff  - index of the first peak in the spectrum whose mass is larger than half the precursor
	// minPathScore - minimum score for any reported de novo sequence
	// minPercScore - percentage of minPathScore that must be met for a path to be reported
	// peakTol    - peak mass tolerance (in Daltons)
	// pmOffset   - symmetric peak's masses should add up to spec.parentMass-1+pmOffset
	// pmTol      - parent mass tolerance (in Daltons)
	// cumSequenceOffset - cumulative mass offset between initial peak masses and the theoretical jump mass (for the jump crossing (M+H)/2)
	// jumps      - valid amino acid mass jumps between peaks (exact masses used to compute parent mass offset)
	// topSequenceOnly - if set to true, return only the maximal-scoring de novo sequence
	// seqsPeaks  - peak masses/indices (when seqsExactMasses = true/false) in the current de novo sequences (set to empty if just initialized)
	// seqsScores - scores for the current de novo sequences (set to empty if just initialized)
	// seqsExactMasses - Determines whether to return de novo sequences as sets of peptide peaks (false) or
	//                    sequences of exact amino acid masses (true). Note that the latter may potentially generate very large numbers of solutions per spectrum!
	float denovo_buildSeqs_enforcePM(Spectrum &spec, vector<vector<list<TwoValues<unsigned int> > > > &predPairs,
												vector<vector<vector<float> > > &maxScores,
												TwoValues<unsigned int> &endsPair, int newPeak,
												unsigned int firstSuff, int tolBinsOffset, float minPathScore, float minPercScore, float peakTol, float pmOffset,
												float pmTol, float resolution, float cumSequenceOffset, AAJumps &jumps, bool topSequenceOnly,
												list<DenovoSeq> &seqsPeaks, bool seqsExactMasses) {
		list<TwoValues<unsigned int> >::iterator predIter;
		TwoValues<float> tmpPair;
		list<float>::iterator scoresIter;
		float maxSeqScore=0, endsPairScore=0;
		int numAA; if(seqsPeaks.empty() or seqsPeaks.front().empty()) numAA=1; else numAA=seqsPeaks.front().size()-1;
	//	int pmTolRange=(int)round((pmTol+(float)numAA*resolution/2)/resolution), numTolBins=2*tolBinsOffset+1;
		int pmTolRange=(int)round(pmTol/resolution), numTolBins=2*tolBinsOffset+1;

	if(DENOVO_DEBUG) cerr<<"  >>> starting with pair ("<<endsPair[0]<<","<<endsPair[1]<<"), # predecessors = "<<predPairs[endsPair[0]][endsPair[1]-firstSuff].size()<<", cumSequenceOffset = "<<cumSequenceOffset<<"\n"; cerr.flush();

		float pmMassDifference, pmMassThresh;
		if(predPairs[endsPair[0]][endsPair[1]-firstSuff].size()==0) {  // If there are no more predecessors then the recursion ends here
			TwoValues<float> first, last;   float seqScore;
			if(numAA==0) { first.set(0,0); last.set(0,0); pmMassDifference = spec.parentMass-AAJumps::massMH; seqScore=0; } else {
				if(seqsExactMasses) {
					first = seqsPeaks.front().peaks.front(); last = seqsPeaks.front().peaks.back();
					pmMassDifference = fabs(last[0]-first[0]-(spec.parentMass-AAJumps::massMH));
				} else {
					first = spec[(unsigned int)(seqsPeaks.front().peaks.front()[0])];
					last  = spec[(unsigned int)(seqsPeaks.front().peaks.back()[0])];
					pmMassDifference = fabs(last[0]-first[0]+cumSequenceOffset-(spec.parentMass-AAJumps::massMH));
				}
				seqScore = seqsPeaks.front().score;
			}
	if(DENOVO_DEBUG) cerr<<"  >>> No predecessors for peaks ("<<endsPair[0]<<","<<endsPair[1]<<")<->("<<first[0]<<","<<last[0]<<"), cumSequenceOffset = "<<cumSequenceOffset<<", score/minPathScore = "<<seqScore<<"/"<<minPathScore<<", --"<<spec.parentMass<<"--"<<AAJumps::massMH<<"--\n"; cerr.flush();
	//		if(fabs(cumSequenceOffset)>pmTol+0.00001 or seqsScores.front()<=minPathScore+0.00001) {
	//*		float pmMassDifference = fabs(last[0]-first[0]-cumSequenceOffset-(spec.parentMass-AAJumps::massMH));
	//*		if(pmMassDifference>pmTol+0.00001 or fabs(cumSequenceOffset)>=peakTol+0.00001 or seqsScores.front()<=minPathScore+0.00001) {
			pmMassThresh = pmTol+max((float)0,(float)numAA)*resolution/2+0.00001;
	//		if(pmMassDifference>pmTol+0.00001 or seqsScores.empty() or seqsScores.front()<=minPathScore-0.00001) {
			if(pmMassDifference>pmMassThresh or seqScore<minPathScore-0.001) {
	if(pmMassDifference>pmTol+0.00001) module_num_rejected_seqs++;
	if(DENOVO_DEBUG) 	cerr<<"  >>> Rejected sequence with pmMassDifference = "<<pmMassDifference<<"/"<<pmMassThresh<<", cumSequenceOffset = "<<cumSequenceOffset<<", score = "<<seqScore<<"/"<<minPathScore-0.001<<": ("<<(pmMassDifference>pmMassThresh)<<","<<(seqScore<minPathScore-0.001)<<")\n"; cerr.flush();
	//for(list<TwoValues<float> >::iterator iter=seqsPeaks.front().begin(); iter!=seqsPeaks.front().end(); iter++)
	//	cerr<<"  --- ("<<(*iter)[0]<<","<<(*iter)[1]<<")"; cerr<<"\n"; cerr.flush();
				seqsPeaks.clear(); return minPathScore;
			}
			if(seqsExactMasses) {
				float massOffset = seqsPeaks.front().peaks.front()[0];  // Correct eventual mass offsets caused by basing de novo sequence on spec[endsPair[0]][0] (see curSeq initialization below)
				for(list<TwoValues<float> >::iterator massIter=seqsPeaks.front().peaks.begin(); massIter!=seqsPeaks.front().peaks.end(); massIter++)
					(*massIter)[0]-=massOffset;
			}
	if(DENOVO_DEBUG) cerr<<"  >>-> Got final sequence with pmMassDifference = "<<pmMassDifference<<", cumSequenceOffset = "<<cumSequenceOffset<<", score = "<<seqScore<<"\n"; cerr.flush();
	//for(list<TwoValues<float> >::iterator iter=seqsPeaks.front().begin(); iter!=seqsPeaks.front().end(); iter++)
	//	cerr<<"("<<(*iter)[0]<<","<<(*iter)[1]<<")"; cerr<<"\n"; cerr.flush();
			return max(minPercScore*seqScore,minPathScore);
		}

		if(fabs(spec.parentMass+pmOffset-AAJumps::massHion*spec.idDist-spec[endsPair[0]][0]-spec[endsPair[1]][0])<=2*peakTol)
			endsPairScore = max(spec[endsPair[0]][1],spec[endsPair[1]][1]); else endsPairScore = spec[endsPair[0]][1]+spec[endsPair[1]][1];

		DenovoSeq curSeq; // curSeq is the intermediate list of peaks where new peaks are added to
		if(newPeak<0) {
			if(seqsExactMasses) { curSeq.peaks.push_back(spec[endsPair[0]]);   curSeq.peaks.push_back(TwoValues<float>(spec[endsPair[1]][0]+cumSequenceOffset,spec[endsPair[1]][1])); }
			else {
				tmpPair.set(float(endsPair[0]),spec[endsPair[0]][1]);   curSeq.peaks.push_back(tmpPair);
				tmpPair.set(float(endsPair[1]),spec[endsPair[1]][1]);   curSeq.peaks.push_back(tmpPair);
			}
	if(DENOVO_DEBUG) cerr<<" -- denovo_buildSeqs_enforcePM: starting with ("<<spec[endsPair[0]][0]<<","<<spec[endsPair[0]][1]<<"), ("<<spec[endsPair[1]][0]+cumSequenceOffset<<","<<spec[endsPair[1]][1]<<"), endsPairScore="<<endsPairScore<<"\n"; cerr.flush();
			seqsPeaks.clear();   curSeq.score = endsPairScore;    seqsPeaks.push_back(curSeq);
		}

		list<DenovoSeq> expandedSeqs;  // Collection of all expanded de-novo seqs
		list<DenovoSeq> curNewSeqs;    // Intermediate recursion current input seq -> output expanded seqs
		float curJump;                              // Current jump between peaks
		TwoValues<unsigned int> jumpsBounds;      // Possible exact masses for curJump
		float curScore = seqsPeaks.front().score;
		for(predIter=predPairs[endsPair[0]][endsPair[1]-firstSuff].begin(); predIter!=predPairs[endsPair[0]][endsPair[1]-firstSuff].end(); predIter++) {
			unsigned int otherPeak;   // on input, seqsPeaks contains at most one list of peaks
			if(endsPair[0]==(*predIter)[0]) {
				newPeak=(int)(*predIter)[1]; otherPeak=endsPair[0]; curJump=spec[newPeak][0]-spec[endsPair[1]][0];
			} else {
				newPeak=(int)(*predIter)[0]; otherPeak=endsPair[1]; curJump=spec[endsPair[0]][0]-spec[newPeak][0];
			}
			float scoreAdd=0, predsPairScore=0;
			if(fabs(spec.parentMass-AAJumps::massHion*spec.idDist-spec[newPeak][0]-spec[otherPeak][0])<=2*peakTol) {
				predsPairScore=max(spec[newPeak][1],spec[otherPeak][1]); if(spec[newPeak][1]>spec[otherPeak][1]) scoreAdd = spec[newPeak][1]-spec[otherPeak][1]; else scoreAdd=0;
			} else { predsPairScore=spec[newPeak][1]+spec[otherPeak][1]; scoreAdd = spec[newPeak][1]; }

			jumpsBounds = jumps.index[(unsigned int)round(curJump/resolution)];
	if(DENOVO_DEBUG) cerr<<"  ...> expanding predecessor with mass "<<curJump<<" (jump indices ["<<jumpsBounds[0]<<":"<<jumpsBounds[1]<<"]), cumSequenceOffset = "<<cumSequenceOffset<<"\n";

			for(unsigned int jumpIdx=jumpsBounds[0]; jumpIdx<jumpsBounds[1]; jumpIdx++) {
				// Find maximum possible predecessors score within parent mass tolerance
				float maxPredScore=0;   int cumPMerr=(int)round((cumSequenceOffset+jumps[jumpIdx]-curJump)/resolution);
				int endTol = min(tolBinsOffset+pmTolRange-cumPMerr,numTolBins),
						startTol = max(0,tolBinsOffset-pmTolRange-cumPMerr),
						vSuff=(*predIter)[1]-firstSuff;
				for(int tolIdx=startTol; tolIdx<endTol; tolIdx++) if(maxPredScore<maxScores[(*predIter)[0]][vSuff][tolIdx]) maxPredScore=maxScores[(*predIter)[0]][vSuff][tolIdx];
				maxSeqScore = curScore + scoreAdd + maxPredScore - predsPairScore;
	if(DENOVO_DEBUG) cerr<<"  ...> jumps mass "<<jumps[jumpIdx]<<" results in cumMassError = "<<cumSequenceOffset+jumps[jumpIdx]-curJump<<" ("<<cumPMerr<<"), range ["<<startTol<<","<<endTol<<"], maxSeqScore = "<<maxSeqScore<<"\n";

				if(maxSeqScore<minPathScore-0.001) {
	if(DENOVO_DEBUG) cerr<<"  ...> SKIPPED!\n";
					continue;
				}

				curSeq = seqsPeaks.front();
				if(seqsExactMasses) {
					if(newPeak==(int)(*predIter)[1]) curSeq.peaks.push_back(TwoValues<float>(curSeq.peaks.back()[0]+jumps[jumpIdx],spec[newPeak][1]));
					else curSeq.peaks.push_front(TwoValues<float>(curSeq.peaks.front()[0]-jumps[jumpIdx],spec[newPeak][1]));
				} else {
					if(newPeak==(int)(*predIter)[1]) curSeq.peaks.push_back(TwoValues<float>(float(newPeak),spec[newPeak][1]));
					else curSeq.peaks.push_front(TwoValues<float>(float(newPeak),spec[newPeak][1]));
				}
	if(DENOVO_DEBUG) cerr<<"  ...> expanding predecessor of ("<<endsPair[0]<<","<<endsPair[1]<<") with jump mass "<<jumps[jumpIdx]<<", incremental score = "<<curScore+scoreAdd<<"\n";
				curSeq.score = curScore+scoreAdd;   curNewSeqs.clear();   curNewSeqs.push_back(curSeq);   // curNewScores.clear();   curNewScores.push_back(curScore+scoreAdd);

	//			float currPMerr = cumSequenceOffset+(jumps[jumpIdx]-curJump);
				float currPMerr = resolution*round((cumSequenceOffset+(jumps[jumpIdx]-curJump))/resolution);   // Sum of rounded masses guarantees same peptides as when computing DP recursion
					float tmpScore=denovo_buildSeqs_enforcePM(spec,predPairs,maxScores,*predIter,newPeak,firstSuff,tolBinsOffset,minPathScore,minPercScore,peakTol,pmOffset,pmTol,resolution,currPMerr,jumps,topSequenceOnly,curNewSeqs,seqsExactMasses);

				if(curNewSeqs.size()>0) {  // If new sequences were reported
					if(tmpScore>minPathScore) minPathScore = tmpScore;
					expandedSeqs.splice(expandedSeqs.end(),curNewSeqs);
					if(expandedSeqs.size()>SOFT_MAX_NUM_SEQS) unique_peak_lists(expandedSeqs,seqsExactMasses?resolution:0.5);
					if(expandedSeqs.size()>MAX_NUM_SEQS) break;
	//cerr<<"  --- Intermediate count = "<<expandedSeqs.size()<<" seqs\n";
				}
			}  // end for over all exact-mass jumps to the same predecessor
			if(expandedSeqs.size()>MAX_NUM_SEQS) break;
		} // end for over all predecessors

		seqsPeaks.clear();   seqsPeaks.splice(seqsPeaks.end(),expandedSeqs);
		return minPathScore;
	}

	//
	//  denovo - generates all de-novo sequences which explain at least minScorePerc
	//   percentage of the total spectrum score. Handles symmetric peaks by summing up
	//   their scores and not using the summed score twice if both peaks happen to be
	//   used - no penalty / no premium model.
	//   Non-standard amino acid mass jumps (i.e. not in jumps) are only allowed at the
	//    ends of the spectrum.
	//
	//  pmOffset=0 for PRM spectra and pmOffset=2 for MS/MS spectra.
	//  idxMatched - Indices of the matched PRMs (sparse set)
	//  minScorePerc - minimum percentage of the maximum denovo score for a sequence to be returned
	//  seqs         - all de-novo reconstructions with a total score
	//                   larger or equal to minScorePerc * sum_of_spectrum_peak_scores
	//
	//  NOTES: - Make sure that spec.idDist is adequately set.
	//         - Originally copied and adapted from alignment_scoring/getMaxSparseSet
	//
	void denovo(Spectrum &inSpec, AAJumps &jumps, float peakTol, float pmTol, float pmOffset,
									float minScorePerc, list<Spectrum> &seqs, bool enforcePM){

		seqs.clear();   if(inSpec.size()==0) return;
		Spectrum spec = inSpec;
		float aaMass = spec.parentMass+(pmOffset-AAJumps::massHion)*spec.idDist, // mass that symmetric peaks add up to
					middle = aaMass/2,   // The spectrum's axis of symmetry
					massBn = aaMass - AAJumps::massH2O,  // Mass of the last b-ion for this spectrum
					curScore,            // Score of the current path through the graph
					minPathScore=0,      // Minimum acceptable path score (updated as higher scoring paths are found)
					maxPathScore=0,      // Maximum path score found
					ionOffset=pmOffset/2;// Peak offset for the b/y-ion fragment masses

		// Includes peaks at zero/aaMass
	//	spec.addZPMpeaks(peakTol,ionOffset,true);   // spec.maximizeZPMpeaks(peakTol);
		spec.addZPMpeaks(peakTol,ionOffset,false);  // Why add the y-endpoints? The spectrum should be symmetric and
																								 //   we prefer the b-ion sequence anyway
		unsigned int numPeaks = spec.size();
		unsigned int idxPref=0, idxSuff=numPeaks-1;   // Spectrum indices of the current path endpoints

		// Set firstSuff to index of the first peak after the middle of the spectrum
		unsigned int firstSuff = min(idxPref,idxSuff); while(firstSuff<numPeaks and spec[firstSuff][0]<middle) firstSuff++;
		unsigned int suffCount = numPeaks-firstSuff;     // Number of suffix peaks - always >0 due to addZPMpeaks above

		// Declare and initialize the DP structures - 2D arrays are filled from top-right to left-bottom (increasing prefix, decreasing suffix)
		vector<vector<list<TwoValues<unsigned int> > > > predPairs;  // list of predecessor pairs for every prefix/suffix pair - pairs are in _spectrum_ indices, not matrix indices
		vector<vector<float> > maxScores;                // pos [i,j] = maximum de-novo sequence score with a prefix using peak i and a suffix using peak j
		list<TwoValues<unsigned int> > pathStarts;       // List of pairs[i,j] that connect prefixes and suffixes with a total path score larger than minPathScore

		predPairs.resize(firstSuff);   maxScores.resize(firstSuff);   pathStarts.clear();
		for(unsigned int i=0; i<firstSuff; i++) {
			predPairs[i].resize(suffCount);   maxScores[i].resize(suffCount);
			for(unsigned int j=0; j<suffCount; j++) { predPairs[i][j].clear(); maxScores[i][j]=-1; }
		}

	//cerr<<"Starting sequencing, spectrum has "<<spec.size()<<" peaks:\n";
	//for(unsigned int i=0; i<spec.size(); i++) cerr<<i<<": "<<spec[i][0]<<", "<<spec[i][1]<<endl;

		unsigned int nextPref, nextSuff, jPivot, numJumps=jumps.size(), jumpType;   idxPref=0; idxSuff=numPeaks-1;
		float curJump;   bool prefixExtension;
		// Initialize scores and predecessors of peak pairs reachable from endpoints
		for(idxPref=0; idxPref<firstSuff; idxPref++)
			for(idxSuff=numPeaks-1; idxSuff>=firstSuff; idxSuff--) {
				unsigned int vNextSuff, vSuff=idxSuff-firstSuff;  // corrected indices to access DP structures

				// Any peaks within tolerance of b0/bn are eligible starters
				if(fabs(spec[idxPref][0]-ionOffset)<=peakTol and fabs(spec[idxSuff][0]-massBn)<=max(peakTol,pmTol))
					maxScores[idxPref][vSuff] = spec[idxPref][1]+spec[idxSuff][1];

				if(maxScores[idxPref][vSuff]<0) {
	//cerr<<"("<<idxPref<<","<<idxSuff<<"/"<<vSuff<<","<<maxScores[idxPref][vSuff]<<"), skipped...\n";
					continue; // Skip unreachable pairs of peaks
				}

				// Try both prefix and suffix extensions - must jump _over_ symmetric suffix/prefix, respectively
				for(jumpType=0; jumpType<2; jumpType++) {
					prefixExtension = (bool) jumpType;

					// Initialize depending on whether this is a prefix or a suffix extension - jumping forward towards the middle and over the complement of the pair
					if(prefixExtension) {
						nextSuff=idxSuff;  nextPref=idxPref+1;  while(nextPref<firstSuff and spec[nextPref][0]<=aaMass-spec[nextSuff][0]+peakTol) nextPref++;
						curJump=spec[nextPref][0]-spec[idxPref][0];
						vNextSuff = vSuff;
					} else {
						nextPref=idxPref;  nextSuff=idxSuff-1;  while(nextSuff>=firstSuff and spec[nextSuff][0]>=aaMass-spec[nextPref][0]-peakTol) nextSuff--;
						curJump=spec[idxSuff][0]-spec[nextSuff][0];
						if (nextSuff>firstSuff) vNextSuff = nextSuff - firstSuff; else vNextSuff = 0;
					}

					// iterate over all reachable peaks
					jPivot=0;
					while(jPivot<numJumps) {
						if(prefixExtension) {  // Find nextPrev by valid jumps from idxPrev
							while(jPivot<numJumps and (fabs(curJump-jumps[jPivot])>peakTol or (nextPref>=firstSuff and nextPref<idxSuff))) {
								if(fabs(curJump-jumps[jPivot])>peakTol and curJump>jumps[jPivot]) jPivot++; else { if(nextPref<idxSuff){ nextPref++; curJump=spec[nextPref][0]-spec[idxPref][0]; } else break; }
							}
							if(jPivot==numJumps or fabs(curJump-jumps[jPivot])>peakTol) break;
						} else {
							// Find nextSuff by valid jumps from idxSuff
							while(jPivot<numJumps and (fabs(curJump-jumps[jPivot])>peakTol or (nextSuff>idxPref and nextSuff<firstSuff))){
								if(fabs(curJump-jumps[jPivot])>peakTol and curJump>jumps[jPivot]) jPivot++; else { if(nextSuff>idxPref){ nextSuff--; curJump=spec[idxSuff][0]-spec[nextSuff][0]; } else break; }
						}
							if(jPivot==numJumps or fabs(curJump-jumps[jPivot])>peakTol) break;
							if (nextSuff>firstSuff) vNextSuff = nextSuff - firstSuff; else vNextSuff = 0;
						}

	//if(idxPref==37) { cerr<<" --- cur = ("<<idxPref<<","<<idxSuff<<"), next = ("<<nextPref<<","<<nextSuff<<")\n"; cerr.flush(); }

						if((nextPref==idxSuff or nextSuff==idxPref) and jPivot<numJumps and (fabs(curJump-jumps[jPivot])<=peakTol)) { // connecting jump
							if(enforcePM or maxScores[idxPref][vSuff]>=minPathScore) {  // if enforcePM then all possible connecting sequences have to be tried to check whether the parent mass matches
								if(prefixExtension) pathStarts.push_front(TwoValues<unsigned int>(idxPref,idxSuff));  // No need to add same start from prefix AND suffix extensions
	//cerr<<"New pathStart at maxScores["<<nextPref<<","<<nextSuff<<"] = "<<maxScores[idxPref][vSuff]<<" from "<<"["<<idxPref<<","<<idxSuff<<"]\n";
								if(maxScores[idxPref][vSuff]>maxPathScore) {
									maxPathScore = maxScores[idxPref][vSuff];
									minPathScore = minScorePerc*maxPathScore;
	//cerr<<"New maxScore at ["<<nextPref<<","<<nextSuff<<"] = "<<maxPathScore<<" from "<<"["<<idxPref<<","<<idxSuff<<"]\n";
								}
							}
							break;
						}

						// Add [idxPrev,idxSuff] to list of predecessors/predScores
						predPairs[nextPref][vNextSuff].push_front(TwoValues<unsigned int>(idxPref,idxSuff));

						// Adjust maximum score of destination pair
						if(fabs(aaMass-spec[nextPref][0]-spec[nextSuff][0])>2*peakTol)
							curScore=maxScores[idxPref][vSuff] + (prefixExtension?spec[nextPref][1]:spec[nextSuff][1]);
	//					else curScore=maxScores[idxPref][vSuff];   // Don't double-count
						else curScore=maxScores[idxPref][vSuff] - (prefixExtension?spec[nextSuff][1]:spec[nextPref][1]) + (spec[nextSuff][1]+spec[nextPref][1])/2;   // Using symmetric peaks adds the average score of the peak pair. (same as not double-counting for symmetric spectra)
	//					else curScore=maxScores[idxPref][vSuff] - (prefixExtension?spec[nextSuff][1]:spec[nextPref][1]) + max(spec[nextSuff][1],spec[nextPref][1]);   // Using symmetric peaks adds the max score of the peak pair. (same as not double-counting for symmetric spectra)
						if(maxScores[nextPref][vNextSuff]<curScore) {
	//if(nextPref==24 or nextPref==43) cerr<<"Prefix extension to ["<<nextPref<<","<<nextSuff<<"] from ["<<idxPref<<","<<idxSuff<<"], curScore = "<<curScore<<"\n";
							maxScores[nextPref][vNextSuff] = curScore;
						}

	//cerr<<" -> ("<<nextPref<<","<<nextSuff<<"/"<<vNextSuff<<","<<maxScores[nextPref][vNextSuff]<<"), ("
	//    <<spec[nextPref][0]<<","<<spec[nextSuff][0]<<")\n";

						if(prefixExtension) { nextPref++; if(nextPref<=idxSuff) curJump=spec[nextPref][0]-spec[idxPref][0]; else break; }
						else { nextSuff--; if(nextSuff>=idxPref) curJump=spec[idxSuff][0]-spec[nextSuff][0]; else break; }
					}  // End while(jPivot<numJumps)

	//if(idxPref==8 and idxSuff==25) {
	//cerr<<"--->>> ("<<idxPref<<","<<idxSuff<<"/"<<vSuff<<","<<maxScores[idxPref][vSuff]<<"), done\n";
	//}
					// Set predecessor/scores for all other peaks reachable from suffix extensions with large jumps
					//  -- note that this is not acceptable when exact parent masses are required (exact jump mass is unknown)
					if(not enforcePM and not prefixExtension and fabs(spec[idxPref][0]-ionOffset)<=peakTol)
						while(nextSuff>=firstSuff) {
							vNextSuff = nextSuff - firstSuff;
							curScore = maxScores[idxPref][vSuff]+spec[nextSuff][1];
							predPairs[nextPref][vNextSuff].push_front(TwoValues<unsigned int>(idxPref,idxSuff));
							if(curScore>maxScores[nextPref][vNextSuff]) maxScores[nextPref][vNextSuff] = curScore;
							nextSuff--;
						}

					// Set predecessor/scores for all other peaks reachable from prefix extensions with large jumps
					//  -- note that this is not acceptable when exact parent masses are required (exact jump mass is unknown)
					if(not enforcePM and prefixExtension and fabs(spec[idxSuff][0]-massBn)<=pmTol)
						while(nextPref<firstSuff) {
							curScore = maxScores[idxPref][vSuff]+spec[nextPref][1];
							predPairs[nextPref][vNextSuff].push_front(TwoValues<unsigned int>(idxPref,idxSuff));
							if(curScore>maxScores[nextPref][vNextSuff]) maxScores[nextPref][vNextSuff] = curScore;
							nextPref++;
						}
				}
			}

		list<TwoValues<unsigned int> >::iterator iter = pathStarts.begin();
		while (not enforcePM and iter!=pathStarts.end())
			if(maxScores[(*iter)[0]][(*iter)[1]-firstSuff]<minPathScore-0.0001) iter=pathStarts.erase(iter); else iter++;

	//cerr<<"Got "<<pathStarts.size()<<" possible starts. enforcePM = "<<enforcePM<<", pmTol = "<<pmTol<<"\n"; cerr.flush();
		list<list<TwoValues<float> > > curSeqs;   list<float> curScores;   Spectrum cSpec;    cSpec=spec;
		list<list<TwoValues<float> > > allSeqs;   list<float> allScores;
		if(enforcePM and minScorePerc>=1-0.0001) minPathScore=0; else minPathScore=minPathScore*.9999;   // minPathScore*.9999 because of numeric rounding errors
		module_num_rejected_seqs = 0;
		allSeqs.clear();   allScores.clear();
		for(list<TwoValues<unsigned int> >::iterator iter=pathStarts.begin(); iter!=pathStarts.end(); iter++) {
			if(maxScores[(*iter)[0]][(*iter)[1]-firstSuff]<minPathScore) continue;
	//cerr<<"Considering pathStart at ["<<(*iter)[0]<<","<<(*iter)[1]<<"] with maxScore "<<maxScores[(*iter)[0]][(*iter)[1]-firstSuff]<<" (minPathScore="<<minPathScore<<")\n";
			if(enforcePM) {
				TwoValues<unsigned int> jumpBounds;
				jumps.find(spec[(*iter)[1]][0]-spec[(*iter)[0]][0],peakTol,jumpBounds);
	//cerr<<" --- Query mass "<<spec[(*iter)[1]][0]-spec[(*iter)[0]][0]<<" got hits to ["<<jumpBounds[0]<<":"<<jumpBounds[1]<<"].\n";
				for(unsigned int jumpIdx=jumpBounds[0]; jumpIdx<jumpBounds[1]; jumpIdx++) {
	//cerr<<" --- starting with middle aa mass = "<<jumps[jumpIdx]<<", minPathScore = "<<minPathScore<<", cumError = "<<spec[(*iter)[1]][0]-spec[(*iter)[0]][0]-jumps[jumpIdx]<<"\n"; cerr.flush();
					curSeqs.clear();   curScores.clear();
					minPathScore=denovo_buildSeqs_enforcePM(spec,predPairs,maxScores,*iter,-1,firstSuff,minPathScore,peakTol,pmOffset,pmTol,jumps[jumpIdx]-(spec[(*iter)[1]][0]-spec[(*iter)[0]][0]),jumps,minScorePerc>=1-0.0001,curSeqs,curScores);
	//cerr<<" --- minPathScore = "<<minPathScore<<"\n"; cerr.flush();
					allSeqs.splice(allSeqs.end(), curSeqs);    allScores.splice(allScores.end(), curScores);
				}
	//cerr<<" ++ Got "<<curSeqs.size()<<" seqs, total count = "<<curSeqs.size()+allSeqs.size()<<".\n"; cerr.flush();
			} else {
				curSeqs.clear();   curScores.clear();
				denovo_buildSeqs(spec,predPairs,maxScores,*iter,-1,firstSuff,minPathScore,peakTol,pmOffset,curSeqs,curScores);
				allSeqs.splice(allSeqs.end(), curSeqs);    allScores.splice(allScores.end(), curScores);
			}

		}

		// Returns only those sequences observing the final minPathScore threshold
	//cerr<<" >> Got "<<allSeqs.size()<<" seqs.\n"; cerr.flush();
		list<float>::iterator scoresIter = allScores.begin();
	//	for(list<list<unsigned int> >::iterator seqsIter=allSeqs.begin(); seqsIter!=allSeqs.end(); seqsIter++, scoresIter++)
		for(list<list<TwoValues<float> > >::iterator seqsIter=allSeqs.begin(); seqsIter!=allSeqs.end(); seqsIter++, scoresIter++) {
	//cerr<<" --- current score = "<<*scoresIter<<", minPathScore = "<<minPathScore<<" ("<<(*scoresIter>=minPathScore-0.0001)<<"): ";
			if(*scoresIter>=minPathScore-0.0001) {
	//cerr<<" accepted\n";
				cSpec.resize(seqsIter->size());  int peakIdx=0;
	//			for(list<unsigned int>::iterator peakIter=seqsIter->begin(); peakIter!=seqsIter->end(); peakIter++)
				for(list<TwoValues<float> >::iterator peakIter=seqsIter->begin(); peakIter!=seqsIter->end(); peakIter++) {
	//				cSpec[peakIdx++]=spec[*peakIter];
					cSpec[peakIdx++]=*peakIter;
				}
				seqs.push_back(cSpec);
			}
	//else cerr<<" rejected\n";
		}
	//cerr<<" >> Returning "<<seqs.size()<<" seqs.\n"; cerr.flush();
	}

	//
	//  denovo_exactPM - generates all de-novo sequences from peptides within pmTol
	//   Daltons of the experimental parent mass (strictly enforced) which explain at
	//   least minScorePerc percentage of the total spectrum score. Handles symmetric
	//   peaks by summing up their scores and not using the summed score twice if
	//   both peaks happen to be used - no penalty / no premium model.
	//   Non-standard amino acid mass jumps (i.e. not in jumps) are only allowed at the
	//    ends of the spectrum.
	//
	//  pmOffset=0 for PRM spectra and pmOffset=2 for MS/MS spectra.
	//  idxMatched - Indices of the matched PRMs (sparse set)
	//  minScorePerc - minimum percentage of the maximum denovo score for a sequence to be returned
	//  seqs         - all de-novo reconstructions with a total score
	//                   larger or equal to minScorePerc * sum_of_spectrum_peak_scores
	//
	//  NOTES: - Make sure that spec.idDist is adequately set.
	//         - Originally copied and adapted from denovo() above
	//
	float denovo_exactPM(Spectrum &inSpec, AAJumps &jumps, float peakTol, float pmTol, float pmOffset,
												 float resolution, float minScorePerc, list<Spectrum> &seqs,
												 bool seqsExactMasses, float minAbsScore){
		seqs.clear();   if(inSpec.size()==0) return 0;
		Spectrum spec = inSpec;
		float aaMass = spec.parentMass+(pmOffset-AAJumps::massHion)*spec.idDist, // mass that symmetric peaks add up to
					middle = aaMass/2,        // The spectrum's axis of symmetry
					massBn = aaMass - AAJumps::massH2O,  // Mass of the last b-ion for this spectrum
					curScore,                 // Score of the current path through the graph
					minPathScore=minAbsScore, // Minimum acceptable path score (updated as higher scoring paths are found)
					maxPathScore=0,           // Maximum path score found
					ionOffset=pmOffset/2;     // Peak offset for the b/y-ion fragment masses
		const int tolBinsOffset = (int)ceil(max(pmTol,peakTol)/resolution), pmTolRange=(int)round(pmTol/resolution);
		const int numTolBins = 2*tolBinsOffset+1;

		// Includes peaks at zero/aaMass and make sure all have the same score
	//	spec.addZPMpeaks(peakTol,ionOffset,false);
		spec.addZPMpeaks(0,ionOffset,false);  // peakTol tolerance may cause an artificially large cumulative parent mass error and exclude the correct peptide

		if(jumps.index.size()==0) jumps.computeIndex(peakTol,resolution);
		unsigned int numPeaks = spec.size();
		unsigned int idxPref=0, idxSuff=numPeaks-1;   // Spectrum indices of the current path endpoints

		// Set firstSuff to index of the first peak after the middle of the spectrum
		unsigned int firstSuff = min(idxPref,idxSuff); while(firstSuff<numPeaks and spec[firstSuff][0]<middle) firstSuff++;
		unsigned int suffCount = numPeaks-firstSuff;     // Number of suffix peaks - always >0 due to addZPMpeaks above

		// Declare and initialize the DP structures - 2D arrays are filled from top-right to left-bottom (increasing prefix, decreasing suffix)
		vector<vector<list<TwoValues<unsigned int> > > > predPairs;  // list of predecessor pairs for every prefix/suffix pair - pairs are in _spectrum_ indices, not matrix indices
		vector<vector<vector<float> > > maxScores;       // pos [i,j,k] = maximum de-novo sequence score with a prefix ending at peak i, a suffix starting at peak j
																										 //   and with a peptide whose parent mass is k*resolution-max(pmTol,peakTol) Da over the experimental parent mass
		list<TwoValues<unsigned int> > pathStarts;       // List of pairs[i,j] that connect prefixes and suffixes with a total path score larger than minPathScore

		predPairs.resize(firstSuff);   maxScores.resize(firstSuff);   pathStarts.clear();
		for(unsigned int i=0; i<firstSuff; i++) {
			predPairs[i].resize(suffCount);   maxScores[i].resize(suffCount);
			for(unsigned int j=0; j<suffCount; j++) {
				predPairs[i][j].clear();
				maxScores[i][j].resize(numTolBins);
				for(int k=0; k<numTolBins; k++) maxScores[i][j][k]=-1;
			}
		}

	if(DENOVO_DEBUG) {
		cerr<<"Starting sequencing, spectrum has "<<spec.size()<<" peaks, parent mass = "<<spec.parentMass<<" (massBn = "<<massBn<<"):\n";
		for(unsigned int i=0; i<spec.size(); i++) cerr<<i<<": "<<spec[i][0]<<", "<<spec[i][1]<<endl;
	}
		unsigned int nextPref, nextSuff, jPivot, jumpType, vNextSuff, vSuff;   idxPref=0; idxSuff=numPeaks-1;
		int idxPMerr,     // Corresponding index for theoretical peptide mass minus experimental parent mass
				addPMerr;
		int endPMerr, startPMerr;  // Set the limits for iterating over parent mass error bins
		float curJump, maxScore;   TwoValues<unsigned int> curJumpIdx;
		bool prefixExtension;
		// Initialize scores and predecessors of peak pairs reachable from endpoints
		for(idxPref=0; idxPref<firstSuff; idxPref++)
			for(idxSuff=numPeaks-1; idxSuff>=firstSuff; idxSuff--) {
				vNextSuff=0; vSuff=idxSuff-firstSuff;  // corrected indices to access DP structures
				bool reachable = false;                // Indicates whether pair [idxPref,idxSuff] is reachable

				// Any peaks within tolerance of b0/bn are eligible starters
				if(fabs(spec[idxPref][0]-ionOffset)<=peakTol and fabs(spec[idxSuff][0]-massBn)<=max(peakTol,pmTol)) {
					idxPMerr = tolBinsOffset + (int)round((spec[idxSuff][0]-massBn-spec[idxPref][0])/resolution);
					if(idxPMerr>=0 and idxPMerr<numTolBins)
						{ reachable=true; maxScores[idxPref][vSuff][idxPMerr] = spec[idxPref][1]+spec[idxSuff][1]; }
				}

				if(not reachable)
					for(idxPMerr=0;idxPMerr<numTolBins;idxPMerr++)
						if(maxScores[idxPref][vSuff][idxPMerr]>=0) { reachable=true; break; }

				if(not reachable) {
	//cerr<<"("<<idxPref<<","<<idxSuff<<"/"<<vSuff<<"), skipped...\n";
					continue; // Skip unreachable pairs of peaks
				}

				// Try both prefix and suffix extensions - must jump _over_ symmetric suffix/prefix, respectively
				for(jumpType=0; jumpType<2; jumpType++) {
					prefixExtension = (bool) jumpType;

					// Initialize depending on whether this is a prefix or a suffix extension - jumping forward towards the middle and over the complement of the pair
					if(prefixExtension) {
						nextSuff=idxSuff;  nextPref=idxPref+1;  while(nextPref<firstSuff and spec[nextPref][0]<=aaMass-spec[nextSuff][0]+peakTol) nextPref++;
						curJump=spec[nextPref][0]-spec[idxPref][0];
						vNextSuff = vSuff;
					} else {
						nextPref=idxPref;  nextSuff=idxSuff-1;  while(nextSuff>=firstSuff and spec[nextSuff][0]>=aaMass-spec[nextPref][0]-peakTol) nextSuff--;
						curJump=spec[idxSuff][0]-spec[nextSuff][0];
						if (nextSuff>=firstSuff) vNextSuff = nextSuff - firstSuff; else continue;  // Connecting jumps only allowed via prefix extensions to avoid repetition
					}

					// iterate over all eligible extensions
					while(true) {
						jPivot = (unsigned int)round(curJump/resolution);   if(jPivot<jumps.index.size()) curJumpIdx = jumps.index[jPivot];
						if(prefixExtension) {  // Find nextPrev by valid jumps from idxPrev
							while(jPivot<jumps.index.size() and (curJumpIdx[0]==curJumpIdx[1] or (nextPref>=firstSuff and nextPref<idxSuff))) {
								if(nextPref>=idxSuff) break;
								if(nextPref>=firstSuff and nextPref<idxSuff) nextPref=idxSuff; else nextPref++;
								curJump=spec[nextPref][0]-spec[idxPref][0];
								jPivot = (unsigned int)round(curJump/resolution);   if(jPivot<jumps.index.size()) curJumpIdx = jumps.index[jPivot];
							}
							if(jPivot>=jumps.index.size() or curJumpIdx[0]==curJumpIdx[1]) break;
						} else {
							while(jPivot<jumps.index.size() and curJumpIdx[0]==curJumpIdx[1] and nextSuff>firstSuff) {  // Last condition enforces connecting jumps only through prefix extensions (no need to repeat connection)
								nextSuff--;   curJump = spec[idxSuff][0]-spec[nextSuff][0];
								jPivot = (unsigned int)round(curJump/resolution);   if(jPivot<jumps.index.size()) curJumpIdx = jumps.index[jPivot];
							}
							if(jPivot>=jumps.index.size() or curJumpIdx[0]==curJumpIdx[1]) break;
							vNextSuff = nextSuff - firstSuff;
						}

	//cerr<<" -> Attempting extension to ("<<nextPref<<","<<nextSuff<<") with mass difference "<<curJump<<", jPivot = "<<jPivot<<", bounds ["<<curJumpIdx[0]<<","<<curJumpIdx[1]<<"]\n";

						maxScore=-1;   startPMerr=numTolBins;   bool predPairsDone=false;   int maxScorePMidx=-1;
						float scoreIncrement;
						if(fabs(aaMass-spec[nextPref][0]-spec[nextSuff][0])>2*peakTol)
							scoreIncrement = (prefixExtension?spec[nextPref][1]:spec[nextSuff][1]);
						else scoreIncrement = -(prefixExtension?spec[nextSuff][1]:spec[nextPref][1]) + (spec[nextSuff][1]+spec[nextPref][1])/2;   // Using symmetric peaks adds the average score of the peak pair.
	//					else scoreIncrement= - (prefixExtension?spec[nextSuff][1]:spec[nextPref][1]) + max(spec[nextSuff][1],spec[nextPref][1]);   // Using symmetric peaks adds the max score of the peak pair. (same as not double-counting for symmetric spectra)
						for(unsigned int jCur=curJumpIdx[0]; jCur<curJumpIdx[1]; jCur++) {
							addPMerr = (int)round((jumps[jCur]-curJump)/resolution);

	//cerr<<" -->> jumpIdx = "<<jCur<<" (mass "<<jumps[jCur]<<"), addPMerr = "<<addPMerr<<"\n";

							if(nextPref==idxSuff or nextSuff==idxPref) {
								// Find maximum de-novo sequence score with a prefix ending at peak i, a suffix starting at peak j and a peptide mass within tolerance
								endPMerr = (int)min(startPMerr,tolBinsOffset+pmTolRange-addPMerr+1);  // Since curJump is larger than previous jump there is no need to search for a better solution beyond the previous startPMerr
								startPMerr = (int)max(0,tolBinsOffset-pmTolRange-addPMerr);
								for(idxPMerr=startPMerr; idxPMerr<endPMerr; idxPMerr++)
									if(maxScores[idxPref][vSuff][idxPMerr]>maxScore)
										{ maxScore=maxScores[idxPref][vSuff][idxPMerr]; maxScorePMidx=idxPMerr; }
							} else {
								// Find mass tolerance range in destination pair matching addPMerr and update maxScores/predPairs
								endPMerr = numTolBins+(int)min(0,-addPMerr);
								startPMerr = max(0,-addPMerr);

	//cerr<<" -->> Updating destination ("<<nextPref<<","<<nextSuff<<") from ("<<idxPref<<","<<idxSuff<<") for jump ("<<jCur<<","<<jumps[jCur]<<") in ["<<startPMerr<<","<<endPMerr<<"], predPairsDone = "<<predPairsDone<<"\n";

								for(idxPMerr=startPMerr; idxPMerr<endPMerr; idxPMerr++) {
									// Add [idxPrev,idxSuff] to list of predecessors/predScores
									if(not predPairsDone) { predPairs[nextPref][vNextSuff].push_front(TwoValues<unsigned int>(idxPref,idxSuff)); predPairsDone=true; }

									// Adjust maximum score of destination pair
	/*								if(fabs(aaMass-spec[nextPref][0]-spec[nextSuff][0])>2*peakTol)
										curScore=maxScores[idxPref][vSuff][idxPMerr] + (prefixExtension?spec[nextPref][1]:spec[nextSuff][1]);
									else curScore=maxScores[idxPref][vSuff][idxPMerr] - (prefixExtension?spec[nextSuff][1]:spec[nextPref][1]) + (spec[nextSuff][1]+spec[nextPref][1])/2;   // Using symmetric peaks adds the average score of the peak pair.
	//								else curScore=maxScores[idxPref][vSuff] - (prefixExtension?spec[nextSuff][1]:spec[nextPref][1]) + max(spec[nextSuff][1],spec[nextPref][1]);   // Using symmetric peaks adds the max score of the peak pair. (same as not double-counting for symmetric spectra)
	*/
									curScore = maxScores[idxPref][vSuff][idxPMerr] + scoreIncrement;
									if(maxScores[nextPref][vNextSuff][idxPMerr+addPMerr]<curScore) maxScores[nextPref][vNextSuff][idxPMerr+addPMerr] = curScore;
								}
							}
						}

						if((nextPref==idxSuff or nextSuff==idxPref) and maxScore>=minPathScore) { // connecting jump
							if(fabs(spec[idxPref][0]-ionOffset)>peakTol or fabs(spec[idxSuff][0]-massBn)>max(peakTol,pmTol)) {  // Don't connect zero to parent mass (valid but useless)
								pathStarts.push_front(TwoValues<unsigned int>(idxPref,idxSuff));
								if(maxScore>maxPathScore) {
									maxPathScore = maxScore;
									minPathScore = minScorePerc*maxPathScore;
	//if(DENOVO_DEBUG) cerr<<"New maxScore at ["<<nextPref<<","<<nextSuff<<"] = "<<maxPathScore<<" from "<<"["<<idxPref<<","<<idxSuff<<"]\n";
								}
	if(DENOVO_DEBUG) cerr<<"New connecting jump with maxScore["<<nextPref<<","<<nextSuff<<"] = "<<maxPathScore<<" from "<<"["<<idxPref<<","<<idxSuff<<"], addPMerr = "<<addPMerr<<", maxScorePMidx = "<<maxScorePMidx<<"\n";
							}
						}

	//cerr<<" -> ("<<nextPref<<","<<nextSuff<<"/"<<vNextSuff<<","<<maxScores[nextPref][vNextSuff]<<"), ("
	//    <<spec[nextPref][0]<<","<<spec[nextSuff][0]<<")\n";
						if(prefixExtension) { nextPref++; if(nextPref>idxSuff) break; curJump=spec[nextPref][0]-spec[idxPref][0]; }
						else { nextSuff--; if(nextSuff<firstSuff) break; curJump=spec[idxSuff][0]-spec[nextSuff][0]; }
					}  // End while(true)

	//if(idxPref==8 and idxSuff==25) {
	//cerr<<"--->>> ("<<idxPref<<","<<idxSuff<<"/"<<vSuff<<","<<maxScores[idxPref][vSuff]<<"), done\n";
	//}
	/* This code became irrelevant - large jumps are now allowed everywhere
					// Set predecessor/scores for all other peaks reachable from suffix extensions with large jumps
					//  -- note that this is not acceptable when exact parent masses are required (exact jump mass is unknown)
					if(not enforcePM and not prefixExtension and fabs(spec[idxPref][0]-ionOffset)<=peakTol)
						while(nextSuff>=firstSuff) {
							vNextSuff = nextSuff - firstSuff;
							curScore = maxScores[idxPref][vSuff]+spec[nextSuff][1];
							predPairs[nextPref][vNextSuff].push_front(TwoValues<unsigned int>(idxPref,idxSuff));
							if(curScore>maxScores[nextPref][vNextSuff]) maxScores[nextPref][vNextSuff] = curScore;
							nextSuff--;
						}

					// Set predecessor/scores for all other peaks reachable from prefix extensions with large jumps
					//  -- note that this is not acceptable when exact parent masses are required (exact jump mass is unknown)
					if(not enforcePM and prefixExtension and fabs(spec[idxSuff][0]-massBn)<=pmTol)
						while(nextPref<firstSuff) {
							curScore = maxScores[idxPref][vSuff]+spec[nextPref][1];
							predPairs[nextPref][vNextSuff].push_front(TwoValues<unsigned int>(idxPref,idxSuff));
							if(curScore>maxScores[nextPref][vNextSuff]) maxScores[nextPref][vNextSuff] = curScore;
							nextPref++;
						}
	*/
				}
			}

	/*cerr<<"maxScores:\n";
	for(idxPref=0; idxPref<firstSuff; idxPref++) {
		for(idxSuff=0; idxSuff<firstSuff; idxSuff++)
			cerr<<maxScores[idxPref][idxSuff]<<"\t";
		cerr<<endl;
	}*/

	/* Will be checked below - no need to fix/update
		list<TwoValues<unsigned int> >::iterator iter = pathStarts.begin();
		while (iter!=pathStarts.end()) {
			if(maxScores[(*iter)[0]][(*iter)[1]-firstSuff]<minPathScore-0.0001)
				iter=pathStarts.erase(iter); else iter++;
		}
	*/

	if(DENOVO_DEBUG) cerr<<"Got "<<pathStarts.size()<<" possible starts (exactPM), pmTol = "<<pmTol<<", min/maxPathScore = "<<minPathScore<<"/"<<maxPathScore<<"\n"; cerr.flush();
		list<DenovoSeq> curSeqs;   Spectrum cSpec;    cSpec=spec;
		list<DenovoSeq> allSeqs;   allSeqs.clear();
		list<DenovoSeq>::iterator seqsIter;
		module_num_rejected_seqs = 0;
		float newMinPathScore;
		for(list<TwoValues<unsigned int> >::iterator iter=pathStarts.begin(); iter!=pathStarts.end(); iter++) {
				curJumpIdx = jumps.index[(unsigned int)round((spec[(*iter)[1]][0]-spec[(*iter)[0]][0])/resolution)];
	if(DENOVO_DEBUG) cerr<<"-> Query mass "<<spec[(*iter)[1]][0]-spec[(*iter)[0]][0]<<", pair ("<<(*iter)[0]<<","<<(*iter)[1]<<"), got hits to ["<<curJumpIdx[0]<<":"<<curJumpIdx[1]<<"].\n";
				for(unsigned int jumpIdx=curJumpIdx[0]; jumpIdx<curJumpIdx[1]; jumpIdx++) {

					float cumSequenceOffset = jumps[jumpIdx]-(spec[(*iter)[1]][0]-spec[(*iter)[0]][0]);
					addPMerr=(int)round(cumSequenceOffset/resolution);

					endPMerr   = min(tolBinsOffset+pmTolRange-addPMerr,numTolBins);
					startPMerr = max(0,tolBinsOffset-pmTolRange-addPMerr);
					vSuff = (*iter)[1]-firstSuff;
					maxScore=0;	for(int idxPMerr=startPMerr; idxPMerr<endPMerr; idxPMerr++) if(maxScore<maxScores[(*iter)[0]][vSuff][idxPMerr]) maxScore=maxScores[(*iter)[0]][vSuff][idxPMerr];
	if(DENOVO_DEBUG) cerr<<"-->> middle aa mass = "<<jumps[jumpIdx]<<", maxScore = "<<maxScore<<"/"<<minPathScore<<", range = ["<<startPMerr<<","<<endPMerr<<"], cumError = "<<cumSequenceOffset<<" ("<<addPMerr<<")\n"; cerr.flush();
					if(maxScore<minPathScore-0.00001) {
	if(DENOVO_DEBUG) cerr<<"-->> SKIPPED!\n"; cerr.flush();
						continue;
					}

	if(DENOVO_DEBUG) cerr<<"-->> starting with middle aa mass = "<<jumps[jumpIdx]<<", maxScore = "<<maxScore<<"/"<<minPathScore<<", cumError = "<<jumps[jumpIdx]-(spec[(*iter)[1]][0]-spec[(*iter)[0]][0])<<"\n"; cerr.flush();
					curSeqs.clear();
	//				newMinPathScore=denovo_buildSeqs_enforcePM(spec,predPairs,maxScores,*iter,-1,firstSuff,tolBinsOffset,minPathScore,minScorePerc,peakTol,pmOffset,pmTol,resolution,cumSequenceOffset,jumps,minScorePerc>=1-0.0001,curSeqs,curScores);
					newMinPathScore=denovo_buildSeqs_enforcePM(spec,predPairs,maxScores,*iter,-1,firstSuff,tolBinsOffset,minPathScore,minScorePerc,peakTol,pmOffset,pmTol,resolution,resolution*(float)addPMerr,jumps,minScorePerc>=1-0.0001,curSeqs,seqsExactMasses);

	if(DENOVO_DEBUG) cerr<<"++++ Got "<<curSeqs.size()<<" new seqs, newMinPathScore = "<<newMinPathScore<<", minPathScore = "<<minPathScore; cerr.flush();
					allSeqs.splice(allSeqs.end(), curSeqs);
					if(newMinPathScore>minPathScore) {
						seqsIter=allSeqs.begin();
						while(seqsIter!=allSeqs.end())
							if(seqsIter->score>=newMinPathScore) seqsIter++; else seqsIter=allSeqs.erase(seqsIter);
						minPathScore = newMinPathScore;
					}
	if(DENOVO_DEBUG) cerr<<", # seqs = "<<allSeqs.size(); cerr.flush();
					if(allSeqs.size()>SOFT_MAX_NUM_SEQS) unique_peak_lists(allSeqs,seqsExactMasses?resolution:0.5);
	if(DENOVO_DEBUG) cerr<<" ("<<allSeqs.size()<<" unique)\n"; cerr.flush();
					if(allSeqs.size()>MAX_NUM_SEQS) break;
				}
			if(allSeqs.size()>MAX_NUM_SEQS) break;
		}
	//cerr<<" -- rejected "<<module_num_rejected_seqs<<" with parent mass error > "<<pmTol<<"\n"; cerr.flush();

		// Quick-not-so-elegant fix to remove repeated sequences possibly from different peptides within mass tolerance
		//   or generated from symmetric peaks (A,B) via predecessors (I,B), (A,J), both reachable from (I,J)
		unique_peak_lists(allSeqs,seqsExactMasses?resolution:0.5);

	/*if(allSeqs.size()>1) {
		for(unsigned int i=0; i<spec.size(); i++) cerr<<i<<": "<<spec[i][0]<<", "<<spec[i][1]<<endl;
		for(seqsIter=allSeqs.begin(); seqsIter!=allSeqs.end(); seqsIter++) {
			cerr<<"Seq: "; cerr.flush();
			for(list<TwoValues<float> >::iterator peakIter=seqsIter->peaks.begin(); peakIter!=seqsIter->peaks.end(); peakIter++)
				cerr<<(*peakIter)[0]<<" ";
			cerr<<endl; cerr.flush();
		}
	}*/

		// Convert the output to a list of Spectrum
	//cerr<<" >> Got "<<allSeqs.size()<<" seqs.\n"; cerr.flush();
	//	for(list<list<unsigned int> >::iterator seqsIter=allSeqs.begin(); seqsIter!=allSeqs.end(); seqsIter++, scoresIter++)
		for(seqsIter=allSeqs.begin(); seqsIter!=allSeqs.end(); seqsIter++) {
			cSpec.resize(seqsIter->size());  int peakIdx=0;
			for(list<TwoValues<float> >::iterator peakIter=seqsIter->peaks.begin(); peakIter!=seqsIter->peaks.end(); peakIter++)
				if(seqsExactMasses) cSpec[peakIdx++]=*peakIter; else cSpec[peakIdx++]=spec[(int)round((*peakIter)[0])];
			seqs.push_back(cSpec);
		}

	if(DENOVO_DEBUG) cerr<<" >> Returning "<<seqs.size()<<" seqs.\n"; cerr.flush();
		return minPathScore;
	}

	//
	//  denovo - generates all de-novo sequences which explain at least minScorePerc
	//   percentage of the total spectrum score. Handles symmetric peaks by summing up
	//   their scores and not using the summed score twice if both peaks happen to be
	//   used - no penalty / no premium model.
	//   Non-standard amino acid mass jumps (i.e. not in neighsL/neighsR) are only
	//    allowed at the ends of the spectrum.
	//
	//  ionOffset - mass offset from b/y ions in relation to PRM masses: 0 for PRM spectra and AAMasses::massHion for MS/MS spectra.
	//  pmOffset  - parent mass offset such that mass_prefix_peak+mass_suffix_peak = spec.parentMass-AAMasses::massMH
	//                Usually, =0 for PRM spectra and pmOffset=2 for MS/MS spectra.
	//  idxMatched - Indices of the matched PRMs (sparse set)
	//  minScorePerc - minimum percentage of the maximum denovo score for a sequence to be returned
	//  seqs         - all de-novo reconstructions with a total score
	//                   larger or equal to minScorePerc * sum_of_spectrum_peak_scores
	//
	//  NOTES: - Make sure that spec.idDist is adequately set.
	//         - neighsL/neighsR-version of the denovo function above
	//
	void denovo(Spectrum &inSpec, vector<vector<int> > &neighsL, vector<vector<int> > &neighsR,
					float peakTol, float pmTol, float ionOffset, float minScorePerc,
					list<Spectrum> &seqs, list<float> &seqScores, bool ctermH2O, bool enforcePM) {

		seqs.clear();   if(inSpec.size()==0) return;
		Spectrum spec = inSpec;
		float pmOffset = 2*ionOffset,
				aaMass = spec.parentMass+(pmOffset-AAJumps::massHion)*spec.idDist, // mass that symmetric peaks add up to
					middle = aaMass/2,   // The spectrum's axis of symmetry
					massBn = aaMass-ionOffset-(ctermH2O?AAJumps::massH2O:0),  // Mass of the last b-ion for this spectrum
					curScore,            // Score of the current path through the graph
					minPathScore=0,      // Minimum acceptable path score (updated as higher scoring paths are found)
					maxPathScore=0;      // Maximum path score found

		unsigned int numPeaks = spec.size();
		unsigned int idxPref=0, idxSuff=numPeaks-1;   // Spectrum indices of the current path endpoints

		// Set firstSuff to index of the first peak after the middle of the spectrum
		unsigned int firstSuff = min(idxPref,idxSuff); while(firstSuff<numPeaks and spec[firstSuff][0]<middle) firstSuff++;
		unsigned int suffCount = numPeaks-firstSuff;     // Number of suffix peaks - always >0 due to addZPMpeaks above

		// Declare and initialize the DP structures - 2D arrays are filled from top-right to left-bottom (increasing prefix, decreasing suffix)
		vector<vector<list<TwoValues<unsigned int> > > > predPairs;  // list of predecessor pairs for every prefix/suffix pair - pairs are in _spectrum_ indices, not matrix indices
		vector<vector<float> > maxScores;                // pos [i,j] = maximum de-novo sequence score with a prefix using peak i and a suffix using peak j
		list<TwoValues<unsigned int> > pathStarts;       // List of pairs[i,j] that connect prefixes and suffixes with a total path score larger than minPathScore

		predPairs.resize(firstSuff);   maxScores.resize(firstSuff);   pathStarts.clear();
		for(unsigned int i=0; i<firstSuff; i++) {
			predPairs[i].resize(suffCount);   maxScores[i].resize(suffCount);
			for(unsigned int j=0; j<suffCount; j++) { predPairs[i][j].clear(); maxScores[i][j]=-1; }
		}

		unsigned int nextPref, nextSuff, jumpType;    idxPref=0; idxSuff=numPeaks-1;
		unsigned int neighIdx;
		bool prefixExtension;
		// Initialize scores and predecessors of peak pairs reachable from endpoints
		for(idxPref=0; idxPref<firstSuff; idxPref++)
			for(idxSuff=numPeaks-1; idxSuff>=firstSuff; idxSuff--) {
				unsigned int vNextSuff, vSuff=idxSuff-firstSuff;  // corrected indices to access DP structures

				// Any peaks within tolerance of b0/bn are eligible starters
				if(fabs(spec[idxPref][0]-ionOffset)<=peakTol and fabs(spec[idxSuff][0]-massBn)<=max(peakTol,pmTol))
					maxScores[idxPref][vSuff] = spec[idxPref][1]+spec[idxSuff][1];

				if(maxScores[idxPref][vSuff]<0) continue; // Skip unreachable pairs of peaks

				// Try both prefix and suffix extensions - must jump _over_ symmetric suffix/prefix, respectively
				for(jumpType=0; jumpType<2; jumpType++) {
					prefixExtension = (bool) jumpType;    neighIdx=0;

					// Initialize depending on whether this is a prefix or a suffix extension - jumping forward towards the middle and over the complement of the pair
					if(prefixExtension) {
						nextSuff=idxSuff;
						if(neighIdx<neighsR[idxPref].size()) {
							nextPref=neighsR[idxPref][neighIdx];
	// Anti-symmetric allowed
							while(neighIdx<neighsR[idxPref].size() and nextPref<firstSuff and spec[nextPref][0]<=aaMass-spec[nextSuff][0]-2*peakTol)
	// Anti-symmetric not allowed
	//						while(neighIdx<neighsR[idxPref].size() and nextPref<firstSuff and spec[nextPref][0]<=aaMass-spec[nextSuff][0]+2*peakTol) // Anti-symmetric not allowed
								if(++neighIdx<neighsR[idxPref].size()) nextPref=neighsR[idxPref][neighIdx];
							vNextSuff = vSuff;
						}
					} else {
						nextPref=idxPref;
						if(neighIdx<neighsL[idxSuff].size()) {
							nextSuff=neighsL[idxSuff][neighIdx];
	// Anti-symmetric allowed
							while(neighIdx<neighsL[idxSuff].size() and nextSuff>=firstSuff and spec[nextSuff][0]>=aaMass-spec[nextPref][0]+2*peakTol)
	// Anti-symmetric not allowed
	//						while(neighIdx<neighsL[idxSuff].size() and nextSuff>=firstSuff and spec[nextSuff][0]>=aaMass-spec[nextPref][0]-2*peakTol) // Anti-symmetric not allowed
								if(++neighIdx<neighsL[idxSuff].size()) nextSuff=neighsL[idxSuff][neighIdx];
						}
						if (nextSuff>firstSuff) vNextSuff = nextSuff - firstSuff; else vNextSuff = 0;
					}

					// iterate over all reachable peaks
					for(; (prefixExtension and neighIdx<neighsR[idxPref].size()) or (not prefixExtension and neighIdx<neighsL[idxSuff].size()); neighIdx++) {
						if(prefixExtension) {  // Find nextPrev by valid jumps from idxPrev
							nextPref=neighsR[idxPref][neighIdx];   if(nextPref>=firstSuff and nextPref!=idxSuff) continue;
						} else {
							nextSuff=neighsL[idxSuff][neighIdx];   if(nextSuff<firstSuff and nextSuff!=idxPref) continue;
							if (nextSuff>firstSuff) vNextSuff = nextSuff - firstSuff; else vNextSuff = 0;
						}

						if(nextPref==idxSuff or nextSuff==idxPref) { // connecting jump
							if(enforcePM or maxScores[idxPref][vSuff]>=minPathScore) {  // if enforcePM then all possible connecting sequences have to be tried to check whether the parent mass matches
								pathStarts.push_front(TwoValues<unsigned int>(idxPref,idxSuff));
								if(maxScores[idxPref][vSuff]>maxPathScore) {
									maxPathScore = maxScores[idxPref][vSuff];
									minPathScore = minScorePerc*maxPathScore;
	//cerr<<"denovo: updated minPathScore to "<<minPathScore<<" from maxPathScore "<<maxPathScore<<" at pair ("<<idxPref<<","<<idxSuff<<")\n";
								}
							}
							break;
						}

						// Add [idxPrev,idxSuff] to list of predecessors/predScores
						predPairs[nextPref][vNextSuff].push_front(TwoValues<unsigned int>(idxPref,idxSuff));

						// Adjust maximum score of destination pair
						if(fabs(aaMass-spec[nextPref][0]-spec[nextSuff][0])>2*peakTol)
							curScore=maxScores[idxPref][vSuff] + (prefixExtension?spec[nextPref][1]:spec[nextSuff][1]);
	//					else curScore=maxScores[idxPref][vSuff];   // Don't double-count
	//					else curScore=maxScores[idxPref][vSuff] - (prefixExtension?spec[nextSuff][1]:spec[nextPref][1]) + (spec[nextSuff][1]+spec[nextPref][1])/2;   // Using symmetric peaks adds the average score of the peak pair. (same as not double-counting for symmetric spectra)
						else curScore=maxScores[idxPref][vSuff] - (prefixExtension?spec[nextSuff][1]:spec[nextPref][1]) + max(spec[nextSuff][1],spec[nextPref][1]);   // Using symmetric peaks adds the max score of the peak pair. (same as not double-counting for symmetric spectra)
	//if(nextPref==6 and nextSuff==34) cerr<<">>> curScore = "<<curScore<<", maxScore = "<<maxScores[nextPref][vNextSuff]<<"\n";
						if(maxScores[nextPref][vNextSuff]<curScore) {
							maxScores[nextPref][vNextSuff] = curScore;
						}
					} // over all reachable peaks
				} // over prefix/suffix extensions
			} // over all possible suffixes

	//for(unsigned int peakIdx=0; peakIdx<spec.size(); peakIdx++)
	//	cerr<<peakIdx<<": "<<spec[peakIdx][0]<<"\t"<<spec[peakIdx][1]<<"\n";

		list<TwoValues<unsigned int> >::iterator iter, tmpIter;
		for(idxPref=0; idxPref<firstSuff; idxPref++)
			for(idxSuff=numPeaks-1; idxSuff>=firstSuff; idxSuff--) {
				unsigned int vSuff=idxSuff-firstSuff;  // corrected indices to access DP structures
				predPairs[idxPref][vSuff].sort();   iter=predPairs[idxPref][vSuff].begin();
				while(iter!=predPairs[idxPref][vSuff].end()) {
					tmpIter=iter; tmpIter++;
					if(tmpIter!=predPairs[idxPref][vSuff].end() and *iter==*tmpIter) {
						iter=predPairs[idxPref][vSuff].erase(iter);
					} else iter++;
				}
			}

		// Remove repeated pathStarts (if the connecting pair was above minPathScore more than once)
		pathStarts.sort();
		iter = pathStarts.begin(), tmpIter;
		while(iter!=pathStarts.end()) { tmpIter=iter; tmpIter++; if(tmpIter!=pathStarts.end() and *iter==*tmpIter) iter=pathStarts.erase(iter); else iter++; }

	//cerr<<"Final minPathScore = "<<minPathScore<<", connecting jumps/scores:\n";
	//for(iter=pathStarts.begin();iter!=pathStarts.end();iter++) cerr<<"("<<(*iter)[0]<<","<<(*iter)[1]<<"): "<<maxScores[(*iter)[0]][(*iter)[1]-firstSuff]<<"\n";

	/*cerr<<"denovo: spectrum has "<<spec.size()<<" peaks, middle = "<<middle<<", firstSuff = "<<firstSuff<<"\n";
	vector<TwoValues<unsigned int> > dbg(9);
	dbg[0].set(2,61); dbg[1].set(2,56); dbg[2].set(7,56);
	dbg[3].set(7,52); dbg[4].set(7,41); dbg[5].set(22,41);
	dbg[6].set(25,41); dbg[7].set(22,36); dbg[8].set(25,36);
	for(unsigned int i=0;i<dbg.size();i++) {
		if(dbg[i][0]<spec.size() and dbg[i][1]<spec.size()) {
			cerr<<"Pair ("<<dbg[i][0]<<","<<dbg[i][1]<<"), masses("<<spec[dbg[i][0]][0]<<","<<spec[dbg[i][1]][0]<<"), maxScores="<<maxScores[dbg[i][0]][dbg[i][1]-firstSuff]<<": ";
			for(list<TwoValues<unsigned int> >::iterator iter=predPairs[dbg[i][0]][dbg[i][1]-firstSuff].begin();iter!=predPairs[dbg[i][0]][dbg[i][1]-firstSuff].end();iter++)
				cerr<<"("<<(*iter)[0]<<","<<(*iter)[1]<<")";
			cerr<<endl;
		}
	}*/

	//	list<list<TwoValues<float> > > curSeqs;   list<float> curScores;   Spectrum cSpec;    cSpec=spec;
	//	list<list<TwoValues<float> > > allSeqs;   list<float> allScores;
		list<list<unsigned int> > curSeqs;   list<float> curScores;   Spectrum cSpec;    cSpec=spec;
		list<list<unsigned int> > allSeqs;   list<float> allScores;
		if(enforcePM and minScorePerc>=1-0.0001) minPathScore=0; else minPathScore=minPathScore*.9999;   // minPathScore*.9999 because of numeric rounding errors
		module_num_rejected_seqs = 0;
		allSeqs.clear();   allScores.clear();
		for(iter=pathStarts.begin(); iter!=pathStarts.end(); iter++) {
			if(enforcePM) {
				cerr<<"ERROR: version of denovo() using neighsL/neighsR does not yet support peptide sequencing with exact parent masses.\n"; exit(-1);
	/*			TwoValues<unsigned int> jumpBounds;
				jumps.find(spec[(*iter)[1]][0]-spec[(*iter)[0]][0],peakTol,jumpBounds);
				for(unsigned int jumpIdx=jumpBounds[0]; jumpIdx<jumpBounds[1]; jumpIdx++) {
					curSeqs.clear();   curScores.clear();
					minPathScore=denovo_buildSeqs_enforcePM(spec,predPairs,maxScores,*iter,-1,firstSuff,minPathScore,peakTol,pmOffset,pmTol,jumps[jumpIdx]-(spec[(*iter)[1]][0]-spec[(*iter)[0]][0]),jumps,minScorePerc>=1-0.0001,curSeqs,curScores);
					allSeqs.splice(allSeqs.end(), curSeqs);    allScores.splice(allScores.end(), curScores);
				}
	*/		} else {
				if(maxScores[(*iter)[0]][(*iter)[1]-firstSuff]<minPathScore) continue;  // Remove pathStarts that fell below the minPathScore threshold (may happen when min/maxPathScore are updated)
				curSeqs.clear();   curScores.clear();
	//if((*iter)[0]!=14 or (*iter)[1]!=24) continue;
				denovo_buildSeqs_idx(spec,predPairs,maxScores,*iter,-1,firstSuff,minPathScore,peakTol,pmOffset,curSeqs,curScores);
	/*
	list<float>::iterator scoresIter=curScores.begin();
	for(list<list<unsigned int> >::iterator seqIter=curSeqs.begin(); seqIter!=curSeqs.end(); seqIter++, scoresIter++) {
		for(list<unsigned int>::iterator peakIter=seqIter->begin(); peakIter!=seqIter->end(); peakIter++)
			cerr<<*peakIter<<",";
		cerr<<" score: "<<*scoresIter<<"\n";
	}*/
				allSeqs.splice(allSeqs.end(), curSeqs);    allScores.splice(allScores.end(), curScores);
			}
		}

		// Returns only those sequences observing the final minPathScore threshold
		list<float>::iterator scoresIter = allScores.begin();
	//	for(list<list<TwoValues<float> > >::iterator seqsIter=allSeqs.begin(); seqsIter!=allSeqs.end(); seqsIter++, scoresIter++) {
		for(list<list<unsigned int> >::iterator seqsIter=allSeqs.begin(); seqsIter!=allSeqs.end(); seqsIter++, scoresIter++) {
			if(*scoresIter>=minPathScore-0.0001) {
				cSpec.resize(seqsIter->size());  int peakIdx=0;
	//			for(list<TwoValues<float> >::iterator peakIter=seqsIter->begin(); peakIter!=seqsIter->end(); peakIter++)
	//				cSpec[peakIdx++]=*peakIter;
				for(list<unsigned int>::iterator peakIter=seqsIter->begin(); peakIter!=seqsIter->end(); peakIter++)
					cSpec[peakIdx++]=spec[*peakIter];
				seqs.push_back(cSpec);    seqScores.push_back(*scoresIter);
			}
		}
	}


	void getNeighborsL(Spectrum &spec, AAJumps &jumps, float peakTol, vector<vector<int> > &neighs) {
		int peakIdx,neighIdx,neighsCount,jumpIdx;
		float curOffset;
		neighs.resize(spec.size());

		for(peakIdx=0; peakIdx<(int)spec.size(); peakIdx++) {
			neighs[peakIdx].resize(peakIdx);   neighsCount=0;   jumpIdx=0;
			for(neighIdx=peakIdx-1; neighIdx>=0; neighIdx--) {
				curOffset = spec[peakIdx][0]-spec[neighIdx][0];
				while(jumpIdx<(int)jumps.size() and jumps[jumpIdx]<curOffset-peakTol) jumpIdx++;
				if(jumpIdx>=(int)jumps.size()) break;
				if(fabs(curOffset-jumps[jumpIdx])<=peakTol) neighs[peakIdx][neighsCount++]=neighIdx;
			}
			neighs[peakIdx].resize(neighsCount);
		}
	}

	//
	//  getNeighborsL - for each peak p in the spectrum, finds the rightmost peak that is
	//    at least minPeakDist-peakTol Da to the left of p. neighs[p].size() is set to zero
	//    when no such peak exists.
	//
	void getNeighborsL(Spectrum &spec, float minPeakDist, float peakTol, vector<vector<int> > &neighs) {
		int peakIdx,neighIdx;
		neighs.resize(spec.size());
		for(peakIdx=0; peakIdx<(int)spec.size(); peakIdx++) {
			neighs[peakIdx].resize(1);   neighs[peakIdx][0]=-1;
			for(neighIdx=peakIdx-1; neighIdx>=0; neighIdx--)
				if(spec[peakIdx][0]-spec[neighIdx][0]>=minPeakDist-peakTol)
					{ neighs[peakIdx][0]=neighIdx; break; }
			if(neighs[peakIdx][0]==-1) neighs[peakIdx].resize(0);
		}
	}

	//
	//  getNeighborsL - for each peak p in the spectrum, finds all peaks that are
	//    at more than minPeakDist-peakTol Da and less than maxPeakDist+peakTol Da
	//    to the left of p. neighs[p].size() is set to zero when no such peaks exist.
	//
	void getNeighborsL(Spectrum &spec, float minPeakDist, float maxPeakDist, float peakTol, vector<vector<int> > &neighs) {
		int peakIdx,neighIdx,neighsCount;
		neighs.resize(spec.size());
		for(peakIdx=0; peakIdx<(int)spec.size(); peakIdx++) {
			neighs[peakIdx].resize(peakIdx);   neighsCount=0;
			for(neighIdx=peakIdx-1; neighIdx>=0; neighIdx--)
				if(spec[peakIdx][0]-spec[neighIdx][0]>=minPeakDist-peakTol) {
					if(spec[peakIdx][0]-spec[neighIdx][0]<=maxPeakDist+peakTol) neighs[peakIdx][neighsCount++]=neighIdx;
					else break;
				}
			neighs[peakIdx].resize(neighsCount);
		}
	}

	void getNeighborsR(Spectrum &spec, AAJumps &jumps, float peakTol, vector<vector<int> > &neighs) {
		int peakIdx,neighIdx,neighsCount,jumpIdx;
		float curOffset;
		neighs.resize(spec.size());
		for(peakIdx=0; peakIdx<(int)spec.size(); peakIdx++) {
			neighs[peakIdx].resize(spec.size()-peakIdx);   neighsCount=0;   jumpIdx=0;
			for(neighIdx=peakIdx+1; neighIdx<(int)spec.size(); neighIdx++) {
				curOffset = spec[neighIdx][0]-spec[peakIdx][0];
				while(jumpIdx<(int)jumps.size() and jumps[jumpIdx]<curOffset-peakTol) jumpIdx++;
				if(jumpIdx>=(int)jumps.size()) break;
				if(fabs(curOffset-jumps[jumpIdx])<=peakTol) neighs[peakIdx][neighsCount++]=neighIdx;
			}
			neighs[peakIdx].resize(neighsCount);
		}
	}

	//
	//  getNeighborsR - for each peak p in the spectrum, finds all peaks that are
	//    at more than minPeakDist-peakTol Da and less than maxPeakDist+peakTol Da
	//    to the right of p. neighs[p].size() is set to zero when no such peaks exist.
	//
	void getNeighborsR(Spectrum &spec, float minPeakDist, float maxPeakDist, float peakTol, vector<vector<int> > &neighs) {
		int peakIdx,neighIdx,neighsCount;
		neighs.resize(spec.size());
		for(peakIdx=0; peakIdx<(int)spec.size(); peakIdx++) {
			neighs[peakIdx].resize(spec.size()-peakIdx);   neighsCount=0;
			for(neighIdx=peakIdx+1; neighIdx<=(int)spec.size(); neighIdx++)
				if(spec[neighIdx][0]-spec[peakIdx][0]>=minPeakDist-peakTol) {
					if(spec[neighIdx][0]-spec[peakIdx][0]<=maxPeakDist+peakTol) neighs[peakIdx][neighsCount++]=neighIdx;
					else break;
				}
			neighs[peakIdx].resize(neighsCount);
		}
	}

	//
	//  aux_denovo_LtoR - auxiliary function to recursively output all suboptimal sequences
	//   determined by denovo_LtoR. Most parameters are exactly as denovo_LtoR, except for
	//
	//  predCumScore, cumScore - see denovo_LtoR()
	//  curPeakIdx    - Index of the peak currently being added to curPath
	//  curPath       - Current path (not including curPeakIdx)
	//  curPathScore  - Cumulative score of the current path (not including curPeakIdx)
	//
	void aux_denovo_LtoR(Spectrum &spec, vector<vector<int> > &neighs, vector<vector<float> > &predCumScore,
											 vector<float> &cumScore, list<list<int> > &seqs, unsigned int curPeakIdx,
											 list<int> curPath, float cumPathScore, float minAcceptableScore) {
		curPath.push_front(curPeakIdx);    cumPathScore+=spec[curPeakIdx][1];
		if(cumPathScore>=minAcceptableScore) seqs.push_back(curPath);
		for(unsigned int neighIdx=0; neighIdx<neighs[curPeakIdx].size(); neighIdx++)
			if(predCumScore[curPeakIdx][neighIdx]+cumPathScore-spec[curPeakIdx][1] >= minAcceptableScore)
				aux_denovo_LtoR(spec,neighs,predCumScore,cumScore,seqs, neighs[curPeakIdx][neighIdx], curPath, cumPathScore, minAcceptableScore);
	}

	//
	//  denovo_LtoR - Finds the best path from the left of the spectrum to the right.
	//    No peak symmetries (b/y) are considered. Each peak's left neighbors are given
	//    in neighs, as determined by one of the getNeighborsL functions above.
	//
	//float denovo_LtoR(Spectrum &spec, vector<vector<int> > &neighs, list<list<int> > &seqs, float minPercMaxScore) {
	float denovo_LtoR(Spectrum &spec, vector<vector<int> > &neighs, list<Spectrum> &seqs, list<float> &seqScores, float minPercMaxScore) {
		vector<float> cumScore(spec.size());            // Accumulated score at each peak
	//	vector<TwoValues<int> > predInfo(spec.size());  // Per peak: predecessor index (pos.0), num predecessors (pos.1)
		vector<vector<float> > predCumScore(spec.size());  // Per peak: best cumulative score per predecessor: predCumScore[i][j] for neighs[i][j]
		unsigned int peakIdx, neighIdx, bestPathEnd=0;
		float bestPathScore=0;

		seqs.clear(); seqScores.clear(); if(spec.size()==0) return 0;

		predCumScore[0].resize(0);   cumScore[0]=spec[0][1];
		for(peakIdx=1; peakIdx<spec.size(); peakIdx++) {
			predCumScore[peakIdx].resize(neighs[peakIdx].size());   cumScore[peakIdx]=spec[peakIdx][1];
			for(neighIdx=0; neighIdx<neighs[peakIdx].size(); neighIdx++) {
				predCumScore[peakIdx][neighIdx] = spec[peakIdx][1]+cumScore[neighs[peakIdx][neighIdx]];
				if(predCumScore[peakIdx][neighIdx] > cumScore[peakIdx]) {
					cumScore[peakIdx] = predCumScore[peakIdx][neighIdx];
	//cerr<<"denovoLR: peak ("<<peakIdx<<","<<spec[peakIdx][0]<<") changed best predecessor to ("<<neighIdx<<","<<neighs[peakIdx][neighIdx]<<","<<spec[neighs[peakIdx][neighIdx]][0]<<"), cumScore"<<cumScore[peakIdx]<<"\n";
	//				predInfo[peakIdx][0] = neighs[peakIdx][neighIdx];
	//				predInfo[peakIdx][1] = predInfo[neighs[peakIdx][neighIdx]][1]+1;
				}
			}
			if(bestPathScore<cumScore[peakIdx])
				{ bestPathEnd=peakIdx; bestPathScore=cumScore[peakIdx]; }
		}

		float minAcceptableScore = minPercMaxScore * bestPathScore;    list<int> curPath;
		list<list<int> > indexSeqs;
		for(peakIdx=spec.size()-1; peakIdx>0; peakIdx--)
			if(cumScore[peakIdx]>=minAcceptableScore)
				aux_denovo_LtoR(spec,neighs,predCumScore,cumScore,indexSeqs, peakIdx, curPath, 0, minAcceptableScore);

		// Returns only those sequences observing the final minPathScore threshold
		Spectrum cSpec;    cSpec.copyNP(spec);   float score;
		for(list<list<int> >::iterator seqsIter=indexSeqs.begin(); seqsIter!=indexSeqs.end(); seqsIter++) {
			cSpec.resize(seqsIter->size());    peakIdx=0;    score=0;
			for(list<int>::iterator peakIter=seqsIter->begin(); peakIter!=seqsIter->end(); peakIter++)
				{ cSpec[peakIdx++]=spec[*peakIter]; score+=spec[*peakIter][1]; }
			seqs.push_back(cSpec);    seqScores.push_back(score);
		}

		return bestPathScore;
	}

	float denovo_LtoR(Spectrum &spec, vector<vector<int> > &neighs, vector<int> &idxMatched) {
		vector<float> cumScore(spec.size());            // Accumulated score at each peak
		vector<TwoValues<int> > predInfo(spec.size());  // Per peak: predecessor index (pos.0), num predecessors (pos.1)
		unsigned int peakIdx, neighIdx, bestPathEnd=0;
		float bestPathScore=0;

		if(spec.size()==0) { idxMatched.resize(0); return 0; }

		predInfo[0][0]=-1;   predInfo[0][1]=0;   cumScore[0]=spec[0][1];
		for(peakIdx=1; peakIdx<spec.size(); peakIdx++) {
			predInfo[peakIdx][0]=-1;   predInfo[peakIdx][1]=0;   cumScore[peakIdx]=spec[peakIdx][1];
			for(neighIdx=0; neighIdx<neighs[peakIdx].size(); neighIdx++) {
				if(spec[peakIdx][1]+cumScore[neighs[peakIdx][neighIdx]]>cumScore[peakIdx]) {
					cumScore[peakIdx] = spec[peakIdx][1]+cumScore[neighs[peakIdx][neighIdx]];
	//cerr<<"denovoLR: peak ("<<peakIdx<<","<<spec[peakIdx][0]<<") changed best predecessor to ("<<neighIdx<<","<<neighs[peakIdx][neighIdx]<<","<<spec[neighs[peakIdx][neighIdx]][0]<<"), cumScore"<<cumScore[peakIdx]<<"\n";
					predInfo[peakIdx][0] = neighs[peakIdx][neighIdx];
					predInfo[peakIdx][1] = predInfo[neighs[peakIdx][neighIdx]][1]+1;
				}
			}
			if(bestPathScore<cumScore[peakIdx])
				{ bestPathEnd=peakIdx; bestPathScore=cumScore[peakIdx]; }
		}

		idxMatched.resize(predInfo[bestPathEnd][1]+1);   int idx=predInfo[bestPathEnd][1];
		while(predInfo[bestPathEnd][0]>=0)
			{ idxMatched[idx--]=bestPathEnd; bestPathEnd=predInfo[bestPathEnd][0]; }

		return bestPathScore;
	}

	//
	// Code to estimate the distribution of peptide scores over a given spectrum
	//
	// Note: scoresDist could be implemented as a sliding window in case there is
	//  a need to save memory for higher resolutions (only need to look back up to
	//  mass(W)/resolution bins).
	//
	// Note: Current implementation enforces a parent mass tolerance of resolution/2.
	//
	void EstimateScoresDistribution(Spectrum &spec, float peakTol, float pmTol, float resolution,
										vector<unsigned int> &hist, bool addIntensities) {
		unsigned int binRadius = (unsigned int)round(peakTol/resolution);
		const unsigned int infinity =  (unsigned int)(ULONG_MAX>>2); // Maximum number of peptides for any given score (to avoid overflowing)
		const unsigned int numScoreBins = 1+1000;   // 1000 = percentage of max score in 0.1 increments
		AAJumps jumps(1,resolution);                // Amino acid masses used to construct all possible peptides
		unsigned int numMassBins = 1+(unsigned int)round((spec.parentMass-AAJumps::massMH+pmTol)/resolution),
								 idxPeak, idxBin, idxScore, peakBin, startBin, endBin;

		hist.resize(numScoreBins);   for(idxScore=0; idxScore<numScoreBins; idxScore++) hist[idxScore]=0;
		if(spec.parentMass<=AAJumps::massMH) return;
		vector<unsigned int> binScores(numMassBins);

		// Remove negative-score peaks
		float minScore=0;
		for(idxPeak=0; idxPeak<spec.size(); idxPeak++) if(spec[idxPeak][1]<minScore) minScore=spec[idxPeak][1];
		if(minScore<0) for(idxPeak=0; idxPeak<spec.size(); idxPeak++) spec[idxPeak][1]-=minScore;

		// Remove peaks that are not reachable from the start/end of the spectrum
	/*	unsigned int idxNew=0;
		for(idxPeak=0; idxPeak<spec.size(); idxPeak++)
			if(spec[idxPeak][0]>=jumps[0]-peakTol and spec[idxPeak][0]<=spec.parentMass-AAJumps::massMH-jumps[0]+peakTol)
				spec[idxNew++]=spec[idxPeak];
		spec.resize(idxNew);
	*/

		// Normalize peak intensities to sum up to numScoreBins
		float totalScore=0;
		for(idxPeak=0; idxPeak<spec.size(); idxPeak++) totalScore+=spec[idxPeak][1];
		float factor = ((float)numScoreBins)/totalScore;
		for(idxPeak=0; idxPeak<spec.size(); idxPeak++) spec[idxPeak][1]*=factor;

		// Add intensities of peaks within tolerance of each bin
		for(idxBin=0; idxBin<numMassBins; idxBin++) binScores[idxBin]=0;
		for(idxPeak=0; idxPeak<spec.size(); idxPeak++) {
			peakBin = (unsigned int)round(spec[idxPeak][0]/resolution);
			startBin = peakBin>binRadius ? peakBin-binRadius : 0 ;
			endBin = peakBin+binRadius>=numMassBins ? numMassBins : peakBin+binRadius+1 ;
			for(idxBin=startBin; idxBin<endBin; idxBin++)
				if(addIntensities) binScores[idxBin]+=(unsigned int)spec[idxPeak][1];
				else binScores[idxBin]=max(binScores[idxBin],(unsigned int)spec[idxPeak][1]);
		}

		// Create and initialize minmaxScore
		vector<TwoValues<unsigned int> > minmaxScore(numMassBins);
		vector<bool> validPepMass(numMassBins);
		for(idxBin=0; idxBin<numMassBins; idxBin++) { minmaxScore[idxBin].set(numScoreBins,0); validPepMass[idxBin]=false; }

		// Create and initialize scoresDist
		// scoresDist[i][j] - number of peptides of parent mass i*resolution with a
		//   score total_spectrum_intensity * j/numScoreBins
		vector<vector<unsigned int> > scoresDist(numMassBins);
		for(idxBin=0; idxBin<numMassBins; idxBin++) {
			scoresDist[idxBin].resize(numScoreBins);
			for(idxScore=0; idxScore<numScoreBins; idxScore++) scoresDist[idxBin][idxScore]=0;
		}
	//	endBin = (unsigned int)round(peakTol/resolution)+1;   endBin=min(endBin,numMassBins);  // Allow mass tolerance at 0
		endBin = 1;  // Allow only exact mass of zero
		for(idxBin=0; idxBin<endBin; idxBin++) {
			scoresDist[idxBin][binScores[idxBin]]=1;  // Only the empty peptide string starts and ends at mass 0
			minmaxScore[idxBin].set(binScores[idxBin],binScores[idxBin]);
			validPepMass[idxBin]=true;
		}

		// Preprocess possible jumps to integer #bins offsets
		unsigned int idxJump, jumpBin;
		vector<unsigned int> jumpBinsOffset(jumps.size());
		for(idxJump=0; idxJump<jumps.size(); idxJump++)
			jumpBinsOffset[idxJump]=(unsigned int)round(jumps[idxJump]/resolution);

		// Compute distribution of scores
		unsigned int scoreOffset;
		for(idxBin=0; idxBin<numMassBins; idxBin++) {
			scoreOffset = binScores[idxBin];
			for(idxJump=0; idxJump<jumpBinsOffset.size(); idxJump++) {
				// Allow mass tolerance for peptide masses
	//			if(jumpBinsOffset[idxJump]>idxBin+binRadius) continue;  // Can't jump to masses <0
	//			jumpBin = (jumpBinsOffset[idxJump]<idxBin) ? idxBin-jumpBinsOffset[idxJump] : 0;
				// Do _not_ allow mass tolerance for peptide masses, only for observed masses (via tolerance in binScores)
				if(jumpBinsOffset[idxJump]>idxBin) continue; else jumpBin=idxBin-jumpBinsOffset[idxJump];  // Can't jump to masses <0
				if(not validPepMass[jumpBin]) continue; // jumpBin was not reachable by a valid peptide

				validPepMass[idxBin]=true;
				startBin = minmaxScore[jumpBin][0] + scoreOffset;    if(minmaxScore[idxBin][0]>startBin) minmaxScore[idxBin][0]=startBin;
				endBin = min(numScoreBins-1,minmaxScore[jumpBin][1]+scoreOffset);   if(minmaxScore[idxBin][1]<endBin) minmaxScore[idxBin][1]=endBin;
				for(idxScore=startBin; idxScore<=endBin; idxScore++)  { // Collect #peptides for each score
					scoresDist[idxBin][idxScore] += scoresDist[jumpBin][idxScore-scoreOffset];
					if(scoresDist[idxBin][idxScore]>infinity) scoresDist[idxBin][idxScore]=infinity;
				}
			}
		}

		// Copy final histogram to hist
	//	for(idxScore=0; idxScore<numScoreBins; idxScore++) hist[idxScore]=scoresDist[numMassBins-1][idxScore];  // When there is tolerance at mass 0
		for(idxBin=numMassBins-min(numMassBins,(unsigned int)round(pmTol/resolution)-1); idxBin<numMassBins; idxBin++)
			for(idxScore=0; idxScore<numScoreBins; idxScore++) hist[idxScore]+=scoresDist[idxBin][idxScore];
	}

	//
	// peaksToMassBins - Each peak in spec creates (1+peakTol/resolution) bins. Each bin's
	//   intensity is given by the summed intensity of all peaks within tolerance.
	//
	void peaksToMassBins(Spectrum &spec, Spectrum &bins, float peakTol, float resolution) {
		unsigned int binRadius = (unsigned int)round(peakTol/resolution);
		unsigned int maxBin = (unsigned int)round(spec[spec.size()-1][0]/resolution)+binRadius,
								 idxPeak, idxBin, peakBin, startBin, endBin;
		vector<float> binScores(maxBin+1);

		// Add intensities of peaks within tolerance of each bin
		for(idxBin=0; idxBin<maxBin; idxBin++) binScores[idxBin]=0;
		for(idxPeak=0; idxPeak<spec.size(); idxPeak++) {
			peakBin = (unsigned int)round(spec[idxPeak][0]/resolution);
			startBin = peakBin>binRadius ? peakBin-binRadius : 0 ;
			endBin = peakBin+binRadius>maxBin ? maxBin : peakBin+binRadius ;
			for(idxBin=startBin; idxBin<=endBin; idxBin++) binScores[idxBin]+=spec[idxPeak][0];
		}

		// Convert binScores vector to bins' Spectrum format (i.e. list of non-zero intensity bins)
		unsigned int nonZeroBins = 0, outputBin;
		for(idxBin=0; idxBin<maxBin; idxBin++) nonZeroBins+=(binScores[idxBin]==0);
		bins.resize(nonZeroBins);
		for(outputBin=0, idxBin=0; idxBin<maxBin; idxBin++)
			if(binScores[idxBin]>0)
				bins[outputBin++].set(((float)idxBin)*resolution,binScores[idxBin]);
	}
}
