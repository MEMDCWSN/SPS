#include "inputParams.h"
#include "batch.h"
#include "graph.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <algorithm>

using namespace std;

class MassIdxPair {
	public:
	float mass; int idx;
};
bool operator<(const MassIdxPair &a, const MassIdxPair &b) { return ((a.mass<b.mass) or (a.mass==b.mass and a.idx<b.idx)); }

/*
  findconnected - Identifies sets of spectra connected by spectral alignments in
   structures relevant to identification of multiple modification states of the
   same peptide.

  INPUT:
    INPUT_SPECS   - Input PRM-spectra (or star spectra, used only for parent masses)
	INPUT_ALIGNS  - Set of detected spectral alignment pairs
    PATH_LENGTH   - Minimum path length between minimally and maximally modified variants
					  >0  - look for highly modified peptides
                      <=0 - look for 0/+1/+1/+2 modification quadruplets:
                              pairs of spectra A,B such that: C=succ(A)=pred(B),
                              D=succ(A)=pred(B) and PM(A)+PM(B)=PM(C)+PM(D)
	MIN_NUM_PATHS - Minimum number of paths between minimally and maximally modified variants
    TOLERANCE_PM  - Parent mass tolerance (in Daltons, default 1 Da)

  OUTPUT:
    OUTPUT_SETS   - Sets of identified connected components, binListArray format, zero-based indices:
					  PATH_LENGTH >0  - Indices of spectra participating in at least one path between
					    minimally and maximally modified variants

                      PATH_LENGTH <=0 - Every pair A,B generates one entry in vSets where
                        the first two elements are A,B and the remaining elements are all
                        possible C,D vertices.
    OUTPUT_NUMPATHS- Number of found paths per component (binArray format, only output if PATH_LENGTH>0)
*/

int main(int argc, char **argv){

    // Get input parameters and check minimum parameter set
    InputParams params; bool paramsOk;
	if(argc<=1) paramsOk=params.readParams("findconnected.params");	else paramsOk=params.readParams(argv[1]);
	if(!paramsOk) {
		cerr << "Error opening parameters file ";
		if(argc<=1) cerr << "findconnected.params\n"; else cerr << argv[1] << endl;
		return -1;
	}
	vector<char *> paramStrs;   paramStrs.resize(5);
	paramStrs[0] = "INPUT_SPECS";
	paramStrs[1] = "INPUT_ALIGNS";
	paramStrs[2] = "OUTPUT_SETS";
	paramStrs[3] = "PATH_LENGTH";
	paramStrs[4] = "MIN_NUM_PATHS";
	if(!params.confirmParams(paramStrs)) {
		cerr << "ERROR: Parameters file ";
		if(argc==1) cerr<<"findconnected.params"; else cerr<<argv[1];
		cerr << " is incomplete. One of the following is missing: INPUT_SPECS, INPUT_ALIGNS, OUTPUT_SETS, PATH_LENGTH, MIN_NUM_PATHS\n";
		return -1;
	}

	char *specSetFN = params.getValue("INPUT_SPECS");
	char *alignsFN = params.getValue("INPUT_ALIGNS");
	char *resultsFN = params.getValue("OUTPUT_SETS");
	int pathLength = params.getValueInt("PATH_LENGTH");
	int minNumPaths = params.getValueInt("MIN_NUM_PATHS");

	int startBaseIdx = params.paramPresent("IDX_START")?params.getValueInt("IDX_START"):0;
	int endBaseIdx = params.paramPresent("IDX_END")?params.getValueInt("IDX_END"):0;
	float peakTol = params.paramPresent("TOLERANCE_PEAK")?(float) params.getValueDouble("TOLERANCE_PEAK"):0.5;
	float pmTol = params.paramPresent("TOLERANCE_PM")?(float) params.getValueDouble("TOLERANCE_PM"):1;

    SpecSet specSet;
    if (specSet.LoadSpecSet_pklbin(specSetFN)<=0) { cerr << "Error reading "<<specSetFN<<"!\n"; return -1; }
    else  cout << "Loading specs complete. Num specs: " << specSet.size() << "\n";

	vector<Results_ASP> aligns;
	if (Load_results(alignsFN, aligns)<=0) { cerr << "Error reading "<<alignsFN<<"!\n"; return -1; }
    else  cout << "Loading aligns complete. Num pairs: " << aligns.size() << "\n";
    if(endBaseIdx==0) endBaseIdx = aligns.size()-1;

	//
	// Sort spectra by [parent mass, original index] + change pairs in aligns correspondigly
	//
	vector<MassIdxPair> massIdxPairs(specSet.size());
	for(unsigned int i=0; i<specSet.size(); i++) { massIdxPairs[i].mass=specSet[i].parentMass;  massIdxPairs[i].idx=i; }
	sort(massIdxPairs.begin(),massIdxPairs.end());
	// --- vector renaming vertices to sorted order
	vector<int> rename_old2new(specSet.size());
	for(unsigned int i=0; i<specSet.size(); i++) { rename_old2new[massIdxPairs[i].idx]=i; }
	// --- change the vertex numbers in aligns + redirect edges from lower index to higher index
	int a,b;
	for(unsigned int i=0; i<aligns.size(); i++) {
		a = rename_old2new[aligns[i].spec1];   b = rename_old2new[aligns[i].spec2];
		aligns[i].spec1=min(a,b);              aligns[i].spec2=max(a,b);
	}


	// Build the graph g based on these sorted spectra + find the components
	MSGraph g;  g.build(aligns);
	vector<int> numPaths;

ofstream outDebug("debug.txt");
for(unsigned int i=0; i<g.vNext.size(); i++){
	if(g.vNext[i].size()==0) continue;
	outDebug << "Vertex " << i <<" (previously "<<massIdxPairs[i].idx+1<<"): ";
	for(unsigned int j=0; j<g.vNext[i].size(); j++)
		outDebug << "(" << g.edges[g.vNext[i][j]][0] <<","<<g.edges[g.vNext[i][j]][1]<<")";
	outDebug << endl;
}
outDebug.close();
//cerr<<"Edges: ";
//for(unsigned int i=0; i<g.edges.size(); i++) cerr << "(" << g.edges[i][0] <<","<<g.edges[i][1]<<")"; cerr<<endl;

	vector<vector<int> > vSets;                // Used when pathLength>0
	vector<int> vSets_unique;

	vector<int> parentMasses(specSet.size());  // Used when pathLength<=0
	for(unsigned int i=0; i<specSet.size(); i++) parentMasses[i]=(int)round(10*specSet[massIdxPairs[i].idx].parentMass);
	list<vector<int> > vSets2;  list<vector<int> >::iterator vSets_iter;

	if (pathLength>0) {
		g.findConnectedPairs(pathLength, minNumPaths, vSets, numPaths);
		cout << "Number of sets: " << vSets.size() << endl;
	} else {
		g.findModPairs(parentMasses, pmTol, vSets2);
		cout << "Number of sets: " << vSets2.size() << endl;
	}

//    g.output_graphviz("test.txt");

//vSets_unique.resize(vSets.size()); for(unsigned int i=0; i<vSets.size(); i++) vSets_unique[i]=i;
	// Remove subsets
	if(pathLength>0) {
		if(vSets.size()==0) vSets_unique.resize(0); else {
			cout << "Filtering subsets ("<<vSets.size()<<" sets)... "; cout.flush();
			vector<TwoValues<int> > sizes(vSets.size());
			for(unsigned int i=0; i<vSets.size(); i++) sizes[i].set(vSets[i].size(),i);
			sort(sizes.begin(),sizes.end());

			vector<bool> setSpecs(specSet.size());   // Set of spectra in specIdx
			vector<bool> isSubset(vSets.size());     // Indicates whether a set is a subset
			for(unsigned int i=0; i<vSets.size(); i++) isSubset[i]=false;

			for(unsigned int i=vSets.size()-1; i>0; i--) {
				unsigned int setIdx = sizes[i][1];	  if(isSubset[setIdx]) continue;
	//cerr<<"Looking for subsets of set "<<setIdx<<endl;
				for(unsigned int j=0; j<specSet.size(); j++) setSpecs[j]=false;
				for(unsigned int j=0; j<vSets[setIdx].size(); j++) setSpecs[vSets[setIdx][j]]=true;

				for(unsigned int j=0; j<vSets.size(); j++) {
					if(isSubset[j] or j==setIdx or vSets[j].size()>vSets[setIdx].size()) continue;
					isSubset[j]=true;
					for(unsigned int k=0; isSubset[j] and k<vSets[j].size(); k++)
						isSubset[j]=isSubset[j] and setSpecs[vSets[j][k]];
	//if(isSubset[j]) cerr<<"set "<<j<<" is a subset of set "<<i<<endl;
				}
			}

			int kept=0;
			for(unsigned int i=0; i<vSets.size(); i++)
				if(!isSubset[i]) { vSets[kept].assign(vSets[i].begin(),vSets[i].end());  numPaths[kept++] = numPaths[i]; }
			vSets.resize(kept);   numPaths.resize(kept);
			cout << vSets.size() << " sets retained after removing subsets\n";
		}
	}

	// Output results renaming the vertices in the found components to their original numbers
	// File format is [N:number of sets], followed by N lines with [number of paths] [list of set members]
	if(pathLength>0) {
		ofstream outs("debug_output.txt");   outs << vSets_unique.size() << endl;
		for(unsigned int i=0; i<vSets_unique.size(); i++) {
			outs<<numPaths[vSets_unique[i]]<<" ";
//cerr<<"Set "<<vSets_unique[i]<<": ";
//for(unsigned int j=0; j<vSets[vSets_unique[i]].size(); j++) cerr<<vSets[vSets_unique[i]][j]<<" ";
			for(unsigned int j=0; j<vSets[vSets_unique[i]].size(); j++) outs<<massIdxPairs[vSets[vSets_unique[i]][j]].idx+1<<" ";
//cerr<<endl;
			outs<<endl;
		}
		outs.close();
	} else {
		ofstream outs("debug_output.txt");   outs << vSets2.size() << endl;
		for(vSets_iter = vSets2.begin(); vSets_iter!=vSets2.end(); vSets_iter++) {
			outs<<(*vSets_iter).size()<<" ";
			for(unsigned int j=0; j<(*vSets_iter).size(); j++) outs<<massIdxPairs[(*vSets_iter)[j]].idx+1<<" ";
			outs<<endl;
		}
		outs.close();
	}


	//
	//  Output binary results file (after recovering original zero-based spectrum indices)
	//
	if(pathLength>0) {
		for(unsigned int i=0; i<vSets.size(); i++)
			for(unsigned int j=0; j<vSets[i].size(); j++)
				vSets[i][j]=massIdxPairs[vSets[i][j]].idx;
		Save_binListArray<int,vector<int>,vector<int>::iterator>(resultsFN,vSets);
		if(params.paramPresent("OUTPUT_NUMPATHS")) Save_binArray(params.getValue("OUTPUT_NUMPATHS"),numPaths);
	} else {
		for(vSets_iter=vSets2.begin(); vSets_iter!=vSets2.end(); vSets_iter++)
			for(unsigned int j=0; j<vSets_iter->size(); j++)
				(*vSets_iter)[j]=massIdxPairs[(*vSets_iter)[j]].idx;
		Save_binListArray<int,vector<int>,vector<int>::iterator>(resultsFN,vSets2);
	}

	return(0);
}
