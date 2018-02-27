#include "alignment_scoring.h"
#include "aminoacid.h"
#include "spectrum.h"
#include "batch.h"
#include "filters.h"

#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>

using namespace std;

int main(int argc, char **argv){
    int numResults, numSpecs;
    vector<TwoValues<float> > ratios;
    vector<Results_ASP> results;
    float ratio;

    if (argc<7) {
        cerr << "Usage: filterRatio <prmSpecs.pkl> <input_results.txt> <RATIO> <output_results.txt> <output_ratios.txt> <output_means.txt>\n";
        return(-1);
    }

    ratio = atof(argv[3]);
    cout << "Keeping all the aligns with both ratios >= " << ratio << endl;

    clock_t startTime = clock();
    SpecSet specSet;
    numSpecs = specSet.LoadSpecSet_pkl(argv[1]);
    clock_t curTime = clock();
    cout << "Specs loaded. Num entries : " << numSpecs << ". " << (curTime-startTime)/CLOCKS_PER_SEC << " secs\n";

    startTime = curTime;
    numResults = Load_resultsASP(argv[2], results);     curTime = clock();
    cout << "Results loaded. Num entries : " << numResults << ". " << (curTime-startTime)/CLOCKS_PER_SEC << " secs\n";

    startTime = curTime;
    FilterMatchRatioASP(specSet, results, 7, ratios);
    curTime = clock();
    cout << "Done filtering in  " << (curTime-startTime)/CLOCKS_PER_SEC << " secs. Writing output..."; cout.flush();

    int i,count=0;
    for(i=0; i<=numResults; i++) if (ratios[i][0]>=ratio & ratios[i][1]>=ratio) count++;

    startTime = curTime;
    ofstream out(argv[4], ios::binary), out2(argv[5], ios::binary);  out << count << endl;
    for(i=0; i<numResults; i++)
        if (ratios[i][0]>=ratio & ratios[i][1]>=ratio)
            { results[i].output(out,';');   out2 << ratios[i][0] << " " << ratios[i][1] << endl; }
    out.close(); out2.close();
    curTime = clock(); cout << "done in " << (curTime-startTime)/CLOCKS_PER_SEC << " secs\n";

    startTime = curTime;
    cout << "Computing means..."; cout.flush();
    vector<TwoValues<float> > scoreMeans;
    ComputeScoreMeans(specSet, results, scoreMeans);
    curTime = clock(); cout << "done in " << (curTime-startTime)/CLOCKS_PER_SEC << " secs\n";

    cout << "Writing means... "; cout.flush();
    out.open(argv[6], ios::binary);
    for(i=0; i<scoreMeans.size(); i++) out << scoreMeans[i][0] << " " << scoreMeans[i][1] << endl;
    out.close();
    cout << "done!\n";

}
