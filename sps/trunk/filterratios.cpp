#include <cstdio>
#include <iostream>
#include <vector>

using namespace std;

int main(int argc, char **argv) {
	if(argc!=6) {
		cerr << "Usage: "<<argv[0]<<" [alignsFN] [ratiosFN] minRatio [outputAlignsFN] [outputRatiosFN]\n";
		return(-1);
	}
	float minRatio = (float) atof(argv[3]);

	FILE *alignsIn = fopen(argv[1],"rb"), *alignsOut = fopen(argv[4],"wb"),
	     *ratiosIn = fopen(argv[2],"rb"), *ratiosOut = fopen(argv[5],"wb");
	if((int)alignsIn==0) { cerr<<"Error opening "<<alignsIn<<"!\n"; return(-1); }
	if((int)alignsOut==0) { cerr<<"Error opening "<<alignsOut<<"!\n"; return(-1); }
	if((int)ratiosIn==0) { cerr<<"Error opening "<<ratiosIn<<"!\n"; return(-1); }
	if((int)ratiosOut==0) { cerr<<"Error opening "<<ratiosOut<<"!\n"; return(-1); }

	float ratio1,ratio2;
	const unsigned int BUFFERSIZE=2048;
	char buffer[BUFFERSIZE], bufferRatios[BUFFERSIZE];
	fgets(buffer,BUFFERSIZE,alignsIn);
	unsigned int numEntries = (int)atof(buffer);
	cout<<"Number of entries in aligns/ratios files: "<<numEntries<<", filtering by ratio >= "<<minRatio<<endl;

	vector<bool> keep(numEntries);   unsigned int keepCount=0;
	for(unsigned int idx=0; idx<numEntries; idx++) {
		fgets(bufferRatios,BUFFERSIZE,ratiosIn);
		ratio1 = (float)atof(strtok(bufferRatios," ;"));
		ratio2 = (float)atof(strtok(NULL," ;"));
		if(min(ratio1,ratio2)>=minRatio) {
			fprintf(ratiosOut,"%f %f\n",ratio1,ratio2);
			keep[idx]=true;   keepCount++;
		} else keep[idx]=false;
	}
	fclose(ratiosIn); fclose(ratiosOut);
cerr<<"Done selecting ratios. Keeping "<<keepCount<<" pairwise alignments.\n";

	fprintf(alignsOut,"%d\n",keepCount);
	for(unsigned int idx=0; idx<numEntries; idx++) {
		fgets(buffer,BUFFERSIZE,alignsIn);
		if(keep[idx]) fprintf(alignsOut,"%s",buffer);
	}
	fclose(alignsIn); fclose(alignsOut);
	return(0);
}
