#include "inputParams.h"
#include "spectrum.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>

using namespace std;

class IndexFile {
	static const unsigned int bufferSize = 4194304;  // 4 Mb line read buffer - supposed to be larger than will ever be needed
	char *buffer;
	public:
	vector<char *> filenames;
	vector<vector<unsigned int > > indices;

	IndexFile() { filenames.resize(0);   indices.resize(0); buffer = (char *)malloc(sizeof(char)*IndexFile::bufferSize); }
	~IndexFile() { reset(); if(buffer) free(buffer); }

	unsigned int size() { return filenames.size(); }
	unsigned int totalCount() { unsigned int tc=0; for(unsigned int i=0;i<indices.size();i++) tc+=indices[i].size(); return tc; }
	void reset() { for(unsigned int i=0;i<filenames.size();i++) if(filenames[i]) free(filenames[i]); filenames.resize(0);   indices.resize(0); }

	short Load(char *filename);
};

int main(int argc, char **argv){
    // Get input parameters and check minimum parameter set
    InputParams params; bool paramsOk;
	if(argc<=1) paramsOk=params.readParams("extract.params");	else paramsOk=params.readParams(argv[1]);
	if(!paramsOk) {
		cerr << "Error opening parameters file ";
		if(argc<=1) cerr << "extract.params\n"; else cerr << argv[1] << endl;
		return -1;
	}
	vector<char *> paramStrs;   paramStrs.resize(2);
	paramStrs[0] = "INPUT_INDEX";
	paramStrs[1] = "OUTPUT_SPECS";
	if(!params.confirmParams(paramStrs)) {
		cerr << "ERROR: Parameters file ";
		if(argc==1) cerr<<"extract.params"; else cerr<<argv[1];
		cerr << " is incomplete. One of the following is missing: INPUT_INDEX, OUTPUT_SPECS\n";
		return -1;
	}

	unsigned int fileIdx, specIdx;
	SpecSet specs;
    IndexFile indexFile;

	char bufferFN[4096];   int numFiles;
	vector<ofstream *> outputFiles;   vector<unsigned int> outputSpecs;
	ifstream buffersFile(params.getValue(paramStrs[1]), ios::binary);
	if (!buffersFile) { cerr << "ERROR: cannot open " << params.getValue(paramStrs[1]) << "\n";  return -1; }
	buffersFile.getline(bufferFN,4096,'\n'); numFiles = (int)atoi(bufferFN);  if(numFiles<=0) { cerr << "ERROR: " << numFiles << " output files?\n";  return -1; }
	outputFiles.resize(numFiles);    outputSpecs.resize(numFiles);
	for(fileIdx=0; fileIdx<outputFiles.size(); fileIdx++) {
		buffersFile.getline(bufferFN,4096,'\n');
		outputFiles[fileIdx] = new ofstream(bufferFN, ios::trunc | ios::binary);   outputSpecs[fileIdx]=0;
		if (!outputFiles[fileIdx] or !(*outputFiles[fileIdx])) { cerr << "ERROR: cannot open output file " << bufferFN << "!\n";  return -1; }
	}

    if(indexFile.Load(params.getValue(paramStrs[0]))<=0) return -1;

	// Copy the spectra to the output files
	float processed=0;
	for(fileIdx=0; fileIdx<indexFile.size(); fileIdx++) {
		cout<<"Processing file "<<fileIdx+1<<"/"<<indexFile.size()<<" ("<<indexFile.filenames[fileIdx]<<")...";
		if(specs.LoadSpecSet_pklbin(indexFile.filenames[fileIdx])<=0) { cerr << "ERROR: cannot open input file " << indexFile.filenames[fileIdx] << "!\n";  return -1; }
		for(specIdx=0; specIdx<indexFile.indices[fileIdx].size(); specIdx++)
			if(indexFile.indices[fileIdx][specIdx]>0) {
				if(outputSpecs[fileIdx]>0) (*outputFiles[indexFile.indices[fileIdx][specIdx]-1])<<"\n"; outputSpecs[fileIdx]++;
				specs[specIdx].output_ms2(*outputFiles[indexFile.indices[fileIdx][specIdx]-1]);
			}
		processed+=indexFile.indices[fileIdx].size();
		cout<<" done ("<<processed/(float)indexFile.totalCount()<<"% completed)\n";
		remove(indexFile.filenames[fileIdx]);
	}

	// Close the output files
	for(fileIdx=0; fileIdx<outputFiles.size(); fileIdx++) outputFiles[fileIdx]->close();
cerr<<"Output files successfully closed\n"; cerr.flush();
	for(fileIdx=0; fileIdx<outputFiles.size(); fileIdx++) delete outputFiles[fileIdx];

/*
	unsigned int fileIdx, finalIdx=0, specIdx;
    SpecSet final, current;
    IndexFile indexFile;

    if(indexFile.Load(params.getValue(paramStrs[0]))<=0) return -1;
    final.resize(indexFile.totalCount());

    for(fileIdx=0; fileIdx<indexFile.size(); fileIdx++) {
    	if(current.LoadSpecSet_pklbin(indexFile.filenames[fileIdx])<=0)
    		{ cerr<<"ERROR loading "<<indexFile.filenames[fileIdx]<<"!\n"; return -1; }

    	// Move spectra from current to final
    	for(specIdx=0; specIdx<indexFile.indices[fileIdx].size(); specIdx++) {
    		if(indexFile.indices[fileIdx][specIdx]>current.size())
    			{ cerr<<"ERROR: Trying to extract spectrum index "<<indexFile.indices[fileIdx][specIdx]<<" from a set of size "<<current.size()<<"!\n"; return -1; }

			if(finalIdx==final.size()) { cerr<<"ERROR: Trying to output more spectra than expected from "<<params.getValue(paramStrs[0])<<"! Incorrect file format?\n"; return -1; }
			final[finalIdx++] = current[indexFile.indices[fileIdx][specIdx]];
			current[indexFile.indices[fileIdx][specIdx]].resize(0);
    	}
    	current.resize(0);
    }
	if(finalIdx<final.size()) { cerr<<"Warning: output contains "<<finalIdx<<" spectra instead of the expected "<<final.size()<<"?\n"; return -1; }

    final.SaveSpecSet_pklbin(params.getValue(paramStrs[1]));
*/
    return(0);
}

short IndexFile::Load(char *filename) {
	reset();
    if(buffer==(char *)0) { cerr<<"ERROR IndexFile::buffer is invalid - can't process "<<filename<<"!\n"; return -1; }
    FILE *input = fopen(filename,"rb");
    if (!input or ferror(input)) { cerr << "ERROR opening " << filename << "!\n"; return -1; }

	int numFiles, numIndices, readSpecIndex;
	unsigned int fileIdx, specIdx;
	char *curToken;

	if(!fgets(buffer,bufferSize,input)) { cerr << "ERROR reading " << filename << "!\n"; return -1; }
	else buffer[strlen(buffer)]=0;
	numFiles = atoi(buffer);  if(numFiles<=0) { cerr << "ERROR reading " << filename << " - invalid number of files: "<<numFiles<<"!\n"; return -1; }

	filenames.resize(numFiles);   for(fileIdx=0; fileIdx<(unsigned int)numFiles; fileIdx++) filenames[fileIdx]=(char*)0;
	indices.resize(numFiles);
	for(fileIdx=0; fileIdx<(unsigned int)numFiles; fileIdx++) {
		if(!fgets(buffer,bufferSize,input)) { cerr << "ERROR reading " << filename << "!\n"; return -1; }
		else buffer[strlen(buffer)]=0;

		// Get filename
		curToken = strtok(buffer,":"); if(!curToken) { cerr << "ERROR: invalid filename in line " << fileIdx+2 << "!\n"; return -1; }
		filenames[fileIdx] = (char *)malloc(strlen(curToken)+1);
		strcpy(filenames[fileIdx],curToken);

		// Get number of indices
		curToken = strtok(NULL,":"); if(!curToken) { cerr << "ERROR: missing number of spectra in line " << fileIdx+2 << "!\n"; return -1; }
		numIndices = atoi(curToken); if(numIndices<=0) { cerr << "ERROR: invalid number of spectra in line " << fileIdx+2 << "!\n"; return -1; }
		indices[fileIdx].resize(numIndices);

		// Get all spectrum indices
		for(specIdx=0; specIdx<(unsigned int)numIndices; specIdx++) {
			curToken = strtok(NULL,":"); if(!curToken) { cerr << "ERROR: missing spectrum index in line " << fileIdx+2 << "!\n"; return -1; }
			readSpecIndex = atoi(curToken); if(readSpecIndex<0) { cerr << "ERROR: invalid spectrum index in line " << fileIdx+2 << "!\n"; return -1; }
			indices[fileIdx][specIdx] = readSpecIndex;
		}
	}

	fclose(input);
	return(1);
}
