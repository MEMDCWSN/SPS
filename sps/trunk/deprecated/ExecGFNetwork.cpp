/*
 * ExecGFNetwork.cpp
 *
 *  Created on: Feb 18, 2011
 *      Author: cboucher@eng.ucsd.edu
 */
// Header Includes
#include "ExecGFNetwork.h"
#include "stdlib.h"

// Module Includes
#include "Logger.h"
#include "FileUtils.h"
#include "ParameterList.h"
#include "PeptideSpectrumMatchSet.h"
#include "aminoacid.h"
#include <sstream>
#include <math.h>
#include <time.h>
#include <iostream>
#include <algorithm>
#include <vector>

using namespace specnets;

namespace specnets {

const unsigned int AAcount = 20;

const int AAmassesInt[] = {
		71,
		156,
		115,
        114,
        160,
        129,
        128,
        57,
        137,
        113,
        113,
        128,
        131,
        147,
        97,
        87,
        101,
        186,
        163,
        99};

const float AAmasses[] = {
		71.0371137870,
		156.1011110260,
		115.0269430310,
        114.0429274460,
        160.0306482000,
        129.0425930950,
        128.0585775100,
        57.0214637230,
        137.0589118610,
        113.0840639790,
        113.0840639790,
        128.0949630160,
        131.0404846050,
        147.0684139150,
        97.0527638510,
        87.0320284090,
        101.0476784730,
        186.0793129520,
        163.0633285370,
        99.0684139150 };

// -------------------------------------------------------------------------

ExecGFNetwork::ExecGFNetwork(void) : ownInput(false) {
	m_name = "ExecGFNetwork";
	m_type = "ExecGFNetwork";

}

// -------------------------------------------------------------------------
ExecGFNetwork::ExecGFNetwork(SpecSet * inSpectra, int dim) : ownInput(false) {
	m_name = "ExecGFNetwork";
	m_type = "ExecGFNetwork";
	numberSpectra = dim;
	scores = inSpectra;

}

// -------------------------------------------------------------------------

ExecGFNetwork::ExecGFNetwork(const ParameterList & inputParams) :
	ExecBase(inputParams), ownInput(false) {
	m_name = "ExecGFNetwork";
	m_type = "ExecGFNetwork";
}

// -------------------------------------------------------------------------

ExecGFNetwork::~ExecGFNetwork(void) {
}

// -------------------------------------------------------------------------

ExecBase * ExecGFNetwork::clone(const ParameterList & inputParams) const {
	return new ExecGFNetwork(inputParams);
}

// -------------------------------------------------------------------------

bool ExecGFNetwork::saveInputData(std::vector<std::string> & filenames) {
	return true;
}


// -------------------------------------------------------------------------

bool ExecGFNetwork::saveOutputData(void) {
	return true;
}

//-----------------------------------------------------------------------------

bool ExecGFNetwork::loadOutputData(void) {
	return true;
}
// -------------------------------------------------------------------------

vector<ExecBase *> const & ExecGFNetwork::split(int numSplit) {
	m_subModules.resize(0);
	return m_subModules;
}

// -------------------------------------------------------------------------

bool ExecGFNetwork::merge(void) {
	return false;
}

// -------------------------------------------------------------------------

bool ExecGFNetwork::validateParams(std::string & error) {
	return true;
}


//-----------------------------------------------------------------------------

void ExecGFNetwork::computeMassJumps()
{
	int numberJumps = 0;
	int jumpsIndex = 0;

	for(int i = 0; i < numberSpectra - 1; i++)
	{
		if(types[i] == FULL)
			continue;

		vector<char> sequence = vector<char> (peptideStrings[i].begin(), peptideStrings[i].end());
		vector<float> masses;
		getMasses(sequence, masses);
		int mass = 0;
		for (int k = 0; k < masses.size(); k++)
			mass += (int) round(masses[k] * 0.9995);

		if(types[i] == PRE)
		{
			jumps[jumpsIndex] = mass;
			numberJumps++;
		}
		if(types[i] == POST || types[i] == SUB)
		{
				size_t found;
				found=peptideStrings[numberSpectra - 1].find(peptideStrings[i]);
				string peptidePiece = peptideStrings[numberSpectra - 1].substr(0, (int)found);
				DEBUG_MSG(peptidePiece);

				vector<char> sequence = vector<char> (peptidePiece.begin(), peptidePiece.end());
				vector<float> masses;
				getMasses(sequence, masses);
				int mass2 = 0;
				for (int k = 0; k < masses.size(); k++)
					mass2 += (int) round(masses[k] * 0.9995);

				jumps[jumpsIndex] = mass2;
				numberJumps++;

				if(types[i] == SUB)
				{
					numberJumps++;
					jumpsIndex++;
					jumps[jumpsIndex] = mass2 + mass;
				}
		}
		jumpsIndex++;

	}

}

//-----------------------------------------------------------------------------

void ExecGFNetwork::outputForPairs(vector<int> inIndices, vector<spectra_type_t> inTypes, string peptide, vector<int> inJumps, int parentMass) {

	partialAnnotations.resize(2);
	partialPeptideScores.resize(2);
	partialSpecProbabilities.resize(2);

	numberSpectra = 2;
	int offset = 200;
	int rowSize = parentMass + 1;
	int colSize = 450 + offset + 1;

	jumps = new int[2*numberSpectra];
	indices = new int[numberSpectra];
	types = new spectra_type_t[numberSpectra];
	peptideScores = new int[numberSpectra];

	for(int i = 0; i < numberSpectra; i++)
	{
		indices[i] = inIndices[i];
		types[i] = inTypes[i];
	}
	for(int i = 0; i < numberSpectra*2; i++)
		jumps[i] = inJumps[i];

	int jumpIndex = 0;
	double iScore_0 = 0;
	double iScore_1 = 0;
	for(int spectraIndex = 0; spectraIndex < numberSpectra; spectraIndex++)
	{
		string partial_peptide = "";
		if(types[spectraIndex] == FULL)
		{
			partial_peptide = computePartialPeptide(peptide, types[spectraIndex], 0, 0);
			peptideScores[spectraIndex] = scorePeptide(partial_peptide, indices[spectraIndex]) + offset;
		}
		else if (types[spectraIndex] == PRE || types[spectraIndex] == POST)
		{
			partial_peptide = computePartialPeptide(peptide, types[spectraIndex], jumps[jumpIndex], 0);
			peptideScores[spectraIndex] = scorePeptide(partial_peptide, indices[spectraIndex]) + offset;
			jumpIndex++;
		}
		else// if (types[spectraIndex] == SUB)
		{
			partial_peptide = computePartialPeptide(peptide, types[spectraIndex], jumps[jumpIndex], jumps[jumpIndex + 1]);
			peptideScores[spectraIndex] = scorePeptide(partial_peptide, indices[spectraIndex]) + offset;
			jumpIndex += 2;
		}
		cout << partial_peptide;
		if(spectraIndex == 0)
		{
			iScore_0 = computeSpectraProbability(indices[spectraIndex], partial_peptide, offset, colSize);
			cout << ", ";
			partialAnnotations[0] = partial_peptide;
		}
		else
		{
			iScore_1 = computeSpectraProbability(indices[spectraIndex], partial_peptide, offset, colSize);
			partialAnnotations[1] = partial_peptide;
			cout << "\t";
		}
	}

	cout << iScore_0 << ", " << iScore_1 << "\t";
	cout << peptideScores[0] - offset << ", " << peptideScores[1] - offset << "\n";

	partialPeptideScores[0] = peptideScores[0] - offset;
	partialPeptideScores[1] = peptideScores[1] - offset;

	partialSpecProbabilities[0] = iScore_0;
	partialSpecProbabilities[1] = iScore_1;

	delete [] jumps;
	delete [] indices;
	delete [] types;
	delete [] peptideScores;

}


//-----------------------------------------------------------------------------

void ExecGFNetwork::outputPSMTriple(vector<int> inIndices, vector<spectra_type_t> inTypes, string peptide,
										vector<int> inJumps, int parentMass) {

	partialAnnotations.resize(3);
	partialPeptideScores.resize(3);
	partialSpecProbabilities.resize(3);

	int numberSpectra = 3;
	int offset = 200;
	int rowSize = parentMass + 1;
	int colSize = 650 + offset + 1;

	jumps = new int[2*numberSpectra];
	indices = new int[numberSpectra];
	types = new spectra_type_t[numberSpectra];
	peptideScores = new int[numberSpectra];

	for(int i = 0; i < numberSpectra; i++)
	{
		indices[i] = inIndices[i];
		types[i] = inTypes[i];
	}
	for(int i = 0; i < 2*numberSpectra; i++)
		jumps[i] = inJumps[i];

	double iScore_0 = 0;
	double iScore_1 = 0;
	double iScore_2 = 0;

	int jumpIndex = 0;
	for(int spectraIndex = 0; spectraIndex < numberSpectra; spectraIndex++)
	{
			string partial_peptide = "";
			if(types[spectraIndex] == FULL)
			{
				partial_peptide = computePartialPeptide(peptide, types[spectraIndex], 0, 0);
				peptideScores[spectraIndex] = scorePeptide(partial_peptide, indices[spectraIndex]) + offset;
			}
			else if (types[spectraIndex] == PRE || types[spectraIndex] == POST)
			{
				partial_peptide = computePartialPeptide(peptide, types[spectraIndex], jumps[jumpIndex], 0);
				peptideScores[spectraIndex] = scorePeptide(partial_peptide, indices[spectraIndex]) + offset;
				jumpIndex++;
			}
			else// if (types[spectraIndex] == SUB)
			{
				partial_peptide = computePartialPeptide(peptide, types[spectraIndex], jumps[jumpIndex], jumps[jumpIndex + 1]);
				peptideScores[spectraIndex] = scorePeptide(partial_peptide, indices[spectraIndex]) + offset;
				jumpIndex += 2;
			}
			cout << partial_peptide;
			if(spectraIndex == 0)
			{
				iScore_0 = computeSpectraProbability(indices[spectraIndex], partial_peptide, offset, colSize);
				cout << ", ";
				partialAnnotations[0] = partial_peptide;
			}
			else if(spectraIndex == 1)
			{
				iScore_1 = computeSpectraProbability(indices[spectraIndex], partial_peptide, offset, colSize);
				partialAnnotations[1] = partial_peptide;
				cout << ", ";
			}
			else
			{
				iScore_2 = computeSpectraProbability(indices[spectraIndex], partial_peptide, offset, colSize);
				partialAnnotations[2] = partial_peptide;
				cout << "\t";
			}
	}

	cout << iScore_0 << ", " << iScore_1 << ", " << iScore_2 << "\t";
	cout << peptideScores[0] - offset << ", " << peptideScores[1] - offset <<  ", " << peptideScores[2] - offset  << "\n";

	partialPeptideScores[0] = peptideScores[0] - offset;
	partialPeptideScores[1] = peptideScores[1] - offset;
	partialPeptideScores[2] = peptideScores[2] - offset;

	partialSpecProbabilities[0] = iScore_0;
	partialSpecProbabilities[1] = iScore_1;
	partialSpecProbabilities[2] = iScore_2;

	delete [] jumps;
	delete [] indices;
	delete [] types;
	delete [] peptideScores;

}


//-----------------------------------------------------------------------------

vector<double> ExecGFNetwork::scorePSMPairs(vector<int> inIndices, vector<spectra_type_t> inTypes, vector<string> peptides,
										vector<int> inJumps, int parentMass) {

	vector<double> retScores;

	numberSpectra = 2;
	int offset = 200;
	int rowSize = parentMass + 1;
	int colSize = 250 + offset + 1;

	jumps = new int[2*numberSpectra];
	indices = new int[numberSpectra];
	types = new spectra_type_t[numberSpectra];
	peptideScores = new int[numberSpectra];

	for(int i = 0; i < numberSpectra; i++)
	{
		indices[i] = inIndices[i];
		types[i] = inTypes[i];
	}
	for(int i = 0; i < numberSpectra*2; i++)
		jumps[i] = inJumps[i];

	constructPSMMatrix(offset, colSize, parentMass);


	for(int i = 0; i < peptides.size(); i++)
	{
		int jumpIndex = 0;
		for(int spectraIndex = 0; spectraIndex < numberSpectra; spectraIndex++)
		{
			string partial_peptide = "";

			if(types[spectraIndex] == FULL)
			{
				partial_peptide = computePartialPeptide(peptides[i], types[spectraIndex], 0, 0);
				peptideScores[spectraIndex] = scorePeptide(partial_peptide, indices[spectraIndex]) + offset;
			}
			else if (types[spectraIndex] == PRE || types[spectraIndex] == POST)
			{
				partial_peptide = computePartialPeptide(peptides[i], types[spectraIndex], jumps[jumpIndex], 0);
				peptideScores[spectraIndex] = scorePeptide(partial_peptide, indices[spectraIndex]) + offset;
				jumpIndex++;
			}
			else// if (types[spectraIndex] == SUB)
			{
				partial_peptide = computePartialPeptide(peptides[i], types[spectraIndex], jumps[jumpIndex], jumps[jumpIndex + 1]);
				peptideScores[spectraIndex] = scorePeptide(partial_peptide, indices[spectraIndex]) + offset;
				jumpIndex += 2;
			}
	//		DEBUG_MSG(partial_peptide);
		}

		double retSpecProb = 0;
		if(peptideScores[0] >= 0 && peptideScores[1] >= 0)
		{
			for (int t0 = peptideScores[0]; t0 < colSize; t0++)	{
				for (int t1 = peptideScores[1]; t1 < colSize; t1++)	{

					/*if(x[parentMass][t0][t1] > 0)
					{
						DEBUG_VAR(x[parentMass][t0][t1]);
						DEBUG_VAR(t1);
						DEBUG_VAR(t0);
					}*/
					retSpecProb += x[parentMass][t0][t1];
				}
			}
		}
		retScores.push_back( retSpecProb );


	}

	for(int d1 = 0; d1 < rowSize; d1++)
	{
		for(int d2 = 0; d2 < colSize; d2++)
	    {
			delete [] x[d1][d2];
	    }
	    delete [] x[d1];
	}

	delete [] x;

	delete [] jumps;
	delete [] indices;
	delete [] types;
	delete [] peptideScores;

	return retScores;

}


//-----------------------------------------------------------------------------

vector<double> ExecGFNetwork::scorePSMTriple(vector<int> inIndices, vector<spectra_type_t> inTypes, vector<string> peptides,
										vector<int> inJumps, int parentMass) {

	int numberSpectra = 3;
	int offset = 200;
	int rowSize = parentMass + 1;
	int colSize = 650 + offset + 1;

	jumps = new int[2*numberSpectra];
	indices = new int[numberSpectra];
	types = new spectra_type_t[numberSpectra];
	peptideScores = new int[numberSpectra];

	for(int i = 0; i < numberSpectra; i++)
	{
		indices[i] = inIndices[i];
		types[i] = inTypes[i];
	}
	for(int i = 0; i < 2*numberSpectra; i++)
		jumps[i] = inJumps[i];

	vector<double> retScores;

	constructPSMMatrixForTriple(offset, colSize, parentMass);

	for(int i = 0; i < peptides.size(); i++)
	{

		DEBUG_MSG(peptides[i]);

		int jumpIndex = 0;
		for(int spectraIndex = 0; spectraIndex < numberSpectra; spectraIndex++)
		{
			//DEBUG_MSG(spectraIndex);

			//DEBUG_MSG(computeSpectraProbability(indices[spectraIndex], masses, offset, colSize));
			if(types[spectraIndex] == FULL)
			{
				peptideScores[spectraIndex] = scorePeptide(computePartialPeptide(peptides[i], types[spectraIndex], 0, 0), indices[spectraIndex]) + offset;
			}
			else if (types[spectraIndex] == PRE || types[spectraIndex] == POST)
			{
				peptideScores[spectraIndex] = scorePeptide(computePartialPeptide(peptides[i], types[spectraIndex], jumps[jumpIndex], 0), indices[spectraIndex]) + offset;
				jumpIndex++;
			}
			else if (types[spectraIndex] == SUB)
			{
				peptideScores[spectraIndex] = scorePeptide(computePartialPeptide(peptides[i], types[spectraIndex], jumps[jumpIndex], jumps[jumpIndex + 1]), indices[spectraIndex]) + offset;
				jumpIndex += 2;
			}
		}
		DEBUG_MSG("peptideScores[0]: " << peptideScores[0] << " peptideScores[1]: " << peptideScores[1] << " peptideScores[2]: " << peptideScores[2]);


		retScores.push_back(scoreSpectraTriple(peptideScores[0], peptideScores[1], peptideScores[2], colSize, parentMass));
	}

	int ROWS = 187;
	for( int r = 0 ; r < ROWS ; r++ )
			delete [] dpTripleMatrix[r];
	delete [] dpTripleMatrix;

	delete [] jumps;
	delete [] indices;
	delete [] types;
	delete [] peptideScores;

	return retScores;

}

//-----------------------------------------------------------------------------

bool ExecGFNetwork::invoke(void) {

	int offset = 250;
	int lastIndex = numberSpectra - 1;
	clock_t t1,t2;

	jumps = new int[2*numberSpectra];
	peptideStrings = new string[numberSpectra];
	indices = new int[numberSpectra];
	types = new spectra_type_t[numberSpectra];
	peptideScores = new int[numberSpectra];
	singleSpecScore = new float[numberSpectra];


	for(int i = 0; i < peptide_results.size(); i++)
	{
		// Obtain the first peptide/spectra and store it in the last position of the array
		indices[lastIndex] = peptide_results[i]->m_scanNum;
		peptideStrings[lastIndex] = peptide_results[i]->m_annotation;
		types[lastIndex] = FULL;

		// Calculate the mass of the first peptide/spectra
		currentParentMass = 0;
		vector<char> sequence = vector<char> (peptideStrings[lastIndex].begin(), peptideStrings[lastIndex].end());
		vector<float> masses;
		getMasses(sequence, masses);
		for (int k = 0; k < masses.size(); k++)
			currentParentMass += (int) round(masses[k] * 0.9995);

		currentParentMass = (int) (currentParentMass);

		int k = 0;
		int j = i + 1;
		bool compute = false;

		while(j < peptide_results.size())
		{

			string temp = peptide_results[j]->m_annotation;
			size_t found=peptideStrings[lastIndex].find(temp);

			if(peptideStrings[lastIndex].compare(temp) == 0)
			{
				//	DEBUG_MSG("FULL");
			/*	types[k] = FULL;
				indices[k] = peptide_results[j]->m_scanNum;
				peptideStrings[k] = peptide_results[j]->m_annotation;
				k++;
				if(k == numberSpectra - 1)
					compute = true;*/
			}
			else if(found==string::npos)
			{
			//	DEBUG_MSG("NONE");
			}
			else if((int)found == 0 )
			{
			//		DEBUG_MSG("PRE");
			/*		types[k] = PRE;
					indices[k] = peptide_results[j]->m_scanNum;
					peptideStrings[k] = peptide_results[j]->m_annotation;
					k++;
					if(k == numberSpectra - 1 )
						compute = true;*/
			}
			else if(temp.length() + (int)found == peptideStrings[numberSpectra - 1].length())
			{
				//DEBUG_MSG("POST");
		/*		types[k] = POST;
				indices[k] = peptide_results[j]->m_scanNum;
				peptideStrings[k] = peptide_results[j]->m_annotation;
				k++;
				if(k == numberSpectra - 1)
					compute = true;*/
			}
			else
			{
				//DEBUG_MSG("SUB");
				types[k] = SUB;
				indices[k] = peptide_results[j]->m_scanNum;
				peptideStrings[k] = peptide_results[j]->m_annotation;
				k++;
				if(k == numberSpectra - 1 )
					compute = true;
			}

			j++;

			if( compute )
			{
				DEBUG_MSG("Calculate individual scores: ");

				for(int spectraIndex = 0; spectraIndex < numberSpectra; spectraIndex++)
				{

				//	DEBUG_MSG(peptideStrings[spectraIndex]);
					singleSpecScore[spectraIndex] = computeSpectraProbability(indices[spectraIndex], peptideStrings[spectraIndex], offset, 150 + offset + 1);

		//			DEBUG_MSG(singleSpecScore[spectraIndex]);

				}
				computeMassJumps();

				t1=clock();
				DEBUG_MSG(computeMultiSpectraProbability(offset, 250 + offset + 1));
				t2=clock();
				double seconds = ((float)t2-(float)t1) / CLOCKS_PER_SEC;
				DEBUG_MSG("seconds:");
				DEBUG_MSG(seconds);
				k = 0;
				compute = false;
			}

		}
	}
	delete [] jumps;
	delete [] peptideStrings;
	delete [] indices;
	delete [] types;
	delete [] peptideScores;
	delete [] singleSpecScore;

	return true;
}
//-----------------------------------------------------------------------------

bool ExecGFNetwork::loadInputData(void) {

	if (m_params.exists("DIMENSION")) {
		numberSpectra = atoi(m_params.getValue("DIMENSION").c_str());
	} else {
		ERROR_MSG("DIMENSION parameter missing.  Check param file.");
		return false;
	}

	if (m_params.exists("PRM_SPECTRA")) {
		if (!scores->LoadSpecSet_mgf(
				m_params.getValue("PRM_SPECTRA").c_str())) {
			ERROR_MSG("Could not load " << m_params.getValue("PRM_SPECTRA"));
			return false;
		}
	}
	if (scores->size() == 0) {
		ERROR_MSG("Could not load PRM spectra. Filename: [" << m_params.getValue("PRM_SPECTRA") << "]") ;
		return false;
	}

	DEBUG_MSG("Loading specs complete. Num specs [" << scores->size() << "]");


	if (m_params.exists("INSPEC_RESULTS")) {
		if (!peptide_results.loadInspectResultsFile(m_params.getValue(
				"INSPEC_RESULTS").c_str(), false)) {
			ERROR_MSG("ERROR loading " << m_params.getValue("INSPEC_RESULTS"));
			return false;
		}
	}

	if (peptide_results.size() == 0) {
		ERROR_MSG("Set of peptide results has size is 0!");
		return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
// Calculates the score of a single spectra

float ExecGFNetwork::computeSpectraProbability(int i, string peptide, int offset, int threshold) {

	AAJumps myjumps(1);
	vector<float> mymasses;
	float mass_float = myjumps.getPeptideMass(peptide);
	int currentParentMass = (int) round(mass_float*0.9995);

	int ROWS = currentParentMass + 1;
	int COLS = threshold;

	double **x = 0;
	x = new double *[ROWS] ;

	for( int r = 0 ; r < ROWS ; r++ )
		x[r] = new double[COLS];


	for (int r = 0; r < ROWS; r++)
		for (int c = 0; c < COLS; c++)
			x[r][c] = 0;

	x[0][offset] = 1;

	for (int m = 1; m < ROWS; m++) {

		//int currScore = (int)round((*scores)[i][m - 1][1]);
		int currScore = ( m < currentParentMass ) ? (int)(*scores)[i][m - 1][1] : 0;
		int start = (currScore < 0) ? 0 : currScore;
		int finish = min(COLS, threshold + currScore);

		for (int j = 0; j < AAcount; j++) {

			if (m - AAmassesInt[j] >= 0)
			{

				if(m - AAmassesInt[j] >= currentParentMass)
				{
					currScore = 0;
					start = 0;
					finish = threshold;

				}
				int massIndex = m - AAmassesInt[j];

				for (int scoreIndex = start; scoreIndex < finish; scoreIndex++)
				{
					x[m][scoreIndex] += x[m - AAmassesInt[j]][scoreIndex - currScore] *0.05;
				}
			}
		}
	}

	int peptideScore = scorePeptide(peptide, i)  + offset;


	double retSpecProb = 0;

	for (int t = peptideScore; t < threshold; t++)
	{

		retSpecProb += x[currentParentMass][t];

	}

	for( int r = 0 ; r < ROWS ; r++ )
		delete [] x[r] ;
	delete [] x;

	return retSpecProb;
}


//-----------------------------------------------------------------------------
// Calculates the score of two spectra
// finished.
// uses indices, types, jumps, numberSpectra, currentParentMass
void ExecGFNetwork::constructPSMMatrix(int offset, int threshold, int parentMass)
{
		int rowSize = parentMass + 1;
		int colSize = threshold;
		x = new double**[rowSize];

		for( int r = 0 ; r < rowSize ; r++ ) {

			x[r] = new double*[colSize];

			for( int r2 = 0 ; r2 < colSize ; r2++ ) {

				x[r][r2] = new double[colSize];

				for( int r3 = 0 ; r3 < colSize ; r3++ ) {

					x[r][r2][r3] = 0.0;
				}
			}
		}

		x[0][offset][offset] = 1;

		int * currScore = new int[2];
		int * start = new int[2];
		int * finish = new int[2];

		int bounds[numberSpectra][2];
		int jumpIndex = 0;
		for(int n = 0; n < 2; n++)
		{
			if(types[n] == FULL)
			{
				bounds[n][0] = 0;
				bounds[n][1] = parentMass;
				//DEBUG_MSG("FULL");
			}
			else if (types[n] == PRE)
			{
				bounds[n][0] = 0;
				bounds[n][1] = jumps[jumpIndex++];
			//	DEBUG_MSG("PRE");
			}
			else if (types[n] == POST)
			{
				bounds[n][0] = jumps[jumpIndex++];
				bounds[n][1] = parentMass;
		//		DEBUG_MSG("POST");
			}
			else if (types[n] == SUB)
			{
				bounds[n][0] = jumps[jumpIndex++];
				bounds[n][1] = jumps[jumpIndex++];
			}
		}

	/*	for(int i = 0; i < numberSpectra; i++)
		{
			DEBUG_MSG("Bounds " << i << " "<< bounds[i][0] << ", " << bounds[i][1]);
		}*/


		for (int m = 1; m < rowSize; m++) {

			for(int n = 0; n < 2; n++)	{
				int scoreOffset = m - bounds[n][0] - 1;

				//delta[n] = (m - 1 < (*scores)[indices[n]].size()) ? (int)(*scores)[indices[n]][m - 1][1] : 0;
				currScore[n] = (scoreOffset >= 0 && m < bounds[n][1] ) ? (int)(*scores)[indices[n]][scoreOffset][1] : 0;
				//rrScore[n] = (scoreOffset >= 0 && scoreOffset < (*scores)[indices[n]].size() ) ? (int)(*scores)[indices[n]][scoreOffset][1] : 0;
				start[n] = (currScore[n] < 0) ? 0 : currScore[n];
				finish[n] = min(colSize, threshold + currScore[n]);
			}

			for (int j = 0; j < AAcount; j++) {

				if(m - AAmassesInt[j] >= 0 && checkBoundaries(m, AAmassesInt[j]) )
				{

					int massIndex = m - AAmassesInt[j];
					for(int n = 0; n < 2; n++)
					{
						if(massIndex >= bounds[n][1])
						{
							currScore[n] = 0;
							start[n] =  0;
							finish[n] = threshold;
						}
					}

					for (int scoreIndex = start[0]; scoreIndex < finish[0]; scoreIndex++) {
						for (int scoreIndex2 = start[1]; scoreIndex2 < finish[1]; scoreIndex2++) {
							x[m][scoreIndex][scoreIndex2] += x[massIndex][scoreIndex - currScore[0]][scoreIndex2 - currScore[1]] * 0.05;

						}
					}
				}
			}
		}

		delete [] currScore;
		delete [] start;
		delete [] finish;
}


//-----------------------------------------------------------------------------
// Calculates the score of three spectra

void ExecGFNetwork::constructPSMMatrixForTriple(int offset, int threshold, int parentMass)
{
		int ROWS = 187;
		int COLS = threshold*threshold*threshold;
		DEBUG_MSG("ExecGFNetwork::computeMultiSpectraProbability_K3 init locals");

		dpTripleMatrix = new double *[ROWS] ;

		for( int r = 0 ; r < ROWS ; r++ )
			dpTripleMatrix[r] = new double[COLS];

		for (int r = 0; r < ROWS; r++)
			for (int c = 0; c < COLS; c++)
				dpTripleMatrix[r][c] = 0;

		dpTripleMatrix[0][offset * threshold * threshold + offset * threshold + offset] = 1;


		int * delta = new int[numberSpectra];
		int * start = new int[numberSpectra];
		int * finish = new int[numberSpectra];

		int threshold_squared = threshold * threshold;

		DEBUG_MSG("Start ExecGFNetwork::computeMultiSpectraProbability_K3 dynamic programming");

		for (int m = 1; m <= parentMass; m++)
		{
			int currentMassRow = m % ROWS;

			int jumpIndex = 0;
			for(int n = 0; n < numberSpectra; n++)
			{
				if(types[n] == FULL)
				{
					delta[n] = (m - 1 < (*scores)[indices[n]].size()) ? (int)(*scores)[indices[n]][m - 1][1] : 0;
				}
				else if (types[n] == PRE)
				{
					delta[n] = (m - 1 < (*scores)[indices[n]].size()) ? (int)(*scores)[indices[n]][m - 1][1] : 0;
					jumpIndex++;
				}
				else if (types[n] == POST)
				{
					delta[n] = ( (m - 1 - jumps[jumpIndex] < (*scores)[indices[n]].size()) && (m - 1 - jumps[jumpIndex] >= 0)) ?
					(int)(*scores)[indices[n]][m - 1 - jumps[jumpIndex] ][1] : 0;
					jumpIndex++;
				}
				else if (types[n] == SUB)
				{
					delta[n] = ((m - 1 - jumps[jumpIndex] < (*scores)[indices[n]].size()) && (m >= jumps[jumpIndex])) ?
					(int)(*scores)[indices[n]][m - 1 - jumps[jumpIndex] ][1] : 0;
					jumpIndex += 2;
				}
				start[n] = (delta[n] >= 0) ? delta[n] : -1*delta[n];
				finish[n] = (delta[n] >= 0) ? threshold: threshold + delta[n];
			}
			int subtract_component = (finish[1] - start[1])*threshold;

			for (int j = 0; j < AAcount; j++)
			{
				int aa = (int) (round(AAmasses[j]));
				if (m >= aa && checkBoundaries(m, aa))
				{
					int previousMassRow =  (m - aa)%ROWS;
					int o = start[0] * threshold_squared + start[1] * threshold + start[2];
					int adjOffset = o - (delta[0] * threshold_squared + delta[1] * threshold + delta[2]);

					for (int scoreIndex = start[0]; scoreIndex < finish[0]; scoreIndex++)
					{
						for (int scoreIndex2 = start[1]; scoreIndex2 < finish[1]; scoreIndex2++)
						{
							int o_inner = o;
							int adjOffset_inner = adjOffset;

							double * yc = dpTripleMatrix[currentMassRow] + o_inner;
							double * yp = dpTripleMatrix[previousMassRow] + adjOffset_inner;

							for (int scoreIndex3 = start[2]; scoreIndex3 < finish[2]; scoreIndex3++)
							{	*yc += *yp * 0.05;
								yp++;
								yc++;
							}
							o += threshold;
							adjOffset += threshold;

						}
						// substract off (finish[1] - start[1])*threshold
						o = o - subtract_component + threshold_squared;
						adjOffset = adjOffset - subtract_component + threshold_squared;

					}
				}
			}
		}

	DEBUG_MSG("ExecGFNetwork::computeMultiSpectraProbability_K3 finish dynamic programming");

	delete [] delta;
	delete [] start;
	delete [] finish;
}


//---------------------------------------------------------------------------------------------------
double inline ExecGFNetwork::scoreSpectraPair(int score1, int score2, int threshold, int parentMass){

		double retSpecProb = 0;


		//DEBUG_MSG(parentMass);
		//DEBUG_MSG("score 1: " << score1);
		//DEBUG_MSG("score 2: " << score2);

		for (int t0 = score1; t0 < threshold; t0++)	{
			for (int t1 = score2; t1 < threshold; t1++)	{
				retSpecProb += x[parentMass][t0][t1];
				DEBUG_VAR(x[parentMass][t0][t1]);
			}
		}

		return retSpecProb;

}

//---------------------------------------------------------------------------------------------------
double inline ExecGFNetwork::scoreSpectraTriple(int score1, int score2, int score3, int threshold, int parentMass){

	int ROWS = 187;
	int threshold_squared = threshold * threshold;

	double retSpecProb = 0;
	int currentMassRow = parentMass % ROWS;

	int t = score1* threshold_squared + score2 * threshold + score3;
	int subtract_component = (threshold - score2)*threshold;

	for (int t0 = score1; t0 < threshold; t0++)
	{
		for (int t1 = score2; t1 < threshold; t1++)
		{
			double * yc = dpTripleMatrix[currentMassRow] + t;

			for (int t2 = score3; t2 < threshold; t2++)
			{
				retSpecProb += *yc;
				yc++;

			}
			t += threshold;
		}
		t = t - subtract_component + threshold_squared;
	}


	DEBUG_MSG(retSpecProb);

}

//-----------------------------------------------------------------------------
// Calculates the score of two spectra
// finished.
// uses indices, types, jumps, numberSpectra, currentParentMass
double ExecGFNetwork::computeMultiSpectraProbability(int offset, int threshold)
{

		int rowSize = currentParentMass + 1;
		int colSize = threshold;
		int arraySize = rowSize * colSize * colSize;
		double * y = new double[arraySize];
		for (int i = 0; i < arraySize; i++)
			y[i] = 0;

		y[0 * colSize * colSize + offset * colSize + offset] = 1;

		int * delta = new int[numberSpectra];
		int * start = new int[numberSpectra];
		int * finish = new int[numberSpectra];

		int colSize_sqrd = colSize * colSize;

		for (int m = 1; m <=currentParentMass; m++) {

			int mOffset = m * colSize * colSize;

			int jumpIndex = 0;
			for(int n = 0; n < numberSpectra; n++)
			{
				if(types[n] == FULL)
				{
					delta[n] = (m - 1 < (*scores)[indices[n]].size()) ? (int)(*scores)[indices[n]][m - 1][1] : 0;
				}
				else if (types[n] == PRE)
				{
					delta[n] = (m - 1 < (*scores)[indices[n]].size()) ? (int)(*scores)[indices[n]][m - 1][1] : 0;
					jumpIndex++;
				}
				else if (types[n] == POST)
				{
											// m - jumps[jumpIndex]
					delta[n] = ( (m - 1 - jumps[jumpIndex] < (*scores)[indices[n]].size()) && (m - 1 - jumps[jumpIndex] >= 0)) ?
					(int)(*scores)[indices[n]][m - 1 - jumps[jumpIndex] ][1] : 0;
					jumpIndex++;
				}
				else if (types[n] == SUB)
				{
					delta[n] = ((m - 1 - jumps[jumpIndex] < (*scores)[indices[n]].size()) && (m >= jumps[jumpIndex])) ?
					(int)(*scores)[indices[n]][m - 1 - jumps[jumpIndex] ][1] : 0;
					jumpIndex += 2;
				}
				start[n] = (delta[n] >= 0) ? delta[n] : -1*delta[n];
				finish[n] = (delta[n] >= 0) ? threshold: threshold + delta[n];
			}

			for (int j = 0; j < AAcount; j++)
			{
				int aa = (int) (round(AAmasses[j]));
				int massIndex = m - aa;
				if (massIndex >= 0 && checkBoundaries(m, aa))
				{
					int massOffset =  massIndex * colSize_sqrd;
					for (int scoreIndex = start[0]; scoreIndex < finish[0]; scoreIndex++)
					{

						int o = mOffset + scoreIndex * colSize + start[1]; //(m,scoreIndex,start)
						int adjOffset = (scoreIndex-delta[0]) * colSize + massOffset + start[1] - delta[1]; //(masIndex,scoreIndex-delta,start-delta)

						for (int scoreIndex2 = start[1]; scoreIndex2 < finish[1]; scoreIndex2++) {

							y[o] += y[adjOffset] * 0.05;
							o++;
							adjOffset++;
						}

					}
				}
			}
		}


		double retSpecProb = 0;

		for (int t0 = peptideScores[0]; t0 < threshold; t0++)
		{
			int tStart= currentParentMass * colSize_sqrd + (t0 * colSize) + peptideScores[1];
			int tFinish = currentParentMass * colSize_sqrd + (t0 * colSize) + threshold;

			for (int t1 = tStart; t1 < tFinish; t1++)
			{
				retSpecProb += y[t1];
			}
		}


		DEBUG_MSG(retSpecProb);


		delete [] y;
		delete [] delta;
		delete [] start;
		delete [] finish;

		return retSpecProb;

}
//-----------------------------------------------------------------------------
// Calculates the score of three spectra

void ExecGFNetwork::computeMultiSpectraProbability_K3(int offset, int threshold)
{
		int ROWS = 187;
		int COLS = threshold*threshold*threshold;
		DEBUG_MSG("ExecGFNetwork::computeMultiSpectraProbability_K3 init locals");

		double **y = 0;
		y = new double *[ROWS] ;

		for( int r = 0 ; r < ROWS ; r++ )
			y[r] = new double[COLS];

		for (int r = 0; r < ROWS; r++)
			for (int c = 0; c < COLS; c++)
				y[r][c] = 0;
		y[0][offset * threshold * threshold + offset * threshold + offset] = 1;


		int * delta = new int[numberSpectra];
		int * start = new int[numberSpectra];
		int * finish = new int[numberSpectra];
		int * peptideScores = new int[numberSpectra];

		DEBUG_MSG("ExecGFNetwork::computeMultiSpectraProbability_K3 init 4d array");

		/*double ****array_4d;
		int dim1 = ROWS;
		int dim2 = threshold;
		int dim3 = threshold;
		int dim4 = threshold;

		array_4d = new double ***[dim1];
		for(int d1 = 0; d1 < dim1; d1++)
		{
		    array_4d[d1] = new double **[dim2];
		    for(int d2 = 0; d2 < dim2; d2++)
		    {
		        array_4d[d1][d2] = new double *[dim3];
		        for(int d3 = 0; d3 < dim3; d3++)
		        {
		            array_4d[d1][d2][d3] = new double[dim4];
		        }
		    }
		}
		for(int d1 = 0; d1 < dim1; d1++)
			for(int d2 = 0; d2 < dim2; d2++)
				for(int d3 = 0; d3 < dim3; d3++)
					for(int d4 = 0; d4 < dim4; d4++)
						array_4d[d1][d2][d3][d4] = 0;

		array_4d[0][offset][offset][offset] = 1;
*/
		int threshold_squared = threshold * threshold;

		DEBUG_MSG("Start ExecGFNetwork::computeMultiSpectraProbability_K3 dynamic programming");

		for (int m = 1; m <= currentParentMass; m++)
		{
			int currentMassRow = m % ROWS;

			for(int n = 0; n < numberSpectra; n++)
			{
				int jumpIndex = 0;
				if(types[n] == FULL)
				{
					delta[n] = (m - 1 < (*scores)[indices[n]].size()) ? (int)(*scores)[indices[n]][m - 1][1] : 0;
				}
				else if (types[n] == PRE)
				{
					delta[n] = (m - 1 < (*scores)[indices[n]].size()) ? (int)(*scores)[indices[n]][m - 1][1] : 0;
					jumpIndex++;
				}
				else if (types[n] == POST)
				{
					delta[n] = ( (m - 1 - jumps[jumpIndex] < (*scores)[indices[n]].size()) && (m - 1 - jumps[jumpIndex] >= 0)) ?
					(int)(*scores)[indices[n]][m - 1 - jumps[jumpIndex] ][1] : 0;
					jumpIndex++;
				}
				else if (types[n] == SUB)
				{
					delta[n] = ((m - 1 - jumps[jumpIndex] < (*scores)[indices[n]].size()) && (m >= jumps[jumpIndex])) ?
					(int)(*scores)[indices[n]][m - 1 - jumps[jumpIndex] ][1] : 0;
					jumpIndex += 2;
				}
				start[n] = (delta[n] >= 0) ? delta[n] : -1*delta[n];
				finish[n] = (delta[n] >= 0) ? threshold: threshold + delta[n];
			}
			int subtract_component = (finish[1] - start[1])*threshold;

			for (int j = 0; j < AAcount; j++)
			{
				int aa = (int) (round(AAmasses[j]));
				if (m >= aa && checkBoundaries(m, aa))
				{
					int previousMassRow =  (m - aa)%ROWS;
					int o = start[0] * threshold_squared + start[1] * threshold + start[2];
					int adjOffset = o - (delta[0] * threshold_squared + delta[1] * threshold + delta[2]);

					for (int scoreIndex = start[0]; scoreIndex < finish[0]; scoreIndex++)
					{
						for (int scoreIndex2 = start[1]; scoreIndex2 < finish[1]; scoreIndex2++)
						{
							//int o = scoreIndex * threshold * threshold + scoreIndex2 * threshold + start[2];
							//int adjOffset = (scoreIndex-delta[0]) * threshold * threshold + (scoreIndex2 - delta[1]) * threshold + start[2] - delta[2];
							int o_inner = o;
							int adjOffset_inner = adjOffset;

							double * yc = y[currentMassRow] + o_inner;
							double * yp = y[previousMassRow] + adjOffset_inner;

							for (int scoreIndex3 = start[2]; scoreIndex3 < finish[2]; scoreIndex3++)
							{
		//						array_4d[currentMassRow][scoreIndex][scoreIndex2][scoreIndex3] += 0.05*array_4d[previousMassRow][scoreIndex - delta[0]][scoreIndex2 - delta[1]][scoreIndex3 - delta[2]];
		//						y[currentMassRow][o_inner] += y[previousMassRow][adjOffset_inner] * 0.05;
		//						o_inner++;
		//						adjOffset_inner++;
								*yc += *yp * 0.05;
								yp++;
								yc++;
							}
							o += threshold;
							adjOffset += threshold;

						}
						// substract off (finish[1] - start[1])*threshold
						o = o - subtract_component + threshold_squared;
						adjOffset = adjOffset - subtract_component + threshold_squared;

					}
				}
			}
		}

		DEBUG_MSG("ExecGFNetwork::computeMultiSpectraProbability_K3 finish dynamic programming");

		double retSpecProb = 0;
		int currentMassRow = currentParentMass % ROWS;

		int t = peptideScores[0] * threshold_squared + peptideScores[1] * threshold + peptideScores[2];
		int subtract_component = (threshold - peptideScores[1])*threshold;

		for (int t0 = peptideScores[0]; t0 < threshold; t0++)
		{

			for (int t1 = peptideScores[1]; t1 < threshold; t1++)
			{
				double * yc = y[currentMassRow] + t;
			//	int tstart = t0 * threshold_squared + t1 * threshold + peptideScores[2];
			//	int tfinish = t0 * threshold_squared + t1 * threshold +  threshold;

				for (int t2 = peptideScores[2]; t2 < threshold; t2++)
				{
				//	retSpecProb += y[currentMassRow][t2];
					retSpecProb += *yc;
					yc++;

				}
				t += threshold;
			}
			t = t - subtract_component + threshold_squared;
		}


		DEBUG_MSG(retSpecProb);



	for( int r = 0 ; r < ROWS ; r++ )
			delete [] y[r];
	delete [] y;

	delete [] delta;
	delete [] start;
	delete [] finish;
}


}
