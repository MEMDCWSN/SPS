/*
 * ExecSpecNetworkEval.cpp
 *
 *  Created on: Jul 2011
 *      Author: cboucher@ucsd.edu
 */
// Header Includes
#include "ExecSpecNetworkEval.h"
#include "stdlib.h"
#include "FdrPeptide.h"
#include <cmath>

#include <ostream>

// Module Includes
#include "Logger.h"
#include "FileUtils.h"
#include "ParameterList.h"
#include "PeptideSpectrumMatchSet.h"
#include "SpectrumPairSet.h"
#include "SpectrumPair.h"
#include "ExecGFNetwork.h"
#include "spectrum.h"
#include "aminoacid.h"
#include "db_fasta.h"
#include <sstream>
#include <math.h>
#include <time.h>
#include <iostream>
#include <algorithm>
#include <vector>

using namespace specnets;

namespace specnets {


const unsigned int AAcount = 20;

const float AAmasses[] = {71.0371137870, 156.1011110260, 115.0269430310, 114.0429274460,
        				160.0306482000, 129.0425930950, 128.0585775100, 57.0214637230,
        				137.0589118610, 113.0840639790, 113.0840639790, 128.0949630160,
        				131.0404846050, 147.0684139150, 97.0527638510, 87.0320284090,
        				101.0476784730, 186.0793129520, 163.0633285370, 99.0684139150 };

// -------------------------------------------------------------------------

ExecSpecNetworkEval::ExecSpecNetworkEval(void) 
		: ownInput(true), filename(0x0), m_db(0x0), m_spectra(0x0), m_pairs(0x0)
{
	m_name = "ExecSpecNetworkEval";
	m_type = "ExecSpecNetworkEval";
}

// -------------------------------------------------------------------------

ExecSpecNetworkEval::ExecSpecNetworkEval(const ParameterList & inputParams) 
  : ExecBase(inputParams), 
    ownInput(true), filename(0x0), m_db(0x0), m_spectra(0x0), m_pairs(0x0)
{
	m_name = "ExecSpecNetworkEval";
	m_type = "ExecSpecNetworkEval";
}

// -------------------------------------------------------------------------

ExecSpecNetworkEval::~ExecSpecNetworkEval(void) {
	
	delete m_spectra;
	delete m_db;
}

// -------------------------------------------------------------------------

ExecBase * ExecSpecNetworkEval::clone(const ParameterList & inputParams) const {
	return new ExecSpecNetworkEval(inputParams);
}

// -------------------------------------------------------------------------

bool ExecSpecNetworkEval::saveInputData(std::vector<std::string> & filenames) {
	return true;
}

//-----------------------------------------------------------------------------

bool ExecSpecNetworkEval::loadOutputData(void) {
	return true;
}

bool ExecSpecNetworkEval::saveOutputData(void) {
	return true;
}

// -------------------------------------------------------------------------

vector<ExecBase *> const & ExecSpecNetworkEval::split(int numSplit) {
	m_subModules.resize(0);
	return m_subModules;
}

// -------------------------------------------------------------------------

bool ExecSpecNetworkEval::merge(void) {
	return false;
}

// -------------------------------------------------------------------------

bool ExecSpecNetworkEval::validateParams(std::string & error) {
	return true;
}


//-----------------------------------------------------------------------------

vector<int> ExecSpecNetworkEval::processShifts(float shift, int parentMass1, int parentMass2) {

		
	int s1 = 0;
	int s2 = 0; 
	int s3 = 0;
	
	vector<int> retShifts;
	retShifts.resize(4);	
	int jumpsIndex = 0;
	
	for(int i = 0; i < dimension*2; i++)
		spectrumJumps[i] = 0;
	
	int int_shift = (int)round(shift* 0.9995);
	
	if (int_shift == 0) 
	{
		s1 = 0;
		s2 = (parentMass1 < parentMass2) ? parentMass1 : parentMass2;
		s3 = (parentMass1 < parentMass2) ? parentMass2 - s2 : parentMass1 - s2;
		currentParentMass = (parentMass1 >= parentMass2) ? parentMass1 : parentMass2;
		spectrumTypes[0] = (parentMass1 == currentParentMass) ? FULL : PRE;
		spectrumTypes[1] = (parentMass2 == currentParentMass) ? FULL : PRE;    
	
		if(spectrumTypes[0] == PRE || spectrumTypes[1]==PRE) 
		{
			spectrumJumps[jumpsIndex] = s2; 
			jumpsIndex++;
		}
				
	} 
	else if (int_shift > 0) 
	{

		s1 = (int) round(shift * 0.9995);
		s2 = (parentMass1 <= parentMass2 + s1) ? parentMass1 - s1 : parentMass2;
		s3 = (parentMass1 <= parentMass2 + s1) ? parentMass2 - s2 : parentMass1 - s1 - s2;
		currentParentMass = (parentMass1 >= parentMass2 + s1) ? parentMass1 : (parentMass2 + s1);
		spectrumTypes[0] = (parentMass1 < parentMass2 + s1) ? PRE : FULL;
		spectrumTypes[1] = (parentMass1 <= parentMass2 + s1) ? POST : SUB;
		
		if(spectrumTypes[0] == PRE) 
		{
			spectrumJumps[jumpsIndex] = parentMass1; 
			jumpsIndex++;
		}
		if(spectrumTypes[1] == POST) 
		{
			spectrumJumps[jumpsIndex] = s1; 
			jumpsIndex++;
		}
		else if(spectrumTypes[1] == SUB) 
		{
			spectrumJumps[jumpsIndex] = s1;
			jumpsIndex++;
			spectrumJumps[jumpsIndex] = s1 + parentMass2;
			jumpsIndex++;
		}
	}						
	else
	{
		
		s1 = (int) round(-1*shift * 0.9995);
		s2 = (parentMass1 + s1 >= parentMass2) ? parentMass2 - s1 : parentMass1;
		s3 = (parentMass1 + s1 >= parentMass2) ? parentMass1 - s2 : parentMass2 - s2 - s1;
	
		spectrumTypes[0] = (parentMass1 + s1 < parentMass2) ? SUB : POST;
		spectrumTypes[1] = (parentMass1 + s1 <= parentMass2) ? FULL : PRE;
		
		if(spectrumTypes[0] == POST) 
		{
			spectrumJumps[jumpsIndex] = s1; 
			jumpsIndex++;
		}
		else if(spectrumTypes[0] == SUB) 
		{
			spectrumJumps[jumpsIndex] = s1;
			jumpsIndex++;
			spectrumJumps[jumpsIndex] = s1 + parentMass1;
			jumpsIndex++;
		}
		if(spectrumTypes[1] == PRE) 
		{
			spectrumJumps[jumpsIndex] = parentMass2; 
			jumpsIndex++;
		}
						
	}
	
	retShifts[0] = s1;
	retShifts[1] = s2;
	retShifts[2] = s3;
	retShifts[3] = jumpsIndex;
		

	
	return retShifts;
}

//-------------------------------------------------------------------------//
void ExecSpecNetworkEval::outputZeroFDRPSMSet(PeptideSpectrumMatchSet psmSet, PeptideSpectrumMatchSet outputPsmSet)
{
	// Output all PSM information
	for(int i = 0; i < psmSet.size(); i++)
	{
		//PeptideSpectrumMatchNetwork * p = (PeptideSpectrumMatchNetwork *)psmSet[i].get();
		//std::tr1::shared_ptr<PeptideSpectrumMatch> psmPointer(new PeptideSpectrumMatchNetwork); //I AM A PSM POINTER!
		std::tr1::shared_ptr<PeptideSpectrumMatchNetwork> p = std::tr1::dynamic_pointer_cast<PeptideSpectrumMatchNetwork>(psmSet[i]); //I AM A PSMNETWORK POINTER

		if(p)
		{
			if(p->m_isTuple)
			{
				cout << p->m_spectrumIndices[0] << ", " << p->m_spectrumIndices[1] << "\t" << p->m_annotation << "\t" << p->m_score << "\t";
				cout << p->m_parentMasses[0] << ", " << p->m_parentMasses[1] << "\t";
				cout << p->m_partialAnnotations[0] << ", " << p->m_partialAnnotations[1] << "\t";
				cout << p->m_partialSpecProbabilities[0] << ", " << p->m_partialSpecProbabilities[1] << "\t";
				cout << p->m_partialPeptideScores[0] << ", " << p->m_partialPeptideScores[1] << "\t";
				cout << "0\n";
			}
			else
			{
				cout << p->m_spectrumIndices[0] << "\t" << p->m_annotation << "\t" << p->m_score << "\t";
				cout << p->m_parentMasses[0] << "\t";
				cout << ",\t";
				cout << ",\t";
				cout << p->m_individualScore << "\t";
				cout <<  "0 \n";
			}
		}
	}
}


//-------------------------------------------------------------------------//
void ExecSpecNetworkEval::outputPSMSet(PeptideSpectrumMatchSet psmSet, PeptideSpectrumMatchSet outputPsmSet) 
{
	// Output all PSM information
	for(int i = 0; i < psmSet.size(); i++)
	{
		//PeptideSpectrumMatchNetwork * p = (PeptideSpectrumMatchNetwork *)psmSet[i].get();
		//std::tr1::shared_ptr<PeptideSpectrumMatch> psmPointer(new PeptideSpectrumMatchNetwork); //I AM A PSM POINTER!
		std::tr1::shared_ptr<PeptideSpectrumMatchNetwork> p = std::tr1::dynamic_pointer_cast<PeptideSpectrumMatchNetwork>(psmSet[i]); //I AM A PSMNETWORK POINTER

		if(!p)
		{
			if(p->m_isTuple)
			{
				cout << p->m_spectrumIndices[0] << ", " << p->m_spectrumIndices[1] << "\t" << p->m_annotation << "\t" << p->m_score << "\t";
				cout << p->m_parentMasses[0] << ", " << p->m_parentMasses[1] << "\t";
				cout << p->m_partialAnnotations[0] << ", " << p->m_partialAnnotations[1] << "\t";
				cout << p->m_partialSpecProbabilities[0] << ", " << p->m_partialSpecProbabilities[1] << "\t";
				cout << p->m_partialPeptideScores[0] << ", " << p->m_partialPeptideScores[1] << "\t";
				cout << outputPsmSet[i]->m_pValue << "\n";
			}
			else
			{
				cout << p->m_spectrumIndices[0] << "\t" << p->m_annotation << "\t" << p->m_score << "\t";
				cout << p->m_parentMasses[0] << "\t";
				cout << ",\t";
				cout << ",\t";
				cout << p->m_individualScore << "\t";
				cout << outputPsmSet[i]->m_pValue << "\n";
			}
		}
	}
}
//-------------------------------------------------------------------------//
std::tr1::shared_ptr<PeptideSpectrumMatchNetwork> ExecSpecNetworkEval::findInPSMSet(PeptideSpectrumMatchSet psmSet, int index)
{
	// Output all PSM information
	for(int i = 0; i < psmSet.size(); i++)
	{

		psmPtr psm = psmSet[i];
		if(psm)
		{
			std::tr1::shared_ptr<PeptideSpectrumMatchNetwork> p = std::tr1::static_pointer_cast<PeptideSpectrumMatchNetwork>(psm); //I AM A PSMNETWORK POINTER

			//PeptideSpectrumMatchNetwork * p = (PeptideSpectrumMatchNetwork *)psmSet[i].get();

			if(p->m_spectrumIndices[0] == index || p->m_spectrumIndices[1] == index)
				return p;

		}
	}
	return std::tr1::shared_ptr<PeptideSpectrumMatchNetwork>();

}

//-------------------------------------------------------------------------//
bool ExecSpecNetworkEval::annotateByPairs(void) 
{
	dimension = 2;
	spectrumIndices.resize(dimension);
	spectrumTypes.resize(dimension);
	spectrumJumps.resize(dimension*2);
		
	vector<string> emptyPartialAnnotations;
	vector<int> emptyPartialPeptideScores;
	vector<double> emptyPartialSpecProbabilities;
	int parentMass = 0;
	int numberPairs = 0;
	PeptideSpectrumMatchSet psmSet;
			
	for(int i = 0; i < m_pairs->size(); i++)
	{	

		spectrumIndices[0] = (*m_pairs)[i].spec1;
		spectrumIndices[1] = (*m_pairs)[i].spec2;	
		
		cout << spectrumIndices[0] << ", " << spectrumIndices[1] << "\n";

		intVecItr = find (annotatedSpectra.begin(), annotatedSpectra.end(), spectrumIndices[0]); 
		if(intVecItr != annotatedSpectra.end())
			continue;

		intVecItr = find (annotatedSpectra.begin(), annotatedSpectra.end(), spectrumIndices[1]); 
		if(intVecItr != annotatedSpectra.end())
			continue;
		
		Spectrum sp1 = (*m_spectra)[spectrumIndices[0]];
		Spectrum sp2 = (*m_spectra)[spectrumIndices[1]];
		
		int parentMass1 = (int) round(sp1.parentMass * 0.9995) - 19;
		int parentMass2 = (int) round(sp2.parentMass * 0.9995) - 19;
		
		double optimalScore = 99;
		string optimalAnnotation = "";
		int optimalPeptideIndex = -1;
		float optimalShift = 0;
				
		/*
		 * Try annotation with the first shift
		 */
		vector<int> shifts = processShifts((*m_pairs)[i].shift1, parentMass1, parentMass2);
		vector<string> peptides = m_db->findPeptides(shifts[0], shifts[1], shifts[2]);	
				
		int shift1 = (int)round((*m_pairs)[i].shift1* 0.9995);
		int shift2 = (int)round((*m_pairs)[i].shift2* 0.9995);			
		if (shift1 >= 0)
		{
			parentMass = ((parentMass2 + shift1) > parentMass1) ? parentMass2 + shift1 : parentMass1;
		}
		else
		{
			shift1 = (int) round(-1*shift1 * 0.9995);
			parentMass = (parentMass2 >= parentMass1 + shift1) ? parentMass2 : parentMass1 + shift1;
		}

		
		if(peptides.size() > 0 )
		{

			annotatedSpectra.push_back(spectrumIndices[0]);
			annotatedSpectra.push_back(spectrumIndices[1]);
					
			vector<double> specPeptideScores = gfNetwork->scorePSMPairs(spectrumIndices, spectrumTypes, peptides, spectrumJumps, parentMass);
		
			for(int p = 0; p < peptides.size(); p++)
			{
		//		DEBUG_MSG(peptides[p] << " " << specPeptideScores[p]);

				if(specPeptideScores[p] < optimalScore && specPeptideScores[p] != 0)
				{
					optimalScore = specPeptideScores[p];
					optimalAnnotation = peptides[p];
					optimalShift = (*m_pairs)[i].shift1;
					optimalPeptideIndex = m_db->peptideIndices[p];
				}
			}
		}
		
		/*
		 * Try annotation with the second shift
		 */
		shifts = processShifts((*m_pairs)[i].shift2, parentMass1, parentMass2);
		if (shift2 >= 0)
		{
			parentMass = (parentMass2 + shift2 >= parentMass1) ? parentMass2 + shift2: parentMass1;
		}
		else
		{
			parentMass = (parentMass2 >= parentMass1 + shift2) ? parentMass2 : parentMass1 + shift2;
		}
		

		peptides = m_db->findPeptides(shifts[0], shifts[1], shifts[2]);
		if(peptides.size() > 0)
		{
			vector<double> specPeptideScores = gfNetwork->scorePSMPairs(spectrumIndices, spectrumTypes, peptides, spectrumJumps, parentMass);
			for(int p = 0; p < peptides.size(); p++)
			{

	//			DEBUG_MSG(peptides[p] << " " << specPeptideScores[p]);
				
				if(specPeptideScores[p] < optimalScore && specPeptideScores[p] != 0)
				{
					optimalScore = specPeptideScores[p];
					optimalAnnotation = peptides[p];
					optimalShift = (*m_pairs)[i].shift2;
					optimalPeptideIndex = m_db->peptideIndices[p];
				
				}
			}
		}
		if(optimalScore < 99)
		{
			annotatedSpectra.push_back(spectrumIndices[0]);
			annotatedSpectra.push_back(spectrumIndices[1]);
			vector<int> shifts = processShifts(optimalShift, parentMass1, parentMass2);
			vector<string> peptides = m_db->findPeptides(shifts[0], shifts[1], shifts[2]);	
					
			int shift = (int)round(optimalShift* 0.9995);
			if (shift >= 0)
			{
				parentMass = ((parentMass2 + shift) > parentMass1) ? parentMass2 + shift : parentMass1;
			}
			else
			{
				shift = (int) round(-1*shift * 0.9995);
				parentMass = (parentMass2 >= parentMass1 + shift) ? parentMass2 : parentMass1 + shift;
			}
				
			cout << spectrumIndices[0] << ", " << spectrumIndices[1] << "\t" << optimalAnnotation << "\t" << optimalScore << "\t";
			cout << parentMass1 << ", " << parentMass2 << "\t";
						
			gfNetwork->outputForPairs(spectrumIndices, spectrumTypes, optimalAnnotation, spectrumJumps, parentMass);
			
			// Calulate FDR
			psmNetworkPtr psm(new PeptideSpectrumMatchNetwork);
			// Spectra
			psm->m_spectra.push_back(&sp1);
			psm->m_spectra.push_back(&sp2);
			psm->m_spectrumIndices.push_back(spectrumIndices[0]);
			psm->m_spectrumIndices.push_back(spectrumIndices[1]);
			psm->m_isDecoy = m_db->checkIfDecoy(optimalPeptideIndex);
			psm->m_isTuple = true;
		    psm->m_parentMasses.resize(2);
		    psm->m_parentMasses[0] = parentMass1;
		    psm->m_parentMasses[1] = parentMass2;
			//SpecProbability
			psm->m_score = optimalScore;
			psm->m_annotation = optimalAnnotation;
			psm->m_partialAnnotations = gfNetwork->partialAnnotations;
			psm->m_partialPeptideScores = gfNetwork->partialPeptideScores;
			psm->m_partialSpecProbabilities = gfNetwork->partialSpecProbabilities;		
			psmSet.push_back((psmPtr&)psm);
		}
	}


	for(int i = 0; i < m_pairs->size(); i++)
	{

		spectrumIndices[0] = (*m_pairs)[i].spec1;
		Spectrum sp1 = (*m_spectra)[spectrumIndices[0]];
		int parentMass1 = (int) round(sp1.parentMass * 0.9995) - 19;

		double optimalScore = 99;
		bool updatePSMSet = false;
		string optimalAnnotation = "";
		int optimalPeptideIndex = -1;
		int index = spectrumIndices[0];

		intVecItr = find (annotatedSpectra.begin(), annotatedSpectra.end(), spectrumIndices[0]);
		if(intVecItr != annotatedSpectra.end())
		{
			// find the element in the psmSet
			//std::tr1::shared_ptr<PeptideSpectrumMatchNetwork> p = findInPSMSet(psmSet, spectrumIndices[0]); // not sure about this either
			//psmPtr p = findInPSMSet(psmSet, spectrumIndices[0]); // not sure about this either

			bool found = false;

			for(int j = 0; j < psmSet.m_psmSet.size(); j++)
			{
				DEBUG_VAR(j);

				psmPtr psm = psmSet.m_psmSet[j];
				if(psm != NULL)
				{
					std::tr1::shared_ptr<PeptideSpectrumMatchNetwork> ptr = std::tr1::static_pointer_cast<PeptideSpectrumMatchNetwork>(psm); //I AM A PSMNETWORK POINTER

					//PeptideSpectrumMatchNetwork * ptr = (PeptideSpectrumMatchNetwork *)psm.get();
					if( ptr->m_spectrumIndices[0] == index || ptr->m_spectrumIndices[1] == index)
					{
						optimalScore = psm->m_score;
						DEBUG_VAR( psm->m_score);
						found = true;
					}
				}
			}

			/*if(found){
				vector<string> peptides = m_db->findPeptideForSingleMass(parentMass1);

				for(int p = 0; p < peptides.size(); p++)
				{
					vector<char> sequence = vector<char> (peptides[p].begin(), peptides[p].end());
					vector<float> masses;
					getMasses(sequence, masses);
					float specPeptideScore = gfNetwork->computeSpectraProbability(spectrumIndices[0], peptides[p], 200, 450 + 200 + 1);

					if(specPeptideScore < optimalScore && specPeptideScore != 0)
					{
						optimalScore = specPeptideScore;
						optimalAnnotation = peptides[p];
						optimalPeptideIndex = m_db->peptideIndices[p];
						updatePSMSet = true;
					}
				}
				if(updatePSMSet)
				{
					cout << spectrumIndices[0] << "\t" << optimalAnnotation << "\t" << optimalScore << "\t" << parentMass1 << "\n";
			//		psmSet.removePsmSetItem(p);

				psmNetworkPtr psm(new PeptideSpectrumMatchNetwork);
					psm->m_spectra.push_back(&sp1);
					psm->m_spectra.push_back(0);
					psm->m_spectrumIndices.push_back(spectrumIndices[0]);
					psm->m_spectrumIndices.push_back(0);
					psm->m_isDecoy = m_db->checkIfDecoy(optimalPeptideIndex);
					psm->m_isTuple = false;
					psm->m_individualScore = gfNetwork->scorePeptide(optimalAnnotation, spectrumIndices[0]);
					psm->m_score = optimalScore;
					psm->m_annotation = optimalAnnotation;
					psm->m_partialAnnotations = emptyPartialAnnotations;
					psm->m_partialPeptideScores = emptyPartialPeptideScores;
					psm->m_partialSpecProbabilities = emptyPartialSpecProbabilities;
					psm->m_parentMasses.resize(1);
					psm->m_parentMasses[0] = parentMass1;
					psmSet.push_back((psmPtr&)psm);
				}
			}else
			{
				DEBUG_MSG("Pointer not found!");
			}*/
		}
	}

	/*PeptideSpectrumMatchSet outputPsmSet;

	if (!FdrPeptide::calculatePValues(psmSet, outputPsmSet)) {
		outputZeroFDRPSMSet(psmSet, outputPsmSet);
	}
	else {
		outputPSMSet(psmSet, outputPsmSet);
	}*/

	return true;
}
//-----------------------------------------------------------------------------

void ExecSpecNetworkEval::outputType(spectra_type_t t) 
{
	if(t == PRE)
	{
		DEBUG_MSG("PRE");
	}
	else if (t == POST)
	{
		DEBUG_MSG("POST");
	} 
	else if (t == SUB)
	{
		DEBUG_MSG("SUB");
	}
	else 		
	{
		DEBUG_MSG("FULL");
	}
}
//-----------------------------------------------------------------------------

void ExecSpecNetworkEval::annotateSingleTuple(int shift_i_int, int shift_j_int, int parentMass1, int parentMass2, int parentMass3) 
{
	spectrumJumps.resize(6);
	vector<int> tempJumps;
	
	tempJumps.resize(6);
	tempJumps[0] = 0;
	tempJumps[1] = parentMass1;
	tempJumps[2] = shift_i_int;
	tempJumps[3] = shift_i_int + parentMass2;
	tempJumps[4] = shift_j_int;
	tempJumps[5] = shift_j_int + parentMass3;
		
	int M = min( min(0, shift_i_int), shift_j_int);
	int parentMass = 0;
				
	for(int tempIndex = 0; tempIndex < tempJumps.size(); tempIndex++)
	{	
		tempJumps[tempIndex] = tempJumps[tempIndex] + abs(M);
		if(parentMass < tempJumps[tempIndex])
			parentMass = tempJumps[tempIndex];
	}
		
	int t = 0;
	int jumpsIndex = 0;
	for(int tempIndex = 0; tempIndex < tempJumps.size(); tempIndex++)
	{	
		if(tempJumps[tempIndex] == 0)
		{
			tempIndex++;
			if(tempJumps[tempIndex] == parentMass)
			{
				spectrumTypes[t] = FULL;
				t++;
			}
			else
			{
				spectrumTypes[t] = PRE;
				spectrumJumps[jumpsIndex] = tempJumps[tempIndex];
				jumpsIndex++;
				t++;
			}		
		} else
		{
			tempIndex++;
			if(tempJumps[tempIndex] == parentMass)
			{
				spectrumTypes[t] = POST;
				spectrumJumps[jumpsIndex] = tempJumps[tempIndex - 1];
				jumpsIndex++;
				t++;
			}
			else
			{
				spectrumTypes[t] = SUB;
				spectrumJumps[jumpsIndex] = tempJumps[tempIndex - 1];
				jumpsIndex++;
				spectrumJumps[jumpsIndex] = tempJumps[tempIndex];
				jumpsIndex++;
				t++;
			}
		}
	}		
		

	vector<string> peptides;
		
	vector<int> temp_shifts = tempJumps;
	std::sort(temp_shifts.begin(), temp_shifts.end());
	temp_shifts.erase(unique(temp_shifts.begin(), temp_shifts.end()), temp_shifts.end());
			
	if(temp_shifts.size() <= 4)
	{
		vector<int> shifts;
		int prev = 0;		
		for(int ii = 1; ii < temp_shifts.size(); ii++)
		{
			shifts.push_back(temp_shifts[ii] - prev);
			prev = temp_shifts[ii];
		}
		for(int ii = shifts.size(); ii <= 3; ii++)				
			shifts.push_back(0);
			
		DEBUG_MSG("shifts_spec1_spec2: " << shifts[0] << " " << shifts[1] << " " << shifts[2] ); 
		peptides = m_db->findPeptides(shifts[0], shifts[1], shifts[2]);	
		
	} else{	
			
		vector<int> shifts;
		int prev = 0;		
		for(int ii = 1; ii < 3; ii++)
		{
			shifts.push_back(temp_shifts[ii] - prev);
			prev = temp_shifts[ii];
		}
		shifts.push_back(temp_shifts[temp_shifts.size() - 1] - prev);
		prev = 0;
		for(int ii = 2; ii < temp_shifts.size(); ii++)
		{
			shifts.push_back(temp_shifts[ii] - prev);
			prev = temp_shifts[ii];
		}			
			
		for(int ii = shifts.size(); ii <= 6; ii++)
			shifts.push_back(0);
		
		DEBUG_MSG("shifts_spec1_spec2: " << shifts[0] << " " << shifts[1] << " " << shifts[2] ); 
		DEBUG_MSG("shifts_spec1_spec3: " <<  shifts[3] << " " << shifts[4] << " " << shifts[5] ); 
		peptides = m_db->findPeptidesForTriple(shifts[0], shifts[1], shifts[2], shifts[3], shifts[4], shifts[5]);
	}


	DEBUG_MSG("Parent Mass 1: " << parentMass1 << " Parent Mass 2: " << parentMass2 << " parentMass3: " << parentMass3);			
	DEBUG_MSG("Shift 1: " << shift_i_int << " Shift 2: " << shift_j_int );			
	DEBUG_MSG("Parent Mass: " << parentMass);		
	DEBUG_MSG(tempJumps[0] << " " << tempJumps[1] << " " << tempJumps[2] << " " << tempJumps[3] << " " << tempJumps[4] << " " << tempJumps[5]);
	DEBUG_MSG(spectrumJumps[0] << " " << spectrumJumps[1] << " " << spectrumJumps[2] << " " << spectrumJumps[3] << " " << spectrumJumps[4] << " " << spectrumJumps[5]);
					
	outputType(spectrumTypes[0]);
	outputType(spectrumTypes[1]);
	outputType(spectrumTypes[2]);

	if(peptides.size() > 0 )
	{
		vector<double> specPeptideScores = gfNetwork->scorePSMTriple(spectrumIndices, spectrumTypes, peptides, spectrumJumps, parentMass);
			
		for(int p = 0; p < peptides.size(); p++)
		{
			DEBUG_MSG(peptides[p] << " " << specPeptideScores[p]);
					
			if(specPeptideScores[p] < optimalTupleScore && specPeptideScores[p] != 0)
			{
				optimalTupleScore = specPeptideScores[p];
				optimalTupleAnnotation = peptides[p];
				gfNetwork->outputPSMTriple(spectrumIndices, spectrumTypes, peptides[p], spectrumJumps, parentMass);
				optPartialAnnotations = gfNetwork->partialAnnotations;
				optPartialPeptideScores = gfNetwork->partialPeptideScores;
				optPartialSpecProbabilities = gfNetwork->partialSpecProbabilities;	
				optimalTuplePeptideIndex = m_db->peptideIndices[p];	
			}
		}
			
	}
}
//-----------------------------------------------------------------------------

bool ExecSpecNetworkEval::annotateByTuples(void) 
{
	dimension = 3;
	spectrumIndices.resize(dimension);
	spectrumTypes.resize(dimension);
	spectrumJumps.resize(dimension*2);
		
	int parentMass = 0;
	int i, j;
	
	PeptideSpectrumMatchSet psmSet;
	
	for(i = 0; i < m_pairs->size(); i++)
	{
		spectrumIndices[0] = (*m_pairs)[i].spec1;
		spectrumIndices[1] = (*m_pairs)[i].spec2;		
			
		// Determine if the spectra has been previously annotated.
		intVecItr = find (annotatedSpectra.begin(), annotatedSpectra.end(), spectrumIndices[0]); 
		if(intVecItr != annotatedSpectra.end())
			continue;

		intVecItr = find (annotatedSpectra.begin(), annotatedSpectra.end(), spectrumIndices[1]); 
		if(intVecItr != annotatedSpectra.end())
			continue;
		
		// look for a path of length two
		bool foundPath = false;
		
		for(j = i + 1; j < m_pairs->size(); j++)
		{
			if(	(*m_pairs)[i].spec1 == (*m_pairs)[j].spec1 )
			{
				spectrumIndices[2] = (*m_pairs)[j].spec2;
				foundPath = true;
				continue;
			}
		}
		
		if(foundPath)
		{
			// Otherwise, a path of length 2 has been found.
			optimalTupleScore = 99;
			optimalTupleAnnotation = "";
			optimalTuplePeptideIndex = -1;
			if(!optPartialAnnotations.empty())
				optPartialAnnotations.clear();
			if(!optPartialPeptideScores.empty())
				optPartialPeptideScores.clear();
			if(!optPartialSpecProbabilities.empty())
				optPartialSpecProbabilities.clear();
				
			Spectrum sp1 = (*m_spectra)[spectrumIndices[0]];
			Spectrum sp2 = (*m_spectra)[spectrumIndices[1]];
			Spectrum sp3 = (*m_spectra)[spectrumIndices[2]];
			int parentMass1 = (int) round(sp1.parentMass * 0.9995) - 19;
			int parentMass2 = (int) round(sp2.parentMass * 0.9995) - 19;
			int parentMass3 = (int) round(sp3.parentMass * 0.9995) - 19;
				
			float shift_i1 = (*m_pairs)[i].shift1;
			int shift_i_int1 = (int)round(shift_i1* 0.9995);
		
			float shift_i2 = (*m_pairs)[i].shift2;
			int shift_i_int2 = (int)round(shift_i2* 0.9995);
		
			float shift_j1 = (*m_pairs)[j].shift1;
			int shift_j_int1 = (int)round(shift_j1* 0.9995);
				
			float shift_j2 = (*m_pairs)[j].shift2;
			int shift_j_int2 = (int)round(shift_j2* 0.9995);
						
			annotateSingleTuple(shift_i_int1, shift_j_int1, parentMass1, parentMass2, parentMass3);
			annotateSingleTuple(shift_i_int1, shift_j_int2, parentMass1, parentMass2, parentMass3);
			annotateSingleTuple(shift_i_int2, shift_j_int1, parentMass1, parentMass2, parentMass3);
			annotateSingleTuple(shift_i_int2, shift_j_int2, parentMass1, parentMass2, parentMass3);
	
			if(optimalTupleScore < 99)
			{
				annotatedSpectra.push_back(spectrumIndices[0]);
				annotatedSpectra.push_back(spectrumIndices[1]);
				annotatedSpectra.push_back(spectrumIndices[2]);
				cout << spectrumIndices[0] << ", " << spectrumIndices[1] << ", " << spectrumIndices[2] << "\t" << optimalTupleAnnotation << "\t" << optimalTupleScore << "\t";
				cout << parentMass1 << ", " << parentMass2 << ", " << parentMass3 << "\t";
				gfNetwork->outputPSMTriple(spectrumIndices, spectrumTypes, optimalTupleAnnotation, spectrumJumps, parentMass);

			}
			/*
				// Calulate FDR
				psmNetworkPtr psm(new PeptideSpectrumMatchNetwork);
					
				// Spectra
				psm->m_spectra.push_back(&sp1);
				psm->m_spectra.push_back(&sp2);
				psm->m_spectra.push_back(&sp3);
				psm->m_spectrumIndices.push_back(spectrumIndices[0]);
				psm->m_spectrumIndices.push_back(spectrumIndices[1]);
				psm->m_spectrumIndices.push_back(spectrumIndices[2]);

				psm->m_isDecoy = m_db->checkIfDecoy(optimalTuplePeptideIndex);
						
				//SpecProbability
				psm->m_score = optimalScore;
				psm->m_annotation = optimalAnnotation;
				psm->m_partialAnnotations = optPartialAnnotations;
				psm->m_partialPeptideScores = optPartialPeptideScores;
				psm->m_partialSpecProbabilities = optPartialSpecProbabilities;

				psmSet.push_back((psmPtr&)psm);
			}*/
		}
	}
	return true;
}
//-----------------------------------------------------------------------------

bool ExecSpecNetworkEval::invoke(void) {

	// Sort the spectrum pairs according to the sum of the alignment scores
	m_pairs->sort_pairs();

	DEBUG_MSG("Adding reversed...");
	m_db->addDecoyReversed();
	DEBUG_MSG("\t\tFinished adding reversed...");

	// Construct the Index for the database
	m_db->populateIndex(false);
	DEBUG_MSG("\t\tFinished populating...");

	//108, 109	DKKIVPR	4.46875e-08	836, 836	DKKIVPR, DKKIVPR	1.875e-08, 9.08829e-06	62, -15
/*	AAJumps myjumps(1);
	string pep = "SKAEAESLYQSK";
	float mass_float = myjumps.getPeptideMass(pep);
	int currentParentMass = (int) round(mass_float*0.9995);

	spectrumIndices.resize(2);
	spectrumTypes.resize(2);
	spectrumJumps.resize(4);

	vector<string> peptides;
	peptides.push_back(pep);

	spectrumIndices[0] = 450;
	spectrumIndices[1] = 703;
	spectrumTypes[0] = POST;
	spectrumTypes[1] = FULL;
	spectrumJumps[0] = 215;
	spectrumJumps[1] = 0;
	spectrumJumps[2] = 0;
	spectrumJumps[3] = 0;

	//clock_t t1,t2;
	//t1=clock();
	vector<double> specPeptideScores = gfNetwork->scorePSMPairs(spectrumIndices, spectrumTypes, peptides, spectrumJumps, currentParentMass);
	DEBUG_MSG(specPeptideScores[0]);
	//t2=clock();
	//double seconds = ((float)t2-(float)t1) / CLOCKS_PER_SEC;
	//DEBUG_MSG("seconds: " << seconds);

		spectrumTypes[0] = PRE;
		spectrumTypes[1] = FULL;
		spectrumJumps[0] = 1106;
		spectrumJumps[1] = 0;
		spectrumJumps[2] = 0;
		spectrumJumps[3] = 0;

		vector<double> specPeptideScores2 = gfNetwork->scorePSMPairs(spectrumIndices, spectrumTypes, peptides, spectrumJumps, currentParentMass);
		DEBUG_MSG(specPeptideScores2[0]);*/


	return annotateByPairs();
	return annotateByTuples();


}
//-----------------------------------------------------------------------------

bool ExecSpecNetworkEval::loadInputData(void) {
	ownInput = false;


	m_spectra = new SpecSet;
	m_db = new DB_fasta;
	m_pairs = new SpectrumPairSet;
	m_pairs->resize(0);
	m_spectra->resize(0);


	if (m_params.exists("PRM_SPECTRA")) {
		if (!m_spectra->LoadSpecSet_mgf(
				m_params.getValue("PRM_SPECTRA").c_str())) {
			ERROR_MSG("Could not load " << m_params.getValue("PRM_SPECTRA"));
			return false;
		}
	}
	if (m_spectra->size() == 0) {
		ERROR_MSG("Could not load PRM Spectra. Filename: [" << m_params.getValue("PRM_SPECTRA") << "]") ;
		return false;
	}


	// Load in spectra pairs
	if (m_params.exists("SPEC_PAIRS")) {
			if (!m_pairs->loadFromBinaryFile(
					m_params.getValue("SPEC_PAIRS").c_str())) {
				ERROR_MSG("Could not load " << m_params.getValue("SPEC_PAIRS"));
				return false;
			}
		}
		if (m_pairs->size() == 0) {
			ERROR_MSG("Could not load spectrum pairs. Filename: [" << m_params.getValue("SPEC_PAIRS") << "]") ;
			return false;
		}

		DEBUG_MSG("Loading spectral pairs complete. Num pairs [" << m_pairs->size() << "]");


	// load in database
    if (!m_params.exists("INPUT_FASTA")) {
			ERROR_MSG("Parameters are incomplete. INPUT_FASTA is missing.");
      return false;
    } else if (m_db->Load(m_params.getValue("INPUT_FASTA").c_str())<=0) {
			ERROR_MSG("Error reading database sequences from "<<m_params.getValue("INPUT_FASTA"));
			return false;
		}


//	DEBUG_MSG("Loading peptides complete. Num peptides[" << m_db->size() << "]");

	gfNetwork = new ExecGFNetwork(m_spectra, dimension);

	return true;
}

//-----------------------------------------------------------------------------

int ExecSpecNetworkEval::processRemainingShifts(	int shift1, int shift2, int shift3, 
													int PM1, int PM2, int PM3, int parentMassSpec1Spec2, int jumpsIndex) {
	int retParentMass = parentMassSpec1Spec2;
	vector<int> temp;
	
	if(shift2 == 0 && shift3 == 0)
	{
		spectrumTypes[2] = (parentMassSpec1Spec2 <= PM3) ? FULL : PRE;
		spectrumJumps[jumpsIndex] = PM3;
		retParentMass = (parentMassSpec1Spec2 <= PM3) ? PM3 : parentMassSpec1Spec2;
	}
	else if(shift2 == 0 && shift3 < 0)
	{
		spectrumTypes[2] = (parentMassSpec1Spec2 <= PM3) ? FULL : PRE;
		spectrumJumps[jumpsIndex] = PM3;
		retParentMass = (parentMassSpec1Spec2 <= PM3) ? PM3 : parentMassSpec1Spec2;
	}
	else if(shift2 == 0 && shift3 > 0)
	{	
		spectrumTypes[2] = (parentMassSpec1Spec2 <= PM3 + shift3) ? POST : SUB;
		retParentMass = (parentMassSpec1Spec2 <= PM3 + shift3) ? PM3 + shift3 : parentMassSpec1Spec2;
		spectrumJumps[jumpsIndex] = shift3;
		jumpsIndex++;
		spectrumJumps[jumpsIndex] = (spectrumTypes[2]  == SUB) ?  PM3 + shift3 : 0;
					
		if(parentMassSpec1Spec2 < PM3 + shift3)
		{				
			temp = spectrumJumps;
			int tempIndex = 0;
			if(spectrumTypes[0] == POST)
			{
				spectrumTypes[0] = SUB;
				spectrumJumps[0] = temp[0];
				spectrumJumps[1] = parentMassSpec1Spec2;
				jumpsIndex = 2;
				tempIndex = 1;
			}
			else if(spectrumTypes[0] == FULL)
			{
				spectrumTypes[0] = PRE;
				spectrumJumps[0] = parentMassSpec1Spec2;
				jumpsIndex = 1;
			}
			else if(spectrumTypes[0] == PRE)
			{
				spectrumJumps[0] = temp[0];
				tempIndex = 1;
				jumpsIndex = 1;
			}
			else if(spectrumTypes[0] == SUB)
			{
				spectrumJumps[0] = temp[0];
				spectrumJumps[1] = temp[1];
				jumpsIndex = 2;
				tempIndex = 2;
			}
						
			if(spectrumTypes[1] == POST)
			{
				DEBUG_MSG("Should not go here: spectrumTypes[1] == POST")
				spectrumTypes[1] = SUB;
				spectrumJumps[jumpsIndex] = temp[tempIndex];
				jumpsIndex++;
				spectrumJumps[jumpsIndex] = PM2;
				jumpsIndex++;
			}
			else if(spectrumTypes[1] == FULL)
			{
				spectrumTypes[1] = PRE;
				spectrumJumps[jumpsIndex] = parentMassSpec1Spec2;
				jumpsIndex++;
			}
			else if(spectrumTypes[1] == PRE)
			{
				spectrumTypes[1] = SUB;
				spectrumJumps[jumpsIndex] = temp[tempIndex];
				jumpsIndex++;
				spectrumJumps[jumpsIndex] = PM2;
				jumpsIndex++;
			}
			else if(spectrumTypes[1] == SUB)
			{
				spectrumJumps[jumpsIndex] = temp[tempIndex];
				jumpsIndex++;
				tempIndex++;
				spectrumJumps[jumpsIndex] = temp[tempIndex];
				jumpsIndex++;
			}
			spectrumJumps[jumpsIndex] = shift3;
			jumpsIndex++;
			spectrumJumps[jumpsIndex] = (spectrumTypes[2]  == SUB) ?  PM3 + shift3 : 0;
			
		}
	}
	else if(shift3 == 0 && shift2 < 0)
	{
		spectrumTypes[2] = (parentMassSpec1Spec2 <= PM3) ? FULL : PRE;
		spectrumJumps[jumpsIndex] = PM3;
		retParentMass = (parentMassSpec1Spec2 <= PM3) ? PM3 : parentMassSpec1Spec2;
	}
	else if(shift3 == 0 && shift2 > 0)
	{
		spectrumTypes[2] = (parentMassSpec1Spec2 <= PM3 + shift2) ? POST : SUB;
		retParentMass = (parentMassSpec1Spec2 <= PM3 + shift2) ? PM3 + shift2 : parentMassSpec1Spec2;
		spectrumJumps[jumpsIndex] = shift2;
		jumpsIndex++;
		spectrumJumps[jumpsIndex] = (spectrumTypes[2]  == SUB) ?  PM3 + shift2 : 0;
						
		if(parentMassSpec1Spec2 < PM3 + shift2)
		{				
			temp = spectrumJumps;
			int tempIndex = 0;
			if(spectrumTypes[0] == POST)
			{
				DEBUG_MSG("Should not go here: spectrumTypes[0] == POST");
			
			}
			else if(spectrumTypes[0] == FULL)
			{
				spectrumTypes[0] = PRE;
				spectrumJumps[0] = PM1;
				jumpsIndex = 1;
			}
			else if(spectrumTypes[0] == PRE)
			{
				spectrumJumps[0] = temp[0];
				tempIndex = 1;
				jumpsIndex = 1;
			}
			else if(spectrumTypes[0] == SUB)
			{
				spectrumJumps[0] = temp[0];
				spectrumJumps[1] = temp[1];
				jumpsIndex = 2;
				tempIndex = 2;
			}
							
			if(spectrumTypes[1] == POST)
			{
				spectrumTypes[1] = SUB;
				spectrumJumps[jumpsIndex] = temp[tempIndex];
				jumpsIndex++;
				spectrumJumps[jumpsIndex] = PM2;
				jumpsIndex++;
			}
			else if(spectrumTypes[1] == FULL)
			{
				DEBUG_MSG("Should not go here: spectrumTypes[1] == FULL");
			}
			else if(spectrumTypes[1] == PRE)
			{
				DEBUG_MSG("Should not go here: spectrumTypes[1] == PRE");
			}
			else if(spectrumTypes[1] == SUB)
			{
				spectrumJumps[jumpsIndex] = temp[tempIndex];
				jumpsIndex++;
				tempIndex++;
				spectrumJumps[jumpsIndex] = temp[tempIndex];
				jumpsIndex++;
			}
			spectrumJumps[jumpsIndex] = shift3;
			jumpsIndex++;
			spectrumJumps[jumpsIndex] = (spectrumTypes[2]  == SUB) ?  PM3 + shift2 : 0;		
		}	
	}
	else if(shift2 < 0 && shift3 < 0)
	{
		if(parentMassSpec1Spec2 + abs(shift2) <= PM3 && parentMassSpec1Spec2 + abs(shift3) <= PM3)
		{
			spectrumTypes[2] = FULL;
			retParentMass = PM3;
		}
		else
		{
			spectrumTypes[2] = PRE;
			retParentMass = (parentMassSpec1Spec2 + abs(shift2) <= parentMassSpec1Spec2 + abs(shift3) ) ? parentMassSpec1Spec2 + abs(shift2) : parentMassSpec1Spec2 + abs(shift3);
		}
		vector<int> temp = spectrumJumps;
		pair<int, int> p = processhelper_negShift(0, 0, retParentMass, parentMassSpec1Spec2, abs(shift2),temp, 0);
		pair<int, int> p2 = processhelper_negShift(1, p.first, retParentMass, parentMassSpec1Spec2, abs(shift3), temp, p.second);																
		spectrumJumps[p2.first] = (spectrumTypes[2] == PRE) ? PM3 : 0;	
	}
	else if(shift2 < 0 && shift3 > 0)
	{
		if(PM3 + shift3 > parentMassSpec1Spec2)
		{
			spectrumTypes[2] = POST;
			retParentMass = PM3 + shift3;
			
			temp = spectrumJumps;
			int tempIndex = 0;
			int jumpsIndex = 0;
			if(spectrumTypes[0] == POST)
			{
				spectrumTypes[0] = SUB;
				spectrumJumps[jumpsIndex] = temp[tempIndex];
				jumpsIndex++;
				spectrumJumps[jumpsIndex] = temp[tempIndex] + PM1;
				tempIndex++;
				jumpsIndex++;				
			}
			else if(spectrumTypes[0] == FULL)
			{
				spectrumTypes[0] = PRE;
				spectrumJumps[jumpsIndex] = PM1;
				jumpsIndex++;
			}
			else if(spectrumTypes[0] == SUB)
			{
				spectrumTypes[0] = SUB;
				spectrumJumps[jumpsIndex] = temp[tempIndex];
				tempIndex++;
				jumpsIndex++;
				spectrumJumps[jumpsIndex] = temp[tempIndex];
				tempIndex++;
				jumpsIndex++;				
			}
			else
			{
				DEBUG_MSG("Should not go here");
			}
			
			if(spectrumTypes[1] == PRE)
			{
				spectrumTypes[1] = PRE;
				spectrumJumps[jumpsIndex] = PM2;
				jumpsIndex++;
				tempIndex++;
			}
			else if(spectrumTypes[1] == FULL)
			{
				spectrumTypes[1] = PRE;
				spectrumJumps[jumpsIndex] = PM2;
				jumpsIndex++;
				tempIndex++;
			}
			else
			{
				DEBUG_MSG("Should not go here");
			}
	
			spectrumJumps[jumpsIndex] = shift3;	
		}
		else
		{
			spectrumTypes[2] = SUB;
			spectrumJumps[jumpsIndex] = shift3;
			jumpsIndex++;
			spectrumJumps[jumpsIndex] = abs(shift2) + PM3;
		}	
	}
	else if(shift2 > 0 && shift3 < 0)
	{
		if(PM3 + shift2 > parentMassSpec1Spec2)
		{
			spectrumTypes[2] = POST;
			retParentMass = PM3 + shift2;
			
			temp = spectrumJumps;
			int tempIndex = 0;
			if(spectrumTypes[0] == PRE)
			{
				spectrumTypes[0] = PRE;
				spectrumJumps[0] = temp[tempIndex];
				jumpsIndex = 1;
				tempIndex++;
			}
			else if(spectrumTypes[0] == FULL)
			{
				spectrumTypes[0] = PRE;
				spectrumJumps[0] = PM1;
				jumpsIndex = 1;
				tempIndex = 0;
			}
			else
			{
				DEBUG_MSG("Should not go here");
			}
				
			if(spectrumTypes[1] == POST)
			{
				spectrumTypes[1] = SUB;
				spectrumJumps[jumpsIndex] = temp[tempIndex];
				jumpsIndex++;
				tempIndex++;
				spectrumJumps[jumpsIndex] = temp[tempIndex] + PM2;
				tempIndex++;
				jumpsIndex++;				
			}
			else if(spectrumTypes[1] == FULL)
			{
				spectrumTypes[1] = PRE;
				spectrumJumps[jumpsIndex] = PM1;
				jumpsIndex++;
			}
			else if(spectrumTypes[1] == SUB)
			{
				spectrumTypes[1] = SUB;
				spectrumJumps[jumpsIndex] = temp[tempIndex];
				tempIndex++;
				jumpsIndex++;
				spectrumJumps[jumpsIndex] = temp[tempIndex];
				tempIndex++;
				jumpsIndex++;				
			}
			else
			{
				DEBUG_MSG("Should not go here");
			}
				spectrumJumps[jumpsIndex] = shift2;	
		}
		else
		{
			spectrumTypes[2] = SUB;
			spectrumJumps[jumpsIndex] = shift2;
			jumpsIndex++;
			spectrumJumps[jumpsIndex] = abs(shift3) + PM3;
		}			
	}
	else if(shift2 > 0 && shift3 > 0)
	{
			if(shift2 > shift3)
			{
				if(PM3 + shift2 > parentMassSpec1Spec2)
				{
					spectrumTypes[2] = POST;
					retParentMass = PM3 + shift2;	
					temp = spectrumJumps;
					int tempIndex = 0;
					if(spectrumTypes[0] == PRE)
					{
						spectrumTypes[0] = PRE;
						spectrumJumps[0] = temp[tempIndex];
						jumpsIndex = 1;
						tempIndex++;						
					}
					else if(spectrumTypes[0] == FULL)
					{
						spectrumTypes[0] = PRE;
						spectrumJumps[0] = PM1;
						jumpsIndex = 1;
						tempIndex = 0;
					}
					else
					{
						DEBUG_MSG("Should not go here");
					}
									
					if(spectrumTypes[1] == POST)
					{
						spectrumTypes[1] = SUB;
						spectrumJumps[jumpsIndex] = temp[tempIndex];
						jumpsIndex++;
						tempIndex++;
						spectrumJumps[jumpsIndex] = temp[tempIndex] + PM2;
						tempIndex++;
						jumpsIndex++;				
					}
					else if(spectrumTypes[1] == SUB)
					{
						spectrumTypes[1] = SUB;
						spectrumJumps[jumpsIndex] = temp[tempIndex];
						tempIndex++;
						jumpsIndex++;
						spectrumJumps[jumpsIndex] = temp[tempIndex];
						tempIndex++;
						jumpsIndex++;				
					}
								
					spectrumJumps[jumpsIndex] = shift2 + shift3;
				}
				else // PM3 + shift2 <= parentMassSpec1Spec2
				{
					spectrumTypes[2] = SUB;
					spectrumJumps[jumpsIndex] = shift2;
					jumpsIndex++;
					spectrumJumps[jumpsIndex] = parentMassSpec1Spec2 - (PM3 + shift2);
				}		
			}
			else //  shift2 <= shift3)
			{
				if(PM3 + shift3 > parentMassSpec1Spec2)
				{
					spectrumTypes[2] = POST;
					retParentMass = PM3 + shift3;	
					temp = spectrumJumps;
					int tempIndex = 0;
					if(spectrumTypes[0] == POST)
					{
						spectrumTypes[0] = SUB;
						spectrumJumps[0] = temp[tempIndex];
						jumpsIndex = 1;
						spectrumJumps[jumpsIndex] = temp[tempIndex] + PM1;
						tempIndex++;
						jumpsIndex++;
					}
					else if(spectrumTypes[0] == FULL)
					{
						spectrumTypes[0] = PRE;
						spectrumJumps[0] = PM1;
						jumpsIndex = 1;
						tempIndex = 0;
					}
					else if(spectrumTypes[0] == SUB)
					{
						spectrumTypes[0] = SUB;
						jumpsIndex = 0;
						spectrumJumps[jumpsIndex] = temp[tempIndex];
						tempIndex++;
						jumpsIndex++;
						spectrumJumps[jumpsIndex] = temp[tempIndex];
						tempIndex++;
						jumpsIndex++;				
					}
					else
					{
						DEBUG_MSG("Should not go here");
					}
													
					if(spectrumTypes[1] == PRE)
					{
						spectrumTypes[1] = PRE;
						spectrumJumps[jumpsIndex] = temp[tempIndex];
						jumpsIndex++;
						tempIndex++;
					}
					else if(spectrumTypes[1] == FULL)
					{
						spectrumTypes[1] = PRE;
						spectrumJumps[jumpsIndex] = PM2;
						jumpsIndex++;
					}
					else if(spectrumTypes[1] == SUB)
					{
						spectrumTypes[1] = SUB;
						spectrumJumps[jumpsIndex] = temp[tempIndex];
						tempIndex++;
						jumpsIndex++;
						spectrumJumps[jumpsIndex] = temp[tempIndex];
						tempIndex++;
						jumpsIndex++;		
					}
												
					spectrumJumps[jumpsIndex] = shift2 + shift3;
				}
				else // PM3 + shift3 <= parentMassSpec1Spec2
				{	
					spectrumTypes[2] = SUB;
					spectrumJumps[jumpsIndex] = shift3;
					jumpsIndex++;
					spectrumJumps[jumpsIndex] = parentMassSpec1Spec2 - (PM3 + shift3);					
				}
			}
	}
	else
	{	
		DEBUG_MSG("Should not go here");
	}
	return retParentMass;
	
}
//---------------------------------------------------------------------------
pair<int, int> ExecSpecNetworkEval::processhelper_negShift(int i, int jumpsIndex, int newParentMass, int oldParentMass, int shift,
												vector<int> temp, int tempIndex) 
{

	if(spectrumTypes[i] == PRE)
	{
		spectrumTypes[i] = SUB;
		spectrumJumps[jumpsIndex] = shift; 
		jumpsIndex++;
		spectrumJumps[jumpsIndex] = shift + temp[tempIndex];
		jumpsIndex++;
		tempIndex++;
	}
	else if(spectrumTypes[i] == POST && newParentMass - shift == oldParentMass)
	{
		spectrumJumps[i] = POST;
		spectrumJumps[jumpsIndex] = shift + temp[tempIndex];
		jumpsIndex++;
		tempIndex++;
	}
	else if(spectrumTypes[i] == POST)
	{
		spectrumJumps[i] = SUB;
		spectrumJumps[jumpsIndex] = shift; 
		jumpsIndex++;
		spectrumJumps[jumpsIndex] = shift + temp[tempIndex];
		jumpsIndex++;
		tempIndex++;
	}
	else if(spectrumTypes[i] == FULL && newParentMass - shift == oldParentMass)
	{
		spectrumJumps[i] = POST;
		spectrumJumps[jumpsIndex] = shift + temp[tempIndex];
		jumpsIndex++;
		tempIndex++;
	}
	else if(spectrumTypes[i] == FULL)
	{
		spectrumJumps[i] = SUB;
		spectrumJumps[jumpsIndex] = shift; 
		jumpsIndex++;
		spectrumJumps[jumpsIndex] = shift + oldParentMass;
		jumpsIndex++;
		tempIndex++;
	}
	else if(spectrumTypes[i] == SUB)
	{
		spectrumJumps[i] = SUB;
		spectrumJumps[jumpsIndex] = temp[tempIndex];
		tempIndex++;
		jumpsIndex++;
		spectrumJumps[jumpsIndex] = temp[tempIndex];
		tempIndex++;
		jumpsIndex++;
				
	}
	else 
	{				
		DEBUG_MSG("Should not go here");		
	}
	pair<int, int> p(jumpsIndex, tempIndex);
	return p;
}


}
