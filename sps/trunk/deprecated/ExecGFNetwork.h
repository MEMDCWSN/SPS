/*
 * ExecGFNetwork.h
 *
 *  Created on: Feb 18, 2011
 *      Author: cabouche
 */

#ifndef EXECGFNETWORK_H_
#define EXECGFNETWORK_H_

// Module Includes
#include "ExecBase.h"
#include "SpectrumPairSet.h"

// External Includes
#include "spectrum.h"
#include "aminoacid.h"
#include "PeptideSpectrumMatchSet.h"


namespace specnets {

enum spectra_type_t {PRE, POST, SUB, FULL};

class ExecGFNetwork: public ExecBase {
public:

		ExecGFNetwork(void);
		ExecGFNetwork(const ParameterList & inputParams);
		ExecGFNetwork(SpecSet *, int);
		virtual ~ExecGFNetwork(void);
		virtual bool invoke(void);
		virtual bool loadInputData(void);
		virtual ExecBase * clone(const ParameterList & input_params) const;
		virtual bool saveOutputData(void);
		virtual bool saveInputData(std::vector<std::string> & filenames);
		virtual bool loadOutputData(void);
		virtual bool merge(void);
		virtual bool validateParams(std::string & error);
		virtual std::vector<ExecBase *> const & split(int numSplit);

		vector<double> scorePSMPairs(vector<int>, vector<spectra_type_t>, vector<string>, vector<int>, int);
		vector<double> scorePSMTriple(vector<int>, vector<spectra_type_t>, vector<string>, vector<int>, int);


		void outputForPairs(vector<int>, vector<spectra_type_t>, string, vector<int>, int);
		void outputPSMTriple(vector<int>, vector<spectra_type_t>, string, vector<int>, int);

		float computeSpectraProbability(int, string, int, int);

		vector<string> partialAnnotations;
		vector<int> partialPeptideScores;
		vector<double> partialSpecProbabilities;



		/*************************************************************************************/
		inline int scorePeptide(string peptide, int index) {

			AAJumps myjumps(1);
			vector<float> mymasses;
			float mass_float = myjumps.getPeptideMass(peptide);
			int mass_int = (int)round(mass_float*0.9995);

			vector<char> sequence = vector<char> (peptide.begin(), peptide.end());
			vector<float> masses;
			getMasses(sequence, masses);

			double partialPeptideMass = 0;
			float peptideScore = 0;

			for (int i = 0; i < masses.size(); i++) {
				partialPeptideMass += masses[i] * 0.9995;
				int m = (int) (partialPeptideMass + 0.5);

				if (m - 1 < (*scores)[index].size() && m < mass_int) {
					peptideScore += (*scores)[index][m - 1][1];

			//		DEBUG_MSG("mass: " << m << ", Score: " <<  peptideScore);

				}
			}
			return (int)round(peptideScore);
		}


private:
		const static float constantParentMass=0.9995;

		double computeMultiSpectraProbability(int, int);

		void computeMultiSpectraProbability_K3(int, int);

		int setDeltaPosition(vector<int>);
		void computeMassJumps();

		/*************************************************************************************/
		inline string computePartialPeptide(string peptide, spectra_type_t type, int massJump1, int massJump2) {

			AAJumps myjumps(1);
			vector<float> mymasses;
			float mass_float = myjumps.getPeptideMass(peptide);
			int PM = (int) round(mass_float*0.9995);

			if(type == FULL)
			{
			///	DEBUG_MSG("FULL");
				return peptide;
			}
			if(type == PRE)
			{

			//	DEBUG_MSG("PRE");
				  for(int p = 1; p < peptide.length(); p++ )
				  {
					  AAJumps myjumps(1);
					  vector<float> mymasses;
					  string prefix = peptide.substr(0, p);
					  float mass = myjumps.getPeptideMass(prefix);
					  int mass_int = (int) round(mass * 0.9995);
					  if(mass_int == (massJump1))
					  {
						  return prefix;
					  }
				  }

			}
			if(type == POST)
			{

			//	DEBUG_MSG("POST " << massJump1);
				 for(int p = 1; p < peptide.length(); p++ )
				 {
					 AAJumps myjumps(1);
					 vector<float> mymasses;
					 string suffix = peptide.substr(p, peptide.length() - p);
					 float mass = myjumps.getPeptideMass(suffix);
					 int mass_int = (int) round(mass * 0.9995);

		//			 DEBUG_MSG(suffix << " " << mass_int);
					 if(PM - mass_int == massJump1)
					 {
						return suffix;
					 }
				}
			}
			if(type == SUB)
			{

			//	DEBUG_MSG("SUB");
				 int p;
				 for(p = 1; p < peptide.length(); p++ )
				 {
					 AAJumps myjumps(1);
					 vector<float> mymasses;
					 string prefix = peptide.substr(0, p);
					 float mass = myjumps.getPeptideMass(prefix);
					 int mass_int = (int) round(mass * 0.9995);
					 if(mass_int == massJump1)
					 {
						continue;
					 }
				}

				for(int len = 1; len <= peptide.length() - p; len++ )
				{
					AAJumps myjumps(1);
					vector<float> mymasses;
					string suffix = peptide.substr(p, len);
					float mass = myjumps.getPeptideMass(suffix);
					int mass_int = (int) round(mass * 0.9995);
					if(mass_int == (massJump2 - massJump1))
					{
						return suffix;
					}
				}
			}
			else
			{
				DEBUG_MSG("Problem with computePartialPeptide");
			}
			return "";
		}



		/*************************************************************************************/
		inline int scorePeptide(vector<float> masses, int index) {

			double partialPeptideMass = 0;
			int peptideScore = 0;

			int m;
			for (int i = 0; i < masses.size(); i++) {
				partialPeptideMass += masses[i] * 0.9995;
				m = (int) round(partialPeptideMass);

				if (m - 1 < (*scores)[index].size()) {
					peptideScore += (int) (*scores)[index][m - 1][1];


				}
			}
		//	DEBUG_MSG("Index: " << index << ", Score: " <<  peptideScore);
			return peptideScore;
		}


		/*************************************************************************************/
		inline bool checkBoundaries(int mass, int aa) {


			int jumpIndex = 0;
			for(int n = 0; n < numberSpectra; n++)
			{
				if (types[n] == PRE || types[n] == POST)
				{
					if(mass - aa < jumps[jumpIndex] && mass > jumps[jumpIndex])
						return false;
					else
						jumpIndex++;

				}
				else if (types[n] == SUB)
				{
					if(mass - aa < jumps[jumpIndex] && mass > jumps[jumpIndex])
						return false;
					jumpIndex++;
					if(mass - aa < jumps[jumpIndex] && mass > jumps[jumpIndex])
						return false;
					jumpIndex++;
				}
			}
			return true;
		}

		void constructPSMMatrix(int, int, int);
		void constructPSMMatrixForTriple(int, int, int);

		double scoreSpectraPair(int, int, int, int);
		double scoreSpectraTriple(int, int, int, int, int);



		/*************************************************************************************/

		PeptideSpectrumMatchSet peptide_results;

		int currentParentMass;
		int numberSpectra;
		bool ownInput;

		spectra_type_t* types;
		int *indices;
	    int *jumps;
		int *peptideScores;
		float *singleSpecScore;

		double * dpMatrix;
		double *** x;
		double ** dpTripleMatrix;
		string * peptideStrings;

	    SpecSet * scores;

	};

} // namespace specnets

#endif /* EXECGFNETWORK_H_ */
