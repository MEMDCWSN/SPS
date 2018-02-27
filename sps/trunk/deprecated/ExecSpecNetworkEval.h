/*
 * ExecSpecNetworkEval.h
 *
 *  Created on: Jul 2011
 *      Author: cboucher@ucsd.edu
 */

#ifndef EXECSPECNETWORKEVAL_H_
#define EXECSPECNETWORKEVAL_H_


// Module Includes
#include "ExecBase.h"
#include "ExecGFNetwork.h"
#include "SpectrumPairSet.h"
#include "PeptideSpectrumMatchNetwork.h"

// External Includes
#include "spectrum.h"
#include "db_fasta.h"
#include "aminoacid.h"
#include "SpectrumPairSet.h"

namespace specnets {

typedef std::tr1::shared_ptr<PeptideSpectrumMatchNetwork> psmNetworkPtr;


class ExecSpecNetworkEval: public ExecBase {
public:

		ExecSpecNetworkEval(void);
		ExecSpecNetworkEval(const ParameterList & inputParams);
		virtual ~ExecSpecNetworkEval(void);
		virtual bool invoke(void);
		virtual bool loadInputData(void);
		virtual ExecBase * clone(const ParameterList & input_params) const;
		virtual bool saveOutputData(void);
		virtual bool saveInputData(std::vector<std::string> & filenames);
		virtual bool loadOutputData(void);
		virtual bool merge(void);
		virtual bool validateParams(std::string & error);
		virtual std::vector<ExecBase *> const & split(int numSplit);
		
		void outputPSMSet(PeptideSpectrumMatchSet, PeptideSpectrumMatchSet); 
		void outputZeroFDRPSMSet(PeptideSpectrumMatchSet, PeptideSpectrumMatchSet);

		
private:
		
		std::tr1::shared_ptr<PeptideSpectrumMatchNetwork>  findInPSMSet(PeptideSpectrumMatchSet, int);

		vector<int> processShifts(float, int, int); 
		bool annotateByPairs(void); 
		bool annotateByTuples(void); 
		
		pair<int, int> processhelper_negShift(int i, int jumpsIndex, int newParentMass, int oldParentMass, int shift,
														vector<int> temp, int tempIndex);
		
		int processRemainingShifts(	int shift1, int shift2, int shift3, int PM1, int PM2, int PM3, int parentMassSpec1Spec2, int jumpsIndex); 

		void annotateSingleTuple(int, int, int, int, int); 
		void outputType(spectra_type_t);
		
		bool ownInput;

		ExecGFNetwork *gfNetwork;
		
		char *filename;
		DB_fasta * m_db; //! Input database of protein sequences
		SpecSet *m_spectra; //! Input spectra
		SpectrumPairSet * m_pairs;  //! Input set of pairs usable for propagations/networks
		
		vector<int> spectrumIndices;
		vector<spectra_type_t> spectrumTypes;
		vector<int> spectrumJumps;
		
		vector<int> annotatedSpectra;
		vector<int>::iterator intVecItr;


		int optimalTuplePeptideIndex;
		double optimalTupleScore;
		string optimalTupleAnnotation;
		vector<string> optPartialAnnotations;
		vector<int> optPartialPeptideScores;
		vector<double> optPartialSpecProbabilities;

		
		
		int currentParentMass;
			

		int dimension;
	};

} // namespace specnets

#endif /* EXECSPECNETWORKEVAL_H_ */
