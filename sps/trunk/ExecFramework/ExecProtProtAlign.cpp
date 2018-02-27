// Header Includes
#include "ExecProtProtAlign.h"

// Module Includes
#include "Logger.h"
//#include "FileUtils.h"

// SpecNets Includes

// Other Includes
#include "Specific.h"

// standard Includes
#include <string>
#include <stdio.h>

using namespace std;
using namespace specnets;

namespace specnets
{
	ExecProtProtAlign::ExecProtProtAlign(void) :
		m_db(0x0), ownInput(true),
		m_alnFileNames(0x0), ownOutput(true)
	{
		m_name = "ExecProtProtAlign";
		m_type = "ExecProtProtAlign";
	}

	// -------------------------------------------------------------------------

	ExecProtProtAlign::ExecProtProtAlign(const ParameterList & inputParams) :
		ExecBase(inputParams), m_db(0x0), ownInput(true),
		m_alnFileNames(0x0), ownOutput(true)
	{
		m_name = "ExecProtProtAlign";
		m_type = "ExecProtProtAlign";
	}

	// -------------------------------------------------------------------------

	ExecProtProtAlign::ExecProtProtAlign(const ParameterList & inputParams,
                                 			 DB_fasta * db,
                                  		 set<unsigned int> * dbIndexes,
                                  		 vector<string> * alnFileNames) :
	  ExecBase(inputParams),
		  m_db(db),
		  m_dbIndexes(dbIndexes),
		  ownInput(false),
		  m_alnFileNames(alnFileNames),
		  ownOutput(false)
	{
		m_name = "ExecProtProtAlign";
		m_type = "ExecProtProtAlign";
	}

	// -------------------------------------------------------------------------

	ExecProtProtAlign::~ExecProtProtAlign(void)
	{
		if (ownInput)
		{
			if(m_db) delete m_db;
		}
		if (ownOutput)
		{
			if(m_alnFileNames) delete m_alnFileNames;
		}
	}

	// -------------------------------------------------------------------------

	ExecBase * ExecProtProtAlign::clone(const ParameterList & inputParams) const
	{
		return new ExecProtProtAlign(inputParams);
	}

	// -------------------------------------------------------------------------

	bool ExecProtProtAlign::invoke(void)
	{
	  DEBUG_TRACE;
		if(!m_db or !m_alnFileNames) return false;
		if(m_params.getValue("CLUSTALW_EXE_DIR").size()==0) return false;
    unsigned int refProtIdx = (unsigned int) m_params.getValueInt("REFERENCE_PROTEIN_IDX");
		if(refProtIdx < 0 or refProtIdx >= m_db->size()) {
			ERROR_MSG("ERROR: Reference Protein Index out of bounds [" << refProtIdx << "]  DB Size = [" << m_db->size() << "]");
			return false;
    }

	  DEBUG_TRACE;
		if(!m_db or m_db->size()==0) {
			ERROR_MSG("ERROR: empty database");
			return false;
		}

	  DEBUG_TRACE;
    float minScore = (float) m_params.getValueDouble("CLUSTALW_MINSCORE");
    if (minScore == -1.0) {
      DEBUG_MSG("CLUSTALW_MINSCORE = -1.0. Exiting.");
      return true;
    }

    string cmdBase = m_params.getValue("CLUSTALW_EXE_DIR");
    if(cmdBase[cmdBase.size()-1] != '/') cmdBase = cmdBase + "/";
    cmdBase = cmdBase + "clustalw ";

    string filenamePrefix;
    if(m_params.exists("CLUSTALW_FILENAME_PREFIX"))
    	filenamePrefix = m_params.getValue("CLUSTALW_FILENAME_PREFIX");
    else filenamePrefix = "clustalw_";

    m_alnFileNames->resize(m_db->size());
	  DEBUG_VAR(m_alnFileNames->size());
    m_numAligns = 0;  // Number of clustalw alignments with a score of at least CLUSTALW_MINSCORE

    if (m_dbIndexes->size() == 0) {
      for (unsigned int i = 0; i < m_db->size(); i++) {
       	if (i != refProtIdx) {
       	  if (!score(filenamePrefix, cmdBase, i, refProtIdx, minScore)) {
       	    return false;
       	  }
      	}
      }
    } else {
      set<unsigned int>::iterator itr = m_dbIndexes->begin();
      set<unsigned int>::iterator itr_end = m_dbIndexes->end();
      for (; itr != itr_end; itr++) {
       	if (*itr != refProtIdx) {
       	  if (!score(filenamePrefix, cmdBase, *itr, refProtIdx, minScore)) {
       	    return false;
       	  }
      	}
      }
    }

    m_alnFileNames->resize(m_numAligns);
	  DEBUG_VAR(m_alnFileNames->size());

    return true;
	}

	// -------------------------------------------------------------------------

	bool ExecProtProtAlign::loadInputData(void)
	{
    if(ownInput) {
    	if(!m_db) m_db = new DB_fasta;
    }

  	if(ownOutput) {
    	if(!m_alnFileNames) m_alnFileNames = new vector<string>;
    }
  	m_alnFileNames->resize(0);

    if (!m_params.exists("INPUT_FASTA")) {
			ERROR_MSG("Parameters are incomplete. INPUT_FASTA is missing.");
      return false;
    } else if (m_db->Load(m_params.getValue("INPUT_FASTA").c_str())<=0) {
			ERROR_MSG("Error reading database sequences from "<<m_params.getValue("INPUT_FASTA"));
			return false;
		}

		return true;
	}

	// -------------------------------------------------------------------------

	bool ExecProtProtAlign::saveOutputData(void)
	{
    if (m_alnFileNames and m_params.exists("CLUSTALW_INDEX")) {
   		FILE *fileID = fopen(m_params.getValue("CLUSTALW_INDEX").c_str(),"wb");
   		if(fileID <= 0) {
   			return false;
   		}
   		for(unsigned int i=0; i<m_alnFileNames->size(); i++)
   			fprintf(fileID,"%s\n", (*m_alnFileNames)[i].c_str());
   		fclose(fileID);
    }

		return true;
	}

	// -------------------------------------------------------------------------

	bool ExecProtProtAlign::saveInputData(std::vector<std::string> & filenames)
	{
		return false;
	}

	// -------------------------------------------------------------------------

	bool ExecProtProtAlign::loadOutputData(void)
	{
		return false;
	}

	// -------------------------------------------------------------------------

  vector<ExecBase*> const & ExecProtProtAlign::split(int numSplit)
	{
		m_subModules.resize(0);
		return m_subModules;
	}

	// -------------------------------------------------------------------------

	bool ExecProtProtAlign::merge(void)
	{
		return false;
	}

	// -------------------------------------------------------------------------

	bool ExecProtProtAlign::validateParams(std::string & error)
	{
    m_isValid = false;

    VALIDATE_PARAM_EXIST("CLUSTALW_EXE_DIR");
    VALIDATE_PARAM_EXIST("REFERENCE_PROTEIN_IDX");
    VALIDATE_PARAM_EXIST("CLUSTALW_MINSCORE");

    m_isValid = true;
    return true;
	}

	// -------------------------------------------------------------------------

  int ExecProtProtAlign::getAlignmentScore(string outputFileName) {
  	BufferedLineReader blr;

  	if(blr.Load(outputFileName.c_str())<0) return -1;
  	for(unsigned int lineIdx=0; lineIdx<blr.size(); lineIdx++)
  		if(strncmp(blr.getline(lineIdx),"Alignment Score ",16)==0)
  			return (int)atoi(&blr.getline(lineIdx)[16]);
  	return -1;
  }

	// -------------------------------------------------------------------------
	bool ExecProtProtAlign::score(string & filenamePrefix,
	                              string & cmdBase,
	                              unsigned int protIdx,
	                              unsigned int refProtIdx,
	                              float minScore)
	{
	  DEBUG_VAR(protIdx);
	  DEBUG_VAR(m_db->getSequence(protIdx));
	  DEBUG_VAR(refProtIdx);
	  DEBUG_VAR(m_db->getSequence(refProtIdx));

    ostringstream fileName;
	  fileName.str("");
	  fileName << filenamePrefix << refProtIdx << "_" << protIdx;

    string cmd = cmdBase + fileName.str() + ".fasta > " + fileName.str() + "_output.txt";
	  fileName << ".fasta";

	  FILE *fileID = fopen(fileName.str().c_str(),"wb");
	  if(fileID <= 0) {
		  ERROR_MSG("ERROR: Can not open " << fileName);
		  return false;
	  }
	  fprintf(fileID,
	          "> (Reference protein index %d) %s %s\n%s\n\n",
	          refProtIdx,
	          m_db->getID(refProtIdx),
	          m_db->getDesc(refProtIdx),
	          m_db->getSequence(refProtIdx));
	  fprintf(fileID,
	          "> (Protein index %d) %s %s\n%s\n\n",
	          protIdx,
	          m_db->getID(protIdx),
	          m_db->getDesc(protIdx),
	          m_db->getSequence(protIdx));
	  fclose(fileID);

    int status = spsSystem(cmd.c_str());
    if (status != 0) {
      ERROR_MSG("Error executing: " << cmd);
      return false;
    }

	  fileName.str("");
	  fileName << filenamePrefix << refProtIdx << "_" << protIdx;
    int curScore = getAlignmentScore(fileName.str()+"_output.txt");
	  if(curScore >= minScore) {
   		fileName.str("");
   		fileName << refProtIdx << ";" << protIdx << ";" << filenamePrefix << refProtIdx
   		         << "_" << protIdx << ".aln";
		  (*m_alnFileNames)[m_numAligns++] = fileName.str();
	  }
	  return true;
  }

}
