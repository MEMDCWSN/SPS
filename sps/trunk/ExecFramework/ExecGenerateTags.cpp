//
// Created by Benjamin Pullman on 7/11/16.
//

/*
 * ExecKLFilter.cpp
 *
 *  Created on: Oct 20, 2015
 *      Author: isna
 */

#include "ExecGenerateTags.h"
#include "ClusterSet.h"


namespace specnets {

    // -------------------------------------------------------------------------
    ExecGenerateTags::ExecGenerateTags(void) :
            ExecBase(), m_inputSpectra(0x0), ownInput(true), ownOutput(true)
    {
        m_inputSpectra = new vector<SpecSet> ();
        m_name = "ExecGenerateTags";
        m_type = "ExecGenerateTags";
    }

    // -------------------------------------------------------------------------
    ExecGenerateTags::ExecGenerateTags(const ParameterList & inputParams) :
            ExecBase(inputParams), m_inputSpectra(0x0), ownInput(true), ownOutput(true)
    {
        m_inputSpectra = new vector<SpecSet> ();
        m_name = "ExecGenerateTags";
        m_type = "ExecGenerateTags";
    }

    // -------------------------------------------------------------------------
//    ExecGenerateTags::ExecGenerateTags(const ParameterList & inputParams, vector<SpecSet> * inputSpectra) :
//            ExecBase(inputParams), m_inputSpectra(inputSpectra), m_filteredSpectra(0x0), ownInput(false), ownOutput(true)
//    {
//        m_name = "ExecGenerateTags";
//        m_type = "ExecGenerateTags";
//    }

    // -------------------------------------------------------------------------
    ExecGenerateTags::~ExecGenerateTags(void) {
        if( ownInput ){
            delete m_inputSpectra;
        }
    }

    // -------------------------------------------------------------------------
    ExecBase * ExecGenerateTags::clone(const ParameterList & inputParams) const {
        return new ExecGenerateTags(inputParams);
    }

    // -------------------------------------------------------------------------
    bool ExecGenerateTags::invoke(void) {

        return true;
    }

    // -------------------------------------------------------------------------
    bool ExecGenerateTags::loadInputData(void) {

      if( m_inputSpectra == 0x0 ) {
          ownInput = true;
          m_inputSpectra = new vector<SpecSet> ();
      }

      if (m_params.exists("INPUT_SPECS")) {
          short int maxMs = 1;
          vector<string> filenames;
          vector<string>::iterator filenames_iter;
          stringSplit(m_params.getValue("INPUT_SPECS"), filenames, ";");
//            SpecSet *all_input = new SpecSet;
//            SpecSet *all_input = new SpecSet;)
          for(filenames_iter = filenames.begin();filenames_iter != filenames.end(); filenames_iter++){
              std::vector<short> msLevel;
              SpecSet *temp = new SpecSet;
              if (filenames_iter->find(".mzXML") != -1) {
                  DEBUG_MSG("Trying to open mzXMLLLLLLLLLLLLLLLLLL");
                  LoadMzxml((*filenames_iter).c_str(),*filenames_iter,*temp,&msLevel,maxMs);
              } else {
                  temp->Load((*filenames_iter).c_str());
              }
              DEBUG_MSG("Loaded Spectra: " << temp->size());
              m_inputSpectra->push_back(*temp);
              DEBUG_MSG("Loaded Spectra: " << temp->size());
          }
      }

      return true;
    }

    // -------------------------------------------------------------------------
    bool ExecGenerateTags::saveOutputData(void) {

      if (m_params.exists("OUTPUT_FOLDER")) {
          vector<SpecSet>::iterator spectrum_iter;
          for(spectrum_iter = m_inputSpectra->begin();spectrum_iter != m_inputSpectra->end(); spectrum_iter++) {

            if( spectrum_iter->size() == 0 ) continue;

              string fileName = (spectrum_iter->begin())->fileName;
              string outputName = m_params.getValue("OUTPUT_FOLDER");
              string outputTag = m_params.getValue("OUTPUT_TAG",".tag");
              DEBUG_MSG("Saving mgf: " << fileName);
              if (fileName.find(".mzXML") != -1) {
                  fileName.replace(fileName.find(".mzXML"),6,outputTag);
              }
              if (fileName.find(".mgf") != -1) {
                  fileName.replace(fileName.find(".mgf"),4,outputTag);
              }
              outputName += "/";
              outputName += fileName;
              spectrum_iter->saveSNRRankedTags2(outputName.c_str(), 4, 0, 20, 0.02, true);

              if (m_params.exists("OUTPUT_FOLDER_DECONV")) {
                  string fileNameDeconv = (spectrum_iter->begin())->fileName;
                  string outputNameDeconv = m_params.getValue("OUTPUT_FOLDER_DECONV");
                  DEBUG_MSG("Saving mgf: " << fileNameDeconv);
                  if (fileNameDeconv.find(".mzXML") != -1) {
                      fileNameDeconv.replace(fileNameDeconv.find(".mzXML"),6,".mgf");
                  }
                  if (fileNameDeconv.find(".mzML") != -1) {
                      fileNameDeconv.replace(fileNameDeconv.find(".mzML"),5,".mgf");
                  }
                  outputNameDeconv += "/";
                  outputNameDeconv += fileNameDeconv;
                  spectrum_iter->SaveSpecSet_mgf(outputNameDeconv.c_str(), true);
              }
          }
      }
      return true;
    }

    // -------------------------------------------------------------------------
    bool ExecGenerateTags::saveInputData(std::vector<std::string> & filenames) {
        return true;
    }

    // -------------------------------------------------------------------------
    bool ExecGenerateTags::loadOutputData(void) {
        return true;
    }

    // -------------------------------------------------------------------------
    vector<ExecBase*> const & ExecGenerateTags::split(int numSplit) {
        return m_subModules;
    }

    // -------------------------------------------------------------------------
    bool ExecGenerateTags::merge(void) {
        return true;
    }

    // -------------------------------------------------------------------------
    bool ExecGenerateTags::validateParams(std::string & error) {

        return true;

    }


} /* namespace specnets */
