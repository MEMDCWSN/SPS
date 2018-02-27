//
// Created by Benjamin Pullman on 7/11/16.
//

/*
 * ExecKLFilter.cpp
 *
 *  Created on: Oct 20, 2015
 *      Author: isna
 */

#include "ExecKLFilter.h"
#include "ClusterSet.h"
#include <iostream>
#include <string>
#include <libgen.h>

namespace specnets {

    // -------------------------------------------------------------------------
    ExecKLFilter::ExecKLFilter(void) :
            ExecBase(), m_inputSpectra(0x0), m_filteredSpectra(0x0), ownInput(true), ownOutput(true)
    {
        m_filteredSpectra = new vector<SpecSet> ();
        m_name = "ExecKLFilter";
        m_type = "ExecKLFilter";
    }

    // -------------------------------------------------------------------------
    ExecKLFilter::ExecKLFilter(const ParameterList & inputParams) :
            ExecBase(inputParams), m_inputSpectra(0x0), m_filteredSpectra(0x0), ownInput(true), ownOutput(true)
    {
        m_filteredSpectra = new vector<SpecSet> ();
        m_name = "ExecKLFilter";
        m_type = "ExecKLFilter";
    }

    // -------------------------------------------------------------------------
//    ExecKLFilter::ExecKLFilter(const ParameterList & inputParams, vector<SpecSet> * inputSpectra) :
//            ExecBase(inputParams), m_inputSpectra(inputSpectra), m_filteredSpectra(0x0), ownInput(false), ownOutput(true)
//    {
//        m_name = "ExecKLFilter";
//        m_type = "ExecKLFilter";
//    }

    // -------------------------------------------------------------------------
    ExecKLFilter::~ExecKLFilter(void) {
        if( ownInput ){
            delete m_inputSpectra;
        }

        if( ownOutput ){
            delete m_filteredSpectra;
        }

    }

    // -------------------------------------------------------------------------
    ExecBase * ExecKLFilter::clone(const ParameterList & inputParams) const {
        return new ExecKLFilter(inputParams);
    }

    // -------------------------------------------------------------------------
    bool ExecKLFilter::invoke(void) {



    	float isotope_space = 1.00235;

        string 	isotable_name           = m_params.getValue("ISOTOPE_TABLE");
        float 	peakTol                 = m_params.getValueFloat("PEAK_TOLERANCE",0.05);
        float 	kl_filter               = m_params.getValueFloat("STRICT_FILTER",0.0);
        int 	max_charge              = m_params.getValueInt("MAX_CHARGE",5);
        int 	min_correction          = m_params.getValueInt("MIN_CORRECTION",0);
        int 	max_correction          = m_params.getValueInt("MAX_CORRECTION",3);\
        int     low_mass                = m_params.getValueInt("LOW_MZ",0);
        int     remove_precursors       = m_params.getValueInt("REMOVE_PRECURSORS",0);
        int     snf                     = m_params.getValueFloat("SNF",0);
        int     only_ms2                = m_params.getValueInt("ONLY_MS2",0);

        specnets::IsoEnvelope* isoenv_strict = new specnets::IsoEnvelope();
        specnets::IsoEnvelope* isoenv = new specnets::IsoEnvelope();
        isoenv_strict->LoadModel(isotable_name.c_str());
        isoenv->LoadModel(isotable_name.c_str());
        vector<SpecSet>::iterator input_spectrum_iter;

        for(input_spectrum_iter = m_inputSpectra->begin(); input_spectrum_iter != m_inputSpectra->end(); input_spectrum_iter++) {

            SpecSet* filtered_specs = new SpecSet();

            int precursor_scan = 0;
            int corrected_charge = 0;
            int corrected_pm = 0;

            if (input_spectrum_iter->size()>0){

                string filename = ((*input_spectrum_iter)[0].fileName);

                DEBUG_VAR(filename);

                vector<specnets::Spectrum>::iterator ms2_iter;

                if (only_ms2 == 0) {

                    specnets::Spectrum ms1;

                    for (ms2_iter = input_spectrum_iter->begin(); ms2_iter != input_spectrum_iter->end(); ms2_iter++) {

                        if ((*ms2_iter).msLevel == 1) {
                            ms1 = *ms2_iter;
                            precursor_scan = ms2_iter->scan;
                        }

                        if ((*ms2_iter).msLevel == 2) {

                            unsigned short charge = (unsigned short) ms2_iter->parentCharge;
                            float parent_mass = ms2_iter->parentMass;

                            unsigned short best_charge = charge;

                            // set score to effective max

                            float best_score = 1000;

                            // guess charge if 0

                            int c_min = (charge != 0) ? charge : 1;
                            int c_max = (charge != 0) ? charge : max_charge;

                            float corrected_parent_mass = 0;
                            float corrected_precursor_mass = 0;

                            for (short iso = (short) min_correction; iso <= max_correction; iso++) {
                                float iso_shifted_mass = isotope_space * iso;
                                float temp_parent_mass = parent_mass - iso_shifted_mass;
                                for (unsigned short c = c_min; c <= c_max; c++) {
                                    float temp_precursor_mass =
                                            (float) (temp_parent_mass + (c - 1) * specnets::AAJumps::massHion) / c;
                                    std::vector<float> mass_env;
                                    isoenv->ExtractEnvelope(temp_precursor_mass, c, ms1, peakTol, mass_env, false);
                                    isoenv->normalize(mass_env, true, 0);
                                    float temp_score = isoenv_strict->ScoreEnvelope(temp_parent_mass, mass_env, false);
                                    if (temp_score < best_score) {
                                        corrected_parent_mass = temp_parent_mass;
                                        corrected_precursor_mass = temp_precursor_mass;
                                        best_charge = c;
                                        best_score = temp_score;
                                    }
                                }
                            }

                            ms2_iter->setParentMass(corrected_parent_mass);

                            std::vector<float> mass_env_strict;
                            isoenv_strict->ExtractEnvelope(corrected_precursor_mass, best_charge, ms1, peakTol,
                                                           mass_env_strict, true);

                            isoenv_strict->normalize(mass_env_strict, true, 0);

                            float strict_score = isoenv_strict->ScoreEnvelope(corrected_parent_mass, mass_env_strict,
                                                                              true);

                            if (charge != best_charge) {
                                corrected_charge++;
                            }

                            if (parent_mass != corrected_parent_mass) {
                                corrected_pm++;
                            }

                            ms2_iter->parentCharge = best_charge;
                            ms2_iter->precursor_kl = strict_score;
                            ms2_iter->precursor_scan = precursor_scan;

                            if (low_mass > 0) {
                                ms2_iter->filterLowMassPeaks(low_mass);
                            }

                            if (remove_precursors == 1) {
                                ms2_iter->removeChargeReducedPrecursors(10,false,0);
                            }

                            if (snf > 0.0) {
                                ms2_iter->filterLowSNR(snf);
                            }

                            if (kl_filter == 0) {
                                filtered_specs->push_back(*ms2_iter);
                            } else if (strict_score < kl_filter) {
                                filtered_specs->push_back(*ms2_iter);
                            }

                        }
                    }
                } else {
                    for (ms2_iter = input_spectrum_iter->begin(); ms2_iter != input_spectrum_iter->end(); ms2_iter++) {
                        if (low_mass > 0) {
                            ms2_iter->filterLowMassPeaks(low_mass);
                        }

                       if (remove_precursors == 1) {
                            ms2_iter->removeChargeReducedPrecursors(10,false,0);
                        }

                        if (snf > 0.0) {
                            ms2_iter->filterLowSNR(snf);
                        }

                        filtered_specs->push_back(*ms2_iter);

                    }
                }

            }
            DEBUG_VAR(filtered_specs->size());
            DEBUG_VAR(corrected_pm);
            DEBUG_VAR(corrected_charge);
            m_filteredSpectra->push_back(*filtered_specs);
        }

        return true;
    }

    // -------------------------------------------------------------------------
    bool ExecKLFilter::loadInputData(void) {

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
                    DEBUG_MSG("Trying to open mzXML");
                    LoadMzxml((*filenames_iter).c_str(),*filenames_iter,*temp,&msLevel,maxMs);
                } else if (filenames_iter->find(".mzML") != -1) {
                    temp->Load((*filenames_iter).c_str());
                } else {
//                    string mgf_name = basename((*filenames_iter).c_str());
                    temp->Load((*filenames_iter).c_str());
//                    temp->setFilename(mgf_name.c_str());
                }
                DEBUG_MSG("Loaded Spectra: " << temp->size());
                m_inputSpectra->push_back(*temp);
            }
        }

        return true;
    }

    // -------------------------------------------------------------------------
    bool ExecKLFilter::saveOutputData(void) {

        int min_peaks               = m_params.getValueFloat("MIN_PEAKS",0);

        if (m_params.exists("OUTPUT_FOLDER")) {
            vector<SpecSet>::iterator filtered_spectrum_iter;
            size_t current_file = 0;
            vector<string> filenames;
            stringSplit(m_params.getValue("INPUT_SPECS"), filenames, ";");
            for(filtered_spectrum_iter = m_filteredSpectra->begin();filtered_spectrum_iter != m_filteredSpectra->end(); filtered_spectrum_iter++) {

            	if( filtered_spectrum_iter->size() <= min_peaks ) continue;

                string fileName = (filtered_spectrum_iter->begin())->fileName;
                vector<string> filename_path;
                stringSplit(filenames[current_file], filename_path, "/");
                string outputName = m_params.getValue("OUTPUT_FOLDER");
                string outputTag = m_params.getValue("OUTPUT_TAG","_filtered") + ".mgf";
                if (fileName.find(".mzXML") != -1) {
                    fileName.replace(fileName.find(".mzXML"),6,outputTag);
                }
                if (fileName.find(".mzML") != -1) {
                    fileName.replace(fileName.find(".mzML"),5,outputTag);
                }
                if (filenames[current_file].find(".mgf") != -1) {
                    fileName = filename_path.back();
                    fileName.replace(fileName.find(".mgf"),4,outputTag);
                }
                outputName += "/";
                outputName += fileName;
                DEBUG_MSG("Saving mgf: " << outputName);
                filtered_spectrum_iter->SaveSpecSet_mgf(outputName.c_str(),true);
                if (m_params.exists("OUTPUT_FOLDER_PKL")) {
                    string fileNamePkl = (filtered_spectrum_iter->begin())->fileName;
                    string outputNamePkl = m_params.getValue("OUTPUT_FOLDER_PKL");
                    string outputTagPkl = m_params.getValue("OUTPUT_TAG","_filtered") + ".pklbin";
                    if (fileNamePkl.find(".mzXML") != -1) {
                        fileNamePkl.replace(fileNamePkl.find(".mzXML"),6,outputTagPkl);
                    }
                    if (filenames[current_file].find(".mgf") != -1) {
                        fileNamePkl = filename_path.back();
                        fileNamePkl.replace(fileNamePkl.find(".mgf"),4,outputTagPkl);
                    }
                    if (fileNamePkl.find(".mzML") != -1) {
                        fileNamePkl.replace(fileNamePkl.find(".mzML"),5,outputTagPkl);
                    }
                    outputNamePkl += "/";
                    outputNamePkl += fileNamePkl;
                    DEBUG_MSG("Saving pickle bin: " << outputNamePkl);
                    filtered_spectrum_iter->SaveSpecSet(outputNamePkl.c_str());
                }
                current_file++;
            }
        };

        if (m_params.exists("KL_OUTPUT")) {
            vector<SpecSet>::iterator spectrum_iter;
            for(spectrum_iter = m_inputSpectra->begin();spectrum_iter != m_inputSpectra->end(); spectrum_iter++) {
                std::ofstream out_summary;
                string fileName = (spectrum_iter->begin())->fileName;
                out_summary.open(m_params.getValue("KL_OUTPUT").c_str(),std::ios_base::app);
                vector<specnets::Spectrum>::iterator spec_iter;
                for (spec_iter = spectrum_iter->begin(); spec_iter != spectrum_iter->end(); spec_iter++) {
                    if (spec_iter->msLevel == 2) {

                        unsigned short charge = (unsigned short) spec_iter->parentCharge;
                        float parent_mass = spec_iter->parentMass;
                        float precursor_mass = (float) (parent_mass + (charge - 1) * specnets::AAJumps::massHion) / charge;

                        out_summary << fileName << "\t";
                        out_summary << spec_iter->scan << "\t";
                        out_summary << spec_iter->precursor_scan << "\t";
                        out_summary << precursor_mass << "\t";
                        out_summary << spec_iter->precursor_kl << "\t";
                        out_summary << spec_iter->size() << "\n";

                    }
                }
            }
        }
        return true;
    }

    // -------------------------------------------------------------------------
    bool ExecKLFilter::saveInputData(std::vector<std::string> & filenames) {
        return true;
    }

    // -------------------------------------------------------------------------
    bool ExecKLFilter::loadOutputData(void) {
        return true;
    }

    // -------------------------------------------------------------------------
    vector<ExecBase*> const & ExecKLFilter::split(int numSplit) {
        return m_subModules;
    }

    // -------------------------------------------------------------------------
    bool ExecKLFilter::merge(void) {
        return true;
    }

    // -------------------------------------------------------------------------
    bool ExecKLFilter::validateParams(std::string & error) {

        m_isValid = false;
        VALIDATE_PARAM_EXIST("ISOTOPE_TABLE");
        VALIDATE_PARAM_EXIST("PEAK_TOLERANCE");
//        VALIDATE_PARAM_EXIST("OUTPUT_FOLDER");
        m_isValid = true;
        return true;

    }


} /* namespace specnets */
