//Module Includes
#include "SpectralLibrary.h"
#include "utils.h"
#include "projectionutils.h"

// Header Include
#include "FirstCorrectTags.h"
#include "ExecFdrPeptide.h"

// System Include
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

using namespace specnets;
using namespace std;


namespace specnets
{

  FirstCorrectTags::FirstCorrectTags(void){
    DEBUG_MSG("1 ------------- START!")

    m_name = "FirstCorrectTags";
    m_type = "FirstCorrectTags";
  }


  FirstCorrectTags::FirstCorrectTags(const ParameterList & inputParams){
    DEBUG_MSG("2 ------------- START!")

    m_name = "FirstCorrectTags";
    m_type = "FirstCorrectTags";
  }


  FirstCorrectTags::~FirstCorrectTags(void){
      DEBUG_MSG("3 ------------- START!")
  }


  ExecBase * FirstCorrectTags::clone(const ParameterList & inputParams) const{
    DEBUG_MSG("4 ------------- START!")

    return new FirstCorrectTags(inputParams);
  }

  bool FirstCorrectTags::invoke(void){
    DEBUG_MSG("5 ------------- START!")
    return true;
  }

  // original first correct tag version
  void get_first_correct_ranks(SpectralLibrary spec_lib)
  {
      //DW
        float peakTol = 0.1;//0.02;
        float comparTol = 0.1;
        int tagLen = 3;
        int gap = 0;
        int tagNum = 20;        //D get top 20 tags from each search spectrum
        int TAG_COMMON = 1;
        bool TAG_FILTERING = true;

        AAJumps jumps(1);
        //D
        FILE *fo = fopen("/media/duong/NewVolume/SPS/code/PRM_first_correct_tags", "w");
        //DFILE *fo = fopen("/media/duong/NewVolume/SPS/code/MSMS_first_correct_tags", "w");

        for (int i = 0; i < spec_lib.size(); i++)
        {
            //D if (spec_lib[i].scan != 557) continue;
            for (int j = 0; j < spec_lib[i].size(); j++)
                cout << spec_lib[i][j][0] << " " << spec_lib[i][j][1] << "; " << endl;
            int j, k, l;
            string annotation = spec_lib[i].psmList.front()->m_annotation;

            /* can uncomment if only considering unmod spectra
            for (j = 0; j < annotation.length(); j++)
                if ((annotation[j] < 'A') or (annotation[j] > 'Z'))
                    break;
            if (j < annotation.length())
                continue;
            */

            list<Tag> tags2;
            //D cout << endl << "spec " << i << " scan " << (*spec_lib)[i].scan << " ";
            //D cout << ExtractCorrectTags((*spec_lib)[i], tags, comparTol, tagLen, false) << endl;
            ExtractCorrectTags(spec_lib[i], tags2, comparTol, tagLen, false, true);     //D no filtering for extracting correct tags
            //D cout << tags.size() << endl;

            list< pair<string, float> > temp2;       //D <tag sequence, tag prefix>
            for (list<Tag>::iterator it = tags2.begin(); it != tags2.end(); it++)
            {
                string strSequence;
                for (int c = 0; c < it->sequence.size(); c++)
                    strSequence = strSequence + jumps.aaLetters[it->sequence[c]];
                //D cout << strSequence.c_str() << " " << strSequence.length() << " " << it->score << endl;

                pair<string, float> p_temp;
                p_temp.first = strSequence;
                p_temp.second = it->flankingPrefix;
                temp2.push_back(p_temp);     //D no consider the score of the tag

                cout << p_temp.first << " " << p_temp.second << "; ";
            }
            cout << endl << "denovo begin" << endl;

            fprintf(fo, "%ld;\t", spec_lib[i].scan);
            list<Tag> tags;
            //D ExtractTags(search_lib[i], tags, peakTol, tagLen, gap, tagNum);
            vector<unsigned int> *f_heap_times = &(ExtractDenovoTags(spec_lib[i], tags, peakTol, tagLen, tagNum, true));      //D do filtering for search spectra
            //D cout << tags.size() << endl;

            list< pair<string, float> > temp;
            long f_rank = -1;       //D rank of first correct tag
            j = 0;
            for (list<Tag>::iterator it = tags.begin(); it != tags.end(); it++)
            {
                string strSequence;
                for (int c = 0; c < it->sequence.size(); c++)
                    strSequence = strSequence + jumps.aaLetters[it->sequence[c]];
                pair<string, float> p_temp;
                p_temp.first = strSequence;
                p_temp.second = it->flankingPrefix;
                temp.push_back(p_temp);

                cout << strSequence.c_str() << " " << strSequence.length() << " " << it->score << ";";
                fprintf(fo, "%s %lf %lf ", (p_temp.first).c_str(), p_temp.second, it->score);

                j++;        //D current rank in the denovo tag set
                if (f_rank < 0)
                    for (list< pair<string, float> >::iterator it2 = temp2.begin(); it2 != temp2.end(); it2++)
                        if (((p_temp.first).compare((*it2).first) == 0) and (fabs(p_temp.second-(*it2).second) < 0.1))
                            f_rank = j;
            }

            cout << endl;
            cout << spec_lib[i].scan << ";" << spec_lib[i].psmList.front()->m_annotation << endl;

            fprintf(fo, ";\t%ld", f_rank);
            fprintf(fo, "\n");
        }
        fclose(fo);
  }



//  void get_first_correct_ranks(SpecSet spec_lib, bool is_PRM_specs, char* FILE_OUT)
//  {
//      //DW
//        float peakTol = 0.1;//0.02;
//        float comparTol = 3.0;//0.1;
//        int tagLen = 3;
//        int gap = 0;
//        int tagNum = 20;        //D get top 20 tags from each search spectrum
//        int TAG_COMMON = 1;
//        bool TAG_FILTERING = true;
//
//        AAJumps jumps(1);
//        //D
//        FILE *fo = fopen(FILE_OUT, "w");
//
//        for (int i = 0; i < spec_lib.size(); i++)
//        {
//            //D if (spec_lib[i].scan != 557) continue;
//            for (int j = 0; j < spec_lib[i].size(); j++)
//                cout << spec_lib[i][j][0] << " " << spec_lib[i][j][1] << "; " << endl;
//            int j, k, l;
//
//            cout << endl << "denovo begin" << endl;
//
//            list<Tag> tags;
//            //D ExtractTags(search_lib[i], tags, peakTol, tagLen, gap, tagNum);
//            vector<unsigned int> *f_heap_times = &(ExtractDenovoTags(spec_lib[i], tags, peakTol, tagLen, tagNum, true));      //D do filtering for search spectra
//            //D cout << tags.size() << endl;
//            fprintf(fo, "%ld;%ld;", spec_lib[i].scan, tags.size());
//
//            list< pair<string, float> > temp;
//            long f_rank = -1;       //D rank of first correct tag
//            j = 0;
//            for (list<Tag>::iterator it = tags.begin(); it != tags.end(); it++)
//            {
//                string strSequence;
//                for (int c = 0; c < it->sequence.size(); c++)
//                    strSequence = strSequence + jumps.aaLetters[it->sequence[c]];
//
//                pair<string, float> p_temp;
//                p_temp.first = strSequence;
//                p_temp.second = it->flankingPrefix;
//                if (is_PRM_specs)
//                    p_temp.second = p_temp.second + AAJumps::massHion;      //D shifting by 1 for PRM tags
//                temp.push_back(p_temp);
//
//                cout << strSequence.c_str() << " " << strSequence.length() << " " << it->score << ";";
//                fprintf(fo, "%s %lf %lf,", (p_temp.first).c_str(), p_temp.second, it->score);
//            }
//
//            cout << endl;
//            //D cout << spec_lib[i].scan << ";" << spec_lib[i].psmList.front()->m_annotation << endl;
//
//            fprintf(fo, ";%ld", f_rank);
//            fprintf(fo, "\n");
//        }
//        fclose(fo);
//  }


  bool FirstCorrectTags::loadInputData(void){
    DEBUG_MSG("6 ------------- START!")

    DEBUG_MSG("bool FirstCorrectTags::loadInputData(void) - START!")

    //Loading the PSMs
    PeptideSpectrumMatchSetSpectralLibraryLoader m_training_psms;

    //D string annotation_file_name = "/media/duong/NewVolume/SPS/code/fix_MSGFDB-94bc69f5-group_by_spectrum_old-main.tsv";
    string annotation_file_name = "/media/duong/NewVolume/SPS/code/fix_MSPLIT_NEW-a859c5af-group_by_spectrum_old-main44.tsv";
    PeptideSpectrumMatchSetSpectralLibraryLoader psm_set;
    psm_set.loadMSGFDBResultsFile(annotation_file_name.c_str());

    for(int i = 0; i < psm_set.size(); i++){
        m_training_psms.push_back(psm_set[i]);
    }

    //Grouping PSMs per file
    map<string, PeptideSpectrumMatchSet> file_psm_map;
    for(int psm_idx = 0; psm_idx < m_training_psms.size(); psm_idx++){
        file_psm_map[m_training_psms[psm_idx]->m_spectrumFile].push_back(m_training_psms[psm_idx]);
    }

    int real_total_spectra = 0;

    //Determining file type by file extension, only supporting plkbin and mzxml
    string spectra_file_name = "/media/duong/NewVolume/SPS/code/PRM_b1925_293T_proteinID_05A_QE3_122212.mzXML";
    //Dstring spectra_file_name = "/media/duong/NewVolume/SPS/code/MSMS_b1925_293T_proteinID_05A_QE3_122212.mzXML";
    int dotpos = spectra_file_name.find_last_of('.');
    int slashpos = spectra_file_name.find_last_of('/');
    string extension = spectra_file_name.substr(dotpos+1);
    string filename_trunc = spectra_file_name;
    if(slashpos != string::npos)
        filename_trunc = spectra_file_name.substr(slashpos+1);
    DEBUG_MSG("filename_trunc"<<filename_trunc<<endl);

    if(strcmp(extension.c_str(), "mzXML") == 0 || strcmp(extension.c_str(), "mzxml") == 0){
        SpecSet tempspecs;
        cout<<"LOADING\t"<<spectra_file_name<<"\t as an mzxml file...\t";
        DEBUG_MSG("LOADING\t"<<spectra_file_name<<"\t as an mzxml file...\t" << endl << "tempspecs.size() = " << tempspecs.size() << endl)
        //LoadMzxml(spectra_file_name.c_str(), tempspecs, NULL, 0);
        //D tempspecs.loadPklBin((spectra_file_name.substr(0,dotpos+1)+"pklbin").c_str(), "/media/duong/NewVolume/SPS/code/MSGFDB-c37c798c-group_by_spectrum_old-main.tsv");
        //D tempspecs.loadPklBin((spectra_file_name.substr(0,dotpos+1)+"pklbin").c_str(), "/Users/dtn074/SPS/MSGFDB-c37c798c-group_by_spectrum_old-main.tsv");
        //D tempspecs.loadPklBin((spectra_file_name.substr(0,dotpos+1)+"pklbin").c_str(), "/media/duong/NewVolume/SPS/code/fix_MSGFDB-94bc69f5-group_by_spectrum_old-main.tsv");
        tempspecs.loadPklBin((spectra_file_name.substr(0,dotpos+1)+"pklbin").c_str(), "/media/duong/NewVolume/SPS/code/fix_MSPLIT_NEW-a859c5af-group_by_spectrum_old-main44.tsv");

        for(int i = 0; i < tempspecs.size(); i++){
            tempspecs[i].fileName = filename_trunc;
        }

        real_total_spectra += tempspecs.size();

        for(int psm_idx = 0; psm_idx < file_psm_map[filename_trunc].size(); psm_idx++){
            int scannum = file_psm_map[filename_trunc][psm_idx]->m_scanNum;
            DEBUG_MSG("BEFORE\t"<<m_training_library.size()<<"\t"<<scannum<<endl);
            m_training_library.insert(m_training_library.end(), tempspecs.begin() + scannum - 1, tempspecs.begin() + scannum);

            file_psm_map[filename_trunc][psm_idx]->m_spectrum = &m_training_library[m_training_library.size()-1];
            m_training_library[m_training_library.size()-1].psmList.push_back(file_psm_map[filename_trunc][psm_idx]);
            DEBUG_MSG("AFTER\t"<<m_training_library.size()<<endl);

        }

        cout<<"DONE with "<<m_training_library.size()<< " out of a total of "<<real_total_spectra<<" spectra"<<endl;
        DEBUG_MSG("DONE with "<<m_training_library.size()<< " out of a total of "<<real_total_spectra<<" spectra"<<endl)

        //D
        cout << "tempspecs.size() = " << tempspecs.size() << endl;
        //Dget_first_correct_ranks(tempspecs, false, "/media/duong/NewVolume/SPS/code/MSMS_tags");
        //Dget_first_correct_ranks(tempspecs, true, "/media/duong/NewVolume/SPS/code/PRM_tags");
    }

    get_first_correct_ranks(m_training_library);

    DEBUG_MSG("bool FirstCorrectTags::loadInputData(void) - END!")
    return true;
  }


  bool FirstCorrectTags::saveOutputData(void){
        DEBUG_MSG("7 ------------- START!")
        return true;
  }


  bool FirstCorrectTags::saveInputData(std::vector<std::string> & filenames){
    DEBUG_MSG("8 ------------- START!")
    return true;
  }

  bool FirstCorrectTags::loadOutputData(void){
    DEBUG_MSG("9 ------------- START!")
    return true;
  }

  std::vector<ExecBase *> const & FirstCorrectTags::split(int numSplit){
    DEBUG_MSG("10 ------------- START!")
    return m_subModules;
  }

  bool FirstCorrectTags::merge(void){
    DEBUG_MSG("11 ------------- START!")
    return true;
  }


  bool FirstCorrectTags::validateParams(std::string & error){
    DEBUG_MSG("12 ------------- START!")
    return true;
  }

  void FirstCorrectTags::load_aminoacid_masses(){
      DEBUG_MSG("13 ------------- START!")
  }

}
