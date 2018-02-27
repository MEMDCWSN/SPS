#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <cstring>
#include <list>
#include <vector>

#include "spectrum.h"
#include "inputParams.h"
#include "mzxml.h"
#include "utils.h"
#include "abruijn.h"

using namespace std;
using namespace specnets;

const char* in_pklbin = "INPUT_PKLBIN (*:*)";
const char* in_mzxml = "INPUT_MZXML (*:*)";
const char* in_mgf = "INPUT_MGF (*:*)";
const char* in_pkl = "INPUT_PKL (*:*)";
const char* in_pkl_dir = "INPUT_PKL_DIR";
const char* in_mic_dir = "INPUT_MIC_DIR";
const char* in_ms2 = "INPUT_MS2 (*:*)";
const char* in_prms = "INPUT_PRMS (*:*)";
const char* in_prmsv3 = "INPUT_PRMSV3 (*:*)";
const char* in_abinfo = "INPUT_ABINFO (*:*)";
const char* in_overlaps = "INPUT_OVERLAPS (*_midx.pklbin) (*:*)";
const char* in_prot_match = "INPUT_MATCHED_PROTS (*_mp.bin) (*:*)";

const char* activ = "ACTIVATION";
const char* fix_charge = "FIX_PRECURSOR_CHARGE";
const char* merge_prec = "MERGE_SPECS_SAME_PREC";
const char* charge_filt = "CHARGE_FILTER";
const char* merge_type = "MERGE_TYPE";
const char* peak_tol_ppm = "TOLERANCE_PEAK_PPM";
const char* intensity_filt = "INTENSITY_FILTER";
const char* merge_triples = "MERGE_TRIPLES";
const char* rank_filt = "RANK_FILTER";
const char* rank_filt_radius = "RANK_FILTER_RADIUS";
const char* conv_charge = "CONVERT_CHARGE";
const char* offset_peaks = "OFFSET_PEAKS_DA";
const char* out_pklbin = "OUTPUT_PKLBIN";
const char* out_mgf = "OUTPUT_MGF";
const char* out_pkl = "OUPUT_PKL";
const char* out_ms2 = "OUTPUT_MS2";
const char* out_abinfo = "OUTPUT_ABINFO";
const char* out_overlaps = "OUTPUT_OVERLAPS (*_midx.pklbin)";
const char* out_prot_match = "OUTPUT_MATCHED_PROTS (*_mp.bin)";

int main(int argc, char *argv[], char **envp)
{
  InputParams params;
  bool paramsOk;
  const char* inputParamFile;
  if (argc <= 1)
  {
    paramsOk = params.readParams("mergeConvert.params");
    inputParamFile = "mergeConvert.params";
  }
  else
  {
    paramsOk = params.readParams(argv[1]);
    inputParamFile = argv[1];
  }
  if (!paramsOk)
  {
    cerr << "Error opening parameters file " << inputParamFile << "\n";
    return -1;
  }

  SpecSet allSpecs;
  abinfo_t all_abinfo;
  vector<int> in_spec_sizes;
  vector<int> local_spec_sizes;
  vector<string> all_in_filenames;
  SpecSet overlaps;
  vector<vector<int> > prot_match;

  ostringstream outs;
  outs << inputParamFile << ".info";

  FILE* output = fopen(outs.str().c_str(), "wb");

  fprintf(output, "Spectra Merge Order:\n<filename>\t<# spectra>\n");
  // if (! Load_binArray(argv[5], prot_match1)) cerr << "ERROR: Failed to load " << argv[5] << "\n";

  int idxUse = 0;
  int idxOv = 0;
  int idxPr = 0;
  int prevSize = 0;
  vector<string> filenames;

  if (params.paramPresent(in_pklbin))
  {
    cout << "Loading pklbin spectra:\n";
    splitText(params.getValue(in_pklbin), filenames, ":");

    for (int i = 0; i < filenames.size(); i++)
    {
      SpecSet in_specs;
      if (!in_specs.loadPklBin(filenames[i].c_str()))
      {
        cerr << "ERROR: Failed to load <" << filenames[i] << ">\n";
        return -1;
      }
      all_in_filenames.push_back(filenames[i]);
      in_spec_sizes.push_back(prevSize);
      prevSize += in_specs.size();
      local_spec_sizes.push_back(in_specs.size());
      cout << in_specs.size() << " from " << filenames[i] << ", ";
      cout.flush();
      fprintf(output, "%s\t%d\n", filenames[i].c_str(), in_specs.size());
      allSpecs.resize(allSpecs.size() + in_specs.size());
      for (int j = 0; j < in_specs.size(); j++)
      {
        allSpecs[idxUse] = in_specs[j];
        idxUse++;
      }
    }
    cout << "\n\n";
  }

  vector<short> msLevel;
  if (params.paramPresent(in_mzxml))
  {
    cout << "Loading mzxml spectra:\n";
    splitText(params.getValue(in_mzxml), filenames, ":");
    for (int i = 0; i < filenames.size(); i++)
    {
      SpecSet in_specs;
      const char* mzxmlFile = filenames[i].c_str();
      if (LoadMzxml(mzxmlFile, in_specs, &msLevel, 2) == 0)
      {
        cerr << "ERROR: Failed to load <" << filenames[i] << ">\n";
        return -1;
      }
      all_in_filenames.push_back(filenames[i]);
      in_spec_sizes.push_back(prevSize);
      prevSize += in_specs.size();
      local_spec_sizes.push_back(in_specs.size());
      cout << in_specs.size() << " from " << filenames[i] << ", ";
      cout.flush();
      fprintf(output, "%s\t%d\n", filenames[i].c_str(), in_specs.size());
      allSpecs.resize(allSpecs.size() + in_specs.size());
      for (int j = 0; j < in_specs.size(); j++)
      {
        allSpecs[idxUse] = in_specs[j];
        idxUse++;
      }
    }
    cout << "\n\n";
  }

  filenames.resize(0);
  if (params.paramPresent(in_mgf))
  {
    cout << "Loading mgf spectra:\n";
    splitText(params.getValue(in_mgf), filenames, ":");
    for (int i = 0; i < filenames.size(); i++)
    {
      SpecSet in_specs;
      if (!in_specs.LoadSpecSet_mgf(filenames[i].c_str()))
      {
        cerr << "ERROR: Failed to load <" << filenames[i] << ">\n";
        return -1;
      }
      all_in_filenames.push_back(filenames[i]);
      in_spec_sizes.push_back(prevSize);
      prevSize += in_specs.size();
      local_spec_sizes.push_back(in_specs.size());
      cout << in_specs.size() << " from " << filenames[i] << ", ";
      cout.flush();
      fprintf(output, "%s\t%d\n", filenames[i].c_str(), in_specs.size());
      allSpecs.resize(allSpecs.size() + in_specs.size());
      for (int j = 0; j < in_specs.size(); j++)
      {
        allSpecs[idxUse] = in_specs[j];
        idxUse++;
      }
    }
    cout << "\n\n";
  }

  filenames.resize(0);
  if (params.paramPresent(in_pkl))
  {
    cout << "Loading pkl spectra:\n";
    splitText(params.getValue(in_pkl), filenames, ":");
    for (int i = 0; i < filenames.size(); i++)
    {
      SpecSet in_specs;
      if (!in_specs.LoadSpecSet_pkl(filenames[i].c_str()))
      {
        cerr << "ERROR: Failed to load <" << filenames[i] << ">\n";
        return -1;
      }
      all_in_filenames.push_back(filenames[i]);
      in_spec_sizes.push_back(prevSize);
      prevSize += in_specs.size();
      local_spec_sizes.push_back(in_specs.size());
      cout << in_specs.size() << " from " << filenames[i] << ", ";
      cout.flush();
      fprintf(output, "%s\t%d\n", filenames[i].c_str(), in_specs.size());
      allSpecs.resize(allSpecs.size() + in_specs.size());
      for (int j = 0; j < in_specs.size(); j++)
      {
        allSpecs[idxUse] = in_specs[j];
        idxUse++;
      }
    }
    cout << "\n\n";
  }

  if (params.paramPresent(in_pkl_dir))
  {
    list<string> dir_files;
    cout << "Loading pkl spectra from directory "
        << params.getValue(in_pkl_dir) << ":\n";
    if (!getdir(params.getValue(in_pkl_dir), dir_files))
    {
      cerr << "ERROR: Failed to load <" << params.getValue(in_pkl_dir) << ">\n";
      return -1;
    }
    int specCount = 0;
    for (list<string>::iterator dirIt = dir_files.begin(); dirIt
        != dir_files.end(); dirIt++)
    {
      SpecSet in_specs;
      string filename = *dirIt;
      if (filename.find(".pkl") == string::npos)
        continue;
      string root = params.getValue(in_pkl_dir);
      root.append("/");
      root.append(filename);
      if (!in_specs.LoadSpecSet_pkl_mic(root.c_str()))
      {
        cerr << "ERROR: Failed to load <" << root << ">\n";
        return -1;
      }
      fprintf(output, "%s\t%d\n", root.c_str(), in_specs.size());
      allSpecs.resize(allSpecs.size() + in_specs.size());
      specCount += in_specs.size();
      for (int j = 0; j < in_specs.size(); j++)
      {
        allSpecs[idxUse] = in_specs[j];
        idxUse++;
      }
    }
    all_in_filenames.push_back(in_pkl_dir);
    in_spec_sizes.push_back(prevSize);
    prevSize += specCount;
    local_spec_sizes.push_back(specCount);
    cout << specCount << " from .pkl files in " << params.getValue(in_pkl_dir)
        << "\n\n";
    cout.flush();
  }

  if (params.paramPresent(in_mic_dir))
  {
    list<string> dir_files;
    cout << "Loading mic spectra from directory "
        << params.getValue(in_mic_dir) << ":\n";
    if (!getdir(params.getValue(in_mic_dir), dir_files))
    {
      cerr << "ERROR: Failed to load <" << params.getValue(in_mic_dir) << ">\n";
      return -1;
    }
    int specCount = 0;
    for (list<string>::iterator dirIt = dir_files.begin(); dirIt
        != dir_files.end(); dirIt++)
    {
      SpecSet in_specs;
      string filename = *dirIt;
      if (filename.find(".mic") == string::npos)
        continue;
      string root = params.getValue(in_mic_dir);
      root.append("/");
      root.append(filename);
      if (!in_specs.LoadSpecSet_pkl_mic(root.c_str()))
      {
        cerr << "ERROR: Failed to load <" << root << ">\n";
        return -1;
      }
      fprintf(output, "%s\t%d\n", root.c_str(), in_specs.size());
      allSpecs.resize(allSpecs.size() + in_specs.size());
      specCount += in_specs.size();
      for (int j = 0; j < in_specs.size(); j++)
      {
        allSpecs[idxUse] = in_specs[j];
        idxUse++;
      }
    }
    all_in_filenames.push_back(in_mic_dir);
    in_spec_sizes.push_back(prevSize);
    prevSize += specCount;
    local_spec_sizes.push_back(specCount);
    cout << specCount << " from .mic files in " << params.getValue(in_mic_dir)
        << "\n\n";
    cout.flush();
  }

  filenames.resize(0);
  if (params.paramPresent(in_ms2))
  {
    cout << "Loading ms2 spectra:\n";
    splitText(params.getValue(in_ms2), filenames, ":");
    for (int i = 0; i < filenames.size(); i++)
    {
      SpecSet in_specs;
      if (!in_specs.LoadSpecSet_ms2(filenames[i].c_str()))
      {
        cerr << "ERROR: Failed to load <" << filenames[i] << ">\n";
        return -1;
      }
      all_in_filenames.push_back(filenames[i]);
      in_spec_sizes.push_back(prevSize);
      prevSize += in_specs.size();
      local_spec_sizes.push_back(in_specs.size());
      cout << in_specs.size() << " from " << filenames[i] << ", ";
      cout.flush();
      fprintf(output, "%s\t%d\n", filenames[i].c_str(), in_specs.size());
      allSpecs.resize(allSpecs.size() + in_specs.size());
      for (int j = 0; j < in_specs.size(); j++)
      {
        allSpecs[idxUse] = in_specs[j];
        idxUse++;
      }
    }
    cout << "\n\n";
  }

  filenames.resize(0);
  if (params.paramPresent(in_prms))
  {
    cout << "Loading prms spectra:\n";
    splitText(params.getValue(in_prms), filenames, ":");
    for (int i = 0; i < filenames.size(); i++)
    {
      SpecSet in_specs;
      if (!in_specs.LoadSpecSet_prmsv3(filenames[i].c_str()))
      {
        cerr << "ERROR: Failed to load <" << filenames[i] << ">\n";
        return -1;
      }
      all_in_filenames.push_back(filenames[i]);
      in_spec_sizes.push_back(prevSize);
      prevSize += in_specs.size();
      local_spec_sizes.push_back(in_specs.size());
      cout << in_specs.size() << " from " << filenames[i] << ", ";
      cout.flush();
      fprintf(output, "%s\t%d\n", filenames[i].c_str(), in_specs.size());
      allSpecs.resize(allSpecs.size() + in_specs.size());
      for (int j = 0; j < in_specs.size(); j++)
      {
        allSpecs[idxUse] = in_specs[j];
        idxUse++;
      }
    }
    cout << "\n\n";
  }

  filenames.resize(0);
  if (params.paramPresent(in_prmsv3))
  {
    cout << "Loading prmsv3 spectra:\n";
    splitText(params.getValue(in_prmsv3), filenames, ":");
    for (int i = 0; i < filenames.size(); i++)
    {
      SpecSet in_specs;
      if (!in_specs.LoadSpecSet_prmsv3(filenames[i].c_str()))
      {
        cerr << "ERROR: Failed to load <" << filenames[i] << ">\n";
        return -1;
      }
      all_in_filenames.push_back(filenames[i]);
      in_spec_sizes.push_back(prevSize);
      prevSize += in_specs.size();
      local_spec_sizes.push_back(in_specs.size());
      cout << in_specs.size() << " from " << filenames[i] << ", ";
      cout.flush();
      fprintf(output, "%s\t%d\n", filenames[i].c_str(), in_specs.size());
      allSpecs.resize(allSpecs.size() + in_specs.size());
      for (int j = 0; j < in_specs.size(); j++)
      {
        allSpecs[idxUse] = in_specs[j];
        idxUse++;
      }
    }
    cout << "\n\n";
  }

  filenames.resize(0);
  if (params.paramPresent(in_abinfo))
  {
    fprintf(output, "\nAbinfo Merge Order:\n<filename>\t<# ab components>\n");
    cout << "Loading abinfo:\n";
    splitText(params.getValue(in_abinfo), filenames, ":");
    for (int i = 0; i < filenames.size(); i++)
    {
      if (i >= in_spec_sizes.size())
      {
        cerr
            << "ERROR: Not enough .pklbin files loaded to determine spectrum index offsets for "
            << filenames[i] << "\n";
        return -1;
      }
      abinfo_t in_ab;
      if (!Load_abinfo(filenames[i].c_str(), in_ab))
      {
        cerr << "ERROR: Failed to load <" << filenames[i] << ">\n";
        return -1;
      }
/*
      int idxCheck = 2;
      cout << "\nContig " << idxCheck << ":\n";
      for (int i1 = 0; i1 < in_ab[idxCheck].second.size(); i1++)
      {
        cout << i1 << " :";
        for (int i2 = 0; i2 < in_ab[idxCheck].second[i1].first.size(); i2++)
        {
          cout << "(" << in_ab[idxCheck].second[i1].first[i2] << " , "
              << in_ab[idxCheck].second[i1].second[i2] << "), ";
        }
        cout << "\n";
      }
      idxCheck = 3;
      cout << "\nContig " << idxCheck << ":\n";
      for (int i1 = 0; i1 < in_ab[idxCheck].second.size(); i1++)
      {
        cout << i1 << " :";
        for (int i2 = 0; i2 < in_ab[idxCheck].second[i1].first.size(); i2++)
        {
          cout << "(" << in_ab[idxCheck].second[i1].first[i2] << " , "
              << in_ab[idxCheck].second[i1].second[i2] << "), ";
        }
        cout << "\n";
      }
      idxCheck = 4;
      cout << "\nContig " << idxCheck << ":\n";
      for (int i1 = 0; i1 < in_ab[idxCheck].second.size(); i1++)
      {
        cout << i1 << " :";
        for (int i2 = 0; i2 < in_ab[idxCheck].second[i1].first.size(); i2++)
        {
          cout << "(" << in_ab[idxCheck].second[i1].first[i2] << " , "
              << in_ab[idxCheck].second[i1].second[i2] << "), ";
        }
        cout << "\n";
      }

      for (idxCheck = 0; idxCheck < in_ab.size(); idxCheck++)
      {
        bool found = false;
        for (int i1 = 0; i1 < in_ab[idxCheck].first.first.size(); i1++)
        {
          if (in_ab[idxCheck].first.first[i1] == 336
              || in_ab[idxCheck].first.first[i1] == 335)
          {
            found = true;
            break;
          }
        }
        if (!found)
        {
          continue;
        }
        cout << "\nContig " << idxCheck << ":\n";
        for (int i1 = 0; i1 < in_ab[idxCheck].second.size(); i1++)
        {
          cout << i1 << " :";
          for (int i2 = 0; i2 < in_ab[idxCheck].second[i1].first.size(); i2++)
          {
            cout << "(" << in_ab[idxCheck].second[i1].first[i2] << " , "
                << in_ab[idxCheck].second[i1].second[i2] << "), ";
          }
          cout << "\n";
        }
      }
*/
      cout << in_ab.size() << " from " << filenames[i] << ", ";
      fprintf(output,
              "%s\t%lu\n",
              filenames[i].c_str(),
              (unsigned long)in_ab.size());
      int size_in = all_abinfo.size();
      for (abinfo_t::iterator compit = in_ab.begin(); compit != in_ab.end(); compit++)
      {
        pair<pair<vector<int> , vector<int> > , vector<pair<vector<int> ,
            vector<double> > > > second = compit->second;
        for (int j = 0; j < second.first.first.size(); j++)
        {
          second.first.first[j] += in_spec_sizes[i];
        }
        for (int j = 0; j < second.second.size(); j++)
        {
          pair<vector<int> , vector<double> > vert = second.second[j];
          for (int k = 0; k < vert.first.size(); k++)
            second.second[j].first[k] += in_spec_sizes[i];
        }
        all_abinfo[compit->first + size_in] = second;
      }
    }
    cout << "\n\n";
  }

  filenames.resize(0);
  if (params.paramPresent(in_overlaps))
  {
    fprintf(output, "\nOverlaps Merge Order:\n<filename>\t<# overlaps>\n");
    cout << "Loading overlaps:\n";
    splitText(params.getValue(in_overlaps), filenames, ":");
    for (int i = 0; i < filenames.size(); i++)
    {
      SpecSet in_specs;
      if (!in_specs.loadPklBin(filenames[i].c_str()))
      {
        cerr << "ERROR: Failed to load <" << filenames[i] << ">\n";
        return -1;
      }
      cout << in_specs.size() << " from " << filenames[i] << ", ";
      fprintf(output, "%s\t%d\n", filenames[i].c_str(), in_specs.size());
      overlaps.resize(overlaps.size() + in_specs.size());
      for (int j = 0; j < in_specs.size(); j++)
      {
        overlaps[idxOv] = in_specs[j];
        idxOv++;
      }
    }
    cout << "\n\n";
  }

  filenames.resize(0);
  if (params.paramPresent(in_prot_match))
  {
    fprintf(output,
            "\nMatched Protein Order:\n<filename>\t<# matched spectra>\n");
    cout << "Loading matched proteins:\n";
    splitText(params.getValue(in_prot_match), filenames, ":");
    for (int i = 0; i < filenames.size(); i++)
    {
      vector<vector<int> > prot_m;
      if (!Load_binArray(filenames[i].c_str(), prot_m))
      {
        cerr << "ERROR: Failed to load <" << filenames[i] << ">\n";
        return -1;
      }
      cout << prot_m.size() << " from " << filenames[i] << ", ";
      fprintf(output,
              "%s\t%lu\n",
              filenames[i].c_str(),
              (unsigned long)prot_m.size());
      prot_match.resize(prot_match.size() + prot_m.size());
      for (int j = 0; j < prot_m.size(); j++)
      {
        prot_match[idxPr] = prot_m[j];
        idxPr++;
      }
    }
    cout << "\n\n";
  }

  if (params.paramPresent(fix_charge))
  {
    for (int i = 0; i < allSpecs.size(); i++)
    {
      if (allSpecs[i].parentCharge == 0 && allSpecs[i].parentMass
          + AAJumps::massMH < allSpecs[i][allSpecs[i].size() - 1][0])
      {

        while (allSpecs[i].parentMass + AAJumps::massMH
            < allSpecs[i][allSpecs[i].size() - 1][0])
        {
          allSpecs[i].parentCharge += 1;
          allSpecs[i].parentMass = (allSpecs[i].parentMass
              * (float)allSpecs[i].parentCharge) - (AAJumps::massHion
              * ((float)allSpecs[i].parentCharge - 1.0));
        }

        cout << "assigned charge " << allSpecs[i].parentCharge
            << " to spectrum " << i << "\n";

      }
    }

  }

  if (params.paramPresent(charge_filt))
  {
    vector<string> chargeRanges;
    splitText(params.getValue(charge_filt), chargeRanges, ":");
    set<short> chargeRangeVals;
    cout << "Extracting parent charges: ";
    for (int i = 0; i < chargeRanges.size(); i++)
    {
      if (chargeRanges[i].length() == 0)
      {
        continue;
      }
      if (!getRanges(chargeRanges[i].c_str(), chargeRangeVals))
      {
        return -1;
      }
      cout << chargeRanges[i] << " from " << all_in_filenames[i] << ", ";
      for (int j = in_spec_sizes[i]; j < in_spec_sizes[i] + local_spec_sizes[i]; j++)
      {
        if (chargeRangeVals.count(allSpecs[j].parentCharge) == 0)
        {
          allSpecs[j].resize(0);
        }
      }
    }
    cout << "\n\n";
  }

  if (params.paramPresent(offset_peaks))
  {
    float peakOffset = params.getValueDouble(offset_peaks);

    cout << "Adding " << peakOffset << " to every peak\n";
    for (int i = 0; i < allSpecs.size(); i++)
    {
      for (int j = 0; j < allSpecs[i].size(); j++)
      {
        allSpecs[i][j][0] += peakOffset;
      }
    }
  }

  if (params.paramPresent(out_mgf) || params.paramPresent(out_pklbin)
      || params.paramPresent(out_pkl) || params.paramPresent(out_ms2))
  {
    fprintf(output, "\nOutput Spectra:\n<filename>\t<# spectra>\n");
  }

  if (params.paramPresent(out_mgf))
  {
    cout << "Saving " << allSpecs.size() << " spectra to <"
        << params.getValue(out_mgf) << "> ...\n";
    const char* activation;
    if (params.paramPresent(activ))
    {
      activation = params.getValue(activ);
    }
    else
    {
      activation = 0;
    }
    if (!allSpecs.SaveSpecSet_mgf(params.getValue(out_mgf), activation))
    {
      cerr << "ERROR: Failed to save to <" << params.getValue(out_mgf) << ">\n";
      return -1;
    }
    fprintf(output, "%s\t%d\n", params.getValue(out_mgf), allSpecs.size());
  }

  if (params.paramPresent(out_pklbin))
  {
    cout << "Saving " << allSpecs.size() << " spectra to <"
        << params.getValue(out_pklbin) << "> ...\n";
    if (!allSpecs.SaveSpecSet_pklbin(params.getValue(out_pklbin)))
    {
      cerr << "ERROR: Failed to save to <" << params.getValue(out_pklbin)
          << ">\n";
      return -1;
    }
    fprintf(output, "%s\t%d\n", params.getValue(out_pklbin), allSpecs.size());
  }

  if (params.paramPresent(out_pkl))
  {
    cout << "Saving " << allSpecs.size() << " spectra to <"
        << params.getValue(out_pkl) << "> ...\n";
    if (!allSpecs.SaveSpecSet_pkl(params.getValue(out_pkl)))
    {
      cerr << "ERROR: Failed to save to <" << params.getValue(out_pkl) << ">\n";
      return -1;
    }
    fprintf(output, "%s\t%d\n", params.getValue(out_pkl), allSpecs.size());
  }

  if (params.paramPresent(out_ms2))
  {
    cout << "Saving " << allSpecs.size() << " spectra to <"
        << params.getValue(out_ms2) << "> ...\n";
    if (!allSpecs.SaveSpecSet_ms2(params.getValue(out_ms2)))
    {
      cerr << "ERROR: Failed to save to <" << params.getValue(out_ms2) << ">\n";
      return -1;
    }
    fprintf(output, "%s\t%d\n", params.getValue(out_ms2), allSpecs.size());
  }

  if (params.paramPresent(out_abinfo))
  {
    fprintf(output, "\nOutput Abinfo:\n<filename>\t<# ab components>\n");
    cout << "Saving " << all_abinfo.size() << " ab components to <"
        << params.getValue(out_abinfo) << "> ...\n";
    if (!Save_abinfo_v1_0(params.getValue(out_abinfo), all_abinfo))
    {
      cerr << "ERROR: Failed to save to <" << params.getValue(out_abinfo)
          << ">\n";
      return -1;
    }
    fprintf(output,
            "%s\t%lu\n",
            params.getValue(out_abinfo),
            (unsigned long)all_abinfo.size());
  }

  if (params.paramPresent(out_overlaps))
  {
    fprintf(output, "\nOutput Overlaps:\n<filename>\t<# overlaps>\n");
    cout << "Saving " << overlaps.size() << " overlaps to <"
        << params.getValue(out_overlaps) << "> ...\n";
    if (!overlaps.SaveSpecSet_pklbin(params.getValue(out_overlaps)))
    {
      cerr << "ERROR: Failed to save to <" << params.getValue(out_overlaps)
          << ">\n";
      return -1;
    }
    fprintf(output, "%s\t%d\n", params.getValue(out_overlaps), overlaps.size());
  }

  if (params.paramPresent(out_prot_match))
  {
    fprintf(output,
            "\nOutput Matched Proteins:\n<filename>\t<# spectra matched>\n");
    cout << "Saving " << prot_match.size() << " protein matched spectra to <"
        << params.getValue(out_prot_match) << "> ...\n";
    if (!Save_binArray(params.getValue(out_prot_match), prot_match))
    {
      cerr << "ERROR: Failed to save to <" << params.getValue(out_prot_match)
          << ">\n";
      return -1;
    }
    fprintf(output,
            "%s\t%lu\n",
            params.getValue(out_prot_match),
            (unsigned long)prot_match.size());
  }
  fclose(output);
  return 0;
}

