/*
 * MappedSpectrum.cpp
 *
 *  Created on: Mar 3, 2011
 *      Author: aguthals
 */

#include "MappedSpectrum.h"

namespace specnets
{

  /**
   * Maps spectrum to peptide and/or proteins by filling in all MappedSpectrum and
   *   MappedPeak class variables
   * @param parentSpecs SpecSet containing spectrum
   * @param idx index of spectrum
   * @param _peptide specnets-style peptide annotation ("" if no annotation)
   * @param proteins sequence of target proteins
   * @param _reverse whether this spectrum should be reversed before annotating
   * @param inputIonTypes needed for Spectrum.annotate
   * @param peakTol Da peak tolerance
   * @return
   */
  void MappedSpectrum::mapProt(SpecSet* parentSpecs,
                               int idx,
                               string _peptide,
                               vector<string>* proteins,
                               bool _reverse,
                               MS2ScoringModel& inputIonTypes,
                               float peakTol)
  {
    ((Spectrum &)(*this)) = (*parentSpecs)[idx];
    this->setPeakTolerance(peakTol);
    reset();
    specIdx = idx;

    residueIdxs.clear();
    reversed = _reverse;
    peptide = _peptide;
    identified = peptide.length() > 0;

    modified = identified && (peptide.find('(') != string::npos);
    numPRMs = 0;
    numSRMs = 0;
    mapped = false;
    residueIdxs.clear();

    if (identified)
    {
      const char* id = peptide.c_str();
      char just_letters[strlen(id) + 1];

      //get rid of special characters and modifications so we can match to protein
      getPepSeq(id, &just_letters[0]);

      AAPeptideSeq = &just_letters[0];

      string::size_type start_pos;
      int residueIdx;
      list<int> starts;
      for (int j = 0; j < proteins->size(); j++)
      {
        if ((*proteins)[j].length() == 0)
        {
          continue;
        }
        //look for matching sequence in protein
        start_pos = (*proteins)[j].find(AAPeptideSeq);

        // keep matching this protein
        while (start_pos != string::npos)
        {
          mapped = true;

          residueIdx = ((int)start_pos);

          if (residueIdxs.count(j) > 0)
          {
            residueIdxs[j].push_back(residueIdx);
          }
          else
          {
            starts.clear();
            starts.push_back(residueIdx);
            residueIdxs[j] = starts;
          }

          start_pos = (*proteins)[j].find(AAPeptideSeq, start_pos + 1);
        }
      }
    }

    if (size() == 0)
    {
      return;
    }

    //cout << "mapping " << idx << ": " << identified << "\n"; cout.flush();

    if (reversed)
    {
      /*if (idx == 11366) {
       DEBUG_MSG("reversing " << idx << ":");
       output(cout);
       }*/
      reverse(0);
      /*if (idx == 11366) {
       DEBUG_MSG("now:");
       output(cout);
       }*/
    }

    addZPMpeaks(-1.0, 0, true);

    massToPeak.clear();
    for (int i = 0; i < size(); i++)
    {
      MZRange peakR(peakList[i][0], peakList[i][1], 0.01);
      massToPeak[peakR] = i;
    }

    map<int, int> BpeakToRes;
    map<int, int> YpeakToRes;

    if (identified)
    {

      AAJumps jumps(1);
      vector<float> prmMasses;
      jumps.getPRMMasses(peptide, prmMasses, 0, (vector<string>*)0, true);
      vector<float> unModPrmMasses;
      jumps.getPRMMasses(AAPeptideSeq,
                         unModPrmMasses,
                         0,
                         (vector<string>*)0,
                         true);

      if (prmMasses.size() != unModPrmMasses.size())
      {
        ERROR_MSG("Modified \"" << peptide << "\" (" << prmMasses.size() << ") and unmodified \"" << AAPeptideSeq << "\" (" << unModPrmMasses.size() << ") annotations have different sizes");
      }

      int numMasses = prmMasses.size();

      Spectrum prmMassSpec;
      prmMassSpec.resize(numMasses);
      prmMassSpec.parentMass = prmMasses[numMasses - 1] + AAJumps::massMH;
      for (int i = 0; i < numMasses; i++)
      {
        prmMassSpec[i][0] = prmMasses[i];
        prmMassSpec[i][1] = 1.0;
      }
      Spectrum unModPrmMassSpec;
      unModPrmMassSpec.resize(numMasses);
      unModPrmMassSpec.parentMass = prmMassSpec.parentMass;
      for (int i = 0; i < numMasses; i++)
      {
        unModPrmMassSpec[i][0] = unModPrmMasses[i];
        unModPrmMassSpec[i][1] = 1.0;
      }

      Spectrum srmMassSpec;
      srmMassSpec.parentMass = prmMassSpec.parentMass;
      prmMassSpec.reverse(0, &srmMassSpec);
      Spectrum unModSrmMassSpec;
      unModSrmMassSpec.parentMass = prmMassSpec.parentMass;
      unModPrmMassSpec.reverse(0, &unModSrmMassSpec);

      /*if (idx == 6213) {
       cout << idx << " - me: mod " << peptide << "(" << numMasses
       << ") and unmod " << AAPeptideSeq << "\n";
       output(cout);

       cout << "\nprms:\n";
       prmMassSpec.output(cout);

       cout << "\nunmod prms:\n";
       unModPrmMassSpec.output(cout);

       cout << "\nsrms:\n";
       srmMassSpec.output(cout);
       }*/

      MZRange peakB;
      MZRange peakY;
      int b_idx, b_idxS, unMod_b_idx, unMod_b_idxS;
      int y_idx, y_idxS, unMod_y_idx, unMod_y_idxS;

      for (int i = 0; i < size(); i++)
      {
        peakB.set(peakList[i][0], 1.0, peakTol);
        peakY.set(parentMass - AAJumps::massHion - peakList[i][0],
                  1.0,
                  peakTol);

        b_idx = prmMassSpec.findClosest(peakB.getMass());
        y_idx = prmMassSpec.findClosest(peakY.getMass());

        b_idxS = srmMassSpec.findClosest(peakB.getMass());
        y_idxS = srmMassSpec.findClosest(peakY.getMass());

        unMod_b_idx = unModPrmMassSpec.findClosest(peakB.getMass());
        unMod_y_idx = unModPrmMassSpec.findClosest(peakY.getMass());

        unMod_b_idxS = unModSrmMassSpec.findClosest(peakB.getMass());
        unMod_y_idxS = unModSrmMassSpec.findClosest(peakY.getMass());

        /*if (idx == 13786) {
         DEBUG_MSG("vert " << i << " bmass=" << peakB.getMass() << " prmmass=" << prmMassSpec[b_idx][0] << "unmod prmmass=" << unModPrmMassSpec[unMod_b_idx][0]);
         }*/

        if (peakB == prmMassSpec[b_idx][0])
        {
          BpeakToRes[i] = b_idx;
        }
        else if (peakY == srmMassSpec[y_idxS][0])
        {
          BpeakToRes[i] = numMasses - 1 - y_idxS;
        }
        else if (peakB == unModPrmMassSpec[unMod_b_idx][0])
        {
          BpeakToRes[i] = unMod_b_idx;
        }
        else if (peakY == unModSrmMassSpec[unMod_y_idxS][0])
        {
          BpeakToRes[i] = numMasses - 1 - unMod_y_idxS;
        }

        if (peakB == srmMassSpec[b_idxS][0])
        {
          YpeakToRes[i] = numMasses - 1 - b_idxS;
        }
        else if (peakY == prmMassSpec[y_idx][0])
        {
          YpeakToRes[i] = y_idx;
        }
        else if (peakB == unModSrmMassSpec[unMod_b_idxS][0])
        {
          YpeakToRes[i] = numMasses - 1 - unMod_b_idxS;
        }
        else if (peakY == unModPrmMassSpec[unMod_y_idx][0])
        {
          YpeakToRes[i] = unMod_y_idx;
        }
      }
    }

    float bEndpt1 = 0;
    float yEndpt1 = parentMass - AAJumps::massHion;

    float bEndpt2 = parentMass - AAJumps::massMH;
    float yEndpt2 = AAJumps::massH2O;

    float mass, intensity;

    peaks.resize(size());
    for (int i = 0; i < peaks.size(); i++)
    {
      peaks[i].reset();

      mass = peakList[i][0];
      intensity = peakList[i][1];
      peaks[i].set(mass, intensity, peakTol);
      peaks[i].peakIdx = i;
      peaks[i].specIdx = idx;
      peaks[i].mapped = false;

      if (BpeakToRes.count(i) > 0)
      {
        peaks[i].mapped = true;
        peaks[i].annotation = "b";
        peaks[i].BpeptideIdx = BpeakToRes[i];
        numPRMs++;
      }
      if (YpeakToRes.count(i) > 0)
      {
        peaks[i].mapped = true;
        peaks[i].annotation = "y";
        peaks[i].YpeptideIdx = YpeakToRes[i];
        numSRMs++;
      }
      if (BpeakToRes.count(i) > 0 && YpeakToRes.count(i) > 0)
      {
        peaks[i].annotation = "b;y";
      }

      if (peaks[i] == bEndpt1 || peaks[i] == bEndpt2 || peaks[i] == yEndpt1
          || peaks[i] == yEndpt2)
      {
        peaks[i].endpt = true;
      }
    }
  }
}
