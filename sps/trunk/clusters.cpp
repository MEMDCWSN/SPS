#include "clusters.h"
#include "aminoacid.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

namespace specnets
{
	// -------------------------------------------------------------------------

	Clusters::Clusters() {
		consensus.resize(0);
		specIdx.resize(0);
		shifts.resize(0);
		endpoints.resize(0);
	}

	// -------------------------------------------------------------------------

	Clusters::Clusters(SpecSet &spectra) {
			consensus = spectra;
			specIdx.resize(spectra.size());
			shifts.resize(spectra.size());
			endpoints.resize(spectra.size());

			for(unsigned int i=0; i<consensus.size(); i++) {
				specIdx[i].resize(0);
				shifts[i].resize(0);
				endpoints[i].resize(0);
			}
	}

	Clusters& Clusters::operator =(const Clusters &other) {
		if (this == &other) {
			return *this;
		}
		consensus = other.consensus;
		specIdx = other.specIdx;
		shifts = other.shifts;
		endpoints = other.endpoints;

		return (*this);
	}

	// -------------------------------------------------------------------------

	void Clusters::setResolution(float newResolution, bool enforceRounding) {
		unsigned int s,i;
		if (enforceRounding)
			for(s=0; s<consensus.size(); s++) {
				consensus[s].setResolution(newResolution,enforceRounding);
				for(i=0; i<shifts[s].size(); i++) shifts[s][i][0]=round(shifts[s][i][0]/newResolution);
			}
		else
			for(s=0; s<consensus.size(); s++) {
				consensus[s].setResolution(newResolution,enforceRounding);
				for(i=0; i<shifts[s].size(); i++) shifts[s][i][0]=shifts[s][i][0]/newResolution;
			}
	}

	/*
	 * /**
     * Indices of the spectra in each cluster (there should be some companion
     * .pkl spectra file).

    vector<vector<int> > specIdx;

    /**
     * Shifts (left and right) of each spectrum in each cluster.

    vector<vector<TwoValues<float> > > shifts;

    /**
     * Consensus spectrum information per cluster - Consensus PRMs.

    SpecSet consensus;

    /**
     * Endpoint types: 0=internal, 1=b0/bN, 2=y0/yN

    vector<vector<short> > endpoints;
    *
    */

	void Clusters::initializeFromAbinfo(SpecSet& contigSpectra,
                                      abinfo_t& componentSets,
                                      SpecSet& starSpectra,
                                      float pkTol,
                                      float pmTol)
  {

    consensus.resize(contigSpectra.size());
    specIdx.resize(consensus.size());
    shifts.resize(consensus.size());
    endpoints.resize(consensus.size());

    starSpectra.addZPMpeaks(pkTol, 0, true);
    vector<bool> stars_reversed(starSpectra.size(), false);

    //map<int, int> specIdxLoc;
    vector<list<TwoValues<int> > > abVertices;

    for (int i = 0; i < contigSpectra.size(); i++) {
      //specIdxLoc.clear();
      if (contigSpectra[i].size() == 0 || componentSets.count(i) == 0
                || componentSets[i].second.size() == 0) {
        consensus[i].resize(0);
        endpoints[i].resize(0);
        shifts[i].resize(0);
        specIdx[i].resize(0);
        continue;
      }

      consensus[i] = contigSpectra[i];
      specIdx[i] = componentSets[i].first.first;
      shifts[i].resize(specIdx[i].size(), TwoValues<float>(0, 0));

      for (int j = 0; j < componentSets[i].first.first.size(); j++) {
        int specIdx = componentSets[i].first.first[j];
        //specIdxLoc[specIdx] = i;
        if (componentSets[i].first.second[j] == 1) {
          if (stars_reversed[specIdx]) {
            ERROR_MSG("Having to reverse spectrum " << specIdx << " twice!")
          }
          starSpectra[specIdx].reverse(0);
          stars_reversed[specIdx] = !stars_reversed[specIdx];
        }
      }

      abVertices.resize(componentSets[i].second.size());

      for (int j = 0; j < componentSets[i].second.size(); j++) {
        pair<vector<int> , vector<double> >* vert =
            &componentSets[i].second[j];
        abVertices[j].clear();
        for (int k = 0; k < vert->first.size(); k++) {
          int specIdx = (*vert).first[k];
          float mass = (*vert).second[k];

          /*if (specIdx == 16299) {
           DEBUG_MSG("Contig " << seqIt->first << " peak=" << j << " spectrum " << specIdx << " w/ mass " << mass);
           }*/
          int peakIdx = starSpectra[specIdx].findClosest(mass);

          if (!MZRange::EqualWithinRange(mass,
                                         starSpectra[specIdx][peakIdx][0],
                                         0.01)) {

            // If we couldn't find the mass, try reversing the spectrum
            starSpectra[specIdx].reverse(0);
            stars_reversed[specIdx] = !stars_reversed[specIdx];
            peakIdx = starSpectra[specIdx].findClosest(mass);

            if (!MZRange::EqualWithinRange(mass,
                                           starSpectra[specIdx][peakIdx][0],
                                           0.01)) {
              ERROR_MSG("Expected mass " << mass << " not equal to closest mass " << starSpectra[specIdx][peakIdx][0] << " in spectrum " << specIdx << " from contig " << i);
            }
          }

          abVertices[j].push_back(TwoValues<int> (specIdx, peakIdx));
        }
      }

      set_endpoints(i, starSpectra, abVertices, pkTol, pmTol);
    }

  }

	void Clusters::set_endpoints(int idx, SpecSet &specs, vector<list<TwoValues<int> > > &abVertices,
			float peakTol, float pmTol) {
		float peakMass, scoreBep, scoreYep;
		list<TwoValues<int> >::iterator vertexIter;
		unsigned int setIdx, specIdx, peakIdx, endpointIdx;

		if(consensus[idx].size()!=abVertices.size()) {
			cerr<<"[Warning in Clusters::set_endpoints()] Cannot set endpoint masses because non-matching number of consensus peaks and ABruijn vertices.";
			endpoints[idx].resize(0);
			return;
		}

		endpoints[idx].resize(consensus[idx].size());   endpointIdx = 0;
		for(setIdx=0; setIdx<abVertices.size(); setIdx++) {
			scoreBep = 0;   scoreYep = 0;
			for(vertexIter=abVertices[setIdx].begin(); vertexIter!=abVertices[setIdx].end(); vertexIter++) {
				specIdx = (*vertexIter)[0];   peakIdx = (*vertexIter)[1];   peakMass = specs[specIdx][peakIdx][0];
				if(abs(peakMass)<=peakTol or abs(peakMass+AAJumps::massMH-specs[specIdx].parentMass)<2*pmTol) scoreBep+=specs[specIdx][peakIdx][1];
				if(abs(peakMass-AAJumps::massH2O)<=peakTol or abs(peakMass+AAJumps::massHion-specs[specIdx].parentMass)<2*pmTol) scoreYep+=specs[specIdx][peakIdx][1];
			}
			if(scoreBep+scoreYep<0.00001) endpoints[idx][setIdx]=0;
			else { if(scoreBep>=scoreYep) endpoints[idx][setIdx]=1; else endpoints[idx][setIdx]=2; }
		}
	}

	void Clusters::getSpecIfB(int idx, Spectrum &putHere, float peakTol, float ionOffset) {
		float deltaY0=0;  // Mass added to all peaks when the first peak is y0 (shift right by AAJumps::massH2O)

		// Adjust masses for peaks that are endpoint masses in the assembled spectra
		putHere = consensus[idx];
		if(not endpoints[idx].empty())
			for(unsigned int peakIdx=0; peakIdx<putHere.size(); peakIdx++) {
				putHere[peakIdx][0]+=deltaY0;
				if(endpoints[idx][peakIdx]==2) {
					if(peakIdx==0) {
						deltaY0 = AAJumps::massH2O;
						putHere.parentMass+=AAJumps::massH2O;  // First mass was y0, need to adjust parent mass (for addZPMpeaks)
					} else putHere[peakIdx][0]-=AAJumps::massH2O;
					if(peakIdx==putHere.size()-1) putHere.parentMass-=AAJumps::massH2O;  // Last mass was yN, need to adjust parent mass (for addZPMpeaks)
				}
			}

		putHere.addZPMpeaks(peakTol,ionOffset,false); // Guarantee peaks at b0/bN with default intensity
	}

	void Clusters::getSpecIfY(int idx, Spectrum &putHere, float peakTol, float ionOffset) {
		Spectrum tmp = consensus[idx];
		float aaMass=tmp.parentMass-AAJumps::massMH;

		// Adjust masses for peaks that are endpoint masses in the assembled spectra
		if(not endpoints[idx].empty())
			for(unsigned int peakIdx=0; peakIdx<tmp.size(); peakIdx++) {
				if(endpoints[idx][peakIdx]==1) {
					tmp[peakIdx][0]+=AAJumps::massH2O;
					if(peakIdx==0) tmp.parentMass-=AAJumps::massH2O;  // First mass was b0, need to adjust parent mass (for addZPMpeaks)
					if(peakIdx==putHere.size()-1) {
						tmp.parentMass+=AAJumps::massH2O;  // Last mass was bN, need to adjust parent mass (for addZPMpeaks)
						aaMass+=AAJumps::massH2O;  // Last mass was bN, need to adjust symmetry
					}
				}
			}

		// Reverse merged spectrum
		int topIdx=tmp.size()-1;  while(topIdx>=0 && tmp[topIdx][0]>aaMass) topIdx--;
		putHere.copyNP(tmp);
		putHere.resize(topIdx+1);    putHere.parentMass = tmp.parentMass;
		for(int i=0; i<=topIdx; i++) putHere[topIdx-i].set(aaMass-tmp[i][0],tmp[i][1]);

		putHere.addZPMpeaks(peakTol,ionOffset,false); // Guarantee peaks at b0/bN
	}


	//
	//  File format is:
	//  <numClusters>
	//  Repeat per cluster:
	//    <numEntries>
	//    Repeat per spectrum
	//        <specIdx> <shift1> <shift2>
	//    <num consensus PRMs> <mhMass>
	//    Repeat per consensus PRM
	//        <mass> <score>
	//    <numEndpoints>
	//    Repeat per endpoint
	//        <endpoint types> (0=internal, 1=b0/bN, 2=y0/yN)
	//        <mass> <score>   <<--- discontinued old format

	int Clusters::Load(const char *contigFilename,
	                   const char *spectrumPklbinFilename /* = 0x0 */,
	                   const char *psmFilename /* = 0x0 */)
  {
		ifstream input(contigFilename, ios::binary);
		if (!input) {
       ERROR_MSG("Can not open  [" << contigFilename << "]" );
		  specIdx.resize(0);
		  return -1;
		}

		unsigned int count, cluster, i;

		input >> count;  specIdx.resize(count);  shifts.resize(count);  consensus.resize(count);  endpoints.resize(count);

		DEBUG_VAR(contigFilename);
		DEBUG_VAR(count);
		for(cluster=0; cluster<specIdx.size(); cluster++) {
				// Read spectrum indices and shifts
				input >> count;   specIdx[cluster].resize(count);   shifts[cluster].resize(count);
				for(i=0; i<count; i++) input >> specIdx[cluster][i] >> shifts[cluster][i][0] >> shifts[cluster][i][1];

				// Read consensus spectrum PRMs
				input >> count >> consensus[cluster].parentMass;   consensus[cluster].resize(count);
				for(i=0; i<count; i++) input >> consensus[cluster][i][0] >> consensus[cluster][i][1];

				// Read endpoints
				input >> count;   endpoints[cluster].resize(count);
//        for(i=0; i<count; i++) input >> endpoints[cluster][i][0] >> endpoints[cluster][i][1];
				for(i=0; i<count; i++) input >> endpoints[cluster][i];
		}
		input.close();

    if (spectrumPklbinFilename) {
      if (Load_pklbin(spectrumPklbinFilename, psmFilename) <= 0) {
        return -2;
      }
    }

		return(specIdx.size());
	}

	int Clusters::Load_pklbin(const char *filename,
	                          const char *psmFilename /* = 0x0 */)
       {
			unsigned int count, cluster, i;

			if(not consensus.loadPklBin(filename)) { specIdx.resize(0); shifts.resize(0); endpoints.resize(0); return -1; }

			specIdx.resize(consensus.size());  shifts.resize(consensus.size());  endpoints.resize(consensus.size());
			for(cluster=0; cluster<consensus.size(); cluster++) {
				specIdx[cluster].resize(0);
				shifts[cluster].resize(0);
				endpoints[cluster].resize(0);
			}

     if (psmFilename) {
       PeptideSpectrumMatchSet psmSetTemp;
       if (psmSetTemp.loadFromFile(psmFilename)) {
         DEBUG_VAR(psmSetTemp.size());
         psmSetTemp.addSpectra(&consensus, true);
       } else {
         ERROR_MSG("Problem loading PSM file [" << psmFilename << "]" );
         return -2;
       }
     }


		return(specIdx.size());
	}

	int Clusters::Save(const char *contigFilename,
	                   const char *spectrumPklbinFilename /* = 0x0 */,
	                   const char *psmFilename /* = 0x0 */)
       {
         ofstream output(contigFilename, ios::binary);
         if (!output) { cerr << "ERROR: cannot open " << contigFilename << "\n";   return -1; }

         unsigned int count, cluster, i;

         output << specIdx.size() << endl;
         for(cluster=0; cluster<specIdx.size(); cluster++) {
           output << specIdx[cluster].size() << endl;
           for(i=0; i<specIdx[cluster].size(); i++) output << specIdx[cluster][i] << " " << shifts[cluster][i][0] << " " << shifts[cluster][i][1] << endl;
           output << consensus[cluster].size() << " " << consensus[cluster].parentMass << endl;
           for(i=0; i<consensus[cluster].size(); i++) output << consensus[cluster][i][0] << " " << consensus[cluster][i][1] << endl;
           output << endpoints[cluster].size() << endl;
//         for(i=0; i<endpoints[cluster].size(); i++) output << endpoints[cluster][i][0] << " " << endpoints[cluster][i][1] << endl;
           for(i=0; i<endpoints[cluster].size(); i++) output << endpoints[cluster][i] << endl;
         }
         output.close();

         if (spectrumPklbinFilename) {
           if (consensus.savePklBin(spectrumPklbinFilename, psmFilename) <= 0) {
             ERROR_MSG("Problem saving spectrum file [" << spectrumPklbinFilename << "]" );
             return -2;
           }
         }

         return(0);
	}

}
