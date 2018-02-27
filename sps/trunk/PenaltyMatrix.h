#ifndef __PenaltyMatrix_H__
#define __PenaltyMatrix_H__

#include "aminoacid.h"

// System Includes
#include <map>
#include <set>
#include <string>
#include <vector>


namespace specnets
{
  const unsigned int CLEAVAGE_START = 0;
  const unsigned int CLEAVAGE_INTERNAL = 1;
  const unsigned int CLEAVAGE_END = 2;

  /*! \brief Holds a penalty matrix.

   A penalty matrix is a two-dimensional matrix indexed by Amino Acid Sequence
   and mass <string, float>.
  */
  class PenaltyMatrix
  {
  public:
    //! \name CONSTRUCTORS
    //@{
    /*! \brief The default constructor.
    */
    PenaltyMatrix(AAJumps & jumps, 
                  float resolution = 1.0, 
                  float knownModPenalty = 1.0, 
                  float unknownPenalty = 1.0, 
                  float unknownMultiplier = 2.0,
                  float minModMass = 100.0,
                  float maxModMass = 100.0);

    //! \name DESTRUCTOR
    //@{
    virtual ~PenaltyMatrix(void);
    //@}

    //! \name ACCESSORS
    //@{
    char ntermChar(void);
    string ntermString(void);
    char ctermChar(void);
    string ctermString(void);
    bool existsNterm(string & seq);
    bool existsCterm(string & seq);
    bool isInMatrix(string & seq, float mass);
    bool isKnown(string & seq, float mass);
    bool isNterm(float mass);
    bool isNterm(string & seq, float mass);
    bool isCterm(float mass);
    bool isCterm(string & seq, float mass);
    float operator()(string & seq, float mass, float averagePeakIntensity = 0.0);
    float operator()(char & aa, float mass, float averagePeakIntensity = 0.0);

    void  getPenalties(string & seq, std::map<float, float> & penaltyMap, float averagePeakIntensity = 0.0);
    float getMass(string aa);
    float getUnknownPenalty(float averagePeakIntensity = 0.0);
    float getKnownPenalty(float averagePeakIntensity = 0.0);
    const map<string, set<float> > & getKnownMods(void);
    const set<float> & getNtermMods(void);
    const set<float> & getNtermMods(string & strAA);
    const map<string, set<float> > & getAllNtermMods(void);
    const set<float> & getCtermMods(void);
    const set<float> & getCtermMods(string & strAA);
    const map<string, set<float> > & getAllCtermMods(void);

    //  Location 0 = before first AA, 1 = internal AA, 2 = end AA
    float getCleavagePenalty(char c, int location);

    bool saveMatrix(std::string & filename);
    bool saveCleavagePenalties(std::string & filename);
    bool saveKnownMods(string & filename);
    bool saveAminoAcids(string & filename);

    void getAminoAcids(vector<string> & aaVec);

    void getSpecProbMods(vector<pair<unsigned int, bool> > & ntermMods,
                         vector<unsigned int> & mods);
    //@}

    //! \name MODIFIERS
    //@{
    bool load(std::string & modFileName,
              std::string & knowmModsFileName,
              std::string & cleavagePenaltiesFileName);
              
    bool loadAminoAcids(string & filename);

    bool loadFromBlosum(std::string & filename, float peakEquivalents);

    bool createFromModificationFreqs(map<float, float> & modFreq,
                                     float minPeakEquivalents,
                                     float maxPeakEquivalents,
                                     float minFrequency,
                                     float averagePeakIntensity = 1.0);

    bool loadKnownModifications(std::string & filename);

    bool loadCleavagePenalties(std::string & filename);

    void addKnownModification(string strAA, float mass);
    void addKnownNtermModification(string strAA, float mass);
    void addKnownCtermModification(string strAA, float mass);
    void addCleavagePenalty(int location, char c, float penalty);
    //@}

    void debug(void);


  protected:

    std::map<std::string, std::map<float, float> > penalties;
    std::map<std::string, float> mapCharMods;
    std::map<float, std::string> mapModChars;
    float m_resolution;
    float m_knownModPenalty;
    float m_unknownPenalty;
    float m_unknownMultiplier;
    float m_allSpectraAveragePeakIntensity;
    float m_minModMass;
    float m_maxModMass;
    
  private:
    //! \name NOT IMPLEMENTED
    //@{
    PenaltyMatrix(const PenaltyMatrix & that);
    //@}

    float roundMass(float mass);

    map<string, set<float> > m_knownMods;
    map<string, set<float> > m_knownNtermMods;
    map<string, set<float> > m_knownCtermMods;
    vector<map<char, float> > m_cleavagePenalties;
  };

} // namespace specnets

#endif // __PenaltyMatrix_H__

