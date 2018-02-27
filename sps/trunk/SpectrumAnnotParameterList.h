/**
 * Parsing for spectrum annotation statistics
 */

#ifndef SPECTRUM_ANNOT_STATISTICS_PARSE_H_
#define SPECTRUM_ANNOT_STATISTICS_PARSE_H_

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

namespace specnets
{
  class SpectrumAnnotParameter
  {
  public:
    /** Spectrum annotation parameter  parsing
     *	@param header - string vector of column names for input
     *	@param fields - string vector of current line split by "\t"
     */
    bool parseFromFields(const vector<string>& header,
                         const vector<string>& fields);

    /**
     * Gets the ion names assigned to a specific fragmentation type
     * @param fragType
     * @return ion names
     */
    string getFragSpecificIonNames(const string& fragType) const;

    string ionNames;
    string statistic;
    string statisticName;
  };

  /*! \brief Class to load in statistics configuration files.

   Statistics configuration files start with a line indiciating the file version
   followed by a header "ion types statistic" and each line following lists the ion types
   and the statistic being considered. As an example:

   Spectrum Statistics V1
   #ion type names should match those in the ms2model file used in SpecNets, except for "all" which
   #accepts all possible types and "na", which indicates that it's not applicable.
   #allowed statistic types: %explained intensity, %explained peaks, %observed ions, total peaks
   # parent mass error ppm, parent mass error da, %observed breaks
   ion types statistic
   all %explained intensity
   all  %explained peaks
   b %observed ions
   y %observed ions
   b-iso,b-NH3,b-H2O,b-H2O-H2O,b-H2O-NH3,a,a-H2O,a-NH3,y-iso,y-NH3,y-H2O,y-H2O-NH3,y-H2O-H2O %observed ions
   b++,b++-NH3,b++-H2O,b++-H2O-H20,b++-H2O-NH3,y++,y++-NH3,y++-H2O,y++-H2O-NH3,y++-H2O-H2O,P++-H2O,P++-NH3,P++ %observed ions
   b %explained intensity
   y %explained intensity
   b-iso,b-NH3,b-H2O,b-H2O-H2O,b-H2O-NH3,a,a-H2O,a-NH3,y-iso,y-NH3,y-H2O,y-H2O-NH3,y-H2O-H2O %explained intensity
   b++,b++-NH3,b++-H2O,b++-H2O-H20,b++-H2O-NH3,y++,y++-NH3,y++-H2O,y++-H2O-NH3,y++-H2O-H2O,P++-H2O,P++-NH3,P++ %explained intensity
   na  total peaks
   na  parent mass error ppm
   na  parent mass error da
   y,b,y++,b++ %observed breaks
   */

  class SpectrumAnnotParameterList
  {
  public:

    /*! \brief Load spectrum statistics configuration file.
     *
     * @param spectrum_annot_file Path to configuration file.
     */
    bool loadSpectrumAnnotFile(const char *spectrum_annot_file);

    /*! \brief Return associated m_params SpectrumAnnotParameter for
     that vector position.
     */
    SpectrumAnnotParameter & operator[](unsigned int i);

    /*! \brief Return associated m_params SpectrumAnnotParameter
     for that vector position.
     */
    const SpectrumAnnotParameter & operator[](unsigned int i) const;

    /*! \brief Set value of one ParameterList to another

     */
    SpectrumAnnotParameterList & operator=(SpectrumAnnotParameterList &other);

    /*! \brief Returns size of m_params vector

     */
    unsigned int size();

    /*! \brief Resizes m_params vector

    @param newSize the new size of the parameter vector
     */
    unsigned int resize(unsigned int newSize);

    vector<SpectrumAnnotParameter> m_params;
  };
}
#endif /* SPECTRUM_STATISTICS_H_ */
