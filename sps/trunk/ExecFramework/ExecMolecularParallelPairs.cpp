// Header Include
#include "ExecMolecularParallelPairs.h"
#include "ExecMergeConvert.h"
#include "ExecTagFilterPairs.h"

// Module Includes
#include "AlignmentUtils.h"
#include "Logger.h"
#include "Filters.h"
#include "FileUtils.h"

// External Includes
#include "utils.h"  // for Save_binArray only
// System Includes
#include "stdlib.h"

static bool DEBUG_FILTERPAIRS_SPLIT = false;

using namespace specnets;
using namespace std;

namespace specnets
{

// -------------------------------------------------------------------------
ExecMolecularParallelPairs::ExecMolecularParallelPairs(void) :
        ExecBase()
{
        m_name = "ExecMolecularParallelPairs";
        m_type = "ExecMolecularParallelPairs";

        m_inputSpectraMS2 = new SpecSet;
        m_inputSpectraMS2_normalized = new SpecSet;
}

// -------------------------------------------------------------------------
ExecMolecularParallelPairs::ExecMolecularParallelPairs(const ParameterList & inputParams) :
        ExecBase(inputParams)
{
        m_name = "ExecMolecularParallelPairs";
        m_type = "ExecMolecularParallelPairs";

        m_inputSpectraMS2 = new SpecSet;
        m_inputSpectraMS2_normalized = new SpecSet;

        m_filteredPairs = new SpectrumPairSet();

}



// -------------------------------------------------------------------------
ExecBase * ExecMolecularParallelPairs::clone(const ParameterList & inputParams) const
{
    return new ExecMolecularParallelPairs(inputParams);
}

// -------------------------------------------------------------------------
bool ExecMolecularParallelPairs::invoke(void)
{
    DEBUG_MSG("NORMALIZING PEAKS");
    this->normalizeInputSpectra();

    int startBaseIdx;
    if (m_params.exists("IDX_START"))
    {
            startBaseIdx = max(0, m_params.getValueInt("IDX_START"));
    }
    else
    {
            startBaseIdx = 0;
    }
    int endBaseIdx;
    if (m_params.exists("IDX_END"))
    {
            endBaseIdx = max(0, m_params.getValueInt("IDX_END"));
    }
    else
    {
            endBaseIdx = -1;
    }


    short minNumMatchedPeaks = (short)m_params.getValueInt("MIN_MATCHED_PEAKS", 6);
    float minRatio = m_params.getValueDouble("MIN_RATIO", 0);
    float peakTol = m_params.getValueDouble("TOLERANCE_PEAK", 0.5);
    float pmTol = m_params.getValueDouble("TOLERANCE_PM", 1);
    DEBUG_VAR(minNumMatchedPeaks);
    DEBUG_VAR(minRatio);
    DEBUG_VAR(peakTol);
    DEBUG_VAR(pmTol);

    if (endBaseIdx < 0)
    {
            endBaseIdx = m_inputSpectraMS2->size() - 1;
    }

    DEBUG_VAR(startBaseIdx);
    DEBUG_VAR(endBaseIdx);


    if (startBaseIdx < 0 || startBaseIdx >= m_inputSpectraMS2->size())
    {
            DEBUG_MSG("Invalid start index [" << startBaseIdx << "] (" << m_inputSpectraMS2->size() << " spectra)");
            return false;
    }

    if (endBaseIdx >= m_inputSpectraMS2->size())
    {
            WARN_MSG("IDX_END (" << endBaseIdx << ") was lowered to the index of the last spectrum in the dataset (" << m_inputSpectraMS2->size() - 1 << ")");
    }
    endBaseIdx = min(m_inputSpectraMS2->size() - 1, (unsigned int)endBaseIdx);

    DEBUG_MSG("Computing alignments between spectra " << startBaseIdx << ":" << endBaseIdx << " and all others in the dataset.");

    vector<int> baseSpecIdx(endBaseIdx - startBaseIdx + 1);
    for (int i = 0; i < endBaseIdx - startBaseIdx + 1; i++)
    {
            baseSpecIdx[i] = startBaseIdx + i;
    }

    float minCosine = (float) m_params.getValueDouble("PAIRS_MIN_COSINE", 0.0);
    bool PROJECTED_COSINE = m_params.getValueBool("PAIRS_PROJECTED_COSINE", false);
    double maxShift = m_params.getValueDouble("MAX_SHIFT", 100);
    double minShift = m_params.getValueDouble("MIN_SHIFT", 0);

    float resolution = m_params.getValueDouble("RESOLUTION", 0.1);

    ExecTagFilterPairs tag_filter;
	if( m_params.getValueBool("TAG_FILTER", false) && m_params.exists("INPUT_TAGS") )
	{
		if( m_params.getValueBool("TAG_FILTER", false) && m_params.exists("INPUT_TAGS") ) {

			if( !tag_filter.prepare(m_inputSpectraMS2,
									m_params.getValue("INPUT_TAGS", ""),
									m_params.getValueInt("MAX_TAG_SIZE", 50),
									true) ){
				ERROR_MSG("ExecTagFilterPairs Error!");
				return false;
			}
		}
	}

    vector<TwoValues<float> > ratios;
    vector<TwoValues<float> > means;
    vector<float> varTerms;

    getPairCosines(*m_inputSpectraMS2_normalized,
                   baseSpecIdx,
                   maxShift, // maximum parent mass difference
                   pmTol,
                   peakTol,
                   minCosine,
                   minRatio,
                   minNumMatchedPeaks,
                   *m_filteredPairs,
                   ratios,
                   means,
                   varTerms,
                   &tag_filter,
                   resolution,
                   0,
                   PROJECTED_COSINE );

    return true;
}

// -------------------------------------------------------------------------
bool ExecMolecularParallelPairs::loadInputData(void)
{

    if (m_params.exists("INPUT_SPECTRA_MS2"))
    {
            DEBUG_MSG("Loading MS2 spectra.. " << m_params.getValue("INPUT_SPECTRA_MS2"));
            if (!m_inputSpectraMS2->Load(m_params.getValue("INPUT_SPECTRA_MS2").c_str())) {
                    DEBUG_MSG("Could not load " << m_params.getValue("INPUT_SPECTRA_MS2"));
                    return false;
            }
            DEBUG_VAR(m_inputSpectraMS2->size());
    }

    startBaseIdx = 0;
    endBaseIdx = m_inputSpectraMS2->size() - 1;

    DEBUG_TRACE;
    if (m_inputSpectraMS2->size() == 0) {
            ERROR_MSG("Input spectra size is 0!");
            return false;
    }

    return true;
}

// -------------------------------------------------------------------------
bool ExecMolecularParallelPairs::saveOutputData(void)
{
    string dataDir = "";

    if (m_params.exists("OUTPUT_ALIGNS"))
    {
            string fileName = dataDir + m_params.getValue("OUTPUT_ALIGNS");
            DEBUG_VAR(m_filteredPairs->size());
            DEBUG_MSG("Outputting aligns..." << fileName);
            m_filteredPairs->saveToBinaryFile(fileName);
    }

    return true;
}

// -------------------------------------------------------------------------
bool ExecMolecularParallelPairs::saveInputData(std::vector<std::string> & filenames)
{
    string paramDir = m_params.getValue("GRID_DATA_DIR_PARAMS", m_params.getValue("GRID_DATA_DIR"));
    if (paramDir.empty()) {
            paramDir = ".";
    }
    string baseDirectory = paramDir + "/";
    string baseFilename = baseDirectory + getName();

    string paramFilename = baseFilename + ".params";
    m_params.writeToFile(paramFilename);

    return true;
}

// -------------------------------------------------------------------------
bool ExecMolecularParallelPairs::loadOutputData(void)
{
    string outDir = m_params.getValue("GRID_DATA_DIR_INTERMEDIATE");
    if (outDir.empty()) {
            outDir = ".";
    }

    if (m_params.exists("OUTPUT_ALIGNS"))
    {
            DEBUG_MSG("Loading alignments...");
            string fileName = outDir + "/" + m_params.getValue("OUTPUT_ALIGNS");
            if (m_filteredPairs->loadFromBinaryFile(fileName) < 0)
            {
                    ERROR_MSG("Could not load: " << fileName);
                    return false;
            }
            DEBUG_MSG("Loaded " << m_filteredPairs->size() << " alignments");
    }

    return true;
}

// ------------------------------------------------------------------------

// -------------------------------------------------------------------------

vector<ExecBase*> const & ExecMolecularParallelPairs::split(int numSplit)
{
    DEBUG_VAR(numSplit);

    if (numSplit < 2)
    {
            DEBUG_MSG("Number split [" << numSplit << "] must be at least 2");
            return m_subModules;
    }

    int startBaseIdx;
    if (m_params.exists("IDX_START"))
    {
            startBaseIdx = max(0, m_params.getValueInt("IDX_START"));
    }
    else
    {
            startBaseIdx = 0;
    }
    DEBUG_VAR(startBaseIdx);
    DEBUG_VAR(m_inputSpectraMS2->size());


    DEBUG_TRACE;

    string dataDirIn = m_params.getValue("GRID_DATA_DIR_IN");
    if (dataDirIn.empty())
    {
            dataDirIn = ".";
    }
    string dataDirIntermediate = m_params.getValue("GRID_DATA_DIR_INTERMEDIATE");
    if (dataDirIntermediate.empty())
    {
            dataDirIntermediate = ".";
    }

    int endBaseIdx;
    if (m_params.exists("IDX_END"))
    {
            endBaseIdx = max(0, m_params.getValueInt("IDX_END"));
            endBaseIdx = min(endBaseIdx, (int)m_inputSpectraMS2->size() - 1);
    }
    else
    {
            endBaseIdx = m_inputSpectraMS2->size() - 1;
    }
    DEBUG_VAR(endBaseIdx);
    if (startBaseIdx == endBaseIdx)
    {
            DEBUG_MSG("Must have more than one row");
            return m_subModules;
    }

    int actual_split = numSplit;
    if(m_inputSpectraMS2->size() < numSplit) {
            actual_split = 1;
    }

    DEBUG_VAR(actual_split);

    int total_nonempty_spectra = 0;
    //Lets count the number of spectra with non zero peaks
    for(int i = startBaseIdx; i < endBaseIdx; i++) {
            if((*m_inputSpectraMS2)[i].size() > 5) {
                    total_nonempty_spectra++;
            }
    }
    DEBUG_VAR(total_nonempty_spectra);
    int running_index_start_idx = startBaseIdx;
    for(int nSplit = 0; nSplit < actual_split; nSplit++) {
        int total_spectra_size = endBaseIdx - startBaseIdx + 1;
        int division_size = total_spectra_size/actual_split;
        int startIndex = division_size * nSplit + startBaseIdx;
        int endIndex = division_size * (nSplit + 1) - 1 + startBaseIdx;
        if(endIndex < 0) {
                endIndex = 0;
        }

        if(nSplit == actual_split - 1) {
                endIndex = endBaseIdx;
        }

        int split_nonempty_spectra = 0;
        for(int i = startIndex; i <= endIndex; i++) {
                if((*m_inputSpectraMS2)[i].size() > 5) {
                        split_nonempty_spectra++;
                }
        }

        DEBUG_MSG("Splitting Cosine: " << nSplit <<"\t" << startIndex<<","<<endIndex<<" with "<<split_nonempty_spectra<<" spectra");

        ParameterList childParams(m_params);

        childParams.removeParam("GRID_EXECUTION"); // necessary for Proteosafe
        childParams.setValue("GRID_DATA_DIR_OUT", dataDirIntermediate);

        // But set the start and end indices
        char buf[128];
        sprintf(buf, "%d", startIndex);
        childParams.setValue("IDX_START", buf);
        sprintf(buf, "%d", endIndex);
        childParams.setValue("IDX_END", buf);

        //DEBUG_MSG("Start [" << startIndex << "] End [" << i << "] Split [" << numSplit << "]");
        ExecBase * theClone = new ExecMolecularParallelPairs(childParams);

        theClone->setName(makeName(m_name, nSplit));

        // Have to set up the output files also so the params will be correct on reload
        string baseName;
        if (m_params.exists("RESULTS_DIR")) {
                baseName = m_params.getValue("RESULTS_DIR", "") + "/";
        }
        baseName = baseName + theClone->getName();

        //DEBUG_VAR(baseName);

        theClone->m_params.setValue("OUTPUT_ALIGNS",
                                    baseName + "_aligns" + ".bin");

        std::string suffix("");
        char bufSplit[128];
        sprintf(bufSplit, "%d", nSplit + 1);
        theClone->m_params.setValue("NUM_SPLIT", bufSplit);

        m_subModules.push_back(theClone);
    }

    DEBUG_VAR(m_subModules.size());
    return m_subModules;
}

// -------------------------------------------------------------------------
bool ExecMolecularParallelPairs::merge(void)
{

    unsigned int startIdx = (unsigned int)m_params.getValueInt("IDX_START", 0);
    unsigned int endIdx = (unsigned int)m_params.getValueInt("IDX_END", 0);
    bool alignPA = m_params.getValueBool("PARTIAL_OVERLAPS", false);
    float pmTol = m_params.getValueDouble("TOLERANCE_PM", 1);
    float minRatio = m_params.getValueDouble("MIN_RATIO", 0);
    float minPValue = m_params.getValueDouble("MAX_PVALUE", 0.05);

    // Find out how many total pairs from all the children
    DEBUG_VAR(m_subModules.size());
    int totalPairs = 0;
    for (int i = 0; i < m_subModules.size(); i++)
    {
        ExecMolecularParallelPairs * fpe = (ExecMolecularParallelPairs *)m_subModules[i];
        DEBUG_VAR(fpe);
        DEBUG_VAR(fpe->m_filteredPairs->size());
        totalPairs += fpe->m_filteredPairs->size();
    }

    DEBUG_VAR(totalPairs);
    if (totalPairs <= 0)
    {
        DEBUG_MSG("Can not merge/filter with " << totalPairs << " input pairs!");
        return false;
    }

    // Resize our array to hold all the pairs
    DEBUG_VAR(m_filteredPairs);
    m_filteredPairs->resize(totalPairs);

    DEBUG_VAR(m_filteredPairs->size());

    // Copy all the result pairs from the children into our array
    int pairCount = 0;
    for (int i = 0; i < m_subModules.size(); i++)
    {
        ExecMolecularParallelPairs * fpe = (ExecMolecularParallelPairs *)m_subModules[i];
        SpectrumPairSet * sps = fpe->m_filteredPairs;
        for (int k = 0; k < sps->size(); k++)
        {
                //DEBUG_VAR(k);
                (*m_filteredPairs)[pairCount++] = sps->operator[](k);
        }
    }
    DEBUG_MSG("Merged " << totalPairs << " pairs");

    return true;
}

// -------------------------------------------------------------------------
bool ExecMolecularParallelPairs::validateParams(std::string & error)
{
        return true;
}

void ExecMolecularParallelPairs::normalizeInputSpectra()
{
    //Probably should do this efficiently to only the things we need
    if (m_inputSpectraMS2_normalized == 0x0)
    {
            m_inputSpectraMS2_normalized = new SpecSet;
    }

    m_inputSpectraMS2_normalized->operator =(*m_inputSpectraMS2);

    for (int i = 0; i < m_inputSpectraMS2_normalized->size(); i++)
    {
        for (unsigned int j = 0; j < (*m_inputSpectraMS2_normalized)[i].size();
             j++)
                (*m_inputSpectraMS2_normalized)[i][j][1] =
                        sqrt((*m_inputSpectraMS2_normalized)[i][j][1]);
        (*m_inputSpectraMS2_normalized)[i].normalize2();
    }
}

} // namespace specnets
