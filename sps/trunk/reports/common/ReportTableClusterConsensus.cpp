////////////////////////////////////////////////////////////////////////////////
#include "ReportTableClusterConsensus.h"

////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
ReportTableClusterConsensus::ReportTableClusterConsensus(const string &projectPath, const string &tableFilename, int columnFilter)
: ReportTableBase(projectPath, tableFilename)
{
  // set the filter column
  setFilterColumn(columnFilter);
  // define view (cluster list)
  defineView();
  // set sort column
  setIdColumn(0);
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableClusterConsensus::defineView(void)
{
  // clear view
  clearView();

  // auxiliary object for column type holding
  ReportColumnTypeString        *auxS;
  ReportColumnTypeImageOnDemand *auxI;
  ReportColumnTypeBox           *auxB;

// colTypes[0] -> (CTstring) Spectrum index
//   --> text = index in specset = "<0>"
//   --> columnLabel = "Scan"
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Cluster<br>Index";
  auxS->text         = "<0>"; //"«0»";
  m_colTypes.push_back(auxS);

// colTypes[7] -> (CTstring) Model %
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Model";
  auxS->text         = "<14>";
  m_colTypes.push_back(auxS);

  auxI = new ReportColumnTypeImageOnDemand();
  auxI->link              = "cluster.<0>.0.html";
  auxI->columnLabel       = "Spectrum";
  auxI->label             = "";
  auxI->iconRenderer      = "specplot";
  auxI->iconDisplayLevel  = 3;
  auxI->alt               = "link";
  //auxI->iconParams      = "--pklbin <projectdir>/spectra/specs_ms.pklbin --spectrum <0> --peptide <4> --target internal --encoding uu64 --zoom 0.3 --notitle";
  auxI->validator       = "<3|4|5>";

  auxI->files.push_back( ReportParamsFiles(1,   SPS_FILE_CONSENSUS_SPECTRA, "",   ""    ) );

  auxI->iconParams.push_back( ReportParamsOption("--spectrum",        "<0>",                                    ""    ) );
  auxI->iconParams.push_back( ReportParamsOption("--peptide",         "<3|4|5>",                                ""    ) );
  auxI->iconParams.push_back( ReportParamsOption("--target",          "internal",                               ""    ) );
  auxI->iconParams.push_back( ReportParamsOption("--aa",              "<aas_file>",                             "<aas_file>" ) );
  auxI->iconParams.push_back( ReportParamsOption("--encoding",        "uu64",                                   ""    ) );
  auxI->iconParams.push_back( ReportParamsOption("--zoom",            "0.35",                                   ""    ) );
//  auxI->iconParams.push_back( ReportParamsOption("--annotation-model","<15>",                                   "<15>") );
//  auxI->iconParams.push_back( ReportParamsOption("--shift-value",     "<16>",                                   "<16>") );
  auxI->iconParams.push_back( ReportParamsOption("--notitle",         "",                                       ""    ) );

  m_colTypes.push_back(auxI);


  // colTypes[2] -> (CTseqsBox)
  auxB = new ReportColumnTypeBox();
  auxB->columnLabel  = "Peptide";

  auxI = new ReportColumnTypeImageOnDemand();
  auxI->columnLabel       = "Reference";
  auxI->label             = "<3>";
  auxI->renderer          = "specplot";
  auxI->linkDisplayLevel  = 3;
  //auxI->params            = "--pklbin <projectdir>/spectra/specs_ms.pklbin --spectrum <0> --peptide <3> --title \"Consensus Spectrum <0> (contig <1>)\" --target internal --encoding uu64";
  auxI->validator         = "<3>";
  auxI->splitLabel        = true;

  auxI->files.push_back( ReportParamsFiles(1,   SPS_FILE_CONSENSUS_SPECTRA, "",   ""    ) );

  //auxI->params.push_back( ReportParamsOption("--pklbin",        "<projectdir>/spectra/specs_ms.pklbin",   ""   ) );
  auxI->params.push_back( ReportParamsOption("--title",           "Cluster Consensus Spectrum <0> (contig <1>)",    ""    ) );
  auxI->params.push_back( ReportParamsOption("--spectrum",        "<0>",                                    ""    ) );
  auxI->params.push_back( ReportParamsOption("--peptide",         "<3>",                                    ""    ) );
  auxI->params.push_back( ReportParamsOption("--aa",              "<aas_file>",                             "<aas_file>" ) );
  auxI->params.push_back( ReportParamsOption("--target",          "internal",                               ""    ) );
//  auxI->params.push_back( ReportParamsOption("--annotation-model","<15>",                                   "<15>") );
//  auxI->params.push_back( ReportParamsOption("--shift-value",     "<16>",                                   "<16>") );
  auxI->params.push_back( ReportParamsOption("--encoding",        "uu64",                                   ""    ) );

  auxB->sequences.push_back(auxI);


  auxI = new ReportColumnTypeImageOnDemand();
  auxI->columnLabel       = "Homolog";
  auxI->label             = "<4>";
  auxI->renderer          = "specplot";
  auxI->linkDisplayLevel  = 3;
  //auxI->params            = "--pklbin <projectdir>/spectra/specs_ms.pklbin --spectrum <0> --peptide <4> --title \"Consensus Spectrum <0> (contig <1>)\" --target internal --encoding uu64";
  auxI->validator         = "<4>";
  auxI->splitLabel        = true;

  auxI->files.push_back( ReportParamsFiles(1,   SPS_FILE_CONSENSUS_SPECTRA, "",   ""    ) );

//  auxI->params.push_back( ReportParamsOption("--pklbin",        "<projectdir>/spectra/specs_ms.pklbin",     ""    ) );
  auxI->params.push_back( ReportParamsOption("--title",           "Cluster Consensus Spectrum <0> (contig <1>)",      ""    ) );
  auxI->params.push_back( ReportParamsOption("--spectrum",        "<0>",                                      ""    ) );
  auxI->params.push_back( ReportParamsOption("--peptide",         "<4>",                                      ""    ) );
  auxI->params.push_back( ReportParamsOption("--aa",              "<aas_file>",                               "<aas_file>" ) );
  auxI->params.push_back( ReportParamsOption("--target",          "internal",                                 ""    ) );
//  auxI->params.push_back( ReportParamsOption("--annotation-model","<15>",                                     "<15>") );
//  auxI->params.push_back( ReportParamsOption("--shift-value",     "<16>",                                     "<16>") );
  auxI->params.push_back( ReportParamsOption("--encoding",        "uu64",                                     ""    ) );

  auxB->sequences.push_back(auxI);

  auxI = new ReportColumnTypeImageOnDemand();
  auxI->columnLabel       = "de Novo";
  auxI->label             = "<5>";
  auxI->renderer          = "specplot";
  auxI->linkDisplayLevel  = 3;
  //auxI->params            = "--pklbin <projectdir>/spectra/specs_ms.pklbin --spectrum <0> --peptide <5> --title \"Consensus Spectrum <0> (contig <1>)\" --target internal --encoding uu64";
  auxI->validator         = "<5>";
  auxI->splitLabel        = true;

  auxI->files.push_back( ReportParamsFiles(1,   SPS_FILE_CONSENSUS_SPECTRA, "",   ""    ) );

//  auxI->params.push_back( ReportParamsOption("--pklbin",        "<projectdir>/spectra/specs_ms.pklbin",   ""    ) );
  auxI->params.push_back( ReportParamsOption("--title",           "Cluster Consensus Spectrum <0> (contig <1>)",    ""    ) );
  auxI->params.push_back( ReportParamsOption("--spectrum",        "<0>",                                    ""    ) );
  auxI->params.push_back( ReportParamsOption("--peptide",         "<5>",                                    ""    ) );
  auxI->params.push_back( ReportParamsOption("--aa",              "<aas_file>",                             "<aas_file>" ) );
  auxI->params.push_back( ReportParamsOption("--target",          "internal",                               ""    ) );
//  auxI->params.push_back( ReportParamsOption("--annotation-model","<15>",                                   "<15>") );
//  auxI->params.push_back( ReportParamsOption("--shift-value",     "<16>",                                   "<16>") );
  auxI->params.push_back( ReportParamsOption("--encoding",        "uu64",                                   ""    ) );

  auxB->sequences.push_back(auxI);


  auxI = new ReportColumnTypeImageOnDemand();
  auxI->dynamic           = true;
  auxI->columnLabel       = "User";
  auxI->label             = "<6>";
  auxI->renderer          = "specplot";
  auxI->linkDisplayLevel  = 3;
  //auxI->params            = "--pklbin <projectdir>/spectra/specs_ms.pklbin --spectrum <0> --peptide <6> --title \"Consensus Spectrum <0> (contig <1>)\" --target internal --encoding uu64";
  auxI->validator         = "<6>";
  auxI->splitLabel        = true;

  auxI->files.push_back( ReportParamsFiles(1,   SPS_FILE_CONSENSUS_SPECTRA, "",   ""    ) );

//  auxI->params.push_back( ReportParamsOption("--pklbin",        "<projectdir>/spectra/specs_ms.pklbin",   ""    ) );
  auxI->params.push_back( ReportParamsOption("--title",           "Cluster Consensus Spectrum <0> (contig <1>)",    ""    ) );
  auxI->params.push_back( ReportParamsOption("--spectrum",        "<0>",                                    ""    ) );
  auxI->params.push_back( ReportParamsOption("--peptide",         "<6>",                                    ""    ) );
  auxI->params.push_back( ReportParamsOption("--aa",              "<aas_file>",                             "<aas_file>" ) );
  auxI->params.push_back( ReportParamsOption("--target",          "internal",                               ""    ) );
//  auxI->params.push_back( ReportParamsOption("--annotation-model","<15>",                                   "<15>") );
//  auxI->params.push_back( ReportParamsOption("--shift-value",     "<16>",                                   "<16>") );
  auxI->params.push_back( ReportParamsOption("--encoding",        "uu64",                                   ""    ) );

  auxB->sequences.push_back(auxI);


  auxS = new ReportColumnTypeString();
  auxS->isInput      = true;
  auxS->dynamic      = true;
  auxS->id           = "input_<row>_<col>";
  auxB->sequences.push_back(auxS);

  auxS = new ReportColumnTypeString();
  auxS->isButton     = true;
  auxS->dynamic      = true;
  auxS->text         = "Update";
  auxS->link         = "--update <0> --pklbin <1> --peptide <6>";
  auxS->onClick      = "javascript:DoOnCick(\"input_<row>_<col>\", n, 6);";
  auxB->sequences.push_back(auxS);

  m_colTypes.push_back(auxB);


// colTypes[2] -> (CTstring) mass
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Mass (m)";
  auxS->text         = "<7>";
  m_colTypes.push_back(auxS);

// colTypes[3] -> (CTstring) charge
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Charge (z)";
  auxS->text         = "<8>";
  m_colTypes.push_back(auxS);

// colTypes[4] -> (CTstring) B%
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "B (%)";
  auxS->text         = "<9>";
  m_colTypes.push_back(auxS);

// colTypes[5] -> (CTstring) Y%
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Y (%)";
  auxS->text         = "<10>";
  m_colTypes.push_back(auxS);

// colTypes[6] -> (CTstring) BY intensity %
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "BY Intensity (%)";
  auxS->text         = "<11>";
  m_colTypes.push_back(auxS);

  // colTypes[2] -> (CTstring) mass
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Tool";
  auxS->text         = "<13>";
  m_colTypes.push_back(auxS);
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableClusterConsensus::defineView2(void)
{
  // clear view
  clearView();

  // auxiliary object for column type holding
  ReportColumnTypeString        *auxS;
  ReportColumnTypeImageOnDemand *auxI;
  ReportColumnTypeBox           *auxB;


  auxB = new ReportColumnTypeBox();
  auxB->columnLabel  = "Sequences";

  auxI = new ReportColumnTypeImageOnDemand();
  auxI->columnLabel       = "";
  auxI->label             = "";
  auxI->iconRenderer      = "specplot";
  auxI->iconDisplayLevel  = 3;
  //auxI->iconParams      = "--pklbin <projectdir>/spectra/specs_ms.pklbin --spectrum <0> --peptide <4> --title \"Consensus Spectrum <0> (contig <1>)\" --target internal --encoding uu64";

  auxI->files.push_back( ReportParamsFiles(1,   SPS_FILE_CONSENSUS_SPECTRA, "",   ""    ) );

//  auxI->iconParams.push_back( ReportParamsOption("--pklbin",        "<projectdir>/spectra/specs_ms.pklbin",   ""    ) );
  auxI->iconParams.push_back( ReportParamsOption("--title",           "Cluster Consensus Spectrum <0> (contig <1>)",    ""    ) );
  auxI->iconParams.push_back( ReportParamsOption("--spectrum",        "<0>",                                    ""    ) );
  auxI->iconParams.push_back( ReportParamsOption("--peptide",         "<3|4|5>",                                ""    ) );
  auxI->iconParams.push_back( ReportParamsOption("--aa",              "<aas_file>",                             "<aas_file>" ) );
  auxI->iconParams.push_back( ReportParamsOption("--target",          "internal",                               ""    ) );
//  auxI->iconParams.push_back( ReportParamsOption("--annotation-model","<15>",                                   "<15>") );
//  auxI->iconParams.push_back( ReportParamsOption("--shift-value",     "<16>",                                   "<16>") );
  auxI->iconParams.push_back( ReportParamsOption("--encoding",        "uu64",                                   ""    ) );

  auxB->sequences.push_back(auxI);


  auxS = new ReportColumnTypeString();
  auxS->isInput      = true;
  auxS->dynamic      = true;
  auxS->id           = "input_<row>_<col>";
  auxB->sequences.push_back(auxS);

  auxS = new ReportColumnTypeString();
  auxS->isButton     = true;
  auxS->dynamic      = true;
  auxS->text         = "Update";
  auxS->link         = "--update <0> --pklbin <1> --peptide <6>";
  auxS->onClick      = "javascript:DoOnCick(\"input_<row>_<col>\", n, 6);";
  auxB->sequences.push_back(auxS);

  m_colTypes.push_back(auxB);
}
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
void ReportTableClusterConsensus::defineViewImages(void)
{
  // clear view
  clearView();

  // auxiliary object for column type holding
  ReportColumnTypeString        *auxS;
  ReportColumnTypeImageOnDemand *auxI;
  ReportColumnTypeBox           *auxB;


  // the image
  auxI = new ReportColumnTypeImageOnDemand();
  auxI->id                = "ls<0>";
  auxI->iconRenderer      = "specplot";
  auxI->iconDisplayLevel  = 3;
  auxI->validator       = "<3|4|5>";
  auxI->files.push_back( ReportParamsFiles(1,   SPS_FILE_CONSENSUS_SPECTRA, "",   ""    ) );
  auxI->iconParams.push_back( ReportParamsOption("--spectrum",        "<0>",                                    ""    ) );
  auxI->iconParams.push_back( ReportParamsOption("--peptide",         "<3|4|5>",                                ""    ) );
  auxI->iconParams.push_back( ReportParamsOption("--aa",              "<aas_file>",                             "<aas_file>" ) );
  auxI->iconParams.push_back( ReportParamsOption("--target",          "internal",                               ""    ) );
  auxI->iconParams.push_back( ReportParamsOption("--encoding",        "uu64",                                   ""    ) );
//  auxI->iconParams.push_back( ReportParamsOption("--annotation-model","<15>",                                   "<15>") );
//  auxI->iconParams.push_back( ReportParamsOption("--shift-value",     "<16>",                                   "<16>") );
  auxI->iconParams.push_back( ReportParamsOption("--zoom",            "0.35",                                   ""    ) );
  auxI->iconParams.push_back( ReportParamsOption("--notitle",         "",                                       ""    ) );
  m_colTypes.push_back(auxI);

  auxI = new ReportColumnTypeImageOnDemand();
  auxI->id                = "lR<0>";
  auxI->renderer          = "specplot";
  auxI->linkDisplayLevel  = 3;
  auxI->validator         = "<3>";
  auxI->files.push_back( ReportParamsFiles(1,   SPS_FILE_CONSENSUS_SPECTRA, "",   ""    ) );
  auxI->params.push_back( ReportParamsOption("--title",           "Cluster Consensus Spectrum <0> (contig <1>)",    ""    ) );
  auxI->params.push_back( ReportParamsOption("--spectrum",        "<0>",                                    ""    ) );
  auxI->params.push_back( ReportParamsOption("--peptide",         "<3>",                                    ""    ) );
  auxI->params.push_back( ReportParamsOption("--aa",              "<aas_file>",                             "<aas_file>" ) );
  auxI->params.push_back( ReportParamsOption("--target",          "internal",                               ""    ) );
//  auxI->params.push_back( ReportParamsOption("--annotation-model","<15>",                                   "<15>") );
//  auxI->params.push_back( ReportParamsOption("--shift-value",     "<16>",                                   "<16>") );
  auxI->params.push_back( ReportParamsOption("--encoding",        "uu64",                                   ""    ) );
  m_colTypes.push_back(auxI);

  auxI = new ReportColumnTypeImageOnDemand();
  auxI->id                = "lH<0>";
  auxI->renderer          = "specplot";
  auxI->linkDisplayLevel  = 3;
  auxI->validator         = "<4>";
  auxI->files.push_back( ReportParamsFiles(1,   SPS_FILE_CONSENSUS_SPECTRA, "",   ""    ) );
  auxI->params.push_back( ReportParamsOption("--title",           "Cluster Consensus Spectrum <0> (contig <1>)",    ""    ) );
  auxI->params.push_back( ReportParamsOption("--spectrum",        "<0>",                                    ""    ) );
  auxI->params.push_back( ReportParamsOption("--peptide",         "<4>",                                    ""    ) );
  auxI->params.push_back( ReportParamsOption("--aa",              "<aas_file>",                             "<aas_file>" ) );
  auxI->params.push_back( ReportParamsOption("--target",          "internal",                               ""    ) );
//  auxI->params.push_back( ReportParamsOption("--annotation-model","<15>",                                   "<15>") );
//  auxI->params.push_back( ReportParamsOption("--shift-value",     "<16>",                                   "<16>") );
  auxI->params.push_back( ReportParamsOption("--encoding",        "uu64",                                   ""    ) );
  m_colTypes.push_back(auxI);

  auxI = new ReportColumnTypeImageOnDemand();
  auxI->id                = "lN<0>";
  auxI->renderer          = "specplot";
  auxI->linkDisplayLevel  = 3;
  auxI->validator         = "<5>";
  auxI->files.push_back( ReportParamsFiles(1,   SPS_FILE_CONSENSUS_SPECTRA, "",   ""    ) );
  auxI->params.push_back( ReportParamsOption("--title",           "Cluster Consensus Spectrum <0> (contig <1>)",    ""    ) );
  auxI->params.push_back( ReportParamsOption("--spectrum",        "<0>",                                    ""    ) );
  auxI->params.push_back( ReportParamsOption("--peptide",         "<5>",                                    ""    ) );
  auxI->params.push_back( ReportParamsOption("--aa",              "<aas_file>",                             "<aas_file>" ) );
  auxI->params.push_back( ReportParamsOption("--target",          "internal",                               ""    ) );
//  auxI->params.push_back( ReportParamsOption("--annotation-model","<15>",                                   "<15>") );
//  auxI->params.push_back( ReportParamsOption("--shift-value",     "<16>",                                   "<16>") );
  auxI->params.push_back( ReportParamsOption("--encoding",        "uu64",                                   ""    ) );
  m_colTypes.push_back(auxI);

  auxI = new ReportColumnTypeImageOnDemand();
  auxI->id                = "l<0>";
  auxI->iconRenderer      = "specplot";
  auxI->iconDisplayLevel  = 3;
  auxI->files.push_back( ReportParamsFiles(1,   SPS_FILE_CONSENSUS_SPECTRA, "",   ""    ) );
  auxI->iconParams.push_back( ReportParamsOption("--title",           "Cluster Consensus Spectrum <0> (contig <1>)",    ""    ) );
  auxI->iconParams.push_back( ReportParamsOption("--spectrum",        "<0>",                                    ""    ) );
  auxI->iconParams.push_back( ReportParamsOption("--peptide",         "<3|4|5>",                                ""    ) );
  auxI->iconParams.push_back( ReportParamsOption("--aa",              "<aas_file>",                             "<aas_file>" ) );
//  auxI->iconParams.push_back( ReportParamsOption("--annotation-model","<15>",                                   "<15>") );
//  auxI->iconParams.push_back( ReportParamsOption("--shift-value",     "<16>",                                   "<16>") );
  auxI->iconParams.push_back( ReportParamsOption("--target",          "internal",                               ""    ) );
  auxI->iconParams.push_back( ReportParamsOption("--encoding",        "uu64",                                   ""    ) );
  m_colTypes.push_back(auxI);
}
////////////////////////////////////////////////////////////////////////////////
}; // namespace
////////////////////////////////////////////////////////////////////////////////
