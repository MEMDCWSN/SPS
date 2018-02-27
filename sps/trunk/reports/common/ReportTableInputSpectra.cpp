////////////////////////////////////////////////////////////////////////////////
#include "ReportTableInputSpectra.h"

////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
ReportTableInputSpectra::ReportTableInputSpectra(const string &projectPath, const string &tableFilename, int columnFilter)
: ReportTableBase(projectPath, tableFilename)
{
  // set the filter column
  setFilterColumn(columnFilter);
  // define view (cluster list)
  defineView();
  // set the sort column
  setIdColumn(2);
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableInputSpectra::populateIDs(void)
{
  for(int i = 0 ; i < m_cells.size() ; i++)
    m_cells[i][0] = parseInt(i);
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableInputSpectra::defineView(void)
{
  // clear view
  clearView();

  // auxiliary object for column type holding
  ReportColumnTypeString        *auxS;
  ReportColumnTypeImageOnDemand *auxI;
  ReportColumnTypeBox           *auxB;

  // colTypes[0] -> (CTstring) Spectrum index
  auxS = new ReportColumnTypeString();
  auxS->columnLabel   = "Spectrum<br>Index";
  auxS->text          = "<1>";
  m_colTypes.push_back(auxS);

  auxS = new ReportColumnTypeString();
  auxS->columnLabel   = "Spectrum<br>Scan";
  auxS->text          = "<2>";
  m_colTypes.push_back(auxS);

// colTypes[8] -> (CTstring) Model %
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Model";
  auxS->text         = "<20>";
  m_colTypes.push_back(auxS);

  // colTypes[1] -> (CTstring) Spectrum file name
  //auxS = new ReportColumnTypeString();
  //auxS->columnLabel   = "File name";
  //auxS->text          = "<6>";
  //auxS->link          = "cluster.<3>.0.html";
  //m_colTypes.push_back(auxS);

  // colTypes[1] -> (CTseqsBox)
  auxB = new ReportColumnTypeBox();
  auxB->columnLabel   = "Spectrum";
  //auxB->link          = "cluster.<3>.0.html";

  auxI = new ReportColumnTypeImageOnDemand();
  auxI->iconRenderer      = "specplot";
  //auxI->iconParams        = "--pklbin <projectdir>/<6>  --spectrum <1> --peptide <8> --notitle --zoom 0.3 --target internal --encoding uu64";
  auxI->iconDisplayLevel  = 3;
  auxI->alt               = "N/A";
  auxI->validator         = "<7|8|9>";

  auxI->files.push_back( ReportParamsFiles(1,   SPS_FILE_SPECSET, "<5>",   ""    ) );

  //auxI->iconParams.push_back( ReportParamsOption("--pklbin",        "<projectdir>/<6>",                         ""   ) );
  //auxI->iconParams.push_back( ReportParamsOption("--pklbin",        "<6>",                         ""   ) );
  auxI->iconParams.push_back( ReportParamsOption("--spectrum",      "<1>",                                      ""   ) );
  auxI->iconParams.push_back( ReportParamsOption("--peptide",       "<7|8|9>",                                  ""   ) );
  auxI->iconParams.push_back( ReportParamsOption("--aa",            "<aas_file>",                               "<aas_file>" ) );
  auxI->iconParams.push_back( ReportParamsOption("--target",        "internal",                                 ""   ) );
  auxI->iconParams.push_back( ReportParamsOption("--encoding",      "uu64",                                     ""   ) );
  auxI->iconParams.push_back( ReportParamsOption("--zoom",          "0.35",                                     ""   ) );
  auxI->iconParams.push_back( ReportParamsOption("--notitle",       "",                                         ""   ) );

  auxB->sequences.push_back(auxI);


  auxS = new ReportColumnTypeString();
  auxS->text          = "<16>";
  auxS->splitText     = true;
  auxB->sequences.push_back(auxS);

  m_colTypes.push_back(auxB);


  // colTypes[2] -> (CTseqsBox)
  auxB = new ReportColumnTypeBox();
  auxB->columnLabel   = "Sequences";

  auxI = new ReportColumnTypeImageOnDemand();
  auxI->columnLabel  = "Reference";
  auxI->label        = "<7>";
  auxI->renderer     = "specplot";
  auxI->linkDisplayLevel  = 3;
  //auxI->params       = "--pklbin <projectdir>/<6>  --spectrum <1> --peptide <7> --title \"Spectrum Scan <2>\" --target internal --encoding uu64";
  auxI->validator    = "<7>";
  auxI->splitLabel   = true;

  auxI->files.push_back( ReportParamsFiles(1,   SPS_FILE_SPECSET, "<5>",   ""    ) );

  //auxI->params.push_back( ReportParamsOption("--pklbin",        "<projectdir>/<6>",                         ""   ) );
  //auxI->params.push_back( ReportParamsOption("--pklbin",        "<6>",                                      ""   ) );
  auxI->params.push_back( ReportParamsOption("--spectrum",      "<1>",                                      ""   ) );
  auxI->params.push_back( ReportParamsOption("--peptide",       "<7>",                                      ""   ) );
  auxI->params.push_back( ReportParamsOption("--aa",            "<aas_file>",                               "<aas_file>" ) );
  auxI->params.push_back( ReportParamsOption("--target",        "internal",                                 ""   ) );
  auxI->params.push_back( ReportParamsOption("--encoding",      "uu64",                                     ""   ) );
  auxI->params.push_back( ReportParamsOption("--title",         "Spectrum Scan <2>",                        ""   ) );

  auxB->sequences.push_back(auxI);


  auxI = new ReportColumnTypeImageOnDemand();
  auxI->columnLabel  = "Homolog";
  auxI->label        = "<8>";
  auxI->renderer     = "specplot";
  auxI->linkDisplayLevel  = 3;
  //auxI->params       = "--pklbin <projectdir>/<6>  --spectrum <1> --peptide <8> --title \"Spectrum Scan <2>\" --target internal --encoding uu64";
  auxI->validator    = "<8>";
  auxI->splitLabel   = true;

  auxI->files.push_back( ReportParamsFiles(1,   SPS_FILE_SPECSET, "<5>",   ""    ) );

  //auxI->params.push_back( ReportParamsOption("--pklbin",        "<projectdir>/<6>",                         ""   ) );
  //auxI->params.push_back( ReportParamsOption("--pklbin",        "<6>",                         ""   ) );
  auxI->params.push_back( ReportParamsOption("--spectrum",      "<1>",                                      ""   ) );
  auxI->params.push_back( ReportParamsOption("--peptide",       "<8>",                                      ""   ) );
  auxI->params.push_back( ReportParamsOption("--aa",            "<aas_file>",                               "<aas_file>" ) );
  auxI->params.push_back( ReportParamsOption("--target",        "internal",                                 ""   ) );
  auxI->params.push_back( ReportParamsOption("--encoding",      "uu64",                                     ""   ) );
  auxI->params.push_back( ReportParamsOption("--title",         "Spectrum Scan <2>",                        ""   ) );

  auxB->sequences.push_back(auxI);


  auxI = new ReportColumnTypeImageOnDemand();
  auxI->columnLabel  = "de Novo";
  auxI->label        = "<9>";
  auxI->renderer     = "specplot";
  auxI->linkDisplayLevel  = 3;
  //auxI->params       = "--pklbin <projectdir>/<6>  --spectrum <1> --peptide <9> --title \"Spectrum Scan <2>\" --target internal --encoding uu64";
  auxI->validator    = "<9>";
  auxI->splitLabel   = true;

  auxI->files.push_back( ReportParamsFiles(1,   SPS_FILE_SPECSET, "<5>",   ""    ) );

  //auxI->params.push_back( ReportParamsOption("--pklbin",        "<projectdir>/<6>",                         ""   ) );
  //auxI->params.push_back( ReportParamsOption("--pklbin",        "<6>",                         ""   ) );
  auxI->params.push_back( ReportParamsOption("--spectrum",      "<1>",                                      ""   ) );
  auxI->params.push_back( ReportParamsOption("--peptide",       "<9>",                                      ""   ) );
  auxI->params.push_back( ReportParamsOption("--aa",            "<aas_file>",                               "<aas_file>" ) );
  auxI->params.push_back( ReportParamsOption("--target",        "internal",                                 ""   ) );
  auxI->params.push_back( ReportParamsOption("--encoding",      "uu64",                                     ""   ) );
  auxI->params.push_back( ReportParamsOption("--title",         "Spectrum Scan <2>",                        ""   ) );

  auxB->sequences.push_back(auxI);


  auxI = new ReportColumnTypeImageOnDemand();
  auxI->dynamic      = true;
  auxI->columnLabel  = "User";
  auxI->label        = "<10>";
  auxI->renderer     = "specplot";
  auxI->linkDisplayLevel  = 3;
  //auxI->params       = "--pklbin <projectdir>/<6>  --spectrum <1> --peptide <10> --title \"Spectrum Scan <2>\" --target internal --encoding uu64";
  auxI->validator    = "<10>";
  auxI->splitLabel   = true;

  auxI->files.push_back( ReportParamsFiles(1,   SPS_FILE_SPECSET, "<5>",   ""    ) );

  //auxI->params.push_back( ReportParamsOption("--pklbin",        "<projectdir>/<6>",                         ""   ) );
  //auxI->params.push_back( ReportParamsOption("--pklbin",        "<6>",                         ""   ) );
  auxI->params.push_back( ReportParamsOption("--spectrum",      "<1>",                                      ""   ) );
  auxI->params.push_back( ReportParamsOption("--peptide",       "<10>",                                     ""   ) );
  auxI->params.push_back( ReportParamsOption("--aa",            "<aas_file>",                               "<aas_file>" ) );
  auxI->params.push_back( ReportParamsOption("--target",        "internal",                                 ""   ) );
  auxI->params.push_back( ReportParamsOption("--encoding",      "uu64",                                     ""   ) );
  auxI->params.push_back( ReportParamsOption("--title",         "Spectrum Scan <2>",                        ""   ) );

  auxB->sequences.push_back(auxI);


  auxS = new ReportColumnTypeString();
  auxS->dynamic      = true;
  auxS->isInput      = true;
  auxS->id           = "input_<row>_<col>";
  auxB->sequences.push_back(auxS);

  auxS = new ReportColumnTypeString();
  auxS->dynamic      = true;
  auxS->isButton     = true;
  auxS->text         = "Update";
  auxS->link         = "--update <0> --column 10 --value <10>";
  auxS->onClick      = "javascript:DoOnCick(\"input_<row>_<col>\", n, 10);";
  auxB->sequences.push_back(auxS);

  m_colTypes.push_back(auxB);


  // colTypes[3] -> (CTstring) mass
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Protein";
  auxS->link         = "protein.<18>.0.html";
  auxS->text         = "<4>";
  m_colTypes.push_back(auxS);

  // colTypes[3] -> (CTstring) mass
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Contig";
  auxS->link         = "contig.<17>.0.html";
  auxS->text         = "<17>";
  m_colTypes.push_back(auxS);

  // colTypes[3] -> (CTstring) mass
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Cluster";
  auxS->link         = "cluster.<3>.0.html";
  auxS->text         = "<3>";
  m_colTypes.push_back(auxS);

  // colTypes[3] -> (CTstring) mass
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Mass (m)";
  auxS->text         = "<11>";
  m_colTypes.push_back(auxS);

  // colTypes[4] -> (CTstring) charge
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Charge (z)";
  auxS->text         = "<12>";
  m_colTypes.push_back(auxS);

  // colTypes[5] -> (CTstring) B%
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "B (%)";
  auxS->text         = "<13>";
  m_colTypes.push_back(auxS);

  // colTypes[6] -> (CTstring) Y%
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Y (%)";
  auxS->text         = "<14>";
  m_colTypes.push_back(auxS);

  // colTypes[7] -> (CTstring) BY intensity %
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "BY intensity (%)";
  auxS->text         = "<15>";
  m_colTypes.push_back(auxS);

  // colTypes[3] -> (CTstring) mass
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Tool";
  auxS->text         = "<19>";
  m_colTypes.push_back(auxS);
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableInputSpectra::defineView2(void)
{
  // clear view
  clearView();

  // auxiliary object for column type holding
  ReportColumnTypeString        *auxS;
  ReportColumnTypeImageOnDemand *auxI;
  ReportColumnTypeBox           *auxB;

  // colTypes[0] -> (CTstring) Spectrum index
  auxS = new ReportColumnTypeString();
  auxS->columnLabel   = "Spectrum<br>Index";
  auxS->text          = "<1>";
  m_colTypes.push_back(auxS);

  auxS = new ReportColumnTypeString();
  auxS->columnLabel   = "Spectrum<br>Scan";
  auxS->text          = "<2>";
  m_colTypes.push_back(auxS);

// colTypes[8] -> (CTstring) Model %
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Model";
  auxS->text         = "<20>";
  m_colTypes.push_back(auxS);

  // colTypes[1] -> (CTstring) Spectrum file name
  //auxS = new ReportColumnTypeString();
  //auxS->columnLabel   = "File name";
  //auxS->text          = "<6>";
  //auxS->link          = "cluster.<3>.0.html";
  //m_colTypes.push_back(auxS);

  // colTypes[1] -> (CTseqsBox)
  auxB = new ReportColumnTypeBox();
  auxB->columnLabel   = "Spectrum";
  //auxB->link          = "cluster.<3>.0.html";

  auxI = new ReportColumnTypeImageOnDemand();
  auxI->iconRenderer      = "specplot";
  auxI->iconDisplayLevel  = 3;
  auxI->alt               = "N/A";
  //auxI->iconParams        = "--pklbin <projectdir>/<6>  --spectrum <1> --peptide <8> --notitle --zoom 0.3 --target internal --encoding uu64";
  auxI->validator         = "<7|8|9>";

  auxI->files.push_back( ReportParamsFiles(1,   SPS_FILE_SPECSET, "<5>",   ""    ) );

  //auxI->iconParams.push_back( ReportParamsOption("--pklbin",        "<projectdir>/<6>",                         ""   ) );
  //auxI->iconParams.push_back( ReportParamsOption("--pklbin",        "<6>",                         ""   ) );
  auxI->iconParams.push_back( ReportParamsOption("--spectrum",      "<1>",                                      ""   ) );
  auxI->iconParams.push_back( ReportParamsOption("--peptide",       "<7|8|9>",                                  ""   ) );
  auxI->iconParams.push_back( ReportParamsOption("--aa",            "<aas_file>",                               "<aas_file>" ) );
  auxI->iconParams.push_back( ReportParamsOption("--target",        "internal",                                 ""   ) );
  auxI->iconParams.push_back( ReportParamsOption("--encoding",      "uu64",                                     ""   ) );
  auxI->iconParams.push_back( ReportParamsOption("--zoom",          "0.35",                                     ""   ) );
  auxI->iconParams.push_back( ReportParamsOption("--notitle",       "",                                         ""   ) );

  auxB->sequences.push_back(auxI);


  auxS = new ReportColumnTypeString();
  auxS->text          = "<16>";
  auxS->splitText     = true;
  auxB->sequences.push_back(auxS);

  m_colTypes.push_back(auxB);


  // colTypes[2] -> (CTseqsBox)
  auxB = new ReportColumnTypeBox();
  auxB->columnLabel   = "Sequences";

  auxI = new ReportColumnTypeImageOnDemand();
  auxI->columnLabel       = "Reference";
  auxI->label             = "<7>";
  auxI->renderer          = "specplot";
  auxI->linkDisplayLevel  = 3;
  //auxI->params            = "--pklbin <projectdir>/<6>  --spectrum <1> --peptide <7> --title \"Spectrum Scan <2>\" --target internal --encoding uu64";
  auxI->validator         = "<7>";
  auxI->splitLabel        = true;

  auxI->files.push_back( ReportParamsFiles(1,   SPS_FILE_SPECSET, "<5>",   ""    ) );

  //auxI->params.push_back( ReportParamsOption("--pklbin",        "<projectdir>/<6>",                         ""   ) );
  //auxI->params.push_back( ReportParamsOption("--pklbin",        "<6>",                         ""   ) );
  auxI->params.push_back( ReportParamsOption("--spectrum",      "<1>",                                      ""   ) );
  auxI->params.push_back( ReportParamsOption("--peptide",       "<7>",                                      ""   ) );
  auxI->params.push_back( ReportParamsOption("--aa",            "<aas_file>",                               "<aas_file>" ) );
  auxI->params.push_back( ReportParamsOption("--target",        "internal",                                 ""   ) );
  auxI->params.push_back( ReportParamsOption("--encoding",      "uu64",                                     ""   ) );
  auxI->params.push_back( ReportParamsOption("--title",         "Spectrum Scan <2>",                        ""   ) );

  auxB->sequences.push_back(auxI);


  auxI = new ReportColumnTypeImageOnDemand();
  auxI->columnLabel       = "Homolog";
  auxI->label             = "<8>";
  auxI->renderer          = "specplot";
  auxI->linkDisplayLevel  = 3;
  //auxI->params            = "--pklbin <projectdir>/<6>  --spectrum <1> --peptide <8> --title \"Spectrum Scan <2>\" --target internal --encoding uu64";
  auxI->validator         = "<8>";
  auxI->splitLabel        = true;

  auxI->files.push_back( ReportParamsFiles(1,   SPS_FILE_SPECSET, "<5>",   ""    ) );

  //auxI->params.push_back( ReportParamsOption("--pklbin",        "<projectdir>/<6>",                         ""   ) );
  //auxI->params.push_back( ReportParamsOption("--pklbin",        "<6>",                         ""   ) );
  auxI->params.push_back( ReportParamsOption("--spectrum",      "<1>",                                      ""   ) );
  auxI->params.push_back( ReportParamsOption("--peptide",       "<8>",                                      ""   ) );
  auxI->params.push_back( ReportParamsOption("--aa",            "<aas_file>",                               "<aas_file>" ) );
  auxI->params.push_back( ReportParamsOption("--target",        "internal",                                 ""   ) );
  auxI->params.push_back( ReportParamsOption("--encoding",      "uu64",                                     ""   ) );
  auxI->params.push_back( ReportParamsOption("--title",         "Spectrum Scan <2>",                        ""   ) );

  auxB->sequences.push_back(auxI);


  auxI = new ReportColumnTypeImageOnDemand();
  auxI->columnLabel       = "de Novo";
  auxI->label             = "<9>";
  auxI->renderer          = "specplot";
  auxI->linkDisplayLevel  = 3;
  //auxI->params            = "--pklbin <projectdir>/<6>  --spectrum <1> --peptide <9> --title \"Spectrum Scan <2>\" --target internal --encoding uu64";
  auxI->validator         = "<9>";
  auxI->splitLabel        = true;

  auxI->files.push_back( ReportParamsFiles(1,   SPS_FILE_SPECSET, "<5>",   ""    ) );

  //auxI->params.push_back( ReportParamsOption("--pklbin",        "<projectdir>/<6>",                         ""   ) );
  //auxI->params.push_back( ReportParamsOption("--pklbin",        "<6>",                         ""   ) );
  auxI->params.push_back( ReportParamsOption("--spectrum",      "<1>",                                      ""   ) );
  auxI->params.push_back( ReportParamsOption("--peptide",       "<9>",                                      ""   ) );
  auxI->params.push_back( ReportParamsOption("--aa",            "<aas_file>",                               "<aas_file>" ) );
  auxI->params.push_back( ReportParamsOption("--target",        "internal",                                 ""   ) );
  auxI->params.push_back( ReportParamsOption("--encoding",      "uu64",                                     ""   ) );
  auxI->params.push_back( ReportParamsOption("--title",         "Spectrum Scan <2>",                        ""   ) );

  auxB->sequences.push_back(auxI);


  auxI = new ReportColumnTypeImageOnDemand();
  auxI->dynamic           = true;
  auxI->columnLabel       = "User";
  auxI->label             = "<10>";
  auxI->renderer          = "specplot";
  auxI->linkDisplayLevel  = 3;
  //auxI->params            = "--pklbin <projectdir>/<6>  --spectrum <1> --peptide <10> --title \"Spectrum Scan <2>\" --target internal --encoding uu64";
  auxI->validator         = "<10>";
  auxI->splitLabel        = true;

  auxI->files.push_back( ReportParamsFiles(1,   SPS_FILE_SPECSET, "<5>",   ""    ) );

  //auxI->params.push_back( ReportParamsOption("--pklbin",        "<projectdir>/<6>",                         ""   ) );
  //auxI->params.push_back( ReportParamsOption("--pklbin",        "<6>",                         ""   ) );
  auxI->params.push_back( ReportParamsOption("--spectrum",      "<1>",                                      ""   ) );
  auxI->params.push_back( ReportParamsOption("--peptide",       "<10>",                                     ""   ) );
  auxI->params.push_back( ReportParamsOption("--aa",            "<aas_file>",                               "<aas_file>" ) );
  auxI->params.push_back( ReportParamsOption("--target",        "internal",                                 ""   ) );
  auxI->params.push_back( ReportParamsOption("--encoding",      "uu64",                                     ""   ) );
  auxI->params.push_back( ReportParamsOption("--title",         "Spectrum Scan <2>",                        ""   ) );

  auxB->sequences.push_back(auxI);


  auxS = new ReportColumnTypeString();
  auxS->dynamic      = true;
  auxS->isInput      = true;
  auxS->id           = "input_<row>_<col>";
  auxB->sequences.push_back(auxS);

  auxS = new ReportColumnTypeString();
  auxS->dynamic      = true;
  auxS->isButton     = true;
  auxS->text         = "Update";
  auxS->link         = "--update <0> --column 10 --value <10>";
  auxS->onClick      = "javascript:DoOnCick(\"input_<row>_<col>\", n, 10);";
  auxB->sequences.push_back(auxS);

  m_colTypes.push_back(auxB);


  // colTypes[3] -> (CTstring) protein name
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Protein";
  auxS->link         = "protein.<18>.0.html";
  auxS->text         = "<4>";
  m_colTypes.push_back(auxS);

  // colTypes[3] -> (CTstring) mass
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Contig";
  auxS->link         = "contig.<17>.0.html";
  auxS->text         = "<17>";
  m_colTypes.push_back(auxS);

  // colTypes[3] -> (CTstring) mass
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Cluster";
  auxS->link         = "cluster.<3>.0.html";
  auxS->text         = "<3>";
  m_colTypes.push_back(auxS);

  // colTypes[3] -> (CTstring) mass
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Mass (m)";
  auxS->text         = "<11>";
  m_colTypes.push_back(auxS);

  // colTypes[4] -> (CTstring) charge
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Charge (z)";
  auxS->text         = "<12>";
  m_colTypes.push_back(auxS);

  // colTypes[5] -> (CTstring) B%
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "B (%)";
  auxS->text         = "<13>";
  m_colTypes.push_back(auxS);

  // colTypes[6] -> (CTstring) Y%
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Y (%)";
  auxS->text         = "<14>";
  m_colTypes.push_back(auxS);

  // colTypes[7] -> (CTstring) BY intensity %
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "BY intensity (%)";
  auxS->text         = "<15>";
  m_colTypes.push_back(auxS);

  // colTypes[3] -> (CTstring) mass
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Tool";
  auxS->text         = "<19>";
  m_colTypes.push_back(auxS);
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableInputSpectra::defineViewNoClusters(void)
{
  // clear view
  clearView();

  // auxiliary object for column type holding
  ReportColumnTypeString        *auxS;
  ReportColumnTypeImageOnDemand *auxI;
  ReportColumnTypeBox           *auxB;

  // colTypes[0] -> (CTstring) Spectrum index
  auxS = new ReportColumnTypeString();
  auxS->columnLabel   = "Spectrum<br>Index";
  auxS->text          = "<1>";
  m_colTypes.push_back(auxS);

  auxS = new ReportColumnTypeString();
  auxS->columnLabel   = "Spectrum<br>Scan";
  auxS->text          = "<2>";
  m_colTypes.push_back(auxS);

// colTypes[8] -> (CTstring) Model %
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Model";
  auxS->text         = "<20>";
  m_colTypes.push_back(auxS);

  // colTypes[1] -> (CTseqsBox)
  auxB = new ReportColumnTypeBox();
  auxB->columnLabel   = "Spectrum";
  //auxB->link          = "cluster.<3>.0.html";

  auxI = new ReportColumnTypeImageOnDemand();
  auxI->iconRenderer      = "specplot";
  auxI->iconDisplayLevel  = 3;
  auxI->alt               = "N/A";
  //auxI->iconParams        = "--pklbin <projectdir>/<6>  --spectrum <1> --peptide <8> --notitle --zoom 0.3 --target internal --encoding uu64";
  auxI->validator         = "<7|8|9>";

  auxI->files.push_back( ReportParamsFiles(1,   SPS_FILE_SPECSET, "<5>",   ""    ) );

  //auxI->iconParams.push_back( ReportParamsOption("--pklbin",        "<projectdir>/<6>",                         ""   ) );
  //auxI->iconParams.push_back( ReportParamsOption("--pklbin",        "<6>",                         ""   ) );
  auxI->iconParams.push_back( ReportParamsOption("--spectrum",      "<1>",                                      ""   ) );
  auxI->iconParams.push_back( ReportParamsOption("--peptide",       "<7|8|9>",                                  ""   ) );
  auxI->iconParams.push_back( ReportParamsOption("--aa",            "<aas_file>",                               "<aas_file>" ) );
  auxI->iconParams.push_back( ReportParamsOption("--target",        "internal",                                 ""   ) );
  auxI->iconParams.push_back( ReportParamsOption("--encoding",      "uu64",                                     ""   ) );
  auxI->iconParams.push_back( ReportParamsOption("--zoom",          "0.35",                                     ""   ) );
  auxI->iconParams.push_back( ReportParamsOption("--notitle",       "",                                         ""   ) );

  auxB->sequences.push_back(auxI);


  auxS = new ReportColumnTypeString();
  auxS->text          = "<16>";
  auxS->splitText     = true;
  auxB->sequences.push_back(auxS);

  m_colTypes.push_back(auxB);


  // colTypes[2] -> (CTseqsBox)
  auxB = new ReportColumnTypeBox();
  auxB->columnLabel   = "Sequences";

  auxI = new ReportColumnTypeImageOnDemand();
  auxI->columnLabel       = "Reference";
  auxI->label             = "<7>";
  auxI->renderer          = "specplot";
  auxI->linkDisplayLevel  = 3;
  //auxI->params            = "--pklbin <projectdir>/<6>  --spectrum <1> --peptide <7> --title \"Spectrum Scan <2>\" --target internal --encoding uu64";
  auxI->validator         = "<7>";
  auxI->splitLabel        = true;

  auxI->files.push_back( ReportParamsFiles(1,   SPS_FILE_SPECSET, "<5>",   ""    ) );

  //auxI->params.push_back( ReportParamsOption("--pklbin",        "<projectdir>/<6>",                         ""   ) );
  //auxI->params.push_back( ReportParamsOption("--pklbin",        "<6>",                         ""   ) );
  auxI->params.push_back( ReportParamsOption("--spectrum",      "<1>",                                      ""   ) );
  auxI->params.push_back( ReportParamsOption("--peptide",       "<7>",                                      ""   ) );
  auxI->params.push_back( ReportParamsOption("--aa",            "<aas_file>",                               "<aas_file>" ) );
  auxI->params.push_back( ReportParamsOption("--target",        "internal",                                 ""   ) );
  auxI->params.push_back( ReportParamsOption("--encoding",      "uu64",                                     ""   ) );
  auxI->params.push_back( ReportParamsOption("--title",         "Spectrum Scan <2>",                        ""   ) );

  auxB->sequences.push_back(auxI);


  auxI = new ReportColumnTypeImageOnDemand();
  auxI->columnLabel       = "Homolog";
  auxI->label             = "<8>";
  auxI->renderer          = "specplot";
  auxI->linkDisplayLevel  = 3;
  //auxI->params            = "--pklbin <projectdir>/<6>  --spectrum <1> --peptide <8> --title \"Spectrum Scan <2>\" --target internal --encoding uu64";
  auxI->validator         = "<8>";
  auxI->splitLabel        = true;

  auxI->files.push_back( ReportParamsFiles(1,   SPS_FILE_SPECSET, "<5>",   ""    ) );

  //auxI->params.push_back( ReportParamsOption("--pklbin",        "<projectdir>/<6>",                         ""   ) );
  //auxI->params.push_back( ReportParamsOption("--pklbin",        "<6>",                         ""   ) );
  auxI->params.push_back( ReportParamsOption("--spectrum",      "<1>",                                      ""   ) );
  auxI->params.push_back( ReportParamsOption("--peptide",       "<8>",                                      ""   ) );
  auxI->params.push_back( ReportParamsOption("--aa",            "<aas_file>",                               "<aas_file>" ) );
  auxI->params.push_back( ReportParamsOption("--target",        "internal",                                 ""   ) );
  auxI->params.push_back( ReportParamsOption("--encoding",      "uu64",                                     ""   ) );
  auxI->params.push_back( ReportParamsOption("--title",         "Spectrum Scan <2>",                        ""   ) );

  auxB->sequences.push_back(auxI);


  auxI = new ReportColumnTypeImageOnDemand();
  auxI->columnLabel       = "de Novo";
  auxI->label             = "<9>";
  auxI->renderer          = "specplot";
  auxI->linkDisplayLevel  = 3;
  //auxI->params            = "--pklbin <projectdir>/<6>  --spectrum <1> --peptide <9> --title \"Spectrum Scan <2>\" --target internal --encoding uu64";
  auxI->validator         = "<9>";
  auxI->splitLabel        = true;

  auxI->files.push_back( ReportParamsFiles(1,   SPS_FILE_SPECSET, "<5>",   ""    ) );

  //auxI->params.push_back( ReportParamsOption("--pklbin",        "<projectdir>/<6>",                         ""   ) );
  //auxI->params.push_back( ReportParamsOption("--pklbin",        "<6>",                         ""   ) );
  auxI->params.push_back( ReportParamsOption("--spectrum",      "<1>",                                      ""   ) );
  auxI->params.push_back( ReportParamsOption("--peptide",       "<9>",                                      ""   ) );
  auxI->params.push_back( ReportParamsOption("--aa",            "<aas_file>",                               "<aas_file>" ) );
  auxI->params.push_back( ReportParamsOption("--target",        "internal",                                 ""   ) );
  auxI->params.push_back( ReportParamsOption("--encoding",      "uu64",                                     ""   ) );
  auxI->params.push_back( ReportParamsOption("--title",         "Spectrum Scan <2>",                        ""   ) );

  auxB->sequences.push_back(auxI);


  auxI = new ReportColumnTypeImageOnDemand();
  auxI->dynamic           = true;
  auxI->columnLabel       = "User";
  auxI->label             = "<10>";
  auxI->renderer          = "specplot";
  auxI->linkDisplayLevel  = 3;
  //auxI->params            = "--pklbin <projectdir>/<6>  --spectrum <1> --peptide <10> --title \"Spectrum Scan <2>\" --target internal --encoding uu64";
  auxI->validator         = "<10>";
  auxI->splitLabel        = true;

  auxI->files.push_back( ReportParamsFiles(1,   SPS_FILE_SPECSET, "<5>",   ""    ) );

  //auxI->params.push_back( ReportParamsOption("--pklbin",        "<projectdir>/<6>",                         ""   ) );
  //auxI->params.push_back( ReportParamsOption("--pklbin",        "<6>",                         ""   ) );
  auxI->params.push_back( ReportParamsOption("--spectrum",      "<1>",                                      ""   ) );
  auxI->params.push_back( ReportParamsOption("--peptide",       "<10>",                                     ""   ) );
  auxI->params.push_back( ReportParamsOption("--aa",            "<aas_file>",                               "<aas_file>" ) );
  auxI->params.push_back( ReportParamsOption("--target",        "internal",                                 ""   ) );
  auxI->params.push_back( ReportParamsOption("--encoding",      "uu64",                                     ""   ) );
  auxI->params.push_back( ReportParamsOption("--title",         "Spectrum Scan <2>",                        ""   ) );

  auxB->sequences.push_back(auxI);


  auxS = new ReportColumnTypeString();
  auxS->dynamic      = true;
  auxS->isInput      = true;
  auxS->id           = "input_<row>_<col>";
  auxB->sequences.push_back(auxS);

  auxS = new ReportColumnTypeString();
  auxS->dynamic      = true;
  auxS->isButton     = true;
  auxS->text         = "Update";
  auxS->link         = "--update <0> --column 10 --value <10>";
  auxS->onClick      = "javascript:DoOnCick(\"input_<row>_<col>\", n, 10);";
  auxB->sequences.push_back(auxS);

  m_colTypes.push_back(auxB);


  // colTypes[3] -> (CTstring) protein name
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Protein";
  auxS->link         = "protein.<18>.0.html";
  auxS->text         = "<4>";
  m_colTypes.push_back(auxS);

  // colTypes[3] -> (CTstring) mass
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Contig";
  auxS->link         = "contig.<17>.0.html";
  auxS->text         = "<17>";
  m_colTypes.push_back(auxS);

  // colTypes[3] -> (CTstring) mass
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Cluster";
  auxS->link         = "cluster.<3>.0.html";
  auxS->text         = "<3>";
  m_colTypes.push_back(auxS);

  // colTypes[3] -> (CTstring) mass
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Mass (m)";
  auxS->text         = "<11>";
  m_colTypes.push_back(auxS);

  // colTypes[4] -> (CTstring) charge
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Charge (z)";
  auxS->text         = "<12>";
  m_colTypes.push_back(auxS);

  // colTypes[5] -> (CTstring) B%
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "B (%)";
  auxS->text         = "<13>";
  m_colTypes.push_back(auxS);

  // colTypes[6] -> (CTstring) Y%
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Y (%)";
  auxS->text         = "<14>";
  m_colTypes.push_back(auxS);

  // colTypes[7] -> (CTstring) BY intensity %
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "BY intensity (%)";
  auxS->text         = "<15>";
  m_colTypes.push_back(auxS);

  // colTypes[3] -> (CTstring) mass
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Tool";
  auxS->text         = "<19>";
  m_colTypes.push_back(auxS);
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableInputSpectra::defineViewImages(void)
{
  // clear view
  clearView();

  // auxiliary object for column type holding
  ReportColumnTypeString        *auxS;
  ReportColumnTypeImageOnDemand *auxI;
  ReportColumnTypeBox           *auxB;

  auxI = new ReportColumnTypeImageOnDemand();
  auxI->id                = "ss<0>";
  auxI->iconRenderer      = "specplot";
  auxI->iconDisplayLevel  = 3;
  auxI->validator         = "<7|8|9>";
  auxI->files.push_back( ReportParamsFiles(1,   SPS_FILE_SPECSET, "<5>",   ""    ) );
  auxI->iconParams.push_back( ReportParamsOption("--spectrum",      "<1>",                                      ""   ) );
  auxI->iconParams.push_back( ReportParamsOption("--peptide",       "<7|8|9>",                                  ""   ) );
  auxI->iconParams.push_back( ReportParamsOption("--aa",            "<aas_file>",                               "<aas_file>" ) );
  auxI->iconParams.push_back( ReportParamsOption("--target",        "internal",                                 ""   ) );
  auxI->iconParams.push_back( ReportParamsOption("--encoding",      "uu64",                                     ""   ) );
  auxI->iconParams.push_back( ReportParamsOption("--zoom",          "0.35",                                     ""   ) );
  auxI->iconParams.push_back( ReportParamsOption("--notitle",       "",                                         ""   ) );
  m_colTypes.push_back(auxI);

  auxI = new ReportColumnTypeImageOnDemand();
  auxI->id                = "sR<0>";
  auxI->renderer          = "specplot";
  auxI->linkDisplayLevel  = 3;
  auxI->validator         = "<7>";
  auxI->files.push_back( ReportParamsFiles(1,   SPS_FILE_SPECSET, "<5>",   ""    ) );
  auxI->params.push_back( ReportParamsOption("--spectrum",      "<1>",                                      ""   ) );
  auxI->params.push_back( ReportParamsOption("--peptide",       "<7>",                                      ""   ) );
  auxI->params.push_back( ReportParamsOption("--aa",            "<aas_file>",                               "<aas_file>" ) );
  auxI->params.push_back( ReportParamsOption("--target",        "internal",                                 ""   ) );
  auxI->params.push_back( ReportParamsOption("--encoding",      "uu64",                                     ""   ) );
  auxI->params.push_back( ReportParamsOption("--title",         "Spectrum Scan <2>",                        ""   ) );
  m_colTypes.push_back(auxI);

  auxI = new ReportColumnTypeImageOnDemand();
  auxI->id                = "sH<0>";
  auxI->renderer          = "specplot";
  auxI->linkDisplayLevel  = 3;
  auxI->validator         = "<8>";
  auxI->files.push_back( ReportParamsFiles(1,   SPS_FILE_SPECSET, "<5>",   ""    ) );
  auxI->params.push_back( ReportParamsOption("--spectrum",      "<1>",                                      ""   ) );
  auxI->params.push_back( ReportParamsOption("--peptide",       "<8>",                                      ""   ) );
  auxI->params.push_back( ReportParamsOption("--aa",            "<aas_file>",                               "<aas_file>" ) );
  auxI->params.push_back( ReportParamsOption("--target",        "internal",                                 ""   ) );
  auxI->params.push_back( ReportParamsOption("--encoding",      "uu64",                                     ""   ) );
  auxI->params.push_back( ReportParamsOption("--title",         "Spectrum Scan <2>",                        ""   ) );
  m_colTypes.push_back(auxI);

  auxI = new ReportColumnTypeImageOnDemand();
  auxI->id                = "sN<0>";
  auxI->renderer          = "specplot";
  auxI->linkDisplayLevel  = 3;
  auxI->validator         = "<9>";
  auxI->files.push_back( ReportParamsFiles(1,   SPS_FILE_SPECSET, "<5>",   ""    ) );
  auxI->params.push_back( ReportParamsOption("--spectrum",      "<1>",                                      ""   ) );
  auxI->params.push_back( ReportParamsOption("--peptide",       "<9>",                                      ""   ) );
  auxI->params.push_back( ReportParamsOption("--aa",            "<aas_file>",                               "<aas_file>" ) );
  auxI->params.push_back( ReportParamsOption("--target",        "internal",                                 ""   ) );
  auxI->params.push_back( ReportParamsOption("--encoding",      "uu64",                                     ""   ) );
  auxI->params.push_back( ReportParamsOption("--title",         "Spectrum Scan <2>",                        ""   ) );
  m_colTypes.push_back(auxI);
}
////////////////////////////////////////////////////////////////////////////////
}; // namespace
////////////////////////////////////////////////////////////////////////////////
