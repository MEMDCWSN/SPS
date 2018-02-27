//
//  main_specdump - stand alone executable for dumping spectrum data
//
#include "Logger.h"
#include "db_fasta.h"
#include <stdlib.h>

using namespace specnets;
using namespace std;

const size_t HISTO_SIZE = 750;

// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));

  DEBUG_TRACE;

  if (argc != 4)
  {
    cerr << "Usage: main_createdecoydb input_db output_db type" << endl;
    cerr << "       types are: reverse, shuffled shuffledinplace" << endl;
    return -1;
  }

  DB_fasta db;
  db.Load(argv[1]);

  if (strcmp(argv[3], "reverse") == 0) {
    DEBUG_TRACE;
    db.addDecoyReversed();
  } else if (strcmp(argv[3], "shuffled") == 0) {
    DEBUG_TRACE;
    db.addDecoyShuffled();
  } else if (strcmp(argv[3], "shuffledinplace") == 0) {
    DEBUG_TRACE;
    db.replaceDecoyShuffled();
  } else {
    cerr << "Unknown type: " << argv[3] << endl;
    return -1;
  }

  DEBUG_TRACE;
  db.Save(argv[2]);
  
  return 0;
}

