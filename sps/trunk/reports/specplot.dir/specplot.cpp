///////////////////////////////////////////////////////////////////////////////
#include "specplot.h"

#include "SpecplotInterface2.h"
#include "PlotSpectrum.h"

#include "Timer.h"
///////////////////////////////////////////////////////////////////////////////
using namespace specnets;
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  //Timer_c t0;
  // Build PlotSpec object
  SpecplotInterface2  specplotInterface;

  // Process options
  int ret = specplotInterface.processOptions(argc, argv);
  // test return value
  if(ret != 0)
    return 0;
  // generate the image
  //cout << "Time to process parameters: " << t0.stop() << endl;
  //Timer_c t1;

  specplotInterface.plot();
  //cout << "Time to generate the image: " << t1.stop() << endl;
  // exit
  return 0;
}
///////////////////////////////////////////////////////////////////////////////
