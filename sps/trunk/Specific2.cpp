////////////////////////////////////////////////////////////////////////////////
// This file contains all functions that are, for some reason, system dependent
////////////////////////////////////////////////////////////////////////////////
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <dirent.h>
#include <iostream>
#include <fstream>

#include "Specific.h"
#include "Logger.h"
#include "StatusFile.h"

////////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace specnets;
////////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// Catch signals, like segfault
//-----------------------------------------------------------------------------
#if defined(__linux__) || defined(__CYGWIN__)
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void segfault_sigaction2(int signal, siginfo_t *si, void *arg)
{
  cerr << "Unexpected error. Aborting." << endl;
  //ERROR_MSG("Caught segfault at address " << si->si_addr);
  exit(-1);
}

#endif
//-----------------------------------------------------------------------------
// The segfault divert function: to where the segfaults are diverted to
//-----------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
int addSegFaultDivert2(void)
{
#if defined(__linux__) || defined(__CYGWIN__)
  // catch signals (redirect)
  struct sigaction sa;

  memset(&sa, 0, sizeof(struct sigaction));
  sigemptyset(&sa.sa_mask);
  sa.sa_sigaction = segfault_sigaction2;
  sa.sa_flags   = SA_SIGINFO;

  sigaction(SIGSEGV, &sa, NULL);
  return 1;
  // end catch signals (redirect)
#endif
  return 0;
}
////////////////////////////////////////////////////////////////////////////////
