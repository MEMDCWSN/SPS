#include "dbg_print.h"

#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>
#include <execinfo.h>
#include <signal.h>
#include <bfd.h>
#include <unistd.h>
#include <cxxabi.h>

#include <iostream>


using namespace std;


namespace
{

const unsigned MAX_FRAMES = 20;

/* globals retained across calls to resolve. */
bfd* abfd = 0;
asymbol **psym = 0;
asection *ptext = 0;

void resolve(char * paddress)
{
  if (! abfd)
  {
    char sname[1024];
    int nl = readlink("/proc/self/exe", sname, sizeof(sname));

    if (nl == -1)
    {
      perror("failed to find executable\n");
      return;
    }
    sname[nl] = 0;

    bfd_init();

    abfd = bfd_openr(sname, 0);
    if (! abfd)
    {
      perror("bfd_openr failed: ");
      return;
    }
    /* oddly, this is required for it to work... */
    bfd_check_format(abfd,bfd_object);

    unsigned storage_needed = bfd_get_symtab_upper_bound(abfd);
    psym = (asymbol **) malloc(storage_needed);
    unsigned cSymbols = bfd_canonicalize_symtab(abfd, psym);

    ptext = bfd_get_section_by_name(abfd, ".text");
  }

  long noffset = ((long)paddress) - ptext->vma;

  if (noffset > 0)
  {
    const char *sfile;
    const char *smangled;
    unsigned nline;

    if (bfd_find_nearest_line(abfd, ptext, psym, noffset, &sfile, &smangled, &nline) && sfile)
    {
      int nstatus;
      char * sdemangled = abi::__cxa_demangle(smangled, 0, 0, & nstatus);

      if (nstatus)
        note(cerr) << sfile << ": " << nline << ": " << smangled << endl;
      else
        note(cerr) << sfile << ": " << nline << ": " << sdemangled << endl;

      free(sdemangled);
    }
  }
}

} // namespace


void dbg_handler(int nsig)
{
  void * parray[MAX_FRAMES];
  size_t nsize;
  size_t i;
  void * pend = (void*) ((128+100) * 2<<20);

  nsize = backtrace(parray, MAX_FRAMES);

  for (i = 2; i < nsize; i++)
    if (parray[i] < pend)
      resolve((char*)parray[i]);
}
