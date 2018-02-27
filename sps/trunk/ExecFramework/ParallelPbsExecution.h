#ifndef __ParallelPbsExecution_H__
#define __ParallelPbsExecution_H__

// Module Includes
#include "ParallelExecution.h"

// System Includes
#include <pthread.h>
#include <vector>

namespace specnets
{
  /*! \brief Executes a specnets module in parallel using the Sun Grid Engine

   */
  class ParallelPbsExecution : public ParallelExecution
  {
  public:

    //! \name CONSTRUCTORS
    //@{
    /*! \brief The constructor

     An execution module is passed to the constructor and it is this module that
     will be executed when invoke is called.
     @sa invoke()
     @param moduleExec the module that will be executed
     */
    ParallelPbsExecution(ExecBase * moduleExec,
                         bool gridExecution,
                         bool manual,
                         bool resume);
    //@}

    //! \name DESTRUCTOR
    //@{
    ~ParallelPbsExecution(void);
    //@}

    //! \name MODIFIERS
    //@{
    /*! \brief Executes the module in parallel on separate grid nodes

     Calls methods on the original module to run portions of that module
     in parallel on separate grid nodes. The number of grid nodes used will be
     equal to the numSplits param. The split() method of the original module
     will be used to obtain the sub-modules for execution. Since all operations
     are performed in memory no calls to the original module to save or load data
     to files are necessary in this implementation. When all threads have completed
     and all batches have been run, the merge() method of the original module
     is called to merge all the sub-module back into a choesive result set.

     @param numSplits The number of sub-modules that will be used to execute the module
     @return True if execution of all sub-modules finished without error, false otherwise.
     */
    virtual bool invoke(int numSplits);
    //@}

  protected:

  private:
    bool m_gridExecution; //! Indicates execution will be on the actual grid
    bool m_manual; //! Indicates whether or not files will be copied manually to the SGE grid
    bool m_resume; //! Indicates (if manual) if resuming after manual SGE execution

    //! \name NOT IMPLEMENTED
    //@{
    ParallelPbsExecution(void);
    ParallelPbsExecution(const ParallelPbsExecution & that);
    //@}

  };

} // namespace specnets

#endif // __ParallelPbsExecution_H__
