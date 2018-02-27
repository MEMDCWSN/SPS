#ifndef __ParallelProcessExecution_H__
#define __ParallelProcessExecution_H__

// Module Includes
#include "ParallelExecution.h"

// System Includes
#include <pthread.h>
#include <vector>

namespace specnets
{
  /*! \brief Executes a specnets module in parallel using processes

   */
  class ParallelProcessExecution : public ParallelExecution
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
    ParallelProcessExecution(ExecBase * moduleExec);
    //@}

    //! \name DESTRUCTOR
    //@{
    ~ParallelProcessExecution(void);
    //@}

    //! \name MODIFIERS
    //@{
    /*! \brief Executes the module in parallel on separate threads

     Calls methods on the original module to run portions of that module
     in parallel on separate threads. The number of threads used will be
     equal to the numCPUs param while the overall execution will be run in
     numNodes batches. The split() method of the original module will be used
     to obtain the sub-modules for execution passing numNodes, and numCPUs
     directly to the module. Since all operations are performed in memory
     no calls to the original module to save or load data to files are
     necessary in this implementation. When all threads have completed and
     all batches have been run, the merge() method of the original module
     is called to merge all the sub-module back into a choesive result set.

     @param numSplits The number of sub-modules that will be used to execute the module
     @return True if execution of all sub-modules finished without error, false otherwise.
     */
    virtual bool invoke(int numSplits);
    //@}

  protected:

  private:
    //! \name NOT IMPLEMENTED
    //@{
    ParallelProcessExecution(void);
    ParallelProcessExecution(const ParallelProcessExecution & that);
    //@}

  };

} // namespace specnets

#endif // __ParallelProcessExecution_H__
