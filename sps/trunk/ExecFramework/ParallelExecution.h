#ifndef __ParallelExecution_H__
#define __ParallelExecution_H__

// Module Includes
#include "ExecBase.h"

namespace specnets
{
  /*! \brief Base class for running execution modules in parallel

   */
  class ParallelExecution
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
    ParallelExecution(ExecBase * moduleExec);
    //@}

    //! \name DESTRUCTOR
    //@{
    virtual ~ParallelExecution(void);
    //@}

    //! \name MODIFIERS
    //@{
    /*! \brief Virtual method that executes the module in parallel

     Performs all operations necessary to run sub-portions of the module in
     parallel. The specifics of this are left to the actual implementation in
     the derived class. However, the basic implementation should call the
     split() method of the module to divide the module into sub-units that
     may be executed in parallel. If necessary (as in the case where execution
     is performed by a separate process) data will be saved. All the sub-modules
     are run and then the results are meged back to the original module using
     the merge() method of the execution module (see ExecBase() for methods).

     @param numSplits The number of sub-modules that will be used to execute the module
     @return True if execution of all sub-modules finished without error, false otherwise.
     */
    virtual bool invoke(int numSplits) = 0;
    //@}

  protected:
    ExecBase * m_moduleExec;

  private:
    //! \name NOT IMPLEMENTED
    //@{
    ParallelExecution(void);
    ParallelExecution(const ParallelExecution & that);
    //@}
  };

} // namespace specnets

#endif // __ParallelExecution_H__
