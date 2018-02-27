///////////////////////////////////////////////////////////////////////////////
#ifndef __SPS_TIMER_H__
#define __SPS_TIMER_H__
///////////////////////////////////////////////////////////////////////////////

#include <ctime>

///////////////////////////////////////////////////////////////////////////////

/*! \brief Timer class

  Provides timer object to measure time passed between code locations.

 */
class Timer_c {

/*! \brief t1 is the initial time; t2 the final time.
 */
  clock_t t1, t2;

/*! \brief True if the clock is running
 */
  bool running;

 public:

    //! \name CONSTRUCTORS
    //@{
    /*! \brief The exemplar constructor.

      Default contructor
     */
    Timer_c() : running(true) {start();};
    //@}

    //! \name DESTRUCTOR
    //@{
    ~Timer_c() {};
    //@}

    /*! \brief Starts timer.

     Starts the timer.
     */
  void   start(void)   {t1 = clock(); running = true;};

    /*! \brief Restarts timer.

     Restarts timer, returning current time.

     @return Time elapsed.
     */
  double restart(void) {double ret = get(); t1 = clock(); return ret;};

    /*! \brief Stops the timer.

     Stops timer, returning current time.

     @return Time elapsed.
     */
  double stop(void)    {double ret = get(); running = false; return ret;};

    /*! \brief Get time.

     Returns current time.

     @return Time elapsed.
     */
 double get(void)     {t2 = clock(); if(running) return ((double)(t2 - t1) / CLOCKS_PER_SEC);return 0.0;};

};
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
