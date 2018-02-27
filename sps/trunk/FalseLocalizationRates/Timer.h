/*
 * Timer.h
 *
 *  Created on: Jul 23, 2012
 *      Author: jsnedecor
 */

#ifndef TIMER_H_
#define TIMER_H_

#include <cstdlib>
#include <sys/time.h>


class Timer
{
    timeval timer[2];

  public:

    timeval start(void);

    timeval stop(void);

    int duration() const;
};
#endif /* TIMER_H_ */
