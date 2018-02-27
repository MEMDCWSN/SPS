/*
 * Timer.cpp
 *
 *  Created on: Jul 23, 2012
 *      Author: jsnedecor
 */

#include "Timer.h"

timeval Timer::start()
{
  gettimeofday(&this->timer[0], NULL);
  return this->timer[0];
}

timeval Timer::stop()
{
  gettimeofday(&this->timer[1], NULL);
  return this->timer[1];
}

int Timer::duration() const
{
  int secs(this->timer[1].tv_sec - this->timer[0].tv_sec);
  int usecs(this->timer[1].tv_usec - this->timer[0].tv_usec);

  if (usecs < 0)
  {
    --secs;
    usecs += 1000000;
  }

  return static_cast<int> (secs * 1000 + usecs / 1000.0 + 0.5);
}
