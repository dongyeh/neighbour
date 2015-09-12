#ifndef TIMER_H
#define TIMER_H

#include <boost/date_time/posix_time/posix_time.hpp>

class Timer
{
public:
    Timer(){restart();}
    ~Timer(){}

    void restart(){_start_time=boost::posix_time::microsec_clock::universal_time();}
    double elapsed()
    {
        boost::posix_time::time_duration td=boost::posix_time::microsec_clock::universal_time()-_start_time;
        return static_cast<double>(td.total_milliseconds())/1000;}
private:
    boost::posix_time::ptime _start_time;
};

#endif // TIMER_H
