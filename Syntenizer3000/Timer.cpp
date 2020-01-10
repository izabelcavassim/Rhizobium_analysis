//
//  Timer.cpp
//  Syntenizer3000
//
//  Created by Camous Moslemi on 13/09/2017.
//
//

#include "Timer.hpp"

void Timer::Start()
{
    t = clock();
}

double Timer::Elapsed()
{
    return difftime(clock(), t)/(double) CLOCKS_PER_SEC;
}
