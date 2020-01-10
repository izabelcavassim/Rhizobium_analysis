//
//  Timer.hpp
//  Syntenizer3000
//
//  Created by Camous Moslemi on 13/09/2017.
//
//

#ifndef Timer_hpp
#define Timer_hpp

#include <tgmath.h>
class Timer
{
public:
    void Start();
    
    double Elapsed();
private:
    clock_t t;
};
#endif /* Timer_hpp */
