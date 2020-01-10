//
//  Progress.hpp
//  Syntenizer3000
//
//  Created by Camous Moslemi on 15/09/2017.
//
//

#ifndef Progress_hpp
#define Progress_hpp

#include <stdio.h>
#include <iterator>

#define STANDARD_PROGRESS_STEPS     10
using namespace std;

class Progress
{
    public:
    Progress(int end, int steps = STANDARD_PROGRESS_STEPS)
    {
        progress = -1;
        this->end = end;
        this->steps = steps;
    }

    void Update(int current)
    {
        current = (current * steps) / end;
        if (current > progress)
        {
            progress = current;
            
            if (progress == 0)
                printf("Progress [");
            printf("%i%% ", progress * (100 / steps) );
            if (progress == steps)
                printf("]\n");
            fflush(stdout);            
        }
    }
    
private:
    int progress;
    int end;
    int steps;
    
};

#endif /* Progress_hpp */
