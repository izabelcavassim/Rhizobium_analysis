//
//  Score.hpp
//  Syntenizer3000
//
//  Created by Camous Moslemi on 19/09/2017.
//
//

#ifndef Score_hpp
#define Score_hpp

//#define GAP_PENALTY     -1.0f
//#define MATCH_SCORE     1.0f
//#define GAP_SCORE       -0.25f

#include <stdio.h>
#include <map>
#include <unordered_map>
#include <tuple>
#include <utility>
#include "Gene.hpp"
//enum Type{FAST, SIMPLE, NAIVE, SOPHISTICATED, ADJUSTED, OLD};

class Score
{
public:
    static int Fast(Gene *, Gene *);
    static double Simple(Gene *, Gene *);
    static double Naive(Gene *, Gene *);
    static double Sophisticated(Gene *, Gene *);
    static double Adjusted(Gene *, Gene *);
    static double Old(Gene *, Gene *);
private:
};

#endif /* Score_hpp */
