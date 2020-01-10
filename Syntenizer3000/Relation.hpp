//
//  Relation.hpp
//  Syntenizer3000
//
//  Created by Camous Moslemi on 23/04/2017.
//
//

#ifndef Relation_hpp
#define Relation_hpp

class Gene;
class Strain;
class Group;
class Relation;

#include <stdio.h>
#include <iostream>
#include "Gene.hpp"

class Relation
{
public:
    Gene *gene;
    double score;
    
    Relation(Gene *g, double s);
};

#endif /* Relation_hpp */
