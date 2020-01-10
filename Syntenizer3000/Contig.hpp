//
//  Contig.hpp
//  Syntenizer3000
//
//  Created by Camous Moslemi on 13/10/2017.
//

#ifndef Contig_hpp
#define Contig_hpp

class Gene;
class Strain;
class Group;
class Relation;

#include <stdio.h>
#include <vector>
#include <string>
#include "Gene.hpp"
#include "Strain.hpp"


using namespace std;

class Contig
{
public:
    string id;
    int length;
    Strain *strain;
    vector<Gene*> genes;
    //vector<bool> orientations;
    //int *coverage;
    Contig();
};

#endif /* Contig_hpp */
