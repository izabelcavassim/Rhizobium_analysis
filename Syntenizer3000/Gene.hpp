//
//  Gene.hpp
//  Syntenizer3000
//
//  Created by Camous Moslemi on 23/04/2017.
//
//

#ifndef Gene_hpp
#define Gene_hpp

class Gene;
class Strain;
class Group;
class Relation;
class Contig;

#include <stdio.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include "Strain.hpp"
#include "Group.hpp"
#include "Relation.hpp"
#include "Score.hpp"
#include "Dataset.hpp"
#include "Contig.hpp"

using namespace std;

class Gene
{
public:
    bool masked;
    string id;
    string product;
    string sequence;
    string colour;
    Contig *contig;
    Strain *strain;
    //double coverage;
    double GC3s;
    int start;
    int end;
    bool orientation;
    bool paralog;
    //int copies;
    Group *group;
    Gene* neighbours[SIZE_OF_NEIGHTBOURHOOD];
    Gene *left_neighbour;
    Gene *right_neighbour;
    //vector<Relation*> relations;
    //vector<double> scoresMatrix;
    int length;
    double score;

    unordered_map<Group*, double> groupscores;
    unordered_map<Gene*, double> scores;

    void CompareNeighbours(Gene *gene);
    double Match(Group *);
    double CalculateGC3s();
};

#endif /* Gene_hpp */
