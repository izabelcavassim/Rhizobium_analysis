//
//  Group.hpp
//  Syntenizer3000
//
//  Created by Camous Moslemi on 23/04/2017.
//
//

#ifndef Group_hpp
#define Group_hpp

class Gene;
class Strain;
class Group;
class Relation;

#include <stdio.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include "Gene.hpp"
//#include "Needleman-Wunsch.hpp"
#include "Score.hpp"
using namespace std;

class Group
{
public:
    string id;
    int coorthologs;
    //bool discordant;
    double algebraicConnectivity;
    double syntenyScoreSophisticated;
    double syntenyScoreOld;
    double syntenyScoreSimple;
    double syntenyScoreAdjusted;
    double syntenyScoreTest;
    double syntenyScoreFast;
    vector<Gene*> genes;
    //vector<double> scoresSum;
    //vector<vector<double>> scoresMatrix;
    bool HasGene(Gene* g);
    void RemoveGene(Gene* g);
    void InsertGenes(vector<Gene*> newGenes, unordered_map<Gene *, Group *> *grouppool);
    int SharedGenes(vector<Gene*> newGenes);
    double SyntenizeSophisticated();
    double SyntenizeOld();
    double SyntenizeSimple();
    double SyntenizeAdjusted();
    double SyntenizeTest();
    double SyntenizeFast();
    double SyntenizeSimpleAgainstGene(Gene *);
    double SyntenizeAgainstGene(Gene *);
    int CountUniqueStrains();
    int CountParalogs();
    int CountOrthologs();
    void Prune();
    bool HasStrain(Strain *);
    void GenerateSyntenyChart(string);
    void GenerateSyntenyHistogram(string OUTPUT_DIRECTORY );

    vector<Group *> Disambiguate();
    vector<Group *> Decluster(double threshold = 0);
    Group(void);
    //~Group(void);
    
    int orphans;
};

#endif /* Group_hpp */
