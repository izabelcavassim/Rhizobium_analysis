//
//  Dataset.hpp
//  Syntenizer3000
//
//  Created by Camous Moslemi on 11/09/2017.
//
//

#ifndef Dataset_hpp
#define Dataset_hpp
class Gene;
class Strain;
class Group;
class Relation;
class Dataset;
class Database;

#include <stdio.h>
#include <vector>
#include <unordered_map>
#include <string>
#include <cmath>
#include "Gene.hpp"
#include "Group.hpp"
#include "Strain.hpp"
#include "Database.hpp"
#include "Progress.hpp"
#include "Timer.hpp"
#include "Score.hpp"

using namespace std;

class Dataset
{
public:
    int genecount;
    double syntenyScoreSophisticated;
    double syntenyScoreOld;
    double syntenyScoreSimple;
    double syntenyScoreAdjusted;
    double syntenyScoreTest;
    double syntenyScoreFast;
    Database *database;
    vector<Group*> groups;
    unordered_map<Gene *, Group *> grouppool;
    unordered_map<string, Group *> grouppool2;
    Dataset(Database *);
    void ScoreSynteny(string);
    void ExportFNA(string);
    //void PrintStats();
};

#endif /* Dataset_hpp */
