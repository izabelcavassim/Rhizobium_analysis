//
//  Database.hpp
//  Syntenizer3000
//
//  Created by Camous Moslemi on 01/09/2017.
//
//

#ifndef Database_hpp
#define Database_hpp

class Gene;
class Strain;
class Group;
class Relation;
class Database;

#include <stdio.h>
#include <vector>
#include <unordered_map>
#include <utility>
#include "Gene.hpp"
#include "Group.hpp"
#include "Strain.hpp"

class Database
{
public:
    unordered_map<string, Gene *> genepool;
    unordered_map<string, Strain *> strainpool;
    vector<pair<string, string>> contigcolours;

    vector<Strain*> strains;
    Database();
};

#endif /* Database_hpp */
