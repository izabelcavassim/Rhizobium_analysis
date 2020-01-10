//
//  Parse.hpp
//  Syntenizer3000
//
//  Created by Camous Moslemi on 01/09/2017.
//
//

#ifndef Parse_hpp
#define Parse_hpp

#include <stdio.h>
#include <string>
#include "Dataset.hpp"
#include "Database.hpp"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <tgmath.h>
#include <iostream>
#include <functional>
#include <unordered_map>
#include "Settings.hpp"
#include "Gene.hpp"
#include "Relation.hpp"
#include "Strain.hpp"
#include "Group.hpp"
//#include "Needleman-Wunsch.hpp"
#include "Timer.hpp"
#include <string>
#include <iostream>
#include <dirent.h>
#include <math.h>
#define GENE_ID_ITEM_NUMBER  1

using namespace std;

class Parse
{
public:
    //static int Groups(string INPUT_FILE, Dataset *dataset,  bool setGeneGroup = true);
    //static int ProteinOrtho(string INPUT_FILE, Dataset *dataset,  bool setGeneGroup = true);
    static int GeneGroups(string INPUT_FILE, Dataset *dataset);
    static int Strains(string INPUT_DIRECTORY, Database *database);
    //static void Scores(string INPUT_DIRECTORY, Dataset *dataset);
    static void Sequences(string INPUT_DIRECTORY, Database *database);
    //static vector<Strain *> Dendro(string INPUT_FILE, Database *database);
    //static void Coverage(string INPUT_DIRECTORY, Database *database);
    static void GroupColours(string INPUT_FILE, Dataset *dataset);
    static void ContigColours(string INPUT_FILE, Database *database);
};

#endif /* Parse_hpp */
