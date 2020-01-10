//
//  Settings.hpp
//  Syntenizer3000
//
//  Created by Camous Moslemi on 13/06/2017.
//
//

#ifndef Settings_hpp
#define Settings_hpp

//#define THRESHOLD                   1e-5
#define TYPE_GFF                    "GFF"
#define TYPE_TAB                    "TAB"
#define TYPE_FAA                    "FAA"
//#define SIMILAR_SEQUENCE_FILE       "similarSequences.txt"
//#define RUN_LENGTH                  10
//#define GAP_PENALTY                 1.0
#define SIZE_OF_NEIGHTBOURHOOD      40

#include <string>
//#include <unordered_map>
//#include <vector>
using namespace std;
/*

class Strain;
class Group;
class Relation;
class Gene;
class Settings
{
public:

    //bool SUBGROUPING_ENABLED = false;
    //int SCORING_SYSTEM = 2;
    //int SIZE_OF_NEIGHTBOURHOOD = 4;
    bool scoreSynteny;
    bool generateStrainChart;
    bool generateSyntenyChart;
    string GROUP_FILE;
    string INPUT_DIRECTORY;
    string OUTPUT_DIRECTORY;
    string PLASMID_IDENTIFIER;
    //string TYPE;
    //unordered_map<pair<int64_t, int64_t>, double> similaritypool;
};
*/
struct
{
    bool scoreSynteny = false;
    bool classifyFragments = false;
    bool disambiguateGroups = false;
    bool generatePresenceAbsenceMatrix = false;
    bool generateGAPIT = false;
    bool generateGRM = false;
    bool generateStrainCharts = false;
    bool generateSyntenyCharts = false;
    bool generateFNAGroups = false;
    double declusterGroups = -1;
    string GROUP_FILE;
    string GROUP_COLOUR_FILE;
    string GENE_PAIR_FILE;
    string GFF_DIRECTORY;
    string SEQUENCE_DIRECTORY;
    string OUTPUT_DIRECTORY;
    string PLASMID_IDENTIFIER;
    string CHROMOSOME_IDENTIFIER;
    string FRAGMENT_IDENTIFIER;
    string CONTIG_COLOUR_FILE;
} settings;


#endif /* Settings_hpp */
