//
//  main.cpp
//  Syntenizer
//
//  Created by Camous.M. on 18/07/16.
//
//
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <tgmath.h>
#include <iostream>
#include <string>
#include <functional>
#include <unordered_map>
#include <set>
#include <iomanip>
#include "Settings.hpp"
#include "Gene.hpp"
#include "Relation.hpp"
#include "Strain.hpp"
#include "Group.hpp"
#include "Parse.hpp"
#include "Database.hpp"
#include "Timer.hpp"

#if defined(__APPLE__) || defined(__linux__) || defined(__unix__)
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#endif

using namespace std;

void GenerateGeneticRelationMatrix(Dataset *dataset, string OUTPUT_FILE, bool header = true)
{
    Timer timer;
    timer.Start();
    printf("Generating genetic relation matrix into:\t%s\n", OUTPUT_FILE.c_str());

    double **presence_absence_matrix = new double*[dataset->groups.size()];
    double **genetic_relation_matrix = new double*[dataset->database->strains.size()];

    for (int i = 0; i < dataset->groups.size(); i++)
        presence_absence_matrix[i] = new double[dataset->database->strains.size()];
    for (int i = 0; i < dataset->database->strains.size(); i++)
        genetic_relation_matrix[i] = new double[dataset->database->strains.size()];

    for (auto strain = dataset->database->strains.begin(); strain != dataset->database->strains.end(); strain++)
    {
        for (auto group = dataset->groups.begin(); group != dataset->groups.end(); group++)
        {
            bool ingroup = false;
            for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
                if ((*gene)->strain == (*strain))
                    ingroup = true;
            //out << (ingroup ? '1' : '0');
            //printf("%i,%i ", distance(dataset->groups.begin(), group), distance(dataset->database->strains.begin(), strain));
            presence_absence_matrix[distance(dataset->groups.begin(), group)][distance(dataset->database->strains.begin(), strain)] = (ingroup ? 2 : 0);
        }
    }
    
    int count = 0;
    double diagonal_mean = 0;

    for (int i = 0; i < dataset->database->strains.size(); i++)
    {
        for (int j = 0; j < dataset->database->strains.size(); j++)
        {
            genetic_relation_matrix[i][j] = 0;
            for (int k = 0; k < dataset->groups.size(); k++)
                genetic_relation_matrix[i][j] += presence_absence_matrix[k][i]*presence_absence_matrix[k][j];
            
            if (i == j)
            {
                count++;
                diagonal_mean += genetic_relation_matrix[i][j];
            }
        }
    }
    
    diagonal_mean /= count;

    ofstream out( OUTPUT_FILE, ifstream::out );
    
    if (header)
        for (auto strain = dataset->database->strains.begin(); strain != dataset->database->strains.end(); strain++)
            out << ';' << (*strain)->id;

    for (int i = 0; i <  dataset->database->strains.size(); i++)
    {
        if (header || (!header && i != 0 ) )
            out << '\n';

        if (header)
            out << dataset->database->strains[i]->id;

        for (int j = 0; j <  dataset->database->strains.size(); j++)
        {
            if (header || (!header && j != 0) )
                out << ';';
            out << genetic_relation_matrix[i][j] / diagonal_mean;
        }
    }

    out.close();

    printf("Generating genetic relation matrix took:\t%i seconds\n\n", (int)timer.Elapsed());
}

void GeneratePresenceAbsenceMatrix(Dataset *dataset, string OUTPUT_FILE, bool header = true)
{
    Timer timer;
    timer.Start();
    
    printf("Generating presence/absence matrix into:\t%s\n", OUTPUT_FILE.c_str());
    
    ofstream out( OUTPUT_FILE, ifstream::out );
    
    if (header)
    {
        out << "Strain";
        for (auto strain = dataset->database->strains.begin(); strain != dataset->database->strains.end(); strain++)
            out << ';' << (*strain)->id;
    }
    for (auto group = dataset->groups.begin(); group != dataset->groups.end(); group++)
    {
        if (header || (!header && group != dataset->groups.begin()) )
            out << '\n';
        if (header)
            out << (*group)->id;

        for (auto strain = dataset->database->strains.begin(); strain != dataset->database->strains.end(); strain++)
        {
            bool ingroup = false;
            for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
                if ((*gene)->strain == (*strain))
                    ingroup = true;
            if (header || (!header && strain != dataset->database->strains.begin()) )
                out << ';';
            out << (ingroup ? '1' : '0');
        }
    }
    
    out << '\n';
    out.close();
    
    printf("Generating presence/absence matrix took:\t%i seconds\n\n", (int)timer.Elapsed());
}

void ExportGroups(string OUTPUT_FILE, Dataset *dataset)
{
    int groups = 1;
    Timer timer;
    timer.Start();
    printf("Exporting gene groups into:\t%s\n", OUTPUT_FILE.c_str());
    ofstream out( OUTPUT_FILE, ifstream::out );
    for (auto group = dataset->groups.begin(); group != dataset->groups.end(); group++, groups++)
    {
        //if (renumerate)
        //    out << "group" + to_string(groups);
        //else
        out << (*group)->id;

        //if (!std::isnan((*group)->algebraicConnectivity))
        //    out << '|' << (*group)->algebraicConnectivity;
        out << ':';
        for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
            out << ' ' << (*gene)->strain->id + '|' + (*gene)->id;
        out << "\n";
    }
    out.close();
    printf("Exporting gene groups took:\t%i seconds\n\n", (int)timer.Elapsed());
}

int AnalyzeParalogs(Dataset *dataset, string OUTPUT_FILE)
{
    int pure_paralog_groups = 0;
    int unambiguous_groups = 0;
    int ambiguous_groups = 0;
    int split_groups = 0;
    int two_paralogs = 0;
    int three_paralogs = 0;
    int four_paralogs = 0;
    int five_paralogs = 0;
    int six_paralogs = 0;
    int sevenplus_paralogs = 0;
    Timer timer;
    timer.Start();
    printf("Exporing paralog info into:\t%s\n", OUTPUT_FILE.c_str());
    ofstream out_analysis( OUTPUT_FILE , ifstream::out );

    for (auto group = dataset->groups.begin(); group != dataset->groups.end(); group++)
    {
        bool has_paralogs = false;
        unordered_map<Strain *, vector<Gene*>*> strains;

        for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
        {
            try
            {
                auto v = strains.at((*gene)->strain);
                v->push_back(*gene);
                has_paralogs = true;
            }
            catch (const out_of_range & e)
            {
                strains[(*gene)->strain] = new vector<Gene*>;
                strains[(*gene)->strain]->push_back(*gene);
            }
        }
        
        if (has_paralogs)
        {
            out_analysis << (*group)->id + ": (";

            vector <Gene*> paralogs;
            vector <Gene*> homologs;
            
            for (auto strain = strains.begin(); strain != strains.end(); strain++)
            {
                if ((*strain).second->size() == 2)
                    two_paralogs++;
                if ((*strain).second->size() == 3)
                    three_paralogs++;
                if ((*strain).second->size() == 4)
                    four_paralogs++;
                if ((*strain).second->size() == 5)
                    five_paralogs++;
                if ((*strain).second->size() == 6)
                    six_paralogs++;
                if ((*strain).second->size() >= 7)
                    sevenplus_paralogs++;

                if ((*strain).second->size() > 1)
                    paralogs.insert(paralogs.end(), (*strain).second->begin(), (*strain).second->end());
                else
                    homologs.insert(homologs.end(), (*strain).second->begin(), (*strain).second->end());
            }
            
            if (homologs.size() == 0)
            {
                pure_paralog_groups++;
                continue;
            }
            
            Strain *prev_strain = (*paralogs.begin())->strain;
            
            bool split = false;
            int total_scoresum = 0;
            int max_scoresum = 0;
            for (auto g1 = paralogs.begin(); g1 != paralogs.end(); g1++)
            {
                double scoresum = 0.0f;
                (*g1)->paralog = true;

                for (auto g2 = homologs.begin(); g2 != homologs.end(); g2++)
                    scoresum += Score::Sophisticated(*g1, *g2);
                
                if (prev_strain != (*g1)->strain)
                {
                    out_analysis << ")\t\t(";
                    prev_strain = (*g1)->strain;
                }
                scoresum /= (double)homologs.size();
                
                for (auto g2 = paralogs.begin(); g2 != paralogs.end(); g2++)
                    if (g1 != g2 && scoresum > 0 && Score::Sophisticated(*g1, *g2) == 0)
                        split = true;
                out_analysis << ' ' + (*g1)->strain->id + '|' + (*g1)->id + " - " + to_string(scoresum);
                
                total_scoresum += scoresum;
                if (scoresum > max_scoresum)
                    max_scoresum = scoresum;
            }
            if (split)
                split_groups++;
            out_analysis << " )\n";

            if (max_scoresum > 0 && total_scoresum == max_scoresum)
                unambiguous_groups++;
            else
                ambiguous_groups++;
        }
    }
    //out_groups.close();
    out_analysis.close();
    printf("Pure paralog groups: %i\n", pure_paralog_groups);
    printf("Split groups: %i\n", split_groups);
    printf("Unambiguous paralog groups: %i\n", unambiguous_groups);
    printf("Ambiguous paralog groups: %i\n", ambiguous_groups);
    printf("two paralogs: %i\n", two_paralogs);
    printf("three paralogs: %i\n", three_paralogs);
    printf("four paralogs: %i\n", four_paralogs);
    printf("five paralogs: %i\n", five_paralogs);
    printf("six paralogs: %i\n", six_paralogs);
    printf("seven or more paralogs: %i\n", sevenplus_paralogs);
    printf("Exporing paralog info took:\t%i seconds\n\n", (int)timer.Elapsed());
    return 0;
}

int SyntenyScoreSubset(Dataset *dataset, string INPUT_FILE)
{
    int counter = 0;
    double scoresumSophisticated = 0.0f;
    //double scoresumOld = 0.0f;
    double scoresumSimple = 0.0f;
    //double scoresumAdjusted = 0.0f;
    //double scoresumFast = 0.0f;
    double syntenyScoreSophisticated;
    //double syntenyScoreOld;
    double syntenyScoreSimple;
    //double syntenyScoreAdjusted;
    //double syntenyScoreFast;
    Timer timer;
    
    ifstream input = ifstream( INPUT_FILE);
    if (!input.is_open())
    {
        printf("\nError: Could not open file %s\n", INPUT_FILE.c_str() );
        exit(1);
    }

    timer.Start();
    string line;
    while (getline(input, line))
    {
        counter++;
        stringstream sstream(line);
        
        if (line.empty() || counter == 1)
            continue;
        
        //if (line[0] != 'g' || line[1] != 'r' || line[2] != 'o' || line[3] != 'u' || line[4] != 'p')
        //    continue;
        
        string item;
        vector<string> items;
        while (getline(sstream, item, ';'))
            items.push_back(item);

        if (items[0] == "")
            continue;

        Group *group = NULL;

        for (auto g = dataset->groups.begin(); g != dataset->groups.end(); g++)
        {
            if ((*g)->id == items[0])
            {
                group = *g;
                break;
            }
        }
        
        if (group == NULL)
        {
            printf("\nError: Group %s not found.\n", items[0].c_str());
            exit(1);
        }
        else
        {
            scoresumSophisticated += group->syntenyScoreSophisticated;
            //scoresumOld += group->syntenyScoreOld;
            scoresumSimple += group->syntenyScoreSimple;
            //scoresumAdjusted += group->syntenyScoreAdjusted;
            //scoresumFast += group->syntenyScoreFast;
        }
    }
    
    syntenyScoreSophisticated = scoresumSophisticated / (double)(counter-1);
    //syntenyScoreOld = scoresumOld / (double)(counter-1);
    syntenyScoreSimple = scoresumSimple / (double)(counter-1);
    //syntenyScoreAdjusted = scoresumAdjusted / (double)(counter-1);
    //syntenyScoreFast = scoresumFast / (double)(counter-1);

    printf("Average Synteny Score Sophisticated:\t[%.2f]\n", syntenyScoreSophisticated);
    //printf("Average Synteny Score Old:\t[%.2f]\n", syntenyScoreOld);
    printf("Average Synteny Score Simple:\t[%.2f]\n", syntenyScoreSimple);
    //printf("Average Synteny Score Adjusted:\t[%.2f]\n", syntenyScoreAdjusted);
    //printf("Average Synteny Score Fast:\t[%.2f]\n", syntenyScoreFast);

    return 0;
}

int SyntenizeGenePairs(Dataset *dataset, string INPUT_FILE, string OUTPUT_FILE)
{
    int counter = 0;
    //double syntenyScoreSophisticated;
    //double syntenyScoreOld;
    //double syntenyScoreSimple;
    //double syntenyScoreAdjusted;
    //double syntenyScoreFast;
    Timer timer;
    
    printf("Scoring synteny for gene pairs from \"%s\":\n", INPUT_FILE.c_str());

    ifstream input = ifstream( INPUT_FILE);
    if (!input.is_open())
    {
        printf("\nError: Could not open file %s\n", INPUT_FILE.c_str() );
        exit(1);
    }
    
    ofstream output( OUTPUT_FILE, ifstream::out );
    
    output << "Strain1;Gene1;Strain2:Gene2;Synteny Simple;Synteny Sophisticated\n";
    timer.Start();
    string line;
    printf("Processing [");
    while (getline(input, line))
    {
        counter++;
        stringstream sstream(line);
        
        if (line.empty() || counter == 1)
            continue;
        
        string item;
        vector<string> items;
        while (getline(sstream, item, ';'))
            items.push_back(item);
        
        if (items[0] == "")
            continue;
        
        Gene *g1;
        try
        {
            g1 = dataset->database->genepool.at(items[0] + '|' + items[1]);
        }
        catch (const out_of_range & e)
        {
            printf("\nError: Gene %s in Strain %s was not found.\n", items[1].c_str(), items[0].c_str());
            exit(1);
        }
        Gene *g2;
        try
        {
            g2 = dataset->database->genepool.at(items[2] + '|' + items[3]);
        }
        catch (const out_of_range & e)
        {
            printf("\nError: Gene %s in Strain %s was not found.\n", items[3].c_str(), items[2].c_str());
            exit(1);
        }

        output << g1->strain->id << ';' << g1->id << g2->strain->id << ';' << g2->id << Score::Simple(g1, g2) << ';' << Score::Sophisticated(g1, g2) << '\n';
        printf("*");
    }
    input.close();
    output.close();
    printf("]\n");

    printf("Scoring synteny for %i gene pairs took:\t%i seconds\n\n", counter, (int)timer.Elapsed());

    return 0;
}

void Disambiguate(string OUTPUT_FILE, Dataset *dataset)
{
    int groups = 0;
    Timer timer;
    timer.Start();
    printf("Disambigating groups into:\t%s\n", OUTPUT_FILE.c_str());

    ofstream out( OUTPUT_FILE, ifstream::out );
    
    //set<Strain *> orphanstrains;
    unordered_map<Strain *, int> orphanstrains;
    unordered_map<string, int> orphancontigs;

    int orphans = 0;
    for (auto group = dataset->groups.begin(); group != dataset->groups.end(); group++)
        if ((*group)->CountParalogs() > 0)
        {
            groups++;
            vector<Group*> splitgroups = (*group)->Disambiguate();
            orphans += (*group)->orphans;
            
            for (auto group2 = splitgroups.begin(); group2 != splitgroups.end(); group2++)
            {
                if ((*group2)->genes.size() == 1)
                {
                    if (orphancontigs.count((*group2)->genes[0]->contig->id))
                        orphancontigs[(*group2)->genes[0]->contig->id]++;
                    else
                        orphancontigs[(*group2)->genes[0]->contig->id] = 1;
                    if (orphanstrains.count((*group2)->genes[0]->strain))
                        orphanstrains[(*group2)->genes[0]->strain]++;
                    else
                        orphanstrains[(*group2)->genes[0]->strain] = 1;

                    continue;
                }
                out << (*group2)->id << ':';
                //out << (*group2)->id;
                //if ((*group2)->algebraicConnectivity != NAN)
                //    out << '|' << (*group2)->algebraicConnectivity;
                //out << ':';
                for (auto gene = (*group2)->genes.begin(); gene != (*group2)->genes.end(); gene++)
                    out << ' ' << (*gene)->strain->id << '|' << (*gene)->id;
                out << '\n';
            }
        }
        else
        {
            out << (*group)->id;
            //if (!std::isnan((*group)->algebraicConnectivity))
            //    out << '|' << (*group)->algebraicConnectivity;
            out << ':';
            for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
                out << ' ' << (*gene)->strain->id << '|' << (*gene)->id;
            out << '\n';
        }
    printf("Orphaned genes %i\n", orphans);
    printf("Orphaned strains %lu\n", orphanstrains.size());

    for (auto orphanedstrain = orphanstrains.begin(); orphanedstrain != orphanstrains.end(); orphanedstrain++)
        printf("%s\t%i\n", orphanedstrain->first->id.c_str(), orphanedstrain->second);
    printf("\n");
    for (auto orphancontig = orphancontigs.begin(); orphancontig != orphancontigs.end(); orphancontig++)
        printf("%s\t%i\n", orphancontig->first.c_str(), orphancontig->second);
    out.close();
    printf("Disambiguating %i groups took:\t%i seconds\n\n", groups, (int)timer.Elapsed());
}

void Decluster(string OUTPUT_FILE, Dataset *dataset, double threshold)
{
    int groups = 0;
    Timer timer;
    timer.Start();
    printf("Declustering groups into:\t%s\n", OUTPUT_FILE.c_str());

    ofstream out( OUTPUT_FILE, ifstream::out );

    //set<Strain *> orphanstrains;
    unordered_map<Strain *, int> orphanstrains;
    unordered_map<string, int> orphancontigs;

    int orphans = 0;
    int erased = 0;

    Progress p(dataset->groups.size()-1);

    for (auto group = dataset->groups.begin(); group != dataset->groups.end(); group++)
    {
        erased += (*group)->genes.size();

        groups++;
        vector<Group*> splitgroups = (*group)->Decluster(threshold);

        int splitgroupgenes = 0;
        for (auto group2 = splitgroups.begin(); group2 != splitgroups.end(); group2++)
        {
            splitgroupgenes += (*group2)->genes.size();
            if ((*group2)->genes.size() == 1)
            {
                if (orphancontigs.count((*group2)->genes[0]->contig->id))
                    orphancontigs[(*group2)->genes[0]->contig->id]++;
                else
                    orphancontigs[(*group2)->genes[0]->contig->id] = 1;
                if (orphanstrains.count((*group2)->genes[0]->strain))
                    orphanstrains[(*group2)->genes[0]->strain]++;
                else
                    orphanstrains[(*group2)->genes[0]->strain] = 1;

                continue;
            }
            if (splitgroups.size() == 1)
                out << (*group)->id << ':';
            else
                out << (*group2)->id << ':';

            for (auto gene = (*group2)->genes.begin(); gene != (*group2)->genes.end(); gene++)
                out << ' ' << (*gene)->strain->id << '|' << (*gene)->id;
            out << '\n';
        }
        erased -= splitgroupgenes;
        p.Update(distance(dataset->groups.begin(), group));
    }
    printf("Erased genes %i\n", erased);
    printf("Orphaned genes %i\n", orphans);
    printf("Orphaned strains %lu\n", orphanstrains.size());

    for (auto orphanedstrain = orphanstrains.begin(); orphanedstrain != orphanstrains.end(); orphanedstrain++)
        printf("%s\t%i\n", orphanedstrain->first->id.c_str(), orphanedstrain->second);
    printf("\n");
    for (auto orphancontig = orphancontigs.begin(); orphancontig != orphancontigs.end(); orphancontig++)
        printf("%s\t%i\n", orphancontig->first.c_str(), orphancontig->second);
    out.close();
    printf("Declustering %i groups took:\t%i seconds\n\n", groups, (int)timer.Elapsed());
}

double CalculateAverageParalogSynteny(Dataset *dataset)
{
    double counterParalogs = 0;
    double scoresumParalogs = 0;
    double counterOrthologs = 0;
    double scoresumOrthologs = 0;
    //double scoresum = 0;
    int ps = 0;
    int os = 0;
    Timer timer;
    timer.Start();
    printf("Calculating average synteny between paralog pairs:\n");
    
    Progress p(dataset->groups.size()-1);

    for (auto group = dataset->groups.begin(); group != dataset->groups.end(); group++)
    {
        if ((*group)->CountParalogs() > 0)
        {
            vector<Gene *> orthologs;
            vector<Gene *> paralogs;
            set<Strain *> strains;
            set<Strain *> paralogstrains;
            
            for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
                if (strains.count((*gene)->strain))
                    paralogstrains.insert((*gene)->strain);
                else
                    strains.insert((*gene)->strain);
            
            for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
                if (paralogstrains.count((*gene)->strain))
                    paralogs.push_back(*gene);
                else
                    orthologs.push_back(*gene);
            
            ps += paralogs.size();
            os += orthologs.size();

            for (auto g1 = paralogs.begin(); g1 != paralogs.end(); g1++)
                for (auto g2 = g1; g2 != paralogs.end(); g2++)
                    if (g1 != g2 && (*g1)->strain == (*g2)->strain)
                    {
                        counterParalogs++;
                        scoresumParalogs += Score::Sophisticated(*g1, *g2);
                        //scoresum += Score::Simple(*g1, *g2);
                    }
            
            for (auto g1 = orthologs.begin(); g1 != orthologs.end(); g1++)
                for (auto g2 = g1; g2 != orthologs.end(); g2++)
                    if (g1 != g2)
                    {
                        counterOrthologs++;
                        scoresumOrthologs += Score::Sophisticated(*g1, *g2);
                        //scoresum += Score::Simple(*g1, *g2);
                    }
        }
        p.Update(distance(dataset->groups.begin(), group));
    }
    printf("Calculating a synteny average of %.2f for %i paralogs and %.2f for %i orthologs took:\t%i seconds\n\n", (scoresumParalogs /counterParalogs), (int)ps,(scoresumOrthologs / counterOrthologs), (int)os,(int)timer.Elapsed());
    return (scoresumParalogs / counterParalogs);
}

void ExportStrainsChartData(string OUTPUT_DIRECTORY, Database *db)
{
    Timer timer;
    timer.Start();
    printf("Exporting gene chart data into:\t%s\n", OUTPUT_DIRECTORY.c_str());
    
    ofstream out( OUTPUT_DIRECTORY + "strains.csv", ifstream::out );
    out << "ID\n";
    
    printf("Processing [");
    for (auto strain = db->strains.begin(); strain != db->strains.end(); strain++)
    {
        out << (*strain)->id << '\n';
        (*strain)->ExportChartData(OUTPUT_DIRECTORY, db);
        printf("*");
    }
    printf("]\n");
    out.close();
    printf("Exporting gene chart data took:\t%i seconds\n\n", (int)timer.Elapsed());
}

void ExportSyntenyChartData(string OUTPUT_DIRECTORY, Dataset *dataset)
{
    Timer timer;
    timer.Start();
    printf("Generating synteny chart data for gene groups into:\t%s\n", OUTPUT_DIRECTORY.c_str());

    Progress p(dataset->groups.size()-1);
    for (auto group = dataset->groups.begin(); group != dataset->groups.end(); group++)
    {
        (*group)->GenerateSyntenyChart(settings.OUTPUT_DIRECTORY);
        p.Update(distance(dataset->groups.begin(), group));
    }
    printf("Generating synteny chart data for %lu gene groups took:\t%i seconds\n\n", dataset->groups.size(), (int)timer.Elapsed());
}

void ExportStatistics(Dataset *dataset, string OUTPUT_FILE)
{
    Timer timer;
    timer.Start();
    printf("Exporting strain statistics into:\t%s\n", OUTPUT_FILE.c_str());
    
    ofstream out( OUTPUT_FILE, ifstream::out );

    out << "Strain ID;Contigs;bp;ubp;GC;CDS;Orthologs;Paralogs;Uniques;tRNA;rRNA;tmRNA;0-10k;10-250k;250k-500k;500k-750k;750k-1000k;1000-5000k;>5000k\n";

    printf("Processing [");
    //int counter = 0;
    for (auto strain = dataset->database->strains.begin(); strain != dataset->database->strains.end(); strain++)
    {
        //int CDS = 0;

        uint64_t orthologs = 0;
        uint64_t paralogs = 0;
        uint64_t uniques = 0;
        uint64_t contigsize10k = 0;
        uint64_t contigsize250k = 0;
        uint64_t contigsize500k = 0;
        uint64_t contigsize750k = 0;
        uint64_t contigsize1000k = 0;
        uint64_t contigsize5000k = 0;
        uint64_t contigsize5000kplus = 0;

        for (auto contig = (*strain)->contigs.begin(); contig != (*strain)->contigs.end(); contig++)
        {
            if ((*contig)->length < 0)
            {
                printf("\nError: Strain %s had contig %s with unknown size!\n", (*strain)->id.c_str(), (*contig)->id.c_str());
                exit(1);
            }
                
            if ((*contig)->length < 10000)
                contigsize10k+= (*contig)->length;
                //contigsize10k++;
            else if ((*contig)->length < 250000)
                contigsize250k+= (*contig)->length;
                //contigsize250k++;
            else if ((*contig)->length < 500000)
                contigsize500k+= (*contig)->length;
                //contigsize500k++;
            else if ((*contig)->length < 750000)
                contigsize750k+= (*contig)->length;
               //contigsize750k++;
            else if ((*contig)->length < 1000000)
                contigsize1000k+= (*contig)->length;
                //contigsize1000k++;
            else if ((*contig)->length < 5000000)
                contigsize5000k+= (*contig)->length;
                //contigsize5000k++;
            else
                contigsize5000kplus+= (*contig)->length;
                //contigsize5000kplus++;
            for (auto g1 = (*contig)->genes.begin(); g1 != (*contig)->genes.end(); g1++)
            {
                if ((*g1)->group == NULL)
                    uniques++;
                else
                {
                    int count = 0;
                    for (auto g2 = (*g1)->group->genes.begin(); g2 != (*g1)->group->genes.end(); g2++)
                        if ((*g2)->strain == *strain)
                            count++;

                    if (count == 1)
                        orthologs++;
                    else
                        paralogs++;
                }
            }
        }
        
        if ((*strain)->CDS != orthologs + paralogs + uniques)
            printf("Warning: CDS (%i) != orthologs (%i) + paralogs (%i) + uniques (%i)\n", (*strain)->CDS,orthologs,paralogs, uniques);

        cout.setf(ios::fixed);

        out << (*strain)->id << ';';
        out << (*strain)->contigs.size() << ';';
        out << (*strain)->bp << ';';
        out << (*strain)->ubp << ';';
        out << setprecision(4) << (*strain)->GC << ';';
        out << (*strain)->CDS << ';';
        out << orthologs << ';';
        out << paralogs << ';';
        out << uniques << ';';
        out << (*strain)->tRNA << ';';
        out << (*strain)->rRNA << ';';
        out << (*strain)->tmRNA << ';';

        if ((*strain)->bp > 0)
        {
            out << setprecision(3) << (contigsize10k / (double)(*strain)->bp)*100.0f << "%;";
            out << setprecision(3) << (contigsize250k / (double)(*strain)->bp)*100.0f << "%;";
            out << setprecision(3) << (contigsize500k / (double)(*strain)->bp)*100.0f << "%;";
            out << setprecision(3) << (contigsize750k / (double)(*strain)->bp)*100.0f << "%;";
            out << setprecision(3) << (contigsize1000k / (double)(*strain)->bp)*100.0f << "%;";
            out << setprecision(3) << (contigsize5000k / (double)(*strain)->bp)*100.0f << "%;";
            out << setprecision(3) << (contigsize5000kplus / (double)(*strain)->bp)*100.0f << "%\n";
        }
        else
            out << "NA;NA;NA;NA;NA;NA;NA\n";

        printf("*");
    }
    printf("]\n");

    out.close();
    printf("Exporting statistics for %lu strains and %lu gene groups took:\t%i seconds\n\n", dataset->database->strains.size(), dataset->groups.size(), (int)timer.Elapsed());
}

void GeneratePresenceAbsenceGAPIT(Dataset *dataset, string OUTPUT_FILE)
{
    Timer timer;
    timer.Start();
    printf("Generating gene presence/absence GAPIT format data into:\t%s\n", OUTPUT_FILE.c_str());

    ofstream out( OUTPUT_FILE, ifstream::out );
    out << "rs;Allele;chrom;Position;Strand;assembly;center;portLSID;assayLSID;panel;Qccode";
    
    for (auto strain = dataset->database->strains.begin(); strain != dataset->database->strains.end(); strain++)
        out << ';' << (*strain)->id;
    out << '\n';
    
    int position = 100;
    vector<int> positions;
    for (int i = 0; i <= dataset->database->strains.size(); i++)
        positions.push_back(1000);

    for (auto group = dataset->groups.begin(); group != dataset->groups.end(); group++)
    {
        set<Strain *> strains;
         for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
             strains.insert((*gene)->strain);

        out << (*group)->id << ';';
        out << "A/G" << ';';
        out << '1' << ';';
        out << position << ';';
        out << "NA;NA;NA;NA;NA;NA;NA";
        for (auto strain = dataset->database->strains.begin(); strain != dataset->database->strains.end(); strain++)
            out << ';' << ((*group)->HasStrain(*strain) ? "AA" : "GG");
        out << '\n';
        positions[strains.size()] += 1000;
        position += 100;
    }

    out.close();
    printf("Generating gene presence/absence GAPIT format data took:\t%i seconds\n\n", (int)timer.Elapsed());
}

void ExportGeneScores(Dataset *dataset, string OUTPUT_FILE)
{
    Timer timer;
    timer.Start();
    printf("Exporting gene scores into:\t%s\n", OUTPUT_FILE.c_str());

    ofstream out( OUTPUT_FILE, ifstream::out );

    out << "Gene1;Gene2;Sophisticated\n";

    for (auto group = dataset->groups.begin(); group != dataset->groups.end(); group++)
    {
        for (auto gene1 = (*group)->genes.begin(); gene1 != (*group)->genes.end(); gene1++)
        {
            for (auto score = (*gene1)->scores.begin();  score != (*gene1)->scores.end(); score++ )
            {
                Gene *gene2 = (*score).first;
                out << (*gene1)->id << ';';
                out << gene2->id << ';';
                out << (*score).second << '\n';
            }
        }
    }
    out.close();
    printf("Exporting gene scores took:\t%i seconds\n\n", (int)timer.Elapsed());
}
/*
void ClassifyFragments(string OUTPUT_FILE, Database *db)
{
    Timer timer;
    timer.Start();
    printf("Classifying contig fragments into:\t%s\n", OUTPUT_FILE.c_str());
    
    ofstream out( OUTPUT_FILE, ifstream::out );
    out << "Strain;Contig;Genes;Chromosomal Percentage;Plasmid Percentage;Chromosomal Synteny;Plasmid Synteny;Coverage;Syneny;Magnitude;Confidence;Classification;Metadata\n";

    printf("Processing [");
    for (auto strain = db->strains.begin(); strain != db->strains.end(); strain++)
    {
        for (auto contig = (*strain)->contigs.begin(); contig != (*strain)->contigs.end(); contig++)
        {
            if ((*contig)->genes.size() < 10)
                continue;
            out << (*strain)->id << ';' << (*contig)->id << ';' << (*contig)->genes.size() <<';';
            double chromosomalGene = 0;
            double plasmidGene = 0;
            double totalGenes = 0;

            double chromosomalSynteny = 0;
            double plasmidSynteny = 0;
            double plasmidTypeSynteny[40] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            double plasmidTypeCount[40] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            double plasmidTypeGenes = 0;

            for (auto g1 = (*contig)->genes.begin(); g1 != (*contig)->genes.end(); g1++)
            {
                if ((*g1)->group == NULL || (*g1)->masked)
                    continue;

                for (auto g2 = (*g1)->group->genes.begin(); g2 != (*g1)->group->genes.end(); g2++)
                {
                    totalGenes++;
                    if ((*g2)->contig->id.find(settings.CHROMOSOME_IDENTIFIER) != string::npos)
                    {
                        chromosomalGene++;
                        //chromosomalSynteny += Score::Sophisticated(*g1, *g2);
                        chromosomalSynteny += Score::Simple(*g1, *g2);
                    }
                    if ((*g2)->contig->id.find(settings.PLASMID_IDENTIFIER) != string::npos)
                    {
                        //double score = Score::Sophisticated(*g1, *g2);
                        double score = Score::Simple(*g1, *g2);
                        int type;
                        size_t loc;
                        string id = (*g2)->contig->id;
                        while ((loc = id.find(settings.PLASMID_IDENTIFIER)) != string::npos)
                        {
                            type = stoi(id.substr(loc+2, 2));
                            plasmidTypeSynteny[type] += score;
                            plasmidTypeCount[type]++;
                            plasmidTypeGenes++;
                            id = id.substr(loc+4);
                        }

                        plasmidSynteny += score;
                        plasmidGene++;
                    }
                }
            }

            double A = 0;
            double B = 0;

            if (totalGenes > 0)
               A = (chromosomalGene / totalGenes) - (plasmidGene / totalGenes);
            if (chromosomalGene > 0)
                B += (chromosomalSynteny / chromosomalGene)/SIZE_OF_NEIGHTBOURHOOD;
            if (plasmidGene > 0)
                B -= (plasmidSynteny / plasmidGene)/SIZE_OF_NEIGHTBOURHOOD;

            double length = 0;

            if ((A <= 0 && B <= 0) || (A <= -0.2 && B <= 0.1) || (A <= 0.1 && B <= -0.2))
                length = -sqrt(A*A+B*B);
            if ((A >= 0 && B >= 0) || (A >= -0.1 && B >= 0.2) || (A >= 0.2 && B >= -0.1))
                length = sqrt(A*A+B*B);

            string confidence = "NA";

            if (abs(length) >= 0.2)
                confidence = "LOW";
            if (abs(length) >= 0.4)
                confidence = "MEDIUM";
            if (abs(length) >= 0.6)
                confidence = "HIGH";
            if (abs(length) >= 0.8)
                confidence = "IDEAL";
            
            string classification = "Ambiguous";
            string metadata;

            if (length >= 0.2)
                classification = "Chromosome";

            if (length <= -0.2)
            {
                //double bestSynteny = -1;
                double bestMatchScore = -1;
                int bestPlasmid = -1;
                for (int i = 0; i < 40; i++)
                {
                    if (plasmidTypeCount[i] > 0)
                    {
                        plasmidTypeSynteny[i] /= plasmidTypeCount[i];
                        double score = (plasmidTypeSynteny[i] / SIZE_OF_NEIGHTBOURHOOD) * 0.75 + (plasmidTypeCount[i] / plasmidTypeGenes) * 0.25;
                        if (score > bestMatchScore)
                        {
                            bestMatchScore = score;
                            bestPlasmid = i;

                            if (score >= 0.0)
                                confidence = "LOW";
                            if (score >= 0.3)
                                confidence = "MEDIUM";
                            if (score >= 0.5)
                                confidence = "HIGH";
                            if (score >= 0.7)
                                confidence = "IDEAL";
                        }

                        metadata += settings.PLASMID_IDENTIFIER + to_string(i) + "|" + to_string(int(plasmidTypeCount[i])) + "|" + to_string(int(plasmidTypeSynteny[i])) + "|" + to_string(score) +" ";
                    }
                }

                if (bestPlasmid < 10)
                    classification = settings.PLASMID_IDENTIFIER + "0" + to_string(bestPlasmid);
                else
                    classification = settings.PLASMID_IDENTIFIER + to_string(bestPlasmid);
            }

            out << setprecision(2) << (chromosomalGene / totalGenes) << ';';
            out << setprecision(2) << (plasmidGene / totalGenes);
            out << ';' << setprecision(2) << (chromosomalGene > 0 ? (chromosomalSynteny / chromosomalGene) : 0.0f) << ';';
            out << setprecision(2) << (plasmidGene > 0 ? (plasmidSynteny / plasmidGene) : 0.0f) << ';';
            out << setprecision(2) << A << ';';
            out << setprecision(2) << B << ';';
            out << length << ';';
            out << confidence << ';';
            out << classification << ';';
            out << metadata << '\n';
        }
        printf("*");
    }
    printf("]\n");
    out.close();
    printf("Classifying contig fragments took:\t%i seconds\n\n", (int)timer.Elapsed());
}
*/

void ClassifyFragments(string OUTPUT_FILE, Database *db)
{
    Timer timer;
    timer.Start();
    printf("Classifying contig fragments into:\t%s\n", OUTPUT_FILE.c_str());

    ofstream out( OUTPUT_FILE, ifstream::out );
    out << "Strain;Contig;Genes;Chromosomal Percentage;Plasmid Percentage;Chromosomal Synteny;Plasmid Synteny;Coverage;Syneny;Magnitude;Confidence;Classification;Metadata\n";

    printf("Processing [");
    for (auto strain = db->strains.begin(); strain != db->strains.end(); strain++)
    {
        for (auto contig = (*strain)->contigs.begin(); contig != (*strain)->contigs.end(); contig++)
        {
            if (/*(*contig)->id_full.find("fragment") != string::npos ||*/ (*contig)->genes.size() < 10)
                continue;
            out << (*strain)->id << ';' << (*contig)->id << ';' << (*contig)->genes.size() <<';';
            double chromosomalGene = 0;
            double plasmidGene = 0;
            double totalGenes = 0;

            double chromosomalSynteny = 0;
            double plasmidSynteny = 0;
            double plasmidTypeSynteny[100] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            double plasmidTypeCount[100] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            double plasmidTypeGenes = 0;

            for (auto g1 = (*contig)->genes.begin(); g1 != (*contig)->genes.end(); g1++)
            {
                if ((*g1)->group == NULL || (*g1)->masked)
                    continue;

                for (auto g2 = (*g1)->group->genes.begin(); g2 != (*g1)->group->genes.end(); g2++)
                {
                    totalGenes++;
                    if ((*g2)->contig->id.find(settings.CHROMOSOME_IDENTIFIER) != string::npos)
                    {
                        chromosomalGene++;
                        //chromosomalSynteny += Score::Sophisticated(*g1, *g2);
                        chromosomalSynteny += Score::Simple(*g1, *g2);
                    }
                    if ((*g2)->contig->id.find(settings.PLASMID_IDENTIFIER) != string::npos)
                    {
                        //double score = Score::Sophisticated(*g1, *g2);
                        double score = Score::Simple(*g1, *g2);
                        int type;
                        size_t loc;
                        string id = (*g2)->contig->id;
                        while ((loc = id.find(settings.PLASMID_IDENTIFIER)) != string::npos)
                        {
                            try
                            {
                                type = stoi(id.substr(loc+2, 2));
                            }
                            catch (exception e)
                            {
                                printf("\nError: Plasmid ID \"%s\" does not confirm to expected format of an identifier character pattern followed by two digits.\n", id.c_str());
                                exit(1);

                            }
                            plasmidTypeSynteny[type] += score;
                            plasmidTypeCount[type]++;
                            plasmidTypeGenes++;
                            id = id.substr(loc+4);
                        }

                        plasmidSynteny += score;
                        plasmidGene++;
                    }
                }
            }

            double A = 0;
            double B = 0;

            if (totalGenes > 0)
                A = (chromosomalGene / totalGenes) - (plasmidGene / totalGenes);
            if (chromosomalGene > 0)
                B += (chromosomalSynteny / chromosomalGene)/SIZE_OF_NEIGHTBOURHOOD;
            if (plasmidGene > 0)
                B -= (plasmidSynteny / plasmidGene)/SIZE_OF_NEIGHTBOURHOOD;

            double magnitude = A + B;

            string confidence = "NA";

            if (abs(magnitude) >= 0.4)
                confidence = "LOW";
            if (abs(magnitude) >= 0.8)
                confidence = "MEDIUM";
            if (abs(magnitude) >= 1.2)
                confidence = "HIGH";
            if (abs(magnitude) >= 1.6)
                confidence = "IDEAL";

            string classification = "Ambiguous";
            string metadata;

            if (magnitude >= 0.4)
                classification = "Chromosome";

            if (magnitude <= -0.4)
            {
                //double bestSynteny = -1;
                double bestMatchScore = -1;
                int bestPlasmid = -1;
                for (int i = 0; i < 100; i++)
                {
                    if (plasmidTypeCount[i] > 0)
                    {
                        plasmidTypeSynteny[i] /= plasmidTypeCount[i];
                        double score = (plasmidTypeSynteny[i] / SIZE_OF_NEIGHTBOURHOOD) * 0.75 + (plasmidTypeCount[i] / plasmidTypeGenes) * 0.25;
                        if (score > bestMatchScore)
                        {
                            bestMatchScore = score;
                            bestPlasmid = i;

                            if (score >= 0.0)
                                confidence = "LOW";
                            if (score >= 0.3)
                                confidence = "MEDIUM";
                            if (score >= 0.5)
                                confidence = "HIGH";
                            if (score >= 0.7)
                                confidence = "IDEAL";
                        }

                        metadata += settings.PLASMID_IDENTIFIER + to_string(i) + "|" + to_string(int(plasmidTypeCount[i])) + "|" + to_string(int(plasmidTypeSynteny[i])) + "|" + to_string(score) +" ";
                    }
                }

                if (bestPlasmid < 10)
                    classification = settings.PLASMID_IDENTIFIER + "0" + to_string(bestPlasmid);
                else
                    classification = settings.PLASMID_IDENTIFIER + to_string(bestPlasmid);
            }

            out << setprecision(2) << (chromosomalGene / totalGenes) << ';';
            out << setprecision(2) << (plasmidGene / totalGenes);
            out << ';' << setprecision(2) << (chromosomalGene > 0 ? (chromosomalSynteny / chromosomalGene) / SIZE_OF_NEIGHTBOURHOOD : 0.0f) << ';';
            out << setprecision(2) << (plasmidGene > 0 ? (plasmidSynteny / plasmidGene) / SIZE_OF_NEIGHTBOURHOOD : 0.0f) << ';';
            out << setprecision(2) << A << ';';
            out << setprecision(2) << B << ';';
            out << magnitude << ';';
            out << confidence << ';';
            out << classification << ';';
            out << metadata << '\n';
        }
        printf("*");
    }
    printf("]\n");
    out.close();
    printf("Classifying contig fragments took:\t%i seconds\n\n", (int)timer.Elapsed());
}
/*
int main(int argc, const char * argv[])
{
    #if defined(__APPLE__) || defined(__linux__) || defined(__unix__)
    int filedesc = open("render_charts.sh", O_RDWR|O_CREAT|O_TRUNC, 755);
    if(filedesc >= 0)
    {
        string script = "Rscript " + ;
        write(filedesc, script.c_str(), script.length());
    }

    #endif
    
    #if defined(_WIN32) || defined(WIN32)
    
    #endif

    return(0);
}
*/

int main(int argc, const char * argv[])
{
    Database *db = new Database();
    Dataset *dataset = new Dataset(db);

    setbuf(stdout, NULL);

    printf("\nSyntenizer 3000 version 1.0 by Camous Moslemi.\n");

    if (argc <= 1)
    {
        printf("Usage:\tsyntenizer --gffdir=/data/ --outdir=/out/ -genegroups=/myproject.poff\n");
        printf("\nOPTIONS\n");
        printf("\t--gffdir=dir\n\t\tSet the directory where strain gff files are located.\n\n");
        printf("\t--outdir=dir\n\t\tSet the output directory for all generated files.\n\n");
        printf("\t--genegroups=file\n\t\tSet the file specifying gene groups.\n\n");
        printf("\t--ffndir=dir\n\t\tSet the directory where gene sequence ffn files are located.\n\n");
        printf("\t--generatefnagroups\n\t\tGenerate gene group fna files.\n\n");
        printf("\t--syntenizegroups\n\t\tGenerate synteny score for gene groups.\n\n");
        printf("\t--syntenizegenepairs=file\n\t\tSet the file specifying gene pairs to be syntenized.\n\n");
        printf("\t--disambiguategroups\n\t\tGenerate gene groups devoid of paralogs.\n\n");
        printf("\t--declustergroups\n\t\tSplit gene groups into subgroups if a collection of genes share a synteny score of 0.\n\n");
        printf("\t--declustergroups=threshold\n\t\tSplit gene groups into subgroups if a collection of genes share a synteny score below provided threshold.\n\n");
        //printf("\t--presenceabsencematrix\n\t\tGenerate presence/absence matrix based on gene groups.\n\n");
        //printf("\t--gapitpresenceabsencematrix\n\t\tGenerate GAPIT formatted presence/absence matrix based on gene groups.\n\n");
        //printf("\t--generelationmatrix\n\t\tGenerate gene relation matrix based on gene groups.\n\n");
        printf("\t--generatestraincharts\n\t\tGenerate data files for strain chart rendering.\n\n");
        printf("\t--generatesyntenycharts\n\t\tGenerate data files for group synteny chart rendering.\n\n");
        printf("\t--contigcolours=file\n\t\tSet the file specifying the contig colours.\n\n");
        printf("\t--chromosomeidentifier=pattern\n\t\tSet the contig name character pattern that identifies a chromosomal contig.\n\n");
        printf("\t--plasmididentifier=pattern\n\t\tSet the contig name character pattern that identifies a plasmid contig.\n\n");
        printf("\t--fragmentidentifier=pattern\n\t\tSet the contig name character pattern that identifies an unclassified contig.\n\n");
        printf("\t--genegrouphighlightcolours=file\n\t\tSet the file specifying gene group highlight colours on strain charts.\n\n");
        printf("\t--classifyfragments\n\t\tGenerate fragment classifications based on synteny and gene group information.\n\n");

        exit(0);
    }
    else
    {
        for (int arg = 1; arg < argc; arg++)
        {
            string argument = argv[arg];

            //printf("Argument %i: %s\n", arg, argument.c_str());

            string::size_type loc;
            if ((loc = argument.find("--gffdir=")) != string::npos)
                settings.GFF_DIRECTORY = argument.substr(argument.find("=")+1);
            if ((loc = argument.find("--outdir=")) != string::npos)
                settings.OUTPUT_DIRECTORY = argument.substr(argument.find("=")+1);
            if ((loc = argument.find("--genegroups=")) != string::npos)
                settings.GROUP_FILE = argument.substr(argument.find("=")+1);
            if ((loc = argument.find("--contigcolours=")) != string::npos)
                settings.CONTIG_COLOUR_FILE = argument.substr(argument.find("=")+1);
            if ((loc = argument.find("--ffndir=")) != string::npos)
                settings.SEQUENCE_DIRECTORY = argument.substr(argument.find("=")+1);
            if ((loc = argument.find("--syntenizegroups")) != string::npos)
                settings.scoreSynteny = true;
            if ((loc = argument.find("--disambiguategroups")) != string::npos)
                settings.disambiguateGroups = true;
            if ((loc = argument.find("--generatefnagroups")) != string::npos)
                settings.generateFNAGroups = true;
            if ((loc = argument.find("--generatestraincharts")) != string::npos)
                settings.generateStrainCharts = true;
            if ((loc = argument.find("--generatesyntenycharts")) != string::npos)
                settings.generateSyntenyCharts = true;
            if ((loc = argument.find("--declustergroups")) != string::npos)
            {
                if (settings.disambiguateGroups)
                {
                    printf("\nError: The declustergroups and disambiguategroups options should not be toggled in the same session. Please toggle each in separate sessions.\n\n");
                    exit(1);
                }
                if ((loc = argument.find("--declustergroups=")) != string::npos)
                    settings.declusterGroups = atof(argument.substr(argument.find("=")+1).c_str());
                else
                    settings.declusterGroups = 0.0f;
            }
            //if ((loc = argument.find("--presenceabsencematrix")) != string::npos)
                settings.generatePresenceAbsenceMatrix = true;
            //if ((loc = argument.find("--gapitpresenceabsencematrix")) != string::npos)
                settings.generateGAPIT = true;
            //if ((loc = argument.find("--generelationmatrix")) != string::npos)
                settings.generateGRM = true;
            if ((loc = argument.find("--chromosomeidentifier=")) != string::npos)
                settings.CHROMOSOME_IDENTIFIER = argument.substr(argument.find("=")+1);
            if ((loc = argument.find("--plasmididentifier=")) != string::npos)
                settings.PLASMID_IDENTIFIER = argument.substr(argument.find("=")+1);
            if ((loc = argument.find("--fragmentidentifier=")) != string::npos)
                settings.FRAGMENT_IDENTIFIER = argument.substr(argument.find("=")+1);
            if ((loc = argument.find("--syntenizegenepairs=")) != string::npos)
                settings.GENE_PAIR_FILE = argument.substr(argument.find("=")+1);
            if ((loc = argument.find("--genegrouphighlightcolours=")) != string::npos)
                settings.GROUP_COLOUR_FILE = argument.substr(argument.find("=")+1);
            if ((loc = argument.find("--classifyfragments")) != string::npos)
                settings.classifyFragments = true;
        }
    }

#if defined(_WIN32) || defined(WIN32)
    char delimiter = '\\';
#else
    char delimiter = '/';
#endif
    if (!settings.GFF_DIRECTORY.empty() && settings.GFF_DIRECTORY[settings.GFF_DIRECTORY.size()-1] != delimiter)
        settings.GFF_DIRECTORY += delimiter;
    if (!settings.OUTPUT_DIRECTORY.empty() && settings.OUTPUT_DIRECTORY[settings.OUTPUT_DIRECTORY.size()-1] != delimiter)
        settings.OUTPUT_DIRECTORY += delimiter;
    if (!settings.SEQUENCE_DIRECTORY.empty() && settings.SEQUENCE_DIRECTORY[settings.SEQUENCE_DIRECTORY.size()-1] != delimiter)
        settings.SEQUENCE_DIRECTORY += delimiter;

    Timer timer;
    timer.Start();

    if (!settings.GFF_DIRECTORY.empty())
    {
        //printf("Attempting to parse strains:\n");
        Parse::Strains(settings.GFF_DIRECTORY, db);
        //printf("Done parsing strains.\n\n");

        if (!settings.SEQUENCE_DIRECTORY.empty())
        {
            //printf("Attempting to parse sequences:\n");
            Parse::Sequences(settings.SEQUENCE_DIRECTORY, db);
            //printf("Done parsing sequences.\n\n");
        }

        if (!settings.GROUP_FILE.empty())
        {
            //printf("Attempting to parse ProteinOrtho groups:\n");
            Parse::GeneGroups(settings.GROUP_FILE, dataset);
            ExportGroups(settings.OUTPUT_DIRECTORY + "gene_groups.csv", dataset);
            //printf("Done parsing ProteinOrtho groups.\n\n");
            
            if (!settings.CONTIG_COLOUR_FILE.empty())
            {
                //printf("Attempting to parse contig colours:\n");
                Parse::ContigColours(settings.CONTIG_COLOUR_FILE, db);
                //printf("Done parsingcontig colours.\n\n");
            }
            if (!settings.GROUP_COLOUR_FILE.empty())
            {
                //printf("Attempting to parse gene group highlight colours:\n");
                Parse::GroupColours(settings.GROUP_COLOUR_FILE, dataset);
                //printf("Done gene group highlight colours.\n\n");
            }
            if (settings.generatePresenceAbsenceMatrix)
            {
                //printf("Attempting to generate Presence Absence Matrix:\n");
                GeneratePresenceAbsenceMatrix(dataset, settings.OUTPUT_DIRECTORY + "presence_absence_matrix.csv", true);
                //printf("Done generating Presence Absence Matrix.\n\n");
            }
            if (settings.generateGAPIT)
            {
                //printf("Attempting to generate GAPIT Presence Absence Matrix:\n");
                GeneratePresenceAbsenceGAPIT(dataset, settings.OUTPUT_DIRECTORY + "presence_absence_gapit.csv");
                //printf("Done generating GAPIT Presence Absence Matrix.\n\n");
            }
            if (settings.generateGRM)
            {
                //printf("Attempting to generate GRM:\n");
                GenerateGeneticRelationMatrix(dataset, settings.OUTPUT_DIRECTORY + "gene_relation_matrix.csv", true);
                //printf("Done generating GRM.\n\n");
            }
            if (settings.scoreSynteny)
            {
                //printf("Attempting to score synteny for gene groups:\n");
                dataset->ScoreSynteny(settings.OUTPUT_DIRECTORY + "groups_synteny.csv");
                //printf("Done scoring synteny for gene groups.\n\n");
            }
            if (!settings.GENE_PAIR_FILE.empty())
            {
                //printf("Attempting to score synteny for gene pairs:\n");
                SyntenizeGenePairs(dataset, settings.GENE_PAIR_FILE, settings.OUTPUT_DIRECTORY + "gene_pairs_synteny.csv");
                //printf("Done scoring synteny for gene pairs.\n\n");
            }
            if (settings.classifyFragments && !settings.CHROMOSOME_IDENTIFIER.empty() && !settings.PLASMID_IDENTIFIER.empty())
            {
                //printf("Attempting to classify fragment contigs:\n");
                ClassifyFragments(settings.OUTPUT_DIRECTORY + "classified_fragments.csv", db);
                //printf("Done classifying fragment contigs.\n\n");
            }
            
            if (settings.generateStrainCharts)
            {
                //printf("Attempting to generate strains chart data:\n");
                ExportStrainsChartData(settings.OUTPUT_DIRECTORY, db);
                //printf("Done exporting generate strains chart data.\n\n");
            }

            if (settings.generateSyntenyCharts)
            {
                //printf("Attempting to generate gene groups synteny chart data:\n");
                ExportSyntenyChartData(settings.OUTPUT_DIRECTORY, dataset);
                //printf("Done generating gene groups synteny chart data.\n\n");
            }
            
            if (settings.disambiguateGroups)
            {
                CalculateAverageParalogSynteny(dataset);
                //printf("Attempting to disambiguate gene groups:\n");
                Disambiguate(settings.OUTPUT_DIRECTORY + "disambiguated_gene_groups.csv", dataset);
                //printf("Done disambiguating gene groups.\n\n");
            }

            if (settings.declusterGroups > -1)
            {
                //printf("Attempting to decluster gene groups:\n");
                Decluster(settings.OUTPUT_DIRECTORY + "declustered_gene_groups.csv", dataset, settings.declusterGroups);
                //printf("Done declustering gene groups.\n\n");
            }

            if (settings.generateFNAGroups)
            {
                //printf("Attempting to generate group fna files:\n");
                dataset->ExportFNA(settings.OUTPUT_DIRECTORY);
                //printf("Done generating group fna files.\n\n");
            }
        }

        //printf("Attempting to export strain metrics:\n");
        ExportStatistics(dataset, settings.OUTPUT_DIRECTORY + "strain_metrics.csv");
        //printf("Done exporting strain metrics.\n\n");

        /*
        printf("Attempting export all calculated synteny scores:\n");
        ExportGeneScores(dataset, settings.OUTPUT_DIRECTORY + ".synteny_scores.csv");
        printf("Done exporting all calculated synteny scores.\n\n");
        */
    }

    printf("Everything took:\t%i seconds\n", (int)timer.Elapsed());

    return 0;
}
