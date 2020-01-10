//
//  Dataset.cpp
//  Syntenizer3000
//
//  Created by Camous Moslemi on 11/09/2017.
//
//

#include "Dataset.hpp"

Dataset::Dataset(Database *db)
{
    database = db;
    genecount = 0;
    syntenyScoreSophisticated = NAN;
    syntenyScoreOld = NAN;
    syntenyScoreSimple = NAN;
    syntenyScoreAdjusted = NAN;
    syntenyScoreTest = NAN;
    syntenyScoreFast = NAN;
}

void Dataset::ScoreSynteny(string OUTPUT_FILE)
{
    Timer timer;

    int counter = 0;
    double scoresumSophisticated = 0.0f;
    //double scoresumOld = 0.0f;
    double scoresumSimple = 0.0f;
    //double scoresumAdjusted = 0.0f;
    //double scoresumFast = 0.0f;
    timer.Start();
    printf("Generating gene group synteny scores into:\t%s\n", OUTPUT_FILE.c_str());

    ofstream out( OUTPUT_FILE, ifstream::out );
    Progress p(groups.size()-1);

    //out << "Group\tConnectivity\tGenes\tStrains\tOrthologs\tParalogs\tChromosome Genes\tPlasmid Genes\tFragment Genes\tProtein Products\tSynteny Score Old\tSynteny Score Adjusted\tSynteny Score Fast\tSynteny Score Simple\tSynteny Score Sophisticated\n";
    out << "Group;Connectivity;Genes;Strains;Orthologs;Paralogs;Chromosome Genes;Plasmid Genes;Fragment Genes;Protein Products;Synteny Score SIMPLE;Synteny Score SOPHISTICATED\n";
    for (auto group = groups.begin(); group != groups.end(); group++)
    {
        set<Group*> groupset;
        if ((*group)->genes.size() > 1)
        {
            scoresumSophisticated += (*group)->SyntenizeSophisticated();
            //scoresumOld += (*group)->SyntenizeOld();
            scoresumSimple += (*group)->SyntenizeSimple();
            //scoresumAdjusted += (*group)->SyntenizeAdjusted();
            //scoresumFast += (*group)->SyntenizeFast();
            counter++;
            set<string> types;

            int chromosome = 0;
            int plasmid = 0;
            int fragment = 0;

            for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
            {
                if (!settings.CHROMOSOME_IDENTIFIER.empty() && (*gene)->contig->id.find(settings.CHROMOSOME_IDENTIFIER) != string::npos)
                    chromosome++;
                if (!settings.PLASMID_IDENTIFIER.empty() && (*gene)->contig->id.find(settings.PLASMID_IDENTIFIER) != string::npos)
                    plasmid++;
                if (!settings.FRAGMENT_IDENTIFIER.empty() && (*gene)->contig->id.find(settings.FRAGMENT_IDENTIFIER) != string::npos)
                    fragment++;

                types.insert((*gene)->product);
            }

            out << (*group)->id.c_str() << ';';
            out << (*group)->algebraicConnectivity << ';';
            out << (*group)->genes.size() << ';';
            out << (*group)->CountUniqueStrains() << ';';
            out << (*group)->CountOrthologs() << ';';
            out << (*group)->CountParalogs() << ';';

            out << (!settings.CHROMOSOME_IDENTIFIER.empty() ? chromosome : '?') << ';';
            out << (!settings.PLASMID_IDENTIFIER.empty() ? plasmid : '?') << ';';
            out << (!settings.FRAGMENT_IDENTIFIER.empty() ? fragment : '?') << ';';

            int counter = 1;
            for (auto type = types.begin(); type != types.end(); type++, counter++)
            {
                out << (*type);
                if (counter < types.size())
                    out << '|';
                else
                    out << ';';
            }

            //out << (*group)->syntenyScoreOld << ';';
            //out << (*group)->syntenyScoreAdjusted << ';';
            //out << (*group)->syntenyScoreFast << ';';
            out << (*group)->syntenyScoreSimple << ';';
            out << (*group)->syntenyScoreSophisticated << "\n";

            p.Update(distance(groups.begin(), group));
        }
    }

    syntenyScoreSophisticated = scoresumSophisticated / (double)counter;
    //syntenyScoreOld = scoresumOld / (double)counter;
    syntenyScoreSimple = scoresumSimple / (double)counter;
    //syntenyScoreAdjusted = scoresumAdjusted / (double)counter;
    //syntenyScoreFast = scoresumFast / (double)counter;

    //printf("\nAverage Synteny Score Old:\t[%.2f]", syntenyScoreOld);
    //printf("\nAverage Synteny Score Adjusted:\t[%.2f]", syntenyScoreAdjusted);
    //printf("\nAverage Synteny Score Fast:\t[%.2f]", syntenyScoreFast);
    printf("Average gene group SIMPLE synteny score:\t%.2f\n", syntenyScoreSimple);
    printf("Average gene groups SOPHISTICATED synteny Score:\t%.2f\n", syntenyScoreSophisticated);
    printf("Generating synteny scores for %lu gene groups took:\t%i seconds\n\n", groups.size(), (int)timer.Elapsed());

    /*
    out << "\n\nAverage Synteny Score Sophisticated\t" << syntenyScoreSophisticated;
    out << "\n\nAverage Synteny Score Old\t" << syntenyScoreOld;
    out << "\n\nAverage Synteny Score Simple\t" << syntenyScoreSimple;
    out << "\n\nAverage Synteny Score Adjusted\t" << syntenyScoreAdjusted;
    out << "\n\nAverage Synteny Score Fast\t" << syntenyScoreFast;
    out << "\n";
    */
    out.close();
}

void Dataset::ExportFNA(string OUTPUT_DIRECTORY)
{
    Timer timer;
    timer.Start();
    printf("Exporting fna gene groups into:\t%s\n", OUTPUT_DIRECTORY.c_str());
    for (auto group = groups.begin(); group != groups.end(); group++)
    {
        ofstream out( OUTPUT_DIRECTORY + (*group)->id + ".fna", ifstream::out );
        for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
        {
            out << '>' << (*group)->id << '|' << (*gene)->strain->id << '|' << (*gene)->id << '\n';
            out << (*gene)->sequence << '\n';
        }
        out.close();
    }
    printf("Exporting %lu fna gene groups took:\t%i seconds\n\n", groups.size(), (int)timer.Elapsed());
}
/*
void Dataset::PrintStats()
{
    //database->strains
    //unordered_map<int, Strain *>
    for (int strainsize = database->strains.size(); strainsize > 1; strainsize-- )
    {
        vector<Strain *> strains;
        for (auto group = groups.begin() ; group < groups.end(); group++ )
        {
            set<Strain *> strainset;
            for (auto gene = (*group)->genes.begin(); gene != (*group)->genes.end(); gene++)
                strainset.insert((*gene)->strain);
            
            if (strainset.size() != strainsize)
                continue;

            vector<Strain *> missing;
        }
    }
}
*/
