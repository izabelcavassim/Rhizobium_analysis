//
//  Group.cpp
//  Syntenizer3000
//
//  Created by Camous Moslemi on 23/04/2017.
//
//

#include "Group.hpp"

//using namespace std;

bool Group::HasGene(Gene* g)
{
    for (auto gene = genes.begin(); gene != genes.end(); gene++)
        if (*gene == g)
            return true;
    
    return false;
}

void Group::RemoveGene(Gene* g)
{
    for (auto gene = genes.begin(); gene != genes.end(); gene++)
        if (*gene == g)
            genes.erase(gene);
}

void Group::InsertGenes(vector<Gene*> newGenes, unordered_map<Gene *, Group *> *grouppool)
{
    for (auto gene = newGenes.begin(); gene != newGenes.end(); gene++)
        if (!this->HasGene(*gene))
        {
            //printf("moving gene %s into %s.\n", (*gene)->id.c_str(), this->id.c_str());
            genes.push_back(*gene);
            (*gene)->group = this;
            (*grouppool)[*gene] = this;
        }
}

int Group::SharedGenes(vector<Gene*> newGenes)
{
    int shared = 0;
    vector<Gene*>::iterator gene;
    
    for (gene = newGenes.begin(); gene != newGenes.end(); gene++)
        if (!this->HasGene(*gene))
            shared++;
    
    return shared;
}

double Group::SyntenizeSophisticated()
{
    int count = 0;
    double scoresum = 0.0f;
    
    for (auto g1 = genes.begin(); g1 != genes.end(); g1++)
        for (auto g2 = g1+1; g2 != genes.end(); g2++)
            if (*g1 != *g2)
            {
                count++;
                double score = Score::Sophisticated(*g1, *g2);

                (*g1)->score += score;
                (*g2)->score += score;
                scoresum += score;
            }

    syntenyScoreSophisticated = scoresum / (double)count;
    return syntenyScoreSophisticated;
}

double Group::SyntenizeOld()
{
    int count = 0;
    double scoresum = 0.0f;
    for (auto g1 = genes.begin(); g1 != genes.end(); g1++)
    {
        for (auto g2 = g1+1; g2 != genes.end(); g2++)
        {
            if (*g1 != *g2)
            {
                count++;
                scoresum += Score::Old(*g1, *g2);
            }
        }
    }
    
    syntenyScoreOld = (double)scoresum / (double)count;
    
    return syntenyScoreOld;
    
}

double Group::SyntenizeSimple()
{
    int count = 0;
    double scoresum = 0.0f;
    for (auto g1 = genes.begin(); g1 != genes.end(); g1++)
    {
        for (auto g2 = g1+1; g2 != genes.end(); g2++)
        {
            if (*g1 != *g2)
            {
                count++;
                scoresum += Score::Simple(*g1, *g2);
            }
        }
    }
    
    syntenyScoreSimple = (double)scoresum / (double)count;
    
    return syntenyScoreSimple;
}

double Group::SyntenizeAdjusted()
{
    int count = 0;
    double scoresum = 0.0f;
    for (auto g1 = genes.begin(); g1 != genes.end(); g1++)
    {
        for (auto g2 = g1+1; g2 != genes.end(); g2++)
        {
            if (*g1 != *g2)
            {
                count++;
                scoresum += Score::Adjusted(*g1, *g2);
            }
        }
    }
    
    syntenyScoreAdjusted = (double)scoresum / (double)count;
    
    return syntenyScoreAdjusted;
    
}

double Group::SyntenizeFast()
{
    int count = 0;
    double scoresum = 0.0f;
    
    for (auto g1 = genes.begin(); g1 != genes.end(); g1++)
        for (auto g2 = g1+1; g2 != genes.end(); g2++)
            if (*g1 != *g2)
            {
                count++;
                double score = Score::Fast(*g1, *g2);
                //(*g1)->score += score;
                //(*g2)->score += score;
                scoresum += score;
            }
    
    syntenyScoreFast = scoresum / (double)count;
    return syntenyScoreFast;
}

double Group::SyntenizeSimpleAgainstGene(Gene *g1)
{
    double count = 0;
    double scoresum = 0;

    for (auto g2 = genes.begin(); g2 != genes.end(); g2++)
        if (g1 != *g2)
        {
            count++;
            scoresum += Score::Simple(g1, *g2);
        }

    return(scoresum / count);
}

double Group::SyntenizeAgainstGene(Gene *g1)
{
    double count = 0;
    double scoresum = 0;

    for (auto g2 = genes.begin(); g2 != genes.end(); g2++)
        if (g1 != *g2)
        {
            count++;
            scoresum += Score::Sophisticated(g1, *g2);
            //scoresum += Score::Simple(g1, *g2);
        }

    return(scoresum / count);
}

int Group::CountUniqueStrains()
{
    set<Strain *> uniquestrains;
    for (auto gene = this->genes.begin(); gene != this->genes.end(); gene++)
        uniquestrains.insert((*gene)->strain);

    return uniquestrains.size();
}

int Group::CountOrthologs()
{
    int orthologs = 0;
    vector<Strain *> strains;
    
    for (auto gene = this->genes.begin(); gene != this->genes.end(); gene++)
        strains.push_back((*gene)->strain);

    for (auto gene = this->genes.begin(); gene != this->genes.end(); gene++)
    {
        if (count(strains.begin(), strains.end(), (*gene)->strain) == 1)
            orthologs++;
    }

    return orthologs;
}

int Group::CountParalogs()
{
    int paralogs = 0;
    vector<Strain *> strains;
    
    for (auto gene = this->genes.begin(); gene != this->genes.end(); gene++)
        strains.push_back((*gene)->strain);
    
    for (auto gene = this->genes.begin(); gene != this->genes.end(); gene++)
    {
        if (count(strains.begin(), strains.end(), (*gene)->strain) > 1)
            paralogs++;
    }
    
    return paralogs;
    /*
    int paralogs = 0;
    set<Strain *> strains;
    set<Strain *> paralogstrains;
    for (auto gene = this->genes.begin(); gene != this->genes.end(); gene++)
    {
        if (strains.count((*gene)->strain))
            paralogstrains.insert((*gene)->strain);
        else
            strains.insert((*gene)->strain);
    }
    
    for (auto gene = this->genes.begin(); gene != this->genes.end(); gene++)
        if (paralogstrains.count((*gene)->strain))
            paralogs++;
    
    return paralogs;
     */
}

void Group::Prune()
{
    int pure_paralog_groups = 0;
    int unambiguous_groups = 0;
    int ambiguous_groups = 0;
    int split_groups = 0;

    bool has_paralogs = false;
    unordered_map<Strain *, vector<Gene*>*> strains;

    for (auto gene = this->genes.begin(); gene != this->genes.end(); gene++)
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
        vector <Gene*> paralogs;
        vector <Gene*> homologs;
        
        for (auto strain = strains.begin(); strain != strains.end(); strain++)
        {
            if ((*strain).second->size() > 1)
                paralogs.insert(paralogs.end(), (*strain).second->begin(), (*strain).second->end());
            else
                homologs.insert(homologs.end(), (*strain).second->begin(), (*strain).second->end());
        }
        
        if (homologs.size() == 0)
        {
            pure_paralog_groups++;
            return;
        }
        
        Strain *prev_strain = (*paralogs.begin())->strain;
        
        bool split = false;
        int total_scoresum = 0;
        int max_scoresum = 0;
        for (auto g1 = paralogs.begin(); g1 != paralogs.end(); g1++)
        {
            double scoresum = 0.0f;
            
            for (auto g2 = homologs.begin(); g2 != homologs.end(); g2++)
                scoresum += Score::Sophisticated(*g1, *g2);
            
            if (prev_strain != (*g1)->strain)
            {
                prev_strain = (*g1)->strain;
            }
            scoresum /= (double)homologs.size();
            
            for (auto g2 = paralogs.begin(); g2 != paralogs.end(); g2++)
                if (g1 != g2 && scoresum > 0 && Score::Sophisticated(*g1, *g2) == 0)
                    split = true;
            
            total_scoresum += scoresum;
            if (scoresum > max_scoresum)
                max_scoresum = scoresum;
        }
        if (split)
            split_groups++;
        
        if (max_scoresum > 0 && total_scoresum == max_scoresum)
            unambiguous_groups++;
        else
            ambiguous_groups++;
    }
}

bool Group::HasStrain(Strain *strain)
{
    for (auto gene = genes.begin(); gene != genes.end(); gene++)
        if ((*gene)->strain == strain)
            return true;

    return false;
}

vector<Group *> Group::Decluster(double threshold)
{
    vector<Group *> splitgroups;

    do
    {
        for (auto g1 = genes.begin(); g1 != genes.end(); g1++)
            (*g1)->score = 0;
        
        for (auto g1 = genes.begin(); g1 != genes.end(); g1++)
        {
            for (auto g2 = g1+1; g2 != genes.end(); g2++)
            {
                double score = Score::Sophisticated(*g1, *g2);
//                double score = Score::Simple(*g1, *g2);

                (*g1)->score += score;
                (*g2)->score += score;
            }
        }

        vector<Gene *> erase;
        Gene *bestcandidate = NULL;
        double bestscore = 0;
        for (auto g1 = genes.begin(); g1 != genes.end(); g1++)
        {
            if ((*g1)->score == 0)
            {
                (*g1)->group = NULL;
                erase.push_back(*g1);
            }
            if ((*g1)->score > bestscore)
            {
                bestscore = (*g1)->score;
                bestcandidate = *g1;
            }
        }

        for (auto gene = erase.begin(); gene != erase.end(); gene++)
            genes.erase(find(genes.begin(), genes.end(), *gene));

        if (bestscore == 0)
            return splitgroups;

        genes.erase(find(genes.begin(), genes.end(), bestcandidate));
        
        Group *newgroup = new Group;
//        newgroup->id = id;
//        if (splitgroups.size() > 0)
//            newgroup->id += "_split" + to_string(splitgroups.size());
        newgroup->id = id + "_split" + to_string(splitgroups.size()+1);

        bestcandidate->group = newgroup;
        bestcandidate->group->genes.push_back(bestcandidate);
        splitgroups.push_back(newgroup);

        bool inserted;
        do
        {
            inserted = false;
            Gene *bestcandidate = NULL;
            double bestscore = 0;
            for (auto g1 = genes.begin(); g1 != genes.end(); g1++)
            {
                double score = newgroup->SyntenizeAgainstGene(*g1);
//                double score = newgroup->SyntenizeSimpleAgainstGene(*g1);
                if (score > bestscore)
                {
                    bestscore = score;
                    bestcandidate = *g1;
                }
            }

            if (bestscore > threshold)
            {
                inserted = true;
                bestcandidate->group = newgroup;
                newgroup->genes.push_back(bestcandidate);
                genes.erase(find(genes.begin(), genes.end(), bestcandidate));
            }
        }
        while (inserted);
    }
    while (genes.size() > 0);

    return splitgroups;
}

vector<Group *> Group::Disambiguate()
{
    vector<Group *> splitgroups;
    vector<Gene *> candidates;
    vector<Gene *> paralogs;
    //vector<Gene *> orthologs;
    set<Strain *> strains;
    set<Strain *> paralogstrains;
    Group *orthologs = new Group;

    vector<double> orthologscores;
    vector<double> paralogscores;

    //printf("%s [%lu] (%i):", id.c_str(), genes.size(), CountParalogs());
    //printf("%s:", id.c_str());
    //size_t loc;
    //if ((loc = id.find(' ')) != string::npos)
    //    id.erase(id.begin()+loc,id.end());

    printf("%s\t%lu\t%.2f\t%i\t", id.c_str(), genes.size(), SyntenizeSophisticated(), CountParalogs());

    orthologs->id = id;
    splitgroups.push_back(orthologs);

    //Make two groups, one of only strains that have a single gene in the group (orthologs), and one with more than one gene (paralogs)
    for (auto gene = genes.begin(); gene != genes.end(); gene++)
        if (strains.count((*gene)->strain))
            paralogstrains.insert((*gene)->strain);
        else
            strains.insert((*gene)->strain);

    for (auto gene = genes.begin(); gene != genes.end(); gene++)
        if (paralogstrains.count((*gene)->strain))
            paralogs.push_back(*gene);
        else
            orthologs->genes.push_back(*gene);

    if (orthologs->genes.size() > 0)
    {
        for (auto strain = paralogstrains.begin(); strain != paralogstrains.end(); strain++)
        {
            Gene *bestcandidate = NULL;
            Gene *nextbestcandidate = NULL;
            for (auto gene = paralogs.begin(); gene != paralogs.end(); gene++)
            {
                if ((*gene)->strain != *strain)
                    continue;

                (*gene)->score = orthologs->SyntenizeAgainstGene(*gene);

                paralogscores.push_back( (*gene)->score);

                if (bestcandidate == NULL || (*gene)->score > bestcandidate->score)
                {
                    if (bestcandidate != NULL)
                        nextbestcandidate = bestcandidate;
                    bestcandidate = *gene;
                }
            }

            if (bestcandidate == NULL)
            {
                printf("Error!\n");
                exit(1);
            }

            if (nextbestcandidate != NULL && bestcandidate->score == nextbestcandidate->score)
            {
                printf("double jeopardy\n!");
            }

            orthologscores.push_back(bestcandidate->score);
            paralogscores.erase(find(paralogscores.begin(), paralogscores.end(), bestcandidate->score));

            orthologs->genes.push_back(bestcandidate);
            paralogs.erase(find(paralogs.begin(), paralogs.end(), bestcandidate));
        }
    }

    do
    {
        // Find the paralog with the highest average synteny with other paralogs and form a group around it.
        Gene *candidate = NULL;
        for (auto g1 = paralogs.begin(); g1 != paralogs.end(); g1++)
        {
            for (auto g2 = g1+1; g2 != paralogs.end(); g2++)
            {
                double score = Score::Sophisticated(*g1, *g2);
                //double score = Score::Simple(*g1, *g2);
                (*g1)->score += score;
                (*g2)->score += score;

                if (candidate == NULL || candidate->score < (*g1)->score)
                    candidate = *g1;

                if (candidate == NULL || candidate->score < (*g2)->score)
                    candidate = *g2;
            }
        }

        if (paralogs.size() == 1)
            candidate = paralogs[0];

        if (orthologs->genes.size() > 0)
        {
            candidate->group = new Group;
            candidate->group->id = id + "-" + to_string(splitgroups.size()+1);
            candidate->group->genes.push_back(candidate);
            splitgroups.push_back(candidate->group);
        }
        else
        {
            orthologs->genes.push_back(candidate);
            candidate->group = orthologs;
        }
        candidates.push_back(candidate);
        paralogs.erase(find(paralogs.begin(), paralogs.end(), candidate));

        // Put any gene scoring above 0 synteny into the newly formed group
        bool restart;
        do
        {
            restart = false;
            for (auto g1 = paralogs.begin(); g1 != paralogs.end(); g1++)
            {
                if (!candidate->group->HasStrain((*g1)->strain) && candidate->group->SyntenizeAgainstGene(*g1) > 0 )
                {
                    orthologscores.push_back(candidate->group->SyntenizeAgainstGene(*g1));
                    candidate->group->genes.push_back(*g1);
                    paralogs.erase(g1);
                    restart = true;
                    break;
                }
                else if (candidate->group->HasStrain((*g1)->strain))
                    paralogscores.push_back(candidate->group->SyntenizeAgainstGene(*g1));
            }
        }
        while(restart);

        // Repeat until no paralogs remain
    }
    while (paralogs.size() > 0);

    //Take all single gene groups and try to make them part of the largest of the new groups that does not have a gene of the same strain in it.
    vector<Gene *> orphans;
    for (auto group = splitgroups.begin(); group != splitgroups.end(); group++)
    {
        if ((*group)->genes.size() == 1)
            orphans.push_back((*group)->genes[0]);
    }

    //splitgroups.erase(remove_if(splitgroups.begin(), splitgroups.end(), [](Group *group) { return group->genes.size() == 1; }), splitgroups.end());

    int unorphaned = 0;
    for (auto gene = orphans.begin(); gene != orphans.end(); gene++)
    {
        Group *bestgroup = NULL;

        for (auto group = splitgroups.begin(); group != splitgroups.end(); group++)
        {
            if ((*group)->HasStrain((*gene)->strain) || (*gene)->group == *group)
                continue;
            if (bestgroup == NULL || (*group)->genes.size() > bestgroup->genes.size())
                bestgroup = *group;
        }

        if (bestgroup != NULL)
        {
            bestgroup->genes.push_back(*gene);
            unorphaned++;

            splitgroups.erase(find(splitgroups.begin(), splitgroups.end(), (*gene)->group));
        }
    }

    this->orphans = 0;
    for (auto group = splitgroups.begin(); group != splitgroups.end(); group++)
        if ((*group)->genes.size() <= 1)
            this->orphans++;
    /**/
    printf("\t%lu\t", splitgroups.size());
    for (auto group = splitgroups.begin(); group != splitgroups.end(); group++)
        if ((*group)->genes.size() > 1)
            printf("%lu\t%.2f\t", (*group)->genes.size(), (*group)->SyntenizeSophisticated());
        else
        {
            this->orphans++;
//            printf("\t%s\t%.2f", (*group)->genes[0]->strain->id.c_str(), orthologs->SyntenizeAgainstGene((*group)->genes[0]));
        }
    /**/
    this->orphans = orphans.size() - unorphaned;
    /*
    printf("%lu\t%f\t%lu\t%f\n",
           orthologscores.size(),
           accumulate(orthologscores.begin(), orthologscores.end(), 1.0f)/(double)orthologscores.size(),
           paralogscores.size(),
           accumulate(paralogscores.begin(), paralogscores.end(), 1.0f)/(double)paralogscores.size());
     */
    printf("\n");
    return splitgroups;
}

void Group::GenerateSyntenyChart(string OUTPUT_DIRECTORY )
{
    ofstream out( OUTPUT_DIRECTORY + id + "_synteny_chart.csv", ifstream::out );

    //printf("Generating Synteny chart for %s\n", id.c_str());
    SyntenizeSophisticated();

    set<Group *> groups;
    double bestscore = -1;
    Gene *candidate = NULL;
    for (auto gene = genes.begin(); gene != genes.end(); gene++)
    {
        for (int i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
            if ((*gene)->neighbours[i] != NULL )
                groups.insert((*gene)->neighbours[i]->group);

        if (candidate == NULL || (*gene)->score > bestscore)
        {
            bestscore = (*gene)->score;
            candidate = *gene;
        }
    }

    //printf("Best candidate %s: %.2f\n", candidate->id.c_str(), candidate->score );

    unordered_map<Group *, int> groupcolours;
    int colour = 0;
    for (auto group = groups.begin(); group != groups.end(); group++)
        groupcolours[*group] = ++colour;

    vector<Gene *> remaining;
    vector<Gene *> sorted;
    sorted.push_back(candidate);
    for (auto gene = genes.begin(); gene != genes.end(); gene++)
        if ((*gene) != candidate)
            remaining.push_back(*gene);

    do
    {
        bestscore = -1;
        Gene *bestmatch = NULL;
        for (auto gene = remaining.begin(); gene != remaining.end(); gene++)
        {
            //double score = Score::Sophisticated(candidate, *gene);
            double score = Score::Simple(candidate, *gene);
            if (bestmatch == NULL || score > bestscore)
            {
                bestscore = score;
                bestmatch = *gene;
            }
        }
        //printf("Best match %s: %.2f\n", bestmatch->id.c_str(), bestscore );

        sorted.push_back(bestmatch);
        remaining.erase(find(remaining.begin(), remaining.end(), bestmatch));
        candidate = bestmatch;
    }
    while (remaining.size() > 0);

    for (auto gene = sorted.begin(); gene != sorted.end(); gene++)
    {
        out << (*gene)->strain->id;
        if (*gene != sorted.back())
            out << ';';
        else
            out << "\n";
    }

    for (int i = 0; i < sorted.size(); i++)
        out << "0" << (i < sorted.size()-1 ? ";" : "\n");
    for (int i = 0; i < sorted.size(); i++)
        out << SIZE_OF_NEIGHTBOURHOOD+2 << (i < sorted.size()-1 ? ";" : "\n");
    for (int i = 0; i < sorted.size(); i++)
        out << (SIZE_OF_NEIGHTBOURHOOD/2)+1 << (i < sorted.size()-1 ? ";" : "\n");

    for (auto line = groupcolours.begin(); line != groupcolours.end(); line++)
    {
        Group *group = line->first;
        int colour = line->second;

        for (int s = 0; s < sorted.size(); s++)
        {
            Gene * g1 = sorted[s];
            for (int i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
                if (g1->neighbours[i] != NULL && g1->neighbours[i]->group == group)
                {
                    if (i+1 <= (SIZE_OF_NEIGHTBOURHOOD/2))
                        out << i+1;
                    else
                        out << i+2;
                    break;
                }

            if (s < sorted.size()-1)
                out << ';';
            else
                out << '\n';
        }
    }
    out.close();
}

Group::Group(void)
{
    coorthologs = -1;
    algebraicConnectivity = NAN;
    syntenyScoreSophisticated = NAN;
    syntenyScoreOld = NAN;
    syntenyScoreSimple = NAN;
    syntenyScoreAdjusted = NAN;
    syntenyScoreFast = NAN;
    syntenyScoreTest = NAN;
}
