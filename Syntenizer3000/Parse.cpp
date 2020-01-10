//
//  Parse.cpp
//  Syntenizer3000
//
//  Created by Camous Moslemi on 01/09/2017.
//
//

#include "Parse.hpp"

//using namespace std;

int Parse::GeneGroups(string INPUT_FILE, Dataset *dataset)
{
    Timer timer;
    bool proteinorthomode = false;
    int genes = 0;
    int groups = 0;
    int shortContigGenes = 0;
    int shortContigGroups = 0;
    int coorthologousgroups = 0;
    double coorthologs = 0.0f;

    size_t loc;
    printf("Parsing gene groups from:\t%s\n", INPUT_FILE.c_str());
    
    ifstream input = ifstream( INPUT_FILE);
    if (!input.is_open())
    {
        printf("\nError: Could not open file %s\n", INPUT_FILE.c_str() );
        exit(1);
    }

    timer.Start();
    
    string line;
    int linenr = 0;
    vector<string> header;
    while (getline(input, line))
    {
        if (line.empty())
            continue;
        
        if (++linenr == 1 && line[0] == '#')
            proteinorthomode = true;

        if (proteinorthomode)
        {
            stringstream sstream(line);
            string item;
            vector<string> items;
            Group *group;

            while (getline(sstream, item, '\t'))
            items.push_back(item);
            
            if (linenr == 1)
            {
                for (int i = 0; i < items.size(); i++)
                {
                    item = items[i];
                    
                    if (i >= 3 && (loc = item.find('.')) != string::npos)
                        item = item.substr(0, loc);
                    
                    header.push_back(item);
                }
                continue;
            }
            
            int numspecies = stoi(items[0]);
            int numgenes = stoi(items[1]);
            string connectivity = items[2];
            bool shortContigFound = false;
            
            groups++;
            group = new Group;
            group->id = "group" + to_string(groups);// + "|" + connectivity;
            group->coorthologs = numgenes - numspecies;
            group->algebraicConnectivity = stod( items[2] );
            dataset->groups.push_back(group);
            
            for (int i = 3; i < items.size(); i++)
            {
                stringstream sstream2(items[i]);
                
                vector<string> geneIDs;
                while (getline(sstream2, item, ','))
                    geneIDs.push_back(item);
                
                for (int j = 0; j < geneIDs.size(); j++)
                {
                    if (geneIDs[j] == "*")
                        continue;
                    
                    Gene *gene = NULL;
                    
                    string strainID = header[i];
                    string geneID = geneIDs[j];
                    
                    try
                    {
                        gene = dataset->database->genepool.at(strainID + '|' + geneID);
                        if (gene->contig->genes.size() < SIZE_OF_NEIGHTBOURHOOD + 1)
                        {
                            shortContigGenes++;
                            shortContigFound = true;
                        }
                    }
                    catch (const out_of_range & e)
                    {
                        printf("\nError: Gene \"%s\" in Strain \"%s\" not found!\n", geneID.c_str(), strainID.c_str() );
                        exit(1);
                    }
                    gene->group = group;
                    group->genes.push_back(gene);
                    dataset->grouppool[gene] = group;
                    dataset->grouppool2[group->id] = group;
                    genes++;
                }
            }

            if (group->genes.size() <= 1)
            {
                printf("\nError: Group \"%s\" does not have the minimum 2 genes required number of members.\n", group->id.c_str());
                exit(1);
            }
            
            if (group->coorthologs > 0)
                coorthologousgroups++;
            if (shortContigFound)
                shortContigGroups++;
            
            coorthologs += group->coorthologs;
        }
        else
        {
            stringstream sstream(line);
            string item;
            Group *group;
            
            if (line.empty())
                continue;
            
            group = new Group;
            dataset->groups.push_back(group);
            
            bool shortContigFound = false;
            set<Strain*> strainset;
            vector<Strain*> strains;
            while (getline(sstream, item, ' '))
            {
                bool idFound = false;
                size_t loc;
                
                if (item.empty())
                    continue;
                
                if (!idFound)
                {
                    loc = item.find(":");
                    
                    if (loc != string::npos)
                    {
                        idFound = true;
                        group->id = item.substr(0, loc).c_str();
                        if ((loc = group->id.find("|"))  != string::npos)
                        {
                            group->algebraicConnectivity = stod(group->id.substr(loc+1, group->id.size()));
                            group->id = group->id.substr(0, loc);
                        }
                        
                        continue;
                    }
                }
                
                loc = item.find("|");
                
                if (loc != string::npos)
                {
                    Gene *gene = NULL;
                    
                    string strainID = item.substr(0, loc);
                    string geneID = item.substr(loc+1, item.size());

                    try
                    {
                        gene = dataset->database->genepool.at(strainID + '|' + geneID);
                        if (gene->contig->genes.size() < SIZE_OF_NEIGHTBOURHOOD + 1)
                        {
                            shortContigGenes++;
                            shortContigFound = true;
                        }
                    }
                    catch (const out_of_range & e)
                    {
                        printf("\nError: Gene \"%s\" in Strain \"%s\" not found!\n", geneID.c_str(), strainID.c_str() );
                        exit(1);
                    }
                    gene->group = group;
                    group->genes.push_back(gene);
                    dataset->grouppool[gene] = group;
                    dataset->grouppool2[group->id] = group;
                    strainset.insert(gene->strain);
                    strains.push_back(gene->strain);
                    genes++;
                }
            }
            if (group->genes.size() <= 1)
            {
                printf("\nError: Group \"%s\" does not have the minimum 2 genes required number of members.\n", group->id.c_str());
                exit(1);
            }
            groups++;
            group->coorthologs = group->genes.size()-strainset.size();
            if (group->coorthologs > 0)
                coorthologousgroups++;
            if (shortContigFound)
                shortContigGroups++;

            coorthologs += group->coorthologs;
        }
    }

    input.close();
    printf("Short Contig Genes:\t%i\t(%.2f%%)\n", shortContigGenes, (shortContigGenes * 100.0f) / genes );
    printf("Short Contig Groups:\t%i\t(%.2f%%)\n", shortContigGroups, (shortContigGroups * 100.0f) / groups );
    printf("Co-orthologs:\t%i\t(%.2f%%)\n", (int)coorthologs, coorthologs*100.0f / (double)genes);
    printf("Co-orthologous groups:\t%i\t(%.2f%%)\n", coorthologousgroups, (double)coorthologousgroups*100.0f / (double)dataset->groups.size());
    printf("Parsing %i genes in %i groups took:\t%i seconds\n\n", genes, groups, (int)timer.Elapsed());
    
    return dataset->groups.size();
}

int Parse::Strains(string INPUT_DIRECTORY, Database *database)
{
    int genes = 0;
    int strains = 0;
    int contigs = 0;
    int shortcontigs = 0;
    Timer timer;
    DIR *dirp = opendir(INPUT_DIRECTORY.c_str());
    struct dirent *dp;

    timer.Start();
    printf("Parsing strains from:\t%s\n", INPUT_DIRECTORY.c_str());
    
    printf("Processing [");
    while ((dp = readdir(dirp)) != NULL)
    {
        string strainID = string(dp->d_name);
        
        if (strainID.size() < 5)
            continue;

        string type;

        if (!strainID.compare(strainID.size()-4,4,".gff"))
            type = TYPE_GFF;
        else if (!strainID.compare(strainID.size()-4,4,".tab"))
            type = TYPE_TAB;
        else
            continue;
    
        ifstream input = ifstream( INPUT_DIRECTORY + strainID);
        Strain *strain = new Strain();
        strainID.resize(strainID.size()-4);
        strain->id = strainID;

        database->strains.push_back(strain);
        database->strainpool[strain->id] = strain;
        strains++;

        int linenr = 0;
        string line;
        unordered_map<string, int> regionsize;
        while (getline(input, line) && line != "##FASTA")
        {
            int counter = 0;
            linenr++;

            bool issequenceregion = false;
            if (type == TYPE_GFF)
            {
                if (line.find("##sequence-region ") != string::npos)
                    issequenceregion = true;

                if (line[0] == '#' && !issequenceregion)
                    continue;

                stringstream sstream(line);
                string s;

                vector<string> items;
                //items.clear();

                char delimiter = issequenceregion ? ' ' : '\t';

                while (getline(sstream, s, delimiter))
                {
                    items.push_back(s);
                    counter++;
                    if (counter >= 9)
                        break;
                }

                if (issequenceregion)
                {
                    int size = stoi(items[3]);
                    regionsize[items[1]] = size;
                    strain->bp += size;
                    continue;
                }

                if (items[2] == "tmRNA")
                    strain->tmRNA++;
                if (items[2] == "rRNA")
                    strain->rRNA++;
                if (items[2] == "tRNA")
                    strain->tRNA++;

                if (items[2] == "CDS")// || items[2] == "tRNA")
                {
                    Gene *gene = new Gene;
                    size_t loc1 = items[8].find("ID=");
                    if (loc1 == string::npos)
                    {
                        printf("\nError: No ID tag found on line %i of %s.\n", linenr, strainID.c_str());
                        exit(1);
                    }
                    size_t loc2 = items[8].find_first_of(";\n", loc1);

                    gene->id = items[8].substr(loc1+3, loc2-(loc1+3));
                    //printf("strinID: %s\n", strain->id.c_str());
                    //printf("geneID: %s\n", gene->id.c_str());
                    //exit(0);
                    //gene->colour = "black";
                    loc1 = items[8].find("product=");
                    if (loc1 != string::npos)
                    {
                        loc2 = items[8].find_first_of(";\n", loc1);

                        gene->product = items[8].substr(loc1+8, loc2-(loc1+8));
                        gene->masked = (gene->product.find("ransposase") != string::npos) || (gene->product.find("ntegrase") != string::npos);
                        //printf("Product: %s\n", gene->product.c_str());
                    }
                    gene->start = stoi(items[3]);
                    gene->end = stoi(items[4]);
                    gene->length = abs(gene->start - gene->end)+1;
                    gene->orientation = (items[6] != "-");
                    gene->paralog = false;
                    gene->strain = strain;
                    gene->left_neighbour = NULL;
                    gene->right_neighbour = NULL;
                    for (int i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
                        gene->neighbours[i] = NULL;
                    gene->group = NULL;
                    //gene->relations.clear();
                    gene->score = 0.0f;
                    //gene->coverage = 0;
                    gene->GC3s = NAN;
                    strain->CDS++;
                    /*
                    if (items[2] == "CDS")
                        strain->CDS++;
                    if (items[2] == "tRNA")
                    {
                        //gene->type = "tRNA";
                        strain->tRNA++;
                    }
                     */

                    string contigID = items[0];

                    Contig *contig;

                    try
                    {
                        contig = strain->contigpool.at(contigID);
                    }
                    catch (const out_of_range & e)
                    {
                        contigs++;
                        contig = new Contig();
                        contig->id = contigID;

                        contig->strain = strain;
                        if (regionsize.size() > 0)
                        {
                            contig->length = regionsize[items[0]];
                            regionsize.erase(items[0]);
                        }
                        strain->contigs.push_back(contig);
                        strain->contigpool[contigID] = contig;
                    }

                    if (regionsize.size() == 0 && gene->end > contig->length)
                        contig->length = gene->end;

                    contig->genes.push_back(gene);
                    gene->contig = contig;
                    
                    database->genepool[ strain->id + "|" + gene->id] = gene;
                    genes++;
                }
            }
            else if (type == TYPE_TAB)
            {
                if (line.empty())
                    continue;

                stringstream sstream(line);
                string s;

                vector<string> items;
                items.clear();

                while (getline(sstream, s, '\t'))
                {
                    items.push_back(s);
                    counter++;
                    //if (counter >= 4)
                    //    break;
                }
                /*
                if (counter != 10)
                {
                    printf("counter %i\n", counter);
                    printf("Error: Unexpected number of tab delimited items on a line inside: %s.tab\n", strainID.c_str());
                    exit(1);
                }
                */
                //if (items[4] == "OFF")
                //    continue;

                strain->CDS++;
                Gene *gene = new Gene;

                gene->id = items[0];
                //TODO parse product if available
                gene->product = items[8];
                gene->orientation = (items[1] == "1");
                gene->paralog = false;
                gene->start = stoi(items[2]);
                gene->end = stoi(items[3]);
                gene->length = abs(gene->start - gene->end)+1;
                gene->strain = strain;
                gene->left_neighbour = NULL;
                gene->right_neighbour = NULL;
                for (int i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
                    gene->neighbours[i] = NULL;
                gene->group = NULL;
                gene->score = 0.0f;
                //gene->coverage = 0;
                gene->GC3s = NAN;
                string contigID = items[5];
                Contig *contig;
                
                //
                try
                {
                    contig = strain->contigpool.at(contigID);
                }
                catch (const out_of_range & e)
                {
                    contigs++;
                    contig = new Contig();
                    contig->id = contigID;
                    
                    contig->strain = strain;
                    /*
                    if (regionsize.size() > 0)
                    {
                        contig->length = regionsize[items[0]];
                        regionsize.erase(items[0]);
                    }
                     */
                    strain->contigs.push_back(contig);
                    strain->contigpool[contigID] = contig;
                }
                if (gene->end > contig->length)
                    contig->length = gene->end;
                
                contig->genes.push_back(gene);
                gene->contig = contig;
                
                database->genepool[ strain->id + "|" + gene->id] = gene;
                genes++;
            }
        }
        
        if (type == TYPE_GFF)
        {
            if (line != "##FASTA")
            {
                input.close();
                //input = ifstream( INPUT_DIRECTORY + strain->id + ".fna");
            }
            
            if (input.is_open())
            {
                int A = 0;
                int T = 0;
                int C = 0;
                int G = 0;
                int N = 0;
                
                Contig *contig = NULL;
                string line;
                string contigsequence;
                while (getline(input, line))
                {
                    stringstream sstream(line);
                    
                    if (line.empty())
                        continue;
                    
                    //if (line[0] == '>')
                    //    continue;
                    
                    if (line[0] == '>')
                    {
                        if (contig != NULL)
                        {
                            for (auto gene = contig->genes.begin(); gene != contig->genes.end(); gene++ )
                            {
                                (*gene)->sequence = contigsequence.substr((*gene)->start-1, (*gene)->length );
                                if (!(*gene)->orientation)
                                {
                                    reverse((*gene)->sequence.begin(),(*gene)->sequence.end());
                                    for (int n = 0; n < (*gene)->sequence.size(); n++)
                                    {
                                        switch ((*gene)->sequence[n])
                                        {
                                            case 'A':
                                                (*gene)->sequence[n] = 'T';
                                                break;
                                            case 'T':
                                                (*gene)->sequence[n] = 'A';
                                                break;
                                            case 'C':
                                                (*gene)->sequence[n] = 'G';
                                                break;
                                            case 'G':
                                                (*gene)->sequence[n] = 'C';
                                                break;
                                        }
                                    }
                                }
                                
                                (*gene)->GC3s = (*gene)->CalculateGC3s();
                            }
                        }
                        string contigID = line.substr(1);
                        try
                        {
                            contig = strain->contigpool.at(contigID);
                        }
                        catch (const out_of_range & e)
                        {
                            //printf("Warning: Contig %s not found.\n", contigID.c_str());
                            contig = NULL;
                        }
                        contigsequence.clear();
                    }
                    else contigsequence += line;
                    A += count(line.begin(), line.end(), 'A');
                    T += count(line.begin(), line.end(), 'T');
                    C += count(line.begin(), line.end(), 'C');
                    G += count(line.begin(), line.end(), 'G');
                    N += count(line.begin(), line.end(), 'N');
                }
                //if (A+T+C+G+N != strain->bp )
                //    printf("Warning: In %s %i bp parsed from gff does not match %i bp parsed from fna. Difference is %i\n", strain->id.c_str(), strain->bp, A+T+C+G+N,strain->bp-(A+T+C+G));
                if (contig != NULL)
                {
                    for (auto gene = contig->genes.begin(); gene != contig->genes.end(); gene++ )
                    {
                        (*gene)->sequence = contigsequence.substr((*gene)->start-1, (*gene)->length );
                        if (!(*gene)->orientation)
                        {
                            reverse((*gene)->sequence.begin(),(*gene)->sequence.end());
                            for (int n = 0; n < (*gene)->sequence.size(); n++)
                            {
                                switch ((*gene)->sequence[n])
                                {
                                    case 'A':
                                        (*gene)->sequence[n] = 'T';
                                        break;
                                    case 'T':
                                        (*gene)->sequence[n] = 'A';
                                        break;
                                    case 'C':
                                        (*gene)->sequence[n] = 'G';
                                        break;
                                    case 'G':
                                        (*gene)->sequence[n] = 'C';
                                        break;
                                }
                            }
                        }
                        
                        (*gene)->GC3s = (*gene)->CalculateGC3s();
                    }
                }
                strain->GC = double(G+C)/double(A+T+C+G);
                strain->ubp = N;
                input.close();
            }
        }
        
        //for (auto c = strain->contigs.begin(); c != strain->contigs.end(); c++)
        for (auto contig = strain->contigs.begin(); contig != strain->contigs.end(); contig++)
        {
            if (TYPE_TAB)
                strain->bp += (*contig)->length;
            //Contig *contig = (*c).second;
            Gene *previousGene = (*contig)->genes.back();

            for (vector<Gene*>::iterator gene = (*contig)->genes.begin(); gene !=  (*contig)->genes.end(); gene++)
            {
                previousGene->right_neighbour = *gene;
                (*gene)->left_neighbour = previousGene;
                previousGene = *gene;
            }
            
            bool isContigLargeEnough;
            int sizeOfContigNeighbourhood;

            if ((*contig)->genes.size() < SIZE_OF_NEIGHTBOURHOOD + 1)
            {
                shortcontigs++;
                isContigLargeEnough = false;
                sizeOfContigNeighbourhood = (*contig)->genes.size()-1;
            }
            else
            {
                isContigLargeEnough = true;
                sizeOfContigNeighbourhood = SIZE_OF_NEIGHTBOURHOOD;
            }

            for (auto gene = (*contig)->genes.begin(); gene != (*contig)->genes.end(); gene++)
            {
                int n;
                Gene *g;
                g = *gene;

                if ((*gene)->orientation)
                {
                    //linear neighbourhood
                    for (n = 0; n < sizeOfContigNeighbourhood/2; n++)
                        g = g->right_neighbour;
                
                    for (n = 0; n < sizeOfContigNeighbourhood; n++)
                    {
                        if (g == *gene)
                            g = g->left_neighbour;
                    
                        (*gene)->neighbours[n+(SIZE_OF_NEIGHTBOURHOOD-sizeOfContigNeighbourhood)/2] = g;
                        g = g->left_neighbour;
                    }
                }
                else
                {
                    //linear neighbourhood
                    for (n = 0; n < sizeOfContigNeighbourhood/2; n++)
                        g = g->left_neighbour;

                    for (n = 0; n < sizeOfContigNeighbourhood; n++)
                    {
                        if (g == *gene)
                            g = g->right_neighbour;
                        
                        (*gene)->neighbours[n+(SIZE_OF_NEIGHTBOURHOOD-sizeOfContigNeighbourhood)/2] = g;
                        g = g->right_neighbour;
                    }
                }
            }
        }

        if (regionsize.size() > 0)
        {
            for (auto c = regionsize.begin(); c != regionsize.end(); c++)
            {
                Contig *contig = new Contig;
                contig->id = c->first;
                contig->length = c->second;
                contig->strain = strain;
                strain->contigs.push_back(contig);
                strain->contigpool[contig->id] = contig;
            }
        }
        printf("*");
    }
    printf("]\n");

    printf("Short contigs: %i\n", shortcontigs);
    printf("Parsing %i genes in %i contigs in %i strains took:\t%i seconds\n\n", genes, contigs, strains, (int)timer.Elapsed());

    return database->strains.size();
}
/*
void Parse::Scores(string INPUT_FILE, Dataset *dataset)
{
    int counter = 0;
    Timer timer;
    printf("Parsing scores from %s:\n", INPUT_FILE.c_str());
    
    ifstream input = ifstream( INPUT_FILE);
    if (!input.is_open())
    {
        printf("Error: Could not open file %s\n", INPUT_FILE.c_str() );
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

        if (line[0] != 'g' || line[1] != 'r' || line[2] != 'o' || line[3] != 'u' || line[4] != 'p')
            continue;

        string item;
        vector<string> items;
        while (getline(sstream, item, '\t'))
            items.push_back(item);

        string groupID = items[0];

        string::size_type loc;
        if ((loc = groupID.rfind('|')) != string::npos)
            groupID.erase(groupID.begin()+loc,groupID.end());

        Group *group = NULL;
        for (auto g = dataset->groups.begin(); g != dataset->groups.end(); g++)
        {
            if ((*g)->id == groupID)
            {
                group = *g;
                break;
            }
        }

        if (group == NULL)
        {
            printf("Group \"%s\" not found!\n", items[0].c_str());
            exit(1);
        }
        group->syntenyScoreSophisticated = stod( items[2]);
        group->syntenyScoreOld = stod( items[3]);
        group->syntenyScoreSimple = stod( items[4]);
        group->syntenyScoreAdjusted = stod( items[5]);
        group->syntenyScoreFast = stod( items[6]);
    }
    input.close();
    printf("Parsing scores took:\t%i seconds\n\n", (int)timer.Elapsed());
}
*/
void Parse::Sequences(string INPUT_DIRECTORY, Database *database)
{
    int genes = 0;
    int strains = 0;
    int contigs = 0;
    int shortcontigs = 0;
    Timer timer;
    DIR *dirp = opendir(INPUT_DIRECTORY.c_str());
    struct dirent *dp;
    
    timer.Start();
    printf("Parsing gene sequences from:\t%s\n", INPUT_DIRECTORY.c_str());
    printf("Processing [");
    
    while ((dp = readdir(dirp)) != NULL)
    {
        string strainID = string(dp->d_name);
        
        if (strainID.size() < 5)
            continue;
        
        if (strainID.compare(strainID.size()-4,4,".ffn"))
            continue;
        
        ifstream input = ifstream( INPUT_DIRECTORY + strainID);
        strainID.resize(strainID.size()-4);
        string::size_type loc;
        
        strains++;
        
        bool NCDS = false;
        string line;
        unordered_map<string, int> regionsize;
        
        Gene *gene = NULL;
        string sequence;
        while (getline(input, line))
        {
            if (line[0] == '>')
            {
                if (gene != NULL && !NCDS)
                {
                    if (!gene->sequence.empty() && gene->sequence != sequence)
                    {
                        printf("\nError: Newly parsed sequence for gene \"%s|%s\" does not match previous.",  gene->strain->id.c_str(), gene->id.c_str());
                        exit(1);
                    }
                    gene->sequence = sequence;
                    gene->GC3s = gene->CalculateGC3s();
                    //gene->sequence.clear();
                }
                NCDS = false;
                string geneID = line;
                string product;
                if ((loc = line.find_first_of(' ')) != string::npos)
                {
                    geneID = line.substr(1, loc-1);
                    product = line.substr(loc+1, line.size());
                }

                try
                {
                    gene = database->genepool.at(strainID + '|' + geneID);
                    genes++;
                }
                catch (const out_of_range & e)
                {
                    //if (product.find("RNA") == string::npos)
                    //    printf("Warning: Gene \"%s|%s\" not found!\n", strainID.c_str(), geneID.c_str());
                    NCDS = true;
                    //exit(1);
                }
                
                sequence.clear();
            }
            else sequence += line;
        }
        if (gene != NULL && !NCDS)
        {
            if (!gene->sequence.empty() && gene->sequence != sequence)
            {
                printf("\nError: Newly parsed sequence for gene \"%s|%s\" does not match previous.",  gene->strain->id.c_str(), gene->id.c_str());
                exit(1);
            }
            gene->sequence = sequence;
            gene->GC3s = gene->CalculateGC3s();
            //gene->sequence.clear();
        }

        input.close();
        printf("*");
    }
    printf("]\n");

    printf("Parsing %i gene sequences in %i strains took:\t%i seconds\n\n", genes, strains, (int)timer.Elapsed());
}
void Parse::ContigColours(string INPUT_FILE, Database *database)
{
    Timer timer;
    timer.Start();
    printf("Parsing contig colours from:\t%s\n", INPUT_FILE.c_str());

    ifstream input = ifstream( INPUT_FILE);
    if (!input.is_open())
    {
        printf("\nError: Could not open file %s\n", INPUT_FILE.c_str() );
        exit(1);
    }
    
    string line;
    int count = 0;
    while (getline(input, line))
    {
        stringstream sstream(line);
        
        if (line.empty() || ++count == 1)
            continue;
        
        string item;
        vector<string> items;
        while (getline(sstream, item, '\t'))
            items.push_back(item);
        
        database->contigcolours.push_back(make_pair(items[0] , items[1]));
    }
    input.close();
    printf("Parsing contig colours took:\t%i seconds\n\n", (int)timer.Elapsed());
}
/*
void Parse::Coverage(string INPUT_DIRECTORY, Database *database)
{
    Timer timer;
    DIR *dirp = opendir(INPUT_DIRECTORY.c_str());
    struct dirent *dp;
    set<Gene*> duplicatedgenes;

    timer.Start();
    printf("Parsing coverage data from directory %s:\n", INPUT_DIRECTORY.c_str());

    while ((dp = readdir(dirp)) != NULL)
    {
        string filename = string(dp->d_name);
        if (filename.size() < 10)
            continue;

        if (filename.compare(filename.size()-9, 9, ".coverage"))
            continue;

        ifstream input = ifstream( INPUT_DIRECTORY + filename);

        Strain *strain = NULL;

        string strainID = filename;
        strainID.resize(strainID.size()-9);
        string::size_type loc;
        if ((loc = strainID.rfind('-')) != string::npos)
            strainID.erase(strainID.begin()+loc,strainID.end());

        try
        {
            strain = database->strainpool.at(strainID);
        }
        catch (const out_of_range & e)
        {
            printf("Error: Strain %s not found!\n", strainID.c_str());
            exit(1);
        }
        vector<int> coverage;
        string contigID;
        string line;
        while (getline(input, line))
        {
            if (line == "")
                continue;
            
            stringstream sstream(line);
            string s;

            vector<string> items;
            items.clear();

            while (getline(sstream, s, '\t'))
            {
                if (s != "")
                    items.push_back(s);
            }

            if (contigID == items[0])
                coverage.push_back(stoi(items[2]));
            else
            {
                if (contigID != "")
                {
                    try
                    {
                        Contig *contig = strain->contigpool.at(contigID);
                        
                        for (auto gene = contig->genes.begin(); gene != contig->genes.end(); gene++)
                        {
                            (*gene)->coverage = 0.0f;

                            for (int i = (*gene)->start-1; i < (*gene)->end; i++ )
                                (*gene)->coverage += coverage[i];

                            (*gene)->coverage /= (double)(*gene)->length;
                        }
                    }
                    catch (const out_of_range & e)
                    {
                        printf("Error: Contig %s in strain %s not found!\n", contigID.c_str(), strain->id.c_str());
                        exit(1);
                    }
                }

                contigID = items[0];
                coverage.clear();
            }
        }
    }
    printf("Parsing coverage data took:\t%i seconds\n\n", (int)timer.Elapsed());
}
*/
void Parse::GroupColours(string INPUT_FILE, Dataset *dataset)
{
    int counter = 0;
    Timer timer;
    
    printf("Parsing gene group highlight colours from:\t%s\n", INPUT_FILE.c_str());

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
            printf("\nError: Gene group \"%s\" not found for colouring\n", items[0].c_str());
            continue;
            //exit(1);
        }
        
        for (auto  gene = group->genes.begin(); gene != group->genes.end(); gene++)
            (*gene)->colour = items[1];
    }
    printf("Parsing gene group highlight colours took:\t%i seconds\n\n", (int)timer.Elapsed());
}

/*
void Parse::Discordants(string INPUT_FILE, Dataset *dataset)
{
    int counter = 0;
    Timer timer;
    
    ifstream input = ifstream( INPUT_FILE);
    if (!input.is_open())
    {
        printf("Error: Could not open file %s\n", INPUT_FILE.c_str() );
        exit(1);
    }
    
    printf("Parsing discordant group info from %s:\n", INPUT_FILE.c_str());
    timer.Start();
    string line;
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
            printf("Error: Group \"%s\" not found \n", items[0].c_str());
            exit(1);
        }
        else
            group->discordant = true;
    }
    printf("Parsing discordant group info took:\t%i seconds\n\n", (int)timer.Elapsed());
}
*/
/*
 vector<Strain*>::iterator strain;
 vector<Gene*>::iterator gene;
 
 for (strain = strains->begin(); strain != strains->end(); strain++)
 {
 for (gene = (*strain)->genes.begin(); gene != (*strain)->genes.end(); gene++)
 {
 int n;
 Gene *g;
 g = *gene;
 
 //linear neighbourhood
 for (n = 0; n < SIZE_OF_NEIGHTBOURHOOD/2; n++)
 g = g->left_neighbour;
 
 for (n = 0; n < SIZE_OF_NEIGHTBOURHOOD; n++)
 {
 if (g == *gene)
 g = g->right_neighbour;
 
 (*gene)->neighbours[n] = g;
 g = g->right_neighbour;
 
 }

}
}

*/
/*
int Parse::fnaGroups0(string INPUT_DIRECTORY, Dataset *dataset)
{
    Timer timer;
    DIR *dirp = opendir(INPUT_DIRECTORY.c_str());
    struct dirent *dp;
    
    printf("Parsing directory %s:\n", INPUT_DIRECTORY.c_str());
    int genes = 0;
    timer.Start();
    while ((dp = readdir(dirp)) != NULL)
    {
        string groupname = string(dp->d_name);
        
        if (groupname.size() < 5 || groupname.compare(groupname.size()-4,4,".fna"))
            continue;
        
        groupname.resize(groupname.size()-4);
        
        ifstream input = ifstream( INPUT_DIRECTORY + dp->d_name);
        if (!input.is_open())
            printf("Warning: Could not open file %s\n", (INPUT_DIRECTORY + dp->d_name).c_str() );
        
        Group *group = NULL;
        string line1;
        string line2;
        while (getline(input, line1))
        {
            if (line1[0] == '>')
            {
                genes++;
                stringstream sstream(line1);
                string s;
                
                vector<string> items;
                items.clear();
                while (getline(sstream, s, '|'))
                    items.push_back(s);
                
                string::size_type loc;
                string strainID = items[1];
                string geneID = items[2];
                
                if ((loc = strainID.rfind('-')) != string::npos)
                    strainID.erase(strainID.begin()+loc,strainID.end());

                if ((loc = geneID.rfind('_')) != string::npos)
                    geneID.erase(geneID.begin(), geneID.begin()+loc+1);
                
                if (group == NULL)
                {
                    group = new Group;
                    group->id = groupname;
                    dataset->groups.push_back(group);
                }

                Strain *strain = NULL;
                try
                {
                    strain = dataset->database->strainpool.at(strainID);
                }
                catch (const out_of_range & e)
                {
                    strain = new Strain();
                    strain->id = strainID;
                    dataset->database->strainpool[strain->id] = strain;
                }

                Gene *gene = NULL;
                try
                {
                    gene = dataset->database->genepool.at(strainID + '|' + geneID);
                    gene->votes++;
                    
                    if (gene->votes > 1)
                    {
                        printf("Warning: Gene %s reference more than twice!\n", (strainID + '|' + geneID).c_str() );
                    }
                }
                catch (const out_of_range & e)
                {
                    gene = new Gene;
                    gene->id = geneID;
                    gene->strain = strain;
                    //gene->group = group;
                    gene->votes = 0;
                    gene->length = 0;
                    getline(input, line2);
                    for (auto sequence = line2.begin(); sequence != line2.end(); sequence++)
                        if (*sequence != '-')
                            gene->length++;
                    dataset->database->genepool[strainID + '|' + geneID] = gene;
                    strain->genes.push_back(gene);
                }
                group->genes.push_back(gene);
                //dataset->grouppool[strainID + '|' + geneID] = group;
                dataset->grouppool[gene] = group;
            }
        }
        
        input.close();
        
        if (group == NULL)
        {
            printf("Warning: No group formed from file %s\n", (INPUT_DIRECTORY + dp->d_name).c_str() );
        }
        else if (group->genes.empty())
        {
            printf("Warning: Group %s from file %s has no genes!\n", group->id.c_str(), (INPUT_DIRECTORY + dp->d_name).c_str() );
        }
    }
    closedir(dirp);
    dataset->genecount = genes;
    printf("genes:\t%d\n", dataset->genecount);
    printf("strains:\t%lu\n", dataset->database->strainpool.size());
    printf("groups:\t%lu\n", dataset->groups.size());
    //printf("ignored genes: %d\n", ignored);
    //printf("unique genes: %lu\n", database->genepool.size());
    //printf("included genes: %i\n", counter);
    //printf("duplicates: %i\n", duplicates);
    printf("Genes per group:\t%f\n", (double)dataset->database->genepool.size()/dataset->groups.size());
    printf("Genes per strain:\t%f\n", (double)dataset->database->genepool.size()/dataset->database->strainpool.size());
    //printf("Largest gene size: %d\n", longest);
    printf("Parsing took:\t%i seconds\n\n", (int)timer.Elapsed());
    return dataset->genecount;

}
*/
/*
int Parse::fnaGroups(string INPUT_DIRECTORY, Database *database)
{
    DIR *dirp = opendir(INPUT_DIRECTORY.c_str());
    struct dirent *dp;
    
    printf("Parsing directory %s:\n", INPUT_DIRECTORY.c_str());
    int duplicates = 0;
    int counter = 0;
    int ignored = 0;
    int longest = 0;
    
    while ((dp = readdir(dirp)) != NULL)
    {
        string groupname = string(dp->d_name);

        //if (!strcmp(dp->d_name, ".DS_Store") || !strcmp(dp->d_name, ".") || !strcmp(dp->d_name, ".."))
        //    continue;
        
        if (groupname.size() < 5 || groupname.compare(groupname.size()-4,4,".fna"))
            continue;
        
        groupname.resize(groupname.size()-4);
        
        //if (!groupname.compare("group10641"))
        //    printf("!\n");
        //groupname.erase(groupname.size()-4, groupname.back());
        
        //int counter = 0;

        //printf("%s\n", dp->d_name);
        ifstream input = ifstream( INPUT_DIRECTORY + dp->d_name);
        
        if (!input.is_open())
            printf("Warning: Could not open file %s\n", (INPUT_DIRECTORY + dp->d_name).c_str() );
        
        Group *group = NULL;
        string line1;
        string line2;
        while (getline(input, line1))
        {
            if (line1[0] == '>')
            {
                counter++;
                stringstream sstream(line1);
                string s;

                vector<string> items;
                items.clear();
                while (getline(sstream, s, '|'))
                    items.push_back(s);

                getline(input, line2);
                line2.erase(remove(line2.begin(), line2.end(), '-'), line2.end());
                
                if (group == NULL)
                {
                    try
                    {
                        group = database->grouppool2.at(line2);
                        //printf("Warning: identical gene in both %s and %s, ignoring the latter\n", group->id.c_str(), groupname.c_str());
                        group = NULL;
                        counter--;
                        ignored++;
                        continue;
                    }
                    catch (const out_of_range & e){}

                    group = new Group;
                    //int counter = 0;
                    //while (dp->d_name[++counter] != '\0'){}
                    //dp->d_name[counter-4] = '\0';
                    //group->id = dp->d_name;
                    group->id = groupname;
                    database->groups.push_back(group);
                }

                Strain *strain = NULL;
                try
                {
                    strain = database->strainpool.at(items[1]);
                }
                catch (const out_of_range & e)
                {
                    strain = new Strain();
                    strain->id = items[1];
                    database->strainpool[strain->id] = strain;
                }
                
                Gene *gene = NULL;
                try
                {
                    //database->grouppool2.at(line2);
                    gene = database->genepool.at(line2);
                    
                    if (gene->id.compare(line2))
                        printf("Hash mismatch:\n%s\n%s\n", gene->id.c_str(), line2.c_str());
                    gene->votes++;
                    gene->antiVotes++;
                    
                    //database->genepool[line2]->votes++;
                    duplicates++;
                }
                catch (const out_of_range & e)
                {
                    //test
                    try
                    {
                        database->grouppool2.at(line2);
                        printf("!!\n");
                        database->genepool.at(line2);
                        printf("!!!\n");
                    }
                    catch (const out_of_range & e) {}

                    gene = new Gene;
                    gene->id = line2;
                    gene->votes = 1;
                    gene->antiVotes = 1;
                    group->genes.push_back(gene);
                    database->genepool[gene->id] = gene;
                    database->grouppool2[gene->id] = group;
                    
                    if (longest < gene->id.size())
                        longest = gene->id.size();
                }
                strain->genes.push_back(gene);
                gene->strains.insert(strain);
            }
        }
        
        input.close();
        
        if (group == NULL)
        {
            //printf("Warning: No group formed from file %s\n", (INPUT_DIRECTORY + dp->d_name).c_str() );
        }
        else if (group->genes.empty())
        {
            printf("Warning: Group %s from file %s has no genes!\n", group->id.c_str(), (INPUT_DIRECTORY + dp->d_name).c_str() );
        }
    }
    closedir(dirp);
    
    printf("strains: %lu\n", database->strainpool.size());
    printf("groups: %lu\n", database->groups.size());
    printf("ignored genes: %d\n", ignored);
    printf("unique genes: %lu\n", database->genepool.size());
    printf("included genes: %i\n", counter);
    printf("duplicates: %i\n", duplicates);
    printf("Genes per group: %f\n", (double)counter/database->groups.size());
    printf("Genes per strain: %f\n", (double)counter/database->strainpool.size());
    printf("Largest gene size: %d\n", longest);

    if (duplicates+database->genepool.size() != counter)
        printf("Warning: Gene count does not add up!\n");
    
    return counter;
}

void Parse::fnaGroups2(string INPUT_DIRECTORY, Database *database)
{
    DIR *dirp = opendir(INPUT_DIRECTORY.c_str());
    struct dirent *dp;
    
    printf("Parsing directory %s:\n", INPUT_DIRECTORY.c_str());

    int coorthologues = 0;
    int totalgenes = 0;
    while ((dp = readdir(dirp)) != NULL)
    {
        string groupname = string(dp->d_name);
        
        if (groupname.size() < 5 || groupname.compare(groupname.size()-4,4,".fna"))
            continue;
        
        groupname.resize(groupname.size()-4);
        ifstream input = ifstream( INPUT_DIRECTORY + dp->d_name);
        
        Group *group = NULL;
        string line1;
        string line2;
        set<string> strainset;
        strainset.clear();
        while (getline(input, line1))
        {
            if (line1[0] == '>')
            {
                stringstream sstream(line1);
                string s;
                
                vector<string> items;
                items.clear();
                while (getline(sstream, s, '|'))
                    items.push_back(s);

                getline(input, line2);
                line2.erase(remove(line2.begin(), line2.end(), '-'), line2.end());

                if (group == NULL)
                {
                    group = new Group;
                    group->id = groupname;
                    database->groups.push_back(group);
                }

                Strain *strain = NULL;
                try
                {
                    strain = database->strainpool.at(items[1]);
                }
                catch (const out_of_range & e)
                {
                    strain = new Strain();
                    strain->id = items[1];
                    database->strainpool[strain->id] = strain;
                }

                Gene *gene = new Gene;
                gene->id = line2;
                group->genes.push_back(gene);
                strain->genes.push_back(gene);
                totalgenes++;

                if (strainset.find(line2) != strainset.end())
                    coorthologues++;
                else
                    strainset.insert(line2);
            }
        }
        
        input.close();
    }
    closedir(dirp);

    printf("number of genes: %d\n", totalgenes);
    printf("number of coorthologues: %d\n", coorthologues);
    printf("percent coorthologues: %f\n", (double) coorthologues / totalgenes);
    printf("number of coorthologues per group: %f\n", (double) coorthologues / database->groups.size());
    
    printf("strains: %lu\n", database->strainpool.size());
    printf("groups: %lu\n", database->groups.size());
    printf("ignored genes: %d\n", ignored);
    printf("unique genes: %lu\n", database->genepool.size());
    printf("actual genes: %i\n", counter);
    printf("duplicates: %i\n", duplicates);
    printf("Genes per group: %f\n", (double)counter/database->groups.size());
    printf("Genes per strain: %f\n", (double)counter/database->strainpool.size());
    printf("Largest gene size: %d\n", longest);
    
    if (duplicates+database->genepool.size() != counter)
        printf("Warning: Gene count does not add up!\n");
}
*/
