//
//  Strain.cpp
//  Syntenizer3000
//
//  Created by Camous Moslemi on 23/04/2017.
//
//

#include "Strain.hpp"

Strain::Strain()
{
    id.clear();
    contigs.clear();
    bp = 0;
    ubp = 0;
    CDS = 0;
    tRNA = 0;
    tmRNA = 0;
    rRNA = 0;
    GC = NAN;
}

void Strain::ExportChartData(string OUTPUT_DIRECTORY, Database *db)
{
    ofstream out0( OUTPUT_DIRECTORY + id + "_contigs.csv", ifstream::out );
    out0 << "ID;Bound;Colour\n";

    for (auto contig = contigs.begin(); contig != contigs.end(); contig++)
    {
        string colour = "black";

        for (auto cmap = db->contigcolours.begin(); cmap != db->contigcolours.end(); cmap++ )
        {
            //printf("%s\n", cmap->first.c_str());
            if ((*contig)->id.find(cmap->first) != string::npos)
            {
                colour = cmap->second;
                break;
            }
        }

        out0 << (*contig)->id << ";0;" << colour << '\n';
        out0 << (*contig)->id << ';' << (*contig)->length << ';' << colour << '\n';
    }
    out0.close();

    ofstream out1( OUTPUT_DIRECTORY + id + "_lanes.csv", ifstream::out );
    //out1 << "ID;Location;Synteny;Abundance;Coverage;GC3s;Type\n";
    out1 << "ID;Location;Synteny;Abundance;GC3s;Colour\n";

    for (auto contig = contigs.begin(); contig != contigs.end(); contig++)
    {
        for (auto gene = (*contig)->genes.begin(); gene != (*contig)->genes.end(); gene++)
            out1 << (*contig)->id << ";"<< (*gene)->start + ((*gene)->length / 2) << ';' << ((*gene)->group != NULL ? (*gene)->group->SyntenizeAgainstGene((*gene))/(double)SIZE_OF_NEIGHTBOURHOOD : 0.0f) << ';' << ((*gene)->group != NULL ? (double)(*gene)->group->CountUniqueStrains() / (double)db->strains.size() : 0.0f / (double)db->strains.size()) << ';' /*<< (*gene)->coverage << ';'*/ << (*gene)->GC3s << ';' << ((*gene)->colour.empty() ? "NA" : (*gene)->colour) << '\n';
            //out1 << (*contig)->id << ";"<< (*gene)->start + ((*gene)->length / 2) << ';' << 0.0f << ';' << ((*gene)->group != NULL ? (double)(*gene)->group->CountUniqueStrains() / (double)db->strains.size() : -0.01f / (double)db->strains.size()) << ';' << (*gene)->coverage << ';' << (*gene)->GC3s << ';' << (*gene)->type << '\n';
    }
    out1.close();

    unordered_map<Group *, set<Gene *>> links;

    for (auto contig = contigs.begin(); contig != contigs.end(); contig++)
        for (auto gene = (*contig)->genes.begin(); gene != (*contig)->genes.end(); gene++)
            if ((*gene)->group != NULL)
                links[(*gene)->group].insert(*gene);

    ofstream out2( OUTPUT_DIRECTORY + id + "_paralogs.csv", ifstream::out );
    out2 << "ID1;Location1;ID2;Location2\n";

    for (auto link = links.begin(); link != links.end(); link++)
        if (link->second.size() > 1)
            for (auto g1 = link->second.begin(); g1 != link->second.end(); g1++)
                for (auto g2 = g1; g2 != link->second.end(); g2++)
                    if (*g1 != *g2)
                        out2 << (*g1)->contig->id << ";" << (*g1)->start + ((*g1)->length / 2)  << ";" << (*g2)->contig->id << ";" << (*g2)->start + ((*g2)->length / 2) << '\n';

    out2.close();
}
