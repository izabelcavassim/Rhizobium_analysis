//
//  Gene.cpp
//  Syntenizer3000
//
//  Created by Camous Moslemi on 23/04/2017.
//
//

#include "Gene.hpp"

double Gene::Match(Group *group)
{
    double matchscore = 0.0f;
    int count = 0;
    
    for (auto gene = group->genes.begin(); gene != group->genes.end(); gene++)
    {
        if (*gene == this)
            continue;

        count++;
        matchscore += Score::Sophisticated(this, *gene);
    }
    matchscore /= (double)count;

    return(matchscore/(double) SIZE_OF_NEIGHTBOURHOOD);
}

double Gene::CalculateGC3s()
{
    double TA = 0;
    double GC = 0;
    for (int codon = 0; codon < sequence.length(); codon+=3)
    {
        if (sequence[codon] == 'T' && (sequence[codon+1] == 'G' || sequence[codon+1] == 'A') && (sequence[codon+2] == 'G' || sequence[codon+2] == 'A'))
            continue;

        if (sequence[codon] == 'A' && sequence[codon+1] == 'T' && sequence[codon+2] == 'G')
            continue;

        if (sequence[codon+2] == 'G' || sequence[codon+2] == 'C')
            GC++;
        if (sequence[codon+2] == 'T' || sequence[codon+2] == 'A')
            TA++;
    }

    return GC / (GC + TA);
}
