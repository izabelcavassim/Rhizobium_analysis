//
//  Score.cpp
//  Syntenizer3000
//
//  Created by Camous Moslemi on 19/09/2017.
//
//

#include "Score.hpp"

int Score::Fast(Gene *gene1, Gene *gene2)
{
    int score = 0;
    
    for (int i = (SIZE_OF_NEIGHTBOURHOOD/2)-1; i >= 0; i--)
        if (gene1->neighbours[i] == NULL || gene2->neighbours[i] == NULL || gene1->neighbours[i]->group == NULL || gene2->neighbours[i]->group == NULL || gene1->neighbours[i]->group != gene2->neighbours[i]->group)
            break;
        else
            score++;
    
    for (int i = (SIZE_OF_NEIGHTBOURHOOD/2); i < SIZE_OF_NEIGHTBOURHOOD; i++)
        if (gene1->neighbours[i] == NULL || gene2->neighbours[i] == NULL || gene1->neighbours[i]->group == NULL || gene2->neighbours[i]->group == NULL || gene1->neighbours[i]->group != gene2->neighbours[i]->group)
            break;
        else
            score++;
    
    return score;
}

double Score::Sophisticated(Gene *gene1, Gene *gene2)
{
    int        fU, fD, fL;
    char       ptr;
    int        i = 0, j = 0;
    double     score = 0.0f;
    int        d = 2;
    int        max;
    int        F[SIZE_OF_NEIGHTBOURHOOD+1][SIZE_OF_NEIGHTBOURHOOD+1];
    char       traceback[SIZE_OF_NEIGHTBOURHOOD+1][SIZE_OF_NEIGHTBOURHOOD+1];
    
    try
    {
        return gene1->scores.at(gene2);
    }
    catch (exception e){}

    try
    {
        return gene2->scores.at(gene1);
    }
    catch (exception e){}

    F[ 0 ][ 0 ] =  0 ;
    traceback[ 0 ][ 0 ] = 'n' ;
    
    for( j = 1; j <= SIZE_OF_NEIGHTBOURHOOD; j++ )
    {
        F[ 0 ][ j ] =  -j * d ;
        traceback[ 0 ][ j ] =  '-' ;
    }
    for( i = 1; i <= SIZE_OF_NEIGHTBOURHOOD; i++ )
    {
        F[ i ][ 0 ] =  -i * d ;
        traceback[ i ][ 0 ] =  '|' ;
    }
    
    for( i = 1; i <= SIZE_OF_NEIGHTBOURHOOD; i++ )
    {
        for( j = 1; j <= SIZE_OF_NEIGHTBOURHOOD; j++ )
        {
            fU = F[ i-1 ][ j ] - d ;
            fD = F[ i-1 ][ j-1 ] + (gene1->neighbours[ j-1 ] != NULL && gene2->neighbours[ i-1 ] != NULL && gene1->neighbours[ j-1 ]->group != NULL && gene2->neighbours[ i-1 ]->group != NULL && gene1->neighbours[ j-1 ]->group == gene2->neighbours[ i-1 ]->group ? 2 : -1);
            fL = F[ i ][ j-1 ] - d ;
            
            if( fU >= fD && fU >= fL )
            {
                max = fU ;
                ptr = '|' ;
            }
            else if( fD > fL )
            {
                max = fD ;
                ptr = '\\' ;
            }
            else
            {
                max = fL ;
                ptr = '-' ;
            }
            F[ i ][ j ] = max;
            traceback[ i ][ j ] = ptr ;
        }
    }
    i-- ; j-- ;
    
    while( i > 0 || j > 0 )
    {
        switch( traceback[ i ][ j ] )
        {
            case '|' :
                score -= 0.5f;
                i-- ;
                break ;
                
            case '\\':
                score += (gene1->neighbours[ j-1 ] != NULL && gene2->neighbours[ i-1 ] != NULL && gene1->neighbours[ j-1 ]->group != NULL && gene2->neighbours[ i-1 ]->group != NULL && gene1->neighbours[ j-1 ]->group == gene2->neighbours[ i-1 ]->group) ? 1 : 0;
                i-- ;  j-- ;
                break ;
                
            case '-' :
                score -= 0.5f;
                j-- ;
        }
    }

    gene1->scores[gene2] = score;
    return score;
}

double Score::Naive(Gene *gene1, Gene *gene2)
{
    double score = 0.0f;
    
    for (int i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
        for (int j = 0; j < SIZE_OF_NEIGHTBOURHOOD; j++)
            if (gene1->neighbours[i] != NULL && gene2->neighbours[j] != NULL && gene1->neighbours[i]->group != NULL && gene2->neighbours[j]->group != NULL && gene1->neighbours[i]->group == gene2->neighbours[j]->group)
                score++;
    
    return score;
}

double Score::Simple(Gene *gene1, Gene *gene2)
{
    double score = 0.0f;
    
    for (int i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
    {
        bool found = false;
        for (int j = 0; j < SIZE_OF_NEIGHTBOURHOOD; j++)
        {
            if (gene1->neighbours[i] != NULL && gene2->neighbours[j] != NULL && gene1->neighbours[i]->group != NULL && gene2->neighbours[j]->group != NULL && gene1->neighbours[i]->group == gene2->neighbours[j]->group)
            {
                found = true;
                break;
            }
        }
        if (found)
            score++;
    }
    
    return score;
}

double Score::Adjusted(Gene *gene1, Gene *gene2)
{
    int i, j;
    double matches = 0.0f;
    double factor = 0.0f;
    
    if (SIZE_OF_NEIGHTBOURHOOD == 100)
        factor = 3.7726f;
    
    if (SIZE_OF_NEIGHTBOURHOOD == 40)
        factor = 2.6517f;
    
    for (i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
    {
        double maxpartialmatch = 0.0;
        double partialmatch;
        
        for (j = 0; j < SIZE_OF_NEIGHTBOURHOOD; j++)
        {
            if (gene1->neighbours[i] != NULL && gene2->neighbours[j] != NULL && gene1->neighbours[i]->group != NULL && gene2->neighbours[j]->group != NULL && gene1->neighbours[i]->group == gene2->neighbours[j]->group)
            {
                partialmatch = 1.0f / (double)(1 + abs(i - j));
                
                if (partialmatch > maxpartialmatch)
                    maxpartialmatch = partialmatch;
            }
        }
        
        double distance = (i < (SIZE_OF_NEIGHTBOURHOOD / 2)) ? (SIZE_OF_NEIGHTBOURHOOD / 2) - i : i+1 - (SIZE_OF_NEIGHTBOURHOOD / 2);
        matches += (SIZE_OF_NEIGHTBOURHOOD/(distance*distance+SIZE_OF_NEIGHTBOURHOOD)) * factor * maxpartialmatch;
    }
    
    return matches;
}

double Score::Old(Gene *gene1, Gene *gene2)
{
    int i, j;
    double matches = 0.0f;
    
    for (i = 0; i < SIZE_OF_NEIGHTBOURHOOD; i++)
    {
        if (gene1->neighbours[i] != NULL &&  gene2->neighbours[i] != NULL && gene1->neighbours[i]->group != NULL &&  gene2->neighbours[i]->group != NULL && gene1->neighbours[i]->group == gene2->neighbours[i]->group)
            matches++;
        else
        {
            double maxpartialmatch = 0.0;
            double partialmatch;
            
            for (j = 0; j < SIZE_OF_NEIGHTBOURHOOD; j++)
            {
                if (gene1->neighbours[i] != NULL && gene2->neighbours[j] != NULL && gene1->neighbours[i]->group != NULL && gene2->neighbours[j]->group != NULL && gene1->neighbours[i]->group == gene2->neighbours[j]->group)
                {
                    partialmatch = 1.0f / double(abs(i-j)+1.0f);
                    
                    if (partialmatch > maxpartialmatch)
                        maxpartialmatch = partialmatch;
                }
            }
            
            matches += maxpartialmatch;
        }
    }
    return matches;
}

