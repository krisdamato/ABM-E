#pragma once

#include "Helpers.h"

namespace ABME
{
    struct GeneticCode
    {
        inline double GetFlipMutationRate() const
        {
            return double(FlipMutationRate) / 255;
        }

        
        inline double GetInsertionMutationRate() const
        {
            return double(InsertionMutationRate) / 255;
        }


        inline double GetDeletionMutationRate() const
        {
            return double(DeletionMutationRate) / 255;
        }


        inline double GetTransMutationRate() const
        {
            return double(TransMutationRate) / 255;
        }


        inline double GetMetaMutationRate() const
        {
            return GlobalSettings::BaseMetaMutationRate + double(MetaMutationRate) / 255;
        }

        inline size_t Length() const
        {
            return ActiveGenes.size();
        }


        Chromosome ActiveGenes;
        uchar FlipMutationRate = 2; // on 255.
        uchar InsertionMutationRate = 2; // on 255.
        uchar DeletionMutationRate = 2; // on 255.
        uchar TransMutationRate = 2; // on 255.
        uchar MetaMutationRate = 1; // on 255.
    };
}