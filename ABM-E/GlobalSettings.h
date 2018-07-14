#pragma once

#include <random>

namespace ABME
{
    class GlobalSettings
    {
    public:
        static void Initialise();
        
        inline static int GetGeneIndex(int i)
        {
            return GeneIndices[i];
        }

        inline static int GetBarcodeRandom()
        {
            auto rand = BarcodeRandoms[NextRandomIndex++];
            if (NextRandomIndex >= RandomsSize) NextRandomIndex = 0;

            return rand;
        }

        static std::mt19937 RNG;
        static const double GeneFlipMutationRate;
        static const double GeneticInsertionRate;
        static const double GeneticDeletionRate;
        static const int NumGenes = 522;
        static const int FixedSeed = 7;
        static const bool Randomise = true;
        static const int NumInteractionUpdates = 10;
        static const int BarcodeSize = 16;
        static const int StartFood = BarcodeSize * BarcodeSize;
        static const int OffspringCost = StartFood / 3;

    protected:
        static std::vector<int> GeneIndices;
        static std::vector<int> BarcodeRandoms;
        static const int RandomsSize = 200000;
        static int NextRandomIndex;
    };
}