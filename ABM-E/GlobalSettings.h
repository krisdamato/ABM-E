#pragma once

#include <random>

namespace ABME
{
    class GlobalSettings
    {
    public:
        static void Initialise(int numThreads);
        
        inline static int GetGeneIndex(int i)
        {
            return GeneIndices[i];
        }


        static std::mt19937 RNG;
        static const int NumGenes = 522;
        static const int FixedSeed = 3;
        static const bool Randomise = false;
        static const int NumInteractionUpdates = 10;
        static const int BarcodeSize = 16;
        static bool ForceEqualChromosomeReproductions;
        static int DistanceStep;
        static bool AllowFreeTileMovement;
        static bool TileDepositsEqualDifference;
        static bool MutationRatesEvolve;
        static double BaseMetaMutationRate;

    protected:
        static std::vector<int> GeneIndices;
        static int NumThreads;
        static int NextRandomIndex;
    };
}