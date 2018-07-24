#pragma once

#include <random>

namespace ABME
{
    class GlobalSettings
    {
    public:
        static void Initialise(int numThreads);

        static std::mt19937 RNG;
        static const int NumGenes = 33554432 + 522; // Includes 5x5 genes...
        static const int FixedSeed = 1;
        static const bool Randomise = false;
        static const int NumInteractionUpdates = 10;
        static const int BarcodeSize = 16;
        static const int CrisisPopulationSize = 500;
        static bool ForceEqualChromosomeReproductions;
        static int DistanceStep;
        static bool AllowFreeTileMovement;
        static bool TileDepositsEqualDifference;
        static bool MutationRatesEvolve;
        static bool UseSingleStructuralMutationRate;
        static double BaseMetaMutationRate;

    protected:
        static int NumThreads;
    };
}