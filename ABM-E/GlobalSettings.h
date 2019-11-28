#pragma once

#include <random>
#include <set>

namespace ABME
{
    class GlobalSettings
    {
    public:
        static std::vector<int> GetGenePatternCollectionSizes();
        static void Initialise(int numThreads);
        static std::set<int> ShuffleIndices(int numIndices, bool useSimplerGenesFirst);

        static std::mt19937 RNG;
        //static const int NumGenes = 33554432 + 522; // Includes 5x5 genes...
        static const int NumGenes = 522; // Up to 3x3 genes
        static const bool Randomise = true;
        static const int NumInteractionUpdates = 10;
        static const int BarcodeSize = 16;
        static const int CrisisPopulationSize = 500;
        static const int BehaviourGenePossibilities = 2;
        static const int InteractionGenePossibilities = 4;
        static const int MaxVitality = 256;
        static const double WorldUpdateProbability;
        static bool ForceEqualChromosomeReproductions;
        static int Seed;
        static int DistanceStep;
        static bool AllowFreeTileMovement;
        static bool TileDepositsEqualDifference;
        static bool MutationRatesEvolve;
        static bool UseSingleStructuralMutationRate;
        static double BaseMetaMutationRate;
		static const int MaxPopulationSize = 300;

    protected:
        static int NumThreads;
    };
}