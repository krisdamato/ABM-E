#pragma once

#include <random>
#include <set>

namespace ABME
{
    class GlobalSettings
    {
    public:
        static std::vector<int> GetGenePatternCollectionSizes(int numGenes);
        static void Initialise(int numThreads);
        static std::set<int> ShuffleIndices(int numIndices, bool useSimplerGenesFirst, int numGenes);

        static std::mt19937 RNG;
		static const int BarcodeSize = 16;
        //static const int NumBehaviourGenes = 33554432 + 522; // Includes 5x5 genes...
        static const int NumBehaviourGenes = 522; // Up to 3x3 genes
        static const bool Randomise = false;
        static const int CrisisPopulationSize = 500;
        static const int BehaviourGenePossibilities = 4;
        static bool ForceEqualChromosomeReproductions;
        static int Seed;
        static int DistanceStep;
        static bool MutationRatesEvolve;
        static bool UseSingleStructuralMutationRate;
        static double BaseMetaMutationRate;
        static double KillActiveMargin;
        static double FoodActiveMargin;
		static const int MaxPopulationSize = 600;
		static const int Patterns1x1 = 1;
		static const int Patterns3x1 = 8;
		static const int Patterns3x3 = 512;
		static const int Patterns5x5 = 33554432;
        static const int VitalityChangePerUpdate = 0;

    protected:
        static int NumThreads;
    };
}