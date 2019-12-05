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
		static const int NumVitalityGenes = BarcodeSize * BarcodeSize;
        static const bool Randomise = true;
        static const int CrisisPopulationSize = 500;
        static const int BehaviourGenePossibilities = 2;
		static const int VitalityGenePossibilities = 2;
        static bool ForceEqualChromosomeReproductions;
        static int Seed;
        static int DistanceStep;
        static bool MutationRatesEvolve;
        static bool UseSingleStructuralMutationRate;
        static double BaseMetaMutationRate;
		static const int MaxPopulationSize = 300;

    protected:
        static int NumThreads;
    };
}