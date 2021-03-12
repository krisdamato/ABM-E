#include "GlobalSettings.h"

#include <algorithm>
#include <numeric>
#include <omp.h>

namespace ABME
{
    std::mt19937 GlobalSettings::RNG;
    int GlobalSettings::NumThreads = 1;
    bool GlobalSettings::ForceEqualChromosomeReproductions = false;
    int GlobalSettings::DistanceStep = 1;
    int GlobalSettings::Seed = 0;
    bool GlobalSettings::UseSingleStructuralMutationRate = false;
    bool GlobalSettings::MutationRatesEvolve = false;
    double GlobalSettings::BaseMetaMutationRate = 1e-6;


    void GlobalSettings::Initialise(int numThreads)
    {
        auto randomDevice = std::random_device();
        if (Randomise) Seed = randomDevice();
        RNG.seed(Seed);

        // Set number of OpenMP threads.
        NumThreads = numThreads;
        omp_set_num_threads(NumThreads);
    }


    std::set<int> GlobalSettings::ShuffleIndices(int numIndices, bool useSimplerGenesFirst, int numGenes)
    {
        auto sizes = GetGenePatternCollectionSizes(numGenes);

        // Chromosomes of up to length 2, choose randomly from the first
        // 2 possible genes. Chromosomes up to length 10 pick both of the
        // first two genes, and randomly from the set of next 8 genes. 
        // Chromosomes up to length 522 pick all genes up to length
        // 10 and randomly from the rest, etc...
        std::set<int> newIndices;
        auto total = 0;
        auto begin = 0;
        if (useSimplerGenesFirst)
        {
            for (auto& s : sizes)
            {
                total += s;
                if (total < numIndices)
                {
                    for (int i = begin; i < total; ++i) newIndices.insert(i);
                }
                else
                {
                    std::uniform_int_distribution<std::mt19937::result_type> dist(begin, total - 1);

                    while (newIndices.size() < numIndices) newIndices.insert(dist(RNG));
                }

                begin += total;
                if (newIndices.size() == numIndices) break;
            }
        }
        else
        {
            std::uniform_int_distribution<std::mt19937::result_type> dist(0, numGenes - 1);

            while (newIndices.size() < numIndices) newIndices.insert(dist(RNG));
        }

        return newIndices;
    }


    std::vector<int> GlobalSettings::GetGenePatternCollectionSizes(int numGenes)
    {
        int total = 0;
        std::vector<int> numTiles{ 1, 3, 9, 25 };
        std::vector<int> sizes;

        for (auto n : numTiles)
        {
            int newSize = std::powl(2, n);
            total += newSize;
            if (total > numGenes) break;
            sizes.push_back(newSize);
        }

        return sizes;
    }
}
