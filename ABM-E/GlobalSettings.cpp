#include "GlobalSettings.h"

#include <numeric>
#include <omp.h>

namespace ABME
{
    std::mt19937 GlobalSettings::RNG;
    const double GlobalSettings::GeneFlipMutationRate = 0.01; // Per gene.
    const double GlobalSettings::GeneticInsertionRate = 0.01; // Per reproduction.
    const double GlobalSettings::GeneticDeletionRate = 0.01; // Per reproduction.
    std::vector<int> GlobalSettings::GeneIndices;
    std::vector<int> GlobalSettings::BarcodeRandoms;
    int GlobalSettings::NextRandomIndex = 0;
    int GlobalSettings::NumThreads = 1;
    bool GlobalSettings::ForceEqualChromosomeReproductions = false;
    int GlobalSettings::DistanceStep = 1;
    bool GlobalSettings::AllowFreeTileMovement = false;


    void GlobalSettings::Initialise(int numThreads)
    {
        auto randomDevice = std::random_device();
        RNG.seed(Randomise ? randomDevice() : FixedSeed);

        // Chromosomes of up to length 2, choose randomly from the first
        // 2 possible genes. Chromosomes up to length 10 pick both of the
        // first two genes, and randomly from the set of next 8 genes. 
        // Finally, chromosomes up to length 522 pick all genes up to length
        // 10 and randomly from the rest.
        GeneIndices = std::vector<int>(NumGenes);
        std::iota(std::begin(GeneIndices), std::end(GeneIndices), 0);

        std::shuffle(GeneIndices.begin(), GeneIndices.begin() + 2, GlobalSettings::RNG);
        std::shuffle(GeneIndices.begin() + 2, GeneIndices.begin() + 10, GlobalSettings::RNG);
        std::shuffle(GeneIndices.begin() + 10, GeneIndices.begin() + 522, GlobalSettings::RNG);

        // Generate a lot of random numbers for use inside barcode updates.
        std::uniform_int_distribution<std::mt19937::result_type> dist(0, BarcodeSize - 1);
        for (auto i = 0; i < RandomsSize; ++i)
        {
            BarcodeRandoms.push_back(dist(RNG));
        }

        // Set number of OpenMP threads.
        NumThreads = numThreads;
        omp_set_num_threads(NumThreads);
    }
}
