#include "GlobalSettings.h"

#include <algorithm>
#include <numeric>
#include <omp.h>

namespace ABME
{
    std::mt19937 GlobalSettings::RNG;
    std::vector<int> GlobalSettings::GeneIndices;
    int GlobalSettings::NextRandomIndex = 0;
    int GlobalSettings::NumThreads = 1;
    bool GlobalSettings::ForceEqualChromosomeReproductions = false;
    int GlobalSettings::DistanceStep = 1;
    bool GlobalSettings::AllowFreeTileMovement = false;
    bool GlobalSettings::TileDepositsEqualDifference = false;
    bool GlobalSettings::MutationRatesEvolve = false;
    double GlobalSettings::BaseMetaMutationRate = 0.001;


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

        // Set number of OpenMP threads.
        NumThreads = numThreads;
        omp_set_num_threads(NumThreads);
    }
}
