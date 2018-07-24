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
    bool GlobalSettings::AllowFreeTileMovement = false;
    bool GlobalSettings::TileDepositsEqualDifference = false;
    bool GlobalSettings::UseSingleStructuralMutationRate = false;
    bool GlobalSettings::MutationRatesEvolve = false;
    double GlobalSettings::BaseMetaMutationRate = 0.001;


    void GlobalSettings::Initialise(int numThreads)
    {
        auto randomDevice = std::random_device();
        RNG.seed(Randomise ? randomDevice() : FixedSeed);

        // Set number of OpenMP threads.
        NumThreads = numThreads;
        omp_set_num_threads(NumThreads);
    }
}
