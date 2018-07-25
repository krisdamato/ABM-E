#pragma once

#include <vector>
#include "GeneticCode.h"

namespace ABME
{
    class Individual;

    /// Interacts multiple individuals for reproduction or death.
    class Interactor
    {
    public:
        static std::vector<Individual*> Interact(std::vector<Individual*> colocations);

    protected:
        static Individual* Interact(Individual& first, Individual& second);
        static Individual* Reproduce(Individual& first, Individual& second);
        
        template <typename T> static Chromosome<T> RecombineChromosomes(Chromosome<T>& first, Chromosome<T>& second, std::uniform_real_distribution<> dist, double metaMutationRate);
    };
}