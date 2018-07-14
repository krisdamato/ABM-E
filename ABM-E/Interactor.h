#pragma once

#include <vector>

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
        static Individual* Reproduce(Individual& first, Individual& second, int food);
        static Individual* ReproduceFixedGenes(Individual& first, Individual& second, int food);
    };
}