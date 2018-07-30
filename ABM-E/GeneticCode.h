#pragma once

#include "Chromosome.h"

namespace ABME
{
    template <typename TParam>
    class GeneticCode
    {
    public:
        inline size_t Length() const
        {
            return BehaviourGenes.Length() + InteractionGenes.Length();
        }


        inline TParam& GetFlipMutationParameter()
        {
            return FlipMutationRate;
        }


        inline TParam& GetMetaMutationParameter()
        {
            return MetaMutationRate;
        }


        inline double GetFlipMutationRate() const
        {
            return double(FlipMutationRate) / TMax;
        }


        inline double GetMetaMutationRate() const
        {
            return GlobalSettings::BaseMetaMutationRate + double(MetaMutationRate) / TMax;
        }


        inline void SetFlipMutationParameter(const TParam& param)
        {
            FlipMutationRate = param;
        }


        inline void SetMetaMutationParameter(const TParam& param)
        {
            MetaMutationRate = param;
        }


        Chromosome<TParam> BehaviourGenes = Chromosome<TParam>(GlobalSettings::BehaviourGenePossibilities); // Behaviour genes can have values in {0, 1}
        Chromosome<TParam> InteractionGenes = Chromosome<TParam>(GlobalSettings::InteractionGenePossibilities); // Interaction genes can have values in {0, ..., 3}

        TParam ProgrammedLifespan = 200;
        TParam ReproductiveAge = 0;

        static const int TMax = std::numeric_limits<TParam>::max();

    protected:
        TParam FlipMutationRate = TParam(0.01 * TMax); // on T_MAX.
        TParam MetaMutationRate = TParam(0.001 * TMax); // on T_MAX.
    };
}