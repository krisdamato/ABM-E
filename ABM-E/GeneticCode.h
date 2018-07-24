#pragma once

#include <limits>
#include "Helpers.h"

namespace ABME
{
    template <typename T>
    class GeneticCode
    {
    public:
        inline double GetFlipMutationRate() const
        {
            return double(FlipMutationRate) / TMax;
        }

        
        inline double GetInsertionMutationRate() const
        {
            return double(InsertionMutationRate) / TMax;
        }


        inline double GetDeletionMutationRate() const
        {
            return GlobalSettings::UseSingleStructuralMutationRate ? GetInsertionMutationRate() : double(DeletionMutationRate) / TMax;
        }


        inline double GetTransMutationRate() const
        {
            return double(TransMutationRate) / TMax;
        }


        inline double GetMetaMutationRate() const
        {
            return GlobalSettings::BaseMetaMutationRate + double(MetaMutationRate) / TMax;
        }

        
        inline T& GetFlipMutationParameter()
        {
            return FlipMutationRate;
        }

        
        inline T& GetInsertionMutationParameter()
        {
            return InsertionMutationRate;
        }
        
        
        inline T& GetDeletionMutationParameter()
        {
            return GlobalSettings::UseSingleStructuralMutationRate ? InsertionMutationRate : DeletionMutationRate;
        }
        
        
        inline T& GetTransMutationParameter()
        {
            return TransMutationRate;
        }
        
        
        inline T& GetMetaMutationParameter()
        {
            return MetaMutationRate;
        }
        
        
        inline size_t Length() const
        {
            return ActiveGenes.size();
        }

        Chromosome ActiveGenes;
        bool HasLargePatterns = false;
        static const int TMax = std::numeric_limits<T>::max();

    private:
        T FlipMutationRate = T(0.01 * TMax); // on T_MAX.
        T InsertionMutationRate = T(0.01 * TMax); // on T_MAX.
        T DeletionMutationRate = T(0.01 * TMax); // on T_MAX.
        T TransMutationRate = T(0.01 * TMax); // on T_MAX.
        T MetaMutationRate = T(0.005 * TMax); // on T_MAX.
    };
}