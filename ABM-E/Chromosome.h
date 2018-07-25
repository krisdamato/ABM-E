#pragma once

#include <limits>
#include "Helpers.h"

namespace ABME
{
    template <typename TParam>
    class Chromosome
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

        
        inline TParam& GetFlipMutationParameter()
        {
            return FlipMutationRate;
        }

        
        inline TParam& GetInsertionMutationParameter()
        {
            return InsertionMutationRate;
        }
        
        
        inline TParam& GetDeletionMutationParameter()
        {
            return GlobalSettings::UseSingleStructuralMutationRate ? InsertionMutationRate : DeletionMutationRate;
        }
        
        
        inline TParam& GetTransMutationParameter()
        {
            return TransMutationRate;
        }


        inline void SetFlipMutationParameter(const TParam& param)
        {
            FlipMutationRate = param;
        }


        inline void SetInsertionMutationParameter(const TParam& param)
        {
            InsertionMutationRate = param;
        }


        inline void SetDeletionMutationParameter(const TParam& param)
        {
            if (GlobalSettings::UseSingleStructuralMutationRate) InsertionMutationRate = param;
            else DeletionMutationRate = param;
        }


        inline void SetTransMutationParameter(const TParam& param)
        {
            TransMutationRate = param;
        }


        GeneSet Genes;

        bool HasLargePatterns = false;

        static const int TMax = std::numeric_limits<TParam>::max();

    private:
        TParam FlipMutationRate = TParam(0.01 * TMax); // on T_MAX.
        TParam InsertionMutationRate = TParam(0.01 * TMax); // on T_MAX.
        TParam DeletionMutationRate = TParam(0.01 * TMax); // on T_MAX.
        TParam TransMutationRate = TParam(0.01 * TMax); // on T_MAX.
    };
}