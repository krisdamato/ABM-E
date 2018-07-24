#include "Interactor.h"

#include "Barcode.h"
#include "Individual.h"

namespace ABME
{
    /// This interacts each subsequent pair, in a non-overlapping way.
    /// Thus, if there are an odd number, one is uninteracted with.
    std::vector<Individual*> Interactor::Interact(std::vector<Individual*> colocations)
    {
        std::vector<Individual*> newIndividuals;
        for (auto i = 0; i < colocations.size() - 1; i += 2)
        {
            auto newIndividual = Interact(*colocations[i], *colocations[i + 1]);
            if (newIndividual != nullptr)
            {
                newIndividuals.push_back(newIndividual);
            }
        }

        return newIndividuals;
    }


    /// This decides the outcome of an interaction:
    /// Either both die, or one dies, or both live and produce an offspring.
    Individual* Interactor::Interact(Individual& first, Individual& second)
    {
        // Clone barcodes.
        auto firstClone = *first.CurrentBarcode;
        auto secondClone = *second.CurrentBarcode;
        auto firstCloneNext = *first.CurrentBarcode;
        auto secondCloneNext = *second.CurrentBarcode;
        auto firstCount = 0;
        auto secondCount = 0;

        for (auto i = 0; i < GlobalSettings::NumInteractionUpdates; ++i)
        {
            // Keep the complement of the intersection.
            firstCloneNext.Subtract(secondClone);
            secondCloneNext.Subtract(firstClone);

            // Update barcodes.
            firstCloneNext.Update(true, first.ItsGeneticCode.HasLargePatterns);
            secondCloneNext.Update(true, second.ItsGeneticCode.HasLargePatterns);

            // Replace barcodes of the next iteration.
            firstClone.SetStringRepresentation(firstCloneNext.GetStringRepresentation());
            secondClone.SetStringRepresentation(secondCloneNext.GetStringRepresentation());

            // Count the number of "live" cells in each.
            firstCount = firstClone.CountLiveCells();
            secondCount = secondClone.CountLiveCells();

            if (firstCount == 0 || firstCount == std::pow(GlobalSettings::BarcodeSize, 2) || secondCount == 0 || secondCount == std::pow(GlobalSettings::BarcodeSize, 2)) break;
        }

        // Kill any that have zero cells.
        if (firstCount == 0 || firstCount == std::pow(GlobalSettings::BarcodeSize, 2)) first.Kill();
        if (secondCount == 0 || secondCount == std::pow(GlobalSettings::BarcodeSize, 2)) second.Kill();

        // If chromosomes have to be equal length, check to make sure.
        if (GlobalSettings::ForceEqualChromosomeReproductions && first.ItsGeneticCode.Length() != second.ItsGeneticCode.Length()) return nullptr;

        // If both are still alive and they are genetically compatible, let's reproduce!
        // Otherwise nothing happens.
        if (first.IsAlive() && second.IsAlive())
        {
            cv::Mat& environment = first.ItsEnvironment.GetMap();
            
            return Reproduce(first, second);
        }

        return nullptr;
    }


    /// Reproduces by randomly selecting each gene from one of the individuals,
    /// inserts/deletes new genes, and applies mutation.
    Individual* Interactor::Reproduce(Individual& first, Individual& second)
    {
        std::uniform_real_distribution<> dist(0.0, 1.0);

        // Create an empty chromosome.
        auto& firstGenetics = first.ItsGeneticCode;
        auto& secondGenetics = second.ItsGeneticCode;
        GeneticCode<ushort> newGeneticCode;

        if (GlobalSettings::MutationRatesEvolve)
        {
            // Crossover metamutation parameters.
            newGeneticCode.GetFlipMutationParameter() = Helpers::Crossover(firstGenetics.GetFlipMutationParameter(), secondGenetics.GetFlipMutationParameter(), dist);
            newGeneticCode.GetInsertionMutationParameter() = Helpers::Crossover(firstGenetics.GetInsertionMutationParameter(), secondGenetics.GetInsertionMutationParameter(), dist);
            if (!GlobalSettings::UseSingleStructuralMutationRate) newGeneticCode.GetDeletionMutationParameter() = Helpers::Crossover(firstGenetics.GetDeletionMutationParameter(), secondGenetics.GetDeletionMutationParameter(), dist);
            newGeneticCode.GetTransMutationParameter() = Helpers::Crossover(firstGenetics.GetTransMutationParameter(), secondGenetics.GetTransMutationParameter(), dist);
            newGeneticCode.GetMetaMutationParameter() = Helpers::Crossover(firstGenetics.GetMetaMutationParameter(), secondGenetics.GetMetaMutationParameter(), dist);

            // Mutate mutation parameters.
            Helpers::BitFlip(newGeneticCode.GetMetaMutationParameter(), dist, newGeneticCode.GetMetaMutationRate());
            Helpers::BitFlip(newGeneticCode.GetFlipMutationParameter(), dist, newGeneticCode.GetMetaMutationRate());
            Helpers::BitFlip(newGeneticCode.GetInsertionMutationParameter(), dist, newGeneticCode.GetMetaMutationRate());
            if (!GlobalSettings::UseSingleStructuralMutationRate) Helpers::BitFlip(newGeneticCode.GetDeletionMutationParameter(), dist, newGeneticCode.GetMetaMutationRate());
            Helpers::BitFlip(newGeneticCode.GetTransMutationParameter(), dist, newGeneticCode.GetMetaMutationRate());
        }
        
        // Pick a random length (from the two).
        std::uniform_int_distribution<std::mt19937::result_type> distLength(std::min(firstGenetics.Length(), secondGenetics.Length()), std::max(firstGenetics.Length(), secondGenetics.Length()));
        int newLength = distLength(GlobalSettings::RNG);

        // Convert chromosomes to vectors of gene indices and values.
        std::vector<int> geneIndices;
        std::vector<uchar> geneValues;
        Helpers::ConvertChromosomeToVectors(firstGenetics.ActiveGenes, geneIndices, geneValues);
        Helpers::ConvertChromosomeToVectors(secondGenetics.ActiveGenes, geneIndices, geneValues);

        std::uniform_int_distribution<std::mt19937::result_type> distIndex(0, firstGenetics.Length() + secondGenetics.Length() - 1);
        std::uniform_int_distribution<std::mt19937::result_type> distGeneIndex(0, GlobalSettings::NumGenes - 1);
        std::uniform_int_distribution<std::mt19937::result_type> distDeleteIndex(0, newLength - 1);

        Chromosome& newChromosome = newGeneticCode.ActiveGenes;

        // Crossover active genes.
        // Pick a gene randomly from the two chromosomes, and ignore it if it already exists.
        for (auto i = 0; i < newLength;)
        {
            int index = distIndex(GlobalSettings::RNG);
            auto geneIndex = geneIndices[index];
            auto geneValue = geneValues[index];
            if (newChromosome.count(geneIndex) == 0)
            {
                newChromosome[geneIndex] = geneValue;
                if (geneIndex >= 522) newGeneticCode.HasLargePatterns = true;
                ++i;
            }
        }

        // Insert mutation.
        if ((dist(GlobalSettings::RNG) < newGeneticCode.GetInsertionMutationRate()) && (newGeneticCode.Length() < GlobalSettings::NumGenes))
        {
            auto geneIndex = -1;
            bool done = false;
            while (!done)
            {
                geneIndex = distGeneIndex(GlobalSettings::RNG);
                done = newChromosome.count(geneIndex) == 0;
            }
            
            uchar geneValue = dist(GlobalSettings::RNG) < 0.5 ? '1' : '0';
            newChromosome[geneIndex] = geneValue;
            if (geneIndex >= 522) newGeneticCode.HasLargePatterns = true;
        }

        // Transmutation (replacement gene with new value).
        if ((dist(GlobalSettings::RNG) < newGeneticCode.GetTransMutationRate()) && (newGeneticCode.Length() < GlobalSettings::NumGenes))
        {
            // Pick a random gene and change its number.
            int changePosition = distDeleteIndex(GlobalSettings::RNG);

            // Remake a new chromosome.
            Chromosome transChromosome;
            int i = 0;
            for (auto&[key, value] : newChromosome)
            {
                if (i == changePosition)
                {
                    bool done = false;
                    auto geneIndex = -1;
                    while (!done)
                    {
                        geneIndex = distGeneIndex(GlobalSettings::RNG);
                        done = newChromosome.count(geneIndex) == 0;
                    }
                    if (geneIndex >= 522) newGeneticCode.HasLargePatterns = true;
                    transChromosome[geneIndex] = dist(GlobalSettings::RNG) < 0.5 ? '1' : '0';
                }
                else
                {
                    transChromosome[key] = value;
                }

                ++i;
            }

            newChromosome = transChromosome;
        }

        // Delete mutation.
        if ((dist(GlobalSettings::RNG) < newGeneticCode.GetDeletionMutationRate()) && (newGeneticCode.Length() >= 2))
        {
            // Remove a random gene.
            int removePosition = distDeleteIndex(GlobalSettings::RNG);
            auto it = newChromosome.begin();
            for (auto i = 0; i < removePosition; ++i, ++it);
            newChromosome.erase(it);
        }

        // Mutate.
        // Note: only gene value is mutated here.
        for (auto&[key, value] : newChromosome)
        {
            bool flip = dist(GlobalSettings::RNG) < newGeneticCode.GetFlipMutationRate();
            if (flip)
            {
                if (value == '1') value = '0';
                if (value == '0') value = '1';
            }
        }

        // Create an individual with this chromosome.
        auto offspring = new Individual(first.ItsEnvironment, newGeneticCode);
        offspring->X = first.X;
        offspring->Y = first.Y;

        // Take a food tile from the environment and update balance.
        if (offspring->BeBorn())
        {
            return offspring;
        }
        
        // If there is no food tile to take, the offspring dies.
        delete offspring;
        return nullptr;
    }
}