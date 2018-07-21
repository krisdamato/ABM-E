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
            firstCloneNext.Update(true);
            secondCloneNext.Update(true);

            // Replace barcodes of the next iteration.
            firstClone.SetStringRepresentation(firstCloneNext.GetStringRepresentation());
            secondClone.SetStringRepresentation(secondCloneNext.GetStringRepresentation());

            // Count the number of "live" cells in each.
            firstCount = firstClone.CountLiveCells();
            secondCount = secondClone.CountLiveCells();

            if (firstCount == 0 || secondCount == 0) break;
        }

        // Kill any that have zero cells.
        if (firstCount == 0) first.Kill();
        if (secondCount == 0) second.Kill();

        // If chromosomes have to be equal length, check to make sure.
        if (GlobalSettings::ForceEqualChromosomeReproductions && first.ItsGeneticCode.Length() != second.ItsGeneticCode.Length()) return nullptr;

        // If both are still alive and they are genetically compatible, let's reproduce!
        // Otherwise nothing happens.
        if (firstCount > 0 && secondCount > 0)
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
        GeneticCode newGeneticCode;

        if (GlobalSettings::MutationRatesEvolve)
        {
            // Crossover metamutation parameters.
            newGeneticCode.FlipMutationRate = Helpers::Crossover(firstGenetics.FlipMutationRate, secondGenetics.FlipMutationRate, dist);
            newGeneticCode.InsertionMutationRate = Helpers::Crossover(firstGenetics.InsertionMutationRate, secondGenetics.InsertionMutationRate, dist);
            newGeneticCode.DeletionMutationRate = Helpers::Crossover(firstGenetics.DeletionMutationRate, secondGenetics.DeletionMutationRate, dist);
            newGeneticCode.MetaMutationRate = Helpers::Crossover(firstGenetics.MetaMutationRate, secondGenetics.MetaMutationRate, dist);

            // Mutate mutation parameters.
            newGeneticCode.MetaMutationRate = Helpers::BitFlip(newGeneticCode.MetaMutationRate, dist, newGeneticCode.GetMetaMutationRate());
            newGeneticCode.FlipMutationRate = Helpers::BitFlip(newGeneticCode.FlipMutationRate, dist, newGeneticCode.GetMetaMutationRate());
            newGeneticCode.InsertionMutationRate = Helpers::BitFlip(newGeneticCode.InsertionMutationRate, dist, newGeneticCode.GetMetaMutationRate());
            newGeneticCode.DeletionMutationRate = Helpers::BitFlip(newGeneticCode.DeletionMutationRate, dist, newGeneticCode.GetMetaMutationRate());
        }
        
        // Pick a random length (from the two).
        int newLength = dist(GlobalSettings::RNG) < 0.5 ? firstGenetics.Length() : secondGenetics.Length();

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