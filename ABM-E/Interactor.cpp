#include "Interactor.h"

#include "Barcode.h"
#include "Individual.h"
#include "GeneticCode.h"

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
        auto firstHasLargePatterns = first.ItsGeneticCode.BehaviourGenes.HasLargePatterns;
        auto secondHasLargePatterns = second.ItsGeneticCode.BehaviourGenes.HasLargePatterns;

        for (auto i = 0; i < GlobalSettings::NumInteractionUpdates; ++i)
        {
            // Keep the complement of the intersection.
            firstCloneNext.Subtract(secondClone);
            secondCloneNext.Subtract(firstClone);

            // Update barcodes.
            firstCloneNext.Update(true, firstHasLargePatterns);
            secondCloneNext.Update(true, secondHasLargePatterns);

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

        // If any of the two are not yet of reproductive age, do nothing.
        if (first.Age < first.ItsGeneticCode.ReproductiveAge || second.Age < second.ItsGeneticCode.ReproductiveAge) return nullptr;

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
            // Crossover parameter mutation rate.
            newGeneticCode.SetFlipMutationParameter(Helpers::Crossover(firstGenetics.GetFlipMutationParameter(), secondGenetics.GetFlipMutationParameter(), dist));
            newGeneticCode.SetMetaMutationParameter(Helpers::Crossover(firstGenetics.GetMetaMutationParameter(), secondGenetics.GetMetaMutationParameter(), dist));

            // Mutate parameter rates.
            Helpers::BitFlip(firstGenetics.GetMetaMutationParameter(), dist, newGeneticCode.GetMetaMutationRate());
            Helpers::BitFlip(firstGenetics.GetFlipMutationParameter(), dist, newGeneticCode.GetMetaMutationRate());

            // Crossover reproductive age and programmed death.
            newGeneticCode.ReproductiveAge = Helpers::Crossover(firstGenetics.ReproductiveAge, secondGenetics.ReproductiveAge, dist);
            newGeneticCode.ProgrammedLifespan = Helpers::Crossover(firstGenetics.ProgrammedLifespan, secondGenetics.ProgrammedLifespan, dist);

            // Mutate reproductive age and programmed death.
            Helpers::BitFlip(newGeneticCode.ReproductiveAge, dist, newGeneticCode.GetFlipMutationRate());
            Helpers::BitFlip(newGeneticCode.ProgrammedLifespan, dist, newGeneticCode.GetFlipMutationRate());
        }
        
        // Recombine both chromosomes.
        newGeneticCode.BehaviourGenes = RecombineChromosomes(firstGenetics.BehaviourGenes, secondGenetics.BehaviourGenes, dist, newGeneticCode.GetMetaMutationRate());
        newGeneticCode.InteractionGenes = RecombineChromosomes(firstGenetics.InteractionGenes, secondGenetics.InteractionGenes, dist, newGeneticCode.GetMetaMutationRate());

        // Create an individual with this chromosome.
        auto offspring = new Individual(first.ItsEnvironment, newGeneticCode);
        offspring->X = first.X;
        offspring->Y = first.Y;

        // Perform the birth sequence.
        if (offspring->BeBorn())
        {
            return offspring;
        }
        
        // If there is no food tile to take, the offspring dies.
        delete offspring;
        return nullptr;
    }


    template <typename TChr>
    Chromosome<TChr> Interactor::RecombineChromosomes(Chromosome<TChr>& first, Chromosome<TChr>& second, std::uniform_real_distribution<> dist, double metaMutationRate)
    {
        Chromosome<TChr> newChromosome(first.MaxGeneValue + 1);

        if (GlobalSettings::MutationRatesEvolve)
        {
            // Crossover metamutation parameters.
            newChromosome.SetFlipMutationParameter(Helpers::Crossover(first.GetFlipMutationParameter(), second.GetFlipMutationParameter(), dist));
            newChromosome.SetInsertionMutationParameter(Helpers::Crossover(first.GetInsertionMutationParameter(), second.GetInsertionMutationParameter(), dist));
            if (!GlobalSettings::UseSingleStructuralMutationRate) newChromosome.SetDeletionMutationParameter(Helpers::Crossover(first.GetDeletionMutationParameter(), second.GetDeletionMutationParameter(), dist));
            newChromosome.SetTransMutationParameter(Helpers::Crossover(first.GetTransMutationParameter(), second.GetTransMutationParameter(), dist));

            // Mutate mutation parameters.
            Helpers::BitFlip(newChromosome.GetFlipMutationParameter(), dist, metaMutationRate);
            Helpers::BitFlip(newChromosome.GetInsertionMutationParameter(), dist, metaMutationRate);
            if (!GlobalSettings::UseSingleStructuralMutationRate) Helpers::BitFlip(newChromosome.GetDeletionMutationParameter(), dist, metaMutationRate);
            Helpers::BitFlip(newChromosome.GetTransMutationParameter(), dist, metaMutationRate);
        }

        // Pick a random length (from the two).
        std::uniform_int_distribution<std::mt19937::result_type> distLength(std::min(first.Genes.size(), second.Genes.size()), std::max(first.Genes.size(), second.Genes.size()));
        int newLength = distLength(GlobalSettings::RNG);

        // Convert chromosomes to vectors of gene indices and values.
        std::vector<int> geneIndices;
        std::vector<uchar> geneValues;
        Helpers::ConvertChromosomeToVectors(first.Genes, geneIndices, geneValues);
        Helpers::ConvertChromosomeToVectors(second.Genes, geneIndices, geneValues);

        std::uniform_int_distribution<std::mt19937::result_type> distIndex(0, first.Genes.size() + second.Genes.size() - 1);
        std::uniform_int_distribution<std::mt19937::result_type> distGeneIndex(0, GlobalSettings::NumGenes - 1);
        std::uniform_int_distribution<std::mt19937::result_type> distDeleteIndex(0, newLength - 1);
        std::uniform_int_distribution<std::mt19937::result_type> distGeneValue(0, newChromosome.MaxGeneValue);

        auto& newGenes = newChromosome.Genes;

        // Crossover active genes.
        // Pick a gene randomly from the two chromosomes, and ignore it if it already exists.
        for (auto i = 0; i < newLength;)
        {
            int index = distIndex(GlobalSettings::RNG);
            auto geneIndex = geneIndices[index];
            auto geneValue = geneValues[index];
            if (newGenes.count(geneIndex) == 0)
            {
                newGenes[geneIndex] = geneValue;
                if (geneIndex >= 522) newChromosome.HasLargePatterns = true;
                ++i;
            }
        }

        // Insert mutation.
        if ((dist(GlobalSettings::RNG) < newChromosome.GetInsertionMutationRate()) && (newLength < GlobalSettings::NumGenes))
        {
            auto geneIndex = -1;
            bool done = false;
            while (!done)
            {
                geneIndex = distGeneIndex(GlobalSettings::RNG);
                done = newGenes.count(geneIndex) == 0;
            }

            uchar geneValue = distGeneValue(GlobalSettings::RNG);
            newGenes[geneIndex] = geneValue;
            if (geneIndex >= 522) newChromosome.HasLargePatterns = true;
        }

        // Transmutation (replacement gene with new value).
        if ((dist(GlobalSettings::RNG) < newChromosome.GetTransMutationRate()) && (newLength < GlobalSettings::NumGenes))
        {
            // Pick a random gene and change its number.
            int changePosition = distDeleteIndex(GlobalSettings::RNG);

            // Remake a new chromosome.
            GeneSet transGenes;
            int i = 0;
            for (auto&[key, value] : newGenes)
            {
                if (i == changePosition)
                {
                    bool done = false;
                    auto geneIndex = -1;
                    while (!done)
                    {
                        geneIndex = distGeneIndex(GlobalSettings::RNG);
                        done = newGenes.count(geneIndex) == 0;
                    }
                    if (geneIndex >= 522) newChromosome.HasLargePatterns = true;
                    transGenes[geneIndex] = distGeneValue(GlobalSettings::RNG);
                }
                else
                {
                    transGenes[key] = value;
                }

                ++i;
            }

            newGenes = transGenes;
        }

        // Delete mutation.
        if ((dist(GlobalSettings::RNG) < newChromosome.GetDeletionMutationRate()) && (newLength >= 2))
        {
            // Remove a random gene.
            int removePosition = distDeleteIndex(GlobalSettings::RNG);
            auto it = newGenes.begin();
            for (auto i = 0; i < removePosition; ++i, ++it);
            newGenes.erase(it);
        }

        // Mutate.
        // Note: only gene value is mutated here.
        for (auto&[key, value] : newGenes)
        {
            bool flip = dist(GlobalSettings::RNG) < newChromosome.GetFlipMutationRate();
            if (flip ) 
            {
                if (newChromosome.MaxGeneValue == 1) value = 1 - value; // Simply flip
                else
                {
                    // Otherwise simple choose from the set of possibilities.
                    value = distGeneValue(GlobalSettings::RNG);
                }
            }
        }

        return newChromosome;
    }
}