#include "Environment.h"

#include <algorithm>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <random>
#include "GlobalSettings.h"
#include "Helpers.h"
#include "Individual.h"
#include "Interactor.h"

namespace ABME
{
    using namespace cv;

    Environment::Environment(int width, int height, float probFood) : Map(height, width, CV_8UC1)
    {
        // Seed RNG.
        auto randomDevice = std::random_device();

        // Generate a map with "food" tiles.
        GenerateRandomFood(probFood);
    }


    int Environment::CountFoodCells() const
    {
        int count = 0;
        for (auto i = 0; i < Map.cols; ++i)
        {
            for (auto j = 0; j < Map.rows; ++j)
            {
                if (Map.at<uchar>(j, i) == 255) ++count;
            }
        }

        for (auto& ind : Individuals)
        {
            count += ind->Balance;
        }

        return count;
    }


    void Environment::Draw(std::string& windowName) const
    {
        // Convert from gayscale to color.
        Mat drawMap = cv::Mat(Map.rows, Map.cols, CV_8UC4);
        cv::cvtColor(Map, drawMap, cv::COLOR_GRAY2BGRA);

        // Get chromosome length range.
        size_t min = 522;
        size_t max = 0;
        for (auto& ind : Individuals)
        {
            min = std::min(min, ind->ItsChromosome.size());
            max = std::max(max, ind->ItsChromosome.size());
        }

        for (auto& ind : Individuals)
        {
            int redLevel = max > min ? 255 * float(ind->ItsChromosome.size() - min) / (max - min) : 128;
            int blueLevel = 255 - redLevel;
            rectangle(drawMap, Rect(ind->X, ind->Y, GlobalSettings::BarcodeSize, GlobalSettings::BarcodeSize), cv::Scalar(blueLevel, 0, redLevel, 64));
        }

        imshow(windowName, drawMap);
    }


    cv::Mat& Environment::GetMap()
    {
        return Map;
    }


    void Environment::Initialise(std::map<int, int> lengthCounts, bool useSameGeneIndices, bool useSimpleGenesFirst)
    {
        std::uniform_int_distribution<std::mt19937::result_type> distWidth(0, Map.cols - GlobalSettings::BarcodeSize);
        std::uniform_int_distribution<std::mt19937::result_type> distHeight(0, Map.rows - GlobalSettings::BarcodeSize);

        // Create individuals.
        for (auto&[length, count] : lengthCounts)
        {
            auto prototype = Helpers::GenerateRandomChromosome(length, useSimpleGenesFirst);
            for (auto i = 0; i < count; ++i)
            {
                auto chromosome = useSameGeneIndices ? Helpers::GenerateRandomChromosome(prototype) : Helpers::GenerateRandomChromosome(length, useSimpleGenesFirst);
                Individuals.push_back(std::make_unique<Individual>(*this, chromosome));
            }
        }

        // Assign random positions and set balances to 1 
        // so that every individual returns one food tile to the
        // environment.
        for (auto& individual : Individuals)
        {
            individual->X = distWidth(GlobalSettings::RNG);
            individual->Y = distHeight(GlobalSettings::RNG);
            individual->Balance++;
        }
    }


    void Environment::RegisterFoodAddition(int numTiles)
    {
        NumFoodCellsToAdd += numTiles;
    }


    void Environment::Update()
    {
        static int i = 0;
        static int born = 0;
        static int killed = 0;
        static int diedNaturally = 0;

        std::uniform_int_distribution<std::mt19937::result_type> dist(0, 2);

        // Take an additive snapshot of the world + barcodes.
        Snapshot = Map.clone();
        for (auto& individual : Individuals)
        {
            BurnBarcode(Snapshot, *individual);
        }

        // Update all individuals.
        for (auto& individual : Individuals)
        {
            individual->Update(Map, Snapshot, Colocations);
            
            // Add random ("Brownian") motion.
            individual->X += dist(GlobalSettings::RNG) - 1;
            individual->Y += dist(GlobalSettings::RNG) - 1;

            individual->X = std::max(0, std::min(individual->X, Map.cols - GlobalSettings::BarcodeSize));
            individual->Y = std::max(0, std::min(individual->Y, Map.rows - GlobalSettings::BarcodeSize));
        }

        // Remove dead individuals.
        for (auto it = Individuals.begin(); it != Individuals.end();)
        {
            if (!(*it)->IsAlive())
            {
                it = Individuals.erase(it);
                ++diedNaturally;
            }
            else
            {
                ++it;
            }
        }

        // Interact colocations.
        for (auto it = Colocations.begin(), end = Colocations.end(); it != end; it = Colocations.upper_bound(it->first))
        {
            auto& location = it->first;
            auto count = Colocations.count(location);
            if (count > 1)
            {
                auto individuals = Colocations.equal_range(location);
                std::vector<Individual*> colocated;
                for (auto it2 = individuals.first; it2 != individuals.second; ++it2)
                {
                    colocated.push_back(it2->second);
                }

                // Interact 'em!
                auto newIndividuals = Interactor::Interact(colocated);
                if (newIndividuals.size() > 0)
                {
                    std::uniform_int_distribution<std::mt19937::result_type> distPos(0, 8);

                    // Add the individuals to our list.
                    born += newIndividuals.size();
                    for (auto& individual : newIndividuals)
                    { 
                        individual->X += distPos(GlobalSettings::RNG) - 4;
                        individual->Y += distPos(GlobalSettings::RNG) - 4;

                        individual->X = std::max(0, std::min(individual->X, Map.cols - GlobalSettings::BarcodeSize));
                        individual->Y = std::max(0, std::min(individual->Y, Map.rows - GlobalSettings::BarcodeSize));

                        Individuals.push_back(std::unique_ptr<Individual>(individual));
                    }
                }
            }
        }

        // Remove more dead individuals.
        for (auto it = Individuals.begin(); it != Individuals.end();)
        {
            if (!(*it)->IsAlive())
            {
                it = Individuals.erase(it);
                ++killed;
            }
            else
            {
                ++it;
            }
        }

        // Generate missing food tiles from dead individuals.
        GenerateRandomFood(NumFoodCellsToAdd);
        NumFoodCellsToAdd = 0;

        if (i % 100 == 0)
        {
            std::cout.precision(3);
            std::cout << i << "] Num. individuals = " << Individuals.size() << "(" << born << " born this cycle, " << killed << " killed, " << diedNaturally << " died naturally)" << std::endl;
            std::map<int, int> genePool;
            for (auto& individual : Individuals)
            {
                genePool[individual->ItsChromosome.size()]++;
            }

            float averageLength = 0.f;
            for (auto&[length, count] : genePool)
            {
                std::cout << "Chromosome length " << length << ": " << count << " individuals.\n";
                averageLength += length * count;
            }
            if (Individuals.size() > 0) averageLength /= Individuals.size();
            std::cout << "\nAvg. chromosome length: " << averageLength << std::endl;

            // Report most popular chromosome.
            std::vector<Chromosome> chromosomes;
            for (auto& ind : Individuals) chromosomes.push_back(ind->ItsChromosome);

            std::pair<Chromosome, int> chrCount = Helpers::MostPopularChromosome(chromosomes);
            std::pair<Chromosome, int> chrTypeCount = Helpers::MostPopularChromosome(chromosomes, true);
            auto geneCountSet = Helpers::GeneStatistics(chromosomes);

            std::cout << "\nMost common chromosome type [" << chrTypeCount.second << "]: \n";
            std::cout << Helpers::ConvertChromosomeToString(chrTypeCount.first, true) << std::endl;
            std::cout << "Most common chromosome [" << chrCount.second << "]: \n";
            std::cout << Helpers::ConvertChromosomeToString(chrCount.first, false) << std::endl;
            std::cout << "Most common genes:\n";
            int j = 0;
            for (auto it = geneCountSet.begin(); it != geneCountSet.end() && j < 5; ++it, ++j)
            {
                std::cout << "[" << it->second << "] " << it->first.first << ":" << it->first.second << std::endl;
            }
            std::cout << std::endl;

            born = 0;
            killed = 0;
            diedNaturally = 0;
        }

        // Clear colocations.
        Colocations.clear();
        i++;
    }


    Individual& Environment::operator[](int index)
    {
        return *Individuals[index];
    }


    void Environment::GenerateRandomFood(float probFood)
    {
        std::uniform_real_distribution<> dist(0.0, 1.0);

        // Generate random "food" tiles.
        for (auto i = 0; i < Map.cols; ++i)
        {
            for (auto j = 0; j < Map.rows; ++j)
            {
                Map.at<uchar>(j, i) = (dist(GlobalSettings::RNG) < probFood) ? 255 : 0;
            }
        }
    }


    /// Generates (or removes) a specific number of foot tiles, depending
    /// on the sign of the argument.
    /// Warning! This will hang the application if it does not find an empty/food tile.
    void Environment::GenerateRandomFood(int numTilesToAdd)
    {
        std::uniform_int_distribution<std::mt19937::result_type> distWidth(0, Map.cols - 1);
        std::uniform_int_distribution<std::mt19937::result_type> distHeight(0, Map.rows - 1);

        while (numTilesToAdd > 0)
        {
            auto x = distWidth(GlobalSettings::RNG);
            auto y = distHeight(GlobalSettings::RNG);
            auto& tile = Map.at<uchar>(y, x);

            if (tile == 255) continue;
            else tile = 255;

            --numTilesToAdd;
        }

        while (numTilesToAdd < 0)
        {
            auto x = distWidth(GlobalSettings::RNG);
            auto y = distHeight(GlobalSettings::RNG);
            auto& tile = Map.at<uchar>(y, x);

            if (tile == 0) continue;
            else tile = 0;

            ++numTilesToAdd;
        }
    }

    void Environment::BurnBarcode(Mat& map, Individual& individual)
    {
        auto& barcode = individual.GetBarcodeString();
        for (auto i = 0; i < barcode.size(); ++i)
        {
            uchar c = barcode[i];
            if (c == '1')
            {
                int x = individual.X + (i % GlobalSettings::BarcodeSize);
                int y = individual.Y + (i / GlobalSettings::BarcodeSize);
                map.at<uchar>(y, x) = 255;
            }
        }
    }
}
