#include "Environment.h"

#include <algorithm>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <random>
#include "GlobalSettings.h"
#include "Helpers.h"
#include "Individual.h"
#include "Interactor.h"
#include "Logger.h"

namespace ABME
{
    using namespace cv;

    Environment::Environment(int width, int height) : Map(height, width, CV_8UC1)
    {
        Map = Scalar(0);

        // Seed RNG.
        auto randomDevice = std::random_device();
    }


    void Environment::AddRegion(cv::Rect region, float probability)
    {
        Regions.push_back(region);
        InitialRegionActiveTiles.push_back(probability * region.area());
        NumActiveTilesToAdd.push_back(0);
    }


    /// Copies the existing population for future release.
    void Environment::CapturePopulation()
    {
        // Clear current contents.
        Captured.clear();

        for (auto& ind : Individuals)
        {
            Captured.push_back(std::unique_ptr<Individual>(ind->Clone(true)));
        }

        PopulationCaptured = true;
    }


    /// Adds or takes away (clamped by the number of tiles that are free at this
    /// point in time) a number of tiles.
    // TODO: Needs fixing.
    int Environment::CauseTileCrisis(int numTilesToAdd)
    {
        //auto activeTiles = CountActiveTiles(false);
        //auto inactiveTiles = Map.cols * Map.rows - activeTiles;

        //// Clamp the number of tiles to add.
        //numTilesToAdd = std::min(inactiveTiles, std::max(numTilesToAdd, -activeTiles));

        //GenerateRandomFood(numTilesToAdd);

        return numTilesToAdd;
    }


    /// Clamps positions to environment maximum dimensions (and barcode margin).
    void Environment::ClampPositions(int& x, int& y) const
    {
        x = std::max(0, std::min(Map.cols - GlobalSettings::BarcodeSize, x));
        y = std::max(0, std::min(Map.rows - GlobalSettings::BarcodeSize, y));
    }


    /// Returns the number of active tiles in the map.
    /// If includeBoundBalance is true, it also adds the number of "tiles" 
    /// bound to the population.
    int Environment::CountActiveTiles(bool includeBoundBalance) const
    {
        int count = 0;
        for (auto i = 0; i < Map.cols; ++i)
        {
            for (auto j = 0; j < Map.rows; ++j)
            {
                if (Map.at<uchar>(j, i) == 255) ++count;
            }
        }

        if (includeBoundBalance)
        {
            for (auto& ind : Individuals)
            {
                for (auto b : ind->Balances) count += b;
            }
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


    std::vector<Rect>& Environment::GetRegions()
    {
        return Regions;
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
            auto newX = GlobalSettings::DistanceStep * (distWidth(GlobalSettings::RNG) / GlobalSettings::DistanceStep);
            auto newY = GlobalSettings::DistanceStep * (distHeight(GlobalSettings::RNG) / GlobalSettings::DistanceStep);

            while (Individual::DetectCollision(Rect(newX, newY, GlobalSettings::BarcodeSize, GlobalSettings::BarcodeSize), Regions))
            {
                newX = GlobalSettings::DistanceStep * (distWidth(GlobalSettings::RNG) / GlobalSettings::DistanceStep);
                newY = GlobalSettings::DistanceStep * (distHeight(GlobalSettings::RNG) / GlobalSettings::DistanceStep);
            }

            individual->X = newX;
            individual->Y = newY;
        }

        // Initialise active tiles.
        InitialiseTiles();
    }


    void Environment::RegisterActiveTileAddition(int regionIndex, int numTiles)
    {
        NumActiveTilesToAdd[regionIndex] += numTiles;
    }


    void Environment::ReleasePopulation()
    {
        // Copy the captured population onto the vector of individuals,
        // but maintain the list for future releases.
        for (auto& ind : Captured)
        {
            Individuals.push_back(std::unique_ptr<Individual>(ind->Clone(true)));
        }
    }


    void Environment::Update()
    {
        static int i = 0;
        static int born = 0;
        static int killed = 0;
        static int diedNaturally = 0;

        std::uniform_int_distribution<std::mt19937::result_type> dist(0, 2);
        auto& logger = Logger::Instance();

        // Take an additive snapshot of the world + barcodes.
        Snapshot = Map.clone();
        for (auto& individual : Individuals)
        {
            BurnBarcode(Snapshot, *individual);
        }

        // Update all individuals.
        for (auto& individual : Individuals)
        {
            individual->Update(Map, Snapshot, Colocations, Regions);
            
            // Add random ("Brownian") motion.
            int newX = individual->X + GlobalSettings::DistanceStep * ((int)dist(GlobalSettings::RNG) - 1);
            int newY = individual->Y + GlobalSettings::DistanceStep * ((int)dist(GlobalSettings::RNG) - 1);
            ClampPositions(newX, newY);

            while (Individual::DetectCollision(Rect(newX, newY, GlobalSettings::BarcodeSize, GlobalSettings::BarcodeSize), Regions))
            {
                newX = individual->X + GlobalSettings::DistanceStep * ((int)dist(GlobalSettings::RNG) - 1);
                newY = individual->Y + GlobalSettings::DistanceStep * ((int)dist(GlobalSettings::RNG) - 1);
                ClampPositions(newX, newY);
            }

            individual->X = newX;
            individual->Y = newY;
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
                        int newX = individual->X + GlobalSettings::DistanceStep * ((int)distPos(GlobalSettings::RNG) - 4);
                        int newY = individual->Y + GlobalSettings::DistanceStep * ((int)distPos(GlobalSettings::RNG) - 4);
                        ClampPositions(newX, newY);

                        while (Individual::DetectCollision(Rect(newX, newY, GlobalSettings::BarcodeSize, GlobalSettings::BarcodeSize), Regions))
                        {
                            newX = individual->X + GlobalSettings::DistanceStep * ((int)dist(GlobalSettings::RNG) - 1);
                            newY = individual->Y + GlobalSettings::DistanceStep * ((int)dist(GlobalSettings::RNG) - 1);
                            ClampPositions(newX, newY);
                        }

                        individual->X = newX;
                        individual->Y = newY;

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
        for (auto i = 0; i < NumActiveTilesToAdd.size(); ++i)
        {
            GenerateRandomFood(Regions[i], NumActiveTilesToAdd[i]);
            NumActiveTilesToAdd[i] = 0;
        }

        if (i % 100 == 0)
        {
            std::setprecision(4);

            std::stringstream log;
            log << i << "] Num. individuals = " << Individuals.size() << "(" << born << " born this cycle, " << killed << " killed, " << diedNaturally << " died naturally)" << std::endl;
            std::map<int, int> genePool;
            
            float averageAge = 0.f;
            for (auto& individual : Individuals)
            {
                genePool[individual->ItsChromosome.size()]++;
                averageAge += individual->Age;
            }
            averageAge /= Individuals.size();

            float averageLength = 0.f;
            for (auto&[length, count] : genePool)
            {
                log << "Chromosome length " << length << ": " << count << " individuals.\n";
                averageLength += length * count;
            }
            if (Individuals.size() > 0) averageLength /= Individuals.size();
            log << "\nAvg. chromosome length: " << averageLength << std::endl;
            log << "\nAvg. age: " << averageAge << std::endl;

            // Report most popular chromosome.
            std::vector<Chromosome> chromosomes;
            for (auto& ind : Individuals) chromosomes.push_back(ind->ItsChromosome);

            std::pair<Chromosome, int> chrCount = Helpers::MostPopularChromosome(chromosomes);
            std::pair<Chromosome, int> chrTypeCount = Helpers::MostPopularChromosome(chromosomes, true);
            auto geneCountSet = Helpers::GeneStatistics(chromosomes);

            log << "\nMost common chromosome type [" << chrTypeCount.second << "]: \n";
            log << Helpers::ConvertChromosomeToString(chrTypeCount.first, true) << std::endl;
            log << "Most common chromosome [" << chrCount.second << "]: \n";
            log << Helpers::ConvertChromosomeToString(chrCount.first, false) << std::endl;
            log << "Most common genes:\n";
            int j = 0;
            for (auto it = geneCountSet.begin(); it != geneCountSet.end() && j < 5; ++it, ++j)
            {
                log << "[" << it->second << "] " << it->first.first << ":" << it->first.second << std::endl;
            }
            log << std::endl;

            if (Captured.size() > 0)
            {
                log << "\nPopulation captured. Press r to release.\n";
            }

            born = 0;
            killed = 0;
            diedNaturally = 0;

            // Log to console and output file.
            logger << log.str();
        }

        // Clear colocations.
        Colocations.clear();
        i++;
    }


    Individual& Environment::operator[](int index)
    {
        return *Individuals[index];
    }


    void Environment::InitialiseTiles()
    {
        std::uniform_real_distribution<> dist(0.0, 1.0);

        for (auto i = 0; i < Regions.size(); ++i)
        {
            GenerateRandomFood(Regions[i], InitialRegionActiveTiles[i]);
        }
    }


    /// Generates (or removes) a specific number of foot tiles, depending
    /// on the sign of the argument.
    /// Warning! This will hang the application if it does not find an empty/food tile.
    void Environment::GenerateRandomFood(cv::Rect& region, int numTilesToAdd)
    {
        std::uniform_int_distribution<std::mt19937::result_type> distWidth(region.x, region.x + region.width - 1);
        std::uniform_int_distribution<std::mt19937::result_type> distHeight(region.y, region.y + region.height - 1);

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
