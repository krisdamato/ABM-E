#include "Individual.h"

#include <opencv2\highgui.hpp>
#include "Barcode.h"
#include "GlobalSettings.h"

namespace ABME
{
    using namespace cv;

    Individual::Individual(Environment& environment, Chromosome chromosome) : ItsEnvironment(environment), ItsChromosome(chromosome)
    {
        CurrentBarcode = std::make_unique<Barcode>(ItsChromosome, GlobalSettings::BarcodeSize, GlobalSettings::BarcodeSize);
    }


    Individual::~Individual()
    {

    }


    void Individual::DrawBarcode(std::string& windowName)
    {
        CurrentBarcode->Draw(windowName);
    }


    /// Adds food tiles as the current active cell positions.
    void Individual::DropFood(cv::Mat& environment, int numFoodTiles)
    {
        Food -= numFoodTiles;
        CurrentBarcode->Output(environment, X, Y, numFoodTiles);
        
        // If some food tiles remain, let the environment to add them randomly.
        ItsEnvironment.RegisterFoodAddition(numFoodTiles);   
    }


    void Individual::ExtractFood(cv::Mat& environment, std::string& representation, int& numTiles)
    {
        // If it's a "sink"...
        if (numTiles > 0)
        {
            for (int i = 0; i < representation.size() && numTiles > 0; ++i)
            {
                if (representation[i] == '1')
                {
                    environment.at<uchar>(Y + i / GlobalSettings::BarcodeSize, X + i % GlobalSettings::BarcodeSize) = 0;
                    --numTiles;
                }
            }
        }
        // Or if it's a "source"...
        else if (numTiles < 0)
        {
            // Replenish food to starting point...
            auto diff = std::min(GlobalSettings::StartFood - Food, -numTiles);
            Food += diff;
            numTiles += diff;

            // ...and return the extra supply to the map.
            for (int i = 0; i < representation.size() && numTiles < 0; ++i)
            {
                if (representation[i] == '0')
                {
                    environment.at<uchar>(Y + i / GlobalSettings::BarcodeSize, X + i % GlobalSettings::BarcodeSize) = 255;
                    ++numTiles;
                }
            }
        }
    }


    const std::string& Individual::GetBarcodeString() const
    {
        return CurrentBarcode->GetStringRepresentation();
    }


    int Individual::InteractWithEnvironment(std::string& representation)
    {
        // Count current active cells.
        auto before = 0;
        for (uchar c : representation)
        {
            if (c == '1') ++before;
        }

        // Apply ruleset to region.
        std::string output(representation.size(), '0');
        CurrentBarcode->UpdateStringRepresentation(representation, output);

        // Count post-update active cells.
        auto after = 0;
        for (uchar c : output)
        {
            if (c == '1') ++after;
        }

        return after - before;
    }


    void Individual::Update(Mat& baseEnvironment, Mat& interactableEnvironment, Environment::ColocationMapType& colocations)
    {
        // Integrate environmental input.
        auto interactionRegion = interactableEnvironment(Rect(X, Y, GlobalSettings::BarcodeSize, GlobalSettings::BarcodeSize));
        auto baseRegion = baseEnvironment(Rect(X, Y, GlobalSettings::BarcodeSize, GlobalSettings::BarcodeSize));
        CurrentBarcode->Input(interactionRegion);

        // Update barcode once.
        CurrentBarcode->Update();

        // Calculate movement and consumption.
        Vec2i movement; int cellsActive = 0;
        CurrentBarcode->ComputeMetrics(movement, cellsActive);

        // Leave or extract some "food".
        auto representation = Helpers::ConvertMatToString(baseRegion);
        int difference = cellsActive - LastCellsActive;
        if (difference != 0) ExtractFood(baseEnvironment, representation, difference);

        // Perform movement.
        X += movement[0];
        Y += movement[1];
        X = std::max(0, std::min(interactableEnvironment.cols - GlobalSettings::BarcodeSize, X));
        Y = std::max(0, std::min(interactableEnvironment.rows - GlobalSettings::BarcodeSize, Y));

        // Only update food if we couldn't extract or deposit new tiles.
        Food -= difference;
        LastCellsActive = cellsActive;

        // Update live status.
        // An individual dies if it has no food or no cell is active (consumption = 0).
        if (Food == 0 || cellsActive == 0)
        {
            Alive = false;

            // The individual still leaves some cells in the environment.
            if (Food + cellsActive > 0) DropFood(baseEnvironment, Food + cellsActive);

            return;
        }

        // Add position to colocations.
        colocations.insert(Environment::ColocationMapType::value_type(Vec2i(X, Y), this));

        // Update individual parameters.
        Age += 1;
    }
}
