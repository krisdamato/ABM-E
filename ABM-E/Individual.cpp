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


    bool Individual::AddDropFood(cv::Mat& environment, std::string& representation, int& difference)
    {
        // If it's a "sink"...
        if (difference > 0)
        {
            int toExtract = 1;
            CurrentBarcode->ExtractTiles(environment, X, Y, toExtract);
            if (toExtract > 0) return false;
        }
        // Or if it's a "source"...
        else if (difference < 0)
        {
            // Return one "food" tile to the environment.
            int toAdd = 1;
            CurrentBarcode->DropTiles(environment, X, Y, toAdd, true);
            if (toAdd > 0) ItsEnvironment.RegisterFoodAddition(toAdd);
        }

        return true;
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
        bool foundFood = true;
        if (difference != 0) foundFood = AddDropFood(baseEnvironment, representation, difference);

        // Perform movement.
        X += movement[0];
        Y += movement[1];
        X = std::max(0, std::min(interactableEnvironment.cols - GlobalSettings::BarcodeSize, X));
        Y = std::max(0, std::min(interactableEnvironment.rows - GlobalSettings::BarcodeSize, Y));

        // Only update food if we couldn't extract or deposit new tiles.
        LastCellsActive = cellsActive;

        // Update live status.
        // An individual dies if it has no food or all or no cell is active.
        if (cellsActive == std::pow(GlobalSettings::BarcodeSize - 2, 2) || cellsActive == 0 || !foundFood)
        {
            Alive = false;
            return;
        }

        // Add position to colocations.
        colocations.insert(Environment::ColocationMapType::value_type(Vec2i(X, Y), this));

        // Update individual parameters.
        Age += 1;
    }
}
