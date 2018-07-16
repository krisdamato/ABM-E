#include "Individual.h"

#include <opencv2/highgui.hpp>
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


    bool Individual::AddDropFood(cv::Mat& environment, int numToTake)
    {
        // If it's a "sink"...
        if (numToTake > 0)
        {
            CurrentBarcode->ExtractTiles(environment, X, Y, numToTake);
            if (numToTake > 0) return false;
        }
        // Or if it's a "source"...
        else if (numToTake < 0)
        {
            // Return "food" tiles to the environment.
            CurrentBarcode->DropTiles(environment, X, Y, numToTake, true);
            if (numToTake < 0) ItsEnvironment.RegisterFoodAddition(-numToTake);
        }

        return true;
    }

    bool Individual::BeBorn()
    {
        if (AddDropFood(ItsEnvironment.GetMap(), 1))
        {
            ++Balance;
            return true;
        }

        return false;
    }


    Individual* Individual::Clone(bool ignoreBalance) const
    {
        auto* individual = new Individual(ItsEnvironment, ItsChromosome);
        individual->Age = Age;
        individual->X = X;
        individual->Y = Y;

        // Copy barcode pattern.
        individual->CurrentBarcode->SetStringRepresentation(CurrentBarcode->GetStringRepresentation());

        // Only copy balance if required. When capturing populations,
        // for instance, balance is best kept 0 so that the overall tile
        // number does not change.
        individual->Balance = ignoreBalance ? 0 : Balance;

        return individual;
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


    void Individual::Kill()
    {
        Alive = false;
        if (Balance > 0) AddDropFood(ItsEnvironment.GetMap(), -Balance);
        else ItsEnvironment.RegisterFoodAddition(Balance);
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
        
        if (difference != 0)
        {
            foundFood = AddDropFood(baseEnvironment, difference > 0 ? 1 : -1);

            // Update consumption balance.
            if (foundFood)
            {
                Balance = difference > 0 ? Balance + 1 : Balance - 1;
            }
        }

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
            Kill();
            return;
        }

        // Add position to colocations.
        colocations.insert(Environment::ColocationMapType::value_type(Vec2i(X, Y), this));

        // Update individual parameters.
        Age += 1;
    }
}
