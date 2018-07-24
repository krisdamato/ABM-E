#include "Individual.h"

#include <opencv2/highgui.hpp>
#include "Barcode.h"
#include "GlobalSettings.h"

namespace ABME
{
    using namespace cv;

    Individual::Individual(Environment& environment, GeneticCode<ushort> geneticCode) : ItsEnvironment(environment), ItsGeneticCode(geneticCode)
    {
        CurrentBarcode = std::make_unique<Barcode>(ItsGeneticCode.ActiveGenes, GlobalSettings::BarcodeSize, GlobalSettings::BarcodeSize);
        Balances = std::vector<int>(ItsEnvironment.GetRegions().size());
    }


    Individual::~Individual()
    {

    }


    void Individual::DrawBarcode(std::string& windowName)
    {
        CurrentBarcode->Draw(windowName);
    }


    bool Individual::AddDropTile(int numToTake)
    {
        // If it's a "sink"...
        if (numToTake > 0)
        {
            CurrentBarcode->ExtractTiles(ItsEnvironment.GetMap(), X, Y, numToTake, ItsEnvironment.GetRegions(), Balances);
            if (numToTake > 0) return false;
        }
        // Or if it's a "source"...
        else if (numToTake < 0)
        {
            // Return "food" tiles to the environment.
            CurrentBarcode->DropTiles(ItsEnvironment.GetMap(), X, Y, numToTake, ItsEnvironment.GetRegions(), Balances, true);
        }

        return true;
    }


    bool Individual::BeBorn()
    {
        if (AddDropTile(1))
        {
            return true;
        }

        return false;
    }


    Individual* Individual::Clone(bool ignoreBalance) const
    {
        auto* individual = new Individual(ItsEnvironment, ItsGeneticCode);
        individual->Age = Age;
        individual->X = X;
        individual->Y = Y;

        // Copy barcode pattern.
        individual->CurrentBarcode->SetStringRepresentation(CurrentBarcode->GetStringRepresentation());

        // Only copy balance if required. When capturing populations,
        // for instance, balance is best kept 0 so that the overall tile
        // number does not change.
        if (!ignoreBalance) individual->Balances = Balances;

        return individual;
    }


    /// Detects whether any of the points of this individual
    /// overlaps with any of the unallowed regions.
    bool Individual::DetectCollision(const cv::Rect& thisRect, std::vector<cv::Rect>& regions)
    {
        const auto area = thisRect.area();
        auto runningArea = 0;
        for (auto& r : regions)
        {
            runningArea += (thisRect & r).area();
            if (runningArea == area) return false;
        }

        return true;
    }


    const std::string& Individual::GetBarcodeString() const
    {
        return CurrentBarcode->GetStringRepresentation();
    }


    void Individual::Kill()
    {
        Alive = false;
        if (GlobalSettings::AllowFreeTileMovement)
        {
            auto count = 0;

            // Sum all additions together.
            for (auto& c : Balances)
            {
                count += c;
            }

            // Drop the tiles where at the individual's active tiles.
            count = -count;
            if (count < 0) CurrentBarcode->DropTiles(ItsEnvironment.GetMap(), X, Y, count, ItsEnvironment.GetRegions(), Balances, true);

            // If there still remain tiles, drop them around the active tiles at the body.
            if (count < 0) CurrentBarcode->DropTiles(ItsEnvironment.GetMap(), X, Y, count, ItsEnvironment.GetRegions(), Balances, false);
            
            // If at this point, there are still some more to add (or take), drop them in the region.
            if (count != 0)
            {
                int index = Helpers::PointInRegionIndex(Point(X, Y), ItsEnvironment.GetRegions());
                ItsEnvironment.RegisterActiveTileAddition(index, -count);
            }
        }
        else
        {
            for (int i = 0; i < Balances.size(); ++i)
            {
                ItsEnvironment.RegisterActiveTileAddition(i, Balances[i]);
            }
        }
    }


    void Individual::Update(Mat& interactableEnvironment, Environment::ColocationMapType& colocations, std::vector<cv::Rect>& offLimitRegions)
    {
        // Integrate environmental input.
        auto interactionRegion = interactableEnvironment(Rect(X, Y, GlobalSettings::BarcodeSize, GlobalSettings::BarcodeSize));
        CurrentBarcode->Input(interactionRegion);

        // Update barcode once.
        CurrentBarcode->Update(true, ItsGeneticCode.HasLargePatterns);

        // Calculate movement and consumption.
        Vec2i movement; int cellsActive = 0;
        CurrentBarcode->ComputeMetrics(movement, cellsActive);

        // Leave or extract some "food".
        int difference = cellsActive - LastCellsActive;
        bool foundTile = true;
        
        if (difference != 0)
        {
            auto numDeposit = GlobalSettings::TileDepositsEqualDifference ? difference : (difference > 0 ? 1 : -1);
            foundTile = AddDropTile(numDeposit);
        }

        // Update live status.
        // An individual dies if it has no food or all or no cell is active.
        if (cellsActive == std::pow(GlobalSettings::BarcodeSize - 2, 2) || cellsActive == 0 || !foundTile)
        {
            Kill();
            return;
        }

        // Compute movement and collisions.
        int newX = X + GlobalSettings::DistanceStep * (movement[0] / GlobalSettings::DistanceStep);
        int newY = Y + GlobalSettings::DistanceStep * (movement[1] / GlobalSettings::DistanceStep);
        ItsEnvironment.ClampPositions(newX, newY);

        const auto greaterMovement = std::max(std::abs(movement[0]), std::abs(movement[1]));
        if (greaterMovement > 0)
        {
            for (auto i = greaterMovement; i > 0;)
            {
                auto thisRect = cv::Rect(newX, newY, GlobalSettings::BarcodeSize, GlobalSettings::BarcodeSize);
                if (!DetectCollision(thisRect, offLimitRegions)) break;

                --i;
                newX = X + GlobalSettings::DistanceStep * (int(i * (float(movement[0]) / greaterMovement) / GlobalSettings::DistanceStep));
                newY = Y + GlobalSettings::DistanceStep * (int(i * (float(movement[1]) / greaterMovement) / GlobalSettings::DistanceStep));
                ItsEnvironment.ClampPositions(newX, newY);
            }
        }
        
        // Perform movement.
        X = newX;
        Y = newY;

        // Only update food if we couldn't extract or deposit new tiles.
        LastCellsActive = cellsActive;

        // Add position to colocations.
        colocations.insert(Environment::ColocationMapType::value_type(Vec2i(X, Y), this));

        // Update individual parameters.
        Age += 1;
    }
}
