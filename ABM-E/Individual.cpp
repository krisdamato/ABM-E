#include "Individual.h"

#include <opencv2/highgui.hpp>
#include <omp.h>
#include "Barcode.h"
#include "GlobalSettings.h"
#include "Helpers.h"

namespace ABME
{
    using namespace cv;

    PatternMap Individual::ShortGenePatternMap = Helpers::GenerateShortPatternMap();
    PatternMap Individual::LongGenePatternMap = Helpers::GenerateLongPatternMap();
	int Individual::nextID = 0;

    Individual::Individual(Environment& environment, GeneticCode<ushort> geneticCode) : ItsEnvironment(environment), ItsGeneticCode(geneticCode)
    {
        CurrentBarcode = std::make_unique<Barcode>(ItsGeneticCode.BehaviourGenes.Genes, GlobalSettings::BarcodeSize, GlobalSettings::BarcodeSize);
    }


    Individual::~Individual()
    {

    }


    void Individual::DrawBarcode(std::string& windowName)
    {
        CurrentBarcode->Draw(windowName);
    }


    bool Individual::BeBorn()
    {
        return true;
    }


    Individual* Individual::Clone(bool ignoreBalance) const
    {
        auto* individual = new Individual(ItsEnvironment, ItsGeneticCode);
        individual->Age = Age;
        individual->X = X;
        individual->Y = Y;

        // Copy barcode pattern.
        individual->CurrentBarcode->SetStringRepresentation(CurrentBarcode->GetStringRepresentation());

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
    }


    /// <summary>
    /// This assumes that world and individual updates are a function of the barcode state. Since there are many, many
    /// barcode states, we will really be collapsing them into more or less independent partitions that lead to different
    /// outcomes.
    /// </summary>
    void Individual::ProcessBarcode(Vec2i& movement)
    {
        // Calculate movement.
        CurrentBarcode->ComputeMetrics(movement);

        // Get active fraction.
        float activeFraction = (float) CurrentBarcode->CountLiveCells() / (GlobalSettings::BarcodeSize * GlobalSettings::BarcodeSize);

        // Map active fraction onto an individual function.
        if ((activeFraction < GlobalSettings::KillActiveMargin) || (activeFraction > (1 - GlobalSettings::KillActiveMargin)))
            Kill();
        else if ((activeFraction < GlobalSettings::FoodActiveMargin) || (activeFraction > (1 - GlobalSettings::FoodActiveMargin)))
        {
            auto cellsToTake = activeFraction < GlobalSettings::FoodActiveMargin ? 1 : -1;
            auto success = ChangeWorld(cellsToTake);
            if (!success) Kill();
        }
    }


    bool Individual::ChangeWorld(int cellsToTake) const
    {
        if (cellsToTake > 0)
            return CurrentBarcode->ExtractTiles(ItsEnvironment.GetMap(), X, Y, cellsToTake);
        else 
            return CurrentBarcode->DropTiles(ItsEnvironment.GetMap(), X, Y, -cellsToTake, true);
    }


    void Individual::Update(Mat& interactableEnvironment, Environment::ColocationMapType& colocations)
    {
        // Integrate environmental input.
        auto interactionRegion = interactableEnvironment(Rect(X, Y, GlobalSettings::BarcodeSize, GlobalSettings::BarcodeSize));
        CurrentBarcode->Input(interactionRegion);

        // Update barcode once.
        CurrentBarcode->Update(true, ItsGeneticCode.BehaviourGenes.HasLargePatterns, Vitality);

        // Update individual and world as a function of barcode.
        Vec2i movement;
        ProcessBarcode(movement);

		// Update individual parameters.
		Age++;
		
        if (Alive == false) return;

        // Clamp movement to avoid wall collisions.
        int newX = X + GlobalSettings::DistanceStep * movement[0];
        int newY = Y + GlobalSettings::DistanceStep * movement[1];
		ItsEnvironment.ClampPositions(newX, newY);

        const auto greaterMovement = std::max(std::abs(movement[0]), std::abs(movement[1]));
        if (greaterMovement > 0)
        {
			auto i = 1;
			for (; i <= greaterMovement; ++i)
			{
				newX = X + GlobalSettings::DistanceStep * (int(i * (float(movement[0]) / greaterMovement)));
				newY = Y + GlobalSettings::DistanceStep * (int(i * (float(movement[1]) / greaterMovement)));

				auto thisRect = cv::Rect(newX, newY, GlobalSettings::BarcodeSize, GlobalSettings::BarcodeSize);
				if (DetectCollision(thisRect, ItsEnvironment.GetRegions())) break;
			}

			--i;
            newX = X + GlobalSettings::DistanceStep * (int(i * (float(movement[0]) / greaterMovement)));
            newY = Y + GlobalSettings::DistanceStep * (int(i * (float(movement[1]) / greaterMovement)));
            ItsEnvironment.ClampPositions(newX, newY);
        }
        
        // Perform movement.
        X = newX;
        Y = newY;

        // Add position to colocations.
        colocations.insert(std::pair(Vec2i(X, Y), this));
    }
}
