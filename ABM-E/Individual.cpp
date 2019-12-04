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


    void Individual::Update(Mat& interactableEnvironment, Environment::ColocationMapType& colocations)
    {
        // Integrate environmental input.
        auto interactionRegion = interactableEnvironment(Rect(X, Y, GlobalSettings::BarcodeSize, GlobalSettings::BarcodeSize));
        CurrentBarcode->Input(interactionRegion);

        // Update barcode once.
        CurrentBarcode->Update(true, ItsGeneticCode.BehaviourGenes.HasLargePatterns, Vitality, ItsGeneticCode.VitalityGenes.Genes);

		// Update individual parameters.
		Age++;
		Vitality--;
		if (Vitality <= 0)
		{
			Alive = false;
			return;
		}

        // Calculate movement and consumption.
        Vec2i movement;
        CurrentBarcode->ComputeMetrics(movement);

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
                if (!DetectCollision(thisRect, ItsEnvironment.GetRegions())) break;

                --i;
                newX = X + GlobalSettings::DistanceStep * (int(i * (float(movement[0]) / greaterMovement) / GlobalSettings::DistanceStep));
                newY = Y + GlobalSettings::DistanceStep * (int(i * (float(movement[1]) / greaterMovement) / GlobalSettings::DistanceStep));
                ItsEnvironment.ClampPositions(newX, newY);
            }
        }
        
        // Perform movement.
        X = newX;
        Y = newY;

        // Add position to colocations.
        colocations.insert(std::pair(Vec2i(X, Y), this));
    }
}
