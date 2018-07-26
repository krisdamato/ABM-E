#include "Individual.h"

#include <opencv2/highgui.hpp>
#include "Barcode.h"
#include "GlobalSettings.h"

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
        individual->Vitality = Vitality;

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
        CurrentBarcode->Update(true, ItsGeneticCode.BehaviourGenes.HasLargePatterns);

        // Update world.
        Vitality += ProcessWorld();

        // Calculate movement and consumption.
        Vec2i movement; int cellsActive = 0;
        CurrentBarcode->ComputeMetrics(movement, cellsActive);

        // Update live status.
        // An individual dies if it has no food or all or no cell is active.
        if (cellsActive == std::pow(GlobalSettings::BarcodeSize - 2, 2) || cellsActive == 0 || Vitality <= 0 || Vitality >= 1000000 || Age >= ItsGeneticCode.ProgrammedLifespan)
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

        // Only update food if we couldn't extract or deposit new tiles.
        LastCellsActive = cellsActive;

        // Add position to colocations.
        colocations.insert(std::pair(Vec2i(X, Y), this));

        // Update individual parameters.
        Age += 1;
    }


    int Individual::ProcessWorld()
    {
        int vitalityUpdate = 0;

        // Convert world at this location to a string.
        auto& wholeMap = ItsEnvironment.GetMap();
        std::string worldString = Helpers::ConvertMatToString(wholeMap(Rect(X, Y, GlobalSettings::BarcodeSize, GlobalSettings::BarcodeSize)));

        // Find pattern matches in this region, and update the map with some small probability.
        auto oldWorldString = worldString;

        // Do the 1D genes first.
        for (auto&[key, val] : ItsGeneticCode.InteractionGenes.Genes)
        {
            if (key >= 10) break;
            std::string pattern = Helpers::GetParentPattern(key);
            vitalityUpdate += UpdateWorld1D(pattern, val, oldWorldString, worldString);
        }

        // Do the 3x3 2D genes next.
        vitalityUpdate += UpdateWorld2D(oldWorldString, 3, worldString);

        // Do the 5x5 2D genes next...
        if (ItsGeneticCode.InteractionGenes.HasLargePatterns) vitalityUpdate += UpdateWorld2D(oldWorldString, 5, worldString);

        // Update the world using the new string.
        for (int i = 0; i < worldString.size(); ++i)
        {
            wholeMap.at<uchar>(Y + i / GlobalSettings::BarcodeSize, X + i % GlobalSettings::BarcodeSize) = worldString[i] == 1 ? 255 : 0;
        }

        return vitalityUpdate;
    }


    int Individual::UpdateWorld1D(std::string& pattern, uchar& geneValue, std::string& oldWorldString, std::string& newWorldString)
    {
        int increment;
        uchar replacement;
        InterpretInteractionGeneValue(geneValue, increment, replacement);
        
        int count = 0;
        std::uniform_real_distribution<> dist(0.0, 1.0);

#pragma omp parallel for
        for (int j = 0; j < GlobalSettings::BarcodeSize; ++j)
        {
            std::string subString = oldWorldString.substr(j * GlobalSettings::BarcodeSize, GlobalSettings::BarcodeSize);
            std::vector<size_t> positions;
            positions.reserve(GlobalSettings::BarcodeSize);

            size_t pos = subString.find(pattern, 0);
            while (pos != std::string::npos)
            {
                count += increment;
                positions.push_back(pos);
                pos = subString.find(pattern, pos + 1);
            }

            // Replace in the new pattern.
            if (pattern.size() == 1)
            {
                for (auto& pos : positions)
                {
                    auto& tile = newWorldString[GlobalSettings::BarcodeSize * j + pos];
                    tile = dist(GlobalSettings::RNG) < GlobalSettings::WorldUpdateProbability ? replacement : tile;
                }
            }
            else if (pattern.size() == 3)
            {
                for (auto& pos : positions)
                {
                    auto& tile = newWorldString[GlobalSettings::BarcodeSize * j + pos + 1];
                    tile = dist(GlobalSettings::RNG) < GlobalSettings::WorldUpdateProbability ? replacement : tile;
                }
            }
        }

        return count;
    }


    int Individual::UpdateWorld2D(std::string& oldWorldString, int patternWidth, std::string& newWorldString)
    {
        auto& map = patternWidth <= 3 ? ShortGenePatternMap : LongGenePatternMap;
        std::uniform_real_distribution<> dist(0.0, 1.0);

        const int edgeLimit = patternWidth - 1;
        const int replaceOffset = (patternWidth - 1) / 2;
        int count = 0;

#pragma omp parallel for
        for (int j = 0; j < GlobalSettings::BarcodeSize - edgeLimit; ++j)
        {
            for (int i = 0; i < GlobalSettings::BarcodeSize - edgeLimit; ++i)
            {
                // Get the pattern at this position of the barcode.
                std::string subString;
                for (auto k = j; k < j + patternWidth; ++k)
                {
                    subString += oldWorldString.substr(k * GlobalSettings::BarcodeSize + i, patternWidth);
                }

                // Find which gene this would require.
                auto& geneIndex = map[subString];

                // Do we have this gene?
                if (ItsGeneticCode.InteractionGenes.Genes.count(geneIndex) == 0) continue;

                // Otherwise get the gene and increment the vitality update by the gene value.
                int increment;
                uchar replacement;
                InterpretInteractionGeneValue(ItsGeneticCode.InteractionGenes.Genes[geneIndex], increment, replacement);
                count += increment;

                // Replace at the right position.
                auto& tile = newWorldString[GlobalSettings::BarcodeSize * (j + replaceOffset) + i + replaceOffset];
                tile = dist(GlobalSettings::RNG) < GlobalSettings::WorldUpdateProbability ? replacement : tile;
            }
        }

        return count;
    }


    void Individual::InterpretInteractionGeneValue(uchar& val, int& vitalityUpdate, uchar& replacement)
    {
        vitalityUpdate = val / 2 == 0 ? -1 : 1;
        replacement = val % 2;
    }
}
