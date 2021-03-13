#include "Barcode.h"

#include <omp.h>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <vector>
#include "Individual.h"

namespace ABME
{
    using namespace cv;


    Barcode::Barcode(GeneSet& behaviourGenes, int width, int height) : behaviourGenes(behaviourGenes), width(width), height(height)
    {
        // Create a string to represent the barcode with the proper length.
        barcode = std::string(width * height, 0);
    }


    Barcode::Barcode(const Barcode& rhs) : behaviourGenes(rhs.behaviourGenes), barcode(rhs.barcode), width(rhs.width), height(rhs.height)
    {
        
    }


    /// Computes movement from current barcode, assuming that
    /// barcode pixels are alternating in the directions (+i, +j, -i, -j).
    void Barcode::ComputeMetrics(Vec2i& movement) const
    {
        float positiveX = 0.f;
        float positiveY = 0.f;

		auto k = 0;
        for (auto i = 0; i < barcode.size(); ++i)
        {
            // Ignore barcode elements on the left/right boundaries.
            auto x = i % width;
			auto y = i / width;
            if (x == 0 || x == width - 1 || y == 0 || y == height - 1) continue;

            // Note: boundary cells are influenced by border effects.
            if (barcode[i] != 0)
            {
                // Normalise movement increments so that updates are fair 
                // with respect to cells available. 
                switch (k % 4)
                {
                case 0:
                    positiveX += 0.5f;
                    break;
                case 1:
                    positiveY += 0.5f;
                    break;
                case 2:
                    positiveX -= 0.5f;
                    break;
                case 3:
                    positiveY -= 0.5f;
                    break;
                }
            }

			++k;
        }

        // Update motion.
        movement[0] = int(positiveX);
        movement[1] = int(positiveY);
    }


    int Barcode::CountLiveCells() const
    {
        auto count = 0;
        for (uchar c : barcode)
        {
            if (c == 1) ++count;
        }

        return count;
    }


    /// Debug draw the current barcode.
    void Barcode::Draw(std::string& windowName) const
    {
        Mat image(height * CellSize, width * CellSize, CV_8UC1);

        for (auto i = 0; i < height; ++i)
        {
            for (auto j = 0; j < width; ++j)
            {
                uchar c = barcode[i * width + j];
                rectangle(image, Point(j * CellSize, i * CellSize), Point(j * CellSize + CellSize, i * CellSize + CellSize), (c == 1 ? 0 : 200), cv::FILLED);
            }
        }

        imshow(windowName, image);
    }

    const std::string& Barcode::GetStringRepresentation() const
    {
        return barcode;
    }


    /// Integrates environmental input into the barcode, additively.
    void Barcode::Input(cv::Mat& environment)
    {
        for (int j = 0; j < environment.rows; j++) 
        {
            for (int i = 0; i < environment.cols; i++)
            {
                if (environment.at<uchar>(j, i) > 0) 
                {
                    barcode[j * width + i] = 1;
                }
            }
        }
    }


    /// Intersects with the rhs barcode, keeping only
    /// those that were set in both this and rhs.
    void Barcode::Intersect(const Barcode& rhs)
    {
        auto rhsBarcode = rhs.barcode;
        for (auto i = 0; i < rhsBarcode.size(); ++i)
        {
            if (rhsBarcode[i] == 0)
            {
                barcode[i] = 0;
            }
        }
    }


    /// Randomly adds food tiles to the environment.
    bool Barcode::DropTiles(cv::Mat& environment, int x, int y, int numToDrop, bool useActiveCells) const
    {
        std::vector<int> relativeIndices;
        for (auto i = 0; i < GlobalSettings::BarcodeSize * GlobalSettings::BarcodeSize; ++i)
        {
            auto& tile = environment.at<uchar>(y + i / GlobalSettings::BarcodeSize, x + i % GlobalSettings::BarcodeSize);
            auto& cell = barcode[i];
            if (tile == 0 && (!useActiveCells || cell == 1)) relativeIndices.push_back(i);
        }

        // Shuffle the indices.
        std::shuffle(relativeIndices.begin(), relativeIndices.end(), GlobalSettings::RNG);

        // Pick spots to deposit.
        for (auto i : relativeIndices)
        {
            auto tileX = x + i % GlobalSettings::BarcodeSize;
            auto tileY = y + i / GlobalSettings::BarcodeSize;

            // Activate tile.
            environment.at<uchar>(tileY, tileX) = 255;
            --numToDrop;
            if (numToDrop <= 0) break;
        }

        return numToDrop == 0;
    }


    bool Barcode::ExtractTiles(cv::Mat& environment, int x, int y, int numToTake) const
    {
        std::vector<int> relativeIndices;
        for (auto i = 0; i < GlobalSettings::BarcodeSize * GlobalSettings::BarcodeSize; ++i)
        {
            auto& tile = environment.at<uchar>(y + i / GlobalSettings::BarcodeSize, x + i % GlobalSettings::BarcodeSize);
            if (tile == 255) relativeIndices.push_back(i);
        }

        // Shuffle the indices.
        std::shuffle(relativeIndices.begin(), relativeIndices.end(), GlobalSettings::RNG);

        // Pick spots to extract from.
        for (auto i : relativeIndices)
        {
            auto tileX = x + i % GlobalSettings::BarcodeSize;
            auto tileY = y + i / GlobalSettings::BarcodeSize;

            // Deactivate tile.
            environment.at<uchar>(tileY, tileX) = 0;
            --numToTake;
            if (numToTake <= 0) break;
        }

        return numToTake == 0;
    }


    /// Replaces the string representation of the current barcode.
    void Barcode::SetStringRepresentation(const std::string& rep)
    {
        barcode = rep;
    }


    /// Subtracts a barcode from this one, not affecting
    /// the rhs.
    void Barcode::Subtract(const Barcode& rhs)
    {
        auto rhsBarcode = rhs.barcode;
        for (auto i = 0; i < rhsBarcode.size(); ++i)
        {
            if (rhsBarcode[i] == 1)
            {
                barcode[i] = 0;
            }
        }
    }

    
    /// Updates the barcode pattern by one step.
    void Barcode::Update(bool usePatternMap, bool useLongPatterns, float& vitality)
    {
        auto oldBarcode = barcode;

        if (!usePatternMap)
        {
            for (auto&[key, val] : behaviourGenes)
            {
                std::string pattern = Helpers::GetParentPattern(key);
                if (pattern.size() <= 3) Update1D(pattern, val, oldBarcode, vitality);
                else if (pattern.size() == 9) Update2D(pattern, val, oldBarcode, vitality);
            }
        }
        else
        {
            // Do the 1D genes first.
            for (auto&[key, val] : behaviourGenes)
            {
                if (key >= 10) break;
                std::string pattern = Helpers::GetParentPattern(key);
                Update1D(pattern, val, oldBarcode, vitality);
            }

            // Do the 3x3 2D genes next.
            Update2DWithPatternMap(oldBarcode, 3, vitality);

            // Do the 5x5 2D genes next...
            if (useLongPatterns) Update2DWithPatternMap(oldBarcode, 5, vitality);
        }

		vitality += GlobalSettings::VitalityChangePerUpdate;
    }


    /// Updates the world map with a small probability.
    /// Only finishes the action if it finds enough tiles to replace the ones added.
    /// Returns whether the update was successful.
    bool Barcode::UpdateWorld(cv::Mat& environment, int x, int y, double probability)
    {
        std::uniform_real_distribution<> dist(0.0, 1.0);
        std::vector<Point> pointsToAdd;
        std::vector<int> removablePoints;

        auto count = 0;
        for (auto i = 0; i < barcode.size(); ++i)
        {
            int tileX = x + i % width;
            int tileY = y + i / width;
            if (barcode[i] == 1 && environment.at<uchar>(tileY, tileX) == 0 && dist(GlobalSettings::RNG) < probability)
            {
                pointsToAdd.push_back(Point(tileX, tileY));
                ++count;
            }

            if (barcode[i] == 0 && environment.at<uchar>(tileY, tileX) == 255) removablePoints.push_back(i);
        }
        
        // We have failed to update if there are fewer tiles to remove.
        if (removablePoints.size() < count) return false;

        // Otherwise just fill the spots immediately.
        for (auto& p : pointsToAdd)
        {
            environment.at<uchar>(p.y, p.x) = 255;
        }

        // Shuffle the removable points.
        std::shuffle(removablePoints.begin(), removablePoints.end(), GlobalSettings::RNG);

        // Pick spots to deposit.
        for (auto i : removablePoints)
        {
            if (count == 0) break;
            auto tileX = x + i % GlobalSettings::BarcodeSize;
            auto tileY = y + i / GlobalSettings::BarcodeSize;

            environment.at<uchar>(tileY, tileX) = 0;
            --count;
        }

        return true;
    }


	inline void Barcode::InterpretGeneValue(uchar& value, int& behaviourEffect, int& vitalityEffect)
	{
		behaviourEffect = value / 2;
		vitalityEffect = 2 * (value % 2) - 1;
	}


	void Barcode::Update1D(std::string& pattern, uchar& replacement, std::string& oldBarcode, float& vitality, std::string* updateInto)
    {
        std::string& update = updateInto == nullptr ? barcode : *updateInto;
		int vitalityEffect = 0;
		int behaviourEffect = 0;
		InterpretGeneValue(replacement, behaviourEffect, vitalityEffect);

		int pVitality = 0;
#pragma omp parallel for reduction(+: pVitality)
        for (int j = 0; j < height; ++j)
        {
            std::string subBarcode = oldBarcode.substr(j * width, width);
            std::vector<size_t> positions;
            positions.reserve(GlobalSettings::BarcodeSize);

            size_t pos = subBarcode.find(pattern, 0);
            while (pos != std::string::npos)
            {
                positions.push_back(pos);
                pos = subBarcode.find(pattern, pos + 1);
            }

            // Replace in the new pattern.
            if (pattern.size() == 1)
            {
                for (auto& pos : positions)
                {
                    update[width * j + pos] = behaviourEffect;
					pVitality += vitalityEffect;
                }
            }
            else if (pattern.size() == 3)
            {
                for (auto& pos : positions)
                {
                    update[width * j + pos + 1] = behaviourEffect;
					pVitality += GlobalSettings::Patterns3x1 * vitalityEffect; // The multiplier normalizes for probability.
                }
            }
        }

		vitality += (float) pVitality / GlobalSettings::Patterns3x3;
    }


    void Barcode::Update2D(std::string& pattern, uchar& replacement, std::string& oldBarcode, float& vitality, std::string* updateInto)
    {
        std::string& update = updateInto == nullptr ? barcode : *updateInto;
		int vitalityEffect = 0;
		int behaviourEffect = 0;
		InterpretGeneValue(replacement, behaviourEffect, vitalityEffect);
		const int patternSize = pattern.size();

		int pVitality = 0;
#pragma omp parallel for reduction(+: pVitality)
        for (int j = 0; j < height - 2; ++j)
        {
            auto subBarcode = oldBarcode.substr(j * width, width);
            std::vector<size_t> firstMatchPositions;
            std::vector<size_t> positions;

            auto firstPatternLine = pattern.substr(0, 3);

            size_t pos = subBarcode.find(firstPatternLine, 0);
            while (pos != std::string::npos)
            {
                firstMatchPositions.push_back(pos);
                pos = subBarcode.find(firstPatternLine, pos + 1);
            }

            // Check that the other lines match too.
            for (auto& pos : firstMatchPositions)
            {
                auto secondLine = oldBarcode.substr((j + 1) * width + pos, 3);
                auto thirdLine = oldBarcode.substr((j + 2) * width + pos, 3);

                if (secondLine == pattern.substr(3, 3) && thirdLine == pattern.substr(6, 3))
                {
                    positions.push_back(pos);
                }
            }

            // Replace in the new pattern.
            for (auto& pos : positions)
            {
                update[width * (j + 1) + pos + 1] = behaviourEffect;
				pVitality += vitalityEffect;
            }
        }

		// This assumes that the patterns are 3x3. Throws an exception otherwise.
        vitality += pVitality ;

		if (patternSize != 9) throw std::exception("Wrong pattern size for the assumptions of this function.");
    }


    void Barcode::Update2DWithPatternMap(std::string& oldBarcode, int patternWidth, float& vitality)
    {
        auto& map = patternWidth <= 3 ? Individual::ShortGenePatternMap : Individual::LongGenePatternMap;

        const int edgeLimit = patternWidth - 1;
        const int replaceOffset = (patternWidth - 1) / 2;

		int pVitality = 0;
#pragma omp parallel for reduction(+: pVitality)
        for (int j = 0; j < height - edgeLimit; ++j)
        {
            for (int i = 0; i < width - edgeLimit; ++i)
            {
                // Get the pattern at this position of the barcode.
                std::string subBarcode;
                for (auto k = j; k < j + patternWidth; ++k)
                {
                    subBarcode += oldBarcode.substr(k * width + i, patternWidth);
                }

                // Find which gene this would require.
                auto& geneIndex = map[subBarcode];

                // Do we have this gene?
                if (behaviourGenes.count(geneIndex) == 0) continue;

				int vitalityEffect = 0;
				int behaviourEffect = 0;
				InterpretGeneValue(behaviourGenes[geneIndex], behaviourEffect, vitalityEffect);

                // Replace at the right position.
                barcode[width * (j + replaceOffset) + i + replaceOffset] = behaviourEffect;
				pVitality += vitalityEffect; // Normalized already.
            }
        }

        vitality += pVitality;

		// This assumes that the patterns are 3x3. Throws an exception otherwise.
		if (patternWidth * patternWidth != 9) throw std::exception("Wrong pattern size for the assumptions of this function.");
    }
}