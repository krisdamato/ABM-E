#include "Barcode.h"

#include <omp.h>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <vector>

namespace ABME
{
    using namespace cv;

    PatternMap Barcode::ShortGenePatternMap = Helpers::GenerateShortPatternMap();
    PatternMap Barcode::LongGenePatternMap = Helpers::GenerateLongPatternMap();


    Barcode::Barcode(Chromosome& chromosome, int width, int height) : chromosome(chromosome), width(width), height(height)
    {
        // Create a string to represent the barcode with the proper length.
        barcode = std::string(width * height, '0');
    }


    Barcode::Barcode(const Barcode& rhs) : chromosome(rhs.chromosome), barcode(rhs.barcode), width(rhs.width), height(rhs.height)
    {
        
    }


    /// Computes movement from current barcode, assuming that
    /// barcode pixels are alternating in the directions (+i, +j, -i, -j).
    /// Also computes current number of active cells.
    /// Only considers the inner 14 x 14 cells (not the boundary cells). 
    void Barcode::ComputeMetrics(Vec2i& movement, int& cellsActive) const
    {
        float positiveX = 0.f;
        float positiveY = 0.f;

        for (auto i = width + 1; i < barcode.size() - (width + 1); ++i)
        {
            // Ignore barcode elements on the left/right boundaries.
            auto x = i % width;
            if (x == 0 || x == 15) continue;

            // Note: boundary cells are influenced by border effects.
            if (barcode[i] != '0')
            {
                // Increment consumption.
                cellsActive += 1;

                // Normalise movement increments so that updates are fair 
                // with respect to cells available. 
                switch (i % 4)
                {
                case 0:
                    positiveX += 1.0f;
                    break;
                case 1:
                    positiveY += 0.75f;
                    break;
                case 2:
                    positiveX -= 0.75f;
                    break;
                case 3:
                    positiveY -= 1.0f;
                    break;
                }
            }
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
            if (c == '1') ++count;
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
                rectangle(image, Point(j * CellSize, i * CellSize), Point(j * CellSize + CellSize, i * CellSize + CellSize), (c == '1' ? 0 : 200), CV_FILLED);
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
                    barcode[j * width + i] = '1';
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
            if (rhsBarcode[i] == '0')
            {
                barcode[i] = '0';
            }
        }
    }


    /// Randomly adds food tiles to the environment.
    /// Note: numToTake must be negative here.
    void Barcode::DropTiles(cv::Mat& environment, int x, int y, int& numToTake, std::vector<Rect>& regions, std::vector<int>& balances, bool useActiveCells) const
    {
        std::vector<int> relativeIndices;
        for (auto i = 0; i < GlobalSettings::BarcodeSize * GlobalSettings::BarcodeSize; ++i)
        {
            auto& tile = environment.at<uchar>(y + i / GlobalSettings::BarcodeSize, x + i % GlobalSettings::BarcodeSize);
            auto& cell = barcode[i];
            if (tile == 0 && (!useActiveCells || cell == '1')) relativeIndices.push_back(i);
        }

        // Shuffle the indices.
        std::shuffle(relativeIndices.begin(), relativeIndices.end(), GlobalSettings::RNG);

        // Pick spots to deposit.
        for (auto i : relativeIndices)
        {
            auto tileX = x + i % GlobalSettings::BarcodeSize;
            auto tileY = y + i / GlobalSettings::BarcodeSize;

            // Update region balance.
            auto index = Helpers::PointInRegionIndex(Point(tileX, tileY), regions);
            --balances[index];

            // Activate tile.
            environment.at<uchar>(tileY, tileX) = 255;
            ++numToTake;
            if (numToTake >= 0) break;
        }
    }


    void Barcode::ExtractTiles(cv::Mat& environment, int x, int y, int& numToTake, std::vector<Rect>& regions, std::vector<int>& balances) const
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

            // Update region balance.
            auto index = Helpers::PointInRegionIndex(Point(tileX, tileY), regions);
            ++balances[index];

            // Deactivate tile.
            environment.at<uchar>(tileY, tileX) = 0;
            --numToTake;
            if (numToTake <= 0) break;
        }
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
            if (rhsBarcode[i] == '1')
            {
                barcode[i] = '0';
            }
        }
    }

    
    /// Updates the barcode pattern by one step.
    /// Note: only allows up to 3x3 patterns.
    void Barcode::Update(bool usePatternMap, bool useLongPatterns)
    {
        auto oldBarcode = barcode;

        if (!usePatternMap)
        {
            for (auto&[key, val] : chromosome)
            {
                std::string pattern = Helpers::GetParentPattern(key);
                if (pattern.size() <= 3) Update1D(pattern, val, oldBarcode);
                else if (pattern.size() == 9) Update2D(pattern, val, oldBarcode);
            }
        }
        else
        {
            // Do the 1D genes first.
            for (auto&[key, val] : chromosome)
            {
                if (key >= 10) break;
                std::string pattern = Helpers::GetParentPattern(key);
                Update1D(pattern, val, oldBarcode);
            }

            // Do the 3x3 2D genes next.
            Update2DWithPatternMap(oldBarcode, 3);

            // Do the 5x5 2D genes next...
            if (useLongPatterns) Update2DWithPatternMap(oldBarcode, 5);
        }
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
            if (barcode[i] == '1' && environment.at<uchar>(tileY, tileX) == 0 && dist(GlobalSettings::RNG) < probability)
            {
                pointsToAdd.push_back(Point(tileX, tileY));
                ++count;
            }

            if (barcode[i] == '0' && environment.at<uchar>(tileY, tileX) == 255) removablePoints.push_back(i);
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


    void Barcode::Update1D(std::string& pattern, uchar& replacement, std::string& oldBarcode, std::string* updateInto)
    {
        std::string& update = updateInto == nullptr ? barcode : *updateInto;

#pragma omp parallel for
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
                    update[width * j + pos] = replacement;
                }
            }
            else if (pattern.size() == 3)
            {
                for (auto& pos : positions)
                {
                    update[width * j + pos + 1] = replacement;
                }
            }
        }
    }


    void Barcode::Update2D(std::string& pattern, uchar& replacement, std::string& oldBarcode, std::string* updateInto)
    {
        std::string& update = updateInto == nullptr ? barcode : *updateInto;

#pragma omp parallel for
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
                update[width * (j + 1) + pos + 1] = replacement;
            }
        }
    }


    void Barcode::Update2DWithPatternMap(std::string& oldBarcode, int patternWidth)
    {
        auto& map = patternWidth <= 3 ? ShortGenePatternMap : LongGenePatternMap;

        const int edgeLimit = patternWidth - 1;
        const int replaceOffset = (patternWidth - 1) / 2;

#pragma omp parallel for
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
                if (chromosome.count(geneIndex) == 0) continue;

                // Replace at the right position.
                barcode[width * (j + replaceOffset) + i + replaceOffset] = chromosome[geneIndex];
            }
        }
    }
}