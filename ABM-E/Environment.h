#pragma once

#include <map>
#include <opencv2/highgui.hpp>
#include "Helpers.h"

namespace ABME
{
    class Individual;

    class Environment
    {
    public:
        using ColocationMapType = std::multimap<cv::Vec2i, Individual*, Vec2iComparator>;

        Environment(int width, int height, float probFood);

        int CountFoodCells() const;
        void Draw(std::string& windowName) const;
        cv::Mat& GetMap();

        void Initialise(std::map<int, int> lengthCounts, bool useSameGeneIndices, bool useSimpleGenesFirst);
        void RegisterFoodAddition(int numTiles);
        void Update();

        Individual& operator[](int index);

    protected:
        void GenerateRandomFood(float probFood);
        void GenerateRandomFood(int numTiles);
        void BurnBarcode(cv::Mat& map, Individual& individual);

        ColocationMapType Colocations;
        cv::Mat Map;
        cv::Mat Snapshot;
        std::vector<std::unique_ptr<Individual>> Individuals;
        int NumFoodCellsToAdd = 0;
    };
}