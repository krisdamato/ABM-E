#pragma once

#include <map>
#include <opencv2/highgui.hpp>
#include "Helpers.h"

namespace ABME
{
    class Individual;

    enum DrawMode
    {
        DrawModeLength,
        DrawModeMutIns,
        DrawModeMutDel,
        DrawModeMutFlip,
        DrawModeMutTrans,
        DrawModeMutMeta,
        DrawModeAge,
        DrawModeBackground,
    };


    class Environment
    {
    public:
        using ColocationMapType = std::multimap<cv::Vec2i, Individual*, Vec2iComparator>;

        Environment(int width, int height);

        void AddPopulation(int numIndividuals, int geneticLength, bool useSameGeneIndices, bool useSimpleGenesFirst);
        void AddRegion(cv::Rect region, float activeProbability);
        void CapturePopulation();
        int CauseTileCrisis(int numTilesToAdd);
        void ClampPositions(int& x, int& y) const;
        int CountActiveTiles() const;
        int CountActiveTiles(int regionIndex) const;
        void Draw(std::string& windowName) const;
        cv::Mat& GetMap();
        std::vector<cv::Rect>& GetRegions();
        void Initialise(std::map<int, int> lengthCounts, bool useSameGeneIndices, bool useSimpleGenesFirst);
        void InitialiseTiles();
        void RegisterActiveTileAddition(int regionIndex, int numTiles);
        void ReleasePopulation();
        void RunMetrics(int& killed, int& born, int& diedNaturally) const;
        void ToggleDrawMode();
        void Update();

        Individual& operator[](int index);

        bool PopulationCaptured = false;

    protected:
        void GenerateRandomTiles(cv::Rect& region, int numTiles);
        void GenerateRandomTiles(int numTiles);
        void BurnBarcode(cv::Mat& map, Individual& individual);
        void ReplenishTiles();

        ColocationMapType Colocations;
        cv::Mat Map;
        cv::Mat Snapshot;
        std::vector<std::unique_ptr<Individual>> Individuals;
        std::vector<std::unique_ptr<Individual>> Captured;
        std::vector<cv::Rect> Regions;
        std::vector<int> InitialRegionActiveTiles;
        std::vector<int> NumActiveTilesToAdd;
        DrawMode drawMode = DrawMode::DrawModeLength;
    };
}
