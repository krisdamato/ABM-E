#pragma once

#include <opencv2/highgui.hpp>
#include <string.h>
#include "Helpers.h"

namespace ABME
{
    class Barcode
    {
    public:
        Barcode(GeneSet& chromosome, int width, int height);
        Barcode(const Barcode& rhs);

        void ComputeMetrics(cv::Vec2i& movement, int& cellsActive) const;
        int CountLiveCells() const;
        void Draw(std::string& windowName) const;
        void DropTiles(cv::Mat& environment, int x, int y, int& numFoodTiles, std::vector<cv::Rect>& regions, std::vector<int>& balances, bool useActiveCells) const;
        void ExtractTiles(cv::Mat& environment, int x, int y, int& numFoodTiles, std::vector<cv::Rect>& regions, std::vector<int>& balances) const;
        const std::string& GetStringRepresentation() const;
        void Input(cv::Mat& environment);
        void Intersect(const Barcode& rhs);
        void SetStringRepresentation(const std::string& rep);
        void Subtract(const Barcode& rhs);
        void Update(bool usePatternMap, bool useLongPatterns);
        bool UpdateWorld(cv::Mat& environment, int x, int y, double probability);

    protected:
        inline void Update1D(std::string& pattern, uchar& replacement, std::string& oldBarcode, std::string* updateInto = nullptr);
        inline void Update2D(std::string& pattern, uchar& replacement, std::string& oldBarcode, std::string* updateInto = nullptr);
        inline void Update2DWithPatternMap(std::string& oldBarcode, int patternWidth);

        GeneSet& behaviourGenes;
        std::string barcode;
        int width;
        int height;

        static PatternMap ShortGenePatternMap;
        static PatternMap LongGenePatternMap;
        static const int CellSize = 16;
    };
}