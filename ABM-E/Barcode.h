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

        void ComputeMetrics(cv::Vec2i& movement) const;
        int CountLiveCells() const;
        void Draw(std::string& windowName) const;
        bool DropTiles(cv::Mat& environment, int x, int y, int numFoodTiles, bool useActiveCells) const;
        bool ExtractTiles(cv::Mat& environment, int x, int y, int numFoodTiles) const;
        const std::string& GetStringRepresentation() const;
        void Input(cv::Mat& environment);
        void Intersect(const Barcode& rhs);
        void SetStringRepresentation(const std::string& rep);
        void Subtract(const Barcode& rhs);
        void Update(bool usePatternMap, bool useLongPatterns, float& vitality);
        bool UpdateWorld(cv::Mat& environment, int x, int y, double probability);

    protected:
		inline void InterpretGeneValue(uchar& value, int& behaviourEffect, int& vitalityEffect);
        inline void Update1D(std::string& pattern, uchar& replacement, std::string& oldBarcode, float& vitality, std::string* updateInto = nullptr);
        inline void Update2D(std::string& pattern, uchar& replacement, std::string& oldBarcode, float& vitality, std::string* updateInto = nullptr);
        inline void Update2DWithPatternMap(std::string& oldBarcode, int patternWidth, float& vitality);

        GeneSet& behaviourGenes;
        std::string barcode;
        int width;
        int height;

        static const int CellSize = 16;
    };
}