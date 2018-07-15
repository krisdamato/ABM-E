#pragma once

#include <opencv2/highgui.hpp>
#include <string.h>
#include "Helpers.h"

namespace ABME
{
    class Barcode
    {
    public:
        Barcode(Chromosome& chromosome, int width, int height);
        Barcode(const Barcode& rhs);

        void ComputeMetrics(cv::Vec2i& movement, int& cellsActive) const;
        int CountLiveCells() const;
        void Draw(std::string& windowName) const;
        void DropTiles(cv::Mat& environment, int x, int y, int& numFoodTiles, bool useActiveCells) const;
        void ExtractTiles(cv::Mat& environment, int x, int y, int& numFoodTiles) const;
        const std::string& GetStringRepresentation() const;
        void Input(cv::Mat& environment);
        void Intersect(const Barcode& rhs);
        void SetStringRepresentation(const std::string& rep);
        void Subtract(const Barcode& rhs);
        void Update();
        void UpdateStringRepresentation(std::string& repr, std::string& output);

    protected:
        inline void Update1D(std::string& pattern, uchar& replacement, std::string& oldBarcode, std::string* updateInto = nullptr);
        inline void Update2D(std::string& pattern, uchar& replacement, std::string& oldBarcode, std::string* updateInto = nullptr);

        Chromosome& chromosome;
        std::string barcode;
        int width;
        int height;

        static const int CellSize = 16;
    };
}