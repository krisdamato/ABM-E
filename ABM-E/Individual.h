#pragma once

#include "Environment.h"
#include "GeneticCode.h"
#include "Helpers.h"

namespace cv
{
    class Mat;
}


namespace ABME
{
    class Barcode;

    class Individual
    {
    public:
        Individual(Environment& environment, GeneticCode chromosome);
        ~Individual();

        bool AddDropTile(int numToTake);
        bool BeBorn();
        Individual* Clone(bool ignoreBalance) const;
        void DrawBarcode(std::string& windowName);
        const std::string& GetBarcodeString() const;
        void Kill();
        void Update(cv::Mat& interactableEnvironment, Environment::ColocationMapType& colocations, std::vector<cv::Rect>& offLimitRegions);

        inline bool IsAlive() const
        {
            return Alive;
        }

        static bool DetectCollision(const cv::Rect& thisRect, std::vector<cv::Rect>& regions);

        Environment& ItsEnvironment;
        GeneticCode ItsGeneticCode;
        std::unique_ptr<Barcode> CurrentBarcode;

        int Age = 0;
        int X = -1;
        int Y = -1;
        int LastCellsActive = 0;
        std::vector<int> Balances;

    protected:
        bool Alive = true;
    };
}