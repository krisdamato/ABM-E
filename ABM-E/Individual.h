#pragma once

#include "Environment.h"
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
        Individual(Environment& environment, Chromosome chromosome);
        ~Individual();

        void DrawBarcode(std::string& windowName);
        bool AddDropFood(cv::Mat& environment, std::string& representation, int& numTiles);
        const std::string& GetBarcodeString() const;
        int InteractWithEnvironment(std::string& regionRepresenation);
        void Update(cv::Mat& baseEnvironment, cv::Mat& interactableEnvironment, Environment::ColocationMapType& colocations);

        Environment& ItsEnvironment;
        Chromosome ItsChromosome;
        std::unique_ptr<Barcode> CurrentBarcode;

        int Age = 0;
        int X = -1;
        int Y = -1;
        int LastCellsActive = 0;
        bool Alive = true;
    };
}