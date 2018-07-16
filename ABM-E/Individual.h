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

        bool AddDropFood(cv::Mat& environment, int numToTake);
        bool BeBorn();
        Individual* Clone(bool ignoreBalance) const;
        void DrawBarcode(std::string& windowName);
        const std::string& GetBarcodeString() const;
        int InteractWithEnvironment(std::string& regionRepresenation);
        void Kill();
        void Update(cv::Mat& baseEnvironment, cv::Mat& interactableEnvironment, Environment::ColocationMapType& colocations);

        inline bool IsAlive() const
        {
            return Alive;
        }

        Environment& ItsEnvironment;
        Chromosome ItsChromosome;
        std::unique_ptr<Barcode> CurrentBarcode;

        int Age = 0;
        int X = -1;
        int Y = -1;
        int LastCellsActive = 0;
        int Balance = 0;

    protected:
        bool Alive = true;
    };
}