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
        Individual(Environment& environment, GeneticCode<ushort> chromosome);
        ~Individual();

        bool BeBorn();
        Individual* Clone(bool ignoreBalance) const;
        void DrawBarcode(std::string& windowName);
        const std::string& GetBarcodeString() const;
        void Kill();
        void Update(cv::Mat& interactableEnvironment, Environment::ColocationMapType& colocations);

        inline bool IsAlive() const
        {
            return Alive;
        }

        static bool DetectCollision(const cv::Rect& thisRect, std::vector<cv::Rect>& regions);

        Environment& ItsEnvironment;
        GeneticCode<ushort> ItsGeneticCode;
        std::unique_ptr<Barcode> CurrentBarcode;

        int Age = 0;
        int X = -1;
        int Y = -1;

        static PatternMap ShortGenePatternMap;
        static PatternMap LongGenePatternMap;

    protected:
		static int nextID;

        bool Alive = true;
		float Vitality = 0.f;
		int ID = nextID++;
    };
}