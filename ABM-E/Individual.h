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

        bool AddDropTile(int numToTake);
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
        int LastCellsActive = 0;
        int Vitality = GlobalSettings::MaxVitality / 2;

        static PatternMap ShortGenePatternMap;
        static PatternMap LongGenePatternMap;

    protected:
        int ProcessWorld();
        int UpdateWorld1D(std::string& pattern, uchar& replacement, std::string& oldWorldString, std::string& update);
        int UpdateWorld2D(std::string& oldBarcode, int patternWidth, std::string& newWorldString);
        static void InterpretInteractionGeneValue(uchar& val, int& vitalityUpdate, uchar& replacement);

        bool Alive = true;
    };
}