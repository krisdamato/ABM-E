#include <ctime>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include "Environment.h"
#include "GlobalSettings.h"
#include "Helpers.h"
#include "Individual.h"
#include "Logger.h"

using namespace ABME;
using namespace cv;

int main(int argc, char** argv)
{
    // Read the number of threads if passed.
    int numThreads = argc > 1 ? std::atoi(argv[1]) : 6;

    // Initialise global parameters.
    GlobalSettings::Initialise(numThreads);
    GlobalSettings::ForceEqualChromosomeReproductions = false;
    GlobalSettings::DistanceStep = 4;
    GlobalSettings::AllowFreeTileMovement = true;
    GlobalSettings::TileDepositsEqualDifference = false;
    GlobalSettings::MutationRatesEvolve = true;
    GlobalSettings::UseSingleStructuralMutationRate = false;

    // Create a logger.
    auto& logger = Logger::Instance();

    // Create an environment and individuals.
    //Environment environment(512, 128);
    //environment.AddRegion(cv::Rect(0, 0, 192, 128), 0.08f);
    //environment.AddRegion(cv::Rect(320, 0, 192, 128), 0.04f);
    //environment.AddRegion(cv::Rect(192, 56, 128, 16), 0.01f);
    
    //Environment environment(256, 128);
    //environment.AddRegion(cv::Rect(0, 0, 120, 128), 0.05f);
    //environment.AddRegion(cv::Rect(136, 0, 120, 128), 0.01f);
    //environment.AddRegion(cv::Rect(120, 56, 16, 16), 0.00f);

    Environment environment(128, 128);
    environment.AddRegion(cv::Rect(0, 0, 56, 128), 0.05f);
    environment.AddRegion(cv::Rect(72, 0, 56, 128), 0.03f);
    environment.AddRegion(cv::Rect(56, 56, 16, 16), 0.00f);

    //Environment environment(128, 128);
    ////environment.AddRegion(cv::Rect(0, 0, 120, 128), 0.08f);
    //environment.AddRegion(cv::Rect(0, 0, 128, 128), 0.01f);

    environment.Initialise({ { 4, 2000 }, { 5, 2000 } }, false, true);

    std::cout << "Starting [" << numThreads << " threads]\n";
    
    clock_t begin = clock();

    bool running = true;
    bool drawEnvironment = true;
    int crisisTiles = 0;
    int intruderGeneticLength = 4;
    
    std::string envWindowName = "ABME - Environment";
    
    if (drawEnvironment) namedWindow(envWindowName, WINDOW_AUTOSIZE);
    
    while (running)
    {
        std::stringstream log;

        environment.Update();
        if (drawEnvironment) environment.Draw(envWindowName);
        auto key = waitKey(1);
        switch (key)
        {
        case 'i':
            drawEnvironment = !drawEnvironment;
            break;
        case 'f':
            std::cout << "Num. active tiles: " << environment.CountActiveTiles() << std::endl;
            break;
        case 'q': 
            running = false;
            break;
        case 'c':
            environment.CapturePopulation();
            log << "Captured the current population.\n";
            Logger::Instance() << log.str();
            break;
        case 'r':
            environment.ReleasePopulation();
            log << "Released the captured population.\n";
            Logger::Instance() << log.str();
            break;
        case '[':
            crisisTiles -= 500;
            std::cout << "Set crisis tiles to " << crisisTiles << std::endl;
            break;
        case ']':
            crisisTiles += 500;
            std::cout << "Set crisis tiles to " << crisisTiles << std::endl;
            break;
        case 't':
            environment.ToggleDrawMode();
            break;
        case '+':
            intruderGeneticLength++;
            std::cout << "Genetic length of intruder population set to " << intruderGeneticLength << std::endl;
            break;
        case '-':
            intruderGeneticLength--;
            std::cout << "Genetic length of intruder population set to " << intruderGeneticLength << std::endl;
            break;
        case '#':
            environment.AddPopulation(GlobalSettings::CrisisPopulationSize, intruderGeneticLength, false, true);
            std::cout << "Added a new population of size " << GlobalSettings::CrisisPopulationSize << " and genetic length " << intruderGeneticLength << std::endl;
            break;
        case 'x':
            int numTiles = environment.CauseTileCrisis(crisisTiles);
            log << "Caused a crisis by adding " << numTiles << " tiles to the map.\n";
            Logger::Instance() << log.str();
            break;
        }
    }

    clock_t end = clock();
    double elapsed = double(end - begin) / CLOCKS_PER_SEC;
    
    std::cout.precision(5);
    std::cout << "Finished in " << elapsed << " s.\n";

    return 0;
}
