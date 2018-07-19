#include <ctime>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include "Environment.h"
#include "GlobalSettings.h"
#include "Helpers.h"
#include "Individual.h"

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

    // Create an environment and individuals.
    //Environment environment(512, 128);
    //environment.AddRegion(cv::Rect(0, 0, 192, 128), 0.08f);
    //environment.AddRegion(cv::Rect(320, 0, 192, 128), 0.04f);
    //environment.AddRegion(cv::Rect(192, 56, 128, 16), 0.01f);
    
    //Environment environment(128, 128);
    //environment.AddRegion(cv::Rect(0, 0, 128, 128), 0.08f);

    Environment environment(512, 128);
    environment.AddRegion(cv::Rect(0, 0, 200, 128), 0.08f);
    environment.AddRegion(cv::Rect(200, 40, 112, 20), 0.01f);
    environment.AddRegion(cv::Rect(312, 0, 200, 128), 0.08f);

    environment.Initialise({ {4, 1000} }, false, true);

    // Get first individual.
    auto& individual = environment[0];

    std::cout << "Starting [" << numThreads << " threads]\n";
    
    std::string envWindowName = "ABME - Environment";
    namedWindow(envWindowName, WINDOW_AUTOSIZE);
    
    std::string indWindowName = "ABME - Barcode";

    clock_t begin = clock();

    bool running = true;
    bool drawEnvironment = true;
    int crisisTiles = 0;

    while (running)
    {
        environment.Update();
        if (drawEnvironment) environment.Draw(envWindowName);
        auto key = waitKey(1);
        switch (key)
        {
        case 'i':
            drawEnvironment = !drawEnvironment;
            break;
        case 'f':
            std::cout << "Num. active tiles: " << environment.CountActiveTiles(true) << std::endl;
            break;
        case 'q': 
            running = false;
            break;
        case 'c':
            environment.CapturePopulation();
            break;
        case 'r':
            environment.ReleasePopulation();
            break;
        case '[':
            crisisTiles -= 500;
            std::cout << "Set crisis tiles to " << crisisTiles << std::endl;
            break;
        case ']':
            crisisTiles += 500;
            std::cout << "Set crisis tiles to " << crisisTiles << std::endl;
            break;
        case 'x':
            int numTiles = environment.CauseTileCrisis(crisisTiles);
            std::cout << "Caused a crisis by adding " << numTiles << " tiles to the map.\n";
            break;
        }
    }

    clock_t end = clock();
    double elapsed = double(end - begin) / CLOCKS_PER_SEC;
    
    std::cout.precision(5);
    std::cout << "Finished in " << elapsed << " s.\n";

    return 0;
}