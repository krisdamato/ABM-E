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

    // Create an environment and individuals.
    Environment environment(512, 512, 0.04f);
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
            std::cout << "Num. food cells: " << environment.CountFoodCells() << std::endl;
            break;
        case 'q': 
            running = false;
            break;
        }
    }

    clock_t end = clock();
    double elapsed = double(end - begin) / CLOCKS_PER_SEC;
    
    std::cout.precision(5);
    std::cout << "Finished in " << elapsed << " s.\n";

    return 0;
}