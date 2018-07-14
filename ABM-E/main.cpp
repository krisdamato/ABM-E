#include "Environment.h"
#include "GlobalSettings.h"
#include "Helpers.h"
#include "Individual.h"

#include <opencv2/core.hpp>
#include <opencv2\highgui.hpp>
#include <opencv2/imgproc.hpp>

using namespace ABME;
using namespace cv;

int main()
{
    GlobalSettings::Initialise();

    // Create an environment and individuals.
    Environment environment(512, 512, 0.15f);
    environment.Initialise({ {4, 1000} }, false, true);

    // Get first individual.
    auto& individual = environment[0];

    std::cout << "Starting...\n";
    
    std::string envWindowName = "ABME - Environment";
    std::string indWindowName = "ABME - Barcode";
    
    namedWindow(envWindowName, WINDOW_AUTOSIZE);
    //namedWindow(indWindowName, WINDOW_AUTOSIZE);
    
    bool running = true;
    bool drawEnvironment = true;
    while (running)
    {
        environment.Update();
        if (drawEnvironment) environment.Draw(envWindowName);
        //individual.DrawBarcode(indWindowName);
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
    std::cout << "Finished.\n";

    return 0;
}