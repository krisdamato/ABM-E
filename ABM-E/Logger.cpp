#include "Logger.h"

#include <ctime>
#include <string>
#include "Helpers.h"

namespace ABME
{
    std::unique_ptr<Logger> Logger::_instance = nullptr;

    Logger::Logger()
    {
        // What is today's date?
        std::string today = Helpers::CurrentTimeString();

        // Create a file with that name and open.
        auto filename = "D:/ABM-E/logs/Log_" + today + ".txt";
        LogFile.open(filename);
        if (!LogFile)
        {
            throw std::exception("Log file couldn't be opened.");
        }
    }


    Logger::~Logger()
    {
        LogFile.close();
    }


    void Logger::operator<<(const std::string& log)
    {
        // Output to screen.
        std::cout << log;

        // Log to file.
        LogFile << log;
    }


    Logger& Logger::Instance()
    {
        if (_instance == nullptr)
        {
            _instance = std::unique_ptr<Logger>(new Logger);
        }

        return *_instance;
    }
}
