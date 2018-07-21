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
        Filename = "D:/ABM-E/logs/Log_" + today + ".txt";
        LogFile.open(Filename);
        if (!LogFile)
        {
            throw std::runtime_error("Log file couldn't be opened.");
        }
	LogFile.close();
    }


    Logger::~Logger()
    {
        LogFile.close();
    }


    void Logger::operator<<(const std::string& log)
    {
	LogFile.open(Filename, std::ios_base::app);

        // Output to screen.
        std::cout << log;

        // Log to file.
        LogFile << log;

	LogFile.close();
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
