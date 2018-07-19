#pragma once

#include <fstream>

namespace ABME
{
    /// Structure that logs to both console and file.
    class Logger
    {
    public:
        Logger();
        ~Logger();

        void operator<<(const std::string& log);

        static Logger& Instance();

    private:
        static std::unique_ptr<Logger> _instance;
        std::ofstream LogFile;
    };
}