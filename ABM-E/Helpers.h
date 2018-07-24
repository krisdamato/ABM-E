#pragma once

#include <exception>
#include <iostream>
#include <iomanip>
#include <opencv2/highgui.hpp>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <sstream>
#include <ctime>
#include <vector>
#include "GlobalSettings.h"

namespace ABME
{
    using Gene = std::pair<int, uchar>;
    using Chromosome = std::map<int, uchar>;
    using PatternMap = std::map<std::string, int>;

    class BadGeneIndexException: std::runtime_error
    {
    public:
        BadGeneIndexException(std::string message = "Bad gene index was passed.") : runtime_error(message.c_str())
        {

        }
    };


    class Vec2iComparator
    {
    public:
        bool operator()(const cv::Vec2i& v1, const cv::Vec2i& v2) const
        {
            if (v1[0] < v2[0]) return true;
            else if (v1[0] == v2[0] && v1[1] < v2[1]) return true;

            return false;
        }
    };


    class GeneCountComparator
    {
    public:
        bool operator()(const std::tuple<int, int, float>& first, const std::tuple<int, int, float>& second) const
        {
            return std::get<1>(first) > std::get<1>(second);
        }
    };


    namespace Helpers
    {
        /// Swaps pointers.
        inline void Swap(void*& first, void*& second)
        {
            void* temp = first;
            first = second;
            second = temp;
        }


        /// Returns a bit pattern as a vector of chars, with 0 being off, and 1 on.
        inline std::string GetParentPattern(int geneIndex)
        {
            auto subtractor = 0;
            auto patternLength = 0;
            if (geneIndex < 2)
            {
                subtractor = 0;
                patternLength = 1;
            }
            else if (geneIndex < 10)
            {
                subtractor = 2;
                patternLength = 3;
            }
            else if (geneIndex < 522)
            {
                subtractor = 10;
                patternLength = 9;
            }
            else if (geneIndex < 33554432 + 522)
            {
                subtractor = 522;
                patternLength = 25;
            }
            else
            {
                std::stringstream warning;
                warning << "Was passed gene index " << geneIndex << ". Cannot exceed 521";
                throw BadGeneIndexException(warning.str());
            }

            geneIndex -= subtractor;

            // Each pattern is now simply the bit representation of the resulting gene index.
            std::string pattern;
            for (auto i = patternLength - 1; i >= 0; --i)
            {
                pattern.push_back((geneIndex & 1 << i) == 1 << i ? '1' : '0');
            }

            return pattern;
        }


        inline PatternMap GenerateShortPatternMap()
        {
            std::cout << "Generating short patterns up to 3x3. ";

            PatternMap map;
            for (int i = 0; i < 522; ++i)
            {
                map[GetParentPattern(i)] = i;
            }

            std::cout << "Finished.\n";

            return map;
        }

        
        inline PatternMap GenerateLongPatternMap()
        {
            std::cout << "Generating 5x5 patterns. This may take a while... ";

            PatternMap map;
            for (int i = 522; i < GlobalSettings::NumGenes; ++i)
            {
                map[GetParentPattern(i)] = i;
            }

            std::cout << "Finished.\n";

            return map;
        }


        inline std::string ConvertMatToString(cv::Mat& region)
        {
            std::string representation(region.cols * region.rows, '0');
            for (int j = 0; j < region.rows; j++)
            {
                for (int i = 0; i < region.cols; i++)
                {
                    if (region.at<uchar>(j, i) > 0)
                    {
                        representation[j * region.cols + i] = '1';
                    }
                }
            }

            return representation;
        }


        /// Prints rules, given a chromosome of genes.
        inline void PrintRulesFromChromosome(Chromosome chromosome)
        {
            for (auto const& [key, val] : chromosome)
            {
                auto pattern = GetParentPattern(key);
                for (auto j = 0; j < pattern.size(); ++j)
                {
                    std::cout << " | " << pattern[j];
                }
                std::cout << " | => " << val << std::endl;
            }
        }


        /// Generates a random chromosome of the requested length.
        inline Chromosome GenerateRandomChromosome(int length, bool simpleGenesFirst)
        {
            std::uniform_int_distribution<std::mt19937::result_type> dist1(0, 1);

            Chromosome chromosome;
            /// Either chose gene indices randomly...
            if (!simpleGenesFirst)
            {
                std::uniform_int_distribution<std::mt19937::result_type> distIndex(0, GlobalSettings::NumGenes - 1);

                for (auto i = 0; i < length;)
                {
                    int geneIndex = distIndex(GlobalSettings::RNG);
                    if (chromosome.count(geneIndex) == 0)
                    {
                        chromosome[geneIndex] = dist1(GlobalSettings::RNG) == 1 ? '1' : '0';
                        ++i;
                    }
                }
            }
            /// ...or start in order of complexity.
            else
            {
                // Chromosomes of up to length 2, choose randomly from the first
                // 2 possible genes. Chromosomes up to length 10 pick both of the
                // first two genes, and randomly from the set of next 8 genes. 
                // Finally, chromosomes up to length 522 pick all genes up to length
                // 10 and randomly from the rest.
                auto geneIndices = std::vector<int>(GlobalSettings::NumGenes);
                std::iota(std::begin(geneIndices), std::end(geneIndices), 0);

                std::shuffle(geneIndices.begin(), geneIndices.begin() + 2, GlobalSettings::RNG);
                std::shuffle(geneIndices.begin() + 2, geneIndices.begin() + 10, GlobalSettings::RNG);
                std::shuffle(geneIndices.begin() + 10, geneIndices.begin() + 522, GlobalSettings::RNG);

                for (auto i = 0; i < length; ++i)
                {
                    int geneIndex = geneIndices[i];
                    chromosome[geneIndex] = dist1(GlobalSettings::RNG) == 1 ? '1' : '0';
                }
            }

            return chromosome;
        }


        /// Copies the chromosome gene indices but changes the values.
        inline Chromosome GenerateRandomChromosome(Chromosome& prototype)
        {
            std::uniform_int_distribution<std::mt19937::result_type> dist1(0, 1);

            Chromosome chromosome;
            for (auto&[index, value] : prototype)
            {
                chromosome[index] = dist1(GlobalSettings::RNG) == 1 ? '1' : '0';
            }

            return chromosome;
        }


        /// Converts a chromosome into two vectors.
        inline void ConvertChromosomeToVectors(Chromosome& chromosome, std::vector<int>& geneIndices, std::vector<uchar>& geneValues)
        {
            for (auto&[index, value] : chromosome)
            {
                geneIndices.push_back(index);
                geneValues.push_back(value);
            }
        }


        /// Counts chromosomes.
        inline std::map<Chromosome, int> ChromosomeCounts(std::vector<Chromosome>& chromosomes)
        {
            std::map<Chromosome, int> counts;
            for (auto& chr : chromosomes)
            {
                ++counts[chr];
            }

            return counts;
        }


        /// Counts chromosomes.
        inline std::map<Chromosome, int> ChromosomeTypeCounts(std::vector<Chromosome>& chromosomes)
        {
            std::map<Chromosome, int> counts;
            for (auto& chr : chromosomes)
            {
                Chromosome emptyChromosome;
                for (auto&[index, count] : chr) emptyChromosome[index] = 0;
                ++counts[emptyChromosome];
            }

            return counts;
        }


        /// Returns the most popular chromosome and its count.
        inline std::pair<Chromosome, int> MostPopularChromosome(std::vector<Chromosome>& chromosomes, bool typeOnly = false)
        {
            // Get counts.
            auto countMap = typeOnly ? ChromosomeTypeCounts(chromosomes) : ChromosomeCounts(chromosomes);

            Chromosome mostPopular;
            int maximum = 0;
            for (auto&[chr, count] : countMap)
            {
                if (count > maximum)
                {
                    maximum = count;
                    mostPopular = chr;
                }
            }

            return std::pair(mostPopular, maximum);
        }


        /// Returns a set (ordered) of tuples with (gene index, gene counts, percentage that are on, i.e. '1')
        inline std::set<std::tuple<int, int, float>, GeneCountComparator> GeneStatistics(std::vector<Chromosome>& chromosomes)
        {
            std::set<std::tuple<int, int, float>, GeneCountComparator> genePool;
            std::map<int, int> onGenes;
            std::map<int, int> offGenes;
            std::set<int> allGenes;

            for (auto& chr : chromosomes)
            {
                for (auto&[index, value] : chr)
                {
                    if (value == '1') ++onGenes[index];
                    else ++offGenes[index];

                    allGenes.insert(index);
                }
            }

            for (auto& geneIndex : allGenes)
            {
                int count = onGenes[geneIndex] + offGenes[geneIndex];
                float onPercentage = 100 * float(onGenes[geneIndex]) / count;
                
                genePool.insert(std::tuple(geneIndex, count, onPercentage));
            }

            return genePool;
        }


        /// Converts a chromosome to a string representation.
        inline std::string ConvertChromosomeToString(Chromosome& chromosome, bool indicesOnly)
        {
            std::stringstream chrString;
            chrString << "|";
            for (auto&[index, value] : chromosome)
            {
                chrString << index << "|";
            }
            chrString << std::endl;
            if (!indicesOnly)
            {
                chrString << "|";
                for (auto&[index, value] : chromosome)
                {
                    chrString << value << "|";
                }
                chrString << std::endl;
            }

            return chrString.str();
        }


        inline bool PointInsideRects(const cv::Point& point, std::vector<cv::Rect>& rects)
        {
            for (auto& rect : rects)
            {
                if (point.inside(rect)) return true;
            }

            return false;
        }


        inline int PointInRegionIndex(const cv::Point& point, std::vector<cv::Rect>& rects)
        {
            for (int i = 0; i < rects.size(); ++i)
            {
                if (point.inside(rects[i])) return i;
            }

            return -1;
        }


        inline std::string CurrentTimeString()
        {
            auto t = std::time(nullptr);
            auto tm = *std::localtime(&t);

            std::ostringstream oss;
            oss << std::put_time(&tm, "%Y-%m-%d_%H-%M-%S");

            return oss.str();
        }


        template <typename T>
        inline T Crossover(T& first, T& second, std::uniform_real_distribution<>& dist)
        {
            T t = 0;
            for (int i = 0; i < 8 * sizeof(T); ++i)
            {
                bool firstBit = (first & 1UL << i) == (1UL << i);
                bool secondBit = (second & 1UL << i) == (1UL << i);

                bool choice = dist(GlobalSettings::RNG) < 0.5 ? firstBit : secondBit;
                t = choice ? t | (1UL << i) : t & ~(1UL << i);
            }

            return t;
        }


        template <typename T>
        inline T BitFlip(T& operand, std::uniform_real_distribution<>& dist, double flipProbability)
        {
            for (int i = 0; i < 8 * sizeof(T); ++i)
            {
                if (dist(GlobalSettings::RNG) < flipProbability)
                {
                    operand ^= (1UL << i);
                }
            }

            return operand;
        }
    }
}
