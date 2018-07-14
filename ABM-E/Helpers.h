#pragma once

#include <exception>
#include <iostream>
#include <opencv2\highgui.hpp>
#include <map>
#include <numeric>
#include <random>
#include <string>
#include <sstream>
#include <vector>
#include "GlobalSettings.h"

namespace ABME
{
    using Chromosome = std::map<int, uchar>;

    class BadGeneIndexException: std::exception
    {
    public:
        BadGeneIndexException(std::string message = "Bad gene index was passed.") : exception(message.c_str())
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
    }
}