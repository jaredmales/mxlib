/** \file generate_fits_status_info.cpp
 * \brief Generates fits_status_info.hpp from the configured CFITSIO status list.
 * \author Jared Males
 */

#include <cerrno>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <format>
#include <algorithm>

int main()
{

    std::ifstream fin;

    std::vector<std::string> ecodes, evals, estrings;

    fin.open("fits_error_codes.txt");
    while(fin)
    {
        std::string line;
        std::getline( fin, line );

        size_t start = line.find_first_not_of( " \t" );

        if( start == std::string::npos )
        {
            continue;
        }

        size_t space1 = line.find( ' ', start );

        if( space1 == std::string::npos )
        {
            continue;
        }

        size_t space2 = space1;
        while(space2 < line.size() && std::isspace(line[space2]))
        {
            ++space2;
        }

        if( space2 == std::string::npos )
        {
            continue;
        }

        size_t space3 = line.find( ' ', space2 );

        if( space3 == std::string::npos )
        {
            continue;
        }



        size_t quotest = line.find( '"', space3 );

        if( quotest == std::string::npos )
        {
            continue;
        }

        ++quotest;

        size_t quoteed = line.find( '"', quotest );

        if( quoteed == std::string::npos )
        {
            continue;
        }

        ecodes.push_back( line.substr( start, space1 - start ) );
        evals.push_back( line.substr( space2, space3 - space2 ));
        estrings.push_back( line.substr( quotest, quoteed - quotest ) );
    }
    fin.close();

    for(size_t n = 0; n < ecodes.size(); ++n)
    {
        evals[n] = ecodes[n];//we'll use the text codes
        std::transform(ecodes[n].begin(), ecodes[n].end(), ecodes[n].begin(),[](unsigned char c){ return std::tolower(c); });

        if(estrings[n].size() == 0)
        {
            estrings[n] = "FITS: " + ecodes[n];
            std::transform(estrings[n].begin(), estrings[n].end(), estrings[n].begin(),[](char c){ if(c=='_') return ' '; else return c; });
        }
        else
        {
            estrings[n] = "FITS: " + estrings[n];
        }

        ecodes[n] = "fits_" + ecodes[n];

    }

    std::ofstream fout;
    fout.open("fits_status_info.hpp");

    fout << "//This is a generated file\n";
    fout << "//Only used during code generation\n";
    fout << "//It should not be installed\n";
    fout << '\n';
    fout << "#include <vector>\n";
    fout << "#include <string>\n\n";
    fout << "void fits_status_info( std::vector<std::string> & ecodes, std::vector<std::string> & evals, std::vector<std::string> & estrings)\n";
    fout << "{\n";
    for(size_t n = 0; n < ecodes.size(); ++n)
    {
        fout << "    ecodes.push_back(\"" + ecodes[n] + "\");" << '\n';
        fout << "    evals.push_back(\"" + evals[n] + "\");" << '\n';
        fout << "    estrings.push_back(\"" + estrings[n] + "\");" << '\n';
        fout << '\n';
    }

    fout << "}\n";

    fout.close();

    return 0;
}
