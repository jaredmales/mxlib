/** \file generate_errno_info.cpp
 * \brief Generates errno_info.hpp from the configured errno name list.
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

    std::vector<std::string> ERRNOS, errnos;
    fin.open("c++11_errnos.txt");
    while(fin)
    {
        std::string EN;
        fin >> EN;
        if(fin)
        {
            ERRNOS.push_back(EN);
            std::string en = EN;
            std::transform(en.begin(), en.end(), en.begin(),[](unsigned char c){ return std::tolower(c); });
            errnos.push_back(en);
        }
    }

    std::ofstream fout;
    fout.open("errno_info.hpp");

    fout << "//This is a generated file\n";
    fout << "//Only used during code generation\n";
    fout << "//It should not be installed\n";
    fout << '\n';
    fout << "#include <vector>\n";
    fout << "#include <string>\n\n";
    fout << "void errno_info( std::vector<std::string> & ERRNOs, std::vector<int> & EVALs, std::vector<std::string> & errnos, std::vector<std::string> & msgs)\n";
    fout << "{\n";
    for(size_t n = 0; n < errnos.size(); ++n)
    {
        fout << "    #ifdef " << ERRNOS[n] << '\n';
        fout << "    ERRNOs.push_back(\"" + ERRNOS[n] + "\");\n";
        fout << "    EVALs.push_back(" + ERRNOS[n] + ");\n";
        fout << "    errnos.push_back(\"" + errnos[n] + "\");\n";
        fout << "    msgs.push_back(strerror(" + ERRNOS[n] + "));\n";
        fout << "    #endif\n\n";
    }

    fout << "}\n";

    fout.close();

    return 0;
}
