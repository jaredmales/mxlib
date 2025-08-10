

#include <cerrno>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <format>
#include <algorithm>

#include "errno_info.hpp"

int main()
{
    std::vector<std::string> error_ts, error_t_msgs;

    std::ifstream fin;
    fin.open( "error_t_list.txt" );

    while( fin )
    {
        std::string line;
        std::getline( fin, line );

        size_t start = line.find_first_not_of( " \t" );

        if( start == std::string::npos )
        {
            continue;
        }

        size_t space = line.find( ' ', start );

        if( space == std::string::npos )
        {
            continue;
        }

        size_t quotest = line.find( '"', space );

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

        error_ts.push_back( line.substr( start, space - start ) );
        error_t_msgs.push_back( line.substr( quotest, quoteed - quotest ) );
    }
    fin.close();

    std::vector<std::string> ERRNOs, errnos, errno_msgs;
    std::vector<int> EVALs;
    errno_info( ERRNOs, EVALs, errnos, errno_msgs );

    // Get max length of a name for formatting
    int maxlen = 0;
    for( auto &en : error_ts )
    {
        if( en.length() > maxlen )
        {
            maxlen = en.length();
        }
    }

    for( auto &en : errnos )
    {
        if( en.length() > maxlen )
        {
            maxlen = en.length();
        }
    }

    std::vector<std::string> uniqueERRNOs, uniqueerrnos;
    for(size_t n = 1; n < ERRNOs.size(); ++n)
    {
        bool unique = true;
        for(size_t m = 0; m < n; ++m)
        {
            if(EVALs[n] == EVALs[m])
            {
                unique = false;
                break;
            }
        }

        if(unique)
        {
            uniqueERRNOs.push_back(ERRNOs[n]);
            uniqueerrnos.push_back(errnos[n]);
        }
    }

    std::cerr << "Unique ECODES: " << uniqueERRNOs.size() << '/' << ERRNOs.size() << '\n';

    std::ofstream fout;
    fout.open( "../error_t.hpp" );

    fout << "/** \\file error_t.hpp" << '\n';
    fout << " *  \\brief The mxlib error_t type and utilities" << '\n';
    fout << " *  \\ingroup error_handling_files" << '\n';
    fout << " */" << '\n';
    fout << '\n';
    fout << "/* **********  THIS FILE IS GENERATED  ********** */" << '\n';
    fout << "/* ********** DO NOT MODIFY OR COMMIT  ********** */" << '\n';
    fout << '\n';
    fout << "#ifndef mx_errno_t_hpp" << '\n';
    fout << "#define mx_errno_t_hpp" << '\n';
    fout << '\n';
    fout << "namespace mx" << '\n';
    fout << "{\n" << '\n';

    fout << "enum class error_t" << '\n';
    fout << "{" << '\n';
    for( size_t n = 0; n < error_ts.size(); ++n )
    {
        std::string symbol = error_ts[n] + ',';
        fout << std::format( "    {:{}} ///< {}\n", symbol, maxlen + 1, error_t_msgs[n] );
    }

    for( size_t n = 0; n < errnos.size() - 1; ++n )
    {
        error_ts.push_back( errnos[n] );
        error_t_msgs.push_back( errno_msgs[n] );

        std::string symbol = errnos[n] + ',';
        fout << std::format( "    {:{}} ///< {} ({})\n", symbol, maxlen+1, errno_msgs[n], ERRNOs[n] );
    }
    error_ts.push_back( errnos.back() );
    error_t_msgs.push_back( errno_msgs.back() );
    fout << std::format( "    {:{}} ///< {} ({})\n", errnos.back(), maxlen+1, errno_msgs.back(), ERRNOs.back() );

    fout << "};" << "\n" << '\n';

    fout << "static constexpr const char * errorName( const error_t & errc)" << '\n';
    fout << "{" << '\n';
    fout << "    switch(errc)" << '\n';
    fout << "    {" << '\n';

    for( size_t n = 0; n < error_ts.size(); ++n )
    {
        fout << "        case error_t::" << error_ts[n] << ":" << '\n';
        fout << "            return \"" << error_ts[n] << "\";" << '\n';
    }
    fout << "        default:" << '\n';
    fout << "            return \"unknown error_t (bug)\";" << '\n';
    fout << "    }" << '\n';
    fout << "}" << '\n';
    fout << '\n';

    fout << "static constexpr const char * errorMessage( const error_t & errc)" << '\n';
    fout << "{" << '\n';
    fout << "    switch(errc)" << '\n';
    fout << "    {" << '\n';

    for( size_t n = 0; n < error_ts.size(); ++n )
    {
        fout << "        case error_t::" << error_ts[n] << ":" << '\n';
        fout << "            return \"" << error_t_msgs[n] << "\";" << '\n';
    }
    fout << "        default:" << '\n';
    fout << "            return \"unknown error_t (bug)\";" << '\n';
    fout << "    }" << '\n';
    fout << "}" << '\n';
    fout << '\n';

    fout << "static constexpr error_t errno2error_t( const int & err)" << '\n';
    fout << "{" << '\n';
    fout << "    switch(err)" << '\n';
    fout << "    {" << '\n';

    for( size_t n = 0; n < uniqueERRNOs.size(); ++n )
    {
        fout << "        case " << uniqueERRNOs[n] << ":" << '\n';
        fout << "            return error_t::" << uniqueerrnos[n] << ";" << '\n';
    }
    fout << "        default:" << '\n';
    fout << "            return error_t::error;" << '\n';
    fout << "    }" << '\n';
    fout << "}" << '\n';

    fout << "} //namespace mx" << '\n';
    fout << "#endif //mx_error_t_hpp" << '\n';

    fout.close();

    /*

        mx::readColumns( "c++11_errnos.txt", errnos );

        for( int i = 0; i < errnos.size(); ++i )
        {
            if( errnos[i] == "ENOTSUP" )
            {
                fout << "      #ifdef " << errnos[i] << "" << '\n';
                fout << "      #if (ENOTSUP != EOPNOTSUPP)" << '\n';
                fout << "      case " << errnos[i] << ":" << '\n';
                fout << "          return \"" << errnos[i] << "\";" << '\n';
                fout << "      #endif" << '\n';
                fout << "      #endif " << '\n';

                continue;
            }

            if( errnos[i] == "EOPNOTSUP" )
            {
                fout << "      #ifdef " << errnos[i] << "" << '\n';
                fout << "      case " << errnos[i] << ":" << '\n';
                fout << "          #if (ENOTSUP == EOPNOTSUPP)" << '\n';
                fout << "             return \"ENOTSUP / EOPNOTSUPP\";" << '\n';
                fout << "          #else" << '\n';
                fout << "             return \"" << errnos[i] << "\";" << '\n';
                fout << "          #endif" << '\n';
                fout << "      #endif " << '\n';

                continue;
            }

            if( errnos[i] == "EWOULDBLOCK" )
            {
                fout << "      #ifdef " << errnos[i] << "" << '\n';
                fout << "      #if (EWOULDBLOCK != EAGAIN)" << '\n';
                fout << "      case " << errnos[i] << ":" << '\n';
                fout << "          return \"" << errnos[i] << "\";" << '\n';
                fout << "      #endif" << '\n';
                fout << "      #endif " << '\n';

                continue;
            }

            if( errnos[i] == "EAGAIN" )
            {
                fout << "      #ifdef " << errnos[i] << "" << '\n';
                fout << "      case " << errnos[i] << ":" << '\n';
                fout << "          #if (EWOULDBLOCK == EAGAIN)" << '\n';
                fout << "             return \"EAGIAN / EWOULDBLOCK\";" << '\n';
                fout << "          #else" << '\n';
                fout << "             return \"" << errnos[i] << "\";" << '\n';
                fout << "          #endif" << '\n';
                fout << "      #endif " << '\n';

                continue;
            }

            fout << "      #ifdef " << errnos[i] << "" << '\n';
            fout << "      case " << errnos[i] << ":" << '\n';
            fout << "          return \"" << errnos[i] << "\";" << '\n';
            fout << "      #endif " << '\n';
        }
            */

    return 0;
}
