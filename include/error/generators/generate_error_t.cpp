

#include <cerrno>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <format>
#include <algorithm>

#include "errno_info.hpp"
#include "fits_status_info.hpp"

#ifdef MXLIB_CUDA
#include "cudaError_info.hpp"
#endif

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

        if( line[start] == '#' )
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

    std::vector<std::string> fits_codes, fits_vals, fits_msgs;
    fits_status_info( fits_codes, fits_vals, fits_msgs );

    #ifdef MXLIB_CUDA
    std::vector<std::string> cuda_codes, cuda_vals, cuda_msgs;
    cuda_info(cuda_vals, cuda_codes, cuda_msgs);
    #endif

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

    for( auto &ft : fits_codes )
    {
        if( ft.length() > maxlen )
        {
            maxlen = ft.length();
        }
    }

    for( auto &ct : cuda_codes )
    {
        if( ct.length() > maxlen )
        {
            maxlen = ct.length();
        }
    }

    std::vector<std::string> uniqueERRNOs, uniqueerrnos;
    for( size_t n = 1; n < ERRNOs.size(); ++n )
    {
        bool unique = true;
        for( size_t m = 0; m < n; ++m )
        {
            if( EVALs[n] == EVALs[m] )
            {
                unique = false;
                break;
            }
        }

        if( unique )
        {
            uniqueERRNOs.push_back( ERRNOs[n] );
            uniqueerrnos.push_back( errnos[n] );
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
    fout << "#include <fitsio.h>" << '\n';
    #ifdef MXLIB_CUDA
    fout << "#include <cuda_runtime.h>" << '\n';
    #endif
    fout << '\n';
    fout << "#ifndef mx_error_t_hpp" << '\n';
    fout << "#define mx_error_t_hpp" << '\n';
    fout << '\n';
    fout << "namespace mx" << '\n';
    fout << "{\n" << '\n';
    fout << '\n';
    fout << "/// The mxlib error codes" << '\n';
    fout << "/** \\ingroup error_handling_codes" << '\n';
    fout << " */" << '\n';
    fout << "enum class error_t" << '\n';
    fout << "{" << '\n';
    for( size_t n = 0; n < error_ts.size(); ++n )
    {
        std::string symbol = error_ts[n] + ',';
        fout << std::format( "    {:{}} ///< {}\n", symbol, maxlen + 1, error_t_msgs[n] );
    }

    for( size_t n = 0; n < errnos.size(); ++n )
    {
        error_ts.push_back( errnos[n] );
        error_t_msgs.push_back( errno_msgs[n] );

        std::string symbol = errnos[n] + ',';
        fout << std::format( "    {:{}} ///< {} ({})\n", symbol, maxlen + 1, errno_msgs[n], ERRNOs[n] );
    }

    for( size_t n = 0; n < fits_codes.size() ; ++n )
    {
        error_ts.push_back( fits_codes[n] );
        error_t_msgs.push_back( fits_msgs[n] );
        std::string symbol = fits_codes[n] + ',';
        fout << std::format( "    {:{}} ///< {}\n", symbol, maxlen + 1, fits_msgs[n] );
    }

    #ifdef MXLIB_CUDA
    for( size_t n = 0; n < cuda_codes.size(); ++n )
    {
        error_ts.push_back( cuda_codes[n] );
        error_t_msgs.push_back( cuda_msgs[n] );
        std::string symbol = cuda_codes[n] + ',';
        fout << std::format( "    {:{}} ///< {}\n", symbol, maxlen + 1, cuda_msgs[n] );
    }
    #endif
    fout << std::format( "    {:{}} ///< {}\n", "__sentinel", maxlen + 1, "do not use" );


    fout << "}; //enum class error_t" << "\n" << '\n';

    fout << "/// Convert a \\ref error_t code to its name" << '\n';
    fout << "/**" << '\n';
    fout << " * \\returns the name of the \\ref error_t code" << '\n';
    fout << " *" << '\n';
    fout << " * \\ingroup error_handling_codes" << '\n';
    fout << " */" << '\n';
    fout << "static constexpr const char * errorName( const error_t & errc /**< [in] the error code to convert*/)"
         << '\n';
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
    fout << "} //errorName" << '\n';
    fout << '\n';

    fout << "/// Get the descriptive message for a \\ref error_t code." << '\n';
    fout << "/**" << '\n';
    fout << " * \\returns the descriptive message corresponding to the \\ref error_t code" << '\n';
    fout << " *" << '\n';
    fout << " * \\ingroup error_handling_codes" << '\n';
    fout << " */" << '\n';
    fout << "static constexpr const char * errorMessage( const error_t & errc /**< [in] the error code for which to "
            "get the message*/)"
         << '\n';
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
    fout << "} //errorMessage" << '\n';
    fout << '\n';

    fout << "/// Convert an errno code to \\ref error_t" << '\n';
    fout << "/**" << '\n';
    fout << " * \\returns the \\ref error_t code corresponding to the errno code" << '\n';
    fout << " *" << '\n';
    fout << " * \\ingroup error_handling_codes" << '\n';
    fout << " */" << '\n';
    fout << "static constexpr error_t errno2error_t( const int & err/**< [in] the errno code to convert*/)" << '\n';
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
    fout << "} //errno2error_t" << '\n';

    fout << '\n';

    fout << "/// Convert a FITS status code to \\ref error_t" << '\n';
    fout << "/**" << '\n';
    fout << " * \\returns the \\ref error_t code corresponding to the FITS status code" << '\n';
    fout << " *" << '\n';
    fout << " * \\ingroup fits_utils" << '\n';
    fout << " */" << '\n';
    fout << "static constexpr error_t fits_status2error_t( const int & err/**< [in] the fits status code to convert*/)"
         << '\n';
    fout << "{" << '\n';
    fout << "    switch(err)" << '\n';
    fout << "    {" << '\n';
    fout << "        case " << "0" << ":" << '\n'; // this is not specified
    fout << "            return error_t::noerror;" << '\n';

    for( size_t n = 0; n < fits_vals.size(); ++n )
    {
        fout << "        case " << fits_vals[n] << ":" << '\n';
        fout << "            return error_t::" << fits_codes[n] << ";" << '\n';
    }
    fout << "        default:" << '\n';
    fout << "            return error_t::error;" << '\n';
    fout << "    }" << '\n';
    fout << "} //fits_status2error_t" << '\n';
    fout << '\n';

    #ifdef MXLIB_CUDA
    fout << "#ifdef MXLIB_CUDA\n";
    fout << "/// Convert a cudaError code to \\ref error_t" << '\n';
    fout << "/**" << '\n';
    fout << " * \\returns the \\ref error_t code corresponding to the cudaError code" << '\n';
    fout << " *" << '\n';
    fout << " * \\ingroup cuda_utils" << '\n';
    fout << " */" << '\n';
    fout << "static constexpr error_t cudaError2error_t( const cudaError_t & err/**< [in] the cudaError code to convert*/)"
         << '\n';
    fout << "{" << '\n';
    fout << "    switch(err)" << '\n';
    fout << "    {" << '\n';
    for( size_t n = 0; n < cuda_vals.size(); ++n )
    {
        fout << "        case " << cuda_vals[n] << ":" << '\n';
        fout << "            return error_t::" << cuda_codes[n] << ";" << '\n';
    }
    fout << "        default:" << '\n';
    fout << "            return error_t::error;" << '\n';
    fout << "    }" << '\n';
    fout << "} //cudaError2error_t" << '\n';
    fout << "#endif //MXLIB_CUDA\n";
    fout << '\n';
    #endif //MXLIB_CUDA
    fout << "} //namespace mx" << '\n';
    fout << "#endif //mx_error_t_hpp" << '\n';

    fout << "\n\n";
    fout << "#ifdef MXLIBTEST_ERROR_T_TESTS" << '\n';
    fout << "#ifndef MXLIBTEST_ERROR_T_TESTS_INC" << '\n';
    fout << "#define MXLIBTEST_ERROR_T_TESTS_INC" << '\n';
    fout << '\n';
    fout << "namespace mx" << '\n';
    fout << "{" << '\n';
    fout << '\n';
    fout << "void error_t_vector( std::vector<error_t> & errcs )" << '\n';
    fout << '{' << '\n';
    fout << "    errcs = { error_t::" << error_ts[0] << ',' << '\n';
    for( size_t n = 1; n < error_ts.size() - 1; ++n )
    {
        fout << "              error_t::" << error_ts[n] << ',' << '\n';
    }
    fout << "              error_t::" << error_ts.back() << "};" << '\n';
    fout << '\n';
    fout << "} //error_t_vector" << '\n';
    fout <<'\n';
    fout << "void errno_vector( std::vector<int> & errnos )" << '\n';
    fout << '{' << '\n';
    fout << "    errnos = { " << uniqueERRNOs[0] << ',' << '\n';
    for( size_t n = 1; n < uniqueERRNOs.size() - 1; ++n )
    {
        fout << "               " << uniqueERRNOs[n] << ',' << '\n';
    }
    fout << "               " << uniqueERRNOs.back() << "};" << '\n';
    fout << '\n';
    fout << "} //errno_vector" << '\n';
    fout << "void fitserr_vector( std::vector<int> & fitserrs )" << '\n';
    fout << '{' << '\n';
    fout << "    fitserrs= { " << fits_vals[0] << ',' << '\n';
    for( size_t n = 1; n < fits_vals.size() - 1; ++n )
    {
        fout << "                " << fits_vals[n] << ',' << '\n';
    }
    fout << "                " << fits_vals.back() << "};" << '\n';
    fout << '\n';
    fout << "} //fitserr_vector" << '\n';
    fout <<'\n';
    fout << "} //namespace mx" << '\n';
    fout << '\n';
    fout << "#endif //MXLIBTEST_ERROR_T_TESTS_INC" << '\n';
    fout << "#endif //MXLIBTEST_ERROR_T_TESTS" << '\n';
    fout.close();

    return 0;
}
