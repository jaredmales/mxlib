

#ifdef MXLIB_CUDA

#include <iostream>
#include <cuda_runtime.h>

#include <type_traits>
#include <limits>
#include <cstring>
#include <vector>
#include <fstream>

int main()
{
    // Some checks to make sure this makes sense

    if( sizeof( cudaError_t ) != sizeof( unsigned int ) )
    {
        std::cerr << "cudaError_t does not appear to be the same size as unsigned int\n";
        return -1;
    }

    if( !std::is_same<unsigned int, std::underlying_type<cudaError_t>::type>() )
    {
        std::cerr << "cudaError_t does not appear to be an unsigned int\n";
        return -1;
    }

    const char *p0 = cudaGetErrorName( static_cast<cudaError_t>( -1 ) );

    if( strcmp( p0, "unrecognized error code" ) != 0 )
    {
        std::cerr << "cudaError_t max value is not unrecognized\n";
        return -1;
    }

    const char *p1 = cudaGetErrorName( static_cast<cudaError_t>( static_cast<cudaError_t>( -1 ) - 1 ) );

    if( strcmp( p1, "unrecognized error code" ) != 0 )
    {
        std::cerr << "cudaError_t max value - 1 is not unrecognized\n";
        return -1;
    }

    if( p0 != p1 )
    {
        std::cerr << "cudaGetErrorName returned different pointers\n";
        return -1;
    }

    int ncodes = 0;
    int lastcode = -1;

    std::vector<cudaError_t> errcs;

    for( unsigned int n = 0; n < cudaErrorApiFailureBase * 2; ++n )
    {
        const char *p = cudaGetErrorName( static_cast<cudaError_t>( n ) );

        if( p != p0 )
        {
            errcs.push_back( static_cast<cudaError_t>( n ) );
        }
    }

    std::ofstream fout;
    fout.open( "cudaError_info.hpp" );

    fout << "//This is a generated file\n";
    fout << "//Only used during code generation\n";
    fout << "//It should not be installed\n";
    fout << '\n';
    fout << "#include <vector>\n";
    fout << "#include <string>\n";
    fout << "#include <cuda_runtime.h>\n\n";
    fout << "void cuda_info( std::vector<std::string> & errVals, std::vector<std::string> & errCodes, "
            "std::vector<std::string> & errMsgs)\n";
    fout << "{\n";
    for( size_t n = 0; n < errcs.size(); ++n )
    {
        std::string EN = cudaGetErrorName( errcs[n] );
        fout << "    errVals.push_back(\"" << EN << "\");\n";

        std::string en = "success";
        if( errcs[n] != cudaSuccess )
        {
            en = EN.substr( 9 );
            for( size_t m = 1; m < en.size(); ++m )
            {
                if(en[m] != tolower( en[m] ))
                {
                    en.insert(m, "_");
                    ++m;
                }
            }

            for( size_t m = 0; m < en.size(); ++m )
            {
                en[m] = tolower( en[m] );
            }
        }

        fout << "    errCodes.push_back(\"cuda_" << en << "\");\n";

        fout << "    errMsgs.push_back(\"CUDA: " << cudaGetErrorString( errcs[n] ) << "\");\n";
    }
    fout << "}\n";

    fout.close();

    return 0;
}
#else //MXLIB_CUDA

int main()
{
    return 0;
}

#endif //MXLIB_CUDA
