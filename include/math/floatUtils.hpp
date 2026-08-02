/** \file floatUtils.hpp
 * \author Jared R. Males
 * \brief Floating-point classification utilities that remain reliable under fast-math optimization.
 * \ingroup gen_math_files
 */

#ifndef math_floatUtils_hpp
#define math_floatUtils_hpp

#include <bit>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace mx
{
namespace math
{

namespace floatUtils_detail
{

/// Convert an extended floating-point value to a classifiable double without overflowing finite values.
template <typename realT>
double normalizedDouble( realT value /**< [in] floating-point value to normalize */ )
{
    return static_cast<double>( value / std::numeric_limits<realT>::max() );
}

} // namespace floatUtils_detail

/// Test whether a floating-point value is NaN, including under finite-math-only optimization.
/**
 * \returns true if value is a quiet or signaling NaN, otherwise false.
 *
 * \ingroup gen_math
 */
template <typename realT>
bool isNan( realT value /**< [in] floating-point value to test */ )
{
    static_assert( std::is_floating_point_v<realT>, "isNan requires a floating-point type" );

    if constexpr( std::numeric_limits<realT>::is_iec559 && sizeof( realT ) == sizeof( std::uint32_t ) )
    {
        constexpr std::uint32_t exponentMask = 0x7f800000U;
        constexpr std::uint32_t mantissaMask = 0x007fffffU;
        const std::uint32_t bits = std::bit_cast<std::uint32_t>( value );
        return ( bits & exponentMask ) == exponentMask && ( bits & mantissaMask ) != 0;
    }
    else if constexpr( std::numeric_limits<realT>::is_iec559 && sizeof( realT ) == sizeof( std::uint64_t ) )
    {
        constexpr std::uint64_t exponentMask = 0x7ff0000000000000ULL;
        constexpr std::uint64_t mantissaMask = 0x000fffffffffffffULL;
        const std::uint64_t bits = std::bit_cast<std::uint64_t>( value );
        return ( bits & exponentMask ) == exponentMask && ( bits & mantissaMask ) != 0;
    }
    else
    {
        const double normalized = floatUtils_detail::normalizedDouble( value );
        return isNan( normalized );
    }
}

/// Test whether a floating-point value is finite, including under finite-math-only optimization.
/**
 * \returns true if value is neither infinite nor NaN, otherwise false.
 *
 * \ingroup gen_math
 */
template <typename realT>
bool isFinite( realT value /**< [in] floating-point value to test */ )
{
    static_assert( std::is_floating_point_v<realT>, "isFinite requires a floating-point type" );

    if constexpr( std::numeric_limits<realT>::is_iec559 && sizeof( realT ) == sizeof( std::uint32_t ) )
    {
        constexpr std::uint32_t exponentMask = 0x7f800000U;
        return ( std::bit_cast<std::uint32_t>( value ) & exponentMask ) != exponentMask;
    }
    else if constexpr( std::numeric_limits<realT>::is_iec559 && sizeof( realT ) == sizeof( std::uint64_t ) )
    {
        constexpr std::uint64_t exponentMask = 0x7ff0000000000000ULL;
        return ( std::bit_cast<std::uint64_t>( value ) & exponentMask ) != exponentMask;
    }
    else
    {
        const double normalized = floatUtils_detail::normalizedDouble( value );
        return isFinite( normalized );
    }
}

} // namespace math
} // namespace mx

#endif // math_floatUtils_hpp
