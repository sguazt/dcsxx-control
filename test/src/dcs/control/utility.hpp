#ifndef DCS_CONTROL_TEST_UTILITY_HPP
#define DCS_CONTROL_TEST_UTILITY_HPP

namespace detail { namespace /*<unnamed>*/ {

template <typename T>
struct complex_cmp
{
    bool operator()(::std::complex<T> const& a, ::std::complex<T> const& b) const
    {
        return a.real() < b.real() || (a.real() == b.real() && a.imag() < b.imag());
    }
};

}} // Namespace detail::<unnamed>*/

#endif // DCS_CONTROL_TEST_UTILITY_HPP
