#pragma once

#include "field.hpp"

namespace fields
{

// R, approximated by double
using Real = double;

template <> struct FieldTraits<double>
{
    static constexpr uint64_t characteristic = 0;
    static constexpr uint64_t cardinality = 0;
    static constexpr bool algebraically_closed = false;

    static std::string name()
    {
        return "R";
    }
};

} // namespace fields
