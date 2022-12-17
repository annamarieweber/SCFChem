#ifndef EXCEPTION_CONSTANTS_H
#define EXCEPTION_CONSTANTS_H

#include <vector>
#include <string>
#include <map>

using std::map;
using std::string;

namespace exception_constants
{
    enum Error_Type
    {
        INVALID_ATOM_IDENTIFIER,
        ATOM_COUNT_MISMATCH,
        INVALID_ATOM_POSITION,
    };

    const map<Error_Type, string> ERROR_MESSAGE = {
        {INVALID_ATOM_IDENTIFIER, "InvalidAtomIdentifier: One or more atom identifier is not a valid atom id."},
        {ATOM_COUNT_MISMATCH, "AtomCountMismatch: Atom Count provided in input file must match number of atoms listed"},
        {INVALID_ATOM_POSITION, "InvalidAtomPosition: Distance between atoms must be > 0. Check coordinates to ensure coordinates are valid."}};
}
#endif