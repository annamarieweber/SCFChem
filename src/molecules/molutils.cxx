#include <armadillo>
#include "constants.h"
#include "baseexception.h"
#include "molutils.h"

using namespace constants;
using namespace arma;

int countBasisFunctions(mat atomMatrix)
{
    int sum = 0;
    vec atoms = atomMatrix.col(0);
    for (int i = 0; i < atoms.size(); i++)
    {
        sum += ATOM_BASIS_FN_MAP[atoms[i] - 1];
    }
    return sum;
}

int countElectronPairs(mat atomMatrix, bool expectEvenElectrons = false)
{
    int n = countBasisFunctions(atomMatrix);
    if (n / 2 * 2 != n)
    {
        if (expectEvenElectrons)
        {
            throw BaseException("InvalidConfig: number of electron pairs must be an integer");
        }
        return n / 2;
    }
    else
    {
        return n / 2;
    }
}

mat broadcastToOrbitals(mat x, mat atomMatrix, int numValenceElectrons)
{

    vec atoms = atomMatrix.col(0);

    mat orbitalRepresentation(numValenceElectrons, numValenceElectrons);
    int k = 0;
    int l = 0;
    for (int i = 0; i < x.n_rows; i++)
    {
        l = 0;

        for (int j = 0; j < x.n_cols; j++)
        {
            orbitalRepresentation.submat(k, l, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1, l + VALENCE_ATOMIC_NUM[atoms(j) - 1] - 1) = mat(VALENCE_ATOMIC_NUM[atoms(i) - 1], VALENCE_ATOMIC_NUM[atoms(j) - 1]).fill(x(i, j));
            l += VALENCE_ATOMIC_NUM[atoms(j) - 1];
        }
        k += VALENCE_ATOMIC_NUM[atoms(i) - 1];
    }
    return orbitalRepresentation;
}

mat broadcastToOrbitalsRA(mat x, mat atomMatrix, int numValenceElectrons)
{
    vec atoms = atomMatrix.col(0);

    mat xOrbitalRepresentation(numValenceElectrons, numValenceElectrons);
    mat yOrbitalRepresentation(numValenceElectrons, numValenceElectrons);
    mat zOrbitalRepresentation(numValenceElectrons, numValenceElectrons);

    mat orbitalRepresentation(3, pow(numValenceElectrons, 2));

    int k = 0;
    int l = 0;

    for (int i = 0; i < sqrt(x.n_cols); i++)
    {
        l = 0;

        for (int j = 0; j < sqrt(x.n_cols); j++)
        {
            xOrbitalRepresentation.submat(k, l, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1, l + VALENCE_ATOMIC_NUM[atoms(j) - 1] - 1) = mat(VALENCE_ATOMIC_NUM[atoms(i) - 1], VALENCE_ATOMIC_NUM[atoms(j) - 1]).fill(x(0, i * sqrt(numValenceElectrons) + j));
            yOrbitalRepresentation.submat(k, l, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1, l + VALENCE_ATOMIC_NUM[atoms(j) - 1] - 1) = mat(VALENCE_ATOMIC_NUM[atoms(i) - 1], VALENCE_ATOMIC_NUM[atoms(j) - 1]).fill(x(1, i * sqrt(numValenceElectrons) + j));
            zOrbitalRepresentation.submat(k, l, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1, l + VALENCE_ATOMIC_NUM[atoms(j) - 1] - 1) = mat(VALENCE_ATOMIC_NUM[atoms(i) - 1], VALENCE_ATOMIC_NUM[atoms(j) - 1]).fill(x(2, i * sqrt(numValenceElectrons) + j));

            l += VALENCE_ATOMIC_NUM[atoms(j) - 1];
        }
        k += VALENCE_ATOMIC_NUM[atoms(i) - 1];
    }

    std::cout << xOrbitalRepresentation << std::endl;
    orbitalRepresentation.row(0) = vectorise(xOrbitalRepresentation, 1);
    orbitalRepresentation.row(1) = vectorise(yOrbitalRepresentation, 1);
    orbitalRepresentation.row(2) = vectorise(zOrbitalRepresentation, 1);

    return orbitalRepresentation;
}

double calcDistance(mat a1, mat a2)
{
    double distance = 0;
    rowvec coordDists = square(a2.row(0) - a1.row(0));
    return sqrt(sum(coordDists));
}

vec valenceElectronCounts(mat atomMatrix)
{
    mat valenceElectrons = atomMatrix.col(0);
    valenceElectrons.for_each([](vec::elem_type &val)
                              { val = VALENCE_ATOMIC_NUM[val - 1]; });
    return valenceElectrons;
}
