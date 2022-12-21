#include <armadillo>
#include "cluster.h"
#include "cndo2.h"
#include "molutils.h"
#include "algoutils.h"

using namespace arma;

algo::cndo2::CNDO2Result algo::cndo2::fockMatrix(mat p_a, mat p_b, Cluster cluster)
{
    mat gamma = algo::utils::gammaMatrix(cluster);

    mat p_tot = p_a + p_b;

    mat gamma_expanded = broadcastToOrbitals(gamma, cluster.atomMatrix, cluster.numValenceElectrons);

    vec gamma_z_tot = (gamma - diagmat(diagvec(gamma))) * cluster.valenceElectronCountsVec;

    vec diag_gamma_z = (cluster.valenceElectronCountsVec.col(0) - (vec(cluster.valenceElectronCountsVec.n_rows, fill::ones) / 2.0)) % diagvec(gamma);

    vec expanded_diag_gamma_z(cluster.numValenceElectrons);
    vec expanded_gamma_z_tot(cluster.numValenceElectrons);
    vec ion_energy_electron_affinity(cluster.numValenceElectrons);
    vec atoms = cluster.atomMatrix.col(0);
    mat h_off_diag(cluster.numValenceElectrons, cluster.numValenceElectrons, fill::ones);
    h_off_diag -= diagmat(vec(cluster.numValenceElectrons, fill::ones));
    vec density_tot(cluster.numValenceElectrons);
    mat g_off_diag(cluster.numValenceElectrons, cluster.numValenceElectrons, fill::ones);
    g_off_diag -= diagmat(vec(cluster.numValenceElectrons, fill::ones));
    vec density_gamma_off_diag(cluster.numValenceElectrons);
    vec density_tot_by_atom(atoms.n_elem);

    mat bonding_params = mat(cluster.numValenceElectrons, cluster.numValenceElectrons, fill::zeros);

    int k = 0;
    for (int i = 0; i < cluster.atomMatrix.n_rows; i++)
    {
        ion_energy_electron_affinity.subvec(k, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1) = vec(ATOM_TO_IONIZATION_ENERGY_ELECTRON_AFFINITY_PARAMS_MAPPING[atoms(i) - 1]);
        expanded_diag_gamma_z.subvec(k, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1).fill(diag_gamma_z(i));
        expanded_gamma_z_tot.subvec(k, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1).fill(gamma_z_tot(i));
        bonding_params.rows(k, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1) += mat(VALENCE_ATOMIC_NUM[atoms(i) - 1], cluster.numValenceElectrons).fill(ATOMIC_BONDING_PARAMETERS[atoms(i) - 1]);
        bonding_params.cols(k, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1) += mat(cluster.numValenceElectrons, VALENCE_ATOMIC_NUM[atoms(i) - 1]).fill(ATOMIC_BONDING_PARAMETERS[atoms(i) - 1]);
        float density_a = sum(diagvec(p_tot.submat(k, k, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1)));
        density_tot.subvec(k, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1).fill(density_a);
        density_tot_by_atom(i) = density_a;

        k += VALENCE_ATOMIC_NUM[atoms(i) - 1];
    }

    vec gamma_density_sum = (gamma - diagmat(diagvec(gamma))) * density_tot_by_atom;
    vec g_diag = (density_tot - diagvec(p_a)) % diagvec(gamma_expanded) + diagvec(broadcastToOrbitals(diagmat(gamma_density_sum), cluster.atomMatrix, cluster.numValenceElectrons));
    g_off_diag = g_off_diag % (-1.0 * p_a % gamma_expanded);

    h_off_diag %= (1.0 / 2.0 * bonding_params % algo::utils::overlapMatrix(cluster));

    vec h_diag = -1 * ion_energy_electron_affinity - expanded_diag_gamma_z - expanded_gamma_z_tot;
    struct algo::cndo2::CNDO2Result res;
    res.fock = diagmat(h_diag) + diagmat(g_diag) + h_off_diag + g_off_diag;
    res.bonding_param_tot = bonding_params;

    return res;
}