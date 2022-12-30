#include <armadillo>
#include "cluster.h"
#include "cndo2.h"
#include "molutils.h"
#include "algoutils.h"
#include "timedfunctional.h"

using namespace arma;

algo::cndo2::CNDO2Result algo::cndo2::fockMatrix(mat p_a, mat p_b, Cluster cluster)
{
    // create timed call to gamma matrix generation
    util::timing::TimedFunctional<std::function<mat(Cluster)>, Cluster> timedGammaMatrix(algo::utils::gammaMatrix);
    // create timed call to broadcast to orbitals generation
    util::timing::TimedFunctional<std::function<mat(mat, mat, int)>, mat, mat, int> timedBroadcastToOrbitals(broadcastToOrbitals);

    mat gamma = timedGammaMatrix(cluster);

    vec idxVec = arma::linspace<arma::vec>(0, cluster.atomMatrix.n_rows - 1, cluster.atomMatrix.n_rows);

    vec valenceElectronCounts = cluster.atomMatrix.col(0);
    valenceElectronCounts.for_each([](mat::elem_type &atom)
                                   { atom = VALENCE_ATOMIC_NUM[atom - 1]; });
    vec endIdx = cumsum(valenceElectronCounts);
    vec startIdx = endIdx - valenceElectronCounts;

    mat p_tot = p_a + p_b;

    mat gamma_expanded = timedBroadcastToOrbitals(gamma, cluster.atomMatrix, cluster.numValenceElectrons);

    vec gamma_z_tot = (gamma - diagmat(diagvec(gamma))) * cluster.valenceElectronCountsVec;

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

    auto fill_params = [&ion_energy_electron_affinity, &expanded_diag_gamma_z, &expanded_gamma_z_tot, &bonding_params, &density_tot_by_atom, &density_tot, p_tot, gamma, gamma_z_tot, startIdx, endIdx, cluster](double &idx)
    {
        ion_energy_electron_affinity.subvec(startIdx[idx], endIdx[idx] - 1) = vec(ATOM_TO_IONIZATION_ENERGY_ELECTRON_AFFINITY_PARAMS_MAPPING[cluster.atomMatrix.col(0)(idx) - 1]);
        expanded_diag_gamma_z.subvec(startIdx[idx], endIdx[idx] - 1).fill(.5);
        expanded_gamma_z_tot.subvec(startIdx[idx], endIdx[idx] - 1).fill(gamma_z_tot(idx));
        bonding_params.rows(startIdx[idx], endIdx[idx] - 1) += mat(VALENCE_ATOMIC_NUM[cluster.atomMatrix.col(0)(idx) - 1], cluster.numValenceElectrons).fill(ATOMIC_BONDING_PARAMETERS[cluster.atomMatrix.col(0)(idx) - 1]);
        bonding_params.cols(startIdx[idx], endIdx[idx] - 1) += mat(cluster.numValenceElectrons, VALENCE_ATOMIC_NUM[cluster.atomMatrix.col(0)(idx) - 1]).fill(ATOMIC_BONDING_PARAMETERS[cluster.atomMatrix.col(0)(idx) - 1]);
        float density_a = sum(diagvec(p_tot.submat(startIdx[idx], startIdx[idx], endIdx[idx] - 1, endIdx[idx] - 1)));
        density_tot.subvec(startIdx[idx], endIdx[idx] - 1).fill(density_a);
        density_tot_by_atom(idx) = density_a;
    };

    idxVec.for_each(fill_params);

    vec gamma_density_sum = (gamma - diagmat(diagvec(gamma))) * density_tot_by_atom;
    vec g_diag = (density_tot - diagvec(p_a)) % diagvec(gamma_expanded) + diagvec(timedBroadcastToOrbitals(diagmat(gamma_density_sum), cluster.atomMatrix, cluster.numValenceElectrons));
    g_off_diag = g_off_diag % (-1.0 * p_a % gamma_expanded);

    h_off_diag %= (1.0 / 2.0 * bonding_params % algo::utils::overlapMatrix(cluster));

    vec h_diag = -1 * ion_energy_electron_affinity - expanded_diag_gamma_z - expanded_gamma_z_tot;
    struct algo::cndo2::CNDO2Result res;
    res.fock = diagmat(h_diag) + diagmat(g_diag) + h_off_diag + g_off_diag;
    res.bonding_param_tot = bonding_params;

    return res;
}