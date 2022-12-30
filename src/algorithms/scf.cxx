#include <armadillo>
#include "scf.h"
#include "cluster.h"
#include "molutils.h"
#include "cndo2.h"
#include "algoutils.h"
#include "timedfunctional.h"
#include <chrono>

using namespace arma;

int algo::scf::run_eig_sym(vec &e_a, mat &c_a, mat &f_a, vec &e_b, mat &c_b, mat &f_b)
{
    eig_sym(e_a, c_a, f_a);
    eig_sym(e_b, c_b, f_b);
    return 1;
}

mat algo::scf::x_component(mat p_a, mat p_b, mat bonding_params)
{
    mat p_tot = p_a + p_b;
    p_tot %= bonding_params;
    return p_tot;
}

void algo::scf::calcSCFEnergy(float threshold, Cluster cluster)
{
    bool mats_eq = false;

    vec atoms = cluster.atomMatrix.col(0);

    // store the alpha and beta density matrices
    mat p_a(cluster.numValenceElectrons, cluster.numValenceElectrons, fill::zeros);
    mat p_b(cluster.numValenceElectrons, cluster.numValenceElectrons, fill::zeros);

    // store the old alpha and beta density matrices
    mat p_a_old(cluster.numValenceElectrons, cluster.numValenceElectrons, fill::zeros);
    mat p_b_old(cluster.numValenceElectrons, cluster.numValenceElectrons, fill::zeros);

    // store the output fo the cndo2 fock matrix calculation
    mat f_a(cluster.numValenceElectrons, cluster.numValenceElectrons, fill::zeros);
    mat f_b(cluster.numValenceElectrons, cluster.numValenceElectrons, fill::zeros);

    // store the eigenvectors
    mat c_a(cluster.numValenceElectrons, cluster.numValenceElectrons, fill::zeros);
    mat c_b(cluster.numValenceElectrons, cluster.numValenceElectrons, fill::zeros);

    // store the eigenvalues
    vec e_a(cluster.numValenceElectrons, fill::zeros);
    vec e_b(cluster.numValenceElectrons, fill::zeros);

    // matrix that stores the bonding param sum
    mat bonding_params(cluster.numValenceElectrons, cluster.numValenceElectrons, fill::zeros);

    // create timed call to fock matrix generation
    util::timing::TimedFunctional<std::function<algo::cndo2::CNDO2Result(arma::mat p_a, arma::mat p_b, Cluster cluster)>, mat, mat, Cluster> timedFockMatrix(algo::cndo2::fockMatrix, "fockMatrix");
    util::timing::TimedFunctional<std::function<int(arma::vec &, arma::mat &, arma::mat &, arma::vec &, arma::mat &, arma::mat &)>, vec &, mat &, mat &, vec &, mat &, mat &> timed_run_eig_sym(algo::scf::run_eig_sym, "eig_sym");
    while (!mats_eq)
    {
        f_a = timedFockMatrix(p_a, p_b, cluster).fock;

        f_b = timedFockMatrix(p_b, p_a, cluster).fock;

        p_a_old = p_a;
        p_b_old = p_b;

        timed_run_eig_sym(e_a, c_a, f_a, e_b, c_b, f_b);

        p_a = (c_a.cols(0, cluster.p - 1) * c_a.cols(0, cluster.p - 1).t());
        p_b = (c_b.cols(0, cluster.q - 1) * c_b.cols(0, cluster.q - 1).t());

        mats_eq = approx_equal(p_a_old, p_a, "absdiff", threshold) && approx_equal(p_b_old, p_b, "absdiff", threshold);

        // std::cout << "p_a diff" << std::endl;
        // std::cout << p_a_old - p_a << std::endl;
        // std::cout << approx_equal(p_a_old, p_a, "absdiff", threshold) << std::endl;
        // std::cout << "p_b diff" << std::endl;
        // std::cout << p_b_old - p_b << std::endl;
        // std::cout << approx_equal(p_b_old, p_b, "absdiff", threshold) << std::endl;
    }

    algo::cndo2::CNDO2Result cndo2Result = algo::cndo2::fockMatrix(mat(cluster.numValenceElectrons, cluster.numValenceElectrons, fill::zeros), mat(cluster.numValenceElectrons, cluster.numValenceElectrons, fill::zeros), cluster);

    mat h_core = cndo2Result.fock;
    bonding_params = cndo2Result.bonding_param_tot;

    float electron_energy = (1.0 / 2.0) * sum(sum((p_a % (h_core + f_a) + p_b % (h_core + f_b))));

    vec p_a_tot(cluster.atomMatrix.n_rows);
    vec p_b_tot(cluster.atomMatrix.n_rows);

    int idx = 0;

    vec valence_atomic_num_vec(cluster.numValenceElectrons);

    for (int i = 0; i < cluster.atomMatrix.n_rows; i++)
    {
        p_a_tot(i) = sum(diagvec(p_a.submat(idx, idx, idx + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1, idx + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1)));
        p_b_tot(i) = sum(diagvec(p_b.submat(idx, idx, idx + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1, idx + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1)));
        valence_atomic_num_vec(i) = VALENCE_ATOMIC_NUM[atoms(i) - 1];

        idx += VALENCE_ATOMIC_NUM[atoms(i) - 1];
    }

    vec p_tot = p_a_tot + p_b_tot;

    mat y(cluster.atomMatrix.n_rows, cluster.atomMatrix.n_rows);

    int a_idx = 0;
    int b_idx = 0;

    for (int i = 0; i < cluster.atomMatrix.n_rows; i++)
    {
        b_idx = 0;
        for (int j = 0; j < cluster.atomMatrix.n_rows; j++)
        {
            y(i, j) = p_tot(i) * p_tot(j) - valence_atomic_num_vec(j) * p_tot(i) - valence_atomic_num_vec(i) * p_tot(j) - accu(pow(p_a.submat(a_idx, b_idx, a_idx + valence_atomic_num_vec(i) - 1, b_idx + valence_atomic_num_vec(j) - 1), 2) + pow(p_a.submat(a_idx, b_idx, a_idx + valence_atomic_num_vec(i) - 1, b_idx + valence_atomic_num_vec(j) - 1), 2));
            b_idx += VALENCE_ATOMIC_NUM[atoms(j) - 1];
        }
        a_idx += VALENCE_ATOMIC_NUM[atoms(i) - 1];
    }

    mat gradient_electron(3, cluster.atomMatrix.n_rows, fill::zeros);
    mat gradient_nuc(3, cluster.atomMatrix.n_rows, fill::zeros);

    mat ra_overlap = utils::overlapMatrixRA(cluster);

    mat ra_gamma = utils::gammaMatrixRA(cluster);

    mat xs_mask = mat(cluster.numValenceElectrons, cluster.numValenceElectrons, fill::ones) - broadcastToOrbitals(diagmat(vec(cluster.atomMatrix.n_rows, fill::ones)), cluster.atomMatrix, cluster.numValenceElectrons);

    mat density_term_expanded_x = xs_mask % x_component(p_a, p_b, bonding_params) % reshape(mat(ra_overlap.row(0)), cluster.numValenceElectrons, cluster.numValenceElectrons);
    mat density_term_x(cluster.atomMatrix.n_rows, cluster.atomMatrix.n_rows);
    mat density_term_expanded_y = xs_mask % x_component(p_a, p_b, bonding_params) % reshape(mat(ra_overlap.row(1)), cluster.numValenceElectrons, cluster.numValenceElectrons);
    mat density_term_y(cluster.atomMatrix.n_rows, cluster.atomMatrix.n_rows);
    mat density_term_expanded_z = xs_mask % x_component(p_a, p_b, bonding_params) % reshape(mat(ra_overlap.row(2)), cluster.numValenceElectrons, cluster.numValenceElectrons);
    mat density_term_z(cluster.atomMatrix.n_rows, cluster.atomMatrix.n_rows);

    int idx_i = 0;
    int idx_j = 0;

    for (int i = 0; i < cluster.atomMatrix.n_rows; i++)
    {
        int idx_j = 0;
        for (int j = 0; j < cluster.atomMatrix.n_rows; j++)
        {
            density_term_x(i, j) = accu(density_term_expanded_x.submat(idx_i, idx_j, idx_i + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1, idx_j + VALENCE_ATOMIC_NUM[atoms(j) - 1] - 1));
            density_term_y(i, j) = accu(density_term_expanded_y.submat(idx_i, idx_j, idx_i + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1, idx_j + VALENCE_ATOMIC_NUM[atoms(j) - 1] - 1));
            density_term_z(i, j) = accu(density_term_expanded_z.submat(idx_i, idx_j, idx_i + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1, idx_j + VALENCE_ATOMIC_NUM[atoms(j) - 1] - 1));

            idx_j += VALENCE_ATOMIC_NUM[atoms(j) - 1];
        }
        idx_i += VALENCE_ATOMIC_NUM[atoms(i) - 1];
    }

    gradient_electron.row(0) = sum(density_term_x) + sum((y - diagmat(diagvec(y))) % reshape(mat(ra_gamma.row(0)), cluster.atomMatrix.n_rows, cluster.atomMatrix.n_rows));
    gradient_electron.row(1) = sum(density_term_y) + sum((y - diagmat(diagvec(y))) % reshape(mat(ra_gamma.row(1)), cluster.atomMatrix.n_rows, cluster.atomMatrix.n_rows));
    gradient_electron.row(2) = sum(density_term_z) + sum((y - diagmat(diagvec(y))) % reshape(mat(ra_gamma.row(2)), cluster.atomMatrix.n_rows, cluster.atomMatrix.n_rows));

    for (int i = 0; i < cluster.atomMatrix.n_rows; i++)
    {
        for (int j = 0; j < cluster.atomMatrix.n_rows; j++)
        {
            if (j != i)
            {
                gradient_nuc.row(0).col(i) += (CONVERSION_FACTOR * -1 * (VALENCE_ATOMIC_NUM[cluster.atomMatrix(i, 0) - 1] * VALENCE_ATOMIC_NUM[cluster.atomMatrix(j, 0) - 1]) * (cluster.atomMatrix.row(i).col(1) - cluster.atomMatrix.row(j).col(1)) * pow(calcDistance(cluster.atomMatrix.row(i).cols(1, 3), cluster.atomMatrix.row(j).cols(1, 3)), -3));
                gradient_nuc.row(1).col(i) += (CONVERSION_FACTOR * -1 * (VALENCE_ATOMIC_NUM[cluster.atomMatrix(i, 0) - 1] * VALENCE_ATOMIC_NUM[cluster.atomMatrix(j, 0) - 1]) * (cluster.atomMatrix.row(i).col(2) - cluster.atomMatrix.row(j).col(2)) * pow(calcDistance(cluster.atomMatrix.row(i).cols(1, 3), cluster.atomMatrix.row(j).cols(1, 3)), -3));
                gradient_nuc.row(2).col(i) += (CONVERSION_FACTOR * -1 * (VALENCE_ATOMIC_NUM[cluster.atomMatrix(i, 0) - 1] * VALENCE_ATOMIC_NUM[cluster.atomMatrix(j, 0) - 1]) * (cluster.atomMatrix.row(i).col(3) - cluster.atomMatrix.row(j).col(3)) * pow(calcDistance(cluster.atomMatrix.row(i).cols(1, 3), cluster.atomMatrix.row(j).cols(1, 3)), -3));
            }
        }
    }

    std::cout << "gradient (Nuclear part)" << std::endl
              << gradient_nuc << std::endl;

    std::cout << "gradient (Electron part)" << std::endl
              << gradient_electron << std::endl;

    std::cout << "gradient" << std::endl
              << gradient_electron + gradient_nuc << std::endl;

    // my original goal here was to eliminate the loop and do a sum of the triangular matrix but i couldnt figure out how to do that with the atom distances and my brain was tired
    float nuc_repulsion = 0;
    for (int i = 0; i < cluster.atomMatrix.n_rows; i++)
    {
        for (int j = 0; j < i; j++)
        {
            nuc_repulsion += (CONVERSION_FACTOR * ((VALENCE_ATOMIC_NUM[cluster.atomMatrix(i, 0) - 1] * VALENCE_ATOMIC_NUM[cluster.atomMatrix(j, 0) - 1]) / calcDistance(cluster.atomMatrix.row(i).cols(1, 3), cluster.atomMatrix.row(j).cols(1, 3))));
        }
    }

    std::cout << "Nuclear Repulsion Energy is " << nuc_repulsion << " eV." << std::endl;
    std::cout << "Electronic Energy is " << electron_energy << " eV." << std::endl;
    std::cout << "Total Energy is " << electron_energy + nuc_repulsion << " eV." << std::endl;
}
