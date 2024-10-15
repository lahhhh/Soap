#include "HarmonyWorker.h"
#include "Custom.h"

arma::mat eigen2arma(const Eigen::MatrixXd& m) {

    int nrow = m.rows();
    int ncol = m.cols();
    std::size_t size = nrow * ncol;

    arma::mat a(nrow, ncol);

    std::memcpy(a.memptr(), m.data(),
        size * sizeof(double));

    return a;
}

Eigen::MatrixXd arma2eigen(const arma::mat& a) {

    std::size_t size = a.n_rows * a.n_cols;

    Eigen::MatrixXd m(a.n_rows, a.n_cols);

    std::memcpy(m.data(), a.memptr(),
        size * sizeof(double));

    return m;
}

arma::mat kmeans_centers(const arma::mat& mat, int k) {
    
    auto [centers, cluster] = custom::kmeans_hartigan_wong_mt(arma2eigen(mat), k);

    return eigen2arma(centers);
}

arma::mat harmony_pow(arma::mat A, const arma::vec& T) {

    for (unsigned c = 0; c < A.n_cols; c++) {
        A.unsafe_col(c) = pow(A.unsafe_col(c), arma::as_scalar(T.row(c)));
    }
    return(A);
}

int my_ceil(float num) {
    int inum = (int)num;
    if (num == (float)inum) {
        return inum;
    }
    return inum + 1;
}

arma::vec find_lambda_cpp(const float alpha, const arma::vec& cluster_E) {
    arma::vec lambda_dym_vec(cluster_E.n_rows + 1, arma::fill::zeros);
    lambda_dym_vec.subvec(1, lambda_dym_vec.n_rows - 1) = cluster_E * alpha;
    return lambda_dym_vec;
}

arma::mat safe_entropy(const arma::mat& X) {
    arma::mat A = X % log(X);
    A.elem(find_nonfinite(A)).zeros();
    return(A);
}

HarmonyWorker::HarmonyWorker(
    const Eigen::MatrixXd& mat,
	const QList<QStringList>& metadata_list
) :
	mat_(mat),
	metadata_list(metadata_list),
	window_size(3),
	block_size(0.05),
    max_iter_harmony(10),
	max_iter_kmeans(20),
    ran_setup(false),
    ran_init(false),
    lambda_estimation(false),
    verbose(false)
{}

void HarmonyWorker::run() {

    try {

        N = metadata_list[0].size();
        int nclust = std::min((int)round(N / 30), 100);
        arma::vec sigma = arma::vec(nclust).fill(0.1);

        //--------- one hot --------------//

        int metadata_number = metadata_list.size();
        int n_sample = metadata_list[0].size();
        QVector<int> levels(metadata_number);
        QVector<QMap<QString, int>> order;
        int level = 0;
        for (int i = 0; i < metadata_number; ++i) {
            const auto& metadata = metadata_list[i];
            QStringList factors = custom::unique(metadata);
            QMap<QString, int> sub_order;
            for (const auto& factor : factors) {
                sub_order[factor] = level++;
            }
            order << sub_order;
        }
        arma::mat Phi = arma::zeros(level, n_sample);
        for (int i = 0; i < metadata_number; ++i) {
            const auto& metadata = metadata_list[i];
            const auto& sub_order = order[i];
            for (int j = 0; j < n_sample; ++j) {
                Phi(sub_order[metadata[j]], j) = 1;
            }
        }
        //-----------------------//

        arma::vec theta = arma::vec(level).fill(2.0);
        arma::vec lambda = arma::vec(level + 1).fill(1.0);
        lambda(0) = 0.0;

        this->setup(
            eigen2arma(this->mat_.transpose()),
            arma::sp_mat(Phi),
            sigma,
            theta,
            lambda,
            0.2,
            0.001,
            0.01,
            nclust
        );
        this->init_cluster_cpp();
        this->harmonize();
        emit x_harmony_ready(arma2eigen(this->Z_corr.t()));
        G_TASK_END;
    }
    catch (std::exception& e) {
        G_TASK_WARN(QString::fromUtf8(e.what()));
        G_TASK_END;
    }
    catch (...) {
        G_TASK_WARN("Meeting Error in computation.");
        G_TASK_END;
    }
};

bool HarmonyWorker::harmonize() {

    for (int i = 0; i < max_iter_harmony; ++i) {
        this->cluster_cpp();
        this->moe_correct_ridge_cpp();
        if (this->check_convergence(1)) {
            G_TASK_LOG("Harmony converged after " + QString::number(i + 1) + " iterations.");
            return true;
        }
    }

    G_TASK_WARN("Harmony did not converge.");

    return true;
};

void HarmonyWorker::setup(
    const arma::mat& __Z, 
    const arma::sp_mat& __Phi,
    const arma::vec& __sigma, 
    const arma::vec& __theta, 
    const arma::vec& __lambda, 
    const float __alpha,
    const float __epsilon_kmeans,
    const float __epsilon_harmony,
    const int __K) 
{

    // Algorithm constants
    N = __Z.n_cols;
    B = __Phi.n_rows;
    d = __Z.n_rows;

    Z_orig = __Z;
    Z_cos = arma::normalise(__Z, 2, 0);
    Z_corr = zeros(size(Z_orig));


    Phi = __Phi;
    Phi_t = Phi.t();

    // Create index
    std::vector<unsigned>counters;
    arma::vec sizes(sum(Phi, 1));
    // std::cout << sizes << std::endl;
    for (unsigned i = 0; i < sizes.n_elem; i++) {
        arma::uvec a(int(sizes(i)));
        index.push_back(a);
        counters.push_back(0);
    }

    arma::sp_mat::const_iterator it = Phi.begin();
    arma::sp_mat::const_iterator it_end = Phi.end();
    for (; it != it_end; ++it)
    {
        unsigned int row_idx = it.row();
        unsigned int col_idx = it.col();
        index[row_idx](counters[row_idx]++) = col_idx;
    }

    Pr_b = sum(Phi, 1) / N;


    epsilon_kmeans = __epsilon_kmeans;
    epsilon_harmony = __epsilon_harmony;

    // Hyperparameters
    K = __K;
    if (__lambda(0) == -1) {
        lambda_estimation = true;
    }
    else {
        lambda = __lambda;
    }
    sigma = __sigma;
    theta = __theta;

    allocate_buffers();
    ran_setup = true;

    alpha = __alpha;

}


void HarmonyWorker::allocate_buffers() {

    _scale_dist = arma::zeros<arma::mat>(K, N);
    dist_mat = arma::zeros<arma::mat>(K, N);
    O = E = arma::zeros<arma::mat>(K, B);

    // Hack: create matrix of ones by creating zeros and then add one!
    arma::sp_mat intcpt = arma::zeros<arma::sp_mat>(1, N);
    intcpt = intcpt + 1;

    Phi_moe = join_cols(intcpt, Phi);
    Phi_moe_t = Phi_moe.t();


    W = arma::zeros<arma::mat>(B + 1, d);
}


void HarmonyWorker::init_cluster_cpp() {

    Y = kmeans_centers(Z_cos, K);

    // Cosine normalization of data centrods
    Y = arma::normalise(Y, 2, 0);

    // (2) ASSIGN CLUSTER PROBABILITIES
    // using a nice property of cosine distance,
    // compute squared distance directly with cross product
    dist_mat = 2 * (1 - Y.t() * Z_cos);

    R = -dist_mat;
    R.each_col() /= sigma;
    R = exp(R);
    R.each_row() /= sum(R, 0);

    // (3) BATCH DIVERSITY STATISTICS
    E = sum(R, 1) * Pr_b.t();
    O = R * Phi_t;

    compute_objective();
    objective_harmony.push_back(objective_kmeans.back());

    dist_mat = 2 * (1 - Y.t() * Z_cos); // Z_cos was changed

    ran_init = true;
}

void HarmonyWorker::compute_objective() {
    const float norm_const = 2000 / ((float)N);
    float kmeans_error = arma::as_scalar(accu(R % dist_mat));
    float _entropy = arma::as_scalar(arma::accu(safe_entropy(R).each_col() % sigma)); // NEW: vector sigma
    float _cross_entropy = arma::as_scalar(
        accu((R.each_col() % sigma) % ((arma::repmat(theta.t(), K, 1) % log((O + E) / E)) * Phi)));

    // Push back the data
    objective_kmeans.push_back((kmeans_error + _entropy + _cross_entropy) * norm_const);
    objective_kmeans_dist.push_back(kmeans_error * norm_const);
    objective_kmeans_entropy.push_back(_entropy * norm_const);
    objective_kmeans_cross.push_back(_cross_entropy * norm_const);
}


bool HarmonyWorker::check_convergence(int type) {
    float obj_new, obj_old;
    switch (type) {
    case 0:
        // Clustering 
        // compute new window mean
        obj_old = 0;
        obj_new = 0;
        for (unsigned i = 0; i < window_size; i++) {
            obj_old += objective_kmeans[objective_kmeans.size() - 2 - i];
            obj_new += objective_kmeans[objective_kmeans.size() - 1 - i];
        }
        if ((obj_old - obj_new) / abs(obj_old) < epsilon_kmeans) {
            return(true);
        }
        else {
            return(false);
        }
    case 1:
        // Harmony
        obj_old = objective_harmony[objective_harmony.size() - 2];
        obj_new = objective_harmony[objective_harmony.size() - 1];
        if ((obj_old - obj_new) / abs(obj_old) < epsilon_harmony) {
            return(true);
        }
        else {
            return(false);
        }
    }

    // gives warning if we don't give default return value
    return(true);
}


int HarmonyWorker::cluster_cpp() {
    unsigned iter;

    // Z_cos has changed
    // R has assumed to not change
    // so update Y to match new integrated data  
    for (iter = 0; iter < max_iter_kmeans; iter++) {

        // STEP 1: Update Y (cluster centroids)
        Y = arma::normalise(Z_cos * R.t(), 2, 0);

        dist_mat = 2 * (1 - Y.t() * Z_cos); // Y was changed


        // STEP 3: Update R    
        update_R();

        // STEP 4: Check for convergence
        compute_objective();

        if (iter > window_size) {
            bool convergence_status = check_convergence(0);
            if (convergence_status) {
                iter++;
                break;
            }
        }
    }

    kmeans_rounds.push_back(iter);
    objective_harmony.push_back(objective_kmeans.back());

    return 0;
}


int HarmonyWorker::update_R() {

    // Generate the 0,N-1 indices
    arma::uvec indices = arma::linspace<arma::uvec>(0, N - 1, N);
    update_order = shuffle(indices);

    // Inverse index
    arma::uvec reverse_index(N, arma::fill::zeros);
    reverse_index.rows(update_order) = indices;

    _scale_dist = -dist_mat; // K x N
    _scale_dist.each_col() /= sigma; // NEW: vector sigma
    _scale_dist = exp(_scale_dist);
    _scale_dist = arma::normalise(_scale_dist, 1, 0);

    // GENERAL CASE: online updates, in blocks of size (N * block_size)
    unsigned n_blocks = (int)(my_ceil(1.0 / block_size));
    unsigned cells_per_block = unsigned(N * block_size);

    // Allocate new matrices
    arma::mat R_randomized = R.cols(update_order);
    arma::sp_mat Phi_randomized(Phi.cols(update_order));
    arma::sp_mat Phi_t_randomized(Phi_randomized.t());
    arma::mat _scale_dist_randomized = _scale_dist.cols(update_order);

    for (unsigned i = 0; i < n_blocks; i++) {
        unsigned idx_min = i * cells_per_block;
        unsigned idx_max = ((i + 1) * cells_per_block) - 1; // - 1 because of submat
        if (i == n_blocks - 1) {
            // we are in the last block, so include everything. Up to 19
            // extra cells.
            idx_max = N - 1;
        }

        auto Rcells = R_randomized.submat(0, idx_min, R_randomized.n_rows - 1, idx_max);
        auto Phicells = Phi_randomized.submat(0, idx_min, Phi_randomized.n_rows - 1, idx_max);
        auto Phi_tcells = Phi_t_randomized.submat(idx_min, 0, idx_max, Phi_t_randomized.n_cols - 1);
        auto _scale_distcells = _scale_dist_randomized.submat(0, idx_min, _scale_dist_randomized.n_rows - 1, idx_max);

        // Step 1: remove cells
        E -= sum(Rcells, 1) * Pr_b.t();
        O -= Rcells * Phi_tcells;

        // Step 2: recompute R for removed cells
        Rcells = _scale_distcells;
        Rcells = Rcells % (harmony_pow(E / (O + E), theta) * Phicells);
        Rcells = normalise(Rcells, 1, 0); // L1 norm columns

        // Step 3: put cells back 
        E += sum(Rcells, 1) * Pr_b.t();
        O += Rcells * Phi_tcells;
    }
    this->R = R_randomized.cols(reverse_index);
    return 0;
}


void HarmonyWorker::moe_correct_ridge_cpp() {

    arma::sp_mat _Rk(N, N);
    arma::sp_mat lambda_mat(B + 1, B + 1);

    if (!lambda_estimation) {
        // Set lambda if we have to
        lambda_mat.diag() = lambda;
    }
    Z_corr = Z_orig;
    for (unsigned k = 0; k < K; k++) {
        if (lambda_estimation) {
            lambda_mat.diag() = find_lambda_cpp(alpha, E.row(k).t());
        }
        _Rk.diag() = R.row(k);
        arma::sp_mat Phi_Rk = Phi_moe * _Rk;
        
        arma::mat inv_cov(eigen2arma(arma2eigen(arma::mat(Phi_Rk * Phi_moe_t + lambda_mat)).inverse()));
        // Calculate R-scaled PCs once
        arma::mat Z_tmp = Z_orig.each_row() % R.row(k);

        // Generate the betas contribution of the intercept using the data
        // This erases whatever was written before in W
        W = inv_cov.unsafe_col(0) * sum(Z_tmp, 1).t();
        // Calculate betas by calculating each batch contribution
        for (unsigned b = 0; b < B; b++) {
            // inv_conv is B+1xB+1 whereas index is B long
            W += inv_cov.unsafe_col(b + 1) * sum(Z_tmp.cols(index[b]), 1).t();
        }

        W.row(0).zeros(); // do not remove the intercept
        Z_corr -= W.t() * Phi_Rk;
    }
    Z_cos = arma::normalise(Z_corr, 2, 0);
}

arma::cube HarmonyWorker::moe_ridge_get_betas_cpp() {
    arma::cube W_cube(B + 1, d, K); // rows, cols, slices

    arma::sp_mat _Rk(N, N);
    arma::sp_mat lambda_mat(B + 1, B + 1);

    if (!lambda_estimation) {
        // Set lambda if we have to
        lambda_mat.diag() = lambda;
    }

    for (unsigned k = 0; k < K; k++) {
        _Rk.diag() = R.row(k);
        if (lambda_estimation) {
            lambda_mat.diag() = find_lambda_cpp(alpha, E.row(k).t());
        }
        arma::sp_mat Phi_Rk = Phi_moe * _Rk;
        W_cube.slice(k) = eigen2arma(arma2eigen(arma::mat(Phi_Rk * Phi_moe_t + lambda_mat)).inverse()) * Phi_Rk * Z_orig.t();
        //W_cube.slice(k) = arma::inv(arma::mat(Phi_Rk * Phi_moe_t + lambda_mat)) * Phi_Rk * Z_orig.t();
    }

    return W_cube;
}
