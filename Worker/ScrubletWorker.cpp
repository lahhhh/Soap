#include "ScrubletWorker.h"

#include "Custom.h"
#include "TruncatedSvd.h"
#include "annoylib.h"
#include "kissrandom.h"


bool check_done(const Eigen::ArrayXd & fsim, const Eigen::MatrixXd & sim, double xatol, double fatol){
    double maxsim = 0, maxfsim = 0 ;
    for(int i = 1; i < sim.rows(); ++i){
        double max1 = (sim.row(i) - sim.row(0)).array().abs().maxCoeff();
        double max2 = std::abs(fsim[i] - fsim[0]);
        maxsim = max1 > maxsim ? max1 : maxsim;
        maxfsim = max2 > maxfsim ? max2 : maxfsim;
    }
    if(maxsim <= xatol && maxfsim <= fatol)return true;
    return false;
}

template <typename T>
Eigen::ArrayXd NelderMead(T const & func, Eigen::ArrayXd x0, int maxiter = 0, int maxfun = 0, double xatol = 0.0001, double fatol = 0.0001, bool adaptive = false){
    double dim = x0.size();
    double rho = 1, chi, psi, sigma;
    int fcall = 1;
    if(adaptive){
        chi = 1 + 2 / dim;
        psi = 0.75 - 1 / (2 * dim);
        sigma = 1 - 1 / dim;
    }else{
        chi = 2;
        psi = 0.5;
        sigma = 0.5;
    }
    double nonzdelt = 0.05, zdelt = 0.00025;

    Eigen::MatrixXd sim = Eigen::MatrixXd::Zero(dim + 1, dim);
    sim.row(0) = x0;
    for(int i = 0; i < dim; ++i){
        Eigen::ArrayXd tmp = x0;
        if(tmp(i) != 0){
            tmp[i] = (1 + nonzdelt) * tmp(i);
        }else{
            tmp[i] = zdelt;
        }
        sim.row(i + 1) = tmp;
    }
    if(!maxiter && !maxfun){
        maxiter = dim * 200;
        maxfun = dim * 200;
    }else if(!maxiter){
        if(maxfun == INT_MAX){
            maxiter = dim * 200;
        }else{
            maxiter = INT_MAX;
        }
    }else if(!maxfun){
        if(maxiter == INT_MAX){
            maxfun = dim * 200;
        }else{
            maxfun = INT_MAX;
        }
    }
    Eigen::ArrayXd fsim = Eigen::ArrayXd::Zero(dim + 1);
    for(int i = 0; i < dim + 1; ++i){
        fsim[i] = func(sim.row(i));
        fcall ++;
    }
    _Cs sort_by_first(fsim, sim);
    int iterations = 1;
    while(fcall < maxfun && iterations < maxiter){
        if(check_done(fsim, sim, xatol, fatol))break;
        Eigen::ArrayXd xbar = (sim.colwise().sum() - sim.row(dim)).array() / dim;
        Eigen::ArrayXd xr = (1 + rho) * xbar - rho * sim.row(dim).array().transpose();
        double fxr = func(xr);
        fcall ++;
        bool doShrink = false;
        if(fxr < fsim[0]){
            Eigen::ArrayXd xe = (1 + rho * chi) * xbar - rho * chi * sim.row(dim).array().transpose();
            double fxe = func(xe);
            fcall ++;
            if(fxe < fxr){
                sim.row(dim) = xe;
                fsim[dim] = fxe;
            }else{
                sim.row(dim) = xr;
                fsim[dim] = fxr;
            }
        }else{
            if(fxr < fsim[dim - 1]){
                sim.row(dim) = xr;
                fsim[dim] = fxr;
            }else{
                if(fxr < fsim[dim]){
                    Eigen::ArrayXd xc = (1 + psi * rho) * xbar - psi * rho * sim.row(dim).array().transpose();
                    double fxc = func(xc);
                    fcall ++;
                    if(fxc < fxr){
                        sim.row(dim) = xc;
                        fsim[dim] = fxc;
                    }else{
                        doShrink = true;
                    }
                }else{
                    Eigen::ArrayXd xcc = (1 - psi) * xbar + psi * sim.row(dim).array().transpose();
                    double fxcc = func(xcc);
                    fcall ++;
                    if(fxcc < fsim[dim]){
                        sim.row(dim) = xcc;
                        fsim[dim] = fxcc;
                    }else{
                        doShrink = true;
                    }
                }
                if(doShrink){
                    for(int i = 1; i < dim + 1; ++i){
                        sim.row(i) = sim.row(0) + sigma * (sim.row(i) - sim.row(0));
                        fsim[i] = func(sim.row(i));
                        fcall ++;
                    }
                }
            }
        }
        _Cs sort_by_first(fsim, sim);
        iterations ++;

    }
    return sim.row(0).array();
}

Eigen::SparseMatrix<double> ScrubletWorker::simulate_doublet(const Eigen::SparseMatrix<double> & mat){
    Eigen::ArrayXi E1(this->n_simulations_), E2(this->n_simulations_);
    std::default_random_engine e;
    e.seed(this->random_state_);
    std::uniform_int_distribution<unsigned> u(0, this->n_cells_ - 1);
    for(int i = 0; i < this->n_simulations_; ++i){
        E1[i] = u(e);
        E2[i] = u(e);
    }
    return _Cs col_reordered(mat, E1) + _Cs col_reordered(mat, E2);
    // to do : syntheticdoubletumisampling < 1
}


ScrubletWorker::ScrubletWorker(
    const Eigen::SparseMatrix<int> & mat, 
    unsigned int random_state, 
    double simulate_doublet_ratio, 
    double expected_doublet_ratio, 
    double stdev_doublet_ratio, 
    double synthetic_doublet_umi_downsampling
) :
    mat_(mat), 
    n_cells_(mat.cols()),
    random_state_(random_state), 
    simulate_doublet_ratio_(simulate_doublet_ratio), 
    expected_doublet_ratio_(expected_doublet_ratio), 
    stdev_doublet_ratio_(stdev_doublet_ratio), 
    synthetic_doublet_umi_downsampling_(synthetic_doublet_umi_downsampling)
{
    this->n_neighbors_ = round(sqrt(0.5 * mat.cols()));
    this->n_simulations_ = round(this->n_cells_ * this->simulate_doublet_ratio_);
    this->k_adjust_ = (int)round( this->n_neighbors_ * (1 + this->n_simulations_ / (double)(this->n_cells_)));
}

Eigen::SparseMatrix<double> ScrubletWorker::normalize(){
    Eigen::SparseMatrix<double> norm = this->mat_.cast<double>();
    Eigen::ArrayXd col_sum = _Cs col_sum(this->mat_).cast<double>();
    col_sum /= col_sum.mean();
    for (int k=0; k < norm.outerSize(); ++k){
        for (Eigen::SparseMatrix<double>::InnerIterator it(norm, k); it; ++it){
            it.valueRef() = it.value() / col_sum[k];
        }
    }
    return norm;
};

void ScrubletWorker::filter_gene( Eigen::SparseMatrix<double> & norm){
    Eigen::ArrayXi gene_ix = _Cs sliced(Eigen::ArrayXi(Eigen::ArrayXi::LinSpaced(this->mat_.rows(), 0, this->mat_.rows() - 1)), _Cs row_sum(norm) > 0);
    Eigen::SparseMatrix<double> sub_norm = _Cs row_sliced(norm, _Cs row_sum(norm) > 0);
    Eigen::ArrayXd mu_gene = _Cs row_mean(sub_norm);
    Eigen::ArrayXd var_gene = _Cs row_sum(Eigen::SparseMatrix<double>(sub_norm.cwiseProduct(sub_norm)));
    var_gene -= pow(mu_gene, 2);
    Eigen::ArrayXd FF_gene = var_gene / mu_gene;
    Eigen::ArrayXd data_x = log(mu_gene);
    Eigen::ArrayXd data_y = log(FF_gene / mu_gene);
    _Cs sort_by_first(data_x, data_y);
    int n_bins = 50;
    double fit_percentile = 0.1;
    double dx = (data_x[data_x.size() - 1] - data_x[0]) / n_bins;
    Eigen::ArrayXd x_out = Eigen::ArrayXd::LinSpaced(n_bins, data_x[0] + dx / 2, data_x[data_x.size() - 1] - dx / 2);
    Eigen::ArrayXd y_out = Eigen::ArrayXd::Zero(x_out.size());
    for(int i = 0; i < x_out.size(); ++i){
        Eigen::ArrayX<bool> ind = (data_x >= (x_out[i] - dx / 2)).cwiseProduct (data_x < (x_out[i] + dx / 2));
        if(ind.count() > 0){
            y_out[i] = _Cs linear_percentile(_Cs sliced(data_y, ind), fit_percentile);
        }
        else{
            y_out[i] = y_out[i - 1];
        }
    }
    auto [h, _, b] = _Cs histogram(Eigen::ArrayXd(log(_Cs sliced(FF_gene, mu_gene > 0))), 200);
    int max_ix;
    h.maxCoeff(&max_ix);
    double c = std::max(exp(b[max_ix]), 1.0);
    Eigen::ArrayXd b0(1);
    b0 << 0.1;
    int err_wt = 1;
    auto g_log = [](const Eigen::ArrayXd& i1, double i2, double i3)->Eigen::ArrayXd{
        return log(exp(- i1) * i2 + i3);
    };
    auto err_fun = [&x_out, &y_out, &g_log, &c, &err_wt](const Eigen::ArrayXd& i)->double{
        double ii = i[0];
        return pow(abs(g_log(x_out, c, ii) - y_out), err_wt).sum();
    };
    double b1 = NelderMead(err_fun, b0)[0];
    int min_counts = 3, min_cells = 3;
    double a = c / (1 + b1) - 1;
    Eigen::ArrayXd v_scores = FF_gene / ((1+a)*(1+b1) + b1 * mu_gene);
    Eigen::ArrayX<bool> ix2 = v_scores > 0;

    gene_ix = _Cs sliced(gene_ix, ix2);
    v_scores = _Cs sliced(v_scores, ix2);
    double min_vscore_pctl = 85;
    double min_vscore = _Cs linear_percentile(v_scores, min_vscore_pctl);
    Eigen::ArrayX<bool> ix = (_Cs row_count<double, true, true>(_Cs row_reordered(norm, gene_ix), min_counts) >= min_cells).cwiseProduct(v_scores > min_vscore);
    gene_ix = _Cs sliced(gene_ix, ix2);
    norm = _Cs row_reordered(norm, gene_ix);
    norm = _Cs row_sliced(norm, _Cs col_var_mt( Eigen::SparseMatrix<double>(norm.transpose())) > 0);
};

Eigen::MatrixXd ScrubletWorker::pca(Eigen::SparseMatrix<double> &norm, Eigen::SparseMatrix<double> &sim){
    norm = _Cs normalize(norm, 100000);
    sim = _Cs normalize(sim, 100000);

    Eigen::MatrixXd scaled_norm = _Cs row_scale_mt(norm);
    Eigen::MatrixXd scaled_sim = _Cs row_scale_mt(sim);
    scaled_norm.transposeInPlace();

    auto [U, S, V] = tsvd(&scaled_norm, 30);

    Eigen::MatrixXd pca_orig = U * S.asDiagonal();
    Eigen::MatrixXd pca_sim = scaled_sim.transpose() * V;

    scaled_norm.resize(0, 0);
    scaled_sim.resize(0, 0);

    int n_all = this->n_simulations_ + this->n_cells_;
    Eigen::MatrixXd pca_all(n_all, pca_orig.cols());
    pca_all << pca_orig, pca_sim;
    return pca_all;
};

void ScrubletWorker::get_nearest_neighbors(Eigen::MatrixXi & index, const Eigen::MatrixXd & pca_all){
    int n_all = index.rows();
    AnnoyIndex<int, double, Euclidean, Kiss64Random, AnnoyIndexSingleThreadedBuildPolicy> ann(pca_all.cols());
    for (int i = 0; i < n_all; ++i){
        std::vector<double> trans(pca_all.row(i).begin(), pca_all.row(i).end());
        ann.add_item(i, trans.data());
    }
    ann.set_seed(this->random_state_);
    ann.build(10);

    const int n_cell_all = pca_all.rows();

#pragma omp parallel for
    for(int i = 0; i < n_cell_all; ++i){
        std::vector<int> result;
        std::vector<double> distances;
        ann.get_nns_by_item(i, this->k_adjust_ + 1, -1, &result, &distances);
        for(int j = 1; j < this->k_adjust_ + 1; ++j){
            index(i, j - 1) = result[j];
        }
    };

    ann.unload();
};

void ScrubletWorker::calculate_score(Eigen::MatrixXi & index){
    Eigen::ArrayXi doublet_labels(this->n_simulations_ + this->n_cells_);
    doublet_labels << Eigen::ArrayXi::Constant(this->n_cells_, 0), Eigen::ArrayXi::Constant(this->n_simulations_, 1);
    for(int i = 0; i < index.rows(); ++i){
        for(int j = 0; j < index.cols(); ++j){
            index(i, j) = doublet_labels[index(i, j)];
        }
    }
    Eigen::ArrayXd n_neighbor_sim = index.rowwise().sum().cast<double>();
    Eigen::ArrayXd n_neighbor_orig = this->k_adjust_ - n_neighbor_sim;
    double rho = this->expected_doublet_ratio_;
    double r = this->n_simulations_ / (double)this->n_cells_;
    double k_adjust_d = this->k_adjust_;
    Eigen::ArrayXd q = (n_neighbor_sim + 1) / (k_adjust_d + 2);
    Eigen::ArrayXd Ld = q * rho / r / (1 - rho - q * (1 - rho - rho / r));
    Eigen::ArrayXd doublet_scores_orig = _Cs sliced(Ld, doublet_labels == 0);
    Eigen::ArrayXd doublet_scores_sim = _Cs sliced(Ld, doublet_labels == 1);
    emit x_scrublet_ready(doublet_scores_orig, doublet_scores_sim);
};

void ScrubletWorker::run(){

    Eigen::SparseMatrix<double> norm = this->normalize();

    this->filter_gene(norm);

    Eigen::SparseMatrix<double> sim = this->simulate_doublet(norm);

    Eigen::MatrixXd pca_all = this->pca(norm, sim);

    Eigen::MatrixXi index(this->n_simulations_ + this->n_cells_, this->k_adjust_);

    this->get_nearest_neighbors(index, pca_all);

    this->calculate_score(index);

    G_TASK_END;
}
