#include "Tsne.h"

#include <random>
#include <QDebug>

#include "sptree.h"
#include "vptree.h"
#include "tsne.h"

Tsne::Tsne(
	unsigned int random_state,
	double perplexity,
	double theta,
	int max_iter,
	int stop_lying_iter,
	int mom_switch_iter,
	double momentum,
	double final_momentum,
	double eta,
	double exaggeration_factor
)
	:
	random_state_(random_state),
	perplexity_(perplexity),
	theta_(theta),
	momentum_(momentum),
	final_momentum_(final_momentum),
	eta_(eta),
	exaggeration_factor_(exaggeration_factor),
	max_iter_(max_iter),
	stop_lying_iter_(stop_lying_iter),
	mom_switch_iter_(mom_switch_iter),
	exact_(theta == 0.0)
{}

// Perform t-SNE
bool Tsne::run(Eigen::MatrixXd* mat, Eigen::MatrixXd& embedding) {

	unsigned int n_vertice = mat->cols();
	if (n_vertice - 1 < 3 * this->perplexity_) {
		return false;
	}

	unsigned int n_dimension = mat->rows();

	// Compute input similarities for exact t-SNE
	if (this->exact_) {
		// Compute similarities
		compute_gaussian_perplexity(*mat);

		// Symmetrize input similarities
		for (int n = 0; n < n_vertice; n++) {
			for (int m = n + 1; m < n_vertice; m++) {
				P(m, n) += P(n, m);
				P(n, m) = P(m, n);
			}
		}

		P /= P.sum();
	}

	// Compute input similarities for approximate t-SNE
	else {

		int K = 3 * this->perplexity_;

		// Compute asymmetric pairwise input similarities
		compute_gaussian_perplexity(*mat, K);

		// Symmetrize input similarities
		symmetrize_matrix(n_vertice);
		double sum_P = 0.0;
		for (unsigned int i = 0; i < row_P[n_vertice]; i++) sum_P += val_P[i];
		for (unsigned int i = 0; i < row_P[n_vertice]; i++) val_P[i] /= sum_P;
	}

	train_iterations(n_vertice, embedding);

	return true;
}

// Perform main training loop
void Tsne::train_iterations(int n_vertice, Eigen::MatrixXd& output) {

	// Initialize solution (randomly)
	srand(this->random_state_);
	output = Eigen::MatrixXd::Random(2, n_vertice);

	// Allocate some memory
	Eigen::ArrayXXd d_output = Eigen::MatrixXd::Zero(2, n_vertice);
	Eigen::ArrayXXd u_output = Eigen::MatrixXd::Zero(2, n_vertice);
	Eigen::ArrayXXd gains = Eigen::MatrixXd::Zero(2, n_vertice);

	// Lie about the P-values
	if (this->exact_) {
		P *= this->exaggeration_factor_;
	}
	else {
		for (int i = 0; i < row_P[n_vertice]; i++) {
			val_P[i] *= this->exaggeration_factor_;
		}
	}

	for (int iter = 0; iter < this->max_iter_; iter++) {

		// Stop lying about the P-values after a while, and switch momentum
		if (iter == this->stop_lying_iter_) {

			if (this->exact_) {
					P /= this->exaggeration_factor_;
			}
			else {
				for (unsigned int i = 0; i < row_P[n_vertice]; i++) {
					val_P[i] /= this->exaggeration_factor_;
				}
			}
		}

		if (iter == this->mom_switch_iter_) {
			this->momentum_ = this->final_momentum_;
		}

		// Compute (approximate) gradient
		if (this->exact_) {
			compute_exact_gradient(output, n_vertice, d_output);
		}
		else {
			compute_gradient(output, d_output);
		}

		// Update gains
		for (unsigned int i = 0; i < n_vertice; i++) {
			gains(0, i) = (sign_tsne(d_output(0, i)) != sign_tsne(u_output(0, i))) ? (gains(0, i) + .2) : (gains(0, i) * .8);
			gains(1, i) = (sign_tsne(d_output(1, i)) != sign_tsne(u_output(1, i))) ? (gains(1, i) + .2) : (gains(1, i) * .8);
		}

		gains = gains.cwiseMax(0.01);

		// Perform gradient update (with momentum and gains)

		u_output = this->momentum_ * u_output - this->eta_ * gains * d_output;

		output.array() += u_output;

		// Make solution zero-mean
		zero_mean(output);
	}

	return;
}

// Compute gradient of the t-SNE cost function (using Barnes-Hut algorithm)
void Tsne::compute_gradient(const Eigen::MatrixXd& embedding, Eigen::ArrayXXd& d_output)
{
	const int n_vertice = embedding.cols();

	// Construct space-partitioning tree on current map
	SPTree* tree = new SPTree(embedding.data(), n_vertice);

	// Compute all terms required for t-SNE gradient
	Eigen::ArrayXXd pos_f = Eigen::ArrayXXd::Zero(2, n_vertice);
	Eigen::ArrayXXd neg_f = Eigen::ArrayXXd::Zero(2, n_vertice);
	//double* neg_f = (double*)calloc(2 * n_vertice, sizeof(double));

	tree->compute_edge_forces(this->row_P.data(), this->col_P.data(), this->val_P.data(), pos_f);

	// Storing the output to sum in single-threaded mode; avoid randomness in rounding errors.
	Eigen::ArrayXd output = Eigen::ArrayXd::Zero(n_vertice);

#pragma omp parallel for 
	for (int n = 0; n < n_vertice; n++) {
		output[n] = tree->compute_non_edge_forces(n, this->theta_, neg_f);
	}

	double sum_Q = output.sum();

	// Compute final t-SNE gradient

	d_output = pos_f - (neg_f / sum_Q);

	delete tree;
}

// Compute gradient of the t-SNE cost function (exact)
void Tsne::compute_exact_gradient(Eigen::MatrixXd& output, int n_vertice, Eigen::ArrayXXd& d_output) {

	// Make sure the current gradient contains zeros
	d_output.setZero();

	// Compute the squared Euclidean distance matrix
	Eigen::ArrayXXd distance = Eigen::ArrayXXd::Zero(n_vertice, n_vertice);
	compute_squared_euclidean_distance(output, distance);

	// Compute Q-matrix and normalization sum
	Eigen::ArrayXXd Q = 1 / (1 + distance);
	double sum_Q = Q.sum() - n_vertice;

	// Perform the computation of the gradient
	for (int n = 0; n < n_vertice; n++) {
		for (int m = 0; m < n_vertice; m++) {
			if (n != m) {
				double mult = (P(m, n) - (Q(m, n) / sum_Q)) * Q(m, n);
				for (int d = 0; d < 2; d++) {
					d_output(d, n) += (output(d, n) - output(d, m)) * mult;
				}
			}
		}
	}
}

// Compute input similarities with a fixed perplexity
void Tsne::compute_gaussian_perplexity(const Eigen::MatrixXd& data) {

	const int n_vertice = data.cols();
	const int n_dimension = data.rows();

	P = Eigen::MatrixXd::Zero(n_vertice, n_vertice);

	// Compute the squared Euclidean distance matrix
	Eigen::ArrayXXd distance = Eigen::ArrayXXd::Zero(n_vertice, n_vertice);
	compute_squared_euclidean_distance(data, distance);

	// Compute the Gaussian kernel row by row
	for (int n = 0; n < n_vertice; n++) {

		// Initialize some variables
		bool found = false;
		double beta = 1.0;
		double min_beta = -DBL_MAX;
		double max_beta = DBL_MAX;
		double tol = 1e-5;
		double sum_P;

		// Iterate until we found a good perplexity
		int iter = 0;
		while (!found && iter < 200) {

			// Compute Gaussian kernel row
			for (int m = 0; m < n_vertice; m++) P(m, n) = exp(-beta * distance(m, n));
			P(n, n) = DBL_MIN;

			// Compute entropy of current row
			sum_P = DBL_MIN;
			sum_P += P.col(n).sum();

			double H = beta * (distance.col(n).array() * P.col(n).array()).sum();
			H = (H / sum_P) + log(sum_P);

			// Evaluate whether the entropy is within the tolerance level
			double Hdiff = H - log(this->perplexity_);
			if (Hdiff < tol && -Hdiff < tol) {
				found = true;
			}
			else {
				if (Hdiff > 0) {
					min_beta = beta;
					if (max_beta == DBL_MAX || max_beta == -DBL_MAX)
						beta *= 2.0;
					else
						beta = (beta + max_beta) / 2.0;
				}
				else {
					max_beta = beta;
					if (min_beta == -DBL_MAX || min_beta == DBL_MAX)
						beta /= 2.0;
					else
						beta = (beta + min_beta) / 2.0;
				}
			}

			// Update iteration counter
			iter++;
		}

		// Row normalize P
		P.col(n).array() /= sum_P;
	}
}


// Compute input similarities with a fixed perplexity using ball trees (this function allocates memory another function should free)

void Tsne::compute_gaussian_perplexity(const Eigen::MatrixXd& data, int K) {

	const int n_vertice = data.cols();
	const int n_dimension = data.rows();

	// Allocate the memory we need
	setup_approximate_memory(n_vertice, K);

	// Build ball tree on data set
	VpTree* tree = new VpTree();
	std::vector<DataPoint> obj_X(n_vertice);
	for (unsigned int n = 0; n < n_vertice; n++) obj_X[n] = DataPoint(n_dimension, n, data);
	tree->create(obj_X);

	// Loop over all points to find nearest neighbors

#pragma omp parallel for 
	for (int n = 0; n < n_vertice; n++) {

		std::vector<DataPoint> indices;
		std::vector<double> distances;
		indices.reserve(K + 1);
		distances.reserve(K + 1);

		// Find nearest neighbors
		tree->search(obj_X[n], K + 1, &indices, &distances);

		double* cur_P = val_P.data() + row_P[n];
		compute_probabilities(this->perplexity_, K, distances.data() + 1, cur_P); // +1 to avoid self.

		unsigned int* cur_col_P = col_P.data() + row_P[n];
		for (int m = 0; m < K; ++m) {
			cur_col_P[m] = indices[m + 1].index(); // +1 to avoid self.
		}
	}

	delete tree;
}

void Tsne::setup_approximate_memory(unsigned int n_vertice, int K) {
	row_P.resize(n_vertice + 1);
	col_P.resize(n_vertice * K);
	val_P.resize(n_vertice * K);
	row_P[0] = 0;
	for (unsigned int n = 0; n < n_vertice; n++) row_P[n + 1] = row_P[n] + K;
	return;
}

void Tsne::compute_probabilities(const double perplexity, const int K, const double* distances, double* cur_P) {

	// Initialize some variables for binary search
	bool found = false;
	double beta = 1.0;
	double min_beta = -DBL_MAX;
	double max_beta = DBL_MAX;
	double tol = 1e-5;

	// Iterate until we found a good perplexity
	int iter = 0; double sum_P;
	while (!found && iter < 200) {

		// Compute Gaussian kernel row
		for (int m = 0; m < K; m++) cur_P[m] = exp(-beta * distances[m] * distances[m]);

		// Compute entropy of current row
		sum_P = DBL_MIN;
		for (int m = 0; m < K; m++) sum_P += cur_P[m];
		double H = 0.0;
		for (int m = 0; m < K; m++) H += beta * (distances[m] * distances[m] * cur_P[m]);
		H = (H / sum_P) + log(sum_P);

		// Evaluate whether the entropy is within the tolerance level
		double Hdiff = H - log(perplexity);
		if (Hdiff < tol && -Hdiff < tol) {
			found = true;
		}
		else {
			if (Hdiff > 0) {
				min_beta = beta;
				if (max_beta == DBL_MAX || max_beta == -DBL_MAX)
					beta *= 2.0;
				else
					beta = (beta + max_beta) / 2.0;
			}
			else {
				max_beta = beta;
				if (min_beta == -DBL_MAX || min_beta == DBL_MAX)
					beta /= 2.0;
				else
					beta = (beta + min_beta) / 2.0;
			}
		}

		// Update iteration counter
		iter++;
	}

	// Row-normalize current row of P.
	for (int m = 0; m < K; m++) {
		cur_P[m] /= sum_P;
	}
	return;
}

void Tsne::symmetrize_matrix(unsigned int n_vertice) {
	// Count number of elements and row counts of symmetric matrix
	Eigen::ArrayXi row_counts = Eigen::ArrayXi::Zero(n_vertice);

	for (unsigned int n = 0; n < n_vertice; n++) {
		for (unsigned int i = row_P[n]; i < row_P[n + 1]; i++) {

			// Check whether element (col_P[i], n) is present
			bool present = false;
			for (unsigned int m = row_P[col_P[i]]; m < row_P[col_P[i] + 1]; m++) {
				if (col_P[m] == n) present = true;
			}
			if (present) row_counts[n]++;
			else {
				row_counts[n]++;
				row_counts[col_P[i]]++;
			}
		}
	}
	int no_elem = 0;
	for (unsigned int n = 0; n < n_vertice; n++) no_elem += row_counts[n];

	// Allocate memory for symmetrized matrix
	std::vector<unsigned int> sym_row_P(n_vertice + 1), sym_col_P(no_elem);
	std::vector<double> sym_val_P(no_elem);

	// Construct new row indices for symmetric matrix
	sym_row_P[0] = 0;
	for (unsigned int n = 0; n < n_vertice; n++) sym_row_P[n + 1] = sym_row_P[n] + row_counts[n];

	// Fill the result matrix
	Eigen::ArrayXi offset = Eigen::ArrayXi::Zero(n_vertice);
	for (unsigned int n = 0; n < n_vertice; n++) {
		for (unsigned int i = row_P[n]; i < row_P[n + 1]; i++) {                                  // considering element(n, col_P[i])

			// Check whether element (col_P[i], n) is present
			bool present = false;
			for (unsigned int m = row_P[col_P[i]]; m < row_P[col_P[i] + 1]; m++) {
				if (col_P[m] == n) {
					present = true;
					if (n <= col_P[i]) {                                                 // make sure we do not add elements twice
						sym_col_P[sym_row_P[n] + offset[n]] = col_P[i];
						sym_col_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = n;
						sym_val_P[sym_row_P[n] + offset[n]] = val_P[i] + val_P[m];
						sym_val_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = val_P[i] + val_P[m];
					}
				}
			}

			// If (col_P[i], n) is not present, there is no addition involved
			if (!present) {
				sym_col_P[sym_row_P[n] + offset[n]] = col_P[i];
				sym_col_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = n;
				sym_val_P[sym_row_P[n] + offset[n]] = val_P[i];
				sym_val_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = val_P[i];
			}

			// Update offsets
			if (!present || (present && n <= col_P[i])) {
				offset[n]++;
				if (col_P[i] != n) offset[col_P[i]]++;
			}
		}
	}

	// Divide the result by two
	for (int i = 0; i < no_elem; i++) sym_val_P[i] /= 2.0;

	// Return symmetrized matrices
	row_P.swap(sym_row_P);
	col_P.swap(sym_col_P);
	val_P.swap(sym_val_P);
}

// Compute squared Euclidean distance matrix (n_verticeo BLAS)
void Tsne::compute_squared_euclidean_distance(const Eigen::MatrixXd& data, Eigen::ArrayXXd& distance) {

	const int n_vertice = data.cols();

#pragma omp parallel for
	for (int i = 0; i < n_vertice; ++i) {
		for (int j = 0; j < i; ++j) {
			distance(i, j) = std::sqrt((data.col(i).array() - data.col(j).array()).square().sum());
		}
	};

	Eigen::ArrayXXd transposed = distance.transpose();
	distance += transposed;
}


// Makes data zero-mean
void Tsne::zero_mean(Eigen::MatrixXd& output) {

	// Compute data mean
	Eigen::ArrayXd mean = output.array().rowwise().mean();

	output.array().colwise() -= mean;

}
