#include "SPTree.h"

#include <QDebug>


// Default constructor for SPTree -- build tree, too!
SPTree::SPTree(const double* inp_data, unsigned int n_dimension)
    :
    data_(inp_data)
{

    // Compute mean, width, and height of current map (boundaries of SPTree)
    int loc = 0;
    double mean_y[2] = {};
    double min_y[2] = { DBL_MAX, DBL_MAX};
    double max_y[2] = { -DBL_MAX, -DBL_MAX};

    for (unsigned int n = 0; n < n_dimension; n++) {

        this->boundary_.corner_1_ += inp_data[n * 2 + 0];
        if (inp_data[loc] < min_y[0]) min_y[0] = inp_data[loc];
        if (inp_data[loc] > max_y[0]) max_y[0] = inp_data[loc];

        this->boundary_.corner_2_ += inp_data[n * 2 + 1];
        if (inp_data[loc + 1] < min_y[1]) min_y[1] = inp_data[loc + 1];
        if (inp_data[loc + 1] > max_y[1]) max_y[1] = inp_data[loc + 1];

        loc += 2;
    }

    this->boundary_.corner_1_ /= (double)n_dimension;
    this->boundary_.corner_2_ /= (double)n_dimension;

    // Construct SPTree    
    this->boundary_.width_1_ = std::max(max_y[0] - mean_y[0], mean_y[0] - min_y[0]) + 1e-5;
    this->boundary_.width_2_ = std::max(max_y[1] - mean_y[1], mean_y[1] - min_y[1]) + 1e-5;
    
    for (unsigned int i = 0; i < n_dimension; i++) {

        this->insert(i);
    }
}


// Constructor for SPTree with particular size and parent (do not fill tree)
SPTree::SPTree(
    SPTree* inp_parent,
    const double* inp_data,
    double corner_1,
    double corner_2,
    double width_1,
    double width_2) 
    :
    parent_(inp_parent),
    data_(inp_data),
    boundary_(corner_1, corner_2, width_1, width_2)
{}


// Destructor for SPTree
SPTree::~SPTree()
{

    for (int i = 0; i < 4; ++i) {

        if (this->children_[i] != nullptr) delete this->children_[i];
    }
}


// Insert a point into the SPTree
bool SPTree::insert(unsigned int new_index)
{
    // Ignore objects which do not belong in this quad tree

    const double* point = this->data_ + new_index * 2;
    if (!this->boundary_.contains_point(point[0], point[1]))
        return false;

    // Online update of cumulative size and center-of-mass
    this->cum_size_++;
    double mult1 = (double)(this->cum_size_ - 1) / (double)this->cum_size_;
    double mult2 = 1.0 / (double)this->cum_size_;

        this->center_of_mass_[0] = this->center_of_mass_[0] * mult1 + mult2 * point[0];
        this->center_of_mass_[1] = this->center_of_mass_[1] * mult1 + mult2 * point[1];
    

    // If there is space in this quad tree and it is a leaf, add the object here
    if (this->is_leaf_ && this->size_ < 1) {
        this->index_ = new_index;
        ++this->size_;
        return true;
    }

    // Don't add duplicates for now (this is not very nice)
    bool any_duplicate = false;
    for (int n = 0; n < this->size_; ++n) {
        
        bool duplicate = true;
        
        for (int d = 0; d < 2; ++d) {
        
            if (point[d] != this->data_[this->index_ * 2 + d]) { 
                duplicate = false; 
                break; 
            }
        }
        
        any_duplicate = any_duplicate || duplicate;
    }
    if (any_duplicate) return true;

    // Otherwise, we need to subdivide the current cell
    if (this->is_leaf_) this->subdivide();

    // Find out where the point can be inserted
    for (int i = 0; i < 4; ++i) {

        if (this->children_[i]->insert(new_index)) {
        
            return true;
        }
    }

    // Otherwise, the point cannot be inserted (this should never happen)
    return false;
}


// Create four children which fully divide this cell into four quads of equal area
void SPTree::subdivide() {

    // Create new children
    double new_corner_1{ 0.0 };
    double new_corner_2{ 0.0 };
    double new_width_1{ 0.0 };
    double new_width_2{ 0.0 };

    for (unsigned int i = 0; i < 4; i++) {

        new_width_1 = 0.5 * this->boundary_.width_1_;        

        if (i % 2 == 1) {
            new_corner_1 = this->boundary_.corner_1_ - 0.5 * this->boundary_.width_1_;            
        }
        else {
            new_corner_1 = this->boundary_.corner_1_ + 0.5 * this->boundary_.width_1_;
        }

        new_width_2 = 0.5 * this->boundary_.width_2_;

        if ((i / 2) % 2 == 1) {
            new_corner_2 = this->boundary_.corner_2_ - 0.5 * this->boundary_.width_2_;
        }
        else {
            new_corner_2 = this->boundary_.corner_2_ + 0.5 * this->boundary_.width_2_;
        }

        this->children_[i] = new SPTree(this, this->data_, new_corner_1, new_corner_2, new_width_1, new_width_2);
    }

    // Move existing points to correct children
    for (unsigned int i = 0; i < this->size_; i++) {

        bool success = false;
        for (unsigned int j = 0; j < 4; j++) {
        
            if (!success) success = this->children_[j]->insert(this->index_);
        }

        this->index_ = -1;
    }

    // Empty parent node
    this->size_ = 0;
    this->is_leaf_ = false;
}


// Checks whether the specified tree is correct
bool SPTree::is_correct()
{

    for (unsigned int n = 0; n < this->size_; n++) {
    
        const double* point = this->data_ + this->index_ * 2;
        if (!this->boundary_.contains_point(point[0], point[1])) {
            return false;
        }
    }

    if (!this->is_leaf_) {
        bool correct = true;
        for (int i = 0; i < 4; i++) {
            correct = correct && this->children_[i]->is_correct();
        }
        return correct;
    }
    else {
        return true;
    }
}


// Compute non-edge forces using Barnes-Hut algorithm
double SPTree::compute_non_edge_forces(unsigned int point_index, double theta, Eigen::ArrayXXd& neg_f) const
{
    double result_sum = 0;
    double buff[2];  // make buff local for parallelization

    // Make sure that we spend no time on empty nodes or self-interactions
    if (this->cum_size_ == 0 || (this->is_leaf_ && this->size_ == 1 && this->index_ == point_index)) return result_sum;

    // Compute distance between point and center-of-mass
    double sqdist = 0.0;
    unsigned int ind = point_index * 2;

    for (unsigned int d = 0; d < 2; d++) {
        buff[d] = this->data_[ind + d] - this->center_of_mass_[d];
        sqdist += buff[d] * buff[d];
    }

    // Check whether we can use this node as a "summary"
    double max_width = std::max(this->boundary_.width_1_, this->boundary_.width_2_);

    if (this->is_leaf_ || max_width / sqrt(sqdist) < theta) {

        // Compute and add t-SNE force between point and current node
        sqdist = 1.0 / (1.0 + sqdist);
        double mult = this->cum_size_ * sqdist;
        result_sum += mult;
        mult *= sqdist;
        for (int d = 0; d < 2; d++) {
            neg_f(d, point_index) += mult * buff[d];
        }
    }
    else {

        // Recursively apply Barnes-Hut to children
        for (unsigned int i = 0; i < 4; i++) {
            result_sum += this->children_[i]->compute_non_edge_forces(point_index, theta, neg_f);
        }
    }
    return result_sum;
}


// Computes edge forces
void SPTree::compute_edge_forces(unsigned int* row_P, unsigned int* col_P, double* val_P, Eigen::ArrayXXd& pos_f) const
{
    const int n_vertice = pos_f.cols();

    // Loop over all edges in the graph

#pragma omp parallel for
    for (int n = 0; n < n_vertice; n++) {
        unsigned int ind1 = n * 2;
        for (unsigned int i = row_P[n]; i < row_P[n + 1]; i++) {

            double buff[2] = {}; // make buff local for parallelization

            // Compute pairwise distance and Q-value
            double sqdist = 1.0;
            unsigned int ind2 = col_P[i] * 2;

            for (unsigned int d = 0; d < 2; d++) {
                buff[d] = this->data_[ind1 + d] - this->data_[ind2 + d];
                sqdist += buff[d] * buff[d];
            }

            sqdist = val_P[i] / sqdist;

            // Sum positive force
            for (int d = 0; d < 2; d++) {
                pos_f(d, n) += sqdist * buff[d];
            }
        }
    }
}