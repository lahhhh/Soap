#pragma once

#include "Identifier.h"

struct Cell {

	Cell() = default;

	Cell(
		double corner_1,
		double corner_2,
		double width_1,
		double width_2) :
		corner_1_(corner_1),
		corner_2_(corner_2),
		width_1_(width_1),
		width_2_(width_2)
	{};

	bool contains_point(double p1, double p2) const {

		if (this->corner_1_ - this->width_1_ > p1 ||
			this->corner_1_ + this->width_1_ < p1 ||
			this->corner_2_ - this->width_2_ > p2 ||
			this->corner_2_ + this->width_2_ < p2
			)
		{
			return false;
		}
		else {
			return true;
		}
	}

	double corner_1_;
	double corner_2_;
	double width_1_;
	double width_2_;
};


class SPTree
{

private:

	// Properties of this node in the tree
	SPTree* parent_{nullptr};
	unsigned int dimension_;
	bool is_leaf_{true};
	unsigned int size_{0};
	unsigned int cum_size_{0};

	// Axis-aligned bounding box stored as a center with half-dimensions to represent the boundaries of this quad tree
	Cell boundary_;

	// Indices in this space-partitioning tree node, corresponding center-of-mass, and list of all children
	const double* data_;
	double center_of_mass_[2]{ 0.0 };
	unsigned int index_;

	// Children
	SPTree* children_[4]{ nullptr };

public:
	SPTree(const double* inp_data, unsigned int n_dimension);
	SPTree(
		SPTree* inp_parent, 
		const double* inp_data, 
		double corner_1,
		double corner_2,
		double width_1,
		double width_2);
	~SPTree();

	bool insert(unsigned int new_index);
	void subdivide();
	bool is_correct();
	double compute_non_edge_forces(unsigned int point_index, double theta, Eigen::ArrayXXd& neg_f) const;
	void compute_edge_forces(unsigned int* row_P, unsigned int* col_P, double* val_P, Eigen::ArrayXXd& pos_f) const;
};

