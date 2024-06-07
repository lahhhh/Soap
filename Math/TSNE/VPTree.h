#pragma once

#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <queue>
#include <limits>
#include <cmath>
#include "DataPoint.h"

double euclidean_distance(const DataPoint& t1, const DataPoint& t2) {
    
    double distance{ 0.0 };
    int n_dimension = t1.dimensionality();

    for (int d = 0; d < n_dimension; ++d) {
        double d1 = t1.data_[d];
        double d2 = t2.data_[d];
        distance += (d1 - d2) * (d1 - d2);
    };

    return std::sqrt(distance);
}


class VpTree
{
public:

    // Default constructor
    VpTree() : _root(0) {}

    // Destructor
    ~VpTree() {
        delete _root;
    }

    // Function to create a new VpTree from data
    void create(const std::vector<DataPoint>& items) {
        delete _root;
        _items = items;
        _root = buildFromPoints(0, items.size());
    }

    // Function that uses the tree to find the k nearest neighbors of target
    void search(const DataPoint& target, int k, std::vector<DataPoint>* results, std::vector<double>* distances)
    {

        // Use a priority queue to store intermediate results on
        std::priority_queue<HeapItem> heap;

        // Variable that tracks the distance to the farthest point in our results
        double tau = DBL_MAX;

        // Perform the search
        this->search(_root, target, k, heap, &tau);

        // Gather final results
        results->clear(); distances->clear();
        while (!heap.empty()) {
            results->push_back(_items[heap.top().index]);
            distances->push_back(heap.top().dist);
            heap.pop();
        }

        // Results are in reverse order
        std::reverse(results->begin(), results->end());
        std::reverse(distances->begin(), distances->end());
    }

private:
    std::vector<DataPoint> _items;

    // Single node of a VP tree (has a point and radius; left children are closer to point than the radius)
    struct Node
    {
        int index;              // index of point in node
        double threshold;       // radius(?)
        Node* left;             // points closer by than threshold
        Node* right;            // points farther away than threshold

        Node() :
            index(0), threshold(0.), left(0), right(0) {}

        ~Node() {               // destructor
            delete left;
            delete right;
        }
    }*_root;


    // An item on the intermediate result queue
    struct HeapItem {
        HeapItem(int index, double dist) :
            index(index), dist(dist) {}
        int index;
        double dist;
        bool operator<(const HeapItem& o) const {
            return dist < o.dist;
        }
    };

    // Distance comparator for use in std::nth_element
    struct DistanceComparator
    {
        const DataPoint& item;
        DistanceComparator(const DataPoint& item) : item(item) {}
        bool operator()(const DataPoint& a, const DataPoint& b) {
            return euclidean_distance(item, a) < euclidean_distance(item, b);
        }
    };

    // Function that (recursively) fills the tree
    Node* buildFromPoints(int lower, int upper)
    {
        if (upper == lower) {     // indicates that we're done here!
            return NULL;
        }

        // Lower index is center of current node
        Node* node = new Node();
        node->index = lower;

        if (upper - lower > 1) {      // if we did not arrive at leaf yet

            // Choose an arbitrary point and move it to the start

            static std::default_random_engine e(0);
            static std::uniform_real_distribution<double> u(0.0, 1.0);

            int i = (int)(u(e) * (upper - lower - 1)) + lower;
            std::swap(_items[lower], _items[i]);

            // Partition around the median distance
            int median = (upper + lower) / 2;
            std::nth_element(_items.begin() + lower + 1,
                _items.begin() + median,
                _items.begin() + upper,
                DistanceComparator(_items[lower]));

            // Threshold of the new node will be the distance to the median
            node->threshold = euclidean_distance(_items[lower], _items[median]);

            // Recursively build tree
            node->index = lower;
            node->left = buildFromPoints(lower + 1, median);
            node->right = buildFromPoints(median, upper);
        }

        // Return result
        return node;
    }

    // Helper function that searches the tree    
    void search(Node* node, const DataPoint& target, unsigned int k, std::priority_queue<HeapItem>& heap, double* ptau)
    {
        if (node == NULL) return;     // indicates that we're done here

        // Compute distance between target and current node
        double dist = euclidean_distance(_items[node->index], target);

        // If current node within radius tau
        if (dist < (*ptau)) {
            if (heap.size() == k) heap.pop();                // remove furthest node from result list (if we already have k results)
            heap.push(HeapItem(node->index, dist));         // add current node to result list
            if (heap.size() == k) *ptau = heap.top().dist;   // update value of tau (farthest point in result list)
        }

        // Return if we arrived at a leaf
        if (node->left == NULL && node->right == NULL) {
            return;
        }

        // If the target lies within the radius of ball
        if (dist < node->threshold) {
            if (dist - (*ptau) <= node->threshold) {         // if there can still be neighbors inside the ball, recursively search left child first
                search(node->left, target, k, heap, ptau);
            }

            if (dist + (*ptau) >= node->threshold) {         // if there can still be neighbors outside the ball, recursively search right child
                search(node->right, target, k, heap, ptau);
            }

            // If the target lies outsize the radius of the ball
        }
        else {
            if (dist + (*ptau) >= node->threshold) {         // if there can still be neighbors outside the ball, recursively search right child first
                search(node->right, target, k, heap, ptau);
            }

            if (dist - (*ptau) <= node->threshold) {         // if there can still be neighbors inside the ball, recursively search left child
                search(node->left, target, k, heap, ptau);
            }
        }
    }
};
