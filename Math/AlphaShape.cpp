#include "AlphaShape.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Alpha_shape_vertex_base_2<K> Vb;
typedef CGAL::Alpha_shape_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds> Triangulation_2;
typedef CGAL::Alpha_shape_2<Triangulation_2> Alpha_shape_2;
typedef Alpha_shape_2::Alpha_shape_edges_iterator Alpha_shape_edges_iterator;

std::pair<QVector<double>, QVector<double>>
alpha_shape(
    const Eigen::ArrayXd& x,
    const Eigen::ArrayXd& y,
    double alpha
) {
    int n_point = x.size();
    std::vector<Point> points(n_point);

    for (int i = 0; i < n_point; ++i) {
        points[i] = Point(x[i], y[i]);
    }

    Alpha_shape_2 A(points.begin(), points.end(), alpha, Alpha_shape_2::GENERAL);

    QVector<double> s_x, s_y;

    for (auto it = A.Alpha_shape_vertices_begin(); it != A.Alpha_shape_vertices_end(); ++it) {
        auto&& p = (*it)->point();
        s_x << p.x();
        s_y << p.y();
    }

    return std::make_pair(s_x, s_y);
};