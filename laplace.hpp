#ifndef LAPLACE_HPP
#define LAPLACE_HPP

#include <vector>

const double R_E = 6371000.0;
std::vector<std::vector<double>> Laplace(
    const std::vector<std::vector<double>>& f,
    const std::vector<std::vector<double>>& cosLat,
    double dlat,
    double dlon
);

#endif // LAPLACE_HPP

