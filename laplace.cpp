#include <vector>
#include <cmath>
const double R_E = 6371000.0;  // I defined it in main code as well

std::vector<std::vector<double>> Laplace(
    const std::vector<std::vector<double>>& f,
    const std::vector<std::vector<double>>& cosLat,
    double dlat,
    double dlon
) {
    int Nlat = f.size();
    int Nlon = f[0].size();
    double R2 = R_E * R_E;

    std::vector<std::vector<double>> lap(Nlat, std::vector<double>(Nlon, 0.0));

    for (int i = 1; i < Nlat - 1; ++i) {
        for (int j = 0; j < Nlon; ++j) {
            int jm1 = (j - 1 + Nlon) % Nlon;  // this is periodic
            int jp1 = (j + 1) % Nlon;

            // second derivative in lon
            double d2fdlon2 = (f[i][jp1] - 2.0 * f[i][j] + f[i][jm1]) / (dlon * dlon);

            // first derivative in latitude of cos(lat) * df/dlat
            double dfdlat_p = (f[i + 1][j] - f[i][j]) / dlat;
            double dfdlat_m = (f[i][j] - f[i - 1][j]) / dlat;
            double term_lat = (cosLat[i + 1][j] * dfdlat_p - cosLat[i - 1][j] * dfdlat_m) / (2.0 * dlat);

            double cos2 = cosLat[i][j] * cosLat[i][j];
            lap[i][j] = (d2fdlon2 / cos2 + term_lat / cosLat[i][j]) / R2;
        }
    }

    return lap;
}
