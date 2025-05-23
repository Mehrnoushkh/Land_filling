#include <vector>
#include <cmath>
#include <iostream>
#include "functions.hpp"
#include "laplace.hpp"
#include <netcdf.h>

const double R_E = 6371000.0;
const double PI = 3.141592653589793;
// padding function
vector<vector<double>> pad_top(const vector<vector<double>>& mat, int rows_to_add, double pad_value = NAN) {
    size_t cols = mat[0].size();
    vector<vector<double>> padded(rows_to_add, vector<double>(cols, pad_value));
    padded.insert(padded.end(), mat.begin(), mat.end());
    return padded;
}

vector<double> concat_vectors(const vector<double>& top, const vector<double>& bottom) {
    vector<double> result = top;
    result.insert(result.end(), bottom.begin(), bottom.end());
    return result;
}

int main() {
    dataset source_data;

    source_data.load_time("time", "SST_1day.nc");
    source_data.load_depth("depth", "SST_1day.nc");
    source_data.load_latitude("latitude", "SST_1day.nc");
    source_data.load_longitude("longitude", "SST_1day.nc");
    size_t Nlat = source_data.Nlat;
    size_t Nlon = source_data.Nlon;

    vector<vector<double>> SST(Nlat, vector<double>(Nlon));
    vector<vector<double>> padded_SST = pad_top(SST, 39);

    vector<double> new_lat(39);
    for (int i = 0; i < 39; ++i)
        new_lat[i] = (-89.75 + i * 0.25) * PI / 180.0;

    // concatenate part, for latitude
    vector<double> extended_lat = concat_vectors(new_lat, source_data.latitude);
    //mask data
    size_t Nlat_extended = padded_SST.size();
    size_t Nlon = padded_SST[0].size();

    vector<vector<bool>> mask(Nlat_extended, vector<bool>(Nlon));

    for (size_t i = 0; i < Nlat_extended; ++i) {
        for (size_t j = 0; j < Nlon; ++j) {
            mask[i][j] = std::isnan(SST[i][j]); 
    }
    }
	
    double dlat = Nlat_extended[1] - Nlat_extended[0];
    double dlon = source_data.longitude[1] - source_data.longitude[0];
    
    vector<vector<double>> LAT(Nlat_extended, vector<double>(Nlon));   
    vector<vector<double>> LON(Nlat_extended, vector<double>(Nlon));
    for (size_t i = 0; i < Nlat_extended; ++i) {
        for (size_t j = 0;j < Nlon; ++j) {
            LAT[i][j] = extended_lat[i];              
            LON[i][j] = source_data.longitude[j];    
         }
     }

	// dA, dx, dy:
    vector<vector<double>> dA(Nlat_extended, vector<double>(Nlon));
    vector<vector<double>> dx(Nlat_extended, vector<double>(Nlon));
    vector<vector<double>> dy(Nlat_extended, vector<double>(Nlon, R_E * dlat));  // dy is constant in x
    for (size_t i = 0; i < Nlat_extended; ++i) {
        for (size_t j = 0; j < Nlon; ++j) {
            double cos_lat = cos(LAT[i][j]);
            dA[i][j] = R_E * R_E * dlat * dlon * cos_lat;
            dx[i][j] = R_E * dlon * cos_lat;
    }
}
    std::vector<double> zonal_mean(Nlat_extended);

    for (size_t i = 0; i < Nlat_extended; ++i) {
        double sum = 0.0;
        int count = 0;
        for (size_t j = 0; j < Nlon; ++j) {
            if (!mask[i][j]) {  // ocean?
               sum += padded_SST[i][j];
               count++;
           }
       }
      zonal_mean[i] = (count > 0) ? (sum / count) : NAN;
      } 
    int SP_cutoff = 0
    for  (size_t i=0; i<Nlat_extended; ++i){
	    bool all_land = True ;
	    for (size_t j = 0; j < Nlon; ++j) {    
	        if (!mask[i][j]){ //ocean areas
			all_land = false;
                        break;
		}
	    }
	    if (all_land){
		    sp_cutoff = i;
		    break;
	    }
    }
    std::vector<double> smooth_mean(Nlat_extended);
    for (int i = 0; i < SP_cutoff; ++i) {
    zonal_mean[i] = zonal_mean[SP_cutoff];
    }
    
    for (size_t i = 0; i < Nlat_extended; ++i) {
	    int LB= std::max(0, static_cast<int> i-8);
            int UB = std::min(static_cast<int>(Nlat_extended - 1), static_cast<int>(i) + 8);
            double sum=0.00;
            int count = 0 ;		    
	    for (int k = LB; k <= UB; ++k) {
		    if (!std::isnan(zonal_mean[k])) {
                       sum += zonal_mean[k];
                       count++;
                    }  
             }
            smooth_mean[i] = (count > 0) ? (sum / count) : NAN;
       }
    std::vector<std::vector<double>> zonal_matrix(Nlat_extended,std::vector<double>(Nlon));

    for (size_t i = 0; i < Nlat_extended; ++i) {
        for (size_t j = 0; j < Nlon; ++j) {
           zonal_matrix[i][j] = smooth_mean[i];
        }
    }
    


    std::vector<std::vector<double>> Filled_SST = padded_SST;
    for (size_t i = 0; i < Nlat_extended; ++i) {
       for (size_t j = 0; j < Nlon; ++j) {
	       if(mask[i][j])
	       {
		       Filled_SST[i][j] = zonal_matrix[i][j];
	       }
	}
    }
    std::vector<std::vector<double>> cosLat(Nlat_extended, std::vector<double>(Nlon));

    for (size_t i = 0; i < Nlat_extended; ++i) {
        double cos_phi = std::cos(extended_lat[i]);
        for (size_t j = 0; j < Nlon; ++j) {
            cosLat[i][j] = cos_phi;  // same value along row
         }
     }

    
    int Nsteps = 3000;
    int Iprint = 300;
    double dx_val = dx[0][0];
    double dy_val = dy[0][0];
    double dt = std::pow(std::min(dx_val, dy_val), 2);

       ////  double dt = std::pow(std::min(dx, dy), 2);
    for (int Istep = 0; Istep <= Nsteps; ++Istep)
    	    double CFL_val = (Istep < 100) ? 0.01 : 0.4;
	    std::vector<std::vector<double>> RHS = Laplace(Filled_SST, cosLat, dlat, dlon);
	    for (size_t i = 0; i < Nlat_extended; ++i) {
     		   for (size_t j = 0; j < Nlon; ++j) {
            		if (!mask[i][j]) {
                		RHS[i][j] = 0.0;
            		}
         	   }
            }
    for (size_t i = 0; i < Nlat_extended; ++i) {
                   for (size_t j = 0; j < Nlon; ++j) {
			   Filled_SST[i][j] += CFL_val * dt * RHS[i][j];
		   }
    }

    for (size_t i = 0; i < Nlat_extended; ++i) {
        for (size_t j = 0; j < Nlon; ++j) {
            if (std::abs(Filled_SST[i][j]) > 100) {
                std::cout << "Unstable value at step " << Istep << ", stopping.\n";
                return 1;
            }
        }
    }

    if (Istep % Iprint == 0) {
        std::cout << "Step " << Istep << " complete.\n";
    }

// write data to output file

    int Nlat = Filled_SST.size();
    int Nlon = Filled_SST[0].size();

    const std::vector<double>& lat_out = extended_lat;
    const std::vector<double>& lon_out = source_data.longitude;

    std::vector<double> flat_data(Nlat * Nlon);
    for (int i = 0; i < Nlat; ++i)
        for (int j = 0; j < Nlon; ++j)
            flat_data[i * Nlon + j] = Filled_SST[i][j];

    int ncid;
    nc_create("filled_lands_file.nc", NC_NETCDF4, &ncid);

    int lat_dimid, lon_dimid;
    nc_def_dim(ncid, "lat", Nlat, &lat_dimid);
    nc_def_dim(ncid, "lon", Nlon, &lon_dimid);

    int lat_varid, lon_varid;
    nc_def_var(ncid, "lat", NC_DOUBLE, 1, &lat_dimid, &lat_varid);
    nc_def_var(ncid, "lon", NC_DOUBLE, 1, &lon_dimid, &lon_varid);
    int dimids[2] = {lat_dimid, lon_dimid};
    int sst_varid;
    nc_def_var(ncid, "Filled_SST", NC_DOUBLE, 2, dimids, &sst_varid);

    nc_enddef(ncid);
    nc_put_var_double(ncid, lat_varid, lat_out.data());
    nc_put_var_double(ncid, lon_varid, lon_out.data());

    nc_put_var_double(ncid, sst_varid, flat_data.data());
    nc_close(ncid);

    return 0
}

