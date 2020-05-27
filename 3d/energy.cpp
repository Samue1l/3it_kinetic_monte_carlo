
#include "energy.h"

std::array<int, 5> get_energy_neighbourhood(int p) {
    /*
    Return the neighbourhood used to calculate the energy of a particule
    */
    std::array<int, 5> r_arr = {p-2,p-1, p, p+1,p+2};
    return r_arr;
}


std::array<int, 3> get_neighbourhood(int p) {
    /*
    Return the general jumping neighbourhood of a particule
    */
    std::array<int, 3> r_arr = {p-1, p, p+1};
    return r_arr;
}


float compute_energy(unsigned char* lat, int x, int y, int z, const int mx, const int my, const int mz) {
    /*
    Compute the energy of a particule
    */
    float e = 0.0f;
    float Eb = 1.0; // eV
    float anisotropy_coeff;
    
    auto x_neighbour = get_energy_neighbourhood(x);
    auto y_neighbour = get_energy_neighbourhood(y);
    auto z_neighbour = get_energy_neighbourhood(z);

    int i =0;
    int v =0;
    for(int ex : x_neighbour){
        ex = (ex+mx)%mx;
        for(int ey : y_neighbour) {
            ey= (ey+my)%my;
            for(int ez : z_neighbour){
                ez = (ez+mz)%mz;
                v = lat[ex*my*mz+ey*mz+ez]>0?1:0;

                e = e + 1-v;
                i = i+v;
            }
            //e = e*(1.2-0.2*x/mx);
        }
    }
    
    return e*Eb;
}
