
#include "energy.h"

std::array<int, 7> get_energy_neighbourhood(int p) {
    /*
    Return the neighbourhood used to calculate the energy of a particule
    */
//    std::array<int, 9> r_arr = {p-4, p-3, p-2, p-1, p, p+1, p+2, p+3, p+4};
    std::array<int, 7> r_arr = {p-3, p-2,p-1, p, p+1, p+2, p+3};
//    std::array<int, 5> r_arr = {p-2,p-1, p, p+1, p+2};
    return r_arr;
}


std::array<int, 3> get_neighbourhood(int p) {
    /*
    Return the general jumping neighbourhood of a particule
    */
    std::array<int, 3> r_arr = {p-1, p, p+1};
    return r_arr;
}
/*
bool is_in_energy_zone(int x, int y) {
    bool is_in;
    int xcenter = 16;
    int ycenter = 16;
    int radius = 4;

    if( pow(x-xcenter,2) + pow(y-ycenter,2) <= pow(radius,2) ){
        is_in=true;
    } else {
        is_in=false;
    }

    return is_in;

}

*/

float compute_energy(unsigned char* lat, int x, int y, int z, const int mx, const int my, const int mz) {
    /*
    Compute the energy of a particule
    */
    float e = 0.0f;
    float Eb = 0.3*1.60217e-19; // eV
    float E0 = 0.25*1.60217e-19;
    float anisotropy_coeff;

    auto x_neighbour = get_energy_neighbourhood(x);

    const float Nmax = pow(x_neighbour.size(), 3)-1;

    auto y_neighbour = get_energy_neighbourhood(y);
    auto z_neighbour = get_energy_neighbourhood(z);

//    const int xcenter = 16;
//    const int ycenter = 16;
//    const int radius = 4;
//
//    int d = sqrt(pow(x-xcenter,2) + pow(y-ycenter,2));
//    if (d >=radius && z>3) {
//        Eb = Eb*(1+2*d*d/(mx*mx));
//    }

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
//    e= 125*e/nmax;
    return 0.5*(124*e*Eb/Nmax-E0);
}
