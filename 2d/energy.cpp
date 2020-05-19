
#include "energy.h"
std::array<int, 5> get_energy_neighbourhood(int p) {
    /*
    Return the neighbourhood used to calculate the energy of a particule
    */
    std::array<int, 5> r_arr = {p-2,p-1, p, p+1,p+2};
    return r_arr;
}
float e_rest = 0.5f;
float e_100 = 1.2f;
float e_111 = 1.2f;
float anisotropy[25] = {
                                        e_100, e_rest, e_111, e_rest, e_100, 
                                        e_rest, e_100, e_111, e_100, e_rest, 
                                        e_111, e_111, 0.0f, e_111, e_111,
                                        e_rest, e_100, e_111, e_100, e_rest,
                                        e_100, e_rest, e_111, e_rest, e_100};

std::array<int, 3> get_neighbourhood(int p) {
    /*
    Return the general jumping neighbourhood of a particule
    */
    std::array<int, 3> r_arr = {p-1, p, p+1};
    return r_arr;
}

float compute_energy(unsigned char* lat, int x, int y, const int mx, const int my) {
    /*
    Compute the energy of a particule
    */
    float e = 0.0f;
    float Eb = 1.0f;
    float anisotropy_coeff;
    
    auto x_neighbour = get_energy_neighbourhood(x);
    auto y_neighbour = get_energy_neighbourhood(y);

    int i =0;
    for(int ex : x_neighbour){
        ex = (ex+mx)%mx;
        for(int ey : y_neighbour) {

            ey= (ey+my)%my;
            
            //anisotropy_coeff = anisotropy[i];
            //e = e + anisotropy_coeff*(1 - lat[ex*my+ey]);
            e = e + (1 - lat[ex*my+ey]);
            i++;
            //e = e*(1.2-0.2*x/mx);
        }
    }
    return e*Eb;
}