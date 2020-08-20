#include "initialization.h"

bool is_cyl(int xc, int yc, int zc, int depth, float radius, int x, int y, int z) {
    bool is_in;
    bool r_check = pow(y-yc,2)+pow(x-xc,2) <= pow(radius,2);

    if(r_check && z< z+depth && z>=zc) {
        return true;
    } else {
        return false;
    }

}

void init_lat(unsigned char* lat, unsigned char* copied_lat, const int mx, const int my, const int mz) {
    /*
    Fill the lat array with the original morphology
    */
    double v;
    unsigned char site;


    const int h_floor = 20;

    int xc1 = 10;
    int yc1 = 10;
    int R1 = 4;

    int xc2 = 22;
    int yc2 = 22;
    int R2 = 4;

    int zc = 4;
    int depth = h_floor-zc;


    for(int i =0; i<mx;i++){
        for(int j=0; j<my;j++){
            for(int k=0; k<mz; k++) {

            bool in1 = is_cyl(xc1,yc1,zc,depth, R1, i,j,k);
            bool in2 = is_cyl(xc2,yc2,zc,depth, R2, i,j,k);

            v= unif(rng);


            if(k<h_floor && !in1 && !in2) {
                site = 1;

            } else {
                site = 0;
            }

            lat[i*my*mz+ j*mz + k] = site;
            copied_lat[i*my*mz+ j*mz + k] = site;
        }
        }
    }
}


void init_dict(std::unordered_map<float, std::vector<int>> &dict, float* rates, const int size) {
    /*
    Fill the unordered map with the pattern
    key=rates (float number)
    value=vector containing the indices (int) of the transitions in the transit array (or arr in functions)
    */
    float r;
    for(int i = 0;i<size;i++) {
        r = rates[i];
        if(dict.find(r) != dict.end()) {
            dict[r].push_back(i);
        } else {
            std::vector<int> v{i};
            dict[r] = v;
        }
    }
}

void calculate_prob(float* arr, std::vector<double> &rates, std::vector<double> &count, int size) {
    /*
    Initialize the rates array and the count array containing the number of occurences of each rates
    */
    for(int i=0; i<size; i++) {

            float value = arr[i];

            if(value !=0){
                std::vector<double>::iterator pt = std::find(std::begin(rates), std::end(rates), value);
                if(pt!= std::end(rates)) {
                    int idx = std::distance(rates.begin(), pt);
                    count[idx]++;
                } else {
                    rates.push_back(value);
                    count.push_back(1.0f);
                }
            }

    }
}

void compute_rates(unsigned char* lat, float* arr, int x, int y, int z, const int mx, const int my, const int mz, const float beta, const int n) {
    /*
    Compute all the jump rates related to one particule (ex: if each particule has 6 neighbours there's 6 different jump rates)
    The rates are calculated both ways so they're divided by 2 in order to leave the total sum inchanged

    A futur improvment would be to change the rates array to calculate only one rates per jump
    */
    unsigned char v = lat[x*my*mz + y*mz + z];
    unsigned char v_reduced = v > 0 ? 1:0;
    const double prefactor = 1e12;

    float rate;

    auto x_neighbour = get_neighbourhood(x);
    auto y_neighbour = get_neighbourhood(y);
    auto z_neighbour = get_neighbourhood(z);

    float e_init = compute_energy(lat, x, y, z, mx, my, mz);
    int sub_idx = 0;

    for(int ex : x_neighbour){
        ex = (ex+mx)%mx;
        for(int ey : y_neighbour) {
            ey= (ey+my)%my;
            for(int ez :z_neighbour) {

                unsigned char atom_reduced = lat[ex*my*mz+ey*mz+ez];
                atom_reduced = atom_reduced >0 ? 1:0; // TODO: Change checks for rate value by using separate function
                if(atom_reduced != v_reduced && ez >= 0 && ez < mz) {

                    float e_final = compute_energy(lat, ex, ey, ez, mx, my, mz);
                    if(v_reduced==1){
                        rate = prefactor*std::min(1.0, 0.5*std::exp(beta*(e_init-e_final))); // Divided by 2 because the rates go both ways
                    } else {
                        rate = prefactor*std::min(1.0, 0.5*std::exp(beta*(e_final-e_init)));
                    }
                    int idx =  n*(x*mz*my+y*mz+z) + sub_idx;
                    arr[idx] = rate;

                    sub_idx++;
                } else if(ex != x || ey != y || ez != z) {
                    rate = 0.0f;
                    int idx =  n*(x*mz*my+y*mz+z) + sub_idx;
                    arr[idx] = rate;

                    sub_idx++;
                }
            }
        }
    }

}


float* init_rates(unsigned char* lat, const int mx, const int my, const int mz, const float beta, const int n) {
    /*
    Compute jump rates for all atoms
    */
    float* arr = new float[mx*my*mz*n + mx*my];

    for(int x=0;x<mx;x++){
        for(int y=0; y<my;y++){
            for(int z=0; z<mz; z++){
                compute_rates(lat, arr, x, y, z, mx, my, mz, beta, n);
            }
        }
    }
    return arr;
}


bool is_in_epi_zone(int x, int y) {
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

void recompute_epitaxy(unsigned char* lat, int* epi, int xo, int yo, const int mx, const int my, const int mz) {
/*Update the epitaxial array by scanning from top to bottom at (xo, yo)
 */
    unsigned char atom;
    bool found_surface = false;
    //    bool is_in = is_in_epi_zone(xo,yo);
    bool is_in = true;

    for(int ez=mz-1;ez>=0;ez--){
        atom = lat[xo*my*mz+mz*yo+ez];
        if( atom >=1) {
            if(found_surface || is_in == false) {
                lat[xo*my*mz+mz*yo+ez] = 1;
            } else {
                lat[xo*my*mz+mz*yo+ez] = 2;
                epi[xo*my+yo] = ez;
                found_surface=true;
            }
        }
    }

}


void init_epitaxy(unsigned char* lat, float* arr, int* epi, float epitaxial_rate, const int mx, const int my, const int mz, const int n) {
    /*
    Compute the rates for the deposition sites
    */
    unsigned char atom;
    for(int x=0;x<mx;x++){
        for(int y=0; y<my;y++){
            arr[n*mx*my*mz + x*my+ y] = 0.0f; //Fill with 0

            recompute_epitaxy(lat, epi, x, y, mx, my, mz);
            for(int z=0; z<mz; z++){

                atom = lat[x*my*mz + y*mz + z];
                if(atom == 2) {
                    epi[x*my+y] = z;
                    arr[n*mx*my*mz + x*my+ y] = epitaxial_rate;
                }
            }

        }
    }
}
