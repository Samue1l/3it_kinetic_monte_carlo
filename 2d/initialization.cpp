#include "initialization.h"
void init_lat(unsigned char* lat, const int mx, const int my) {
    /*
    Fill the lat array with the original morphology
    */
    double v, p;

    for(int i =0; i<mx;i++){
        for(int j=0; j<my;j++){
            v= unif(rng);

            if (i<0.8*mx) {
                p=0.5;
            } else {
                p=0;
            }

            lat[i*my+j] = v < p ? 0:1;
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

void compute_rates(unsigned char* lat, float* arr, int x, int y, const int mx, const int my, const float beta, const int n) {
    /*
    Compute all the jump rates related to one particule (ex: if each particule has 6 neighbours there's 6 different jump rates)
    The rates are calculated both ways so they're divided by 2 in order to leave the total sum inchanged

    A futur improvment would be to change the rates array to calculate only one rates per jump 
    */
    unsigned char v = lat[x*my + y];
    float rate;

    auto x_neighbour = get_neighbourhood(x);
    auto y_neighbour = get_neighbourhood(y);

    float e_init = compute_energy(lat, x, y, mx, my);
    int sub_idx = 0;
    for(int ex : x_neighbour){
        ex = (ex+mx)%mx;
        for(int ey : y_neighbour) {
            ey= (ey+my)%my;

            if(lat[ex*my+ey] != v ) {

                float e_final = compute_energy(lat, ex, ey, mx, my);
                if(v==1){
                    rate = std::min(1.0, 0.5*std::exp(beta*(e_init-e_final))); // Divided by 2 because the rates go both ways
                } else {
                    rate = std::min(1.0, 0.5*std::exp(beta*(e_final-e_init)));
                }
                int idx = n*(x*my+y) + sub_idx;
                arr[idx] = rate;

                sub_idx++;
            } else if(ex != x || ey != y) {
                rate = 0.0f;
                int idx = n*(x*my+y) + sub_idx;
                arr[idx] = rate;

                sub_idx++;
            }
        }
    }
}

float* init_rates(unsigned char* lat, float* arr, const int mx, const int my, const float beta, const int n) {
    /*
    Compute jump rates for all atoms
    */

    for(int x=0;x<mx;x++){
        for(int y=0; y<my;y++){

            compute_rates(lat, arr, x, y, mx, my, beta, n);

        }
    }
    return arr;
}
