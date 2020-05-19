
#include "update.h"

void add_rate(float rate, std::vector<double> &rates, std::vector<double> &c_rates) {
/*
Add rate to the rates array and update the count array
*/
    if(rate !=0){
                std::vector<double>::iterator pt = std::find(std::begin(rates), std::end(rates), rate);
                if(pt != std::end(rates)) {
                    int idx = std::distance(rates.begin(), pt);
                    c_rates[idx]++;
                } else {
                    rates.push_back(rate);
                    c_rates.push_back(1.0f);
                }
            }
}

void decrease_rate(float rate,  std::vector<double> &rates, std::vector<double> &c_rates){
/*
Remove rate to the rates array and update the count array
*/
    if(rate != 0) {
                std::vector<double>::iterator pt = std::find(std::begin(rates), std::end(rates), rate);
                if(pt!= std::end(rates)) {
                    int idx = std::distance(rates.begin(), pt);
                    c_rates[idx]--;
                } else {
                    throw "Rate not found";
                }
            }
}
void recompute_rates(unsigned char* lat, float* arr, int x, int y, std::vector<double> &rates, std::vector<double> &c_rates, std::unordered_map<float, std::vector<int>> &dict,
                        const int mx, const int my, const float beta, const int n) {
/*
Update the transition array and mapping for an atom at position x,y,z
*/
    unsigned char v = lat[x*my + y];
    float rate, old_rate;
    std::vector<int>* trans_idx;

    std::array<int, 3> x_neighbour = get_neighbourhood(x);
    std::array<int, 3> y_neighbour = get_neighbourhood(y);

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
                old_rate = arr[idx];
                arr[idx] = rate;

                decrease_rate(old_rate, rates, c_rates);
                if(old_rate!=0){
                    trans_idx = &dict[old_rate];
                    (*trans_idx).erase(std::remove((*trans_idx).begin(), (*trans_idx).end(), idx));
                }
                add_rate(rate, rates, c_rates);
                dict[rate].push_back(idx);


                int other_idx = n*(ex*my+ey) + n-1-sub_idx;
                old_rate = arr[other_idx];
                arr[other_idx]=0.0f;
                decrease_rate(old_rate, rates, c_rates);
                if(old_rate!=0){
                    trans_idx = &dict[old_rate];
                    (*trans_idx).erase(std::remove((*trans_idx).begin(), (*trans_idx).end(), other_idx));
                }

                sub_idx = sub_idx + 1;
            } else if(ex != x || ey != y) {
                rate = 0.0f;
                int idx = n*(x*mx+y) + sub_idx;
                int other_idx = n*(ex*my+ey) + n-1-sub_idx;

                old_rate = arr[idx];
                arr[idx] = rate;
                decrease_rate(old_rate, rates, c_rates);
                if(old_rate!=0){
                    trans_idx = &dict[old_rate];
                    (*trans_idx).erase(std::remove((*trans_idx).begin(), (*trans_idx).end(), idx));
                }

                old_rate = arr[other_idx];
                arr[other_idx]=0.0f;
                decrease_rate(old_rate, rates, c_rates);
                if(old_rate!=0){
                    trans_idx = &dict[old_rate];
                    (*trans_idx).erase(std::remove((*trans_idx).begin(), (*trans_idx).end(), other_idx));
                }

                sub_idx = sub_idx + 1;
            }
        }
    }
}



float get_choice(std::vector<double> &l, std::vector<double> &weights) {
/*Choose a rate among the rate array weighted by their occurence
*/
    double r = unif(rng);
    for(int i=0;i<weights.size();i++){
        if(r<weights[i] ){
            return l[i];
        }
    }
    //throw "Choice cannot be made";
    return l[weights.size()-1];
}


std::vector<int> get_transitions(float* arr, float rate, int size) {
/*
Get all the transitions corresponding to one rate
*/
    std::vector<int> r_arr;
    for(int i=0;i<size;i++) {

            if(arr[i] == rate){
                r_arr.push_back(i);
            }

    }
    return r_arr;
}

void element_wise_product(std::vector<double> &v1,std::vector<double> &v2,std::vector<double> &out) {
/*
Compute the element wise product of 2 vectors of double
 */
    for(int i=0;i<v1.size();i++){
        out.push_back(v1[i]*v2[i]);
    }
}



double update_lattice(unsigned char* lat, float* arr,std::vector<double> &rates, std::vector<double> &c_rates,std::unordered_map<float, std::vector<int>> &dict,
                             const int mx, const int my, const int n, const float beta, int i) {
/* Update the lattice by choosing a transition and realizing it
*/
    std::vector<double> l_rates, weights;
    std::vector<int>* possible_transitions;
    int xo,yo,xf,yf;
    int chosen_transition;
    double R, chosen_rate;
    double dt;
    int size = mx*my*n;
    unsigned char lxy;

    element_wise_product(rates,c_rates,l_rates);


    R = std::accumulate(l_rates.begin(), l_rates.end(),0.0);

    weights.push_back(l_rates[0]/R);
    for(int i=1;i<l_rates.size();i++){
        weights.push_back(weights[i-1]+l_rates[i]/R);
    }


    chosen_rate = get_choice(rates, weights);
    possible_transitions = &dict[chosen_rate];
    //int r_idx = *select_randomly(possible_transitions.begin(), possible_transitions.end());
    int r_idx = round(unif(rng)*((*possible_transitions).size()-1));

    chosen_transition = (*possible_transitions)[r_idx];

    int o_i = (int) chosen_transition/n;
    int sub_i = chosen_transition % n;

    xo = (int) o_i/my;
    yo = o_i % my;

    switch (sub_i)
    {
    case 0:
         xf = (xo==0) ? mx-1 : xo-1;
         yf = (yo==0) ? my-1 : yo-1;
        break;
    case 1:
         xf = (xo==0) ? mx-1 : xo-1;
         yf = yo;
        break;
    case 2:
         xf = (xo==0) ? mx-1 : xo-1;
         yf = (yo==my-1) ? 0 : yo+1;
        break;
    case 3:
         xf = xo;
         yf = (yo==0) ? my-1 : yo-1;
        break;
    case 4:
         xf = xo;
         yf = (yo==my-1) ? 0 : yo+1;;
        break;
    case 5:
         xf = (xo==mx-1) ? 0 : xo+1;
         yf = (yo==0) ? my-1 : yo-1;
        break;
    case 6:
         xf = (xo==mx-1) ? 0 : xo+1;
         yf = yo;
        break;
    case 7:
         xf = (xo==mx-1) ? 0 : xo+1;
         yf = (yo==my-1) ? 0 : yo+1;
        break;
    }

    lxy = lat[xo*my+yo];
    lat[xo*my+yo] = lat[xf*my+yf];
    lat[xf*my+yf] = lxy;

    recompute_rates(lat,arr,xo,yo,rates,c_rates,dict,mx,my,beta,n);
    recompute_rates(lat,arr,xf,yf,rates,c_rates,dict,mx,my,beta,n);

    dt = -1.0*std::log(unif(rng))/R;

    return dt;
}
