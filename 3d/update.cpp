
#include "update.h"

void add_rate(float rate, std::vector<double> &rates, std::vector<double> &c_rates) {
/*
Add rate to the rates array and update the count array
*/
    if(rate !=0){
                std::vector<double>::iterator pt = std::find(std::begin(rates), std::end(rates), rate);
                if(pt!= std::end(rates)) {
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
void recompute_rates(unsigned char* lat, float* arr, int x, int y,int z, std::vector<double> &rates, std::vector<double> &c_rates, std::unordered_map<float, std::vector<int>> &dict,
                        const int mx, const int my, const int mz, const float beta, const int n) {
/*
Update the transition array and mapping for an atom at position x,y,z
*/
    unsigned char v = lat[x*my*mz + y*mz + z];
    unsigned char v_reduced = v > 0 ? 1:0;
    float rate, old_rate;
    std::vector<int>* trans_idx;

    auto x_neighbour = get_neighbourhood(x);
    auto y_neighbour = get_neighbourhood(y);
    auto z_neighbour = get_neighbourhood(z);

    float e_init = compute_energy(lat, x, y, z, mx, my, mz);
    int sub_idx = 0;
    for(int ex : x_neighbour){
        ex = (ex+mx)%mx;
        for(int ey : y_neighbour) {
            ey= (ey+my)%my;
            for(int ez : z_neighbour) {
                //ez = (ez+mz)%mz;
                unsigned char atom = lat[ex*my*mz+ey*mz+ez];
                unsigned char atom_reduced = atom >0 ?1:0;
                if( ez >=0 && ez < mz){
                    if(atom_reduced != v_reduced) {

                        float e_final = compute_energy(lat, ex, ey, ez, mx, my, mz);
                        if(v_reduced==1){
                            rate = std::min(1.0, 0.5*std::exp(beta*(e_init-e_final))); // Divided by 2 because the rates go both ways
                        } else {
                            rate = std::min(1.0, 0.5*std::exp(beta*(e_final-e_init)));
                        }
                        int idx = n*(x*mz*my+y*mz+z) + sub_idx;
                        old_rate = arr[idx];
                        arr[idx] = rate;

                        decrease_rate(old_rate, rates, c_rates);
                        if(old_rate!=0){
                            trans_idx = &dict[old_rate];
                            (*trans_idx).erase(std::remove((*trans_idx).begin(), (*trans_idx).end(), idx));
                        }
                        add_rate(rate, rates, c_rates);
                        dict[rate].push_back(idx);


                        int other_idx = n*(ex*my*mz+ey*mz+ez) + n-1-sub_idx;
                        old_rate = arr[other_idx];
                        arr[other_idx]=0.0f;
                        decrease_rate(old_rate, rates, c_rates);
                        if(old_rate!=0){
                            trans_idx = &dict[old_rate];
                            (*trans_idx).erase(std::remove((*trans_idx).begin(), (*trans_idx).end(), other_idx));
                        }

                        //sub_idx = sub_idx + 1;
                    } else if(ex != x || ey != y || ez != z) {
                        rate = 0.0f;
                        int idx = n*(x*mz*my+y*mz+z) + sub_idx;
                        int other_idx = n*(ex*my*mz+ey*mz+ez) + n-1-sub_idx;

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




                    }
                }
                if(ex != x || ey != y || ez != z) sub_idx++;
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
    return -1;
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

void recompute_epitaxy(unsigned char* lat, int* epi, int xo, int yo, int zo, const int mx, const int my, const int mz) {
/*Update the epitaxial array by scanning the whole simulation box
 */
    unsigned char atom;
    bool found_surface = false;
    for(int ez=mz-1;ez>=0;ez--){
        atom = lat[xo*my*mz+mz*yo+ez];
        if( atom >=1) {
            if(found_surface) {
                lat[xo*my*mz+mz*yo+ez] = 1;
            } else {
                lat[xo*my*mz+mz*yo+ez] = 2;
                epi[xo*my+yo] = ez;
                found_surface=true;
            }
        }
    }

}


double update_lattice(unsigned char* lat, int* epi, float* arr,std::vector<double> &rates, std::vector<double> &c_rates,std::unordered_map<float, std::vector<int>> &dict,
                             const int mx, const int my, const int mz, const int n, const float beta) {
/* Update the lattice by choosing a transition and realizing it
*/
    std::vector<double> l_rates, weights;
    std::vector<int>* possible_transitions;
    int xo,yo,zo,xf,yf,zf;
    int chosen_transition;
    double R, chosen_rate;
    double dt;
    int size = mx*my*mz*n;
    unsigned char lxy, fxy;

    l_rates.reserve(rates.size());
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


    if(chosen_transition < n*mx*my*mz) {
        int o_i = (int) chosen_transition/n;
        int sub_i = chosen_transition % n;

        xo = (int) o_i/(my*mz);
        o_i = o_i % (my*mz);
        yo = (int) o_i/mz;
        zo = o_i % mz;

        switch(sub_i) {
    case 0:
    xf = xo - 1;
    yf = yo - 1;
    zf = zo - 1;
   break;
case 1:
    xf = xo - 1;
    yf = yo - 1;
    zf = zo ;
   break;
case 2:
    xf = xo - 1;
    yf = yo - 1;
    zf = zo + 1;
   break;
case 3:
    xf = xo - 1;
    yf = yo ;
    zf = zo - 1;
   break;
case 4:
    xf = xo - 1;
    yf = yo ;
    zf = zo ;
   break;
case 5:
    xf = xo - 1;
    yf = yo ;
    zf = zo + 1;
   break;
case 6:
    xf = xo - 1;
    yf = yo + 1;
    zf = zo - 1;
   break;
case 7:
    xf = xo - 1;
    yf = yo + 1;
    zf = zo ;
   break;
case 8:
    xf = xo - 1;
    yf = yo + 1;
    zf = zo + 1;
   break;
case 9:
    xf = xo ;
    yf = yo - 1;
    zf = zo - 1;
   break;
case 10:
    xf = xo ;
    yf = yo - 1;
    zf = zo ;
   break;
case 11:
    xf = xo ;
    yf = yo - 1;
    zf = zo + 1;
   break;
case 12:
    xf = xo ;
    yf = yo ;
    zf = zo - 1;
   break;
case 13:
    xf = xo ;
    yf = yo ;
    zf = zo + 1;
   break;
case 14:
    xf = xo ;
    yf = yo + 1;
    zf = zo - 1;
   break;
case 15:
    xf = xo ;
    yf = yo + 1;
    zf = zo ;
   break;
case 16:
    xf = xo ;
    yf = yo + 1;
    zf = zo + 1;
   break;
case 17:
    xf = xo + 1;
    yf = yo - 1;
    zf = zo - 1;
   break;
case 18:
    xf = xo + 1;
    yf = yo - 1;
    zf = zo ;
   break;
case 19:
    xf = xo + 1;
    yf = yo - 1;
    zf = zo + 1;
   break;
case 20:
    xf = xo + 1;
    yf = yo ;
    zf = zo - 1;
   break;
case 21:
    xf = xo + 1;
    yf = yo ;
    zf = zo ;
   break;
case 22:
    xf = xo + 1;
    yf = yo ;
    zf = zo + 1;
   break;
case 23:
    xf = xo + 1;
    yf = yo + 1;
    zf = zo - 1;
   break;
case 24:
    xf = xo + 1;
    yf = yo + 1;
    zf = zo ;
   break;
case 25:
    xf = xo + 1;
    yf = yo + 1;
    zf = zo + 1;
   break;

        }

        xf = (xf+mx)%mx;
        yf = (yf+my)%my;
        zf = (zf+mz)%mz;

        lxy = lat[xo*my*mz+yo*mz+zo];
        fxy = lat[xf*my*mz+yf*mz+zf];
        lat[xo*my*mz+yo*mz+zo] = fxy;
        lat[xf*my*mz+yf*mz+zf] = lxy;

        recompute_epitaxy(lat, epi, xo,yo,zo,mx,my,mz);
        recompute_epitaxy(lat, epi, xf,yf,zf,mx,my,mz);
        /*

        if(lxy == 2 ) {
            int zeo = (zo-1+mz)%mz;
            int zef = (zf-1+mz)%mz;
            lat[xo*my*mz+yo*mz+zeo] = 2;
            lat[xf*mz*my+yf*mz+zef] = 1;

            epi[xo*my+yo] = zeo;
            epi[xf*my+yf] = zef;
        } else if (fxy==2) {
            int zeo = (zo-1+mz)%mz;
            int zef = (zf-1+mz)%mz;
            lat[xo*my*mz+yo*mz+zeo] = 1;
            lat[xf*mz*my+yf*mz+zef] = 2;

            epi[xo*my+yo] = zeo;
            epi[xf*my+yf] = zef;
        }
*/
        recompute_rates(lat,arr,xo,yo,zo,rates,c_rates,dict,mx,my,mz,beta,n);
        recompute_rates(lat,arr,xf,yf,zf,rates,c_rates,dict,mx,my,mz,beta,n);

    } else {
        int xo = (int) (chosen_transition-n*mx*my*mz)/my;
        int yo = chosen_transition % my;
        int zo = epi[xo*my + yo];

        int xf = xo;
        int yf = yo;
        int zf = zo + 1;
        if (zf == mz) {
            throw "Top of simulation box overflow";
        }
        lat[xo*my*mz+yo*mz+zo] = 1;
        lat[xf*my*mz+yf*mz+zf] = 2;
        epi[xf*my + yf] = zf;

        recompute_rates(lat,arr,xo,yo,zo,rates,c_rates,dict,mx,my,mz,beta,n);
        recompute_rates(lat,arr,xf,yf,zf,rates,c_rates,dict,mx,my,mz,beta,n);


    }

    dt = -1.0*std::log(time_dist(rng)+1)/R;
    return dt;
}
