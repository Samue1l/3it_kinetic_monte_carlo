
#include "libs.h"
#include "initialization.h"
#include "update.h"

float get_total_energy(unsigned char* lat, const int mx, const int my) {
/*Compute the total energy of the simulation box
*/
    float tot_e = 0.0f;
    for(int x=0;x<mx;x++){
        for(int y=0;y<my;y++){
            if(lat[x*my+y]==1){
                tot_e = tot_e + compute_energy(lat, x, y, mx, my)  ;
            }
        }
    }
    return tot_e;
}

int get_count(unsigned char* lat, const int mx, const int my) {
/*Compute the number of bulk sites
 */
    int count=0;
    for(int x=0;x<mx;x++){
        for(int y=0;y<my;y++){
            if(lat[x*my+y]==1){
                count++;
            }
        }
    }
    return count;
}


void write_lat(unsigned char* lat, std::string filename, const int mx, const int my) {
/*Write the disposition of the lattice to a text file
*/
printf("Writing lattice to %s...\n", filename.c_str());

  std::ofstream f;
  f.open(filename);
  if (f.is_open()) {
    f << mx << " ";
    f << my << " ";
    f << std::endl;
    for (int i = 0; i < mx; i++) {
      for (int j = 0; j < my; j++) {
            f << (int)lat[i * my + j ] << " ";

      }
    }
  }
  f.close();
}


int main(void) {
    const int mx = 128; // Dimension X
    const int my = 128; // Dimension Y
    const int n = 8; // Number of neighbours
    const float beta = 1.5; // Inverse of kT
    int niter = 2000000; // Number of iterations

    double t = 0.0;
    float dt, e;
    int step = (int) round(niter/10)-1;
    int f_step = (int) round(niter/100)-1;
    int tt = step/10;
    int nf =0;
    std::string base_file_name = "./output_files/iter_";

    unsigned char* box = new unsigned char[mx*my];
    float* transit = new float[mx*my*n];

    std::vector<double> rates, c_rates;
    std::unordered_map<float, std::vector<int>> rate_dict;
    
    rng.seed(ss);

    std::string ori = "./output_files/ori.txt";
    std::string final_f = "./output_files/lat.txt";

    printf("Initialising lattice and rates ----\n");
    init_lat(box,mx,my);
    init_rates(box,transit, mx, my, beta, n);
    calculate_prob(transit, rates, c_rates, n*my*mx);
    init_dict(rate_dict, transit, n*mx*my);

    printf("%d sites ----\n", mx*my);
    printf("Energy: %.1f -- Bulk sites: %d\n", get_total_energy(box,mx,my), get_count(box,mx,my));
    printf("%f\n", unif(rng));
    write_lat(box, ori, mx, my);
    auto t0 = std::chrono::high_resolution_clock::now();

    // Update loop
    for(int i=0;i<niter;i++) {
        dt = update_lattice(box, transit, rates, c_rates, rate_dict, mx, my, n, beta, i);
        t = t+dt;
        if(i%step==0){
            e = get_total_energy(box, mx,my);
            printf("t: %.5f / Energy: %.1f\n",t,e );
            std::string if_name = base_file_name + std::to_string(nf) + std::string(".txt");
            write_lat(box, if_name, mx, my);
            nf++;
        }

    }

    auto t1 = std::chrono::high_resolution_clock::now();
    double duration = (double) std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count();

    printf("%d iterations done in %f sec----\n", niter, duration*1e-6);
    write_lat(box, final_f, mx, my);
    delete[] box;
    delete[] transit;

    return 0;
}
