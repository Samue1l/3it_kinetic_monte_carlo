#include "libs.h"
#include "initialization.h"
#include "update.h"


float get_total_energy(unsigned char* lat, const int mx, const int my, const int mz) {
/*Compute the total energy of the simulation box
*/
    float tot_e = 0.0f;
    for(int x=0;x<mx;x++){
        for(int y=0;y<my;y++){
            for(int z=0;z<mz;z++){
                if(lat[x*my*mz+y*mz+z]==1){
                    tot_e = tot_e + compute_energy(lat, x, y,z, mx, my, mz)  ;
                }
            }
        }
    }
    return tot_e;
}

int get_count(unsigned char* lat, const int mx, const int my, const int mz) {
/*Compute the number of bulk sites
 */
    int count=0;
    for(int x=0;x<mx;x++){
        for(int y=0;y<my;y++){
            for(int z=0; z<mz;z++) {
                if(lat[x*my*mz+y*mz+z]>=1){
                    count++;
                }
            }
        }
    }
    return count;
}

void write_lat(unsigned char* lat, std::string filename, const int mx, const int my, const int mz) {
/*Write the disposition of the lattice to a text file
*/
printf("Writing lattice to %s...\n", filename.c_str());

  std::ofstream f;
  f.open(filename);
  if (f.is_open()) {
    f << mx << " ";
    f << my << " ";
    f << mz << " ";
    f << std::endl;
    for (int i = 0; i < mx; i++) {
      for (int j = 0; j < my; j++) {
          for(int k=0; k<mz; k++){
            f << (int)lat[i * my * mz + j * mz + k] << " ";
        }
      }
    }
  }
  f.close();
}

int main(void) {
    const int mx = 32; // Dimension X
    const int my = 32; // Dimension Y
    const int mz = 48; // Dimension Z
    const int n = 26;  // Number of neighbours for annealing
    const float beta = 1.5; // Inverse of kT
    float epi_rate = std::exp(-8); // Epitaxial rate

    int niter = 200000; // Number of iterations
    double t = 0.0;
    float dt, e; 
    int step = (int) round(niter/10)-1;
    int f_step = (int) round(niter/10)-1;
    int nf =0;

    bool wf = false;
    std::string base_file_name = "./output_files/iter_";

    unsigned char* box = new unsigned char[mx*my*mz];
    int* epi_height = new int[mx*my];
    std::vector<double> rates, c_rates;
    std::unordered_map<float, std::vector<int>> rate_dict;
    srand(24101);

    std::string ori = "./output_files/ori.txt";
    std::string final_f = "./output_files/lat_128.txt";

    printf("Initialising lattice and rates ----\n");
    init_lat(box, mx,my,mz);

    float* transit = init_rates(box,mx,my,mz, beta,n);
    init_epitaxy(box, transit, epi_height, epi_rate, mx, my, mz, n);

    const int size = n*my*mx*mz + my*mx;
    calculate_prob(transit, rates, c_rates, size);
    init_dict(rate_dict, transit, size);

    printf("%d sites ----\n", mx*my*mz);
    printf("Energy: %.1f -- Bulk sites: %d\n", get_total_energy(box,mx,my,mz), get_count(box,mx,my,mz));
    write_lat(box, ori, mx, my, mz);
    auto t0 = std::chrono::high_resolution_clock::now();

    // Update loop
    for(int i=0;i<niter;i++) {
        dt = update_lattice(box, epi_height,transit, rates, c_rates, rate_dict, mx, my, mz, n, beta);
        t = t+dt;
        if(i%step==0){
            e = get_total_energy(box, mx,my,mz);
            printf("t: %.5f / Energy: %.1f\n",t,e );
            //write_lat(box, final_f, mx, my, mz);
        }

        if(i%f_step==0 && wf) {
            std::string if_name = base_file_name + std::to_string(nf) + std::string(".txt");
            write_lat(box, if_name, mx, my, mz);
            nf++;
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    double duration = (double) std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count();

    printf("%d iterations done in %f sec----\n", niter, duration*1e-6);
    write_lat(box, final_f, mx, my, mz);

    delete[] box;
    delete[] transit;
    delete[] epi_height;

    return 0;
}
