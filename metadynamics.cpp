/*
C++ code to simulate metadynamics (MTD) and obtain a free energy surface (FES)
from a potential energy surface (PES).
Daniel J. Sharpe
*/

#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <random>
using namespace std;

#define OUTFILE "walker_info"; // root name for walker info output files
typedef vector<pair<vector<double>,vector<double>>> fp_vec;

/* Zwanzig's arbitrary rough potential. Parabolic with many small potential barriers
   superimposed, the amplitude parameter eps determining the roughness.
   The added linear term breaks symmetry, by an extent according to parameter zeta.
   The frequency of the superimposed cosine and sine waves is given by a & b
   parameters, respectively  */
class Pot_funcs {

    private:

    double eps; double zeta; double a; double b;

    public:

    Pot_funcs(vector<double> func_args, int func_id) {
        if (func_id==1) { // initialise Zwanzig potential
            eps = func_args[0];
            zeta = func_args[1];
            a = func_args[2];
            b = func_args[3]; }
    }
    ~Pot_funcs() {}

    double zwanzig_pot (const vector<double> &x) {

        double u = 0.;
        for (auto x_i: x) {
            u += pow(x_i,2.) + (eps*(cos(a*x_i)+sin(b*x_i))) - (zeta*x_i); }
        return u;
    }
};
typedef double (Pot_funcs::*Pes_membfunc)(const vector<double> &x);

class Metadynamics;

/* class containing functions (e.g. integrators) required for MD simulation */
class MD_Simn {

    friend Metadynamics;

    // return a single random number from a normal distribution with zero mean and given std_dev
    static double rand_normal(double std_dev, int seed=17) {
        static std::default_random_engine generator (seed);
        static std::normal_distribution<double> distribution(0.0,std_dev);
        return distribution(generator);
    }

    // numerical derivative by central difference method
    static vector<double> cen_diff(vector<double> xt, Pot_funcs *pot_funcs_obj, \
                                   Pes_membfunc pes_func, double h=0.001) {
        int i=0;
        vector<double> cen_diff_deriv;
        for (auto x_i: xt) {
            vector<double> xt_fwd = xt; vector<double> xt_bwd = xt;
            xt_fwd[i] += h, xt_bwd[i] -= h;
            double fwd_potval = (pot_funcs_obj->*pes_func)(xt_fwd);
            double bwd_potval = (pot_funcs_obj->*pes_func)(xt_bwd);
            cen_diff_deriv.emplace_back(-(fwd_potval-bwd_potval)/(2.*h)); // note neg. of deriv
            i++;
        }
        return cen_diff_deriv;
    }

    /* (first-order) Langevin MD single step (reduced units). Full expression for integrator is:
       r(t+dt) = r(t) + (((F(r(t))/m)*gamma)*dt) + ((((2*kB*T)/(m*gamma))**(1/2))*xi_w*(dt**(1/2)))
       gamma is a frictional coeff, xi_w a random number (representing Weiner process), dt is timestep */
    static vector<double> langevin_md(vector<double> xt, double T, Pot_funcs *pot_funcs_obj, \
                          Pes_membfunc pes_func, double gamma=0.1, double dt=0.001) {
        int i=0;
        vector<double> force_pes = cen_diff(xt,pot_funcs_obj,pes_func);
        vector<double> force_bias(2,0.1); // placeholder
        vector<double> force_tot(force_pes.size(),0.);
        for (int j=0;j<force_tot.size();j++) {
            force_tot[j] = force_pes[j] + force_bias[j]; }
        for (auto x_i: xt) {
            double xi_w = rand_normal(0.1);
            xt[i] += ((force_tot[i]/gamma)*dt) + (pow(2.*T/gamma,0.5)*xi_w*pow(dt,0.5));
            i++;
        }
        return xt;
    }
};


/* class for a single walker in multiple-walker metadynamics (MW-MTD) */
class Walker {

    private:

    int id;
    // the bias potential for a given walker is defined by a set of Gaussians (and also 
    fp_vec gauss_dumped; // vector of std::pair(gauss_mean,gauss_stddev) for dumped Gaussians
    vector<double> gauss_heights; // vector of heights for dumped Gaussians

    public:

    double T; // simulation temperature for this walker
    vector<double> xt; // current coordinates

    Walker(int walker_id, double temp, vector<double> coords);

    void walker_update_gaussians(vector<double> gauss_width, double gauss_height) {
        gauss_dumped.emplace_back(make_pair(xt,gauss_width));
        gauss_heights.emplace_back(gauss_height);
    }

    // print information (temp, bias potential info...) of walker to file
    void walker_print_info() {};
};

Walker::Walker(int walker_id, double temp, vector<double> coords) {
    id = walker_id;
    T = temp;
    xt = coords;
};


/* Metadynamics class sets up the walkers and drives simulation */
class Metadynamics {

    private:

    bool wt; // if true, we use well-tempered MTD (rescale height of Gaussians)
    int N_replicas; // number of walkers (have MW-MTD if N_replicas > 1)
    int n_steps; // total number of MD steps
    int n_dep; // step interval for depositing Gaussians
    double T_lo, T_hi; // low/high temperature distribution of walkers
    double W; // (starting, if using well-tempered MTD) height of Gaussians
    vector<double> sigma; // standard deviation of Gaussians
    vector<double> W_walkers; // if wt; current (rescaled) height of Gaussian for each walker
    double tau_G; // deposition stride
    double omega; // energy rate
    int n_gauss; // total no. of Gaussians to be deposited
    vector<double> xt; // current coordinates
    Pot_funcs *pes_inst; // instance of class containing potentials
    Pes_membfunc pes; // potential energy surface function

    public:

    Metadynamics(const int N_rep, int n_md, int n_intvl1, double T1, double T2, double gauss_h, \
                 vector<double> gauss_w, double dep_stride, vector<double> x0, Pot_funcs *pot_funcs_obj, Pes_membfunc pes_func);
    ~Metadynamics();
    vector<Walker> walkers;

    void drive_simn() {

        double pot_val; double gauss_height;

        for (int i=0;i<n_steps;i++) {
            cout << "Taking a step..." << endl;
            vector<Walker> walkers_temp = walkers; // temporary vector
            int j=0;
            for (auto mtd_walker: walkers) {
                // single step for a single walker
                //cout << "  walker #:" << j+1 << endl;
                //cout << "   coords before: " << mtd_walker.xt[0] << "  " << mtd_walker.xt[1] << endl;
                walkers_temp[j].xt = MD_Simn::langevin_md(mtd_walker.xt,mtd_walker.T,pes_inst,pes);
                //cout << "   coords after:  " << walkers_temp[j].xt[0] << "  " << walkers_temp[j].xt[1] << endl;
                pot_val = (pes_inst->*pes)(mtd_walker.xt);
                //cout << "   pot value: " << pot_val << endl;
                if ((i!=0)&&(i%n_dep==0)) { // deposit a Gaussian
                    if (wt) { // rescale Gaussian height
                    // W_walkers[j] = ...;
                    gauss_height = W_walkers[j];
                    } else { gauss_height = W; }
                    walkers_temp[j].walker_update_gaussians(sigma,gauss_height);
                }
                j++;
            }
            walkers = walkers_temp;
        }
        for (auto mtd_walker: walkers) {
            mtd_walker.walker_print_info(); }
    }
};

Metadynamics::Metadynamics(const int N_rep, int n_md, int n_intvl1, double T1, double T2, double gauss_h, \
                           vector<double> gauss_w, double dep_stride, vector<double> x0, Pot_funcs *pot_funcs_obj, Pes_membfunc pes_func) {

    N_replicas = N_rep; n_steps = n_md; n_dep = n_intvl1;
    T_lo = T1; T_hi = T2;
    W = gauss_h; sigma = gauss_w, tau_G = dep_stride; omega = W/tau_G; n_gauss = n_steps / n_dep;
    xt = x0;
    wt = false;
    pes_inst = pot_funcs_obj; pes = pes_func;
    if (wt) { vector<double> W_walkers(W,N_replicas); }

    // initialise replicas
    for (int i=0;i<N_replicas;i++) {
        walkers.emplace_back(Walker(i+1,T_lo+(double(i)*((T_hi-T_lo)/double(N_replicas-1))),xt));
        cout << "Walker #" << i << " temp: " << walkers[i].T << "\n";
    }
}

Metadynamics::~Metadynamics() {

    cout << "Called destructor\n";
}


int main () {

    // parameters for Zwanzig's potential
    double eps = 0.02; double zeta = 0.05; double a = 167.; double b = 73.;
    vector<double> zwanzig_params = {eps,zeta,a,b};

    vector<double> x0 = {0.5,0.5};

    Pot_funcs zwanzig_inst(zwanzig_params,1);
    Pes_membfunc pes = &Pot_funcs::zwanzig_pot;
    Pot_funcs *pes_inst_ptr;
    pes_inst_ptr = &zwanzig_inst;

    cout << "Calling Zwanzig's potential: " << zwanzig_inst.zwanzig_pot(x0) << "\n";
    cout << "Calling Zwanzig's pot. by member ptr: " << (zwanzig_inst.*pes)(x0) << "\n";
    cout << "And again, with the class instance as a ptr: " << (pes_inst_ptr->*pes)(x0) << "\n";

    vector<double> gauss_width(0.05,2);
    Metadynamics mtd1(5,1000,10,1.0,2.0,1.5,gauss_width,0.2,x0,pes_inst_ptr,pes);
    mtd1.drive_simn();

    return 0;
}
