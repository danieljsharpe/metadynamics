/*
C++ code to simulate metadynamics (MTD) and obtain a free energy surface (FES)
from a potential energy surface (PES).
Daniel J. Sharpe
*/

#include <iostream>
#include <cmath>
#include <vector>
#include <array>
using namespace std;

typedef vector<pair<float,float>> fp_vec;

/* */
class MD_Simn {

    friend Walker;

    // Langevin MD step
    double langevin_md() {}
};

/* class for a single walker in multiple-walker metadynamics (MW-MTD) */
class Walker {

    private:

    fp_vec gauss_dumped; // vector of std::pair(gauss_mean,gauss_stddev) for dumped Gaussians
    vector<double> xt; // current coordinates

    public:

    double T; // simulation temperature for this walker

    Walker(double T, vector<double> xt);
};

Walker::Walker(double temp) { T = temp; };

/* Metadynamics class sets up the walkers and drives simulation */
class Metadynamics {

    private:

    bool wt; // if true, we use well-tempered MTD (rescale height of Gaussians)
    int N_replicas; // number of walkers (have MW-MTD if N_replicas > 1)
    int n_steps; // total number of MD steps
    int n_dep; // step interval for depositing
    double T_lo, T_hi; // low/high temperature distribution of walkers
    double W; // height of Gaussians
    double tau_G; // deposition stride
    int n_gauss; // total no. of Gaussians to be deposited
    vector<double> xt; // current coordinates

    public:

    Metadynamics(const int N_rep, int n_md, int n_intvl1, double T1, double T2, double gauss_h, \
                 double dep_stride, vector<double> x0);
    ~Metadynamics();
    vector<Walker> walkers;
};

Metadynamics::Metadynamics(const int N_rep, int n_md, int n_intvl1, double T1, double T2, double gauss_h, \
                           double dep_stride, vector<double> x0) {

    N_replicas = N_rep;
    n_steps = n_md;
    n_dep = n_intvl1;
    T_lo = T1;
    T_hi = T2;
    W = gauss_h;
    tau_G = dep_stride;
    n_gauss = n_steps / n_dep;
    xt = x0;

    // initialise replicas
    for (int i=0;i<N_replicas;i++) {
        walkers.emplace_back(Walker(T_lo+(double(i)*((T_hi-T_lo)/double(N_replicas-1)))));
        cout << "Walker #" << i << " temp: " << walkers[i].T << "\n";
    }
}

Metadynamics::~Metadynamics() {

    cout << "Called destructor\n";
}

/* Zwanzig's arbitrary rough potential. Parabolic with many small potential barriers
   superimposed, the amplitude parameter eps determining the roughness.
   The added linear term breaks symmetry, by an extent according to parameter zeta.
   The frequency of the superimposed cosine and sine waves is given by a & b
   parameters, respectively  */
double zwanzig_pot (const vector<double> &x, double eps, double zeta, double a, double b) {

    double u = 0.;
    for (auto x_i: x) {
        u += pow(x_i,2.) + (eps*(cos(a*x_i)+sin(b*x_i))) - (zeta*x_i); }
    return u;
}

int main () {

    // parameters for Zwanzig's potential
    double eps = 0.02;
    double zeta = 0.05;
    double a = 167.;
    double b = 73.;

    vector<double> x0 = {0.5,0.5};

    cout << "Calling Zwanzig's potential: " << zwanzig_pot(x0,eps,zeta,a,b) << "\n";

    Metadynamics mtd1(5,5000,100,1.0,2.0,1.5,0.2,x0);

    return 0;
}
