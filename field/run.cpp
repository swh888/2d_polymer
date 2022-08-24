#include <cstdlib>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <string>
#include <map>
#include <random>
#include <iomanip>
#include <algorithm>
#include "assert.h"

using namespace std;

// declaring functions
bool check_occ(int rep_loc,
               int rand_n,
               int poly_l,
               int boxl,
               vector<vector<int>> &poly_coord,
               vector<vector<int>> const &dirvec);

void slither(int rep_loc,
             int rand_n,
             int poly_l,
             vector<vector<int>> &poly_coord,
             vector<vector<int>> const &dirvec);

void calc_energies(int poly_l,
                   int boxl,
                   double amplitude,
                   vector<vector<int>> &poly_coord,
                   vector<vector<int>> &e_lattice,
                   vector<double> &energies);

bool attempt_pivot(int pivot_loc,
                   int sym_op,
                   int pivot_dir,
                   float n_acc,
                   int poly_l,
                   int boxl,
                   double beta,
                   double amplitude,
                   vector<vector<int>> &poly_coord,
                   vector<vector<vector<int>>> &symvec,
                   vector<vector<int>> &e_lattice,
                   vector<double> &energies);

void calc_rg(int poly_l, int poly_center, int box_center, vector<vector<int>> &poly_coord, vector<vector<float>> &rg_vec, bool accept);

float calc_rn(int poly_l, int poly_center, int box_center, vector<vector<int>> &poly_coord);

void center_poly(int poly_l,
                 int poly_center,
                 int box_center,
                 vector<vector<int>> &poly_coord);

int main()
{
    // program options
    bool reptation = false; // pick one of reptation or pivot
    bool pivot = true;
    int n_steps = 7000000; // total number of MC steps
    int equil_steps = 2000000; // steps to discard
    int poly_l = 1000; // degree of polymerization
    int save_freq = 500; // save properties every save_freq steps
    int save_conf = 200000; // save configuration every save_conf steps
    int maxlag = 0;
    bool calc_acf = false;

    double amplitude = 0.05; // lamellar field amplitude
    int width = 20; // lamellar field width
    double beta = 1.0; // inverse temperature

    // basic checks for reasonable parameters
    assert(reptation != pivot); // check for invalid algorithm selection
    assert(n_steps > 0);
    assert(poly_l > 1);
    assert(save_freq > 0);
    assert(n_steps > equil_steps);

    // initialize variables
    int accprob = 0;
    bool accept;
    int rand_n;
    int sym_op;
    int rep_loc;
    int pivot_loc;
    int pivot_dir;
    double n_acc;
    float rg;
    float rn;
    float rg_sum = 0.0;
    float rn_sum = 0.0;
    int boxl = 2 * width * ceil((poly_l + 3)/float(2 * width)); // box length (square lattice)
    cout << "box: " << boxl << "\n";
    int poly_center = floor(poly_l/2); // index of monomer to center in box
    int box_center = floor(boxl/2); // location of box center
    int init_value = 0; // just used to initialize vectors with zeros

    // create vectors for data storage
    vector<vector<int>> poly_coord;
    poly_coord.resize(poly_l, vector<int>(3, init_value)); // matrix of polymer coordinates
    for (int i=0; i<poly_l/2; i++) {
        poly_coord[i][2] = 0;
    }
    for (int i=poly_l/2; i<poly_l; i++) {
        poly_coord[i][2] = 1;
    }

    vector<vector<int>> e_lattice;
    e_lattice.resize(boxl, vector<int>(boxl, init_value)); // lattice coordinates and A and B interactions
    for (int x=0; x<boxl; x++) {
        for (int y=0; y<boxl; y++) {
            e_lattice[x][y] = (x/width)%2;
        }
    }

    vector<double> energies;
    energies.resize(poly_l, 0.0);

    /*
    for (int x=0; x<boxl; x++) {
        for (int y=0; y<boxl; y++) {
            e_lattice[x][y][0] = (x/width)%2;
            if (e_lattice[x][y][0] == 0) {
                e_lattice[x][y][1] = 1;
            }
        }
    }
    */

    vector<vector<int>> dirvec;
    dirvec.resize(4, vector<int>(2, init_value)); // matrix of possible trial directions (reptation)
    dirvec[0][0] = -1;
    dirvec[0][1] = 0;
    dirvec[1][0] = 1;
    dirvec[1][1] = 0;
    dirvec[2][0] = 0;
    dirvec[2][1] = -1;
    dirvec[3][0] = 0;
    dirvec[3][1] = 1;

    vector<vector<vector<int>>> symvec;
    symvec.resize(4, vector<vector<int>>(2, vector<int>(2))); // symmetry function matrices
    // counter-clockwise 90 deg
    symvec[0][0][0] = 0;
    symvec[0][0][1] = -1;
    symvec[0][1][0] = 1;
    symvec[0][1][1] = 0;
    // clockwise 90 deg
    symvec[1][0][0] = 0;
    symvec[1][0][1] = 1;
    symvec[1][1][0] = -1;
    symvec[1][1][1] = 0;
    // reflect across x-axis
    symvec[2][0][0] = 1;
    symvec[2][0][1] = 0;
    symvec[2][1][0] = 0;
    symvec[2][1][1] = -1;
    // reflect across y-axis
    symvec[3][0][0] = -1;
    symvec[3][0][1] = 0;
    symvec[3][1][0] = 0;
    symvec[3][1][1] = 1;

    // output data storage
    vector<vector<float>> rg_vec;
    vector<float> rn_vec;

    ofstream rgout("rg.txt");
    ofstream rnout("rn.txt");
    ofstream confout("conf.txt");

    // initialize random number generators
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> accnumber(0.0, 1.0);
    uniform_int_distribution<> monomer(1, poly_l-2);
    uniform_int_distribution<> direction(0, 3);
    uniform_int_distribution<> head_tail(0, 1);
    uniform_int_distribution<> symmetry(0, 3); // 4 possible symmetry operations

    // initialize polymer as horizontal line in center of box
    for (int i=0; i<poly_l; i++) {
        poly_coord[i][0] = i + 1;
        poly_coord[i][1] = box_center;
    }
    center_poly(poly_l, poly_center, box_center, poly_coord);
    calc_energies(poly_l, boxl, amplitude, poly_coord, e_lattice, energies);

    // calculate radius of gyration and squared end-to-end distance of initial polymer
    accept = true;
    calc_rg(poly_l, poly_center, box_center, poly_coord, rg_vec, accept);
    rn = calc_rn(poly_l, poly_center, box_center, poly_coord);
    rn_vec.push_back(rn);

    // run MC simulation for n_steps
    for (int n=0; n<n_steps; n++) {
        if (reptation) {
            cout << "Pivot implementation incomplete... do not use. Exiting...\n";
            return 1;
            // pick 'head' or 'tail' of polymer to attempt reptation
            rep_loc = head_tail(gen);
            // pick random direction to attempt move in
            rand_n = direction(gen);
            // check if new attempted site is occupied after moving
            n_acc = accnumber(gen);
            accept = check_occ(rep_loc, rand_n, poly_l, boxl, poly_coord, dirvec);
            // if accepted, move polymer along contour
            /*
            if (accept) {
                accprob++;
                slither(rep_loc, rand_n, poly_l, poly_coord, dirvec);
                // shift new polymer to center of box
                //center_poly(poly_l, poly_center, box_center, poly_coord);
                rg = calc_rg(poly_l, poly_center, box_center, poly_coord);
                rn = calc_rn(poly_l, poly_center, box_center, poly_coord);
            }
            */
        }
        else if (pivot) {
            pivot_loc = monomer(gen);
            sym_op = symmetry(gen);
            pivot_dir = head_tail(gen);
            n_acc = accnumber(gen);
            accept = attempt_pivot(pivot_loc, sym_op, pivot_dir, n_acc, poly_l, boxl, beta, amplitude, poly_coord, symvec, e_lattice, energies);
            if (accept) {
                accprob++;
                //center_poly(poly_l, poly_center, box_center, poly_coord);
                if (n%save_freq == 0 && n>equil_steps) {
                    calc_rg(poly_l, poly_center, box_center, poly_coord, rg_vec, accept);
                }
                rn = calc_rn(poly_l, poly_center, box_center, poly_coord);
            }
            //cout << "acc: " << accept << "\n";
        }

        if (n%save_freq == 0 && n>equil_steps) {
            //rg_vec.push_back(rg);
            if (!accept) {
                calc_rg(poly_l, poly_center, box_center, poly_coord, rg_vec, accept);
            }
            rn_vec.push_back(rn);
        }

        if (n%save_conf == 0 && n>equil_steps) {
            for (int j=0; j<poly_l; j++) {
                confout << n << " " << poly_coord[j][0] << " " << poly_coord[j][1] << " " << poly_coord[j][2] << "\n";
            }
            confout << "\n\n";
        }

        if (n%500000 == 0) {
            cout << "step " << n << "\n";
        }
    }
    cout << "Acceptance Probability: " << float(accprob)/float(n_steps) << "\n";

    // write to output files
    rgout << "# step rg^2\n";
    rnout << "# step rn^2\n";

    for (int i=0; i<rg_vec.size(); i++) {
        rgout << i*save_freq << " " << rg_vec[i][0] << " " << rg_vec[i][1] << " " << rg_vec[i][2] << "\n";
    }

    for (auto const& x : rg_vec) {
        //rgout << x << "\n";
        rg_sum += x[0];
    }
    for (auto const& x : rn_vec) {
        rnout << x << "\n";
        rn_sum += x;
    }

    cout << "Average rg^2: " << rg_sum/rg_vec.size() << "\n";
    cout << "Average rn^2: " << rn_sum/rn_vec.size() << "\n";

    // calculate ACF
    if (calc_acf) {
        ofstream rgacfout("rg_acf.txt");
        ofstream rnacfout("rn_acf.txt");
        double denomg = 0.0;
        double denomn = 0.0;
        double numg = 0.0;
        double numn = 0.0;
        vector<double> rkg;
        vector<double> rkn;
        int size = rg_vec.size();
        double avg_rg = rg_sum / double(rg_vec.size());
        double avg_rn = rn_sum / double(rn_vec.size());
        for(int i=0; i<size; i++){
            denomg += (rg_vec[i][0]-avg_rg)*(rg_vec[i][0]-avg_rg);
            denomn += (rn_vec[i]-avg_rn)*(rn_vec[i]-avg_rn);
        }
        for(int k=0; k<maxlag; k++){
            for(int i=0; i<size-k; i++){
                numg += (rg_vec[i][0]-avg_rg)*(rg_vec[i+k][0]-avg_rg);
                numn += (rn_vec[i]-avg_rn)*(rn_vec[i+k]-avg_rn); 
            }
            rkg.push_back(numg/denomg);
            rkn.push_back(numn/denomn);
            numg = 0.0;
            numn = 0.0;
        }

        // write ACF to file
        for (int i=0; i<rkg.size(); i++) {
            rgacfout << i*save_freq << " " << rkg[i] << "\n";
        }
        for (int i=0; i<rkn.size(); i++) {
            rnacfout << i*save_freq << " " << rkn[i] << "\n";
        }
    }

    return 0;
}

// function to check if trial position would be occupied after reptation
bool check_occ(int rep_loc, int rand_n, int poly_l, int boxl, vector<vector<int>> &poly_coord, vector<vector<int>> const &dirvec)
{
    int new_site;
    int x_old;
    int x_new;
    int y_old;
    int y_new;
    int opposite_x;
    int opposite_y;
    int midpoint, mid_x, mid_y;
    int mon_type;
    int opposite_mon_type;
    bool acc;
    double n_acc;
    double old_energy;
    double new_energy;
    if (rep_loc == 0) {
        midpoint = poly_l/2 + 1;
        x_old = poly_coord[0][0];
        y_old = poly_coord[0][1];
        opposite_x = poly_coord[poly_l-1][0];
        opposite_y = poly_coord[poly_l-1][1];
        mid_x = poly_coord[midpoint][0];
        mid_y = poly_coord[midpoint][1];
        x_new = x_old + dirvec[rand_n][0];
        y_new = y_old + dirvec[rand_n][1];
        mon_type = poly_coord[0][2];
        opposite_mon_type = poly_coord[poly_l-1][2];
    }
    else if (rep_loc == 1) {
        midpoint = poly_l/2 + 1;
        x_old = poly_coord[poly_l-1][0];
        y_old = poly_coord[poly_l-1][1];
        opposite_x = poly_coord[0][0];
        opposite_y = poly_coord[0][1];
        mid_x = poly_coord[midpoint][0];
        mid_y = poly_coord[midpoint][1];
        x_new = x_old + dirvec[rand_n][0];
        y_new = y_old + dirvec[rand_n][1];
        mon_type = poly_coord[poly_l-1][2];
        opposite_mon_type = poly_coord[0][2];
    }
    // if new site is occupied, check if it is the opposite end of polymer
    // accept if so since the opposite end will move
    for (int i=0; i<poly_l; i++) {
        if (x_new == poly_coord[i][0] && y_new == poly_coord[i][1]) {
            if (x_new == opposite_x && y_new == opposite_y) {
                break;
            }
            else {
                return false;
            }
        }
    }
    // calculate change in energy
    return true;
}

void calc_rg(int poly_l, int poly_center, int box_center, vector<vector<int>> &poly_coord, vector<vector<float>> &rg_vec, bool accept)
{
    vector<float> rg_frame(3);
    int x_shift;
    int y_shift;
    float dx = 0.0;
    float dy = 0.0;
    float x_sum_sq = 0.0;
    float y_sum_sq = 0.0;
    float x_avg = 0.0;
    float y_avg = 0.0;
    float rg = 0.0;
    if (accept) {
        vector<vector<int>> copy_coord(poly_coord);
        x_shift = box_center - copy_coord[poly_center][0];
        y_shift = box_center - copy_coord[poly_center][1];
        for (int j=0; j<poly_l; j++) {
            copy_coord[j][0] += x_shift;
            copy_coord[j][1] += y_shift;
        }
        for (int j=0; j<poly_l; j++) {
            x_avg += copy_coord[j][0];
            y_avg += copy_coord[j][1];
        }
        x_avg = x_avg/float(poly_l);
        y_avg = y_avg/float(poly_l);
        for (int j=0; j<poly_l; j++) {
            dx = copy_coord[j][0] - x_avg;
            dy = copy_coord[j][1] - y_avg;
            x_sum_sq += dx*dx;
            y_sum_sq += dy*dy;
        }
        rg_frame[0] = (x_sum_sq + y_sum_sq)/float(poly_l);
        rg_frame[1] = x_sum_sq/float(poly_l);
        rg_frame[2] = y_sum_sq/float(poly_l);
    }
    else {
        rg_frame = rg_vec[rg_vec.size()-1];
    }
    rg_vec.push_back(rg_frame);
}

float calc_rn(int poly_l, int poly_center, int box_center, vector<vector<int>> &poly_coord)
{
    int x_shift;
    int y_shift;
    vector<vector<int>> copy_coord(poly_coord);
    x_shift = box_center - copy_coord[poly_center][0];
    y_shift = box_center - copy_coord[poly_center][1];
    for (int j=0; j<poly_l; j++) {
        copy_coord[j][0] += x_shift;
        copy_coord[j][1] += y_shift;
    }

    float rn = 0.0;
    float dx = 0.0;
    float dy = 0.0;
    dx = copy_coord[0][0] - copy_coord[poly_l-1][0];
    dy = copy_coord[0][1] - copy_coord[poly_l-1][1];
    rn = dx*dx + dy*dy;
    return rn;
}

void slither(int rep_loc, int rand_n, int poly_l, vector<vector<int>> &poly_coord, vector<vector<int>> const &dirvec)
{
    int j;
    int k;
    int x_old;
    int x_new;
    int x_shift;
    int y_old;
    int y_new;
    int y_shift;
    int opposite_x;
    int opposite_y;
    if (rep_loc == 0) {
        x_old = poly_coord[0][0];
        y_old = poly_coord[0][1];
        opposite_x = poly_coord[poly_l-1][0];
        opposite_y = poly_coord[poly_l-1][1];
        x_new = x_old + dirvec[rand_n][0];
        y_new = y_old + dirvec[rand_n][1];
        poly_coord[0][0] = x_new;
        poly_coord[0][1] = y_new;
        for (j=1; j<poly_l; j++) {
            x_shift = x_old - poly_coord[j][0];
            y_shift = y_old - poly_coord[j][1];
            x_old = poly_coord[j][0];
            y_old = poly_coord[j][1];
            poly_coord[j][0] = x_old + x_shift;
            poly_coord[j][1] = y_old + y_shift;
        }
    }
    else if (rep_loc == 1) {
        x_old = poly_coord[poly_l-1][0];
        y_old = poly_coord[poly_l-1][1];
        opposite_x = poly_coord[0][0];
        opposite_y = poly_coord[0][1];
        x_new = x_old + dirvec[rand_n][0];
        y_new = y_old + dirvec[rand_n][1];
        poly_coord[poly_l-1][0] = x_new;
        poly_coord[poly_l-1][1] = y_new;
        for (j=poly_l-2; j>=0; j--) {
            x_shift = x_old - poly_coord[j][0];
            y_shift = y_old - poly_coord[j][1];
            x_old = poly_coord[j][0];
            y_old = poly_coord[j][1];
            poly_coord[j][0] = x_old + x_shift;
            poly_coord[j][1] = y_old + y_shift;
        }
    }
}

void center_poly(int poly_l, int poly_center, int box_center, vector<vector<int>> &poly_coord)
{
    int x_shift;
    int y_shift;
    x_shift = box_center - poly_coord[poly_center][0];
    y_shift = box_center - poly_coord[poly_center][1];
    for (int j=0; j<poly_l; j++) {
        poly_coord[j][0] += x_shift;
        poly_coord[j][1] += y_shift;
    }
}

void calc_energies(int poly_l, int boxl, double amplitude, vector<vector<int>> &poly_coord, vector<vector<int>> &e_lattice, vector<double> &energies)
{
    int e_unscaled = 0;
    double energy = 0.0;
    int x, y, mon_type;
    for (int i=0; i<poly_l; i++) {
        x = poly_coord[i][0];
        y = poly_coord[i][1];
        x = x - boxl * round(x/boxl);
        y = y - boxl * round(y/boxl);
        mon_type = poly_coord[i][2];
        if (poly_coord[i][2] != e_lattice[x][y]) {
            e_unscaled = 1;
            energies[i] = amplitude * e_unscaled;
        }
        else {
            energies[i] = 0.0;
        }
    }
}

bool attempt_pivot(int pivot_loc, int sym_op, int pivot_dir, float n_acc, int poly_l, int boxl, double beta, double amplitude, vector<vector<int>> &poly_coord, vector<vector<vector<int>>> &symvec, vector<vector<int>> &e_lattice, vector<double> &energies)
{
    // reverse order if randomly chosen
    if (pivot_dir == 1) {
        pivot_loc = poly_l - pivot_loc - 1;
        reverse(poly_coord.begin(), poly_coord.end());
        reverse(energies.begin(), energies.end());
    }

    int new_size = poly_l - pivot_loc;
    int x;
    int y;
    int mon_type;
    int e_unscaled;
    bool acc;
    double acc_prob;
    double old_energy = 0.0;
    double new_energy = 0.0;
    double delta_e;
    auto start = poly_coord.begin() + pivot_loc;
    auto end = poly_coord.end();
    vector<double> new_energies(new_size);
    vector<vector<int>> coord_copy(end - start);
    vector<vector<int>> modified_matrix;
    modified_matrix.resize(new_size, vector<int>(new_size, 0));
    copy(start, end, coord_copy.begin());
    // shift pivot point to (0,0)
    for (int i=pivot_loc; i<poly_l; i++) {
        coord_copy[i-pivot_loc][0] = poly_coord[i][0] - poly_coord[pivot_loc][0];
        coord_copy[i-pivot_loc][1] = poly_coord[i][1] - poly_coord[pivot_loc][1];
    }
    // multiply by symmetry matrix to transform
    for (int i=0; i<new_size; i++) {
        for (int j=0; j<2; j++) {
            for (int k=0; k<2; k++) {
                modified_matrix[i][j] += symvec[sym_op][j][k] * coord_copy[i][k];
            }
        }
    }
    // shift pivoted section back to original coordinates
    for (int i=pivot_loc; i<poly_l; i++) {
        modified_matrix[i-pivot_loc][0] += poly_coord[pivot_loc][0];
        modified_matrix[i-pivot_loc][1] += poly_coord[pivot_loc][1];
        // check for overlaps
        for (int j=0; j<pivot_loc; j++) {
            if (modified_matrix[i-pivot_loc][0]==poly_coord[j][0] && modified_matrix[i-pivot_loc][1]==poly_coord[j][1]) {
                if (pivot_dir == 1) {
                    reverse(poly_coord.begin(), poly_coord.end());
                    reverse(energies.begin(), energies.end());
                }
                return false;
            }
        }
    }
    // calculate energies of moved points
    for (int i=pivot_loc; i<poly_l; i++) {
        //cout << "i: " << i << "\n";
        x = modified_matrix[i-pivot_loc][0];
        y = modified_matrix[i-pivot_loc][1];
        //cout << x << "\n";
        //cout << floor(x/float(boxl)) << "\n";
        x = x - boxl * floor(x/float(boxl));
        y = y - boxl * floor(y/float(boxl));
        //y = y - boxl * round(y/boxl);
        //cout << x << " " << y << "\n";
        //cout << "before\n";
        //cout << e_lattice[x][y] << "\n";
        //cout << "after\n";
        mon_type = poly_coord[i][2];
        if (mon_type == e_lattice[x][y]) {
            e_unscaled = 0;
        }
        else {
            e_unscaled = 1;
        }
        new_energies[i-pivot_loc] = amplitude * e_unscaled;
        new_energy += amplitude * e_unscaled;
        old_energy += energies[i];
    }
    // calculate acceptance prob
    delta_e = beta * (old_energy - new_energy);
    acc_prob = exp(delta_e);
    if (n_acc >= acc_prob) {
        acc = false;
    }
    else if (n_acc < acc_prob) {
        acc = true;
        // modify original coordinates if accepted
        for (int i=pivot_loc; i<poly_l; i++) {
            poly_coord[i][0] = modified_matrix[i-pivot_loc][0];
            poly_coord[i][1] = modified_matrix[i-pivot_loc][1];
            energies[i] = new_energies[i-pivot_loc];
        }
    }
    // restore to original order if reversed
    if (pivot_dir == 1) {
        reverse(poly_coord.begin(), poly_coord.end());
        reverse(energies.begin(), energies.end());
    }
    return acc;
}

