#include <cstdlib>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <string>
#include <map>
#include <random>
#include <iomanip>
#include <chrono>
#include "assert.h"

using namespace std;

// declaring functions
bool check_occ(int rep_loc,
               int rand_n,
               int poly_l,
               vector<vector<int>> &poly_coord,
               vector<vector<int>> const &dirvec);

void slither(int rep_loc,
             int rand_n,
             int poly_l,
             vector<vector<int>> &poly_coord,
             vector<vector<int>> const &dirvec);

bool attempt_pivot(int pivot_loc,
                   int sym_op,
                   int poly_l,
                   vector<vector<int>> &poly_coord,
                   vector<vector<vector<int>>> &symvec);

float calc_rg(int poly_l, vector<vector<int>> &poly_coord);

float calc_rn(int poly_l, vector<vector<int>> &poly_coord);

void center_poly(int poly_l,
                 int poly_center,
                 int box_center,
                 vector<vector<int>> &poly_coord);

int main()
{
    // program options
    string algorithm = "ALGALG"; // pick one of reptation or pivot
    int n_steps = 1000000; // number of MC steps
    int equil_steps = 500000;
    int poly_l = NNNN; // degree of polymerization
    int save_freq = 100000000; // save properties every save_freq steps
    int save_conf = 100000000; // save configuration every save_conf steps
    int maxlag = 0; // max number of WRITTEN frames (every save_freq steps) to calculate the ACF

    // check algorithm choice
    bool pivot;
    bool reptation;
    if (algorithm == "reptation") {
        reptation = true;
        pivot = false;
    }
    else if (algorithm == "pivot") {
        pivot = true;
        reptation = false;
    }
    else {
        cout << "Invalid algorithm choice--choose reptation or pivot. Exiting...\n";
        return 0;
    }

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
    float rg;
    float rn;
    float rg_sum = 0.0;
    float rn_sum = 0.0;
    int boxl = poly_l + 3; // box length (square lattice)
    int poly_center = floor(poly_l/2); // index of monomer to center in box
    int box_center = floor(boxl/2); // location of box center
    int init_value = 0; // just used to initialize vectors with zeros

    // create vectors for data storage
    vector<vector<int>> poly_coord;
    poly_coord.resize(poly_l, vector<int>(2, init_value)); // matrix of polymer coordinates

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
    vector<float> rg_vec;
    vector<float> rn_vec;

    ofstream rgout("rg.txt");
    ofstream rnout("rn.txt");
    ofstream rgacfout("rg_acf.txt");
    ofstream rnacfout("rn_acf.txt");
    ofstream avrgout("avg_rg.txt");
    ofstream avrnout("avg_rn.txt");
    ofstream confout("conf.txt");

    // initialize random number generators
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> monomer(0, poly_l-2);
    uniform_int_distribution<> direction(0, 3);
    uniform_int_distribution<> head_tail(0, 1);
    uniform_int_distribution<> symmetry(0, 3); // 4 possible symmetry operations

    // initialize polymer as horizontal line in center of box
    for (int i=0; i<poly_l; i++) {
        poly_coord[i][0] = i + 1;
        poly_coord[i][1] = box_center;
    }
    center_poly(poly_l, poly_center, box_center, poly_coord);

    // calculate radius of gyration and squared end-to-end distance of initial polymer
    rg = calc_rg(poly_l, poly_coord);
    rn = calc_rn(poly_l, poly_coord);
    rg_vec.push_back(rg);
    rn_vec.push_back(rn);

    // get the start time
    auto start_time = chrono::steady_clock::now();

    // run MC simulation for n_steps
    for (int n=0; n<n_steps; n++) {
        if (reptation) {
            // pick 'head' or 'tail' of polymer to attempt reptation
            rep_loc = head_tail(gen);
            // pick random direction to attempt move in
            rand_n = direction(gen);
            // check if new attempted site is occupied after moving
            accept = check_occ(rep_loc, rand_n, poly_l, poly_coord, dirvec);
            // if accepted, move polymer along contour
            if (accept) {
                accprob++;
                slither(rep_loc, rand_n, poly_l, poly_coord, dirvec);
                // shift new polymer to center of box
                center_poly(poly_l, poly_center, box_center, poly_coord);
                rg = calc_rg(poly_l, poly_coord);
                rn = calc_rn(poly_l, poly_coord);
            }
        }
        else if (pivot) {
            pivot_loc = monomer(gen);
            sym_op = symmetry(gen);
            pivot_dir = head_tail(gen);
            accept = attempt_pivot(pivot_loc, sym_op, poly_l, poly_coord, symvec);
            if (accept) {
                accprob++;
                center_poly(poly_l, poly_center, box_center, poly_coord);
                rg = calc_rg(poly_l, poly_coord);
                rn = calc_rn(poly_l, poly_coord);
            }
        }

        if (n%save_freq == 0 && n>equil_steps) {
            rg_vec.push_back(rg);
            rn_vec.push_back(rn);
        }

        if (n%save_conf == 0 && n>equil_steps) {
            for (int j=0; j<poly_l; j++) {
                confout << n << " " << poly_coord[j][0] << " " << poly_coord[j][1] << "\n";
            }
            confout << "\n\n";
        }

        if (n%500000 == 0) {
            cout << "step " << n << "\n";
        }
    }
    cout << "Acceptance Probability: " << float(accprob)/float(n_steps) << "\n";

    // output the time to screen
    auto end_time = chrono::steady_clock::now();
    cout << "Time/step(ms) = " << chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count()/(float(n_steps)) << "\n";

    // write to output files
    
    rgout << "# step rg^2\n";
    rnout << "# step rn^2\n";

    for (auto const& x : rg_vec) {
        rgout << x << "\n";
        rg_sum += x;
    }
    for (auto const& x : rn_vec) {
        rnout << x << "\n";
        rn_sum += x;
    }

    // calculate ACF
    if (maxlag > 0) {
        double denomg = 0.0;
        double denomn = 0.0;
        double numg = 0.0;
        double numn = 0.0;
        vector<double> rkg;
        vector<double> rkn;
        int size = rg_vec.size();
        //float avg_rg = reduce(rg_vec.begin(), rg_vec.end()) / float(rg_vec.size());
        //float avg_rn = reduce(rn_vec.begin(), rn_vec.end()) / float(rn_vec.size());
        double avg_rg = rg_sum / double(rg_vec.size());
        double avg_rn = rn_sum / double(rn_vec.size());
        avrgout << avg_rg << "\n";
        avrnout << avg_rn << "\n";
        for(int i=0; i<size; i++){
            denomg += (rg_vec[i]-avg_rg)*(rg_vec[i]-avg_rg);
            denomn += (rn_vec[i]-avg_rn)*(rn_vec[i]-avg_rn);
        }
        for(int k=0; k<maxlag; k++){
            for(int i=0; i<size-k; i++){
                numg += (rg_vec[i]-avg_rg)*(rg_vec[i+k]-avg_rg);
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

    cout << "Average rg^2: " << rg_sum/rg_vec.size() << "\n";
    cout << "Average rn^2: " << rn_sum/rn_vec.size() << "\n";

    return 0;
}

// function to check if trial position would be occupied after reptation
bool check_occ(int rep_loc, int rand_n, int poly_l, vector<vector<int>> &poly_coord, vector<vector<int>> const &dirvec)
{
    int new_site;
    int x_old;
    int x_new;
    int y_old;
    int y_new;
    int opposite_x;
    int opposite_y;
    if (rep_loc == 0) {
        x_old = poly_coord[0][0];
        y_old = poly_coord[0][1];
        opposite_x = poly_coord[poly_l-1][0];
        opposite_y = poly_coord[poly_l-1][1];
        x_new = x_old + dirvec[rand_n][0];
        y_new = y_old + dirvec[rand_n][1];
    }
    else if (rep_loc == 1) {
        x_old = poly_coord[poly_l-1][0];
        y_old = poly_coord[poly_l-1][1];
        opposite_x = poly_coord[0][0];
        opposite_y = poly_coord[0][1];
        x_new = x_old + dirvec[rand_n][0];
        y_new = y_old + dirvec[rand_n][1];
    }
    // if new site is occupied, check if it is the opposite end of polymer
    // accept if so since the opposite end will move
    for (int i=0; i<poly_l; i++) {
        if (x_new == poly_coord[i][0] && y_new == poly_coord[i][1]) {
            if (x_new == opposite_x && y_new == opposite_y) {
                return true;
            }
            else {
                return false;
            }
        }
    }
    return true;
}

float calc_rg(int poly_l, vector<vector<int>> &poly_coord)
{
    int j;
    float dx = 0.0;
    float dy = 0.0;
    float x_sum_sq = 0.0;
    float y_sum_sq = 0.0;
    float x_avg = 0.0;
    float y_avg = 0.0;
    float rg = 0.0;
    for (j=0; j<poly_l; j++) {
        x_avg += poly_coord[j][0];
        y_avg += poly_coord[j][1];
    }
    x_avg = x_avg/float(poly_l);
    y_avg = y_avg/float(poly_l);
    for (j=0; j<poly_l; j++) {
        dx = poly_coord[j][0] - x_avg;
        dy = poly_coord[j][1] - y_avg;
        x_sum_sq += dx*dx;
        y_sum_sq += dy*dy;
    }
    rg = (x_sum_sq + y_sum_sq)/float(poly_l);
    return rg;
}

float calc_rn(int poly_l, vector<vector<int>> &poly_coord)
{
    float rn = 0.0;
    float dx = 0.0;
    float dy = 0.0;
    dx = poly_coord[0][0] - poly_coord[poly_l-1][0];
    dy = poly_coord[0][1] - poly_coord[poly_l-1][1];
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

bool attempt_pivot(int pivot_loc, int sym_op, int poly_l, vector<vector<int>> &poly_coord, vector<vector<vector<int>>> &symvec)
{
    int new_size = poly_l - pivot_loc;
    auto start = poly_coord.begin() + pivot_loc;
    auto end = poly_coord.end();
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
                return false;
            }
        }
    }
    // modify original coordinates if accepted
    for (int i=pivot_loc; i<poly_l; i++) {
        poly_coord[i][0] = modified_matrix[i-pivot_loc][0];
        poly_coord[i][1] = modified_matrix[i-pivot_loc][1];
    }
    return true;
}

