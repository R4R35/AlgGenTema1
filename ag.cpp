#include <iostream>
#include <cmath>
#include <cstring>
#include <vector>
#include <limits>
#include <random>
#include <bitset>
#include <functional>
#include <chrono>

using namespace std;
// =========================== CONSTANTS =======================================


const float PI = 2*asin(1.0);
const float EPSILON = 1e-8;

// ============================ HELPER FUNCTIONS =================================

vector<bool> random_bits(int D,int l){
    int size = D*l;
    vector<bool> r(size,false);
    mt19937 rng;
    rng.seed(time(0));
    for(int i=0;i<size;i++){
      r[i] = rng()%2;
    }
    return r;
}
/*
Nu cred ca sunt necesare implementarile la functiile astea still

vector<bool> cstringToVect(char* s){
    vector<bool>v;
    int l = strlen(s);
    for(int  ii =0 ;ii< l;++ii){
        v.push_back(s[ii] - '0');
    }
    return v;
}


 */
float decode(vector<bool> v,float low,float high){
    unsigned long long int sum = 0;
    unsigned long long int mul = 1;
    for(int ii =v.size()-1;ii>=0;--ii){
        sum += v[ii] * mul;
        mul <<= 1;
    }
    float div = pow( 2,v.size()-1);
    float x = float(sum) / div;
    x= low + x*(high-low);
    return x;
}

vector<double> decode_Bits(vector<bool> bitstring,int D,int l, float low,float high) {
    vector<double> X(D);

    for (int i = 0; i < D; i++) {
    unsigned long long int sum =0;

        //we extract the chromosome and start decoding it
        for (int j = 0; j < l; j++) {
            sum = (sum << 1) | bitstring[i*l +j];
        }

        //we normalize it in [0,1]
        double div = pow(2,l)-1;
        double x = double(sum) / div;

        //now X[i] is within [low,high]
        X[i] = low + x*(high-low);
    }
    return X;
}

//================================= MATH FUNCTIONS ===========================================

double rastrigin(vector<double> X){
    size_t n = X.size();
    double sum = 10*n;
    for(int ii= 0; ii< n; ++ii){
        sum += pow(X[ii],2) -10 *cos(2 * PI * X[ii]);
    }
    return sum;
}

double deJong(vector<double> X){
    size_t n = X.size();
    double sum =0;
    for(int ii = 0; ii < n; ++ii){
        sum +=pow(X[ii],2);
    }
    return sum;
}

double schwefel(vector<double> X){
    float n = X.size();
    double sum =0;
    for(int ii = 0; ii<n;++ii){
        sum +=(-X[ii] * sin((sqrt(fabs(X[ii])))));
    }
    return sum;
}

double michalewicz(vector<double> X){
    float n = X.size();
    double sum = 0;
    for(int ii = 0; ii< n; ++ii){
        sum += sin(X[ii])*(pow(sin(ii*pow(X[ii],2))/PI,20));
    }
    sum*= -1;
    return sum;
}

//========================================= HILL CLIMB ===========================================
struct HCResult {
    vector<bool> X;
    double best_val;
    int nrIterations;
    double time;
};

struct bounds {
    float low,high;
    string name;
};

HCResult HillClimbing(const function<double(const vector<double>&)>& f,vector<bool>vc , int D,int l, bounds &bd, string& variant, int max_iterations = 100000) {

    int pos =-1;
    int iter=0;
    bool improved = true;
    double eval_current = f(decode_Bits(vc,D,l,bd.low,bd.high));
    auto begin  = chrono::high_resolution_clock::now();

    double best_neighbour = numeric_limits<double>::infinity(); //we generate the neighbour. the first one is +infinity;

    while (improved && iter < max_iterations) {
        improved= false;
        iter++;

        for (int ii= 0;ii<vc.size();++ii) {

            vc[ii]=!vc[ii];
            float eval = f((decode_Bits(vc,D,l,bd.low,bd.high)));
            vc[ii]= !vc[ii];//we return it to the current state

            if (variant == "first" && eval < best_neighbour ) {
                eval_current = eval;
                best_neighbour = eval;
                vc[ii]=!vc[ii];
                improved = true;
                break;
            };

            if (variant == "best" && eval < eval_current) {
                best_neighbour = eval;
                pos = ii;
            };

            if (variant == "worst" && eval < eval_current) {
                if (pos ==-1 || eval < best_neighbour) {
                    best_neighbour = eval;
                    pos = ii;
                }
            };
        }
        if (variant != "first" && pos !=-1 && best_neighbour <eval_current ) {
        //if we managed to find a global improvement(for best or worst)
            eval_current = best_neighbour;
            vc[pos]=!vc[pos];
            improved = true;
        }

    }
    auto end = chrono::high_resolution_clock::now();
    float elapsed = chrono::duration_cast<chrono::duration<float>>(end - begin).count();
    return {vc,eval_current,iter,elapsed};
}
//========================================== SIMULATED ANNEALING ===================================
struct sa_result {
    vector<bool> X;
    double best_val;
    int nrIterations;
    double time;
};

sa_result SimulatedAnnealing(const function<double(const vector<double>&)>& f,
                             vector<bool> vc, int D, int l,
                             bounds &bd, int max_iterations = 100000)
{
    float T = 1000.0f;
    const float T_min = 1e-5;
    const float alpha = 0.99f;
    int iters = 0;

    mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
    uniform_real_distribution<float> urd(0.0f, 1.0f);
    uniform_int_distribution<int> uid(0, D * l - 1);

    auto begin = chrono::high_resolution_clock::now();

    float eval_vc = f(decode_Bits(vc, D, l, bd.low, bd.high));
    float best_val = eval_vc;
    vector<bool> best_sol = vc;

    while (T > T_min && iters < max_iterations) {
        vector<bool> neighbour = vc;
        int num_flips = 1 + rng() % 2; // random flip of 1,2 or 3 bits
        for (int i = 0; i < num_flips; i++) {
            int flip_pos = uid(rng);
            neighbour[flip_pos] = !neighbour[flip_pos];
        }

        double eval_neighbour = f(decode_Bits(neighbour, D, l, bd.low, bd.high));
        double delta = eval_neighbour - eval_vc;

        if (delta < 0) {
            vc = neighbour;
            eval_vc = eval_neighbour;
        } else {
            float p = exp(-delta / T);
            if (urd(rng) < p) {
                vc = neighbour;
                eval_vc = eval_neighbour;
            }
        }

        if (eval_vc < best_val) {
            best_val = eval_vc;
            best_sol = vc;
        }

        T *= alpha;
        ++iters;
    }

    auto end = chrono::high_resolution_clock::now();
    double elapsed = chrono::duration<double>(end - begin).count();

    return {best_sol, best_val, iters, elapsed};
}
//========================================== MAIN ====================
int main(){
    int D = 30;
    int l = 20;
    int i = 0;
    vector<bool> X =random_bits(D,l);
    bounds bd={-5.12,5.12};
    vector<string> variant = {"first","best","worst"};

    sa_result saResult = SimulatedAnnealing(rastrigin,X,D,l,bd);
    cout<<"Simulated Annealing: "<<endl;
    cout<<"Minimum Value: " << saResult.best_val<<endl;
    cout<<"Nr of iterations: "<<saResult.nrIterations<<endl;
    cout<<"Time: "<<saResult.time<<endl<<endl;
/*
    HCResult results= HillClimbing(rastrigin,X,D,l,bd,variant[i]);
    cout<<"Hill Climbing with variant : "<<variant[i]<<endl;
    cout<<"Minimum Value: "<< results.best_val<<endl;
    cout<<"Nr of iterations: "<< results.nrIterations << endl;
    cout<<"Average Time: "<< results.time<< endl;
*/

    return 0;
}