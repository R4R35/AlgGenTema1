#include <iostream>
#include <cmath>
#include <cstring>
#include <vector>
#include <limits>
#include <random>
#include <bitset>
#include <functional>

using namespace std;
// =========================== CONSTANTS =======================================


const float PI = 2*asin(1.0);


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

vector<float> decode_Bits(vector<bool> bitstring,int D,int l, float low,float high) {
    vector<float> X(D);

    for (int i = 0; i < D; i++) {
    unsigned long long int sum =0;

        //we extract the chromosome and start decoding it
        for (int j = 0; j < l; j++) {
            sum = (sum << 1) | bitstring[i*l +j];
        }

        //we normalize it in [0,1]
        float div = pow(2,l)-1;
        float x = float(sum) / div;

        //now X[i] is within [low,high]
        X[i] = low + x*(high-low);
    }
    return X;
}

//================================= MATH FUNCTIONS ===========================================

float rastrigin(vector<float> X){
    size_t n = X.size();
    float sum = 10*n;
    for(int ii= 0; ii< n; ++ii){
        sum += pow(X[ii],2) -10 *cos(2 * PI * X[ii]);
    }
    return sum;
}

float deJong(vector<float> X){
    size_t n = X.size();
    float sum =0;
    for(int ii = 0; ii < n; ++ii){
        sum +=pow(X[ii],2);
    }
    return sum;
}

float schwefel(vector<float> X){
    float n = X.size();
    float sum =0;
    for(int ii = 0; ii<n;++ii){
        sum +=(-X[ii] * sin((sqrt(fabs(X[ii])))));
    }
    return sum;
}

float michalewicz(vector<float> X){
    float n = X.size();
    float sum = 0;
    for(int ii = 0; ii< n; ++ii){
        sum += sin(X[ii])*(pow(sin(ii*pow(X[ii],2))/PI,20));
    }
    sum*= -1;
    return sum;
}

//========================================= HILL CLIMB ===========================================
struct HCResult {
    vector<bool> X;
    float best_val;
    int nrIterations;
    //double time;
};

struct bounds {
    float low,high;
    string name;
};

HCResult HillClimbing(const function<float(const vector<float>&)>& f,vector<bool>vc , int D,int l, bounds &bd, string& variant, int max_iterations = 100000) {

    int pos =-1;
    int iter=0;
    bool improved = true;
    float eval_current = f(decode_Bits(vc,D,l,bd.low,bd.high));

    while (improved && iter < max_iterations) {
        improved= false;
        iter++;

        float best_neighbour = numeric_limits<float>::infinity(); //we generate the neighbour. the first one is -infinity;
        for (int ii= 0;ii<vc.size();++ii) {

            vc[ii]=!vc[ii];
            float eval = f((decode_Bits(vc,D,l,bd.low,bd.high)));
            vc[ii]= !vc[ii];//we return it to the current state

            if (variant == "first" && eval < best_neighbour ) {
                eval_current = eval;
                vc[ii]=!vc[ii];
                improved = true;
                break;
            };

            if (variant == "best" && eval < eval_current) {
                best_neighbour = eval;
                pos = ii;
            };

            if (variant == "worst" && eval < eval_current) {
                if (pos ==-1 || eval > best_neighbour) {
                    best_neighbour = eval;
                    pos = ii;
                }
            };
        }
        if (variant != "first") {
        //if we managed to find a global improvement(global or worst)
            if (pos !=-1 && best_neighbour <eval_current ) {
                eval_current = best_neighbour;
                vc[pos]=!vc[pos];
                improved = true;
            }
        }

    }
    return {vc,eval_current,iter};
}


//========================================== MAIN ====================
int main(){
    vector<bool> X =random_bits(10,20);
    bounds bd={-5.12,5.12};
    vector<string> variant = {"first","best","worst"};
    HCResult results= HillClimbing(rastrigin,X,10,20,bd,variant[0]);
    cout<<"Minimum Value: "<< results.best_val<<endl;
    cout<<"Nr of iterations: "<< results.nrIterations << endl;
    return 0;
}