#include <iostream>
#include <cmath>
#include <cstring>
#include <vector>

using namespace std;
float a = -5.12; float b = 5.12;
const float PI = 3.14;

float decode(vector<bool> v){
    unsigned long long int sum = 0;
    unsigned long long int mul = 1;
    for(int ii =v.size()-1;ii>=0;--ii){
        sum += v[ii] * mul;
        mul <<= 1;
    }
    float div = pow( 2,v.size()-1);
    float x = float(sum) / div;
    x= a + x*(b-a);
    return x;
}

vector<bool> cstringToVect(char* s){
    vector<bool>v;
    int l = strlen(s);
    for(int  ii =0 ;ii< l;++ii){
        v.push_back(s[ii] - '0');
    }
    return v;
}

float rastrigin(vector<float> X){
    int n = X.size();
    float sum = 10*n;
    for(int ii= 0; ii< n; ++ii){
        sum += pow(X[ii],2) -10 *cos(2 * PI * X[ii]);
    }
    return sum;
}

float deJong(vector<float> X){
    float n = X.size();
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
        sum +=(-X[ii] * sin((sqrt(abs(X[ii])))));
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

int main(int argc, char* argv[]){
    vector<bool> v = cstringToVect(argv[1]);
    float x = decode(v);
    vector<float> X = {x};
    cout << deJong(X);
}