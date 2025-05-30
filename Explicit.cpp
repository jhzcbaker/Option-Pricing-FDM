#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

vector<vector<double>> explicitScheme(double Vol, double Int_Rate, char PType, double Strike, double Expiration, double S_max, char EType, int NAS) {
    vector<double> S(NAS + 1, 0.0);
    double dS = S_max / NAS;
    double dt = 0.9 / (Vol * Vol * NAS * NAS);
    int NTS = static_cast<int>(Expiration / dt) + 1;
    dt = Expiration / NTS;
    double D = 0;

    vector<vector<double>> V(NAS + 1, vector<double>(NTS + 1, 0.0));
    int q = (PType == 'P') ? -1 : 1;

    for (int i = 0; i <= NAS; ++i) {
        S[i] = i * dS;
        V[i][NTS] = max(q * (S[i] - Strike), 0.0);
    }


    for (int k = NTS - 1; k >= 0; --k) {
        for (int i = 1; i < NAS; ++i) {
            double A = 0.5 * ( Vol * Vol * i * i - i * (Int_Rate - D)) * dt;
            double B = 1 - (Int_Rate + Vol * Vol * i * i) * dt;
            double C = 0.5 * ( Vol * Vol * i * i + i * (Int_Rate - D)) * dt;
            V[i][k] = A * V[i - 1][k + 1] + B * V[i][k + 1] + C * V[i + 1][k + 1];
        }
        V[0][k] = V[0][k + 1] * (1 - Int_Rate * dt);
        V[NAS][k] = 2 * V[NAS - 1][k] - V[NAS - 2][k];
    }
    return V;
}






int main() {
    double S0, E, T, sigma, r;
    cout << "Enter Stock Price (S0): ";
    cin >> S0;
    cout << "Enter Strike Price (E): ";
    cin >> E;
    cout << "Enter Time to Expiry (T in years): ";
    cin >> T;
    cout << "Enter Volatility (sigma as decimal): ";
    cin >> sigma;
    cout << "Enter Risk-free Interest Rate (r as decimal): ";
    cin >> r;
    
    double S_max = 2 * E;
    int NAS = 100;
    
    vector<vector<double>> call_values_explicit = explicitScheme(sigma, r, 'C', E, T, S_max, 'N', NAS);
    vector<vector<double>> put_values_explicit = explicitScheme(sigma, r, 'P', E, T, S_max, 'Y', NAS);
   
    
    int idx = static_cast<int>(S0 / (S_max / NAS));
    cout << "\nCall Option Price (Explicit Scheme): " << call_values_explicit[idx].front() << endl;
    cout << "Put Option Price (Explicit Scheme): " << put_values_explicit[idx].front() << endl;

    
    return 0;
}

    /*for (int k = 1; k <= NTS; ++k) {
        for (int i = 1; i < NAS; ++i) {
            double Delta = (V[i + 1][k - 1] - V[i - 1][k - 1]) / (2 * dS);
            double Gamma = (V[i + 1][k - 1] - 2 * V[i][k - 1] + V[i - 1][k - 1]) / (dS * dS);
            double Theta = -0.5 * Vol * Vol * S[i] * S[i] * Gamma - Int_Rate * S[i] * Delta + Int_Rate * V[i][k - 1];
            V[i][k] = V[i][k - 1] - dt * Theta;
        }
        V[0][k] = V[0][k - 1] * (1 - Int_Rate * dt);
        V[NAS][k] = 2 * V[NAS - 1][k] - V[NAS - 2][k];
    }
    */