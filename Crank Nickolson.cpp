#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;


// Function implementing the Crank-Nicolson scheme
vector<vector<double>> crankNicolsonScheme(double Vol, double Int_Rate, char PType, double Strike, double Expiration, double S_max, int NAS, int NTS) {
    vector<double> S(NAS + 1, 0.0);
    double dS = S_max / NAS; // Step size in asset price space
    double dt = Expiration / NTS; // Step size in time
    double D = 0;

    // Matrix to store option prices
    vector<vector<double>> V(NAS + 1, vector<double>(NTS + 1, 0.0));
    int q = (PType == 'P') ? -1 : 1;

    // Initialize the asset prices and final option values
    for (int i = 0; i <= NAS; ++i) {
        S[i] = i * dS;
        V[i][NTS] = max(q * (S[i] - Strike), 0.0); // Payoff at expiry
    }

    // Coefficients for the tridiagonal matrix
    vector<double> a(NAS + 1, 0.0), b(NAS + 1, 1.0), c(NAS + 1, 0.0), d(NAS + 1, 0.0);
    vector<double> L(NAS + 1, 0.0), U(NAS + 1, 0.0), Y(NAS + 1, 0.0);

    // Iterate backwards in time
    for (int k = NTS - 1; k >= 0; --k) {
        for (int i = 1; i < NAS; ++i) {
            double sigmaSq = Vol * Vol * i * i;
            double rPart = (Int_Rate - D) * i;

            a[i] = -0.25 * dt * (sigmaSq - rPart);
            b[i] = 1 + 0.5 * dt * (sigmaSq + Int_Rate);
            c[i] = -0.25 * dt * (sigmaSq + rPart);
            
            // Right-hand side values using Crank-Nicolson averaging
            d[i] = (1 - 0.5 * dt * (sigmaSq + Int_Rate)) * V[i][k + 1] 
                   + 0.25 * dt * (sigmaSq - rPart) * V[i - 1][k + 1] 
                   + 0.25 * dt * (sigmaSq + rPart) * V[i + 1][k + 1];
        }
        
        // Forward elimination step in the tridiagonal matrix algorithm
        L[1] = b[1];
        U[1] = c[1] / L[1];
        Y[1] = d[1] / L[1];

        for (int i = 2; i < NAS; ++i) {
            L[i] = b[i] - a[i] * U[i - 1];
            U[i] = c[i] / L[i];
            Y[i] = (d[i] - a[i] * Y[i - 1]) / L[i];
        }

        // Back substitution step
        V[NAS - 1][k] = Y[NAS - 1];
        for (int i = NAS - 2; i >= 1; --i) {
            V[i][k] = Y[i] - U[i] * V[i + 1][k];
        }

        // Boundary conditions
        V[0][k] = V[0][k + 1] * (1 - Int_Rate * dt);
        V[NAS][k] = 2 * V[NAS - 1][k] - V[NAS - 2][k];
    }

    return V;
}


vector<vector<double>> crankNicolson(double Vol, double Int_Rate, char PType, double Strike, double Expiration, double S_max, int NAS, int NTS) {
    vector<double> S(NAS + 1, 0.0);
    double dS = S_max / NAS;
    double dt = Expiration / NTS;
    double D = 0;

    vector<vector<double>> V(NAS + 1, vector<double>(NTS + 1, 0.0));
    int q = (PType == 'P') ? -1 : 1;

    for (int i = 0; i <= NAS; ++i) {
        S[i] = i * dS;
        V[i][NTS] = max(q * (S[i] - Strike), 0.0);
    }

    vector<double> a(NAS, 0.0), b(NAS, 1.0), c(NAS, 0.0), d(NAS, 0.0);
    
    for (int k = NTS - 1; k >= 0; --k) {
        for (int i = 1; i < NAS; ++i) {
            double alpha = 0.25 * dt * (Vol * Vol * i * i - Int_Rate * i);
            double beta = -0.5 * dt * (Vol * Vol * i * i + Int_Rate);
            double gamma = 0.25 * dt * (Vol * Vol * i * i + Int_Rate * i);
            
            a[i - 1] = -alpha;
            b[i - 1] = 1 - beta;
            c[i - 1] = -gamma;
            d[i - 1] = alpha * V[i - 1][k + 1] + (1 + beta) * V[i][k + 1] + gamma * V[i + 1][k + 1];
        }
        
        for (int i = 1; i < NAS - 1; ++i) {
            double m = a[i] / b[i - 1];
            b[i] -= m * c[i - 1];
            d[i] -= m * d[i - 1];
        }
        
        V[NAS - 1][k] = d[NAS - 2] / b[NAS - 2];
        for (int i = NAS - 2; i > 0; --i) {
            V[i][k] = (d[i - 1] - c[i - 1] * V[i + 1][k]) / b[i - 1];
        }
    }

    return V;
}




int main() {
    double S0, E, T, sigma, r;
    int NAS = 1000, NTS = 1000;
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
    double S_max = 2 * S0;
    vector<vector<double>> call_values_crank_nicolson = crankNicolsonScheme(sigma, r, 'C', E, T, S_max, NAS, NTS);
    vector<vector<double>> put_values_crank_nicolson = crankNicolsonScheme(sigma, r, 'P', E, T, S_max, NAS, NTS);
  
    int idx = static_cast<int>(S0 / (S_max / NAS));

    cout << "\nPut Option Price (Explicit Scheme): " << put_values_crank_nicolson[idx].front() << endl;
    cout << "Call Option Price (Explicit Scheme): " << call_values_crank_nicolson[idx].front() << endl;


    vector<vector<double>> call_values_crank = crankNicolson(sigma, r, 'C', E, T, S_max, NAS, NTS);
    vector<vector<double>> put_values_crank = crankNicolson(sigma, r, 'P', E, T, S_max, NAS, NTS);

    cout << "\nPut Option Price (Explicit Scheme): " << put_values_crank[idx].front() << endl;
    cout << "Call Option Price (Explicit Scheme): " << call_values_crank[idx].front() << endl;


    return 0;
}
