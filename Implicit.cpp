#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

vector<vector<double>> implicitScheme(double Vol, double Int_Rate, char PType, double Strike, double Expiration, double S_max, int NAS, int NTS) {
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

    vector<double> a(NAS + 1, 0.0), b(NAS + 1, 1.0), c(NAS + 1, 0.0), d(NAS + 1, 0.0);
    
    for (int k = NTS - 1; k >= 0; --k) {
        for (int i = 1; i < NAS; ++i) {
            a[i] = -0.5 * (Vol * Vol * i * i - i * (Int_Rate - D)) * dt;
            b[i] = 1 + (Int_Rate + Vol * Vol * i * i) * dt;
            c[i] = -0.5 * (Vol * Vol * i * i + i * (Int_Rate - D)) * dt;
            d[i] = V[i][k + 1];
        }

        for (int i = 2; i < NAS; ++i) {
            b[i] -= a[i] * c[i - 1] / b[i - 1];
            d[i] -= a[i] * d[i - 1] / b[i - 1];
        }

        V[NAS - 1][k] = d[NAS - 1] / b[NAS - 1];
        for (int i = NAS - 2; i >= 1; --i) {
            V[i][k] = (d[i] - c[i] * V[i + 1][k]) / b[i];
        }

        V[0][k] = V[0][k + 1] * (1 - Int_Rate * dt);
        V[NAS][k] = 2 * V[NAS - 1][k] - V[NAS - 2][k];
    }

    return V;
}

// Function implementing the Jacobi method
vector<vector<double>> jacobiScheme(double Vol, double Int_Rate, char PType, double Strike, double Expiration, double S_max, int NAS, int NTS) {
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

    vector<double> a(NAS + 1, 0.0), b(NAS + 1, 1.0), c(NAS + 1, 0.0), d(NAS + 1, 0.0);
    
    const double tol = 1e-6;
    const int max_iter = 1000;

    for (int k = NTS - 1; k >= 0; --k) {
        for (int i = 1; i < NAS; ++i) {
            a[i] = -0.5 * (Vol * Vol * i * i - i * (Int_Rate - D)) * dt;
            b[i] = 1 + (Int_Rate + Vol * Vol * i * i) * dt;
            c[i] = -0.5 * (Vol * Vol * i * i + i * (Int_Rate - D)) * dt;
            d[i] = V[i][k + 1];
        }

        // Initialize the guess solution with the previous time step
        vector<double> V_new(NAS + 1, 0.0);
        for (int i = 1; i < NAS; ++i) {
            V_new[i] = V[i][k + 1];  // Initial guess
        }

        // Jacobi Iteration
        int iter = 0;
        double error;
        do {
            vector<double> V_old = V_new;
            error = 0.0;

            for (int i = 1; i < NAS; ++i) {
                V_new[i] = (d[i] - a[i] * V_old[i - 1] - c[i] * V_old[i + 1]) / b[i];
                error = max(error, fabs(V_new[i] - V_old[i]));
            }

            iter++;
            if (iter > max_iter) break;  // Avoid infinite loop

        } while (error > tol);

        // Store the result in V
        for (int i = 1; i < NAS; ++i) {
            V[i][k] = V_new[i];
        }

        // Boundary conditions
        V[0][k] = V[0][k + 1] * (1 - Int_Rate * dt);
        V[NAS][k] = 2 * V[NAS - 1][k] - V[NAS - 2][k];
    }

    return V;
}

vector<vector<double>> gaussSeidelScheme(double Vol, double Int_Rate, char PType, double Strike, double Expiration, double S_max, int NAS, int NTS) {
    vector<double> S(NAS + 1, 0.0);
    double dS = S_max / NAS;
    double dt = Expiration / NTS;
    double D = 0; // Dividend yield (set to 0 if not given)

    vector<vector<double>> V(NAS + 1, vector<double>(NTS + 1, 0.0));
    int q = (PType == 'P') ? -1 : 1;  // Put or Call option

    // Initialize stock price grid and terminal condition
    for (int i = 0; i <= NAS; ++i) {
        S[i] = i * dS;
        V[i][NTS] = max(q * (S[i] - Strike), 0.0);
    }

    vector<double> a(NAS + 1, 0.0), b(NAS + 1, 1.0), c(NAS + 1, 0.0), d(NAS + 1, 0.0);
    
    const double tol = 1e-6;
    const int max_iter = 1000;

    for (int k = NTS - 1; k >= 0; --k) {
        // Coefficients for the tridiagonal system
        for (int i = 1; i < NAS; ++i) {
            a[i] = -0.5 * (Vol * Vol * i * i - i * (Int_Rate - D)) * dt;
            b[i] = 1 + (Int_Rate + Vol * Vol * i * i) * dt;
            c[i] = -0.5 * (Vol * Vol * i * i + i * (Int_Rate - D)) * dt;
            d[i] = V[i][k + 1];  // Right-hand side (previous time step)
        }

        // Gauss-Seidel Iteration
        int iter = 0;
        double error;
        do {
            error = 0.0;

            for (int i = 1; i < NAS; ++i) {
                double old_V = V[i][k];

                // Using updated values immediately (Gauss-Seidel formula)
                V[i][k] = (d[i] - a[i] * V[i - 1][k] - c[i] * V[i + 1][k]) / b[i];

                // Compute error to check for convergence
                error = max(error, fabs(V[i][k] - old_V));
            }

            iter++;
            if (iter > max_iter) break;  // Avoid infinite loop

        } while (error > tol);

        // Boundary conditions
        V[0][k] = V[0][k + 1] * (1 - Int_Rate * dt); // Dirichlet condition
        V[NAS][k] = 2 * V[NAS - 1][k] - V[NAS - 2][k]; // Linear extrapolation
    }

    return V;
}

int main() {
    double S0, E, T, sigma, r, tol;
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
    cout << "Enter Convergence Tolerance (tol): ";
    cin >> tol;
    
    double S_max = 2 * S0;
    int NAS = 1000;
    int NTS = 1000;
    int idx = static_cast<int>(S0 / (S_max / NAS));
    
    vector<vector<double>> put_values_jacobi = jacobiScheme(sigma, r, 'P', E, T, S_max, NAS, NTS);
    vector<vector<double>> call_values_jacobi = jacobiScheme(sigma, r, 'C', E, T, S_max, NAS, NTS);

    cout << "Put Option Price (Jacobi Scheme): " << put_values_jacobi[idx].front() << endl;
    cout << "Call Option Price (Jacobi Scheme): " << call_values_jacobi[idx].front() << endl;

    vector<vector<double>> put_values_gauss = gaussSeidelScheme(sigma, r, 'P', E, T, S_max, NAS, NTS);
    vector<vector<double>> call_values_gauss = gaussSeidelScheme(sigma, r, 'C', E, T, S_max, NAS, NTS);

    cout << "Put Option Price (Gauss-Seidel Scheme): " << put_values_gauss[idx].front() << endl;
    cout << "Call Option Price (Gauss-Seidel Scheme): " << call_values_gauss[idx].front() << endl;
    
    vector<vector<double>> put_values_implicit = implicitScheme(sigma, r, 'P', E, T, S_max, NAS, NTS);
    vector<vector<double>> call_values_implicit = implicitScheme(sigma, r, 'C', E, T, S_max, NAS, NTS);

    cout << "Put Option Price (Implicit Scheme): " << put_values_implicit[idx].front() << endl;
    cout << "Call Option Price (Implicit Scheme): " << call_values_implicit[idx].front() << endl;
    
    return 0;
}


/*vector<vector<double>> sorScheme(double Vol, double Int_Rate, char PType, double Strike, double Expiration, double S_max, int NAS, int NTS, double tol) {
    vector<double> S(NAS + 1, 0.0);
    double dS = S_max / NAS;
    double dt = Expiration / NTS;
    double D = 0;
    double omega = 1.0;
    int prevNoIts = 0;
    
    vector<vector<double>> V(NAS + 1, vector<double>(NTS + 1, 0.0));
    int q = (PType == 'P') ? -1 : 1;

    for (int i = 0; i <= NAS; ++i) {
        S[i] = i * dS;
        V[i][NTS] = max(q * (S[i] - Strike), 0.0);
    }

    vector<double> a(NAS + 1, 0.0), b(NAS + 1, 1.0), c(NAS + 1, 0.0), d(NAS + 1, 0.0), temp(NAS + 1, 0.0);
    
    for (int k = NTS - 1; k >= 0; --k) {
        for (int i = 1; i < NAS; ++i) {
            a[i] = -0.5 * (Vol * Vol * i * i - i * (Int_Rate - D)) * dt;
            b[i] = 1 + (Int_Rate + Vol * Vol * i * i) * dt;
            c[i] = -0.5 * (Vol * Vol * i * i + i * (Int_Rate - D)) * dt;
            d[i] = V[i][k + 1];
        }

        double Err = tol + 1;
        int NoIts = 0;
        while (Err > tol) { // Ensure max iterations to prevent infinite loops
            Err = 0;
            for (int i = 1; i < NAS; ++i) {
                double newVal = V[i][k] + omega * (d[i] - c[i] * V[i + 1][k] - b[i] * V[i][k] - a[i] * V[i - 1][k]) / b[i];
                Err += pow(newVal - V[i][k], 2);
                V[i][k] = newVal;
            }
            NoIts++;
        }

        if (k == NTS - 1 || NoIts < prevNoIts) {
            omega += 0.05;
        } else if (NoIts > prevNoIts) {
            omega -= 0.05;
        }
        prevNoIts = NoIts;

        V[0][k] = max(Strike - S[0], 0.0);  // Fixing put option boundary
        V[NAS][k] = 0;  // Deep ITM put option value should be 0 at large S
    }

    return V;
}
    */