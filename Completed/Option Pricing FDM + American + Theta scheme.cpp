#include <iostream>
#include <vector>
#include <cmath>


using namespace std;

// Explicit Scheme
vector<vector<double>> explicitScheme(double Vol, double Int_Rate, char PType, double Strike, double Expiration, double S_max, int NAS, int NTS, double D, char EType) {
        // Define stock price grid
        vector<double> S(NAS + 1, 0.0);
        double dS = S_max / NAS; // Step size for stock prices


        /*
        // This reparametrises dt to ensure stability, method from Wilmotts on Quantitative Finance
        // Compute time step size with stability condition
        double dt = 0.9 / (Vol * Vol * NAS * NAS);
        NTS = static_cast<int>(Expiration / dt) + 1; // Number of time steps
        dt = Expiration / NTS; // Adjust time step size to fit exactly in Expiration
        */


        // Initialize option value matrix
        vector<vector<double>> V(NAS + 1, vector<double>(NTS + 1, 0.0));
        int q = (PType == 'P') ? -1 : 1; // Determine if the option is a Put (-1) or Call (1)

        
        double dt = Expiration / NTS;
        if (dt > (1 / (Vol * Vol * NAS * NAS))) {
            cout << "Stabitity Condition Violated (dt <= 1 / (Sigma^2 * N^2))" << endl;
            return V;
        }
        
    
        // Set up stock prices and terminal condition (payoff at expiration)
        for (int i = 0; i <= NAS; ++i) {
            S[i] = i * dS;
            V[i][NTS] = max(q * (S[i] - Strike), 0.0); // Payoff function
        }
    
        // Iterate backward in time to compute option values
        for (int k = NTS - 1; k >= 0; --k) {
    
            /*      Other Boundary Conditions
            if (PType == 'P') {
                V[0][k] = V[0][k + 1] * (1 - Int_Rate * dt); // Put option boundary condition
                V[NAS][k] = 0.0;
            } else{
                V[0][k] = 0.0; // Call option lower boundary condition
                V[NAS][k] = S_max * exp(D - k * dt) - Strike * exp(-Int_Rate * k * dt);
            }
            

            // Boundary conditions from Wilmotts on Quantitative Finance
            V[0][k] = V[0][k + 1] * (1 - Int_Rate * dt);
            V[NAS][k] = 2 * V[NAS - 1][k] - V[NAS - 2][k];
            
            */
    
            // Boundary conditions
            V[0][k] = V[0][k + 1] * (1 - Int_Rate * dt);
            V[NAS][k] = -(NAS) * (Int_Rate - D) * dt * V[NAS - 1][k + 1] - (1 + (NAS * (Int_Rate - D) - Int_Rate) * dt) * V[NAS][k + 1];
            
    
            for (int i = 1; i < NAS; ++i) {
                // Coefficients for finite difference scheme
                double sigmaSq = Vol * Vol * i * i;
                double rPart = (Int_Rate - D) * i;
                
                double A = 0.5 * dt * (sigmaSq - rPart);
                double B = 1 - dt * (sigmaSq + Int_Rate);
                double C = 0.5 * dt * (sigmaSq + rPart);
                
                // Compute option value at grid point (i, k)
                V[i][k] = A * V[i - 1][k + 1] + B * V[i][k + 1] + C * V[i + 1][k + 1];
            }

        // Early exercise condition for American options
        if (EType == 'Y') {
            for (int i = 1; i < NAS; ++i) {
                V[i][k] = max(V[i][k], V[i][NTS]); // V[i][NTS] is set equal to the payoff funtion above when terminal conditions are set on line 38
            }
        }

    }

    return V;
}

// Implicit Scheme
vector<vector<double>> implicitScheme(double Vol, double Int_Rate, char PType, double Strike, double Expiration, double S_max, int NAS, int NTS, double D, char EType, int maxIter = 1000, double tol = 1e-6){
    vector<double> S(NAS + 1, 0.0);
    double dS = S_max / NAS;
    double dt = Expiration / NTS;

    vector<vector<double>> V(NAS + 1, vector<double>(NTS + 1, 0.0));
    int q = (PType == 'P') ? -1 : 1;  // Put or Call option

    // Initialize stock price grid and terminal condition
    for (int i = 0; i <= NAS; ++i) {
        S[i] = i * dS;
        V[i][NTS] = max(q * (S[i] - Strike), 0.0);
    }

    // Coefficients for the tridiagonal system
    vector<double> a(NAS), b(NAS), c(NAS), d(NAS);

    // Iterate backward in time
    for (int k = NTS - 1; k >= 0; --k) {
        // Boundary conditions
        V[0][k] = V[0][k + 1] * (1 - Int_Rate * dt);
        V[NAS][k] = -(NAS) * (Int_Rate - D) * dt * V[NAS - 1][k + 1] - (1 + (NAS * (Int_Rate - D) - Int_Rate) * dt) * V[NAS][k + 1];


        // Construct the tridiagonal system
        for (int i = 1; i < NAS; ++i) {
            double sigmaSq = Vol * Vol * i * i;
            double rPart = (Int_Rate - D) * i;

            a[i] = -0.5 * dt * (sigmaSq - rPart);
            b[i] = 1 + dt * (sigmaSq + Int_Rate);
            c[i] = -0.5 * dt * (sigmaSq + rPart);
            d[i] = V[i][k + 1];  // Right-hand side (previous time step)
        }

        // Apply Gauss-Seidel Iteration
        double error = 1.0;
        int iter = 0;

        while (error > tol && iter < maxIter) {
            error = 0.0;

            for (int i = 1; i < NAS; ++i) {
                double oldV = V[i][k];

                V[i][k] = (d[i] - a[i] * V[i - 1][k] - c[i] * V[i + 1][k]) / b[i];

                error = max(error, fabs(V[i][k] - oldV)); // Compute max error for convergence check
            }

            iter++;
        }
        
        // Early exercise condition for American options
        if (EType == 'Y') {
            for (int i = 1; i < NAS; ++i) {
                V[i][k] = max(V[i][k], V[i][NTS]); // V[i][NTS] is set equal to the payoff funtion above when terminal conditions are set
            }
        }

    }

    return V;
}



// Crank-Nicolson Scheme
vector<vector<double>> crankNicolson(double Vol, double Int_Rate, char PType, double Strike, double Expiration, double S_max, int NAS, int NTS, double D, char EType, int maxIter = 1000, double tol = 1e-6) {
    vector<double> S(NAS + 1, 0.0);
    double dS = S_max / NAS; // Step size in asset price space
    double dt = Expiration / NTS; // Step size in time

    // Matrix to store option prices
    vector<vector<double>> V(NAS + 1, vector<double>(NTS + 1, 0.0));
    int q = (PType == 'P') ? -1 : 1;

    // Initialize the asset prices and final option values
    for (int i = 0; i <= NAS; ++i) {
        S[i] = i * dS;
        V[i][NTS] = max(q * (S[i] - Strike), 0.0); // Payoff at expiry
    }

    // Coefficients for the Crank-Nicolson scheme
    vector<double> a(NAS), b(NAS), c(NAS), d(NAS);

    // Iterate backwards in time
    for (int k = NTS - 1; k >= 0; --k) {

        // Boundary conditions
        V[0][k] = V[0][k + 1] * (1 - Int_Rate * dt);
        V[NAS][k] = -(NAS) * (Int_Rate - D) * dt * V[NAS - 1][k + 1] - (1 + (NAS * (Int_Rate - D) - Int_Rate) * dt) * V[NAS][k + 1];


        // Compute coefficients for the Crank-Nicolson scheme
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

        // Apply Gauss-Seidel Iteration
        double error = 1.0;
        int iter = 0;

        while (error > tol && iter < maxIter) {
            error = 0.0;

            for (int i = 1; i < NAS; ++i) {
                double oldV = V[i][k];

                V[i][k] = (d[i] - a[i] * V[i - 1][k] - c[i] * V[i + 1][k]) / b[i];

                error = max(error, fabs(V[i][k] - oldV)); // Compute max error for convergence check
            }

            iter++;
        }

        // Early exercise condition for American options
        if (EType == 'Y') {
            for (int i = 1; i < NAS; ++i) {
                V[i][k] = max(V[i][k], V[i][NTS]); // V[i][NTS] is set equal to the payoff funtion above when terminal conditions are set
            }

        }

    }

    return V;
}


// Crank-Nicolson Theta Scheme
vector<vector<double>> crankNicolsonTheta(double Vol, double Int_Rate, char PType, double Strike, double Expiration, double S_max, int NAS, int NTS, double D, char EType, double theta, char thetaType, int maxIter = 1000, double tol = 1e-6) {
    vector<double> S(NAS + 1, 0.0);
    double dS = S_max / NAS; // Step size in asset price space
    double dt = Expiration / NTS; // Step size in time

    // Matrix to store option prices
    vector<vector<double>> V(NAS + 1, vector<double>(NTS + 1, 0.0));
    int q = (PType == 'P') ? -1 : 1;

    // Initialize the asset prices and final option values
    for (int i = 0; i <= NAS; ++i) {
        S[i] = i * dS;
        V[i][NTS] = max(q * (S[i] - Strike), 0.0); // Payoff at expiry
    }

    // Coefficients for the Crank-Nicolson scheme
    vector<double> a(NAS), b(NAS), c(NAS), d(NAS);

    // Iterate backwards in time
    for (int k = NTS - 1; k >= 0; --k) {

        // Boundary conditions
        V[0][k] = V[0][k + 1] * (1 - Int_Rate * dt);
        V[NAS][k] = -(NAS) * (Int_Rate - D) * dt * V[NAS - 1][k + 1] - (1 + (NAS * (Int_Rate - D) - Int_Rate) * dt) * V[NAS][k + 1];


        // Compute coefficients for the Crank-Nicolson scheme
        for (int i = 1; i < NAS; ++i) {
            if (thetaType == 'Y') {
                theta = 0.5 + 1.0 / 12.0 * S[i] * dS * dt; // Calculate special theta at each step, S[i] already == dS * i
            }
            double sigmaSq = Vol * Vol * i * i;
            double rPart = (Int_Rate - D) * i;

            a[i] = -0.5 * theta * dt * (sigmaSq - rPart);
            b[i] = 1 + theta * dt * (sigmaSq + Int_Rate);
            c[i] = -0.5 * theta * dt * (sigmaSq + rPart); // A weighting of theta is applied to the implicit scheme

            // Right-hand side values using Crank-Nicolson averaging
            d[i] = (1 - (1 - theta) * dt * (sigmaSq + Int_Rate)) * V[i][k + 1]
                   + 0.5 * (1 - theta) * dt * (sigmaSq - rPart) * V[i - 1][k + 1]
                   + 0.5 * (1 - theta) * dt * (sigmaSq + rPart) * V[i + 1][k + 1]; // A weighting of (1 - theta) is applied to the explicit component
        }

        // Apply Gauss-Seidel Iteration
        double error = 1.0;
        int iter = 0;

        while (error > tol && iter < maxIter) {
            error = 0.0;

            for (int i = 1; i < NAS; ++i) {
                double oldV = V[i][k];

                V[i][k] = (d[i] - a[i] * V[i - 1][k] - c[i] * V[i + 1][k]) / b[i]; // Here the implicit weighting is taken from the LHS and applied to the RHS to give the final result 

                error = max(error, fabs(V[i][k] - oldV)); // Compute max error for convergence check
            }

            iter++;
        }

        // Early exercise condition for American options
        if (EType == 'Y') {
            for (int i = 1; i < NAS; ++i) {
                V[i][k] = max(V[i][k], V[i][NTS]); // V[i][NTS] is set equal to the payoff funtion above when terminal conditions are set
            }

        }

    }

    return V;
}




double interpolate(double float_idx, int idx, const vector<vector<double>>& values) {
    return (1 - float_idx + idx)  * values[idx].front() + (float_idx - idx) * values[idx + 1].front(); // Interpolate to increase accuracey
}


int main() {

    double S0, E, T, sigma, r, D, theta;
    int NAS, NTS;
    char EType, thetaType;
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
    cout << "Enter Dividend Yeild (D as decimal) :";
    cin >> D;
    cout << "Enter Number of Asset Steps (M) :";
    cin >> NAS;
    cout << "Enter Number of Time Steps (N) :";
    cin >> NTS;
    cout << "Enter If Option is American (Y/N) :";
    cin >> EType;
    cout << "Would you like to use theta special (1/2 + 1/12 * dS * dS * dt) (Y/N) :";
    cin >> thetaType;
    if (thetaType != 'Y'){
        cout << "Enter the value of Theta for Crank-Nicholson (theta as decimal):";
        cin >> theta;
    } else {
        theta == 0;
    }
    /*
    //Test conditions
    double S0 = 100;
    double E = 100;
    double T = 1;
    double sigma = 0.2;
    double r = 0.05;
    double D = 0;
    double theta = 1; 
    int NAS = 100;
    int NTS = 1000;
    char thetaType = 'N';
    char EType = 'N';
    */
    double S_max = 2 * S0;

    double float_idx = S0 / (S_max / NAS);

    int idx = static_cast<int>(float_idx);

    cout << "\n-------------------------------" << endl;
    if (EType == 'Y') {
        cout << "    American Option Prices" << endl;
    } else {
        cout << "    European Option Prices" << endl;
    }
    cout << "-------------------------------" << endl;

    vector<vector<double>> put_values_explicit = explicitScheme(sigma, r, 'P', E, T, S_max, NAS, NTS, D, EType);
    vector<vector<double>> call_values_explicit = explicitScheme(sigma, r, 'C', E, T, S_max, NAS, NTS, D, EType);

    cout << "\nPut Option Price (Explicit Scheme): " << interpolate(float_idx, idx, put_values_explicit) << endl;
    cout << "Call Option Price (Explicit Scheme): " << interpolate(float_idx, idx, call_values_explicit) << endl;

    
    vector<vector<double>> put_values_implicit = implicitScheme(sigma, r, 'P', E, T, S_max, NAS, NTS, D, EType);
    vector<vector<double>> call_values_implicit = implicitScheme(sigma, r, 'C', E, T, S_max, NAS, NTS, D, EType);

    cout << "\nPut Option Price (Implicit Scheme): " << interpolate(float_idx, idx, put_values_implicit) << endl;
    cout << "Call Option Price (Implicit Scheme): " << interpolate(float_idx, idx, call_values_implicit) << endl;


    vector<vector<double>> call_values_crank = crankNicolson(sigma, r, 'C', E, T, S_max, NAS, NTS, D, EType);
    vector<vector<double>> put_values_crank = crankNicolson(sigma, r, 'P', E, T, S_max, NAS, NTS, D, EType);

    cout << "\nPut Option Price (Crank Nicolson Scheme): " << interpolate(float_idx, idx, put_values_crank) << endl;
    cout << "Call Option Price (Crank Nicolson Scheme): " << interpolate(float_idx, idx, call_values_crank) << endl;


    vector<vector<double>> call_values_crank_theta = crankNicolsonTheta(sigma, r, 'C', E, T, S_max, NAS, NTS, D, EType, theta, thetaType);
    vector<vector<double>> put_values_crank_theta = crankNicolsonTheta(sigma, r, 'P', E, T, S_max, NAS, NTS, D, EType, theta, thetaType);

    cout << "\nPut Option Price (Theta Scheme): " << interpolate(float_idx, idx, put_values_crank_theta) << endl;
    cout << "Call Option Price (Theta Scheme): " << interpolate(float_idx, idx, call_values_crank_theta) << endl;
    

    system("pause");
    return 0;
}
/*
The impact of dividends on American options:
1. American Call Options:
    If a stock pays dividends, it is likely that the call option will be exercised early before the dividend is paid.
    This is because after a dividend is paid the stock price often drops making the call option less valuable.
    As such early exercise for in the money calls is often optimal before the ex-dividend date.
    
2. American Put Options:
    If a stock pays dividends, it is likely that the put option will be exercised after the dividend is paid.
    This is because after a dividend is paid the stock price often drops making the put option more valuable.
    Therefore, exercising before the ex-dividend date is often suboptimal even if the put is heavily in the money, as it will likely be worth more after the dividend has been paid.

When seen in my codes results a higher dividend rate leads to put options becoming more expensive and call options becoming less expensive shrinking the gap in their price difference.

Theta Special:
    When theta is 0, we have the explicit scheme, when theta is 1, we have the implicit scheme and with theta = 0.5 we have the Crank Nicholson Scheme

    Thus, by setting theta to theta special the error is improved.

    Theta special modifies the weighting dynamically which improves stability and accuracy.

*/