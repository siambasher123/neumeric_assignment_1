
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

const int MAX = 100;
const double TOL = 1e-6;

void jacobi_iterative_method(const vector<vector<double>>& a, const vector<double>& b, vector<double>& c, int n)
{
    vector<double> arr1(n, 0);
    int iteration = 0;
    bool converged;

    do
    {
        converged = true;
        for (int i = 0; i < n; ++i)
        {
            double sum = 0.0;
            for (int j = 0; j < n; ++j)
            {
                if (i != j)
                {
                    sum += a[i][j] * arr1[j];
                }
            }
            c[i] = (b[i] - sum) / a[i][i];

            if (fabs(c[i] - arr1[i]) > TOL)
                converged = false;
        }
        arr1 = c;
        iteration++;
    }
    while (!converged && iteration < MAX);
}


void gauss_seidel_iterative_method(const vector<vector<double>>& a, const vector<double>& b, vector<double>& c, int n)
{
    int iteration = 0;
    bool converged;

    do
    {
        converged = true;
        for (int i = 0; i < n; ++i)
        {
            double total = 0.0;
            for (int j = 0; j < n; ++j)
            {
                if (i != j)

                    total += a[i][j] * c[j];

            }
            double count1 = (b[i] - total) / a[i][i];

            if (fabs(count1 - c[i]) > TOL)
                converged = false;

            c[i] = count1;
        }
        iteration++;
    }
    while (!converged && iteration < MAX);
}


void gauss_elimination(vector<vector<double>> a, vector<double> b, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int k = i + 1; k < n; k++)
        {
            double factor = a[k][i] / a[i][i];
            for (int j = 0; j < n; j++)
            {
                a[k][j] -= factor * a[i][j];
            }
            b[k] -= factor * b[i];
        }
    }
    vector<double> c(n);
    for (int i = n - 1; i >= 0; i--)
    {
        c[i] = b[i];
        for (int j = i + 1; j < n; j++)
        {
            c[i] -= a[i][j] * c[j];
        }
        c[i] /= a[i][i];
    }
    cout << "Solution by Gauss Elimination: ";
    for (double x : c) cout << x << " ";
    cout << endl;
}
void gauss_jordan_elimination(vector<vector<double>> a, vector<double> b,int n)
{
    for (int i = 0; i < n; i++)
    {
        double pivoting = a[i][i];
        for (int j = 0; j < n; j++)
        {
            a[i][j] /= pivoting;
        }
       b[i] /= pivoting;
        for (int k = 0; k < n; k++)
        {
            if (k != i)
            {
                double f = a[k][i];
                for (int j = 0; j < n; j++)
                {
                    a[k][j] -= f *a[i][j];
                }
               b[k] -= f * b[i];
            }
        }
    }
    cout << "Solution by Gauss-Jordan Elimination: ";
    for (double x : b) cout << x << " ";
    cout << endl;
}

void lower_upper_factorization(const vector<vector<double>>& a, vector<double>  b, int n)
{
    vector<vector<double>> L(n, vector<double>(n, 0));
    vector<vector<double>> U(n, vector<double>(n, 0));
    vector<double> Y(n);
    vector<double> X(n);


    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            U[i][j] = a[i][j];
            for (int k = 0; k < i; k++)
            {
                U[i][j] -= L[i][k] * U[k][j];
            }
        }

        for (int j = i + 1; j < n; j++)
        {
            L[j][i] = a[j][i];
            for (int k = 0; k < i; k++)
            {
                L[j][i] -= L[j][k] * U[k][i];
            }
            L[j][i] /= U[i][i];
        }
        L[i][i] = 1.0;
    }


    cout << "Lower Matrix L:\n";
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << setw(10) << L[i][j] << " ";
        }
        cout << endl;
    }


    cout << "\nUpper Matrix U:\n";
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << setw(10) << U[i][j] << " ";
        }
        cout << endl;
    }

    // Forward substitution for L * Y = B
    for (int i = 0; i < n; i++)
    {
        Y[i] = b[i];
        for (int j = 0; j < i; j++)
        {
            Y[i] -= L[i][j] * Y[j];
        }
    }


    cout << "\nY matrix:\n";
    for (int i = 0; i < n; i++)
    {
        cout << Y[i] << " ";
    }
    cout << endl;

    // Back substitution for U * X = Y
    for (int i = n - 1; i >= 0; i--)
    {
        X[i] = Y[i];
        for (int j = i + 1; j < n; j++)
        {
            X[i] -= U[i][j] * X[j];
        }
        X[i] /= U[i][i];
    }


    cout << "\nSolution Matrix X:\n";
    for (int i = 0; i < n; i++)
    {
        cout << X[i] << " ";
    }
    cout << endl;
}

int main()
{
    int n, z;
    cout << "Enter the number of variables: ";
    cin >> n;

    vector<vector<double>> a(n, vector<double>(n));
    vector<double> b(n);
    vector<double> c(n, 0);

    cout << "Enter the coefficients of the matrix A:\n";
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            cin >> a[i][j];

    cout << "Enter the constants vector B:\n";
    for (int i = 0; i < n; i++)
        cin >> b[i];

    cout << "Choose a method to solve the  linear system:\n";
    cout << "1. Jacobi Iterative Method \n";
    cout << "2. Gauss-Seidel Iterative Method \n";
    cout << "3. Gauss Elimination\n";
    cout << "4. Gauss-Jordan Elimination\n";
    cout << "5. LU Factorization\n";
    cout << "Enter your choice: ";
    cin >> z;

    switch (z)
    {
    case 1:
        jacobi_iterative_method(a,b,c, n);
        cout << "Solution by Jacobi Method: ";
        for (double x : c) cout << x << " ";
        cout << endl;
        break;
    case 2:
         gauss_seidel_iterative_method (a,b,c, n);
        {

            cout << "Solution by Gauss-Seidel Method: ";
            for (double x : c) cout << x << " ";
            cout << endl;
            break;
        case 3:
            gauss_elimination(a,b, n);
            break;
        case 4:
            gauss_jordan_elimination(a,b, n);
            break;
        case 5:
            lower_upper_factorization(a,b, n);
            break;
        default:
            cout << "Invalid choice. Please select a valid method." << endl;
        }

        return 0;
    }
}
