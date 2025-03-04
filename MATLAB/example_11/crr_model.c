#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <cmath>
#include <vector>

#define max(a,b) (a>b?a:b)
#define min(a,b) (a<b?a:b)

namespace py = pybind11;

// Function declarations
static void runge_kutta(int nvar, int nt, py::array_t<double> tout, py::array_t<double> y0, py::dict data, py::dict options, py::array_t<double> y);
static void rk2(py::dict data, int nvar, int s, double t, double h, std::vector<double>& u, std::vector<double>& LTE);
static int fRhs(int s, double t, std::vector<double>& u, std::vector<double>& udot, py::dict data);
static double expFlux(double Sr, double a);
static double exponen(double x);

void runge_kutta(int nvar, int nt, py::array_t<double> tout, py::array_t<double> y0, py::dict data, py::dict options, py::array_t<double> y)
{
    double hin, hmax_, hmin_, reltol, order;
    double h, t, t1, t2, wrms;
    int ns = nt - 1;
    
    // Create references to input and output data
    auto tout_buffer = tout.unchecked<1>();
    auto y0_buffer = y0.unchecked<1>();
    auto y_buffer = y.mutable_unchecked<2>();
    
    // Memory allocation for temporary vectors
    std::vector<double> LTE(nvar), ytmp(nvar), w(nvar);
    
    // Extract integration options
    hin = options["InitialStep"].cast<double>();
    hmax_ = options["MaxStep"].cast<double>();
    hmin_ = options["MinStep"].cast<double>();
    reltol = options["RelTol"].cast<double>();
    order = options["Order"].cast<double>();
    
    // Initialize y with y0
    for (int i = 0; i < nvar; i++) {
        y_buffer(i, 0) = y0_buffer(i);
    }

    // Integrate from tout[0] to tout[nt-1]
    for (int s = 1; s <= ns; s++) {
        t1 = tout_buffer(s - 1);
        t2 = tout_buffer(s);
        h = hin;
        h = max(hmin_, min(h, hmax_));
        h = min(h, t2 - t1);
        
        // Set initial y
        for (int i = 0; i < nvar; i++) {
            ytmp[i] = y_buffer(i, s - 1);
        }
        
        t = t1;
        while (t < t2) {
            // Advance solution by step h
            rk2(data, nvar, s, t, h, ytmp, LTE);
            
            // Check step acceptance
            wrms = 0;
            for (int i = 0; i < nvar; i++) {
                w[i] = 1.0 / (reltol * fabs(ytmp[i]) + options["AbsTol"].cast<py::array_t<double>>()(i));
                wrms += pow(w[i] * LTE[i], 2);
            }
            wrms = sqrt(wrms / nvar);
            
            if (wrms <= 1) {
                for (int i = 0; i < nvar; i++) {
                    ytmp[i] = ytmp[i];
                }
                t += h;
            }
            h = h * max(0.2, min(5.0, 0.9 * pow(wrms, -1.0 / order)));
            h = max(hmin_, min(h, hmax_));
            h = min(h, t2 - t);
        }
    }
}

void rk2(py::dict data, int nvar, int s, double t, double h, std::vector<double>& u, std::vector<double>& LTE)
{
    std::vector<double> udotE(nvar), uE(nvar), udot(nvar);
    
    // Euler method
    fRhs(s, t, u, udotE, data);
    for (int i = 0; i < nvar; i++) {
        uE[i] = u[i] + h * udotE[i];
    }
    
    // Heun method
    fRhs(s, t + h, uE, udot, data);
    for (int i = 0; i < nvar; i++) {
        u[i] = u[i] + 0.5 * h * (udotE[i] + udot[i]);
    }
    
    // Compute LTE
    for (int i = 0; i < nvar; i++) {
        LTE[i] = fabs(uE[i] - u[i]);
    }
}

int fRhs(int s, double t, std::vector<double>& u, std::vector<double>& udot, py::dict data)
{
    // Define parameters from data
    auto P = data["P"].cast<py::array_t<double>>().unchecked<1>();
    auto Ep = data["Ep"].cast<py::array_t<double>>().unchecked<1>();
    double Imax = data["Imax"].cast<double>();
    double Sumax = data["Sumax"].cast<double>();
    double Qsmax = data["Qsmax"].cast<double>();
    double aE = data["aE"].cast<double>();
    double aF = data["aF"].cast<double>();
    double aS = data["aS"].cast<double>();
    double Kf = data["Kf"].cast<double>();
    double Ks = data["Ks"].cast<double>();
    
    // Model calculations
    double Si = u[0], Su = u[1], Sf = u[2], Ss = u[3];
    double Precip = P[s - 1];
    double EvapI, P_e, Ep_e;

    if (Imax > 0.0) {
        EvapI = Ep[s - 1] * expFlux(Si / Imax, 50.0);
        P_e = P[s - 1] * expFlux(Si / Imax, -50.0);
        Ep_e = std::max(0.0, Ep[s - 1] - EvapI);
    } else {
        EvapI = 0.0;
        P_e = P[s - 1];
        Ep_e = Ep[s - 1];
    }

    double Evap = Ep_e * expFlux(Su / Sumax, aE);
    double Perc = Qsmax * expFlux(Su / Sumax, aS);
    double Runoff = P_e * expFlux(Su / Sumax, aF);
    double FastQ = Sf / Kf;
    double SlowQ = Ss / Ks;

    udot[0] = Precip - EvapI - P_e;
    udot[1] = P_e - Evap - Perc - Runoff;
    udot[2] = Runoff - FastQ;
    udot[3] = Perc - SlowQ;

    return 0;
}

double expFlux(double Sr, double a)
{
    Sr = std::max(0.0, std::min(1.0, Sr));
    if (fabs(a) < 1e-6) {
        return Sr;
    } else {
        return (1.0 - exponen(-a * Sr)) / (1.0 - exponen(-a));
    }
}

double exponen(double x)
{
    return exp(std::min(300.0, x));
}

PYBIND11_MODULE(crr_model, m) {
    m.def("runge_kutta", &runge_kutta, "Runge-Kutta method for rainfall-runoff model integration");
}
