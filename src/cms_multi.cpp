#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include "bramauxiliary.h"


using namespace std;


// output file 
string filename("iter_cms_multi");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  

//parameters and arrays to count the numbers
// (and hence frequencies of all genotypes
const int Nmito = 100;
const int Ngeno = 3;
double pi = 0;
double pvecf[Nmito+1][3];
double pvecm[Nmito+1][3];
double pvecfold[Nmito+1][3];
double pvecmold[Nmito+1][3];

double mu_auto[3][3];

// parameters
double mu = 0;
double nu = 0;
int Nbottleneck = 1;
double c = 0;
double Ws = 0;
double alpha_ws = 0;
double alpha_cms = 0;
double fa = 0;
int skip = 1;
int init = 0;

// initialize some parameters from command line arguments
void init_arguments(int argc, char *argv[])
{
    mu = atof(argv[1]);
    nu = atof(argv[2]);
    Nbottleneck = atoi(argv[3]);
    c = atof(argv[4]);
    Ws = atof(argv[5]);
    alpha_ws = atof(argv[6]);
    alpha_cms = atof(argv[7]);
    fa = atof(argv[8]);

    if (Nbottleneck > Nmito)
    {
        Nbottleneck = Nmito;
    }


    mu_auto[0][0] = (1-nu) * (1-nu);
    mu_auto[0][1] = 2 * (1-nu) * nu;
    mu_auto[0][2] = nu * nu;
    mu_auto[1][0] = 2 * (1-nu) * nu;
    mu_auto[1][1] = (1-nu)*(1-nu) + nu * nu;
    mu_auto[1][2] = 2 * (1-nu) * nu;
    mu_auto[2][0] = nu * nu;
    mu_auto[2][1] = 2 * (1-nu) * nu;
    mu_auto[2][2] = (1-nu) * (1-nu);

}

void mutation()
{
    double pvecftplus1[Nmito+1][Ngeno];
    double pvecmtplus1[Nmito+1][Ngeno];


    // loop over the rows of the matrix M
    for (int j = 0; j <= Nmito; ++j)
    {
        for (int k = 0; k < Ngeno; ++k)
        {
            // set t+1 to 0
            pvecftplus1[j][k] = 0;
            pvecmtplus1[j][k] = 0;

            // mutation probability of individuals
            // containing l cms mitochondria
            for (int l = 0; l <= Nmito; ++l)
            {
                double mu_lj = 0;

                int max_x = Nmito - l < j ? Nmito - l : j;

                // x mitochondria mutate from 0 to 1
                //
                // x >= j - l; i.e., the number of mutations from 0 to 1
                // should be equal to or larger than the difference between the number of 1s
                // between t and t+1
                //
                // similarly, x < Nmito - l, j  x should be smaller than the the number of 0s 
                // initially present and smaller than the number of final number of 1s
                //
                // loop through the total number of mitochondria in state 0
                // before mutation occurs
                for (int x = (j - l > 0 ? j - l : 0); x <= max_x; ++x)
                {
                    // we start with j type 1 mitochondria
                    // from the Nmito - j type 0 mitochondria

                    mu_lj += gsl_ran_binomial_pdf(x, mu, Nmito - l) * gsl_ran_binomial_pdf(l - j + x, mu, l);
                }

                for (int h = 0; h < Ngeno; ++h)
                {
                    pvecftplus1[j][k] += pvecf[l][h] * mu_lj * mu_auto[h][k];
                    pvecmtplus1[j][k] += pvecm[l][h] * mu_lj * mu_auto[h][k];
                }
            }
        }
    }

    double sumcheck = 0.0;
    double sumcheckm = 0.0;

    // update frequencies
    for (int i = 0; i <= Nmito; ++i)
    {
        for (int j = 0; j < Ngeno; ++j)
        {
            pvecf[i][j] = pvecftplus1[i][j];
            pvecm[i][j] = pvecmtplus1[i][j];
            sumcheck += pvecf[i][j];
            sumcheckm += pvecm[i][j];
        }
    }

    assert(sumcheck - 1.0 < 1e-10);
    assert(sumcheckm - 1.0 < 1e-10);

}

// bottlenecks
void bottleneck()
{
    double pvecftplus1[Nmito+1][Ngeno];
    double pvecmtplus1[Nmito+1][Ngeno];

    double B_inverse = 1.0 / Nbottleneck;

    if (Nbottleneck == Nmito)
    {
        return;
    }

    // loop through the distribution of i Cm mitochondria
    // after the bottleneck and resampling
    for (int i = 0; i <= Nmito; ++i)
    {
        for (int k = 0; k < Ngeno; ++k)
        {
            pvecftplus1[i][k] = 0;
            pvecmtplus1[i][k] = 0;

            // loop through the distribution of j Cm mitochondria
            // that are starting before the bottleneck
            for (int j = 0; j <= Nmito; ++j)
            {
                double b_ij = 0;

                // loop through the distribution of h1 Cm mitochondria at the end of the bottleneck
                for (int h1 = 0; h1 <= Nbottleneck; ++h1)
                {
                    double prob = h1 == 0 ? 0 : h1 * B_inverse;

                    b_ij += gsl_ran_hypergeometric_pdf(h1, j, Nmito - j, Nbottleneck)  * gsl_ran_binomial_pdf(i, prob, Nmito);
                }

                pvecftplus1[i][k] += b_ij * pvecf[j][k];
                pvecmtplus1[i][k] += b_ij * pvecm[j][k];
            }
        }
    }

    double sumcheck = 0.0;
    for (int i = 0; i <= Nmito; ++i)
    {
        for (int j = 0; j < Ngeno; ++j)
        {
            pvecf[i][j] = pvecftplus1[i][j];
            pvecm[i][j] = pvecmtplus1[i][j];
            sumcheck += pvecm[i][j];
        }
    }

    assert(sumcheck - 1.0 < 1e-10);
}

void meiosis_1_2()
{
    // allocate array for the distribution of mitochondria 
    // after this life cycle stage
    double pvecftplus1[Nmito+1][Ngeno];
    double pvecmtplus1[Nmito+1][Ngeno];

    for (int i = 0; i <= Nmito; ++i)
    {
        for (int j = 0; j < Ngeno; ++j)
        {
            pvecftplus1[i][j] = 0;
        }
    }

    for (int i = 0; i <= Nmito; ++i)
    {
        // first meiotic subdivision
        // implying duplication of the mother cell into two daughter cells
        for (int k = 0; k <= Nmito; ++k)
        {
            // duplication and splitting into daughter cells according to hypergeometric
            // distribution yields a probability that an individual with k mitochondria 
            // produces a daughter cell with i mitochondria
            // notice that in the description of Hadjivasiliou et al 2012
            // summation only starts for k = j/2 -> k = M. However, 
            double e1_ij = gsl_ran_hypergeometric_pdf(i, 2*k, 2*(Nmito-k), Nmito);

            // hence proportion of individuals bearing i mitochondria is cumulatively
            // multiplied by e1_ij
            pvecftplus1[i][0] += e1_ij * (pvecf[k][0] + 1.0 / 6 * pvecf[k][1]);
            pvecmtplus1[i][0] += e1_ij * (pvecm[k][0] + 1.0 / 6 * pvecm[k][1]);
            pvecftplus1[i][1] += 2.0 /3 * e1_ij * pvecf[k][1];
            pvecmtplus1[i][1] += 2.0 /3 * e1_ij * pvecm[k][1];
            pvecftplus1[i][2] += e1_ij * (pvecf[k][2] + 1.0 / 6 * pvecf[k][1]);
            pvecmtplus1[i][2] += e1_ij * (pvecm[k][2] + 1.0 / 6 * pvecm[k][1]);
        }
    }

    for (int i = 0; i <= Nmito; ++i)
    {
        for (int j = 0; j < Ngeno; ++j)
        {
            pvecf[i][j] = pvecftplus1[i][j];
            pvecm[i][j] = pvecmtplus1[i][j];
            pvecftplus1[i][j] = 0;
            pvecmtplus1[i][j] = 0;
        }
    }
    
    // second round of meiosis
    // split Nmito mitochondria in two for each granddaughter cell
    for (int i = 0; i <= Nmito/2; ++i)
    {
        for (int j = i; j <= Nmito; ++j)
        {
            double e2_ij = gsl_ran_hypergeometric_pdf(i, j, Nmito-j, Nmito/2);

            pvecftplus1[i][0] += e2_ij * pvecf[j][0] + .5 * e2_ij * pvecf[j][1];
            pvecmtplus1[i][0] += e2_ij * pvecm[j][0] + .5 * e2_ij * pvecm[j][1];
            pvecftplus1[i][1] += e2_ij * pvecf[j][2] + .5 * e2_ij * pvecf[j][1];
            pvecmtplus1[i][1] += e2_ij * pvecm[j][2] + .5 * e2_ij * pvecm[j][1];
        }
    }

    double wm, wf, wmbar, wfbar;

    wmbar = 0;
    wfbar = 0;


    // now calculate gametic fitness
    for (int i = 0; i <= Nmito; ++i)
    {
        for (int j = 0; j < Ngeno; ++j)
        {
            // male fitness in the presence of a recessive restorer
            wm = j >= 1 ? 1.0 * (1-c) : 1.0 - pow((double) i / Nmito,alpha_cms); 
            wf = 1.0; 
            
            // recessive restorer
            if (j < 1)
            {
                wf *= fa * (1.0 + pow((double) i / Nmito,alpha_cms));
            }

            if (i > 0)
            {
                if (Ws > 0)
                {
                    wf *= Ws;
                }
                else
                {
                   wf *= 1.0 - pow((double) i / Nmito, alpha_ws);
                }
            }

            if (i > 5)
            {
                cout << i << " " << wf << " " << wm << ";" << endl;
            }

            pvecf[i][j] = wf * pvecftplus1[i][j];
            pvecm[i][j] = wm * pvecmtplus1[i][j];
            
            wmbar += pvecm[i][j];
            wfbar += pvecf[i][j];

        }
    }

    // now normalize freqs (dividing by mean fitness)
    for (int i = 0; i <= Nmito; ++i)
    {
        for (int j = 0; j < Ngeno; ++j)
        {
            pvecf[i][j] /= wfbar;
            pvecm[i][j] /= wmbar;
        }
    }
}

// syngamy assuming uniparental inheritance
void syngamy()
{
    double pvectplus1[Nmito+1][Ngeno];

    for (int i = 0; i <= Nmito; ++i)
    {
        pvectplus1[i][0] = 0;
        pvectplus1[i][1] = 0;
        pvectplus1[i][2] = 0;

        // number of Cm mitochondria present in sperm
        for (int rm = 0; rm <= Nmito/2; ++rm)
        {
            // rf = the number of Cm mitochondria present in the egg
            for (int rf = 0; rf <= Nmito/2; ++rf)
            {
                pvectplus1[i][0] += pvecf[rf][0] * pvecm[rm][0] *
                   gsl_ran_binomial_pdf(i, 2.0 * rf / Nmito, Nmito);
                
                pvectplus1[i][1] += (pvecf[rf][0] * pvecm[rm][1] + pvecf[rf][1] * pvecm[rm][0]) *
                   gsl_ran_binomial_pdf(i, 2.0 * rf / Nmito, Nmito);
                
                pvectplus1[i][2] += (pvecf[rf][1] * pvecm[rm][1]) *
                   gsl_ran_binomial_pdf(i, 2.0 * rf / Nmito, Nmito);
            }
        }
        //cout << pvectplus1[i][0] << " " << pvectplus1[i][1] << " " << pvectplus1[i][2] << endl;
    }

    double sumcheck = 0.0;
    for (int i = 0; i <= Nmito; ++i)
    {
        for (int j = 0; j < Ngeno; ++j)
        {
            pvecf[i][j] = pvectplus1[i][j];
            pvecm[i][j] = pvectplus1[i][j];
            sumcheck += pvecm[i][j];
        }
    }

    assert(sumcheck - 1.0 < 1e-10);
}

void write_data(int const time)
{
    DataFile << time << ";";

    double sumpvec = 0;
    double cms_freq = 0;
    double rest_freq = 0;

    for (int i = 0; i <= Nmito; ++i)
    {
        for (int j = 0; j < Ngeno; ++j)
        {
            DataFile << pvecf[i][j] << ";" << pvecm[i][j] << ";";

            sumpvec += pvecf[i][j] + pvecm[i][j];

            cms_freq += pvecf[i][j] * ((double)i/Nmito);
            rest_freq += pvecf[i][j] * .5 * j;
        }
    }

    DataFile << sumpvec << ";" << cms_freq << ";" << rest_freq << ";" << endl;

}

void write_data_headers()
{
    DataFile << "time;";

    for (int i = 0; i <= Nmito; ++i)
    {
        for (int j = 0; j < Ngeno; ++j)
        {
            DataFile << "pf" << i  << j << ";pm" << i << j << ";";   
        }
    }

    DataFile << "sump;cms_freq;rest_freq;" << endl;
}

void write_parameters()
{
    DataFile << endl << endl <<
        "mu;" << mu << endl <<
        "Nbottleneck;" << Nbottleneck << endl <<
        "c;" << c << endl <<
        "Ws;" << Ws << endl <<
        "alpha_ws;" << alpha_ws << endl <<
        "alpha_cms;" << alpha_cms << endl <<
        "fa;" << fa << endl; 
}

bool check_bounds()
{
    bool check = true;

    for (int i = 0; i <= Nmito; ++i)
    {
        for (int j = 0; j < Ngeno; ++j)
        {
            if (fabs(pvecf[i][j] - pvecfold[i][j]) > 1e-9)
            {
                check = false;
                break;
            }

            if (fabs(pvecm[i][j] - pvecmold[i][j]) > 1e-9)
            {
                check = false;
                break;
            }
        }
    }
    return(check);
}


int main(int argc, char **argv)
{
    // initialize all the parameters
    init_arguments(argc, argv);

    // initialize vector pvec
    for (int i = 0; i <= Nmito; ++i)
    {
        for (int j = 0; j < Ngeno; ++j)
        {
            pvecm[i][j] = 0.0;
            pvecf[i][j] = 0.0;
        }
    }

    pvecm[0][0] = 0.9;
    pvecf[0][0] = 0.9;
    pvecm[1][0] = 0.05;
    pvecf[1][0] = 0.05;
    pvecm[0][1] = 0.05;
    pvecf[0][1] = 0.05;

    write_parameters();
    write_data_headers();

    for (int time = 0; time < 1e06; ++time)
    {
        double sumcheck = 0.0;
        for (int i = 0; i <= Nmito; ++i)
        {
            for (int j = 0; j < Ngeno; ++j)
            {
                pvecfold[i][j] = pvecf[i][j];
                pvecmold[i][j] = pvecm[i][j];
                sumcheck += pvecmold[i][j];
            }
        }
        cout << "sumcheck init: " << sumcheck << endl;

        mutation();
        // cout << "mut gehad" << endl;
            // write_data(time);
        bottleneck();
         // cout << "bottleneck gehad" << endl;
            //write_data(time);
        meiosis_1_2();
         // cout << "meiosis gehad" << endl;
          //  write_data(time);
        syngamy();
         // cout << "syngamy gehad" << endl;
           // write_data(time);

        if (check_bounds())
        {
            write_data(time);
            break;
        }

        if (time % skip == 0)
        {
            write_data(time);
        }
    }

}
