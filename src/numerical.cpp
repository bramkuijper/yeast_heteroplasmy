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

int const Nmito = 50;

double p[Nmito+1];

double w[Nmito+1];
double sumw = 0;
double ch = 0;

size_t generation = 0;


// initialize files to write data 
string filename("iter_yeast_drift");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  // output file 


void selection() {

    double sumw = 0;

    for (size_t i = 0; i <= Nmito; ++i)
    {
        w[i] = i < Nmito/2 ? 1.0 - ch * 2 * i / (Nmito) * 2 * i / (Nmito) :
            1.0 - ch * 2 * (Nmito - i) / (Nmito) * 2 * (Nmito - i) / (Nmito);

        p[i] *= w[i];

        sumw += p[i];
    }

    cout << sumw << endl;

    for (size_t i = 0; i <= Nmito; ++i)
    {
        p[i] /= sumw;
    }

}

void replicate() {

    double ptplus1[Nmito+1];
    
    for (size_t i = 0; i <= Nmito; ++i)
    {
        ptplus1[i] = 0;
    }

    double sum = 0;

    for (size_t i = 0; i <= Nmito; ++i)
    {
//        for (size_t k = 0; k <= 2*Nmito; ++k)
//        {
//            for (size_t j = 0; j <= Nmito; ++j)
//            {
//                ptplus1[i] += 
//                    gsl_ran_hypergeometric_pdf(i, k, 2*Nmito - k, Nmito) *
//                    gsl_ran_binomial_pdf(k, j / Nmito, 2 * Nmito) * p[j];
//
//            }
//        }
        for (size_t j = 0; j <= Nmito; ++j)
        {
            ptplus1[i] += 
                gsl_ran_hypergeometric_pdf(i, j, 2*Nmito - j, Nmito) * p[j];

        }

        sum += ptplus1[i];
    }

    cout << "s: " << sum << endl;


    for (size_t i = 0; i <= Nmito; ++i)
    {
        p[i] = ptplus1[i];
    }
}

void ascus()
{
    double ptplus1[4*Nmito+1];
    double ptplus2[2*Nmito+1];
    
    for (size_t i = 0; i <= 2*Nmito; ++i)
    {
        ptplus2[i] = 0;
    }

    for (size_t i = 0; i <= 4*Nmito; ++i)
    {
        // duplication to 4*M
        // when both parents contribute equally (i.e., when 
        // the number of parent-1 mitochondria are M/(2*M)
        //
        ptplus1[i] = gsl_ran_binomial_pdf(i, 0.5, 4 * Nmito);

    }

    // now make two daughter cells
    for (size_t i = 0; i <= 2*Nmito; ++i)
    {
        for (size_t j = 0; j <= 4*Nmito; ++j)
        {
            ptplus2[i] += gsl_ran_hypergeometric_pdf(i, j, 4 * Nmito - j, 2 * Nmito) * ptplus1[j];

        }
    }

    for (size_t i = 0; i <= Nmito;++i)
    {
        p[i] = 0;
    }


    // again divide the two daughter cells into two
    for (size_t i = 0; i <= Nmito; ++i)
    {
        for (size_t j = 0; j <= 2*Nmito; ++j)
        {
            p[i] += gsl_ran_hypergeometric_pdf(i, j, 2 * Nmito - j,  Nmito) * ptplus2[j];
        }
    }
}

void write_data_headers()
{
    DataFile << "generation;mean_p;";

    for (size_t i = 0; i <= Nmito; ++i)
    {
        DataFile << "p_" << i << ";";
    }

    DataFile << endl;
}

void write_data()
{
    double pfreq = 0;

    for (size_t i = 0; i <= Nmito; ++i)
    {
        pfreq += p[i] * ((double)i/Nmito);    
    }

    DataFile << generation << ";" << pfreq << ";";

    for (size_t i = 0; i <= Nmito; ++i)
    {
        DataFile << p[i] << ";";
    }
    DataFile << endl;
}

// initialize arguments
void init_arguments(int argc, char *argv[])
{
    ch = atof(argv[1]);
}

void write_parameters()
{
    DataFile << endl << endl
        << "ch;" << ch << endl;
}



int main(int argc, char **argv)
{
    // initialization
    init_arguments(argc, argv);

    write_data_headers();
    ascus();
    write_data();

    for (generation = 1; generation < 30; ++generation)
    {
        selection();
        replicate();
        write_data();
    }

    write_parameters();
}
