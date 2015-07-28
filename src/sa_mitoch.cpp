// individual-based simulation with multiple mitochondria 


#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "bramauxiliary.h"

using namespace std;

// parameters
const int numgen = 5000; // number of generations 
const int N = 10000; // population size
const int Nmito = 200; // number of mitochondria per individual
const int Nmito_haplo = 100; // number of mitochondria per individual
const int clutch = 5;  // number of offspring each individual produces
int Nbottleneck = 50; // size of the mitochondrial bottleneck
double mu = 0; // mutation probability
double seed = 0; // random seed (initialized in Init)
double sf = 0; // female selection coefficient against Mm
double sm = 0; // male selection coefficient against Mf
double k = 0; // parameter describing dominance / recessivity 
bool biparental = 0; // biparental inheritance or not
double leakage = 0; // proportion of alleles that are leaked

// counters
int generation = 0; // generation number
int NOva = 0; // number of ova in generation t
int NSperm = 0; // number of sperm in generation t
int Nm = 0; // number of males in generation t
int Nf = 0; // number of females in generation t
int Nmsurv = 0; // number of surviving males in generation t
int Nfsurv = 0; // number of surviving females in generation t
double total_time = 0; // keep track of the total time the simulation lasts

// others
gsl_rng_type const * T; // gnu scientific library rng type
gsl_rng *r; // gnu scientific rng 

// the individual
struct Individual
{
    bool mito[Nmito]; // array containing all mitochondria
    double fraction_Mm; // proportion Mm alleles of total
    bool sex; // male or female?

    int id; // for debugging purposes
};

struct Gamete
{
    // each gamete receives half of the number of 
    // mitochondria
    bool haplomito[Nmito_haplo];
};

// make arrays for population and the next generation
typedef Individual Population[N];
typedef Gamete GametePool[N*clutch];
Population Males, Females;
GametePool Ova,Sperm;


// initialize files to write data 
string filename("sim_sa_mito");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  // output file 

// initialize some parameters from command line arguments
void init_arguments(int argc, char *argv[])
{
    mu = atof(argv[1]);
    Nbottleneck = atof(argv[2]);
    sf = atof(argv[3]);
    sm = atof(argv[4]);
    k = atof(argv[5]);
    leakage = atof(argv[6]);
}

// write output data
void write_data()
{
    double mean_prop_mm = 0;
    double ss_prop_mm = 0;

    for (int i = 0; i < Nm; ++i)
    {
        mean_prop_mm += Males[i].fraction_Mm;
        ss_prop_mm +=Males[i].fraction_Mm * Males[i].fraction_Mm;
    }

    for (int i = 0; i < Nf; ++i)
    {
        mean_prop_mm += Females[i].fraction_Mm;
        ss_prop_mm += Females[i].fraction_Mm * Females[i].fraction_Mm;
    }

    mean_prop_mm /= Nf + Nm;

    double var_prop_mm = ss_prop_mm / (Nf + Nm) - mean_prop_mm * mean_prop_mm;

    DataFile << generation << ";" <<
                mean_prop_mm << ";" <<
                var_prop_mm << ";" <<
                Nfsurv << ";" <<
                Nmsurv << ";" << endl;
}

// parameters at the end of the sim
void write_parameters()
{
    total_time = time(NULL) - total_time;
    DataFile << endl << endl << 
                    "inheritance;" << (biparental ? "biparental" : "uniparental") << endl
                    << "k;" << k << endl
                    << "sf;" << sf << endl
                    << "sm;" << sm << endl
                    << "mu;" << mu << endl
                    << "leakage;" << leakage << endl
                    << "seed;" << seed  << endl
                    << "N;" << N << endl
                    << "Nmito;" << Nmito << endl
                    << "clutch;" << clutch << endl
                    << "Nbottleneck;" << Nbottleneck << endl
                    << "total_time;" << total_time << endl;
}

// initialize the population to particular values
// at the start of the simulation
void init_pop()
{ 
    // start the time
    total_time = time(NULL);

    // obtain a seed from current nanosecond count
	seed = get_nanoseconds();

    // set the seed to the random number generator
    // stupidly enough, this can only be done by setting
    // a shell environment parameter
    stringstream s;
    s << "GSL_RNG_SEED=" << setprecision(10) << seed;
    putenv(const_cast<char *>(s.str().c_str()));

    // set up the random number generators
    // (from the gnu gsl library)
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    // initialize all individuals
    // with a random sample of Mm and Mf
    // mitochondria
    for (int i = 0; i < N/2; ++i)
    {
        Males[i].fraction_Mm = 0;
        for (int j = 0; j < Nmito; ++j)
        {
            Males[i].mito[j] = gsl_ran_bernoulli(r, 0.5);
            Males[i].fraction_Mm += Males[i].mito[j];
            Males[i].id = 5;
        }
        Males[i].fraction_Mm /= Nmito;
    }
    
    for (int i = 0; i < N/2; ++i)
    {
        Females[i].fraction_Mm = 0;
        for (int j = 0; j < Nmito; ++j)
        {
            Females[i].mito[j] = gsl_ran_bernoulli(r, 0.5);
            Females[i].fraction_Mm += Females[i].mito[j];
            Females[i].id = 6;
        }
        
        Females[i].fraction_Mm /= Nmito;
    }

    Nf = N/2;
    Nm = N/2;
    NOva = 0;
    NSperm = 0;
}


// Mutate individual mitochondria
void Mutate(Individual &ind)
{
    for (int i = 0; i < Nmito; ++i)
    {
        if (gsl_rng_uniform(r) < mu)
        {
            ind.mito[i] = !ind.mito[i];
        }
    }

}

// create four gametes
void create_gametes(Individual &Parent, bool sex)
{
    // implement bottleneck as follows:
    // 1. random shuffle array with mitochondria
    // 2. pick the first nb mitochondria

    bool mito_bottlenecked[Nbottleneck]; 

    // sample without replacement Nbottleneck mitochondria
    // from the set of mitochonri
    gsl_ran_choose(r, 
                    mito_bottlenecked, 
                    Nbottleneck, 
                    Parent.mito, 
                    Nmito, 
                    sizeof(bool));

    // then sample from the bottlenecked collection
    // to get the original number of Nmito mitochondria
    gsl_ran_sample(r, 
            Parent.mito,
            Nmito,
            mito_bottlenecked,
            Nbottleneck,
            sizeof(bool));

    // allocate space for a duplicated set of mitochondria
    bool mito_dup[Nmito*2];

    // duplicate each mitochondrion
    for (int i = 0; i < Nmito;++i)
    {
        mito_dup[i*2] = Parent.mito[i];
        mito_dup[i*2+1] = Parent.mito[i];
    }

    // run meiosis 1,2
    for (int i = 0; i < rint(0.25 * clutch)/2; ++i)
    {
        // initialize two daughtercells
        Individual daughtercell1;
        Individual daughtercell2;

        // re-shuffle array of all mitochondria that have been
        // previously duplicated in parental (mast)cell
        gsl_ran_shuffle(r, mito_dup, Nmito*2, sizeof(bool));

        // transmit the mitochondria to the daughtercell
        // different ways of doing this, all are actually O(Nmito)
        for (int j = 0; j < Nmito; ++j)
        {
            // even mitochondria go to one daughtercell
            daughtercell1.mito[j] = mito_dup[j*2];
            // the odd ones to the other
            daughtercell2.mito[j] = mito_dup[j*2+1];
        }

        // shuffle mitochondrial ordering in daughtercells to guarantee
        // that mitochondrial transmission to the gametes
        // is random
        gsl_ran_shuffle(r, daughtercell1.mito, Nmito, sizeof(bool));
        gsl_ran_shuffle(r, daughtercell2.mito, Nmito, sizeof(bool));

        // generate the four gametes
        // and add them to the gamete pool
        Gamete gamete1; 
        Gamete gamete2; 
        Gamete gamete3; 
        Gamete gamete4; 

        for (int j = 0; j < rint(0.5 * Nmito); ++j)
        {
            gamete1.haplomito[j] = daughtercell1.mito[j*2];
            gamete2.haplomito[j] = daughtercell1.mito[j*2 + 1];
            gamete3.haplomito[j] = daughtercell2.mito[j*2];
            gamete4.haplomito[j] = daughtercell2.mito[j*2+1];
        }

        if (sex)
        {
            // add gametes to the gamete pool
            Ova[NOva++] = gamete1;
            Ova[NOva++] = gamete2;
            Ova[NOva++] = gamete3;
            Ova[NOva++] = gamete4;
            
            assert(NOva <= N * clutch);
        }
        else
        {
            Sperm[NSperm++] = gamete1;
            Sperm[NSperm++] = gamete2;
            Sperm[NSperm++] = gamete3;
            Sperm[NSperm++] = gamete4;

            assert(NSperm <= N * clutch);
        }
    }
}

// Survival selection 
void survive()
{
    // let females survive
    for (int i = 0; i < Nf; ++i)
    {
        // female dies
        // remove from pool
        if (gsl_rng_uniform(r) > 1 - pow(Females[i].fraction_Mm,k)*sf)
        {
            Females[i] = Females[Nf-1];
            --Nf;
            --i;
        }
    }

    Nfsurv = Nf;
    
    // let males survive
    for (int i = 0; i < Nm; ++i)
    {
        // male dies
        // remove from pool
        if (gsl_rng_uniform(r) > 1 - (1-pow(Males[i].fraction_Mm,k))*sm)
        {
            Males[i] = Males[Nm-1];
            --Nm;
            --i;
        }
    }
    
    Nmsurv = Nm;

    if (Nm == 0 || Nf == 0)
    {
        write_data();
        write_parameters();
        exit(1);
    }
}

// generate the new generation and replace the old
void replace_adults()
{
    NOva = 0;
    NSperm = 0;

    // make ova
    for (int i = 0; i < Nf; ++i)
    {
        create_gametes(Females[i],1);
    }

    assert(NOva < N * clutch);

    // make sperm
    if (biparental || leakage > 0)
    {
        for (int i = 0; i < Nm; ++i)
        {
            create_gametes(Males[i],0);
        }
    }

    Nm = 0;
    Nf = 0;

    // create the next generation
    for (int i = 0; i < N; ++i)
    {
        Individual kid;
        kid.id = 7;
        // randomly choose egg (eggs can, in principle, be used multiple times)
        Gamete ova = Ova[gsl_rng_uniform_int(r, NOva)];

        Gamete sperm;

        // randomly choose sperm
        if (biparental || leakage > 0)
        {
            // sample random spemr
            sperm = Sperm[gsl_rng_uniform_int(r, NSperm)];
        }

        kid.fraction_Mm = 0;
        // transfer mitochondria to kid
        for (int j = 0; j < Nmito/2; ++j)
        {
            kid.mito[j*2] = ova.haplomito[j];

            if (biparental)
            {
                kid.mito[j*2+1] = sperm.haplomito[j];
            }
            else if (leakage > 0)
            {
                if (j < rint(leakage * Nmito))
                {
                    kid.mito[j*2+1] = sperm.haplomito[gsl_rng_uniform_int(r, Nmito/2)];
                }
                else
                {
                    kid.mito[j*2+1] = ova.haplomito[gsl_rng_uniform_int(r, Nmito/2)];
                }
            }
            else
            {
                kid.mito[j*2+1] = ova.haplomito[gsl_rng_uniform_int(r, Nmito/2)];
            }

            kid.fraction_Mm += kid.mito[j*2];
            kid.fraction_Mm += kid.mito[j*2+1];
        }

        kid.fraction_Mm /= Nmito;

        // determine individual's sex
        // assuming a 1:1 sex ratio
        if (gsl_rng_uniform(r) < 0.5)
        {
            Mutate(kid);
            Males[Nm++] = kid;
        }
        else
        {
            Mutate(kid);
            Females[Nf++] = kid;
        }
    }
}

// write the headers of the datafile
void write_data_headers()
{
    DataFile << "generation;prop_mm;var_prop_mm;Nfsurv;Nmsurv;" << endl;
}

// write statistics to file


int main(int argc, char **argv)
{
    // initialize all the parameters
    init_arguments(argc, argv);

    // initialize all members of the population
    init_pop();

    // write headers of the data file
    write_data_headers();

    // skip indicates at which generation 
    // interval data is written to file
    int skip = 1;

    // loop through all the generations
    for (generation = 0; generation < numgen; ++generation)
    {
        survive();

        replace_adults();

        if (generation % skip == 0)
        {
            write_data();
        }
    }

    write_data();
    write_parameters();
}
