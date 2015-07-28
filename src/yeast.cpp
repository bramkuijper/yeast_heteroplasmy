#include <ctime>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include "bramauxiliary.h"
#include "individual.hpp"

using namespace std;
const int numgen = 20;

const int Nreplicate = 10000;

size_t nMito_error = 5;
size_t nMito_min = 20;


vector <Individual> Pop;

double ch = 0;

int generation = 0;
int total_time = 0;
int seed = 0;
int experiment_i = 0;
int replicate = 0;

size_t type = 0;

// initialize random number generators
gsl_rng_type const * T; // gnu scientific library rng type
gsl_rng *r; // gnu scientific rng 

// initialize files to write data 
string filename("sim_yeast");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  // output file 

// comparison function as we need to do some array sorting
static int compare(const void *a, const void *b) {
    return *(int*)a - *(int*)b;
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
}

// make mitochondria both parents
// duplicate them
// distribute them randomly over the 4 cells of the ascus
// later on we add more elaborate mechanisms of ascus production
// where the mixture of mitochondria is not random
void make_tetrad(Individual *ascus, size_t ascus_size)
{
    // allocate array of mitochondria from two parents
    // that is going to be replicated in order to make spores
    bool mito_diploid[4*(nMito_min + nMito_error)];

    // number of mitochondria of parent 1
    size_t nMito_parent_1 = nMito_min 
        + (gsl_rng_uniform(r) < 0.5 ? -1.0 : 1.0) * gsl_rng_uniform_int(r, nMito_error);

    // number of mitochondria of parent 2
    size_t nMito_parent_2 = nMito_min
        + (gsl_rng_uniform(r) < 0.5 ? -1.0 : 1.0) * gsl_rng_uniform_int(r, nMito_error);
 

    size_t nMito_total = nMito_parent_1 + nMito_parent_2; 

    assert(nMito_parent_1 >= 0 && nMito_total < 2 * (nMito_min + nMito_error));
  
    // generate mitochondria parent 1
    for (size_t i = 0; i < nMito_parent_1; ++i)
    {
        mito_diploid[i] = 0;
    }

    // generate mitochondria parent 2
    for (size_t i = nMito_parent_1; i < nMito_total; ++i)
    {
        mito_diploid[i] = 1;
    }

    // duplicate mitochondria by random sampling
    // from the existing mitochondria and duplicating those
    //
    // right before the ascus is going to be produced
    for (size_t i = nMito_total; i < nMito_total * 2; ++i)
    {
        assert(i < 4 * (nMito_min + nMito_error));
        mito_diploid[i] = mito_diploid[gsl_rng_uniform_int(r, nMito_total)];
    }

    nMito_total = nMito_total*2;

    // draw three random numbers, the relative values of those
    // determine the number of mitochondria going to each spore
    // see http://stackoverflow.com/questions/8064629/random-numbers-that-add-to-100-matlab/8068956#8068956 
    int nmito_tetrad[4];

    for (size_t i = 0; i < 3; ++i)
    {
        nmito_tetrad[i] = gsl_rng_uniform_int(r, nMito_total);

    }

    // the last number is the total number of mitochondria that are going to be distributed
    nmito_tetrad[3] = nMito_total;

    // sort the numbers in accumulating order
    qsort(nmito_tetrad, 4, sizeof(int), compare);

    // var to remember the last number of mitochondria to be transmitted
    size_t lastnumber = 0;

    // remember the randomly chosen mitochondrion to go to a spore
    // so that we can delete it from the parent
    size_t random_mitochondrion;

    // now create ascus through random sampling
    for (size_t i = 0; i < ascus_size; ++i)
    {
        // calculate number of mitochondria going
        // to this tetrad
        size_t nmito_to_this_tetrad = nmito_tetrad[i] - lastnumber;
        lastnumber = nmito_tetrad[i];

        // set the begin count of mitochondria to 0
        // for each spore 
        ascus[i].nMito = 0;
        ascus[i].nMitoMut = 0;

        for (size_t j = 0; j < nmito_to_this_tetrad; ++j)
        {
            assert(nMito_total > 0);

            // sample a random mitochondrion
            random_mitochondrion = gsl_rng_uniform_int(r, nMito_total);

            // move mitochondrion to the spore

            ascus[i].nMitoMut += mito_diploid[random_mitochondrion];

            
            ++ascus[i].nMito;

            // remove mitochondrion from the diploid by replace it
            // with the last one of the stack and reducing the count
            mito_diploid[random_mitochondrion] = mito_diploid[--nMito_total];

        }
            
        cout << ascus[i].nMitoMut << endl;
    }
}

// variety of functions possible
double selection(size_t const nMito, size_t const nMut)
{
    double val = 0;

    if (type == 0)
    {
        val = nMut < .5 * nMito ? 
                1.0 - ch * 2.0 * nMut / nMito * 2.0 * nMut / nMito
                :
                1.0 - ch * 2.0 * (nMito - nMut) / nMito * 2.0 * (nMito - nMut) / nMito;
    } else if (type == 1)
    {
        val = 1.0 - ch * pow((double) nMut / nMito,2);
    } else 
    {
        val = nMut < .5 * nMito ? 
                1.0 - ch * (1.0 -  2.0 * nMut / nMito * 2.0 * nMut / nMito)
                :
                1.0 - ch * (1.0 - 2.0 * (nMito - nMut) / nMito * 2.0 * (nMito - nMut) / nMito);
    }

    // selection against heteroplasmy
    return(val);
}

void create_kid(Individual &parent, Individual &kid)
{
    kid.nMito = nMito_min
        + (gsl_rng_uniform(r) < 0.5 ? -1.0 : 1.0) * gsl_rng_uniform_int(r, nMito_error);

    // replicate mitochondria via sampling with replacement
    // in the model by Hadjivasiliou et al, (p4 in supplement)
    // a binomial prob was used, with a total sample of M mitochondria
    // with probability
    // of sampling a mutant mitochondrion
    // given by 2m/M. However, this is a bit odd, as it could potentially lead to situations
    // where a gamete with m mutant mitochondria eventually contributes fewer than 
    // m mitochondria to the zygote. 
    //
    // here we are going to sample only the kids mitochondria from the parental
    // ones 
    size_t nmito_mut_total = parent.nMitoMut 
        + gsl_ran_binomial(r, 
                (double) parent.nMitoMut / parent.nMito, kid.nMito);

    kid.nMitoMut = gsl_ran_hypergeometric(r, 
            nmito_mut_total,
            kid.nMito + parent.nMito - nmito_mut_total,
            kid.nMito
            );

    parent.nMitoMut = nmito_mut_total - kid.nMitoMut;

    assert(parent.nMitoMut >= 0);
}

// write the histogram data to a file
void write_data(
        gsl_histogram const *h, 
        size_t const popsize, 
        size_t const spore_number,
        size_t const generation)
{
    for (size_t i = 0; i < gsl_histogram_bins(h); ++i)
    {
        DataFile << generation << ";" 
        << popsize << ";" 
        << spore_number << ";"
        << (double)i / gsl_histogram_bins(h) << ";"
        << gsl_histogram_get(h, i) / popsize << ";" << endl;
    }
}

// write headers to the file
void init_file(gsl_histogram *h)
{
    DataFile << "generation;N;spore;p_i;freq;" << endl;
}

// cross two parents and grow them
void cross_and_grow()
{
    // generate an ascus of a pairing containing
    // 50% mutant mitochondria from parent 1
    // and 50% wt mitochondria from parent 2
    Individual * ascus = new Individual[4];

    make_tetrad(ascus, 4);

    // counters for stats
    size_t popsize = 0;

    // make a histogram to get an idea of the frequency distribution
    gsl_histogram * h = gsl_histogram_alloc(200);
    gsl_histogram_set_ranges_uniform(h, 0, 1.01);
    
    // initialize the data file
    init_file(h);
    
    // grow colonies from each spore
    for (size_t spore_i = 0; spore_i < 4; ++spore_i)
    {
        cout << "running spore " << spore_i << endl;
        // add the spore as the initial member of the population
        Pop.insert(
                Pop.end(), 
                Individual(ascus[spore_i].nMito, ascus[spore_i].nMitoMut)
                );

        // update the histogram
        gsl_histogram_increment(h, (double)ascus[spore_i].nMitoMut / ascus[spore_i].nMito);

        ++popsize;
        
        write_data(h, popsize, spore_i, 0);

        // now execute duplication (after selection) of each individual cell 
        for (int generation = 1; generation < 26; ++generation)
        {
            cout << "generation: " << generation << " " << popsize << endl;
            size_t popsize_new = popsize;

            for (size_t i = 0; i < popsize; ++i)
            {
                // perform weak selection. If individual survive let it reproduce
                if (gsl_rng_uniform(r) < selection(Pop[i].nMito, Pop[i].nMitoMut))
                {
                    Individual kid;

                    create_kid(Pop[i], kid);

                    // update the histogram
                    gsl_histogram_increment(h, (double)Pop[i].nMitoMut / Pop[i].nMito);

                    Pop.insert(Pop.end(), kid);

                    ++popsize_new;
                }
            }
            
            popsize = popsize_new;

            // output data
            write_data(h, popsize, spore_i, generation);

        }

        // done with one spore, now the other
        // so clear the population
        Pop.clear();
        popsize = 0;

        // and empty the histogram so that we can generate
        // new stats for the next spore
        gsl_histogram_reset(h);
    }

    // remove histogram
    gsl_histogram_free(h);
}


void write_parameters()
{
    DataFile << endl << endl
                << "type;" << type << endl
                << "seed;" << seed << endl
                << "nmito_error;" << nMito_error << endl
                << "nmito_min;" << nMito_min << endl
                << "ch;" << ch << endl;
}


// write the headers of the datafile
// initialize some parameters from command line arguments
void init_arguments(int argc, char *argv[])
{
    ch = atof(argv[1]);
    type = atoi(argv[2]);
    nMito_min = atoi(argv[3]);
    nMito_error = atoi(argv[4]);
}

int main(int argc, char **argv)
{
    // get command line arguments
    init_arguments(argc,argv);

    // initialize all members of the population
    init_pop();

    cross_and_grow();

    write_parameters();
}
