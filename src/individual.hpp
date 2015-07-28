// individual class

#ifndef __INDIVIDUAL_INCLUDED__
#define __INDIVIDUAL_INCLUDED__

#include <cstddef>

class Individual
{

    public:
        size_t nMito;
        size_t nMitoMut;

    Individual(std::size_t const nMito, std::size_t const nMitoMut);
    
    Individual();

    Individual(Individual const &Individual);
};

#endif
