#include "individual.hpp"
#include <cstddef>

Individual::Individual(std::size_t const nMito, std::size_t const nMitoMut) : nMito(nMito), nMitoMut(nMitoMut)
{
}

Individual::Individual() 
{
    nMito = 0;
    nMitoMut = 0;
}

Individual::Individual(Individual const &Individual) {

    nMito = Individual.nMito;
    nMitoMut = Individual.nMitoMut;

}


