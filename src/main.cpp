#include <iostream>
#include <fstream>
#include "Molecule.h"
#include "FileInputParser.h"
#include "CNDO.h"
#include "overlapIntegrals.h"

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <input file>" << std::endl;
        return 1;
    }

    FileInputParser parser;
    Molecule molecule = parser.parseMolecule(argv[1]);

    arma::mat overlapMatrix = calcOverlapMatrix(molecule);

    CNDO cndo(molecule, overlapMatrix);
    cndo.scfCycle();

    arma::mat gradient = cndo.calcTotalEnergyGradMat();

    std::cout << "Nuclear Repulsion Energy: " << cndo.calcNuclearRepulsionEnergy() << " eV." << std::endl;
    std::cout << "Total Energy: " << cndo.calcTotalEnergy() << " eV." << std::endl;
    std::cout << "Gradient: " << std::endl
              << gradient << std::endl;

    return 0;
}
