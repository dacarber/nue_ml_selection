/**
 * @file analysis.C
 * @brief ROOT macro to be used with CAFAna executable to run the selection.
 * @author daniel.carber@colostate.edu
*/

#include "include/analysis.h"
#include "include/container.h"
#include "include/csv_maker.h"
#include "sbnana/CAFAna/Core/Binning.h"

using namespace ana;

/**
 * The main function of the selection. Creates a container for the CAFAna
 * Spectrum objects and populates it with a variety of variables that define
 * the selection.
 * @return none.
*/
void analysis()
{
    /**
     * 1. NuMI Nu neutrino (full flux) + out-of-time cosmics collection charge only(v09_89_01p01).
    */

    SpecContainer spectra("/pnfs/icarus/persistent/users/lkashur/numi_nu_cosmic_v09_89_01p01/output_flat/*.flat.root", "spectra_nucosmics.root", -1, 9.36E19 ); //Collection plane only charge
    
    spectra.add_spectrum1d("sSelected", Binning::Simple(1, 0, 2), kInfoVar);

    spectra.run();
}