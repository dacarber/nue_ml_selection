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
void analysis_csv()
{
    /**
     * 1. NuMI Nu neutrino (full flux) + out-of-time cosmics all plane(v09_89_01p01).   
     * 2. NuMI Nu neutrino (full flux) + out-of-time cosmics collection charge only(v09_89_01p01).
    */

    //SpecContainer spectra("/pnfs/icarus/scratch/users/dcarber/hdf5_files/NuMI_Nu/v09_89_01p01/merged_caf/*.flat.root", "spectra_nucosmics.root", -1, 9.36E19 ); //All plane charge
    SpecContainer spectra("/pnfs/icarus/persistent/users/dcarber/Nue_analysis/flat_cafs/collection_plane_v09_89_01p01/NuMI_nu_corsika_collection_part_*.root", "spectra_nucosmics.root", -1, -1 ); //Collection plane only charge
    spectra.add_spectrum1d("sSelected", Binning::Simple(1, 0, 2), kInfoVar);

    spectra.run();
}