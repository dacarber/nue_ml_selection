/**
 * @file analysis.C
 * @brief ROOT macro to be used with CAFAna executable to run the selection.
 * @author justin.mueller@colostate.edu
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
     * 1. BNB neutrino (full flux) + out-of-time cosmics (v09_63_01).
     * 2. BNB in-time cosmics + out-of-time cosmics (v09_63_01).
    */
    //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/bnb_nucosmics_v6.flat.root", "spectra_nucosmics.root", 1.253e19, 2.5e20);
    //SpecContainer spectra("/exp/icarus/data/users/mueller/mlcafs/bnb_intime.flat.root", "spectra_intime.root", 9070*2.05e14, 2.5e20);
    //SpecContainer spectra("/exp/icarus/data/users/dcarber/numi_nue/sample_numi_nue.flat.root", "spectra_nuecosmics.root", 8.7e19, 8.7e19 );
    //SpecContainer spectra("/pnfs/icarus/scratch/users/dcarber/hdf5_files/numi_nu_corsika_230918/larcv_*.0_mlreco_ana.flat.root", "spectra_nucosmics.root", 1.83e19, 1.83e19 );
    SpecContainer spectra("/pnfs/icarus/scratch/users/dcarber/hdf5_files/numi_intime_cosmics_230918/lite_files/flat_cafs/larcv_72733929.0_mlreco_.flat.root", "spectra_intimecosmics.root", 1.83e19, 1.83e19 );
    /**
     * 3. BNB neutrino (full flux) + out-of-time cosmics *     Central Value    * (v09_82_02_01).
     * 4. BNB neutrino (full flux) + out-of-time cosmics * Coherent Noise +4.5% * (v09_82_02_01).
     * 5. BNB neutrino (full flux) + out-of-time cosmics *  Elli. Recombination * (v09_82_02_01).
     * 6. BNB neutrino (full flux) + out-of-time cosmics * Untuned Signal Shape * (v09_82_02_01).
    */
    //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/ml_hdf5/bnb_nu_sys/systematics_cv_v2.flat.root", "spectra_cv.root", -1, 2.5e20);
    //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/ml_hdf5/bnb_nu_sys/systematics_tpcnoise_coh_p1_v2.flat.root", "spectra_tpcnoise_coh_p1.root", -1, 2.5e20);
    //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/ml_hdf5/bnb_nu_sys/systematics_recombination.flat.root", "spectra_recombination.root", -1, 2.5e20);
    //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/ml_hdf5/bnb_nu_sys/systematics_untunedsigshape.flat.root", "spectra_untunedsigshape.root", -1, 2.5e20);
    
    /**
     * 7. BNB neutrino-only (full flux)  *     Central Value    * (v09_82_02_01).
     * 8. BNB neutrino-only (full flux)  * Coherent Noise +4.5% * (v09_82_02_01).
     * 9. BNB neutrino-only (full flux)  *  Elli. Recombination * (v09_82_02_01).
     * 10. BNB neutrino-only (full flux) * Untuned Signal Shape * (v09_82_02_01).
    */
    //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/systematics/sample_cv.flat.root", "spectra_cv.root", -1, 2.5e20);
    //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/systematics/sample_cohnoise.flat.root", "spectra_tpcnoise_coh_p1.root", -1, 2.5e20);
    //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/systematics/sample_intnoise.flat.root", "spectra_intnoise.root", -1, 2.5e20);
    //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/systematics/sample_recombination.flat.root", "spectra_recombination.root", -1, 2.5e20);
    //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/systematics/sample_sigshape.flat.root", "spectra_untunedsigshape.root", -1, 2.5e20);
    
    /**
     * 11. MPV/MPR sample (v09_63_00).
    */
    //SpecContainer spectra("/exp/icarus/data/users/mueller/mlcafs/mpv_mpr.flat.root", "spectra_mpvmpr.root", 1e19, 2.5e20);

    /**
     * Spectra (1D) for interactions.
    */
    spectra.add_spectrum1d("sCountParticles", Binning::Simple(20, 0, 20), kCountParticles);
    spectra.add_spectrum1d("sCountPrimaries", Binning::Simple(20, 0, 20), kCountPrimaries);
    spectra.add_spectrum1d("sCountParticlesTruth", Binning::Simple(20, 0, 20), kCountParticlesTruth);
    spectra.add_spectrum1d("sCountPrimariesTruth", Binning::Simple(20, 0, 20), kCountPrimariesTruth);
    spectra.add_spectrum1d("sEnergy_1e1p_signal_bias", Binning::Simple(50,-1,1), kEnergy_1e1p_signal_bias);
    spectra.add_spectrum1d("sEnergy_1e1p_othernu_bias", Binning::Simple(50,-1,1), kEnergy_1e1p_othernu_bias);
    spectra.add_spectrum1d("sEnergy_1e1p_cosmic_bias", Binning::Simple(50,-1,1), kEnergy_1e1p_cosmic_bias);
    spectra.add_spectrum1d("sEnergy_1eNp_1p_signal_bias", Binning::Simple(50,-1,1), kEnergy_1eNp_1p_signal_bias);
    spectra.add_spectrum1d("sEnergy_1eNp_Np_signal_bias", Binning::Simple(50,-1,1), kEnergy_1eNp_Np_signal_bias);
    spectra.add_spectrum1d("sEnergy_1eNp_othernu_bias", Binning::Simple(50,-1,1), kEnergy_1eNp_othernu_bias);
    spectra.add_spectrum1d("sEnergy_1eNp_cosmic_bias", Binning::Simple(50,-1,1), kEnergy_1eNp_cosmic_bias);

    spectra.add_spectrum1d("sNuEnergy_1e1p_signal_bias", Binning::Simple(50,-1,1), kNuEnergy_1e1p_signal_bias);
    spectra.add_spectrum1d("sNuEnergy_1e1p_othernu_bias", Binning::Simple(50,-1,1), kNuEnergy_1e1p_othernu_bias);
    spectra.add_spectrum1d("sNuEnergy_1eNp_1p_signal_bias", Binning::Simple(50,-1,1), kNuEnergy_1eNp_1p_signal_bias);
    spectra.add_spectrum1d("sNuEnergy_1eNp_Np_signal_bias", Binning::Simple(50,-1,1), kNuEnergy_1eNp_Np_signal_bias);
    spectra.add_spectrum1d("sNuEnergy_1eNp_othernu_bias", Binning::Simple(50,-1,1), kNuEnergy_1eNp_othernu_bias);
    spectra.add_spectrum1d("sCos_open_angle", Binning::Simple(50,0,3.2), kOpeningAngle_1e1p);
    spectra.add_spectrum1d("sP_momentum_bias", Binning::Simple(50,-1,1), kP_momentum_bias);
    spectra.add_spectrum1d("sAz_proton_angle_bias", Binning::Simple(50,-1,1), kAz_proton_angle_bias);
    spectra.add_spectrum1d("sPol_proton_angle_bias", Binning::Simple(50,-1,1), kPol_proton_angle_bias);
    spectra.add_spectrum1d("sAz_electron_angle_bias", Binning::Simple(50,-1,1), kAz_electron_angle_bias);

    spectra.add_spectrum1d("sNuMI_electron_angle_bias", Binning::Simple(50,0,3.2), kNuMI_electron_angle_bias);

    spectra.add_spectrum1d("sPol_electron_angle_bias", Binning::Simple(50,-1,1), kPol_electron_angle_bias);

    spectra.add_spectrum1d("sAz_proton_angle", Binning::Simple(50,0,3.2), kProtonAz_1e1p);
    spectra.add_spectrum1d("sPol_proton_angle", Binning::Simple(50,0,3.2), kProtonPolar_1e1p);
    spectra.add_spectrum1d("sAz_electron_angle", Binning::Simple(50,0,3.2), kElectronAz_1e1p);
    spectra.add_spectrum1d("sPol_electron_angle", Binning::Simple(50,0,3.2), kElectronPolar_1e1p);
    /**
     * Spectra (2D) for counting selection statistics by interaction categorization (efficiency).
    */
    spectra.add_spectrum2d("sCountTTP_NoCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryTTP_NoCut, kCountTTP_NoCut);
    spectra.add_spectrum2d("sCountTTP_FVCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryTTP_FVCut, kCountTTP_FVCut);
    spectra.add_spectrum2d("sCountTTP_FVConCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryTTP_FVConCut, kCountTTP_FVConCut);
    spectra.add_spectrum2d("sCountTTP_FVConTop1e1pCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryTTP_FVConTop1e1pCut, kCountTTP_FVConTop1e1pCut);
    spectra.add_spectrum2d("sCountTTP_All1e1pCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryTTP_All1e1pCut, kCountTTP_All1e1pCut);
    spectra.add_spectrum2d("sCountTTP_FVConTop1eNpCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryTTP_FVConTop1eNpCut, kCountTTP_FVConTop1eNpCut);
    spectra.add_spectrum2d("sCountTTP_All1eNpCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryTTP_All1eNpCut, kCountTTP_All1eNpCut);
    spectra.add_spectrum2d("sCountTTP_FVConTop1eXCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryTTP_FVConTop1eXCut, kCountTTP_FVConTop1eXCut);
    spectra.add_spectrum2d("sCountTTP_All1eXCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryTTP_All1eXCut, kCountTTP_All1eXCut);

    /**
     * Spectra (2D) for counting selection statistics by interaction categorization (purity).
    */
    spectra.add_spectrum2d("sCountPTT_NoCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryPTT_NoCut, kCountPTT_NoCut);
    spectra.add_spectrum2d("sCountPTT_FVCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryPTT_FVCut, kCountPTT_FVCut);
    spectra.add_spectrum2d("sCountPTT_FVConCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryPTT_FVConCut, kCountPTT_FVConCut);
    spectra.add_spectrum2d("sCountPTT_FVConTop1e1pCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryPTT_FVConTop1e1pCut, kCountPTT_FVConTop1e1pCut);
    spectra.add_spectrum2d("sCountPTT_All1e1pCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryPTT_All1e1pCut, kCountPTT_All1e1pCut);
    spectra.add_spectrum2d("sCountPTT_FVConTop1eNpCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryPTT_FVConTop1eNpCut, kCountPTT_FVConTop1eNpCut);
    spectra.add_spectrum2d("sCountPTT_All1eNpCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryPTT_All1eNpCut, kCountPTT_All1eNpCut);
    spectra.add_spectrum2d("sCountPTT_FVConTop1eXCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryPTT_FVConTop1eXCut, kCountPTT_FVConTop1eXCut);
    spectra.add_spectrum2d("sCountPTT_All1eXCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryPTT_All1eXCut, kCountPTT_All1eXCut);

    /**
     * Spectra (2D) for visible energy.
    */
    spectra.add_spectrum2d("sVisibleEnergyTTP_NoCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryTTP_NoCut, kVisibleEnergyTTP_NoCut);
    spectra.add_spectrum2d("sVisibleEnergyTTP_FVCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryTTP_FVCut, kVisibleEnergyTTP_FVCut);
    spectra.add_spectrum2d("sVisibleEnergyTTP_FVConCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryTTP_FVConCut, kVisibleEnergyTTP_FVConCut);
    spectra.add_spectrum2d("sVisibleEnergyTTP_FVConTop1e1pCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryTTP_FVConTop1e1pCut, kVisibleEnergyTTP_FVConTop1e1pCut);
    spectra.add_spectrum2d("sVisibleEnergyTTP_All1e1pCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryTTP_All1e1pCut, kVisibleEnergyTTP_All1e1pCut);
    spectra.add_spectrum2d("sVisibleEnergyTTP_FVConTop1eNpCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryTTP_FVConTop1eNpCut, kVisibleEnergyTTP_FVConTop1eNpCut);
    spectra.add_spectrum2d("sVisibleEnergyTTP_All1eNpCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryTTP_All1eNpCut, kVisibleEnergyTTP_All1eNpCut);
    spectra.add_spectrum2d("sVisibleEnergyTTP_FVConTop1meXCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryTTP_FVConTop1eXCut, kVisibleEnergyTTP_FVConTop1eXCut);
    spectra.add_spectrum2d("sVisibleEnergyTTP_All1eXCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryTTP_All1eXCut, kVisibleEnergyTTP_All1eXCut);

    /**
     * Spectra (2D) for flash time.
    */
    spectra.add_spectrum2d("sFlashTime_NoCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -1000, 1000), kCategoryPTT_NoCut, kFlashTimePTT_NoCut);
    spectra.add_spectrum2d("sFlashTime_FVCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -1000, 1000), kCategoryPTT_FVCut, kFlashTimePTT_FVCut);
    spectra.add_spectrum2d("sFlashTime_FVConCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -1000, 1000), kCategoryPTT_FVConCut, kFlashTimePTT_FVConCut);
    spectra.add_spectrum2d("sFlashTime_FVConTop1e1pCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -1000, 1000), kCategoryPTT_FVConTop1e1pCut, kFlashTimePTT_FVConTop1e1pCut);
    spectra.add_spectrum2d("sFlashTime_All1e1pCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -1000, 1000), kCategoryPTT_All1e1pCut, kFlashTimePTT_All1e1pCut);
    spectra.add_spectrum2d("sFlashTime_FVConTop1eNpCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -1000, 1000), kCategoryPTT_FVConTop1eNpCut, kFlashTimePTT_FVConTop1eNpCut);
    spectra.add_spectrum2d("sFlashTime_All1eNpCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -1000, 1000), kCategoryPTT_All1eNpCut, kFlashTimePTT_All1eNpCut);
    spectra.add_spectrum2d("sFlashTime_FVConTop1eXCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -1000, 1000), kCategoryPTT_FVConTop1eXCut, kFlashTimePTT_FVConTop1eXCut);
    spectra.add_spectrum2d("sFlashTime_All1eXCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -1000, 1000), kCategoryPTT_All1eXCut, kFlashTimePTT_All1eXCut);
    spectra.add_spectrum2d("sFlashTime_Zoomed_NoCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -4, 10), kCategoryPTT_NoCut, kFlashTimePTT_NoCut);
    spectra.add_spectrum2d("sFlashTime_Zoomed_FVCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -4, 10), kCategoryPTT_FVCut, kFlashTimePTT_FVCut);
    spectra.add_spectrum2d("sFlashTime_Zoomed_FVConCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -4, 10), kCategoryPTT_FVConCut, kFlashTimePTT_FVConCut);
    spectra.add_spectrum2d("sFlashTime_Zoomed_FVConTop1e1pCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -4, 10), kCategoryPTT_FVConTop1e1pCut, kFlashTimePTT_FVConTop1e1pCut);
    spectra.add_spectrum2d("sFlashTime_Zoomed_All1e1pCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -4, 10), kCategoryPTT_All1e1pCut, kFlashTimePTT_All1e1pCut);
    spectra.add_spectrum2d("sFlashTime_Zoomed_FVConTop1eNpCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -4, 10), kCategoryPTT_FVConTop1eNpCut, kFlashTimePTT_FVConTop1eNpCut);
    spectra.add_spectrum2d("sFlashTime_Zoomed_All1eNpCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -4, 10), kCategoryPTT_All1eNpCut, kFlashTimePTT_All1eNpCut);
    spectra.add_spectrum2d("sFlashTime_Zoomed_FVConTop1eXCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -4, 10), kCategoryPTT_FVConTop1eXCut, kFlashTimePTT_FVConTop1eXCut);
    spectra.add_spectrum2d("sFlashTime_Zoomed_All1eXCut", Binning::Simple(10, 0, 10), Binning::Simple(50, -4, 10), kCategoryPTT_All1eXCut, kFlashTimePTT_All1eXCut);

    /**
     * Spectra (2D) for (stacked) reconstructed quantities.
    */
    spectra.add_spectrum2d("sFlashTimePTT_NoCut", Binning::Simple(10, 0, 10), Binning::Simple(60, -4, 10), kCategoryTopologyPTT_NoCut, kFlashTimePTT_NoCut);
    spectra.add_spectrum2d("sVisibleEnergyPTT_Topology_All1e1pCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryTopologyPTT_All1e1pCut, kVisibleEnergyPTT_All1e1pCut);
    spectra.add_spectrum2d("sVisibleEnergyPTT_InteractionMode_All1e1pCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryInteractionModePTT_All1e1pCut, kVisibleEnergyPTT_All1e1pCut);
    spectra.add_spectrum2d("sFlashTimePTT_Topology_All1e1pCut", Binning::Simple(10, 0, 10), Binning::Simple(60, -4, 10), kCategoryTopologyPTT_All1e1pCut, kFlashTimePTT_All1e1pCut);
    spectra.add_spectrum2d("sVisibleEnergyPTT_Topology_All1eNpCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryTopologyPTT_All1eNpCut, kVisibleEnergyPTT_All1eNpCut);
    spectra.add_spectrum2d("sVisibleEnergyPTT_InteractionMode_All1eNpCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryInteractionModePTT_All1eNpCut, kVisibleEnergyPTT_All1eNpCut);
    spectra.add_spectrum2d("sFlashTimePTT_Topology_All1eNpCut", Binning::Simple(10, 0, 10), Binning::Simple(60, -4, 10), kCategoryTopologyPTT_All1eNpCut, kFlashTimePTT_All1eNpCut);
    spectra.add_spectrum2d("sVisibleEnergyPTT_All1eXCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryPTT_All1eXCut, kVisibleEnergyPTT_All1eXCut);
    spectra.add_spectrum2d("sVisibleEnergyPTT_Topology_All1eXCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryTopologyPTT_All1eXCut, kVisibleEnergyPTT_All1eXCut);
    spectra.add_spectrum2d("sVisibleEnergyPTT_InteractionMode_All1eXCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryInteractionModePTT_All1eXCut, kVisibleEnergyPTT_All1eXCut);
    spectra.add_spectrum2d("sFlashTimePTT_Topology_All1eXCut", Binning::Simple(10, 0, 10), Binning::Simple(60, -4, 10), kCategoryTopologyPTT_All1eXCut, kFlashTimePTT_All1eXCut);

    spectra.add_spectrum2d("sFlashTimeTTP_NoCut", Binning::Simple(10, 0, 10), Binning::Simple(60, -4, 10), kCategoryTopologyTTP_NoCut, kFlashTimeTTP_NoCut);
    spectra.add_spectrum2d("sVisibleEnergyTTP_Topology_All1e1pCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryTopologyTTP_All1e1pCut, kVisibleEnergyTTP_All1e1pCut);
    spectra.add_spectrum2d("sVisibleEnergyTTP_InteractionMode_All1e1pCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryInteractionModeTTP_All1e1pCut, kVisibleEnergyTTP_All1e1pCut);
    spectra.add_spectrum2d("sFlashTimeTTP_Topology_All1e1pCut", Binning::Simple(10, 0, 10), Binning::Simple(60, -4, 10), kCategoryTopologyTTP_All1e1pCut, kFlashTimeTTP_All1e1pCut);
    spectra.add_spectrum2d("sVisibleEnergyTTP_Topology_All1eNpCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryTopologyTTP_All1eNpCut, kVisibleEnergyTTP_All1eNpCut);
    spectra.add_spectrum2d("sVisibleEnergyTTP_InteractionMode_All1eNpCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryInteractionModeTTP_All1eNpCut, kVisibleEnergyTTP_All1eNpCut);
    spectra.add_spectrum2d("sFlashTimeTTP_Topology_All1eNpCut", Binning::Simple(10, 0, 10), Binning::Simple(60, -4, 10), kCategoryTopologyTTP_All1eNpCut, kFlashTimeTTP_All1eNpCut);
    spectra.add_spectrum2d("sVisibleEnergyTTP_All1eXCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryTTP_All1eXCut, kVisibleEnergyTTP_All1eXCut);
    spectra.add_spectrum2d("sVisibleEnergyTTP_Topology_All1eXCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryTopologyTTP_All1eXCut, kVisibleEnergyTTP_All1eXCut);
    spectra.add_spectrum2d("sVisibleEnergyTTP_InteractionMode_All1eXCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryInteractionModeTTP_All1eXCut, kVisibleEnergyTTP_All1eXCut);
    spectra.add_spectrum2d("sFlashTimeTTP_Topology_All1eXCut", Binning::Simple(10, 0, 10), Binning::Simple(60, -4, 10), kCategoryTopologyTTP_All1eXCut, kFlashTimeTTP_All1eXCut);

    /**
     * Spectra (2D) for particles.
    */
    spectra.add_spectrum2d("sCal_electron", Binning::Simple(50, 50, 1000), Binning::Simple(50, 50, 1000), kCalTruth_electron, kCal_electron);
    spectra.add_spectrum2d("sCal_electron2electron", Binning::Simple(50, 50, 1000), Binning::Simple(50, 50, 1000), kCalTruth_electron, kCal_electron2electron);
    spectra.add_spectrum2d("sCal_electron_bias2d", Binning::Simple(10, 0, 1000), Binning::Simple(250,-0.25,0.25), kCalTruth_electron, kCal_electron_bias);
    spectra.add_spectrum1d("sCal_electron_bias", Binning::Simple(75,-1,1), kCal_electron_bias);
    spectra.add_spectrum1d("sCal_noncc_electron_bias", Binning::Simple(75,-1,1), kCal_noncc_electron_bias);
    spectra.add_spectrum1d("sCal_wellreco_electron_bias", Binning::Simple(75,-1,1), kCal_wellreco_electron_bias);
    spectra.add_spectrum1d("sCal_wellreco_proton_bias", Binning::Simple(75,-1,1), kCal_wellreco_proton_bias);
    spectra.add_spectrum1d("sCCOverlap", Binning::Simple(50, 0, 1), kCCOverlap);
    spectra.add_spectrum1d("sNonCCOverlap", Binning::Simple(50, 0, 1), kNonCCOverlap);



    /**
     * Spectra (2D) for matched (truth-to-predicted) particles.
    */
    spectra.add_spectrum2d("sPrimary_confusion", Binning::Simple(2,0,2), Binning::Simple(2,0,2), kPrimaryTruth, kPrimary);
    spectra.add_spectrum2d("sPID_confusion", Binning::Simple(5,0,5), Binning::Simple(5,0,5), kPIDTruth, kPID);
    spectra.add_spectrum2d("sPrimaryPID_confusion", Binning::Simple(10,0,10), Binning::Simple(10,0,10), kPrimaryPIDTruth, kPrimaryPID);
    spectra.add_spectrum2d("sPrimary_Neutrino_confusion", Binning::Simple(2,0,2), Binning::Simple(2,0,2), kPrimaryTruth_Neutrino, kPrimary_Neutrino);
    spectra.add_spectrum2d("sPID_Neutrino_confusion", Binning::Simple(5,0,5), Binning::Simple(5,0,5), kPIDTruth_Neutrino, kPID_Neutrino);
    spectra.add_spectrum2d("sPrimaryPID_Neutrino_confusion", Binning::Simple(10,0,10), Binning::Simple(10,0,10), kPrimaryPIDTruth_Neutrino, kPrimaryPID_Neutrino);
    spectra.add_spectrum2d("sPrimary_Cosmic_confusion", Binning::Simple(2,0,2), Binning::Simple(2,0,2), kPrimaryTruth_Cosmic, kPrimary_Cosmic);
    spectra.add_spectrum2d("sPID_Cosmic_confusion", Binning::Simple(5,0,5), Binning::Simple(5,0,5), kPIDTruth_Cosmic, kPID_Cosmic);
    spectra.add_spectrum2d("sPrimaryPID_Cosmic_confusion", Binning::Simple(10,0,10), Binning::Simple(10,0,10), kPrimaryPIDTruth_Cosmic, kPrimaryPID_Cosmic);
    spectra.add_spectrum2d("sLowX", Binning::Simple(100,-400,400), Binning::Simple(100,-400,400), kLowX, kLowXTruth);

    spectra.add_spectrum2d("sPrimaryWellReco_confusion", Binning::Simple(2,0,2), Binning::Simple(2,0,2), kPrimaryWellRecoTruth, kPrimaryWellReco);
    spectra.add_spectrum2d("sPIDWellReco_confusion", Binning::Simple(5,0,5), Binning::Simple(5,0,5), kPIDWellRecoTruth, kPIDWellReco);
    spectra.add_spectrum2d("sPrimaryPIDWellReco_confusion", Binning::Simple(10,0,10), Binning::Simple(10,0,10), kPrimaryPIDWellRecoTruth, kPrimaryPIDWellReco);

    spectra.add_spectrum2d("sPrimaryWellReco_Neutrino_confusion", Binning::Simple(2,0,2), Binning::Simple(2,0,2), kPrimaryWellRecoTruth_Neutrino, kPrimaryWellReco_Neutrino);
    spectra.add_spectrum2d("sPIDWellReco_Neutrino_confusion", Binning::Simple(5,0,5), Binning::Simple(5,0,5), kPIDWellRecoTruth_Neutrino, kPIDWellReco_Neutrino);
    spectra.add_spectrum2d("sPrimaryPIDWellReco_Neutrino_confusion", Binning::Simple(10,0,10), Binning::Simple(10,0,10), kPrimaryPIDWellRecoTruth_Neutrino, kPrimaryPIDWellReco_Neutrino);

    /**
     * Spectra (2D) for correlating truth quantities.
    */
    spectra.add_spectrum2d("sScatteringProtonOverlap", Binning::Simple(50, 0.25, 1), Binning::Simple(25, 0, 1), kProtonScattering, kLeadingProtonOverlap);

    /**
     * Dummy spectra for dumping particle-level information to a CSV log file.
    */
    spectra.add_spectrum1d("sSelected", Binning::Simple(1, 0, 2), kInfoVar);
    //spectra.add_spectrum1d("sSignal", Binning::Simple(1, 0, 2), kSignal);

    spectra.run();
}