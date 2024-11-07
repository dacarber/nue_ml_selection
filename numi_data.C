
/**
 * @file data.C
 * @brief ROOT macro to be used with CAFAna executable to run the selection on
 * data events.
 * @author justin.mueller@colostate.edu
*/

//#include "include/analysis.h"
//#include "include/datalog.h"
#include "include/definitions.h"
#include "include/cuts.h"
#include "include/variables.h"
#include "include/nue_variables.h"
#include "include/container.h"
#include "sbnana/CAFAna/Core/Binning.h"

#include <fstream>

using namespace ana;

#define OUT(STREAM,TAG) STREAM << std::fixed << TAG << ","
#define CSV(VAL) VAL << ","

std::ofstream output("output_data_crtpmt.log");
std::ofstream output_evt("output_evt.log");

/**
 * Writes reconstructed variables selected interactions.
 * @param sr is an SRSpillProxy that attaches to the StandardRecord of the
 * current spill.
 * @param j the reco interaction (selected).
 * @return None.
*/
void write_reco(const caf::SRSpillProxy* sr, const caf::SRInteractionDLPProxy& j)
{
    output  << CSV(sr->hdr.run) << CSV(sr->hdr.evt) << CSV(sr->hdr.subrun)
            << CSV(vars::image_id(j)) << CSV(vars::id(j))
            << CSV(vars::leading_electron_ke(j))
            << CSV(vars::leading_proton_ke(j))
            << CSV(vars::visible_energy(j))
            << CSV(vars::leading_electron_pt(j))
            << CSV(vars::leading_proton_pt(j))
            << CSV(vars::electron_polar_angle(j))
            << CSV(vars::electron_azimuthal_angle(j))
            << CSV(vars::opening_angle(j))
            << CSV(vars::interaction_pt(j))
            << CSV(vars::phiT(j))
            << CSV(vars::alphaT(j))
            << CSV(vars::electron_softmax(j))
            << CSV(vars::proton_softmax(j))
            << CSV(cuts::all_1e1p_cut(j))
            << CSV(cuts::all_1eNp_cut(j))
            << CSV(cuts::all_1eX_cut(j))
            //<< CSV(cuts::crtpmt_veto_data(sr))
            << CSV(j.volume_id)
            << std::endl;
}

/**
 * Writes the reconstructed variables for the selected interactions.
 * @param sr is an SRSpillProxy that attaches to the StandardRecord of the
 * current spill.
 * @return A dummy vector of doubles.
*/
const SpillMultiVar kDataInfo([](const caf::SRSpillProxy* sr)
{
    /**
     * Loop over reconstructed interactions and log interaction-level
     * information. 
    */
    for(auto const & i : sr->dlp)
    {
        if(cuts::all_1eX_cut(i) || cuts::all_1eNp_cut(i) || cuts::all_1e1p_cut(i))
        {
            OUT(output,"DATA");
            write_reco(sr, i);
        }
    }

    output_evt  << CSV(sr->hdr.run) << CSV(sr->hdr.evt) << CSV(sr->hdr.subrun) << std::endl;

    return std::vector<double>{1};
});

/**
 * Enumerates the cut that each interaction passes.
 * @param sr is an SRSpillProxy that attaches to the StandardRecord of the
 * current spill.
 * @return A vector of doubles containing the enumeration of cuts passed
 * by each interaction.
*/
const SpillMultiVar kOffbeam1mu1pCut([](const caf::SRSpillProxy* sr)
{
    std::vector<double> cut_vector;
    for(auto const & i : sr->dlp)
    {
        // No cut: 0, fiducial: 1, contained: 2, topological: 3, flash: 4
        double cut(0);
        if(cuts::fiducial_cut(i))
            cut = 1;
        if(cuts::containment_cut(i))
            cut = 2;
        if(cuts::fiducial_containment_topological_1e1p_cut(i))
            cut = 3;
        if(cuts::all_1e1p_cut(i))
            cut = 4;
        cut_vector.push_back(cut);
    }
    return cut_vector;
});

/**
 * Enumerates the cut that each interaction passes.
 * @param sr is an SRSpillProxy that attaches to the StandardRecord of the
 * current spill.
 * @return A vector of doubles containing the enumeration of cuts passed
 * by each interaction.
*/
const SpillMultiVar kOffbeam1muNpCut([](const caf::SRSpillProxy* sr)
{
    std::vector<double> cut_vector;
    for(auto const & i : sr->dlp)
    {
        // No cut: 0, fiducial: 1, contained: 2, topological: 3, flash: 4
        double cut(0);
        if(cuts::fiducial_cut(i))
            cut = 1;
        if(cuts::containment_cut(i))
            cut = 2;
        if(cuts::fiducial_containment_topological_1eNp_cut(i))
            cut = 3;
        if(cuts::all_1eNp_cut(i))
            cut = 4;
        cut_vector.push_back(cut);
    }
    return cut_vector;
});

/**
 * Enumerates the cut that each interaction passes.
 * @param sr is an SRSpillProxy that attaches to the StandardRecord of the
 * current spill.
 * @return A vector of doubles containing the enumeration of cuts passed
 * by each interaction.
*/
const SpillMultiVar kOffbeam1muXCut([](const caf::SRSpillProxy* sr)
{
    std::vector<double> cut_vector;
    for(auto const & i : sr->dlp)
    {
        // No cut: 0, fiducial: 1, contained: 2, topological: 3, flash: 4
        double cut(0);
        if(cuts::fiducial_cut(i))
            cut = 1;
        if(cuts::containment_cut(i))
            cut = 2;
        if(cuts::fiducial_containment_topological_1eX_cut(i))
            cut = 3;
        if(cuts::all_1eX_cut(i))
            cut = 4;
        cut_vector.push_back(cut);
    }
    return cut_vector;
});

const SpillMultiVar kHandscanInfo([](const caf::SRSpillProxy* sr)
{
    /**
     * Loop over reconstructed interactions and log interaction-level
     * information. No truth information can be used.
    */
    for(auto const & i : sr->dlp)
    {
        if(cuts::topological_1eNp_cut(i))
        {
            size_t leading_eon(vars::leading_particle_index(i, 2));
            size_t leading_proton(vars::leading_particle_index(i, 4));
            OUT(output,"INTERACTION")   << CSV(sr->hdr.run) << CSV(sr->hdr.evt)
                                        << CSV(vars::image_id(i)) << CSV(vars::id(i))
                                        << CSV(vars::cryostat(i)) << CSV(i.is_fiducial)
                                        << CSV(i.is_contained) << CSV(cuts::topology(i))
                                        << CSV(cuts::flash_cut_data(i))
                                        << CSV(i.vertex[0]) << CSV(i.vertex[1]) << CSV(i.vertex[2])
                                        //<< CSV(i.particles[leading_eon].length)
                                        << CSV(vars::leading_electron_ke(i))
                                        << CSV(i.particles[leading_proton].length)
                                        << CSV(vars::leading_proton_ke(i))
                                        << CSV(vars::flash_time(i))
                                        //<< CSV(i.particles[leading_electron].end_point[0])
                                        //<< CSV(i.particles[leading_electron].end_point[1])
                                        //<< CSV(i.particles[leading_electron].end_point[2])
                                        //<< CSV(i.particles[leading_electron].start_dir[0])
                                        //<< CSV(i.particles[leading_electron].start_dir[1])
                                        //<< CSV(i.particles[leading_electron].start_dir[2])
                                        //<< CSV(i.particles[leading_proton].end_point[0])
                                        //<< CSV(i.particles[leading_proton].end_point[1])
                                        //<< CSV(i.particles[leading_proton].end_point[2])
                                        //<< CSV(i.particles[leading_proton].start_dir[0])
                                        //<< CSV(i.particles[leading_proton].start_dir[1])
                                        //<< CSV(i.particles[leading_proton].start_dir[2])
                                        << std::endl;
        }
    }

    return std::vector<double>{1};
});

/**
 * The main function of the selection. Creates a container for the CAFAna
 * Spectrum objects and populates it with a variety of variables that define
 * the selection.
 * @return none.
*/
void numi_data()
{

    SpecContainer spectra("/pnfs/icarus/persistent/users/lkashur/numi_run2_v09_89_01p01/output_flat/*.flat.root", "spectra_data_new_weights.root", -1, -1);



    spectra.add_spectrum1d("sDataInfo", Binning::Simple(1, 0, 2), kDataInfo);

    spectra.run();
}
