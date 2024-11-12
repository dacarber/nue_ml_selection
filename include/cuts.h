/**
 * @file cuts.h
 * @brief Header file for definitions of selection cuts.
 * @author justin.mueller@colostate.edu 
*/
#ifndef CUTS_H
#define CUTS_H

#include <functional>
#include <vector>
#include <string>
#include <sstream>
#include <numeric>

namespace cuts
{
    /**
     * Apply a cut on whether a match exists.
     * @tparam T the object type (true or reco, interaction or particle).
     * @param obj the object to select on.
     * @return true if the object is matched.
    */
    template<class T>
        bool matched(const T & obj) { return obj.is_matched == 1; }

    /**
     * Apply a cut on the quality of the reconstruction using the overlap
     * fraction as a discriminant.
     * @tparam T the object type (true or reco, interaction or particle).
     * @param obj the object to selection on.
     * @return true if the object is well reconstructed.
    */
    template<class T>
        bool wellreco(const T & obj) { return matched(obj) && obj.match_ids[0] > 0.95; }

    /**
     * Apply a cut on the validity of the flash match.
     * @tparam T the type of interaction (true or reco).
     * @param interaction on which to place the flash validity cut.
     * @return true if the interaction is flash matched and the time is valid.
    */
    template<class T>
        bool valid_flashmatch(const T & interaction) { return !std::isnan(interaction.flash_time) && interaction.is_flash_matched == 1; }

    /**
     * Count the primaries of the interaction with cuts applied to each particle.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to find the topology of.
     * @return the count of primaries of each particle type within the
     * interaction.
     */
    template<class T>
        std::vector<uint32_t> count_primaries(const T & interaction)
        {
            std::vector<uint32_t> counts(5, 0);
            for(auto &p : interaction.particles)
            {
                if(p.is_primary)
                {
                    double energy(p.pid > 1 ? p.csda_ke : p.calo_ke);
                    if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                        energy = p.energy_init;
                    if(p.pid == 2 && energy > 25) // Muon greater than 50 cm.
                        ++counts[p.pid];
                    else if((p.pid == 1 && energy >= 70) || (p.pid == 4 && energy > 40) || (p.pid == 3 && energy >25) || (p.pid == 0 && energy >25)) //
                        ++counts[p.pid];
                }
            }
            return counts;
        }

    /**
     * Check if the particle meets final state signal requirements.
     * Particles must be primary and have an energy above threshold.
     * Muons must have a length of at least 50 cm (143.425 MeV), protons
     * must have an energy above 50 MeV, and all other particles must have
     * an energy above 25 MeV.
     * @tparam T the type of particle (true or reco).
     * @param particle to check.
     * @return true if the particle is a final state signal particle.
    */
    template<class T>
        bool final_state_signal(const T & p)
        {
            bool passes(false);
            if(p.is_primary)
            {
                double energy(p.pid > 1 ? p.csda_ke : p.calo_ke);
                if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
                    energy = p.energy_init;

                if((p.pid == 1 && energy > 70) || (p.pid != 1 && p.pid < 4 && energy > 25) || (p.pid == 4 && energy > 40))
                    passes = true;
            }
            return passes;
        }

    /**
     * Find the topology of the interaction with cuts applied to each particle.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to find the topology of.
     * @return the topology of the interaction as a string (e.g 0ph0e1mu0pi1p).
     */
    template<class T>
        std::string topology(const T & interaction)
        {
            std::vector<uint32_t> counts(count_primaries(interaction));
            std::stringstream ss;
            ss  << counts[0] << "ph"
                << counts[1] << "e"
                << counts[2] << "mu"
                << counts[3] << "pi"
                << counts[4] << "p";
            return ss.str();
        }

    /**
     * Apply no cut (all interactions/particles passed).
     * @tparam T the type of object (true or reco, interaction or particle).
     * @param interaction/particle to select on.
     * @return true (always).
     */
    template<class T>
        bool no_cut(const T & obj) { return true; }

    /**
     * Apply a fiducial volume cut. Interaction vertex must be within 25 cm of
     * x and y detector faces, 50 cm of downstream (+) z face, and 30 cm of
     * upstream (-) z face.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the vertex is in the fiducial volume.
     */
    template<class T>
        bool fiducial_cut(const T & interaction) { return interaction.is_fiducial; }
    
    /**
     * Apply a containment volume cut. All points within the interaction must be
     * at least 5 cm from the detector boundaries.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the vertex is contained.
     */
    template<class T>
        bool containment_cut(const T & interaction) {
            bool passes(true); 
            for(const auto & p : interaction.particles){
                if(p.is_primary)
                {
                    if((p.pid > 1 && p.is_contained == false))
                        passes = false;
                }
                
                //return interaction.is_contained; 
            }
            return passes;
        }

    /**
     * Apply a 1mu1p topological cut. The interaction must have a topology
     * matching 1mu1p as defined by the conditions in the topology() function.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction has a 1mu1p topology.
     */
    template<class T>
        bool topological_1e1p_cut(const T & interaction) { return topology<T>(interaction) == "0ph1e0mu0pi1p"; }
    
    /**
     * Apply a 1muNp topological cut. The interaction must have a topology
     * matching 1muNp as defined by the conditions in the count_primaries()
     * function.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction has a 1muNp topology.
     */
    template<class T>
        bool topological_1eNp_cut(const T & interaction)
        {
            std::vector<uint32_t> c(count_primaries(interaction));
            return c[0] == 0 && c[1] == 1 && c[2] == 0 && c[3] == 0 && c[4] >= 1;
        }

    /**
     * Apply a 1muX topological cut. The interaction must have a topology
     * matching 1muX as defined by the conditions in the count_primaries()
     * function.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction has a 1muX topology.
     */
    template<class T>
        bool topological_1eX_cut(const T & interaction)
        {
            std::vector<uint32_t> c(count_primaries(interaction));
            return c[1] == 1;
        }
    
    /**
     * Apply a flash time cut. The interaction must be matched to an in-time
     * flash. The in-time definition is valid for BNB simulation.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction has been matched to an in-time flash.
     */
    template<class T>
        bool flash_cut(const T & interaction)
        {
            if(!valid_flashmatch(interaction))
                return false;
            else
                return (interaction.flash_time >= 0) && (interaction.flash_time <= 9.6);
        }

    /**
     * Apply a flash time cut. The interaction must be matched to an in-time
     * flash. The in-time definition is valid for BNB data.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction has been matched to an in-time flash.
     */
    template<class T>
        bool flash_cut_data(const T & interaction)
        {
            if(!valid_flashmatch(interaction))
                return false;
            else
                return (interaction.flash_time >= 0) && (interaction.flash_time <= 9.6);
        }

    /**
     * Apply a fiducial and containment cut (logical "and" of both).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial and containment cut.
     */
    template<class T>
        bool fiducial_containment_cut(const T & interaction) { return fiducial_cut<T>(interaction) && containment_cut<T>(interaction); }

    /**
     * Apply a fiducial, containment, and topological (1mu1p) cut (logical
     * "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment, and
     * topological cut.
     */
    template<class T>
        bool fiducial_containment_topological_1e1p_cut(const T & interaction) { return fiducial_cut<T>(interaction) && containment_cut<T>(interaction) && topological_1e1p_cut<T>(interaction); }

    /**
     * Apply a fiducial, containment, and topological (1muNp) cut (logical
     * "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment, and
     * topological cut.
     */
    template<class T>
        bool fiducial_containment_topological_1eNp_cut(const T & interaction) { return fiducial_cut<T>(interaction) && containment_cut<T>(interaction) && topological_1eNp_cut<T>(interaction); }

    /**
     * Apply a fiducial, containment, and topological (1muX) cut (logical
     * "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment, and
     * topological cut.
     */
    template<class T>
        bool fiducial_containment_topological_1eX_cut(const T & interaction) { return fiducial_cut<T>(interaction) && containment_cut<T>(interaction) && topological_1eX_cut<T>(interaction); }

    /**
     * Apply a fiducial, containment, topological (1mu1p), and flash time cut
     * (logical "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment,
     * topological, and flash time cut.
     */
    template<class T>
        bool all_1e1p_cut(const T & interaction) { return topological_1e1p_cut<T>(interaction) && flash_cut<T>(interaction) && containment_cut<T>(interaction); } //&& fiducial_cut<T>(interaction) 

    /**
     * Apply a fiducial, containment, topological (1muNp), and flash time cut
     * (logical "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment,
     * topological, and flash time cut.
     */
    template<class T>
        bool all_1eNp_cut(const T & interaction) { return topological_1eNp_cut<T>(interaction) && flash_cut<T>(interaction) && containment_cut<T>(interaction); } //
    
    /**
     * Apply a fiducial, containment, topological (1muX), and flash time cut
     * (logical "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment,
     * topological, and flash time cut.
     */
    template<class T>
        bool all_1eX_cut(const T & interaction) { return topological_1eX_cut<T>(interaction)  && flash_cut<T>(interaction) && containment_cut<T>(interaction); } //&& fiducial_cut<T>(interaction)

    /**
     * Defined the true neutrino interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a neutrino interaction.
    */
    template<class T>
        bool neutrino(const T & interaction) { 
            if(interaction.nu_id >= 0){
                return true; 
            }
            else{
                return false;
            }

        }

    /**
     * Define the true cosmic interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a cosmic.
     */
    template<class T>
        bool cosmic(const T & interaction) { 
                if(interaction.nu_id < 0){
                return true; 
            }
            else{
                return false;
            }

        }

    /**
     * Define the true neutrino interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a neutrino interaction.
    */
    template<class T>
        bool matched_neutrino(const T & interaction) { return interaction.match.size() > 0 && neutrino(interaction); }

    /**
     * Define the true neutrino interaction classification (well-reconstructed).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a well-reconstructed neutrino interaction.
    */
    template<class T>
        bool wellreco_neutrino(const T & interaction) { return wellreco(interaction) && neutrino(interaction); }

    /**
     * Define the true neutrino interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a neutrino interaction.
    */
    template<class T>
        bool matched_cosmic(const T & interaction) { return interaction.match.size() > 0 && cosmic(interaction); }

    /**
     * Define the true 1mu1p interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a 1mu1p neutrino interaction.
     */
    template<class T>
        bool signal_1e1p(const T & interaction) { return topological_1e1p_cut<T>(interaction) && neutrino(interaction); }

    /**
     * Define the true 1muNp interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a 1muNp neutrino interaction.
     */
    template<class T>
        bool signal_1eNp(const T & interaction) { return topological_1eNp_cut<T>(interaction) && neutrino(interaction); }

    /**
     * Define the true 1muNp interaction classification (N > 1 strictly).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a 1muNp (N > 1) neutrino interaction.
     */
    template<class T>
        bool signal_1eNp_Nnot1(const T & interaction) { return !topological_1e1p_cut(interaction) && topological_1eNp_cut<T>(interaction) && neutrino(interaction); }

    /**
     * Define the true 1muX interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a 1muX neutrino interaction.
     */
    template<class T>
        bool signal_1eX(const T & interaction) { return topological_1eX_cut<T>(interaction) && neutrino(interaction); }

    /**
     * Define the true 1muX interaction classification (not 1muNp).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a 1muX (not 1muNp) neutrino interaction.
     */
    template<class T>
        bool signal_1eX_not_1eNp(const T & interaction) { return !topological_1eX_cut(interaction) && !topological_1eNp_cut<T>(interaction) && neutrino(interaction); }

    /**
     * Define the true "other neutrino" interaction classification (1mu1p).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is an "other neutrino" interaction.
     */
    template<class T>
        bool other_nu_1e1p(const T & interaction) { return !topological_1e1p_cut<T>(interaction) && neutrino(interaction); }
    
    /**
     * Define the true "other neutrino" interaction classification (1muNp).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is an "other neutrino" interaction.
     */
    template<class T>
        bool other_nu_1eNp(const T & interaction) { return !topological_1eNp_cut<T>(interaction) && neutrino(interaction); }

    /**
     * Define the true "other neutrino" interaction classification (1muX).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is an "other neutrino" interaction.
     */
    template<class T>
        bool other_nu_1eX(const T & interaction) { return !topological_1eX_cut<T>(interaction) && neutrino(interaction); }

    /**
     * Define true muon particle classification.
     * @tparam T the type of particle (true or reco)
     * @param particle to select on.
     * @return true if the particle is a muon.
    */
    template<class T>
        bool electron(const T & particle) { return particle.pid == 1 && particle.is_contained; }

    /**
     * Define true muon particle classification (matched).
     * @tparam T the type of particle (true or reco)
     * @param particle to select on.
     * @return true if the particle is a muon (matched).
    */
    template<class T>
        bool matched_electron(const T & particle) { return electron(particle) && matched(particle); }

    /**
     * Define true proton particle classification.
     * @tparam T the type of particle (true or reco)
     * @param particle to select on.
     * @return true if the particle is a proton.
    */
    template<class T>
        bool proton(const T & particle) { return particle.pid == 4 && particle.is_contained; }

    /**
     * Define true proton particle classification (matched).
     * @tparam T the type of particle (true or reco)
     * @param particle to select on.
     * @return true if the particle is a proton (matched).
    */
    template<class T>
        bool matched_proton(const T & particle) { return proton(particle) && matched(particle); }

    /**
     * Define cut for particles crossing the cathode.
     * @tparam T the type of particle (true or reco).
     * @param particle to select on.
     * @return true if the particle is cathode-crossing.
    */
    template<class T>
        bool cathode_crossing(const T & particle)
        {
            return particle.is_cathode_crosser;
        }

    /**
     * Define cut for muons crossing the cathode.
     * @tparam T the type of particle (true or reco).
     * @param particle to select on.
     * @return true if the particle is a muon that crosses the cathode.
    */
    template<class T>
        bool cathode_crossing_electron(const T & particle) { return particle.pid == 1 && cathode_crossing(particle); }
        
    /**
     * Define cut for muons not crossing the cathode.
     * @tparam T the type of particle (true or reco).
     * @param particle to select on.
     * @return true if the particle is a muon that does not cross the cathode.
    */
    template<class T>
        bool non_cathode_crossing_electron(const T & particle) { return particle.pid == 1 && !cathode_crossing(particle); }
    
    /**
     * Define muons contained to a single TPC.
     * @tparam T the type of particle (true or reco)
     * @param particle to select on.
     * @return true if the particle is a muon contained to a single TPC.
    */
    template<class T>
        bool contained_tpc_electron(const T & particle) { return electron(particle) && !cathode_crossing(particle); }

    /**
     * Define muons contained that are well-reconstructed.
     * @tparam T the type of particle (true or reco)
     * @param particle to select on.
     * @return true if the particle is a contained muon that is
     * well-reconstructed.
    */
    template<class T>
        bool wellreco_electron(const T & particle) { return electron(particle) && wellreco(particle); }
    /**
     * Define muons contained that are well-reconstructed.
     * @tparam T the type of particle (true or reco)
     * @param particle to select on.
     * @return true if the particle is a contained muon that is
     * well-reconstructed.
    */
    template<class T>
        bool wellreco_proton(const T & particle) { return proton(particle) && wellreco(particle); }
}
#endif