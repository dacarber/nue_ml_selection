/**
 * @file nue_variables.h
 * @brief Header file for definitions of selection variables in the context
 * of the nue analyses.
 * @author justin.mueller@colostate.edu
*/
#ifndef NUMU_VARIABLES_H
#define NUMU_VARIABLES_H

#include "variables.h"
#include "cuts.h"
#include <algorithm>
#include <cmath>
#include <TVector3.h>

namespace vars
{
    /**
     * Helper methods for calculating particle/interaction-level variables.
    */

    /**
     * Variable for the transverse momentum of a particle.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the transverse momentum of the particle
    */
    
    template<class T>
        double momentum(const T & particle)
        {
            if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
                return std::sqrt(std::pow(particle.momentum[0], 2) + std::pow(particle.momentum[1], 2) + std::pow(particle.momentum[2], 2));
            else
                return std::sqrt(std::pow(particle.momentum[0], 2) + std::pow(particle.momentum[1], 2) + std::pow(particle.momentum[2], 2));
        }
    
    /**
     * Variable for the polar angle (w.r.t the z-axis) of the particle.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the polar angle of the particle.
    */
    template<class T>
        double polar_angle(const T & particle)
        {
            if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
                return std::acos(particle.start_dir[2]);
            else
                return std::acos(particle.start_dir[2]);
        }

    /**
     * Variable for the azimuthal angle (w.r.t the z-axis) of the particle.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the azimuthal angle of the particle.
    */
    template<class T>
        double azimuthal_angle(const T & particle)
        {
            if(particle.start_dir[1] >0);
                return std::acos(particle.start_dir[0] / std::sqrt(std::pow(particle.start_dir[0], 2) + std::pow(particle.start_dir[1], 2)));
            else
                return -std::acos(particle.start_dir[0] / std::sqrt(std::pow(particle.start_dir[0], 2) + std::pow(particle.start_dir[1], 2)));
        }


    /**
     * Variable for the azimuthal angle (w.r.t the z-axis) of the particle.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the azimuthal angle of the particle.
    */
    template<class T>
        double NuMI_angle(const T & particle)
        {   
            double x;
            double y;
            double z;
            double r;
            if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>){
                x = (31512.0380) - particle.start_point[0];
                y = (3364.4912) - particle.start_point[1];
                z = (73363.2532) - particle.start_point[2];
                r = std::sqrt(std::pow(x, 2)+std::pow(y, 2)+std::pow(z, 2));
                x = x/r;
                y = y/r;
                z = z/r;
                return std::acos(x *particle.start_dir[0] + y *particle.start_dir[1]+z *particle.start_dir[2]);
            }
            else{
                x = (31512.0380) - particle.start_point[0];
                y = (3364.4912) - particle.start_point[1];
                z = (73363.2532) - particle.start_point[2];
                r = std::sqrt(std::pow(x, 2)+std::pow(y, 2)+std::pow(z, 2));
                x = x/r;
                y = y/r;
                z = z/r;
                return std::acos(x *particle.start_dir[0] + y *particle.start_dir[1]+z *particle.start_dir[2]);
            }
        }
    /**
     * Variable for the azimuthal angle (w.r.t the z-axis) of the particle.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the azimuthal angle of the particle.
    */
    template<class T>
        double NuMI_nu_theta(const T & interaction)
        {   
            double r;
            std::vector<double> dir_vector(3,0);
               
            dir_vector[0] = interaction.vertex[0] - (31512.0380);
            dir_vector[1] = interaction.vertex[1] - (3364.4912);
            dir_vector[2] = interaction.vertex[2] - (73363.2532);
            r = std::sqrt(std::pow(dir_vector[0], 2)+std::pow(dir_vector[1], 2)+std::pow(dir_vector[2], 2));
            dir_vector[0] = dir_vector[0]/r;
            dir_vector[1] = dir_vector[1]/r;
            dir_vector[2] = dir_vector[2]/r;                                      
            return std::acos(dir_vector[2]);
            
        }
        template<class T>
        double NuMI_nu_phi(const T & interaction)
        {   
            double r;
            std::vector<double> dir_vector(3,0);
               
            dir_vector[0] = interaction.vertex[0] - (31512.0380);
            dir_vector[1] = interaction.vertex[1] - (3364.4912);
            dir_vector[2] = interaction.vertex[2] - (73363.2532);
            r = std::sqrt(std::pow(dir_vector[0], 2)+std::pow(dir_vector[1], 2)+std::pow(dir_vector[2], 2));
            dir_vector[0] = dir_vector[0]/r;
            dir_vector[1] = dir_vector[1]/r;
            dir_vector[2] = dir_vector[2]/r;                                      
            if(dir_vector[1] >0);
                return std::acos(dir_vector[0] / std::sqrt(std::pow(dir_vector[0], 2) + std::pow(dir_vector[1], 2)));
            else
                return -std::acos(dir_vector[0] / std::sqrt(std::pow(dir_vector[0], 2) + std::pow(dir_vector[1], 2)));
            
        }
    /**
     * Methods for calculating the reconstructed variables for the numu analyses.
    */

    /**
     * Variable for total visible energy of interaction.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the total visible energy of the interaction.
    */
    /*template<class T>
        double visible_energy(const T & interaction)
        {
            double energy(0);
            for(const auto & p : interaction.particles)
            {
                if(p.is_primary)
                {
                    if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                    {
                        energy += p.energy_deposit;
                    }
                    else
                    {
                        if(p.pid < 2) energy += p.calo_ke;
                        else energy += p.csda_ke;
                    }
                    if(p.pid == 2) energy += MUON_MASS;
                    else if(p.pid == 3) energy += PION_MASS;
                }
            }
            return energy;
        }
    */
    /**
     * Variable for finding the leading muon kinetic energy.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the kinetic energy of the leading muon.
    */
    template<class T>
        double leading_muon_ke(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 2));
            double energy(csda_ke(interaction.particles[i]));
            if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                energy = ke_init(interaction.particles[i]);
            return energy;
        }
    /**
     * Variable for finding the leading muon kinetic energy.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the kinetic energy of the leading muon.
    */
    template<class T>
        int non_primary_number(const T & interaction)
        {
            int num_par(0);
            for(const auto & p : interaction.particles)
            {   double energy(p.pid > 1 ? p.csda_ke : p.calo_ke);
                if(p.is_primary == false)
                {
                    num_par+=1;                    
                }
                else
                {
                    if (p.pid == 0 and energy <= 70) num_par+=1;
                    else if (p.pid == 1 and energy <= 25) num_par+=1;
                    else if (p.pid == 2 and energy <= 25) num_par+=1;
                    else if (p.pid == 3 and energy <= 25) num_par+=1;
                    else if (p.pid == 4 and energy <= 40) num_par+=1;
                }
            }
            return num_par;
        }
    /**
     * Variable for finding the leading proton kinetic energy.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the kinetic energy of the leading muon.
    */
    /*
    template<class T>
        double leading_proton_ke(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 4));
            double energy(csda_ke(interaction.particles[i]));
            if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                energy = ke_init(interaction.particles[i]);
            return energy;
        }
    */
    /**
     * Variable for the transverse momentum of the leading muon.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the transverse momentum of the leading muon.
    */
    /*
    template<class T>
        double leading_electron_pt(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 1));
            return transverse_momentum(interaction.particles[i]);
        }
    */
    /**
     * Variable for the transverse momentum of the leading proton.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the transverse momentum of the leading proton.
    */
    
    template<class T>
        double leading_proton_p(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 4));
            return momentum(interaction.particles[i]);
        }
    /**
     * Variable for the transverse momentum of the leading proton.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the transverse momentum of the leading proton.
    */
    
    template<class T>
        double true_leading_proton_p(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 4));
            return momentum(interaction.truth_particles[i]);
        }
    
    /**
     * Variable for the muon polar angle.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the polar angle of the leading muon.
    */
    template<class T>
        double electron_polar_angle(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 1));
            return polar_angle(interaction.particles[i]);
        }

    /**
     * Variable for the muon azimuthal angle.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the azimuthal angle of the leading muon.
     */
    template<class T>
        double electron_azimuthal_angle(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 1));
            return azimuthal_angle(interaction.particles[i]);
        }
    /**
     * Variable for the muon azimuthal angle.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the azimuthal angle of the leading muon.
     */
    template<class T>
        double electron_NuMI_angle(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 1));
            return NuMI_angle(interaction.particles[i]);
        }
    template<class T>
        double proton_polar_angle(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 4));
            return polar_angle(interaction.particles[i]);
        }

    /**
     * Variable for the muon azimuthal angle.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the azimuthal angle of the leading muon.
     */
    template<class T>
        double proton_azimuthal_angle(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 4));
            return azimuthal_angle(interaction.particles[i]);
        }
    /**
     * Variable for the opening angle between leading muon and
     * proton.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the opening angle between the leading muon and
     * proton.
    */
    template<class T>
        double opening_angle(const T & interaction)
        {
            auto & e(interaction.particles[leading_particle_index(interaction, 1)]);
            auto & p(interaction.particles[leading_particle_index(interaction, 4)]);
            if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                return std::acos(e.start_dir[0] * p.start_dir[0] + e.start_dir[1] * p.start_dir[1] + e.start_dir[2] * p.start_dir[2]);
            else
                return std::acos(e.start_dir[0] * p.start_dir[0] + e.start_dir[1] * p.start_dir[1] + e.start_dir[2] * p.start_dir[2]);
        }


    template <class T> 
        TVector3  electron_transverse_momentum(const T & interaction)
      {
        //TVector3 beamdir(0, 0, 1); // BNB
        std::vector<double> dir_vector(3,0);                    
        dir_vector[0] = interaction.vertex[0] - (31512.0380);
        dir_vector[1] = interaction.vertex[1] - (3364.4912);
        dir_vector[2] = interaction.vertex[2] - (73363.2532);
        r = std::sqrt(std::pow(dir_vector[0], 2)+std::pow(dir_vector[1], 2)+std::pow(dir_vector[2], 2));
        dir_vector[0] = dir_vector[0]/r;
        dir_vector[1] = dir_vector[1]/r;
        dir_vector[2] = dir_vector[2]/r;                                                                                                                                                                                              
        TVector3 beamdir(dir_vector[0], dir_vector[1], dir_vector[2]); // NuMI                                                                                                                                                                                                 

        // Loop over particles                                                                                                                                                                                                                                        
        
        size_t i(leading_particle_index(interaction, 1));
        auto & particle = interaction.particles[i];
            //if(!particle.is_primary) continue;

            // pT = p - pL                                                                                                                                                                                                                                            
            //    = p-(p dot beamdir) * beamdir                                                                                                                                                                                                                         
            TVector3 p;
            TVector3 pL;
            TVector3 pT;

            

            p.SetX(part.momentum[0]);
            p.SetY(part.momentum[1]);
            p.SetZ(part.momentum[2]);

            pL = p.Dot(beamdir) * beamdir;
            pT = p - pL;

            pT0 += pT[0];
            pT1 += pT[1];
            pT2 += pT[2];

                                                                                                                                                                                                                                          
        
        return pT;
      }

    template <class T> 
        TVector3 proton_transverse_momentum(const T & interaction)
      {
        //TVector3 beamdir(0, 0, 1); // BNB
        std::vector<double> dir_vector(3,0);                    
        dir_vector[0] = interaction.vertex[0] - (31512.0380);
        dir_vector[1] = interaction.vertex[1] - (3364.4912);
        dir_vector[2] = interaction.vertex[2] - (73363.2532);
        r = std::sqrt(std::pow(dir_vector[0], 2)+std::pow(dir_vector[1], 2)+std::pow(dir_vector[2], 2));
        dir_vector[0] = dir_vector[0]/r;
        dir_vector[1] = dir_vector[1]/r;
        dir_vector[2] = dir_vector[2]/r;                                                                                                                                                                                              
        TVector3 beamdir(dir_vector[0], dir_vector[1], dir_vector[2]); // NuMI                                                                                                                                                                                                 

        // Output                                                                                                                                                                                                                                                     
        double pT0(0);

        // Loop over particles                                                                                                                                                                                                                                        
        for(auto & part : interaction.particles)
          {

            if(!part.is_primary or part.pid != 4 ) continue;

            // pT = p - pL                                                                                                                                                                                                                                            
            //    = p-(p dot beamdir) * beamdir                                                                                                                                                                                                                         
            TVector3 p;
            TVector3 pL;
            TVector3 pT;

            

            p.SetX(part.momentum[0]);
            p.SetY(part.momentum[1]);
            p.SetZ(part.momentum[2]);

            pL = p.Dot(beamdir) * beamdir;
            pT = p - pL;
            pT0 += pT[coord];
                                                                                                                                                                                                                                          
        }
        return pT;
      }


    template<class T>
        double delta_pt(const T & interaction)
        {
            TVector3 plT(electron_transverse_momentum(interaction));
            TVector3 ppT(proton_transverse_momentum(interaction));
            TVector3 delta_p = plT+ppT;
            return delta_p.Mag;
        }
    template<class T>
        double delta_alphaT(const T & interaction)
        {
            TVector3 plT(electron_transverse_momentum(interaction));
            TVector3 ppT(proton_transverse_momentum(interaction));
            TVector3 delta_p = plT+ppT;
            double delta_a = std::acos(-plT.Dot(delta_p)/(plT.Mag * delta_p.Mag));
            return delta_a;
        }
    template<class T>
        double delta_phiT(const T & interaction)
        {
            TVector3 plT(electron_transverse_momentum(interaction));
            TVector3 ppT(proton_transverse_momentum(interaction));
            double delta_phi = std::acos(-plT.Dot(ppT)/(plT.Mag * ppT.Mag));
            return delta_phi;
        }

    /**
     * Variable for phi_T of the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the phi_T of the interaction.
    */
    template<class T>
        double phiT(const T & interaction)
        {
            double lpx(0), lpy(0), hpx(0), hpy(0);
            for(const auto & p : interaction.particles)
                if(cuts::final_state_signal(p))
                {
                    if(p.pid > 2)
                    {
                        if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                        {
                            std::vector<double> dir_vector(3,0);
                            std::vector<double> z_axis = {0,0,1};
                            dir_vector[0] = particle.start_point[0] - (31512.0380);
                            dir_vector[1] = particle.start_point[1] - (3364.4912);
                            dir_vector[2] = particle.start_point[2] - (73363.2532);
                            r = std::sqrt(std::pow(dir_vector[0], 2)+std::pow(dir_vector[1], 2)+std::pow(dir_vector[2], 2));
                            dir_vector[0] = dir_vector[0]/r;
                            dir_vector[1] = dir_vector[1]/r;
                            dir_vector[2] = dir_vector[2]/r;

                            std::vector<double> cross_product = {-dir_vector[1],-dir_vector[0],0}


                            hpx += p.momentum[0];
                            hpy += p.momentum[1];
                        }
                        else
                        {
                            hpx += p.momentum[0];
                            hpy += p.momentum[1];
                        }
                    }
                    else if(p.pid == 2)
                    {
                        if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                        {
                            lpx += p.momentum[0];
                            lpy += p.momentum[1];
                        }
                        else
                        {
                            lpx += p.momentum[0];
                            lpy += p.momentum[1];
                        }
                    }
                }
            return std::acos((-hpx * lpx - hpy * lpy) / (std::sqrt(std::pow(hpx, 2) + std::pow(hpy, 2)) * std::sqrt(std::pow(lpx, 2) + std::pow(lpy, 2))));
        }

    /**
     * Variable for alpha_T of the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the phi_T of the interaction.
    */
    template<class T>
        double alphaT(const T & interaction)
        {
            double lpx(0), lpy(0), px(0), py(0);
            for(const auto & p : interaction.particles)
                if(cuts::final_state_signal(p))
                {
                    if(p.pid == 1)
                    {
                        if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                        {
                            x = (31512.0380) - particle.start_point[0];
                            y = (3364.4912) - particle.start_point[1];
                            z = (73363.2532) - particle.start_point[2];
                            lpx += p.momentum[0];
                            lpy += p.momentum[1];
                        }
                        else
                        {
                            lpx += p.momentum[0];
                            lpy += p.momentum[1];
                        }
                    }
                    if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                    {
                        px += p.momentum[0];
                        py += p.momentum[1];
                    }
                    else
                    {
                        px += p.momentum[0];
                        py += p.momentum[1];
                    }
                }
            return std::acos((-px * lpx - py * lpy) / (std::sqrt(std::pow(px, 2) + std::pow(py, 2)) * std::sqrt(std::pow(lpx, 2) + std::pow(lpy, 2))));
        }

    /**
     * Variable for the muon softmax score for the leading muon of the
     * interaction.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the muon softmax score of the leading muon.
    */
    
    template<class T>
        double electron_softmax(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 1));
            return interaction.particles[i].pid_scores[1];
        }
    
    /**
     * Variable for the proton softmax score for the leading proton of the
     * interaction.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the proton softmax score of the leading proton.
    */
    
    template<class T>
        double proton_softmax(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 4));
            return interaction.particles[i].pid_scores[4];
        }
    
}
#endif