# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp -p atom_vec_meso.cpp ..
  cp -p pair_sph_heatconduction.cpp ..
  cp -p pair_sph_idealgas.cpp ..
  cp -p pair_sph_lj.cpp ..
  cp -p pair_sph_rhosum.cpp ..
  cp -p pair_sph_taitwater.cpp ..
  cp -p pair_sph_taitwater_morris.cpp ..
  cp -p pair_sph_protondiffusion.cpp ..
  cp -p pair_sph_nafion.cpp ..
  cp -p pair_sph_elecpot.cpp ..
  cp -p compute_meso_e_atom.cpp ..
  cp -p compute_meso_rho_atom.cpp ..
  cp -p compute_meso_t_atom.cpp ..
  cp -p compute_meso_prtn_conc_atom.cpp ..
  cp -p compute_meso_elec_pot_atom.cpp ..
  cp -p fix_meso.cpp ..
  cp -p fix_meso_stationary.cpp ..
  cp -p fix_electrostatic.cpp ..

  cp -p atom_vec_meso.h ..
  cp -p pair_sph_heatconduction.h ..
  cp -p pair_sph_idealgas.h ..
  cp -p pair_sph_lj.h ..
  cp -p pair_sph_rhosum.h ..
  cp -p pair_sph_taitwater.h ..
  cp -p pair_sph_taitwater_morris.h ..
  cp -p pair_sph_protondiffusion.h ..
  cp -p pair_sph_nafion.h ..
  cp -p pair_sph_elecpot.h ..
  cp -p compute_meso_e_atom.h ..
  cp -p compute_meso_rho_atom.h ..
  cp -p compute_meso_t_atom.h ..
  cp -p compute_meso_prtn_conc_atom.h ..
  cp -p compute_meso_elec_pot_atom.h ..
  cp -p fix_meso.h ..
  cp -p fix_meso_stationary.h ..
  cp -p fix_electrostatic.h ..

elif (test $1 = 0) then

  rm -f ../atom_vec_meso.cpp
  rm -f ../pair_sph_heatconduction.cpp
  rm -f ../pair_sph_idealgas.cpp
  rm -f ../pair_sph_lj.cpp
  rm -f ../pair_sph_rhosum.cpp
  rm -f ../pair_sph_taitwater.cpp
  rm -f ../pair_sph_taitwater_morris.cpp
  rm -f ../pair_sph_protondiffusion.cpp
  rm -f ../pair_sph_nafion.cpp
  rm -f ../pair_sph_elecpot.cpp
  rm -f ../compute_meso_e_atom.cpp
  rm -f ../compute_meso_rho_atom.cpp
  rm -f ../compute_meso_t_atom.cpp
  rm -f ../compute_meso_prtn_conc_atom.cpp
  rm -f ../compute_meso_elec_pot_atom.cpp
  rm -f ../fix_meso.cpp
  rm -f ../fix_meso_stationary.cpp
  rm -f ../fix_electrostatic.cpp

  rm -f ../atom_vec_meso.h
  rm -f ../pair_sph_heatconduction.h
  rm -f ../pair_sph_idealgas.h
  rm -f ../pair_sph_lj.h
  rm -f ../pair_sph_rhosum.h
  rm -f ../pair_sph_taitwater.h
  rm -f ../pair_sph_taitwater_morris.h
  rm -f ../pair_sph_protondiffusion.h
  rm -f ../pair_sph_nafion.h
  rm -f ../pair_sph_elecpot.h
  rm -f ../compute_meso_e_atom.h
  rm -f ../compute_meso_rho_atom.h
  rm -f ../compute_meso_t_atom.h
  rm -f ../compute_meso_prtn_conc_atom.h
  rm -f ../compute_meso_elec_pot_atom.h
  rm -f ../fix_meso.h
  rm -f ../fix_meso_stationary.h
  rm -f ../fix_electrostatic.h  

fi
