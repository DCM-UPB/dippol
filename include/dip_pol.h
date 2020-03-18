#ifndef DIP_POL_H
#define DIP_POL_H

#include <loos.hpp>

void val_init(sys_info *sys, const std::vector<loos::AtomicGroup> &mol, std::vector<vect_3d> &dip0, std::vector<vect_3d> &dip0_mol, mat_sym_3d *pol0, mat_sym_3d *pol0_mol, std::ofstream &file_traj, std::ofstream &file_dip0, std::ofstream &file_pol0, const int step);
void comp_dip_pol(sys_info *sys, const std::vector<loos::AtomicGroup> &mol, std::vector<vect_3d> &dip0, std::vector<vect_3d> &dip0_mol, vect_3d *dip, vect_3d *dip_mol, mat_sym_3d *pol0, mat_sym_3d *pol0_mol, mat_3d *pol, mat_3d *pol_mol, vect_3d *e_ind, mat_3d *a_ind, vect_3d *v3_tmp, mat_3d *m3_tmp, mat_sym_3d *Tij_mc, mat_sym_3d *Tij_at, double *fthole, std::ofstream &fileo, std::ofstream &file_Tij_mc, std::ofstream &file_Tij_at, std::ofstream &file_dip, std::ofstream &file_pol, std::ofstream &file_dip_ind, std::ofstream &file_pol_ind, int step);
void sumT_dip(mat_sym_3d *tjk, vect_3d &dip, vect_3d &v3_tmp);
void sumT_pol(mat_sym_3d *tjk, mat_3d *pol, mat_3d *m3_tmp);
void init_dip_pol(std::vector<vect_3d> &dip0, vect_3d *dip, mat_sym_3d *pol0, mat_3d *pol, sys_info *sys);
void calc_tij(sys_info *sys, const std::vector<loos::AtomicGroup> &mol, double *fthole, mat_sym_3d *Tij, int size,std::ofstream &file_Tij, int step);
void scf_protection(sys_info *sys, const std::vector<loos::AtomicGroup> &mol, std::vector<vect_3d> &dip0, vect_3d *dip, mat_sym_3d *pol0, mat_sym_3d *pol0_mol, mat_3d *pol, mat_sym_3d *Tij, vect_3d *e_ind, vect_3d *v3_tmp, mat_3d *m3_tmp, int step);
void thole(double *r, double *fthole, int thole);

#endif 
