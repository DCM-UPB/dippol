/*
  In this file can be found the subroutines
  which are necessary to calculate the dipole moment
  and the polarizability tensor

*/

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <cstdlib>
#include <math.h>
#include "main.h"
#include "dip_pol.h"
#include "mol_atom.h"
#include "my_math.h"
#include "write.h"

using std::cout;
using std::endl;
using std::setw;

void val_init(sys_info *sys, const std::vector<loos::AtomicGroup> &mol, std::vector<vect_3d> &dip0, std::vector<vect_3d> &dip0_mol, mat_sym_3d *pol0, mat_sym_3d *pol0_mol, std::ofstream &file_traj, std::ofstream &file_dip0, std::ofstream &file_pol0, const int step){
  //int i=0;
  double v_oh1[3]={0.}, v_oh2[3]={0.};
  loos::GCoord voh1, voh2;
    
  #pragma omp parallel for
  for(uint mol_i=0 ; mol_i < mol.size(); mol_i++){
      #ifndef YUKI
      const loos::GCoord com = mol[mol_i].centerOfMass();
      file_traj << "X " << std::setw(12) << com.x()  << std::setw(12) << com.y()  << std::setw(12) << com.z() << endl;
      #endif
      
      //NOTE new temporary code, 
      // Should be replaced with proper calculation of dipole and pol
      // from ref molecule
      voh1 = mol[mol_i][1]->coords() - mol[mol_i][0]->coords();
      voh2 = mol[mol_i][2]->coords() - mol[mol_i][0]->coords();
      voh1 /= voh1.length();
      voh2 /= voh2.length();
      
      v_oh1[0] = voh1[0];
      v_oh1[1] = voh1[1];
      v_oh1[2] = voh1[2];
      v_oh2[0] = voh2[0];
      v_oh2[1] = voh2[1];
      v_oh2[2] = voh2[2];
      
      /*==============================================
       * Note:
       * typ_dip 1 = atomic dipole moment
       * typ_pol 1 = atomic polarizabilit
       * 
       * The atomic polarizability can be used even
       * without atomic dipole: 1 on 3 point independents
       * 
       * Therefore there are 4 possibilities:
       * at_dip+ mol_pol is not possible,
       * all three other combinations are possible, i.e.:
       * at_dip + at_pol, mol_dip + at_pol, mol_dip + mol_pol
       = ============*===================================*/
      //NOTE 
      //DANGER omp not tested with anything besides mol/mol!!!
      //NOTE
      if((*sys).typ_dip){ /*Atomic dipole moment*/
          if((*sys).typ_pol){ /* Atomic polarizability */
              dippol0_atat(sys,v_oh1,v_oh2,&(dip0[3*mol_i]),&(pol0[3*mol_i]));
              dip_at2mol(&dip0[3*mol_i],&dip0_mol[mol_i]);
              pol0_at2mol(&(pol0[3*mol_i]),&(pol0_mol[mol_i]));
          }
          else{ /* Molecular polarizability */
              std::cout << "The combination of atomic dipole with molecular polarizability is not possible!" << std::endl;
              std::cout << "Please modify the input accordingly." << std::endl;
              std::cout << "End of program" << std::endl;
              exit(0);
          }
      }
      else{/*Molecular dipole moment */
          if((*sys).typ_pol){ /* Atomic polarizability */
              dippol0_molat(sys,v_oh1,v_oh2,&dip0[mol_i],&(pol0[3*mol_i]));
              dip0_mol[mol_i]=dip0[mol_i];
              pol0_at2mol(&(pol0[3*mol_i]),&(pol0_mol[mol_i]));
          }
          else{ /* Molecular polarizability */ //NOTE I should only focus on this for now
              //std::cout << "Parsing molecule " << mol_i << std::endl;
              dippol0_molmol(sys, mol[mol_i], dip0[mol_i], &(pol0[mol_i]));
              //cout << dip0[mol_i] << endl;
              dip0_mol[mol_i]=dip0[mol_i];
              pol0_mol[mol_i]=pol0[mol_i];
          }
      }
      
  }
  
  
  #ifdef DEBUG
  for(int i=0 ; i<(*sys).nb_dip ; i++){
      file_dip0 <<  setw(10) << step << setw(5) << i << setw(18) << dip0[i].x() 
      << setw(18) << dip0[i].y()  << setw(18) << dip0[i].z()  << endl;
      //fprintf(file_dip0,"%10d %5d %9.4f %9.4f %9.4f\n",step,i,dip0[i].x(),dip0[i].y(),dip0[i].z());
  }
  for(int i=0 ; i<(*sys).nb_pol ; i++){
      file_pol0 << setw(10) << step << setw(5) << i 
      << setw(18) << pol0[i].xx
      << setw(18) << pol0[i].yx
      << setw(18) << pol0[i].yy
      << setw(18) << pol0[i].zx
      << setw(18) << pol0[i].zy
      << setw(18) << pol0[i].zz << endl;
      //fprintf(file_pol0,"%10d %5d %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n",step,i,pol0[i].xx,pol0[i].yx,pol0[i].yy,pol0[i].zx,pol0[i].zy,pol0[i].zz);
  }
  #endif
  
  return;
}



/*
  Complete dipole (permanent + induced) and effective polarizbility
  by Morita and Hynes, J. Phys. Chem. B 2002, 106, 673-685
  and calculated only once (Yuki version)
*/
void comp_dip_pol(sys_info *sys, const std::vector<loos::AtomicGroup> &mol, std::vector<vect_3d> &dip0, std::vector<vect_3d> &dip0_mol, vect_3d *dip, vect_3d *dip_mol, mat_sym_3d *pol0, mat_sym_3d *pol0_mol, mat_3d *pol, mat_3d *pol_mol, vect_3d *e_ind, mat_3d *a_ind, vect_3d *v3_tmp, mat_3d *m3_tmp, mat_sym_3d *Tij_mc, mat_sym_3d *Tij_at, double *fthole, std::ofstream &fileo, std::ofstream &file_Tij_mc, std::ofstream &file_Tij_at, std::ofstream &file_dip, std::ofstream &file_pol, std::ofstream &file_dip_ind, std::ofstream &file_pol_ind, std::ofstream &file_dip_f, std::ofstream &file_pol_f, int step){
  int test_scf=0;
  
  /* Calculation of the dipole field tensor Tij*/
  calc_tij(sys,mol,fthole,Tij_mc,(*sys).nb_mol,file_Tij_mc,step);/*Molecular*/
  if((*sys).typ_dip+(*sys).typ_pol){/*Atomic*/
    calc_tij(sys,mol,fthole,Tij_at,(*sys).nb_point,file_Tij_at,step);
  }
  
  /* Initialization of the total terms*/
  init_dip_pol(dip0,dip,pol0,pol,sys);
  /* Initialization of the induced terms at the order 0*/
  init_dip_pol(dip0,e_ind,pol0,a_ind,sys);
  
  /*========================*/
  /*========SCF loop========*/
  /*========================*/
  test_scf=1;/*To know when we go in or out of the while loop*/
  while(test_scf){
      /*Initialization*/
      test_scf=0;
      #pragma omp parallel for
      for(int i=0;i<(*sys).nb_dip;i++){
          v3_tmp[i].x()=0.; v3_tmp[i].y()=0.; v3_tmp[i].z()=0.;
      }
      
      #pragma omp parallel for
      for(int i=0;i<(*sys).nb_pol;i++){
          m3_tmp[i].xx=0.; m3_tmp[i].xy=0.; m3_tmp[i].xz=0.;
          m3_tmp[i].yx=0.; m3_tmp[i].yy=0.; m3_tmp[i].yz=0.;
          m3_tmp[i].zx=0.; m3_tmp[i].zy=0.; m3_tmp[i].zz=0.;
      }
      
      
      /*-Sum Tij.e_ind[j], -Sum Tij.a_ind[j]*/
      /*-Sum Tji.e_ind[i], -Sum Tji.a_ind[i]*/
      /* Dipole*/
      if((*sys).typ_dip){/*Atomic*/
          int ij=0;
          for(int i=0;i<(*sys).nb_dip;i++){
              for(int j=i+1;j<(*sys).nb_dip;j++){
                  sumT_dip(&(Tij_at[ij]),e_ind[j],v3_tmp[i]);
                  sumT_dip(&(Tij_at[ij]),e_ind[i],v3_tmp[j]);
                  ij++;
              }
          }
      }
      else{/*Molecular*/
          //#pragma omp parallel for
          for(int i=0;i<(*sys).nb_dip;i++){
              int factor =0;
              for (int f_i=0; f_i < i; ++f_i)
                  factor += f_i+1;
              for(int j=i+1;j<(*sys).nb_dip;j++){
                  int ij = i*(*sys).nb_dip - factor + j -i-1;
                  //cout << std::setw(5) << i 
                  //<< std::setw(5) << j
                  //<< std::setw(5) << ij << endl;
                  //ij = index;
                  sumT_dip(&(Tij_mc[ij]),e_ind[j],v3_tmp[i]);
                  sumT_dip(&(Tij_mc[ij]),e_ind[i],v3_tmp[j]);
                  //cout << "Checking " << i << endl;
                  //cout << i << " " << v3_tmp[i] << "  " << j << "  " << v3_tmp[i] << endl;
                  //ij++;
              }
          }
      }

      /*Polarizability*/
      if((*sys).typ_pol){/*Atomic*/
          int ij=0;
          for(int i=0;i<(*sys).nb_pol;i++){
              int factor =0;
              for (int f_i=0; f_i < i; ++f_i)
                  factor += f_i+1;        
              for(int j=i+1;j<(*sys).nb_pol;j++){
                  ij = i*(*sys).nb_pol - factor + j -i-1;
                  sumT_pol(&(Tij_at[ij]),&(a_ind[j]),&(m3_tmp[i]));
                  sumT_pol(&(Tij_at[ij]),&(a_ind[i]),&(m3_tmp[j]));
                  //ij++;
              }
          }
      }
      else{/*Molecular*/
          for(int i=0;i<(*sys).nb_pol;i++){
              int factor =0;
              for (int f_i=0; f_i < i; ++f_i)
                  factor += f_i+1;        
              for(int j=i+1;j<(*sys).nb_pol;j++){
                  int ij = i*(*sys).nb_pol - factor + j -i-1;
                  sumT_pol(&(Tij_mc[ij]),&(a_ind[j]),&(m3_tmp[i]));
                  sumT_pol(&(Tij_mc[ij]),&(a_ind[i]),&(m3_tmp[j]));
                  //ij++;
              }
          }
      }
      
      /*==========================================================*/
      /* Calculation of the induced dipole moment at the order N+1*/
      /*==========================================================*/
      if((*sys).typ_dip){/*Atomic*/
          for(int i=0;i<(*sys).nb_dip;i++){
              /*If atomic dipole, therefore there are atomic polarizabilities
               So it is sure t*hat pol0[i] is define*d*/
              mult_msym_v_3d(&(pol0[i]),v3_tmp[i],e_ind[i]);
          }
      }
      else{/*Molecular*/
          for(int i=0;i<(*sys).nb_dip;i++){
              //cout << m_to_eigen(pol0_mol[i]) << endl;
              //cout << v3_tmp[i] << endl;
              //cout << e_ind[i] << endl;
              mult_msym_v_3d(&(pol0_mol[i]),v3_tmp[i],e_ind[i]);
              //cout << v3_tmp[i] << "   " << e_ind[i] << endl;
          }
      }
      
      for(int i=0;i<(*sys).nb_dip;i++){
          /*Update of the complete dipole at the order N+1*/
          dip[i].x() += e_ind[i].x();
          dip[i].y() += e_ind[i].y();
          dip[i].z() += e_ind[i].z();
          
          /*The norm of the induced dipole is not converged*/
          if(sqrt(pow(e_ind[i].x(),2)+pow(e_ind[i].y(),2)+pow(e_ind[i].z(),2))>DIP_CONV){
              test_scf=1;
          }
          
          #ifdef DEBUG
          file_dip_ind << setw(10) << step << setw(5) << i
          << setw(18) << e_ind[i].x()
          << setw(18) << e_ind[i].y()
          << setw(18) << e_ind[i].z() << endl;
          
          file_dip << setw(10) << step << setw(5) << i 
          << setw(18) << dip[i].x()
          << setw(18) << dip[i].y()
          << setw(18) << dip[i].z() << endl;
          #endif
          
      }
      
      /*===========================================================*/
      /* Calculation of the induced polarizability at the order N+1*/
      /*===========================================================*/
      for(int i=0;i<(*sys).nb_pol;i++){
          mult_msym_m_3d(&(pol0[i]),&(m3_tmp[i]),&(a_ind[i]));
          
          /*Update of the complete polarizability at the order N+1*/
          pol[i].xx += a_ind[i].xx ; pol[i].xy += a_ind[i].xy ; pol[i].xz += a_ind[i].xz ;
          pol[i].yx += a_ind[i].yx ; pol[i].yy += a_ind[i].yy ; pol[i].yz += a_ind[i].yz ;
          pol[i].zx += a_ind[i].zx ; pol[i].zy += a_ind[i].zy ; pol[i].zz += a_ind[i].zz ;
          
          /*The trace of the induced polarizability is not converged*/
          if(sqrt(pow(a_ind[i].xx+a_ind[i].yy+a_ind[i].zz,2))>POL_CONV ){
              test_scf=1;
          }
          
          #ifdef DEBUG
          file_pol_ind << setw(10) << step << setw(5) << i
          << setw(18) << a_ind[i].xx<< setw(18) << a_ind[i].xy<< setw(18) << a_ind[i].xz
          << setw(18) << a_ind[i].yx<< setw(18) << a_ind[i].yy<< setw(18) << a_ind[i].yz
          << setw(18) << a_ind[i].zx<< setw(18) << a_ind[i].zy<< setw(18) << a_ind[i].zz << endl;

          file_pol << setw(10) << step << setw(5) << i
          << setw(18) << pol[i].xx<< setw(18) << pol[i].xy<< setw(18) << pol[i].xz
          << setw(18) << pol[i].yx<< setw(18) << pol[i].yy<< setw(18) << pol[i].yz
          << setw(18) << pol[i].zx<< setw(18) << pol[i].zy<< setw(18) << pol[i].zz << endl;
          #endif
          
          /*Protection against infinite loop
           * If we are out of the limits, only the first order dipole
           * moment and polarizability are written for the whole step*/
          if(pol[i].xx>POL_LIM ||pol[i].xy>POL_LIM ||pol[i].xz>POL_LIM || \
              pol[i].yx>POL_LIM ||pol[i].yy>POL_LIM ||pol[i].yz>POL_LIM ||	\
              pol[i].zx>POL_LIM ||pol[i].zy>POL_LIM ||pol[i].zz>POL_LIM){
              scf_protection(sys,mol,dip0,dip,pol0,pol0_mol,pol,Tij_mc,e_ind,v3_tmp,m3_tmp,step);
          i=(*sys).nb_mol; /*Out of the for loop*/
          test_scf=0; /*Out of the while loop*/
              }
              
      }/* Induced polarizability*/
  }/*SCF*/
  
  //cout << "done molecular5" << endl;
  /*===================================
   *  Calculation of the molecular values
   * ===================================*/
  /* Dipole moment */
  if((*sys).typ_dip){
      for(int i=0;i<(*sys).nb_mol;i++){
          dip_at2mol(&(dip[3*i]),&(dip_mol[i]));
      }
  }
  else{
      dip_mol=dip;
  }
  
  /* Polarizability */
  if((*sys).typ_pol){ 
      for(int i=0;i<(*sys).nb_mol;i++){
          pol_at2mol(&(pol[3*i]),&(pol_mol[i]));
      }
  }
  else{
      pol_mol=pol;
  }
  
  for(int i=0;i<(*sys).nb_dip;i++){
      file_dip_f << setw(10) << step << setw(5) << i 
      << setw(18) << dip[i].x()
      << setw(18) << dip[i].y()
      << setw(18) << dip[i].z() << endl;
  }
  
  for(int i=0;i<(*sys).nb_dip;i++){
      file_pol_f << setw(10) << step << setw(5) << i
      << setw(18) << pol[i].xx<< setw(18) << pol[i].xy<< setw(18) << pol[i].xz
      << setw(18) << pol[i].yx<< setw(18) << pol[i].yy<< setw(18) << pol[i].yz
      << setw(18) << pol[i].zx<< setw(18) << pol[i].zy<< setw(18) << pol[i].zz << endl;
  }
      
  /* We write the info about the dipole moment and the polarizability */
  write_dippol(*sys,dip0_mol,dip_mol,pol0_mol,pol_mol,fileo);
  //cout << "done write" << endl;
  return;
}



/* Calculation of the Sum involving Tjk and dip[k]
   Preliminary step before the calculation of the induced dipole

   Remark1: The sign "-" is already included 
   Remark2: tjk is diagonal (only lower part)

*/
void sumT_dip(mat_sym_3d *tjk, vect_3d &dip, vect_3d &v3_tmp){
  
  /* Contribution of k on the induced field of j */
  v3_tmp.x() -= (*tjk).xx*dip.x() + (*tjk).yx*dip.y() + (*tjk).zx*dip.z();
  v3_tmp.y() -= (*tjk).yx*dip.x() + (*tjk).yy*dip.y() + (*tjk).zy*dip.z();
  v3_tmp.z() -= (*tjk).zx*dip.x() + (*tjk).zy*dip.y() + (*tjk).zz*dip.z();

  return;
}




/* Calculation of the Sum involving Tjk and pol[k]
   Preliminary step before the calculation of the induced polarizability 

   Remark1: The sign "-" is already included 
   Remark2: tjk is diagonal (only lower part)

*/
void sumT_pol(mat_sym_3d *tjk, mat_3d *pol, mat_3d *m3_tmp){
  /*Diagonal terms of j*/
  (*m3_tmp).xx -= (*tjk).xx*(*pol).xx + (*tjk).yx*(*pol).yx + (*tjk).zx*(*pol).zx;
  (*m3_tmp).yy -= (*tjk).yx*(*pol).xy + (*tjk).yy*(*pol).yy + (*tjk).zy*(*pol).zy;
  (*m3_tmp).zz -= (*tjk).zx*(*pol).xz + (*tjk).zy*(*pol).yz + (*tjk).zz*(*pol).zz;

  /*Lower part of j*/
  (*m3_tmp).yx -= (*tjk).yx*(*pol).xx + (*tjk).yy*(*pol).yx + (*tjk).zy*(*pol).zx;
  (*m3_tmp).zx -= (*tjk).zx*(*pol).xx + (*tjk).zy*(*pol).yx + (*tjk).zz*(*pol).zx;
  (*m3_tmp).zy -= (*tjk).zx*(*pol).xy + (*tjk).zy*(*pol).yy + (*tjk).zz*(*pol).zy;

  /*Upper part of j*/
  (*m3_tmp).xy -= (*tjk).xx*(*pol).xy + (*tjk).yx*(*pol).yy + (*tjk).zx*(*pol).zy;
  (*m3_tmp).xz -= (*tjk).xx*(*pol).xz + (*tjk).yx*(*pol).yz + (*tjk).zx*(*pol).zz;
  (*m3_tmp).yz -= (*tjk).yx*(*pol).xz + (*tjk).yy*(*pol).yz + (*tjk).zy*(*pol).zz;

  return;
}

/*Initialization before the calculation of 
  the dipole moment and the polarizability*/
void init_dip_pol(std::vector<vect_3d> &dip0, vect_3d *dip, mat_sym_3d *pol0, mat_3d *pol, sys_info *sys){
    int i=0;
    //cout << "initializing dipoles" << endl;
    #pragma omp parallel for
    for(int i=0;i<(*sys).nb_dip;i++){
        dip[i].x()=dip0[i].x();
        dip[i].y()=dip0[i].y();
        dip[i].z()=dip0[i].z();
    }
    #pragma omp parallel for
    for(int i=0;i<(*sys).nb_pol;i++){
        pol[i].xx=pol0[i].xx; pol[i].xy=pol0[i].yx; pol[i].xz=pol0[i].zx;
        pol[i].yx=pol0[i].yx; pol[i].yy=pol0[i].yy; pol[i].yz=pol0[i].zy;
        pol[i].zx=pol0[i].zx; pol[i].zy=pol0[i].zy; pol[i].zz=pol0[i].zz;
    }
    
    return;
}


/* Calculation of the dipole field tensor Tij
 * Only for k>j because Tjk=Tkj and Tjj=0 */
void calc_tij (sys_info* sys, const std::vector<loos::AtomicGroup> &mol, double* fthole, mat_sym_3d* Tij, int size, std::ofstream &file_Tij, int step )
{
    int i=0, j=0;
    double r3=0., r5=0.;
    
//     /*Calculation of the distance between the atoms i and j*/
//     if(size==(*sys).nb_mol){ /*Mass center*/
//         
//         #pragma omp parallel for schedule(static)
//         for(i=0;i<size;i++){
//             int factor =0;
//             for (int f_i=0; f_i < i; ++f_i)
//                 factor += f_i+1;
//             for(j=i+1;j<size;j++){
//                 const int index = i*size - factor + j -i -1;
//                 loos::AtomicGroup moli = mol[i].copy();
//                 loos::AtomicGroup molj = mol[j].copy();
//                 loos::GCoord distance = mol[i].centerOfMass() - mol[j].centerOfMass();
//                 distance.reimage((*sys).cell);
//                 cout << distance << "  " << index << endl;
//                 x[index] = distance.x();
//                 y[index] = distance.y();
//                 z[index] = distance.z();
//             }
//         }
//     }
//     else{ /*Atomic*/ //NOTE should be killed
//         for(i=0;i<size;i++){
//             for(j=i+1;j<size;j++){
//                 int index = i*size+j;
//                 x[index] = at_coord[3*i+0] - at_coord[3*j+0];
//                 y[index] = at_coord[3*i+1] - at_coord[3*j+1];
//                 z[index] = at_coord[3*i+2] - at_coord[3*j+2];
//                 ij++;
//             }
//         }
//     }
//     cout << "done" << endl;
// 
//     #pragma omp parallel for
//     for(i=0;i<size*(size-1)/2;i++){
//         pbc((*sys).cell,&x[i],&y[i],&z[i]); //NOTE done with reimage above
//         r[i]=sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
//         cout << "Distance: " << r[i] << endl;
//     }
//     
    /* Calculation of the screening functions (Thole)
     * The parameters are already set to 1.d0 if we do not
     * want these screening function
     */
    
    //WARNING remove thole function
    //cannot use atomic code anymore
    //if((*sys).thole){
    //    thole(r,fthole,size);
    //}
    
    /* Calculation of Tij */
    //WARNING now atomic is broken, 
    //I have removed the old tedious code to compute x,y,z,r 
    // without writing code for atomic 
    std::vector<double> x(size*(size-1)/2,0.0);
    std::vector<double> y(size*(size-1)/2,0.0);
    std::vector<double> z(size*(size-1)/2,0.0);
    std::vector<double> r(size*(size-1)/2,0.0);
    //cout << size*(size-1)/2 << endl;
    
    #pragma omp parallel for
    for(int i=0;i<size;i++){
        int factor =0;
        for (int f_i=0; f_i < i; ++f_i)
            factor += f_i+1;
        //cout << i << "  " << factor;
        for(int j=i+1;j<size;j++){
            int ij = i*size - factor +j-i-1;
            //cout << " " << ij << endl;
            loos::GCoord distance = mol[i].centerOfMass() - mol[j].centerOfMass();
            distance.reimage((*sys).cell);
            r[ij] = distance.length();
            x[ij] = distance.x();
            y[ij] = distance.y();
            z[ij] = distance.z();
        }
    }
    //cout <<"current size " << r.size() << endl;
    
    #pragma omp parallel for 
    for(int i=0;i<size;i++){
        int factor =0;
        for (int f_i=0; f_i < i; ++f_i)
            factor += f_i+1;
        for(int j=i+1;j<size;j++){
            int ij = i*size - factor +j-i-1;
            //loos::GCoord distance = mol[i].centerOfMass() - mol[j].centerOfMass();
            //distance.reimage((*sys).cell);
            //double r = distance.length();
            //double x = distance.x();
            //double y = distance.y();
            //double z = distance.z();
            if(r[ij]<=(*sys).cutoff){
                if((*sys).thole){
                    r3=pow(r[ij],-3.)* fthole[2*ij  ]; /*This is 1/r3 with screening function*/
                    r5=3/pow(r[ij],5.)*fthole[2*ij+1]; /*This is 3/r5 with screening function*/
                }
                else{
                    r3=pow(r[ij],-3.); //This is 1/r3
                    r5=3/pow(r[ij],5.); //This is 3/r5
                }
                /*Tij (=Tji) calculation*/
                Tij[ij].xx=r3-r5*x[ij]*x[ij];
                Tij[ij].yx=  -r5*y[ij]*x[ij]; Tij[ij].yy=r3-r5*y[ij]*y[ij];
                Tij[ij].zx=  -r5*z[ij]*x[ij]; Tij[ij].zy=  -r5*z[ij]*y[ij]; Tij[ij].zz=r3-r5*z[ij]*z[ij];
                
            }
            else{
                Tij[ij].xx=0.0;
                Tij[ij].yx=0.0;	Tij[ij].yy=0.0;
                Tij[ij].zx=0.0;	Tij[ij].zy=0.0;	Tij[ij].zz=0.0;
            }
        }
    }

    /*If we do not want intra molecular DID (for 3 points model only)*/
    if(size>(*sys).nb_mol){
        if((*sys).intra==0){
            int ij=0;
            for(int i=0;i<size;i=i+3){/*Only the O are passed*/
                /*No interaction between O and H1*/
                Tij[ij].xx=0.0;
                Tij[ij].yx=0.0;	Tij[ij].yy=0.0;
                Tij[ij].zx=0.0;	Tij[ij].zy=0.0;	Tij[ij].zz=0.0;
                
                /*No interaction between O and H2*/
                ij++;
                Tij[ij].xx=0.0;
                Tij[ij].yx=0.0;	Tij[ij].yy=0.0;
                Tij[ij].zx=0.0;	Tij[ij].zy=0.0;	Tij[ij].zz=0.0;
                
                /*No interaction between H1 and H2*/
                ij=ij+size-i-2;
                Tij[ij].xx=0.0;
                Tij[ij].yx=0.0;	Tij[ij].yy=0.0;
                Tij[ij].zx=0.0;	Tij[ij].zy=0.0;	Tij[ij].zz=0.0;
                
                ij=ij+2*size-2*i-5; /*The intermolecular interactions are passed*/
            }
        }
    }
    
//     #ifdef DEBUG
//     //ij=0; 
//     for(int i=0;i<size;i++){
//         int factor =0;
//         for (int f_i=0; f_i < i; ++f_i)
//             factor += f_i+1;
//         for(int j=i+1;j<size;j++){
//             int ij = i*size - factor +j-i-1;
//             file_Tij << setw(10) << step << setw(10) << i << setw(5) << j
//             << setw(18) << Tij[ij].xx
//             << setw(18) << Tij[ij].yx
//             << setw(18) << Tij[ij].yy
//             << setw(18) << Tij[ij].zx
//             << setw(18) << Tij[ij].zy
//             << setw(18) << Tij[ij].zz << endl;
//         } 
//     } 
//     #endif 
  
    return;
}


/*Protection against infinite loop
 * If we are out of the limits, only the first order dipole
 * moment and polarizability are written for the whole step
 * XXX TRAVAIL ARBEIT WORK
 * Check if it works
 */
void scf_protection(sys_info *sys, const std::vector<loos::AtomicGroup> &mol, std::vector<vect_3d> &dip0, vect_3d *dip, mat_sym_3d *pol0, mat_sym_3d *pol0_mol, mat_3d *pol, mat_sym_3d *Tij_mc, vect_3d *e_ind, vect_3d *v3_tmp, mat_3d *m3_tmp, int step){
    int i=0, j=0, ij=0;
    std::ofstream file_pb;
    //cout << "done5" << endl;
    /*Initialization*/
    init_dip_pol(dip0,dip,pol0,pol,sys);
    //#pragma omp parallel for
    for(int i=0;i<(*sys).nb_mol;i++){
        v3_tmp[i].x()=0.; v3_tmp[i].y()=0.; v3_tmp[i].z()=0.;
    }
    //#pragma omp parallel for
    for(int i=0;i<(*sys).nb_pol;i++){
        m3_tmp[i].xx=0.; m3_tmp[i].xy=0.; m3_tmp[i].xz=0.;
        m3_tmp[i].yx=0.; m3_tmp[i].yy=0.; m3_tmp[i].yz=0.;
        m3_tmp[i].zx=0.; m3_tmp[i].zy=0.; m3_tmp[i].zz=0.;
    }
    
    /*-Sum Tij.dip[j], -Sum Tij.pol[j]*/
    /*-Sum Tji.dip[i], -Sum Tji.pol[i] */
    //cout << "done4" << endl;
    ij=0;
    //#pragma omp parallel for
    for(int i=0;i<(*sys).nb_mol;i++){
        int factor =0;
        for (int f_i=0; f_i < i; ++f_i)
            factor += f_i+1;
        for(int j=i+1;j<(*sys).nb_mol;j++){
            ij = i*(*sys).nb_mol - factor + j -i-1;
            sumT_dip(&(Tij_mc[ij]),dip[j],v3_tmp[i]);
            sumT_dip(&(Tij_mc[ij]),dip[i],v3_tmp[j]);
            //ij++;
        }
    }
    //cout <<"done3" << endl;
    //ij=0;
    for(int i=0;i<(*sys).nb_point;i++){
        int factor =0;
        for (int f_i=0; f_i < i; ++f_i)
            factor += f_i+1;
        for(int j=i+1;j<(*sys).nb_point;j++){
            ij = i*(*sys).nb_point - factor + j -i-1;
            sumT_pol(&(Tij_mc[ij]),(&pol[j]),(&m3_tmp[i]));
            sumT_pol(&(Tij_mc[ij]),(&pol[i]),(&m3_tmp[j]));
            //ij++;
        }
    }
    
    
    /*1-Sum Tij.pol[j]*/
    for(int i=0;i<(*sys).nb_point;i++){
        m3_tmp[i].xx+=1.;
        m3_tmp[i].yy+=1.;
        m3_tmp[i].zz+=1.;
    }
    
    //#pragma omp parallel for
    for(int i=0;i<(*sys).nb_mol;i++){
        mult_msym_v_3d(&(pol0_mol[i]),v3_tmp[i],e_ind[i]);/*INDUCED dipole moment*/
        /*Total dipole*/
        dip[i].x() = dip0[i].x() + e_ind[i].x();
        dip[i].y() = dip0[i].y() + e_ind[i].y();
        dip[i].z() = dip0[i].z() + e_ind[i].z();
        //cout << i << "   " << dip[i] << "   " << e_ind[i] << endl;
    }
    
    //#pragma omp parallel for
    for(int i=0;i<(*sys).nb_point;i++){
        mult_msym_m_3d(&(pol0[i]),&(m3_tmp[i]),&(pol[i]));/*TOTAL polarizability*/
    }
    
    
    /*Record of the positions with a problem of convergence
     *for further debugging*/
    file_pb.open(CONV_PB,std::ofstream::app);
    file_pb << (*sys).nb_mol << endl;
    file_pb << "i = " << step << " , this is the positions which had some convergence problem\n";
    for(int i=0;i<(*sys).nb_mol;i++){
        file_pb << " X  " << setw(15) << mol[i].centerOfMass().x() << setw(15) << mol[i].centerOfMass().y() << setw(15) << mol[i].centerOfMass().z() << endl;
    }
    file_pb.close();
}



/*This subroutine returns the 2 screening functions*/
void thole(double *r, double *fthole, int size){
    int i=0, j=0, ij=0;
    double x=0;
    
    for(int i=0;i<size;i++){
        for(int j=i+1;j<size;j++){
            x=A_THOLE*r[ij];
            
            /* Screening of the term in 1/r3 */
            fthole[2*ij  ]=1-(1 + x + x*x/2.)*exp(-x);
            /* Screening of the term in 1/r5 */
            fthole[2*ij+1]=1-(1 + x + x*x/2. + pow(x,3.)/6)*exp(-x);
            
            ij++;
        }
    }
    
};


