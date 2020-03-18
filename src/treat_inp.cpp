/*
  In this file can be found the subroutines
  which are necessary with the input files
  (Either from the user or from CP2K)

  I suppose that the xyz-file comes from CP2K
  and therefore is formated!!!!!!
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include "main.h"
#include "dip_pol.h"
#include "my_math.h"
#include "treat_inp.h"
#include "tools.hpp"
#include <loos.hpp>
#include <vector>
#include <sstream>


using std::cout;
using std::endl;
using std::setw;

/*==================================*/
/*Reading the input file of the user*/
/*==================================*/
int read_input(input_info *input, sys_info *sys, int argc, char *argv[]){
  int i=0;
  //int first_H=0;
  char symb=0;
  char chain[CHAIN_SIZE]="",p[CHAIN_SIZE]="",q[CHAIN_SIZE]="",r[CHAIN_SIZE]="";
  std::ifstream file;
  
  if(argc!=2){
    printf("You have to specify the name of the input file:\n");
    printf("%s input_file\n",argv[0]);
    printf("\nEnd of program\n");
    exit(0);
  }
  //FOPEN_SAFE(file,argv[1],"r");
  file.open(argv[1]);

  /* Position file */
  //fgets(chain,CHAIN_SIZE,file);
  file.getline(chain,CHAIN_SIZE);
  sscanf(chain,"%s",(*input).namei);

  /* topology file */ //Hossam
  //fgets(chain,CHAIN_SIZE,file); //Hossam
  file.getline(chain,CHAIN_SIZE);
  sscanf(chain,"%s",(*input).namet); //Hossam

  loos::AtomicGroup current_frame = loos::createSystem((*input).namet); //Hossam
  //cout << current_frame.size() << endl; //Hossam
  //cout << current_frame.splitByMolecule().size() << endl; //Hossam
  
  /* Output file (dippol.dat) */
  //fgets(chain,CHAIN_SIZE,file);
  file.getline(chain,CHAIN_SIZE);
  sscanf(chain,"%s",(*input).nameo);

#ifndef YUKI
  /*Output file for the mass center position*/
  //fgets(chain,CHAIN_SIZE,file);
  file.getline(chain,CHAIN_SIZE);
  sscanf(chain,"%s",(*input).name_mass);
#endif
  
  /*Size of the box*/
  //fgets(chain,CHAIN_SIZE,file);
  //sscanf(chain,"%lf %lf %lf",&((*sys).cell.x()),&((*sys).cell.y()),&((*sys).cell.z()) );
  file >> (*sys).cell.x();
  file >> (*sys).cell.y();
  file >> (*sys).cell.z();
  cout << "Box size is: " << (*sys).cell << endl;
  file.ignore(1000,'\n');

  /*Initial / Final step*/
  //fgets(chain,CHAIN_SIZE,file);
  file.getline(chain,CHAIN_SIZE);
  sscanf(chain,"%d %d",&((*input).stepi),&((*input).stepf));

  /*Order of the atoms (only the first character will be read)*/
  //fgets(chain,CHAIN_SIZE,file);
  file.getline(chain,CHAIN_SIZE);
  sscanf(chain,"%d",&((*sys).at_p_mol));/*Atoms per molecule*/
  //MALLOC_SAFE((*sys).order,(*sys).at_p_mol,int);

  //first_H=1;
  for(i=0;i<(*sys).at_p_mol;i++){
    //symb=fgetc(file);
    file.get(symb);
    //file.ignore(1000,'\n');
    //cout << symb << endl;
    //fgets(chain,CHAIN_SIZE,file);
    file.getline(chain,CHAIN_SIZE);

//     if(symb=='O'){(*sys).order[i]=0;}
//     else if(symb=='H'){
//       if(first_H){
// 	(*sys).order[i]=3;/*H1*/
// 	first_H=0;
//       }
//       else{
// 	(*sys).order[i]=6;/*H2*/
//       }
//     }
//     else{(*sys).order[i]=9;}/*Any dummy atom*/
  }

  /*
    Polarization of the beam
    For (*sys).p (*sys).q (*sys).r:
    X=0, Y=1, Z=2
  */
  //fgets(chain,CHAIN_SIZE,file);
  file.getline(chain,CHAIN_SIZE);
  sscanf(chain,"%s %s %s",p,q,r);
  //cout << p << q << r << endl;
  if(strcmp(r,"X")==0){(*sys).r=0;}
  else if(strcmp(r,"Y")==0){(*sys).r=1;}
  else if(strcmp(r,"Z")==0){(*sys).r=2;}
  else{printf("Problem with the polarization %s\nEnd of program.\n",r);exit(0);}

  if(strcmp(q,"X")==0){(*sys).q=0;}
  else if(strcmp(q,"Y")==0){(*sys).q=1;}
  else if(strcmp(q,"Z")==0){(*sys).q=2;}
  else if(strcmp(q,"S")==0){(*sys).q=-(*sys).r-1;}
  else{printf("Problem with the polarization %s\nEnd of program.\n",q);exit(0);}

  if(strcmp(p,"X")==0){(*sys).p=0;}
  else if(strcmp(p,"Y")==0){(*sys).p=1;}
  else if(strcmp(p,"Z")==0){(*sys).p=2;}
  else if(strcmp(p,"S")==0){
    (*sys).p=-(*sys).r-1;
    /* Extra control for SS {X,Y,Z} polarization */
    if(strcmp(p,q)!=0){
      printf("Problem with the polarization\nYou certainly wanted to write \"S S %s\"\nEnd of program.\n",r);exit(0);
    }
  }
  else{printf("Problem with the polarization %s\nEnd of program.\n",p);exit(0);}
  
  /*Cutoff for the induction*/
  //fgets(chain,CHAIN_SIZE,file);
  file.getline(chain,CHAIN_SIZE);
  sscanf(chain,"%lf",&((*sys).cutoff));

  /*Do we use the Thole damping factor for Tij?*/
  //fgets(chain,CHAIN_SIZE,file);
  file.getline(chain,CHAIN_SIZE);
  sscanf(chain,"%d",&((*sys).thole));

  /*Atomic or molecular dipole moment?*/
  //symb=fgetc(file);
  file.get(symb);
  //fgets(chain,CHAIN_SIZE,file);
  file.getline(chain,CHAIN_SIZE);
  if(symb=='a'||symb=='A'){
    (*sys).typ_dip=1;/* 1=Atomic dipole moment */
    //fgets(chain,CHAIN_SIZE,file);
    file.getline(chain,CHAIN_SIZE);
    sscanf(chain,"%lf %lf",&((*sys).MO_z),&((*sys).MH_z));
  }
  else if(symb=='m'||symb=='M'){
    (*sys).typ_dip=0;/* 0=Molecular dipole moment */
    //fgets(chain,CHAIN_SIZE,file);
    file.getline(chain,CHAIN_SIZE);
    sscanf(chain,"%lf %lf %lf",&((*sys).M_x),&((*sys).M_y),&((*sys).M_z));
    //cout << std::setw(12) << (*sys).M_x << (*sys).M_y << (*sys).M_z << endl;
  }
  else{
    printf("Problem with the kind of dipole moment. Please chose \"A(a)tomic\" or \"M(m)olecular\"\nEnd of program.\n");exit(0);
  }

  /*Atomic or molecular polarizability?*/
  //symb=fgetc(file);
  file.get(symb);
  //fgets(chain,CHAIN_SIZE,file);
  file.getline(chain,CHAIN_SIZE);
  if(symb=='a'||symb=='A'){
    (*sys).typ_pol=1;/* 1=Atomic polarizability */
    //fgets(chain,CHAIN_SIZE,file);
    file.getline(chain,CHAIN_SIZE);
    sscanf(chain,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",			\
	   &((*sys).AO_xx),&((*sys).AO_yy),&((*sys).AO_zz),		\
	   &((*sys).AH_xx),&((*sys).AH_yy),&((*sys).AH_zz),		\
	   &((*sys).AH_yx),&((*sys).AH_zx),&((*sys).AH_zy));
    //cout << (*sys).AO_xx << " " << &((*sys).AO_yy) << " " << &((*sys).AO_zz) << endl;
  }
  else if(symb=='m'||symb=='M'){
      //cout << "Using molecular polarizability" << endl;
    if((*sys).typ_dip){
      printf("What you are doing is not possible\n");
      printf("To have the ATOMIC induced dipole moment,\n");
      printf("you need an ATOMIC polarizability\n");
      printf("alpha0i * Tij * mu_j  (with j!=i)\n");
      printf("End of program\n");
      exit(0);
    }
    (*sys).typ_pol=0;/* 0=Molecular polarizability */
    //fgets(chain,CHAIN_SIZE,file);
    //file.ignore(1000,'\n');
    file.getline(chain,CHAIN_SIZE);
    sscanf(chain,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",		\
    &((*sys).A_xx),&((*sys).A_xy),&((*sys).A_xz),	\
    &((*sys).A_yx),&((*sys).A_yy),&((*sys).A_zz),	\
    &((*sys).A_zx),&((*sys).A_zy),&((*sys).A_zz));
    //cout << (*sys).A_xx << " " << (*sys).A_xy << " " << (*sys).A_xz << endl;
    //cout << (*sys).A_yx << " " << (*sys).A_yy << " " << (*sys).A_yz << endl;
    //cout << (*sys).A_zx << " " << (*sys).A_zy << " " << (*sys).A_zz << endl;
  }
  else{
    printf("Problem with the kind of polarizability. Please chose \"A(a)tomic\" or \"M(m)olecular\"\nEnd of program.\n");exit(0);
  }
  //file.ignore(1000,'\n');
  /*Do we take into account the intra molecular DID?*/
  //fgets(chain,CHAIN_SIZE,file);
  file.getline(chain,CHAIN_SIZE);
  sscanf(chain,"%d",&((*sys).intra));
  //cout << (*sys).intra << endl;
  
  //NOTE now new code to parse ref molecules.
  //NOTE should at some point create a map from the first current
  //NOTE linking each molecule id to a pointer to the correct reference mol
    std::string line;
    std::stringstream ss;
    while (getline(file,line))
    {
        if (line[0] == '[')
        {   
            bool found_dip=0, found_alpha=0;
            mol_type new_type;
            std::string name = line;
            std::size_t name_end = name.find("]");
            if (name_end ==std::string::npos)
            {
               cout << "Ill-formated residue name, exiting!" << endl;
               exit(0);
            }
            name = name.substr(1,name_end-1);
            int mol_size;
            file >> mol_size;
            file.ignore(1000,'\n');
            loos::AtomicGroup new_mol(3);
            for (int i=0; i < new_mol.size(); i++)
            {
                std::string at_label;
                double at_x, at_y, at_z;
                file >> at_label >> at_x >> at_y >> at_z;
                file.ignore(1000,'\n');
                new_mol[i]->index(i+1);
                new_mol[i]->id(i);
                new_mol[i]->name(at_label);
                new_mol[i]->coords(loos::GCoord(at_x,at_y,at_z));
            }
            cout << "found molecule with name " << name << endl;
            cout << "deducing masses from atomic labels:" << endl;
            for (int i=0; i<new_mol.size(); i++)
                new_mol[i]->mass(deduceMass(new_mol[i]->name()));
            new_mol.deduceAtomicNumberFromMass(0.1);
            cout << new_mol << endl;
            
            new_type.ref = new_mol.copy();
            new_type.ref.centerAtOrigin();

            //calculate the PAS of the ref once and for all!
            cout << "Computing principal axes of reference molecule " << name << endl;
            
            Eigen::Matrix3d tensor = inertia_tensor(new_type.ref);
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> pas_es(tensor);
            new_type.princ = pas_es.eigenvectors();
            cout << "Principal axes:" << endl;
            cout << new_type.princ << endl << endl;
            cout << "Principal moments: ";
            cout << pas_es.eigenvalues().transpose() << endl;
            
            //cout << "PAS from LOOS" << endl;
            //std::vector<loos::GCoord> pas = new_mol.momentsOfInertia();
            //for (int i=0; i < pas.size(); i++)
            //    cout << pas[i] << endl;
            
            getline(file,line);
            ss.str(line);
            std::string first;
            ss >> first;
            while ( (first != "END") && file.good() )
            {
                if (first == "DIP")
                {
                    found_dip = 1;
                    ss >> new_type.dip[0] >> new_type.dip[1] >> new_type.dip[2];
                    cout << "Dipole: " << setw(10) << new_type.dip[0]
                    << setw(10) << new_type.dip[1]
                    << setw(10) << new_type.dip[2] << endl;
                }
                else if ( first == "ALPHA")
                {
                    found_alpha = 1;
                    ss >> new_type.alpha[0][0] >> new_type.alpha[0][1] >> new_type.alpha[0][2]
                    >> new_type.alpha[1][0] >> new_type.alpha[1][1] >> new_type.alpha[1][2]
                    >> new_type.alpha[2][0] >> new_type.alpha[2][1] >> new_type.alpha[2][2];
                    new_type.alpha_m << new_type.alpha[0][0] , new_type.alpha[0][1] , new_type.alpha[0][2]
                    , new_type.alpha[1][0] , new_type.alpha[1][1] , new_type.alpha[1][2]
                    , new_type.alpha[2][0] , new_type.alpha[2][1] , new_type.alpha[2][2];
                    cout << "Alpha\n";
                }
                else if (first == "BETA" )
                {
                    cout << "found beta\n";
                }
                getline(file,line);
                ss.str(line);
                ss >> first;
            }
            if (!found_dip || !found_alpha)
            {
                cout << "Dipole or alpha tensor not provided, exiting!" << endl;
                exit(0);
            }
            (*sys).all_mol_types[name] = new_type;
            cout << endl;
        }
    }
    
    //cout << (*sys).all_mol_types["WAT"].dip[2] << endl;
    //cout << (*sys).all_mol_types["CHL"].dip[2] << endl;
    file.close();
    return 0;
}



/*==============================*/
/*Reading of the trajectory file*/
/*==============================*/
int trajectory(sys_info *sys, input_info input, char *argv[]) {
    //H FILE *filei=NULL;
    std::ofstream file_debug;
    int step=0, i=0, j=0;
    //double time=0.;
    //char chain[CHAIN_SIZE]="";
    
    //double *x=NULL, *y=NULL, *z=NULL, *r=NULL, *fthole=NULL, *at_coord;
    double *fthole=NULL;
    //double *at_coord;
    vect_3d *dip=NULL, *dip_mol=NULL, *e_ind=NULL, *v3_tmp=NULL;
    //mol_info *mol=NULL;
    mat_sym_3d *pol0=NULL, *pol0_mol=NULL, *Tij_mc=NULL, *Tij_at=NULL;
    mat_3d *pol=NULL, *pol_mol=NULL, *a_ind=NULL, *m3_tmp=NULL;
    std::ofstream fileo;
    std::ofstream file_traj, file_dip0, file_pol0, file_Tij_mc, file_Tij_at;
    std::ofstream file_dip,file_pol,file_dip_ind, file_pol_ind;
    
    /*To facilitate the reading, I extract the info from structure input*/
    //H int stepi=input.stepi, stepf=input.stepf;
    //H char namei[CHAIN_SIZE]="";
    //H strcpy(namei,input.namei);
    
    //Open the traj and topology as loos objects
    loos::AtomicGroup topology = loos::createSystem(input.namet); //Hossam
    topology.periodicBox(sys->cell);
    //cout << topology.size() << endl; //Hossam
    //cout << topology.splitByMolecule().size() << endl; //Hossam
    loos::pTraj traj = createTrajectory(input.namei,topology);//Hossam
    traj->readFrame(0);
    traj->updateGroupCoords(topology);
    std::vector<loos::AtomicGroup> all_mols = topology.splitByMolecule();
    
    /*Preparation of the different values requiered for the allocation*/
    //H (*sys).nb_mol  =(*sys).nb_at/(*sys).at_p_mol;
    //H (*sys).nb_point=(*sys).nb_at/(*sys).at_p_mol;
    (*sys).nb_mol = all_mols.size();
    (*sys).nb_point = (*sys).nb_mol;
    (*sys).nb_dip   = (*sys).nb_mol;
    (*sys).nb_pol   = (*sys).nb_mol;
    
    if((*sys).typ_dip){/*Atomic dipole*/
        //H (*sys).nb_dip  =3*(*sys).nb_mol;
        //H (*sys).nb_point=3*(*sys).nb_mol;
        (*sys).nb_dip   = topology.size();
        (*sys).nb_point = topology.size();
    }
    if((*sys).typ_pol){/*Atomic polarizability*/
        //H (*sys).nb_pol  =3*(*sys).nb_mol;
        //H (*sys).nb_point=3*(*sys).nb_mol;
        (*sys).nb_pol   = topology.size();
        (*sys).nb_point = topology.size();
    }
    
    fileo.open(input.nameo);
    #ifndef YUKI
    file_traj.open(input.name_mass);
    fileo << "This file has been written by " << argv[0] << ". Edit it at your own risk.\n";
    fileo << "#Mass_center_file " << input.name_mass << endl;
    fileo << "#Cell_parameters " << (*sys).cell.x() << " " << (*sys).cell.y() << " " << (*sys).cell.z() << endl;
    fileo << "#Cutoff_Ang " << (*sys).cutoff << endl;
    fileo << "#Dip0 Dip_comp Pol0 Pol_comp\n";
    #endif
    
    
    //NOTE need to carefully check the following allocations
    /*Creation of the different arrays and files*/
    if((*sys).typ_dip + (*sys).typ_pol){ /* Atomic dipole moment or polarizability */
        MALLOC_SAFE(Tij_at,(*sys).nb_point*((*sys).nb_point-1)/2,mat_sym_3d);
    }
    //H MALLOC_SAFE(mol     ,(*sys).nb_mol                        ,mol_info);
    //MALLOC_SAFE(at_coord,3*topology.size()                      ,double);
    //MALLOC_SAFE(dip0    ,(*sys).nb_dip                        ,vect_3d);
    MALLOC_SAFE(dip     ,(*sys).nb_dip                        ,vect_3d);
    //MALLOC_SAFE(dip0_mol,(*sys).nb_mol                        ,vect_3d);
    MALLOC_SAFE(dip_mol ,(*sys).nb_mol                        ,vect_3d);
    MALLOC_SAFE(pol0    ,(*sys).nb_pol                        ,mat_sym_3d);
    MALLOC_SAFE(pol     ,(*sys).nb_pol                        ,mat_3d);
    MALLOC_SAFE(pol0_mol,(*sys).nb_mol                        ,mat_sym_3d);
    MALLOC_SAFE(pol_mol ,(*sys).nb_mol                        ,mat_3d);
    MALLOC_SAFE(e_ind   ,(*sys).nb_dip                        ,vect_3d);
    MALLOC_SAFE(v3_tmp  ,(*sys).nb_dip                        ,vect_3d);
    MALLOC_SAFE(a_ind   ,(*sys).nb_pol                        ,mat_3d);
    MALLOC_SAFE(m3_tmp  ,(*sys).nb_pol                        ,mat_3d);
    MALLOC_SAFE(Tij_mc  ,(*sys).nb_mol*((*sys).nb_mol-1)/2    ,mat_sym_3d);
    MALLOC_SAFE(fthole  ,(*sys).nb_point*((*sys).nb_point-1)  ,double);
    
    /* Initialization of thole screening factor to 1 (No thole factor)*/
    for(i=0;i<(*sys).nb_point*((*sys).nb_point-1);i++){
        fthole[i]=1.;
    }
    
    
    #ifdef DEBUG
    file_debug.open(DEBUG_TRAJ);
    file_dip0.open(DEBUG_DIP0);
    file_pol0.open(DEBUG_POL0);
    file_Tij_mc.open(DEBUG_TIJ_MC);
    file_dip_ind.open(DEBUG_DIP_IND);
    file_dip.open(DEBUG_DIP);
    file_pol.open(DEBUG_POL);
    file_pol_ind.open(DEBUG_POL_IND);
    file_debug << "#Step_number Record_number Mol_number atom_number x y z\n";
    file_dip0 << "#Record_number Mol_number dip0_x dip0_y dip0_z\n";
    file_pol0 << "#Record_number Mol_number pol0_xx pol0_yx pol0_yy pol0_zx pol0_zy pol0_zz\n";
    file_Tij_mc << "#Record_number Mol_i Mol_j Tij_xx Tij_yx Tij_yy Tij_zx Tij_zy Tij_zz\n";
    file_dip_ind << "#Record_number Mol_number e_ind_x e_ind_y e_ind_z at the order N\n";
    file_dip << "#Record_number Mol_number dip_x dip_y dip_z\n";
    file_pol << "#Record_number Mol_number pol_xx pol_xy pol_xz pol_yx pol_yy pol_yz pol_zx pol_zy pol_zz\n";
    file_pol_ind << "#Record_number Mol_number a_ind_xx xy xz yx yy yz zx zy zz at the order N\n";
    if((*sys).typ_dip + (*sys).typ_pol){ /* Atomic dipole moment or polarizability */
        file_Tij_at.open(DEBUG_TIJ_AT);
        file_Tij_at << "#Record_number Mol_i Mol_j Tij_xx Tij_yx Tij_yy Tij_zx Tij_zy Tij_zz\n";
    }
    #endif
    
    #ifndef YUKI
    file_traj << (*sys).nb_mol << endl << endl;
    #endif
    
    //NOW loop over trajectory frames
    for(int frame_i=input.stepi; (frame_i < input.stepf) && (frame_i < traj->nframes()) ; frame_i++)
    {
        traj->readFrame(frame_i);
        traj->updateGroupCoords(topology);
        //NOTE probably the culprit, Kill it!
        //at_coord = topology.coordsAsVector().data();
        //cout << "SIZESIZE " << topology.coordsAsVector().size();
        std::vector<vect_3d> dip0((*sys).nb_dip);
        std::vector<vect_3d> dip0_mol((*sys).nb_mol);
                
        
        //#pragma omp parallel for
        for (int mol_i = 0; mol_i < all_mols.size(); mol_i++)
        {
            all_mols[mol_i].reimage();
        }
        
        /*================================================================================
         * Calculation of* the complete dipole moment / polarizability (permanent + induced)
         *================================================================================*/
        cout << "calling val_init " << endl;
        val_init(sys,all_mols,dip0,dip0_mol,pol0,pol0_mol,file_traj,file_dip0,file_pol0,step);/*dip0, pol0*/
        cout << "done val_init" << endl;
        comp_dip_pol(sys, all_mols, dip0, dip0_mol, dip, dip_mol, pol0, pol0_mol, pol, pol_mol, e_ind, a_ind, v3_tmp, m3_tmp, Tij_mc, Tij_at, fthole, fileo, file_Tij_mc, file_Tij_at, file_dip, file_pol, file_dip_ind, file_pol_ind, step);
        
        
        (*sys).step_max++; /* Real number of steps recorded (can be lower than stepf-stepi+1) */
    }//end loop over trajectory
    
    /* Release the memory */
    free(dip);
    free(dip_mol);
    free(pol0);
    free(pol);
    free(pol0_mol);
    free(pol_mol);
    free(e_ind);
    free(v3_tmp);
    free(m3_tmp);
    free(Tij_mc);
    free(fthole);
    free(a_ind);
    if((*sys).typ_dip + (*sys).typ_pol){free(Tij_at);}
    
    fileo.close();
    #ifndef YUKI
    file_traj.close();
    #endif
    #ifdef DEBUG
    file_dip0.close();
    file_pol0.close();
    file_Tij_mc.close();
    file_debug.close();
    file_dip.close();
    file_pol.close();
    file_dip_ind.close();
    if((*sys).typ_dip + (*sys).typ_pol){file_Tij_at.close();}
    #endif
    
    return 0;
}

