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
#include "../include/main.h"
#include "../include/dip_pol.h"
#include "../include/my_math.h"
#include "../include/treat_inp.h"
#include <loos.hpp>
#include <vector>
#include <sstream>

/*==================================*/
/*Reading the input file of the user*/
/*==================================*/
void read_input(input_info *input, sys_info *sys, int argc, char *argv[]){
  int i=0, first_H=0;
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
  std::cout << current_frame.size() << std::endl; //Hossam
  std::cout << current_frame.splitByMolecule().size() << std::endl; //Hossam
  
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
  std::cout << "Box size is: " << (*sys).cell << std::endl;
  file.ignore(1000,'\n');

  /*Initial / Final step*/
  //fgets(chain,CHAIN_SIZE,file);
  file.getline(chain,CHAIN_SIZE);
  sscanf(chain,"%d %d",&((*input).stepi),&((*input).stepf));

  /*Order of the atoms (only the first character will be read)*/
  //fgets(chain,CHAIN_SIZE,file);
  file.getline(chain,CHAIN_SIZE);
  sscanf(chain,"%d",&((*sys).at_p_mol));/*Atoms per molecule*/
  MALLOC_SAFE((*sys).order,(*sys).at_p_mol,int);

  first_H=1;
  for(i=0;i<(*sys).at_p_mol;i++){
    //symb=fgetc(file);
    file.get(symb);
    //file.ignore(1000,'\n');
    //std::cout << symb << std::endl;
    //fgets(chain,CHAIN_SIZE,file);
    file.getline(chain,CHAIN_SIZE);

    if(symb=='O'){(*sys).order[i]=0;}
    else if(symb=='H'){
      if(first_H){
	(*sys).order[i]=3;/*H1*/
	first_H=0;
      }
      else{
	(*sys).order[i]=6;/*H2*/
      }
    }
    else{(*sys).order[i]=9;}/*Any dummy atom*/
  }

  /*
    Polarization of the beam
    For (*sys).p (*sys).q (*sys).r:
    X=0, Y=1, Z=2
  */
  //fgets(chain,CHAIN_SIZE,file);
  file.getline(chain,CHAIN_SIZE);
  sscanf(chain,"%s %s %s",p,q,r);
  //std::cout << p << q << r << std::endl;
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
    //std::cout << std::setw(12) << (*sys).M_x << (*sys).M_y << (*sys).M_z << std::endl;
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
    //std::cout << (*sys).AO_xx << " " << &((*sys).AO_yy) << " " << &((*sys).AO_zz) << std::endl;
  }
  else if(symb=='m'||symb=='M'){
      //std::cout << "Using molecular polarizability" << std::endl;
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
    //std::cout << (*sys).A_xx << " " << (*sys).A_xy << " " << (*sys).A_xz << std::endl;
    //std::cout << (*sys).A_yx << " " << (*sys).A_yy << " " << (*sys).A_yz << std::endl;
    //std::cout << (*sys).A_zx << " " << (*sys).A_zy << " " << (*sys).A_zz << std::endl;
  }
  else{
    printf("Problem with the kind of polarizability. Please chose \"A(a)tomic\" or \"M(m)olecular\"\nEnd of program.\n");exit(0);
  }
  //file.ignore(1000,'\n');
  /*Do we take into account the intra molecular DID?*/
  //fgets(chain,CHAIN_SIZE,file);
  file.getline(chain,CHAIN_SIZE);
  sscanf(chain,"%d",&((*sys).intra));
  //std::cout << (*sys).intra << std::endl;
  
  //NOTE now new code to parse ref molecules.
  //NOTE should at some point create a map from the first current
  //NOTE linking each molecule id to a pointer to the correct reference mol
    std::string line;
    std::stringstream ss;
    while (getline(file,line))
    {
        if (line[0] == '[')
        {   
            mol_type new_type;
            std::string name = line;
            std::size_t name_end = name.find("]");
            if (name_end ==std::string::npos)
            {
               std::cout << "Ill-formated residue name, exiting!" << std::endl;
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
                new_mol[i]->name(at_label);
                new_mol[i]->coords(loos::GCoord(at_x,at_y,at_z));
            }
            new_type.ref = new_mol.copy();
            std::cout << "found molecule with name " << name << std::endl;
            std::cout << new_mol << std::endl;
            getline(file,line);
            ss.str(line);
            std::string first;
            ss >> first;
            while ( (first != "END") && file.good() )
            {
                if (first == "DIP")
                {
                  ss >> new_type.dip[0] >> new_type.dip[1] >> new_type.dip[2];
                  std::cout << "found dip\n";
                  std::cout << new_type.dip[0] << new_type.dip[1] << new_type.dip[2] << std::endl;
                }
                else if ( first == "ALPHA")
                {
                    ss >> new_type.alpha[0][0] >> new_type.alpha[0][1] >> new_type.alpha[0][2]
                         >> new_type.alpha[1][0] >> new_type.alpha[1][1] >> new_type.alpha[1][2]
                         >> new_type.alpha[2][0] >> new_type.alpha[2][1] >> new_type.alpha[2][2];
                    std::cout << "found alpha\n";
                }
                else if (first == "BETA" )
                {
                   std::cout << "found beta\n";
                }
                getline(file,line);
                ss.str(line);
                ss >> first;
            }
            (*sys).all_mol_types[name] = new_type;
            std::cout << std::endl;
        }
    }

    //std::cout << (*sys).all_mol_types["WAT"].dip[2] << std::endl;
    //std::cout << (*sys).all_mol_types["CHL"].dip[2] << std::endl;
    file.close();
}



/*==============================*/
/*Reading of the trajectory file*/
/*==============================*/
void trajectory(sys_info *sys, input_info input, char *argv[]) {
    //H FILE *filei=NULL;
    FILE *file_debug=NULL;
    int step=0, i=0, j=0, start_read=0, end_read=0;
    double time=0.;
    char chain[CHAIN_SIZE]="";
    
    double *x=NULL, *y=NULL, *z=NULL, *r=NULL, *fthole=NULL, *at_coord;
    vect_3d *dip0=NULL,*dip0_mol=NULL, *dip=NULL, *dip_mol=NULL, *e_ind=NULL, *v3_tmp=NULL;
    //mol_info *mol=NULL;
    mat_sym_3d *pol0=NULL, *pol0_mol=NULL, *Tij_mc=NULL, *Tij_at=NULL;
    mat_3d *pol=NULL, *pol_mol=NULL, *a_ind=NULL, *m3_tmp=NULL;
    FILE *fileo=NULL, *file_traj=NULL,*file_dip0=NULL,*file_pol0=NULL;
    FILE *file_Tij_mc=NULL, *file_Tij_at=NULL, *file_dip=NULL,*file_pol=NULL,*file_dip_ind=NULL, *file_pol_ind=NULL;
    
    /*To facilitate the reading, I extract the info from structure input*/
    //H int stepi=input.stepi, stepf=input.stepf;
    //H char namei[CHAIN_SIZE]="";
    //H strcpy(namei,input.namei);
    
    //Open the traj and topology as loos objects
    loos::AtomicGroup topology = loos::createSystem(input.namet); //Hossam
    std::cout << topology.size() << std::endl; //Hossam
    //std::cout << topology.splitByMolecule().size() << std::endl; //Hossam
    loos::pTraj traj = createTrajectory(input.namei,topology);//Hossam
    
    //H FOPEN_SAFE(filei,namei,"r");
    FOPEN_SAFE(fileo,input.nameo,"w+");
    #ifndef YUKI
    FOPEN_SAFE(file_traj,input.name_mass,"w+");
    fprintf(fileo,"#This file has been written by %s. Edit it at your own risk.\n",argv[0]);
    fprintf(fileo,"#Mass_center_file %s\n",input.name_mass);
    fprintf(fileo,"#Cell_parameters %f %f %f\n",(*sys).cell.x(),(*sys).cell.y(),(*sys).cell.z());
    fprintf(fileo,"#Cutoff_Ang %f\n",(*sys).cutoff);
    fprintf(fileo,"#Dip0 Dip_comp Pol0 Pol_comp\n");
    #endif
    
    /*The trajectory file is read while we have not reached EOF
     or* the maximal value authorized by the user */
    for(int frame_i=input.stepi; (frame_i < input.stepf) && (frame_i < traj->nframes()) ; frame_i++)
    {
        traj->readFrame(frame_i);
        traj->updateGroupCoords(topology);
    //H while(fscanf(filei,"%d",&((*sys).nb_at)) != EOF && !end_read){
    //H    /*Comment line reading*/
    //H    fscanf(filei,"%s",chain);/* Pass the i */
    //H    fscanf(filei,"%s",chain);/* Pass the = */
    //H    fscanf(filei,"%d",&step);
    //H    fscanf(filei,"%s",chain);/*Pass the ,*/
    //H    fscanf(filei,"%s",chain);/*Pass the time*/
    //H    fscanf(filei,"%s",chain);/*Pass the =*/
    //H    fscanf(filei,"%lf",&time);
        if (frame_i == input.stepi)
        {
            
            /*========================*/
            /*First step selected only*/
            /*========================*/
            //H if(step>=input.stepi ){ //NOTE why >=??
            //H start_read=1;
            
            /*Preparation of the different values requiered for the allocation*/
            //H (*sys).nb_mol  =(*sys).nb_at/(*sys).at_p_mol;
            //H (*sys).nb_point=(*sys).nb_at/(*sys).at_p_mol;
            (*sys).nb_mol   = topology.splitByMolecule().size();
            (*sys).nb_point = (*sys).nb_mol;
            (*sys).nb_dip   = (*sys).nb_mol;
            (*sys).nb_pol   = (*sys).nb_mol;
            
            //read the molecules here and populate nb_dip and nb_point
            
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
            
            
            //NOTE need to carefully check the following allocations
            /*Creation of the different arrays and files*/
            if((*sys).typ_dip + (*sys).typ_pol){ /* Atomic dipole moment or polarizability */
                MALLOC_SAFE(Tij_at,(*sys).nb_point*((*sys).nb_point-1)/2,mat_sym_3d);
            }
            //H MALLOC_SAFE(mol     ,(*sys).nb_mol                        ,mol_info);
            MALLOC_SAFE(at_coord,9*(*sys).nb_mol                      ,double);
            MALLOC_SAFE(dip0    ,(*sys).nb_dip                        ,vect_3d);
            MALLOC_SAFE(dip     ,(*sys).nb_dip                        ,vect_3d);
            MALLOC_SAFE(dip0_mol,(*sys).nb_mol                        ,vect_3d);
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
            MALLOC_SAFE(x       ,(*sys).nb_point*((*sys).nb_point-1)/2,double);
            MALLOC_SAFE(y       ,(*sys).nb_point*((*sys).nb_point-1)/2,double);
            MALLOC_SAFE(z       ,(*sys).nb_point*((*sys).nb_point-1)/2,double);
            MALLOC_SAFE(r       ,(*sys).nb_point*((*sys).nb_point-1)/2,double);
            MALLOC_SAFE(fthole  ,(*sys).nb_point*((*sys).nb_point-1)  ,double);
            
            /* Initialization of thole screening factor to 1 (No thole factor)*/
            for(i=0;i<(*sys).nb_point*((*sys).nb_point-1);i++){
                fthole[i]=1.;
            }
            
            
            #ifdef DEBUG
            FOPEN_SAFE(file_debug,DEBUG_TRAJ,"w+");
            FOPEN_SAFE(file_dip0,DEBUG_DIP0,"w+");
            FOPEN_SAFE(file_pol0,DEBUG_POL0,"w+");
            FOPEN_SAFE(file_Tij_mc ,DEBUG_TIJ_MC ,"w+");
            FOPEN_SAFE(file_dip_ind,DEBUG_DIP_IND,"w+");
            FOPEN_SAFE(file_dip    ,DEBUG_DIP    ,"w+");
            FOPEN_SAFE(file_pol    ,DEBUG_POL    ,"w+");
            FOPEN_SAFE(file_pol_ind,DEBUG_POL_IND,"w+");
            fprintf(file_debug,"#Step_number Record_number Mol_number atom_number x y z\n");
            fprintf(file_dip0,"#Record_number Mol_number dip0_x dip0_y dip0_z\n");
            fprintf(file_pol0,"#Record_number Mol_number pol0_xx pol0_yx pol0_yy pol0_zx pol0_zy pol0_zz\n");
            fprintf(file_Tij_mc,"#Record_number Mol_i Mol_j Tij_xx Tij_yx Tij_yy Tij_zx Tij_zy Tij_zz\n");
            fprintf(file_dip_ind,"#Record_number Mol_number e_ind_x e_ind_y e_ind_z at the order N\n");
            fprintf(file_dip,"#Record_number Mol_number dip_x dip_y dip_z\n");
            fprintf(file_pol,"#Record_number Mol_number pol_xx pol_xy pol_xz pol_yx pol_yy pol_yz pol_zx pol_zy pol_zz\n");
            fprintf(file_pol_ind,"#Record_number Mol_number a_ind_xx xy xz yx yy yz zx zy zz at the order N\n");
            if((*sys).typ_dip + (*sys).typ_pol){ /* Atomic dipole moment or polarizability */
                FOPEN_SAFE(file_Tij_at ,DEBUG_TIJ_AT ,"w+");
                fprintf(file_Tij_at,"#Record_number Mol_i Mol_j Tij_xx Tij_yx Tij_yy Tij_zx Tij_zy Tij_zz\n");
            }
            #endif
        //H }     
        }
        
        
        //H if(frame_i>=input.stepf){end_read=1;}
        //H fgets(chain,CHAIN_SIZE,filei);/* Read the end of the line */
        
        std::vector<loos::AtomicGroup> all_mols = topology.splitByMolecule();
        for (int mol_i = 0; mol_i < all_mols.size(); mol_i++)
        {
            all_mols[mol_i].reimage();
            std::cout << "reimaged molecule " << mol_i << " with " << (*sys).cell << std::endl;
        }
 
//H         /*Body read if we are after the starting step*/
//H         if(start_read){
//H             #ifndef YUKI
//H             fprintf(file_traj,"%d\n",(*sys).nb_mol);
//H             fprintf(file_traj,"i = %d , time = %f . Only the mass center position is written.\n",step,time);
//H             #endif
//H             
//H             for(i=0;i<(*sys).nb_mol;i++){
//H                 for(j=0;j<(*sys).at_p_mol;j++){
//H                     fscanf(filei,"%s %lf %lf %lf",chain,				\
//H                     (double *) &(mol[i])+(*sys).order[j],			\
//H                     (double *) &(mol[i])+(*sys).order[j]+1,		\
//H                     (double *) &(mol[i])+(*sys).order[j]+2);
//H                     
//H                 }
//H                 H pbc_mol((*sys).cell,&(mol[i])); /*Mass center + centering*/
//H                 
//H                 #ifdef DEBUG
//H                 fprintf(file_debug,"%10d %10d %5d 1 %f %f %f\n",step,(*sys).step_max,i,mol[i].xO ,mol[i].yO ,mol[i].zO );
//H                 fprintf(file_debug,"%10d %10d %5d 2 %f %f %f\n",step,(*sys).step_max,i,mol[i].xH1,mol[i].yH1,mol[i].zH1);
//H                 fprintf(file_debug,"%10d %10d %5d 3 %f %f %f\n",step,(*sys).step_max,i,mol[i].xH2,mol[i].yH2,mol[i].zH2);
//H                 #endif
//H             }
            
            
            /*=====================================================
             Storage of th*e coordinates into an array of double
             This array will have (sooner or later to replace the
             unneficient structure "mol")
             ==================================================*/
            //H for(i=0;i<(*sys).nb_mol;i++){
            //H   at_coord[9*i+0]=mol[i].xO ; at_coord[9*i+1]=mol[i].yO ; at_coord[9*i+2]=mol[i].zO ; 
            //H   at_coord[9*i+3]=mol[i].xH1; at_coord[9*i+4]=mol[i].yH1; at_coord[9*i+5]=mol[i].zH1;
            //H   at_coord[9*i+6]=mol[i].xH2; at_coord[9*i+7]=mol[i].yH2; at_coord[9*i+8]=mol[i].zH2;
            //H }
            
            //at_coord = &(topology.coordsAsVector()[0]);
            at_coord = topology.coordsAsVector().data();
            
            /*================================================================================
             Calculation of* the complete dipole moment / polarizability (permanent + induced)
             ================================================================================*/
            val_init(sys,all_mols,dip0,dip0_mol,pol0,pol0_mol,file_traj,file_dip0,file_pol0,step);/*dip0, pol0*/
            comp_dip_pol(sys,all_mols,at_coord,dip0,dip0_mol,dip,dip_mol,pol0,pol0_mol,pol,pol_mol,e_ind,a_ind,v3_tmp,m3_tmp,Tij_mc,Tij_at,x,y,z,r,fthole,fileo,file_Tij_mc,file_Tij_at,file_dip,file_pol,file_dip_ind,file_pol_ind,step);/*Induction*/
            
            
            
            (*sys).step_max++; /* Real number of steps recorded (can be lower than stepf-stepi+1) */
        }//end loop over trajectory
        else{/*Passing the lines before stepi*/
            for(i=1;i<=(*sys).nb_at;i++){
                fgets(chain,CHAIN_SIZE,filei);
            }
        }
    //} //NOTE found extra brace here, check if something is wrong
    
    /* Release the memory */
    //free(mol);
    free(at_coord);
    free(dip0);
    free(dip);
    free(dip0_mol);
    free(dip_mol);
    free(pol0);
    free(pol);
    free(pol0_mol);
    free(pol_mol);
    free(e_ind);
    free(v3_tmp);
    free(m3_tmp);
    free(Tij_mc);
    free(x);
    free(y);
    free(z);
    free(r);
    free(fthole);
    if((*sys).typ_dip + (*sys).typ_pol){free(Tij_at);}
    
    fclose(filei);
    fclose(fileo);
    #ifndef YUKI
    fclose(file_traj);
    #endif
    #ifdef DEBUG
    fclose(file_debug);
    fclose(file_dip0);
    fclose(file_pol0);
    fclose(file_Tij_mc);
    fclose(file_dip);
    fclose(file_pol);
    fclose(file_dip_ind);
    if((*sys).typ_dip + (*sys).typ_pol){fclose(file_Tij_at);}
    #endif
    
    
    
    return;
}

