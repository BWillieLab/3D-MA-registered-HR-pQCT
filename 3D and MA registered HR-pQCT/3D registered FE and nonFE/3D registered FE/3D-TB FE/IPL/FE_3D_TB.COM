$!      
$!         _/_/_/  _/_/_/    _/        
$!          _/    _/    _/  _/           Image Processing Language  
$!         _/    _/_/_/    _/ 
$!        _/    _/        _/             (c) SCANCO Medical
$!     _/_/_/  _/        _/_/_/_/				     
$!
$!========================================================================================
$! Developed by Seyedmahdi(Mahdi) Hosseinitabatabaei (mahdi.tabatabaei@mail.mcgill.ca; https://www.researchgate.net/profile/Seyedmahdi-Hosseinitabatabaei)
$! Dr. B.Willie's lab, Shriners Hospital for Children - McGill University, Montreal, Canada
$! Version P1 (07-DEC-2021)
$!========================================================================================
$! Introduction:
$! Another approach to eliminate interpolation error for 3D-registered microFE is to transform the microFE boundary conditions (i.e., compressive loading vector)
$! instead of the images after finding the common volume using 3D-registration, since the loading vectors can be transformed without error.
$! This method is named 3D-TB.
$! For the 3D-TB method, instead of rotating the moving image to align it with the reference image, the compressive displacement vector was inverse transformed
$! to match the orientation of the moving image using the same transformation matrix obtained during 3D-registration
$! (i.e., the resulting vector had components in x, y, and z directions, whereas the initial vector only had z component).
$! This way, we prevented introducing interpolation error into the images. To flatten the volumes, we created a flat common mask in the domain of the reference image,
$! and inverse transformed it to the domain of the moving image. Then, to fill the cropped regions, we added a layer of stiff material with an elastic modulus of 410 GPa and Poisson’s ratio of 0.3.
$! The addition of the stiff material enabled the application of the boundary conditions in IPL.
$! Next, during post-processing in python, reaction forces were projected onto the unit vector representing the direction of applied compressive load and aggregated to obtain the total reaction force.
$! Finally, stiffness and failure load were calculated.
$!========================================================================================
$! Important Parameters:
$!
$!========================================================================================
$! This script:
$!				1. Creatses the flattened common volume in the domain of the 1st image, and inverse transform it to the domain of the 2nd image
$!              2. Reads the segmented images and remove the non matching regions
$!              3. Fills the non matching regions with a stiff material
$!              4. Perform FE on the common region by applying the loading to the stiff material
$
$! INPUTS:
$!              a. The grayscale ".AIM" file of the 1st image
$!              b. The grayscale ".AIM" file of the 2nd image
$!              c. The top slice number (indexed from 1) to be flattened (the slice number at which the flat region starts)
$!              d. The top slice number (indexed from 1) to be flattened (the slice number at which the flat region ends)
$
$! OUTPUTS:
$!              A. Segmented image of the 1st image cropped by the common region, flattened, and filled with stiff material
$!              B. Segmented image of the 1st image cropped by the common region, flattened, and filled with stiff material
$!              C. FE LISTING and POSTLIST files containing all the information (nodal reaction forces, and elemental energy equivalent strains)
$
$! Version history:
$!				N/A
$! call as:
$! @FE_3D_TB.com    c0000001.aim    c0000002.aim		TOP		BOTTOM
$!========================================================================================
$ FILENAME 	 = F$PARSE("''P1'",,,"NAME") - ""
$ DEV      		= F$PARSE("''P1'",,,"DEVICE")
$ DIR_READ   	= F$PARSE("''P1'",,,"DIRECTORY")
$ DIR1   		:== 'DEV''DIR_READ''FILENAME'
$
$ FILENAME2 	 = F$PARSE("''P2'",,,"NAME") - ""
$ DEV2      		= F$PARSE("''P2'",,,"DEVICE")
$ DIR_READ2   	= F$PARSE("''P2'",,,"DIRECTORY")
$ DIR2   		:== 'DEV2''DIR_READ2''FILENAME2'
$
$ ! 3D-reg inputs ===================
$ IPL_GOBJ          :=  'DEV''DIR_READ''FILENAME'.GOBJ
$ IPL_GOBJ2          :=  'DEV2''DIR_READ2''FILENAME2'.GOBJ
$ IPL_GOBJ_3D         :=  'DEV''DIR_READ''FILENAME'_3D.GOBJ
$ IPL_GOBJ_3D2          :=  'DEV2''DIR_READ2''FILENAME2'_3D.GOBJ
$ IPL_FNAME		     :=  'DEV''DIR_READ''FILENAME'_3D_GOBJ.aim
$ IPL_FNAME2		     :=  'DEV2''DIR_READ2''FILENAME2'_3D_GOBJ.aim
$ IPL_SEG		     :=  'DEV''DIR_READ''FILENAME'_SEG.aim
$ IPL_SEG2		     :=  'DEV2''DIR_READ2''FILENAME2'_SEG.aim
$ IPL_FNAMEC		     :=  'DEV''DIR_READ''FILENAME'_3D_GOBJ_C.aim
$ IPL_FNAMEC2		     :=  'DEV2''DIR_READ2''FILENAME2'_3D_GOBJ_C.aim
$ IPL_SEGC		     :=  'DEV''DIR_READ''FILENAME'_SEG_C.aim
$ IPL_SEGC2		     :=  'DEV2''DIR_READ2''FILENAME2'_SEG_C.aim
$ IPL_FE		     :=  'DEV''DIR_READ''FILENAME'_3D_TB_IPL_FLAT.aim
$ IPL_FE2		     :=  'DEV2''DIR_READ2''FILENAME2'_3D_TB_IPL_FLAT.aim
$
$ define comm_file 'DIR1'"_2COMMON1_BLK.AIM"
$ define tmat_file 'DEV2''DIR_READ2''FILENAME'_'FILENAME2'_TMAT.DAT
$
$ top1=F$INTEGER(P3)
$ bottom1=F$INTEGER(P4)
$
$ top=top1-1
$ bottom=bottom1-1
$ dim= bottom-top 
$
$ XCTNO1 = F$EXTRACT(15,8,DIR_READ)
$ MEASNO1 = F$EXTRACT(24,8,DIR_READ)
$ XCTNO2 = F$EXTRACT(15,8,DIR_READ2)
$ MEASNO2 = F$EXTRACT(24,8,DIR_READ2)
$
$  JOBFILE1 :== "DISK1:"'[RESEARCH.RESULTS]'XCTNO1'_'MEASNO1'_FE_3D_TB'
$  JOBFILE2 :== "DISK1:"'[RESEARCH.RESULTS]'XCTNO2'_'MEASNO2'_FE_3D_TB'
$
$ IF F$INTEGER(P5) .EQ. 2
$ THEN
$		full_dim = 168
$		E_mod = 8748
$ ELSE
$		full_dim = 110
$		E_mod = 6829
$ ENDIF
$
$ IPL_BATCH


! First, convert the full gobjs into masks of unregistered scans

/gobj_to_aim
  -gobj_filename             "IPL_GOBJ
  -output                    mask_full
  -peel_iter                 0
  

/gobj_to_aim
  -gobj_filename             "IPL_GOBJ2
  -output                    mask_full2
  -peel_iter                 0
  
  
! Convert the 3D-registered gobj into registered masks

/gobj_to_aim
  -gobj_filename             "IPL_GOBJ_3D
  -output                    mask_3d
  -peel_iter                 0
  

/gobj_to_aim
  -gobj_filename             "IPL_GOBJ_3D2
  -output                    mask_3d2
  -peel_iter                 0



! Flatten the common volume

! ------ read common aim, clean it and create a cortex mask file ---------
/read
  -name                      comm
  -filename                  comm_file

/set_value
  -input                     comm
  -value_object              1
  -value_background          0
  
/sub_get
  -input                     comm
  -output                    comm_flat
  -pos                       -1 -1 "top
  -dim                       -1 -1 "dim
  -global_pos_flag           false
  
! Identify the stiff material region for each image
/transform
	-in 					comm_flat
	-out 					comm_flat2
	-Tmat_file_name 				tmat_file
	-img_interpol_option 				0         ! Nearest neighbor for binary masks to prevent dilation due to inteprolation
	-el_size_mm_out 				-1 -1 -1
	-Tmat_invert				true




/multiply_volumes
  -input1                    mask_full
  -input2                    comm_flat
  -output                    mask_flat
  -common_region_only        false
  -multiply_zero_pixels      true
  
/multiply_volumes
  -input1                    mask_full2
  -input2                    comm_flat2
  -output                    mask_flat2
  -common_region_only        false
  -multiply_zero_pixels      true




! Subtract the registered and unregistered masks to find the region that needs to be converted to a stiff material on the top and bottom

/subtract_aims
  -input1                    mask_full
  -input2                    mask_flat
  -output                    stiff
  
/subtract_aims
  -input1                    mask_full2
  -input2                    mask_flat2
  -output                    stiff2


! Set value of the stiff material to 100

/set_value
  -input                    stiff
  -value_object             100
  -value_background         0
  
/set_value
  -input                    stiff2
  -value_object             100
  -value_background         0

! Crop the segmented images using the 3D-registered mask

/read
  -name                      seg
  -filename                  "IPL_SEG

/cl_ow_rank_extract
	-input_output			seg
	-first_rank				1
	-last_rank				1
	-connect_boundary		false
	-value_in_range			127
	
/gobj_maskaimpeel_ow
  -input_output             seg
  -gobj_filename             "IPL_GOBJ_3D
  -peel_iter                 0
  
  
  
/read
  -name                      seg2
  -filename                  "IPL_SEG2

/cl_ow_rank_extract
	-input_output			seg2
	-first_rank				1
	-last_rank				1
	-connect_boundary		false
	-value_in_range			127
	
/gobj_maskaimpeel_ow
  -input_output             seg2
  -gobj_filename             "IPL_GOBJ_3D2
  -peel_iter                 0


  

/multiply_volumes
  -input1                    seg
  -input2                    comm_flat
  -output                    seg_comm
  -common_region_only        false
  -multiply_zero_pixels      true
  
  
/multiply_volumes
  -input1                    seg2
  -input2                    comm_flat2
  -output                    seg_comm2
  -common_region_only        false
  -multiply_zero_pixels      true


/rename
  -old_name                  seg_comm
  -new_name                  seg
  
/rename
  -old_name                  seg_comm2
  -new_name                  seg2


/set_value
  -input                    seg
  -value_object             127
  -value_background         0
  
  
/set_value
  -input                    seg2
  -value_object             127
  -value_background         0  
  
! Counting the number of voxels in the stiff material and segmented images without the stiff material to update the critical volume for failure load calculation
! This updated critical value compensates for the number of voxels added from the stiff material
/vox_scanco_param
  -input                     stiff
  -region_number             0

/vox_scanco_param
  -input                     stiff2
  -region_number             0
  
/vox_scanco_param
  -input                     seg
  -region_number             0

/vox_scanco_param
  -input                     seg2
  -region_number             0
  
  
! Attach the stiff material to the top and bottom of the segmented bone

/concat
  -input1                    seg
  -input2                    stiff
  -output                    fe_input1
  -common_region_only        false
  -add_not_overlay           false
  -make_edge                 false
  -shift_ofin2               0 0 0
  -turnangle                 0.000000
  -turnpoint_global          -1 -1
  
  
/concat
  -input1                    seg2
  -input2                    stiff2
  -output                    fe_input2
  -common_region_only        false
  -add_not_overlay           false
  -make_edge                 false
  -shift_ofin2               0 0 0
  -turnangle                 0.000000
  -turnpoint_global          -1 -1

/bounding_box_cut
  -input                    fe_input1
  -output                   fe_input1t
  -z_only                   false
  -border                   0 0 0


/bounding_box_cut
  -input                    fe_input2
  -output                   fe_input2t
  -z_only                   false
  -border                   0 0 0

! Write the inputs for FE

/write_v020
  -name                      fe_input1t
  -filename                  "IPL_FE
  -compress_type             bin
  -version_020               false 
  
/write_v020
  -name                      fe_input2t
  -filename                  "IPL_FE2
  -compress_type             bin
  -version_020               false 
..
$! EXIT
$
$IPLFE_BATCH
/read
  -name                      fe_input1
  -filename                  "IPL_FE


/read
  -name                      fe_input2
  -filename                  "IPL_FE2
!!========================================================================================
!! FE-Analysis
!!========================================================================================
!
! Material 1 = Trabecular Bone
! Material 2 = Cortical   Bone
! Material 3 (100) = Stiff material

/fe_solve3
  -in                        fe_input1
  -fea_file_name             "JOBFILE1
  -problem_nr                110
  -scale_factor              1
  -list_option               4
  -tolerance_force           1.000000E-03
  -tolerance_displ           1.000000E-03
  -max_nr_iter               20000
  -restart_option            0
  -comp_val_mat_001          126
  -Ymodulus_mat_001          "E_mod
  -Poissonr_mat_001          3.000000E-01
  -comp_val_mat_002          127
  -Ymodulus_mat_002          "E_mod
  -Poissonr_mat_002          3.000000E-01
  -comp_val_mat_003          100
  -Ymodulus_mat_003          4.100000E+05
  -Poissonr_mat_003          3.000000E-01

!
!!========================================================================================
!! FE Post-processing
!!========================================================================================
!
! VARIABLE_NR=42 only to create Von Mises stress distribution 3D output
! Critical Volume = 70 MPa
! No results written to the database here

! VARIABLE_NR=43 to determine Failure Load after Pistoia´s criterion
! Critical Volume = 0.007 = 7000 microstrain
! Results written to the database
!/db_scanco_activate
!  -write                     true

/fe_post
  -post_file_name            "JOBFILE1
  -output                    out_not_written
  -variable_nr               43
  -loadcase_nr               1
  -list_option               4
  -interpol_option           0
  -averaging_option          0
  -test_dir                  0
  -critical_volume           -0.02
  -critical_value            0.007
  -defplot_scaling           0.000000E+00
  
!/write_v020
!  -name                      OUT1
!  -filename                  "OUT1_FILE 
!  -compress_type             none
!  -version_020               false  
  
  
/fe_solve3
  -in                        fe_input2
  -fea_file_name             "JOBFILE2
  -problem_nr                110
  -scale_factor              1
  -list_option               4
  -tolerance_force           1.000000E-03
  -tolerance_displ           1.000000E-03
  -max_nr_iter               20000
  -restart_option            0
  -comp_val_mat_001          126
  -Ymodulus_mat_001          "E_mod
  -Poissonr_mat_001          3.000000E-01
  -comp_val_mat_002          127
  -Ymodulus_mat_002          "E_mod
  -Poissonr_mat_002          3.000000E-01
  -comp_val_mat_003          100
  -Ymodulus_mat_003          4.100000E+05
  -Poissonr_mat_003          3.000000E-01
  
  
/fe_post
  -post_file_name            "JOBFILE2
  -output                    out_not_written
  -variable_nr               43
  -loadcase_nr               1
  -list_option               4
  -interpol_option           0
  -averaging_option          0
  -test_dir                  0
  -critical_volume           -0.02
  -critical_value            0.007
  -defplot_scaling           0.000000E+00
  
..
$ DIRECT :== 'DEV''DIR_READ'
$ SET DEFAULT 'DIRECT'
$ CREATE CHECK_TB.TXT
$ EXIT