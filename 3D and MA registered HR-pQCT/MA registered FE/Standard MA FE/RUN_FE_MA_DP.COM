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
$! This script:
$!				1. Performs FE analysis on the MA-registered regions of the two input images using the same periosteal mask approach(MA-SP method)
$
$! INPUTS:
$!              a. The grayscale ".AIM" file of the 1st image
$!              b. The grayscale ".AIM" file of the 2nd image
$!              c. Scanner generation (integer)
$
$! OUTPUTS:
$!              A. FE outcomes (AIM, LISTING, POST, POSTLIST) using Pistoia's criterion
$
$! Version history:
$
$! call as:
$! @RUN_FE_MA_DP.com    c0000001.aim    c0000002.aim		Scanner generation
$!========================================================================================
$! Important Parameters:
$!
$!========================================================================================
$
$ FILENAME 	 = F$PARSE("''P1'",,,"NAME") - ""
$ DEV      		= F$PARSE("''P1'",,,"DEVICE")
$ DIR_READ   	= F$PARSE("''P1'",,,"DIRECTORY")
$ DIR1   		:== 'DEV''DIR_READ''FILENAME'
$
$ FILENAME2 	 = F$PARSE("''P2'",,,"NAME") - ""
$ DEV2      		= F$PARSE("''P2'",,,"DEVICE")
$ DIR_READ2   	= F$PARSE("''P2'",,,"DIRECTORY")
$ DIR2   		:== 'DEV2''DIR_READ2''FILENAME2'
$ GOBJT           :== 'DEV2''DIR_READ2''FILENAME2'_TRANSFORMED.GOBJ	
$
$  JOBFILE1 := 'DIR1'_FE_MA_DP_JOB
$  JOBFILE2 := 'DIR2'_FE_MA_DP_JOB
$ FE_FLAT_1   		:== 'DIR1'_FE_MA_DP_FLAT.AIM
$ FE_FLAT_2   		:== 'DIR2'_FE_MA_DP_FLAT.AIM
$
$
$ IF P3 .EQ. 2
$ THEN
$		full_dim = 168
$		E_mod = 8748
$ ELSE
$		full_dim = 110
$		E_mod = 6829
$ ENDIF
$
$ IPLFE_BATCH

! STEP 1:
! Read the flat segmented images to run FE
/read
  -name                      aim1
  -filename                  "FE_FLAT_1
  
/read
  -name                      aim2
  -filename                  "FE_FLAT_2


! STEP 2:
! ----- FE-analysis ------------------------

/fe_solve3
  -in                        aim1
  -fea_file_name             "JOBFILE1
  -problem_nr                33
  -scale_factor              -0.01
  -list_option               1
  -tolerance_force           1.000000E-03
  -tolerance_displ           1.000000E-03
  -max_nr_iter               20000
  -restart_option            0
  -comp_val_mat_001          127
  -Ymodulus_mat_001          "E_mod
  -Poissonr_mat_001          3.000000E-01
  -comp_val_mat_002          126
  -Ymodulus_mat_002          "E_mod
  -Poissonr_mat_002          3.000000E-01
  -comp_val_mat_003          0
  -Ymodulus_mat_003          0
  -Poissonr_mat_003          3.000000E-01

! STEP 3:
! ----- FE-post process --------------------

/fe_post
  -post_file_name            "JOBFILE1
  -output                    out_not_written
  -variable_nr               43
  -loadcase_nr               1
  -list_option               1
  -interpol_option           1
  -averaging_option          0
  -test_dir                  0  
  -critical_volume           -0.02 
  -critical_value            0.007 
  -defplot_scaling           0.000000E+00
  
  
! =================== 2nd image ===============

/fe_solve3
  -in                        aim2
  -fea_file_name             "JOBFILE2
  -problem_nr                33
  -scale_factor              -0.01
  -list_option               1
  -tolerance_force           1.000000E-03
  -tolerance_displ           1.000000E-03
  -max_nr_iter               20000
  -restart_option            0
  -comp_val_mat_001          127
  -Ymodulus_mat_001          "E_mod
  -Poissonr_mat_001          3.000000E-01
  -comp_val_mat_002          126
  -Ymodulus_mat_002          "E_mod
  -Poissonr_mat_002          3.000000E-01
  -comp_val_mat_003          0
  -Ymodulus_mat_003          0
  -Poissonr_mat_003          3.000000E-01


/fe_post
  -post_file_name            "JOBFILE2
  -output                    out_not_written
  -variable_nr               43
  -loadcase_nr               1
  -list_option               1
  -interpol_option           1
  -averaging_option          0
  -test_dir                  0  
  -critical_volume           -0.02 
  -critical_value            0.007 
  -defplot_scaling           0.000000E+00
..
$ EXIT
