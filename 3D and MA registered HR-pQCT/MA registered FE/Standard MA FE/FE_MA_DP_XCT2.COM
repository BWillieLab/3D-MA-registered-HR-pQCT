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
$! Matched-angle (MA) image registration is an alternative approach to the standard 3D registrationto to reduce interpolation error.
$! Image interpolation acts as a low pass filter and blurs (widens) the image.
$! Using MA registration, the moving and reference images are both rotated by the same angle in opposite directions,
$! and the error due to image rotation is propagated to both images, hence reducing the bias in interpolation error.
$!
$! The current script the standard version of MA-registered microFE where each scan is segmented using its own rotated mask is used.
$!
$! For the MA-SP method, we rotated the reference and moving images using rotation angles of half of the original angles, but in opposite directions.
$! This way, the MA method aligns the two images in a middle domain. Briefly, using the registration function in IPL,
$! we divided the rotation angles from registration by half to identify the middle domain.
$! Then, each of the reference and moving images where registered to the middle domain. Finally, we performed microFE on the flat common volume.
$!
$! * The MA registration procedure is developed bsaed on scripts obtained from Dr. Bert van Rietbergen
$!========================================================================================
$! Important Parameters:
$!
$!========================================================================================
$! This script:
$!				1. Aligns the grayscale images of two scans from XCT device in the middle domain
$!              2. Segmentes the images in the common domain
$!              3. Masks the images using a dilated common mask (same periosteal mask: SP)
$!              4. Flattens the segmented images at the top and bottom surfaces for applying FE boundary conditions
$
$! INPUTS:
$!              a. The grayscale ".AIM" file of the 1st image
$!              b. The grayscale ".AIM" file of the 2nd image
$!              c. The top slice number (indexed from 0) to be flattened (the slice number at which the flat region starts)
$!              d. The top slice number (indexed from 0) to be flattened (the slice number at which the flat region ends)
$
$! OUTPUTS:
$!              A. Grayscale, segmented, cortical mask, and trabecular mask of the 1st image cropped by the common region
$!              B. Grayscale, segmented, cortical mask, and trabecular mask of the 2nd image cropped by the common region
$!              C. Periosteal, endosteal, and cortical GOBJs of the 1st image cropped by the common region
$!              D. Periosteal, endosteal, and cortical GOBJs of the 2nd image cropped by the common region
$!              E. Grayscale image and periosteal GOBJ of the reference image in the middle domain
$!              F. The segmented image of the 1st image cropped by the common region and flattened
$!              G. The segmented image of the 2nd image cropped by the common region and flattened
$!              H. An overlay of the segmented images in the middle domain
$!              I. The transformation matrix aligning the 1st image onto the middle domain (initial transformation acquired from registering 2nd image to 1st) (.DAT file)
$!              J. The transformation matrix aligning the 1st image onto the middle domain (.DAT file)
$!              K. The transformation matrix aligning the 2nd image onto the middle domain (.DAT file)
$
$! Version history:
$!				N/A
$! call as:
$! @FE_MA_DP_XCT2.com    c0000001.aim    c0000002.aim		TOP		BOTTOM
$!========================================================================================
$! Defining the files and directories for files
$ FILENAME 	 = F$PARSE("''p1'",,,"NAME") - ""
$ DEV      		= F$PARSE("''p1'",,,"DEVICE")
$ DIR_READ   	= F$PARSE("''p1'",,,"DIRECTORY")
$ FILENAME2  	= F$PARSE("''p2'",,,"NAME") - ""
$ DEV2       		= F$PARSE("''p2'",,,"DEVICE")
$ DIR_READ2  	= F$PARSE("''p2'",,,"DIRECTORY")
$ DIR1   		:== 'DEV''DIR_READ''FILENAME'
$ DIR2   		:== 'DEV2''DIR_READ2''FILENAME2'
$ JOBNAME   		:== 'DEV2''DIR_READ2''FILENAME2'
$ TRANSFORMATIONBL    	:== 'DEV''DIR_READ''FILENAME'_'FILENAME2'_TMATBL.DAT
$ TRANSFORMATIONFU    	:== 'DEV2''DIR_READ2''FILENAME'_'FILENAME2'_TMATFU.DAT
$
$ BL_TMAT_FILE    	:== 'DEV2''DIR_READ2''FILENAME'_'FILENAME2'_MA_TMAT.DAT
$
$ MSKT1   		:== 'DEV''DIR_READ''FILENAME'_MA_MASKT1.AIM
$ MSKT2   		:== 'DEV2''DIR_READ2''FILENAME2'_MA_MASKT2.AIM
$ CTMSKT1   		:== 'DEV''DIR_READ''FILENAME'_MA_CTMASKT1.AIM
$ CTMSKT2   		:== 'DEV2''DIR_READ2''FILENAME2'_MA_CTMASKT2.AIM
$ TBMSKT1   		:== 'DEV''DIR_READ''FILENAME'_MA_TBMASKT1.AIM
$ TBMSKT2   		:== 'DEV2''DIR_READ2''FILENAME2'_MA_TBMASKT2.AIM
$ SEG1   		:== 'DEV''DIR_READ''FILENAME'_MA_SEG1.AIM
$ SEG2   		:== 'DEV2''DIR_READ2''FILENAME2'_MA_SEG2.AIM
$ OVERLAY1_2   		:== 'DEV''DIR_READ''FILENAME'_MA_OVERLAYED.AIM
$
$ FE_FLAT_1   		:== 'DEV''DIR_READ''FILENAME'_FE_MA_DP_FLAT.AIM
$ FE_FLAT_2   		:== 'DEV2''DIR_READ2''FILENAME2'_FE_MA_DP_FLAT.AIM
$
$ GST_REF 		:== 'DEV''DIR_READ''FILENAME'_MA_REF.AIM
$ GOBJT_REF 		:== 'DEV''DIR_READ''FILENAME'_MA_REF.GOBJ
$
$ GST1 		:== 'DEV''DIR_READ''FILENAME'_MA_GST1.AIM
$ GST2 		:== 'DEV2''DIR_READ2''FILENAME2'_MA_GST2.AIM
$
$ GOBJT1           :== 'DEV''DIR_READ''FILENAME'_MA.GOBJ
$ GOBJT2           :== 'DEV2''DIR_READ2''FILENAME2'_MA.GOBJ
$ CTGOBJT1           :== 'DEV''DIR_READ''FILENAME'_CTMA.GOBJ
$ CTGOBJT2           :== 'DEV2''DIR_READ2''FILENAME2'_CTMA.GOBJ
$ TBGOBJT1           :== 'DEV''DIR_READ''FILENAME'_TBMA.GOBJ
$ TBGOBJT2           :== 'DEV2''DIR_READ2''FILENAME2'_TBMA.GOBJ
$
$ 	define gobj1_file 'DIR1'".GOBJ"
$ 	define gobj2_file 'DIR2'".GOBJ"
$
$ 	define ctmask1_file 'DIR1'"_CORT_MASK.AIM" 	! This is different between XCT1 (_BLCK) & XCT2 (_CORT_MASK)
$ 	define ctmask2_file 'DIR2'"_CORT_MASK.AIM"
$
$ 	define tbmask1_file 'DIR1'"_TRAB_MASK.AIM"	! This is different between XCT1 (_TRAB_BLCK) & XCT2 (_TRAB_MASK)
$ 	define tbmask2_file 'DIR2'"_TRAB_MASK.AIM"
$
$ top=F$INTEGER(P3)
$ bottom=F$INTEGER(P4)
$
$ dim= bottom-top 
$!
$ 	ipl_scanco_prog := $um:ipl_scanco_m.exe
$!
$ ON ERROR THEN EXIT
$!
$ 	IPL_BATCH

! Read the input files
/read
	-name 					in1
	-filename					"P1
		
/read
	-name 					in2
	-filename 					"P2

/read
  -name                      ctmsk1
  -filename                  ctmask1_file

/set_value
  -input                     ctmsk1
  -value_object              1
  -value_background          0

  
/read
  -name                      tbmsk1
  -filename                  tbmask1_file
  
/set_value
  -input                     tbmsk1
  -value_object              1
  -value_background          0
  
  
/concat
	-input1 					ctmsk1
	-input2 					tbmsk1
	-output 					msk1
	-common_region_only 			false
	-add_not_overlay 				true
	-make_edge 				false
	-shift_ofin2 				0 0 0
	-turnangle 					0.0
	-turnpoint_global 				-1 -1

/set_value
	-input					msk1
	-value_object				127
	-value_background				0	
	
/read
  -name                      ctmsk2
  -filename                  ctmask2_file

/set_value
  -input                     ctmsk2
  -value_object              1
  -value_background          0

  
/read
  -name                      tbmsk2
  -filename                  tbmask2_file
  
/set_value
  -input                     tbmsk2
  -value_object              1
  -value_background          0
  
/concat
	-input1 					ctmsk2
	-input2 					tbmsk2
	-output 					msk2
	-common_region_only 			false
	-add_not_overlay 				true
	-make_edge 				false
	-shift_ofin2 				0 0 0
	-turnangle 					0.0
	-turnpoint_global 				-1 -1

/set_value
	-input					msk2
	-value_object				127
	-value_background				0	
	
! First, find the transformation matrix that aligns the baseline image to the followup image, divide it by 2 so that it can be used to map everything to the middle domain	
	
/register
  -in1                        in2
  -gobj_filename_in1          gobj2_file
  -in2                        in1
  -gobj_filename_in2          gobj1_file
  -Tmat_file_name             "BL_TMAT_FILE
  -orientation_search         1
  -initial_rotation           0.0  0.0  0.0
  -initial_translation        0.0  0.0  0.0
  -delta_rotation             0.10 0.10 0.50
  -delta_translation          0.10 0.10 0.50
  -resolution_scaling         8 4 2                  ! no need for accurate registration here
  -delta_scaling              1.00 0.10 0.10
  -tolerance                  0.00001
  -min_corr_coef              0.5
  -min_method                 1
  -object_func                1
  -max_nr_iter                1000
  -scale_factor               0.5	
	

! Now that the transformation for the middle domain is found, transform the baseline image to that domain	

/transform
  -in                        in1
  -out                       in1_ref
  -Tmat_file_name            "BL_TMAT_FILE
  -img_interpol_option       2						
  -el_size_mm_out            -1.000 -1.000 -1.000

! Before saving the image, move the image in the middle domain by 500 voxels in z direction so that any negative coordinates will be prevented
/header_geo_set
  -input                     in1_ref
  -off_new                   -1 -1 -1
  -pos_new                   -1 -1 500
  -el_size_mm_new            -1.000 -1.000 -1.000
  

/write_v020
  -name                      in1_ref
  -filename                  "GST_REF
  -compress_type             none
  -version_020               false	


! Transform the mask of the baseline image to the middle domain to create a gobj for the future registrations
/transform
  -in                        msk1
  -out                       tmp
  -Tmat_file_name            "BL_TMAT_FILE
  -img_interpol_option       0
  -el_size_mm_out            -1.000 -1.000 -1.000

! Before saving the gobj, move the image in the middle domain by 500 voxels in z direction so that any negative coordinates will be prevented
/header_geo_set
  -input                     tmp
  -off_new                   -1 -1 -1
  -pos_new                   -1 -1 500
  -el_size_mm_new            -1.000 -1.000 -1.000


/togobj_from_aim
  -input                     tmp
  -gobj_filename             "GOBJT_REF
  -min_elements              0
  -max_elements              0
  -curvature_smooth          1	
	

	
	
! Now, register the baseline image to the middle domain to incorporate the 500 slice offset
/register
	-in1 					in1_ref
	-gobj_filename_in1 				"GOBJT_REF
	-in2 					in1
	-gobj_filename_in2 				gobj1_file
	-Tmat_file_name 				"TRANSFORMATIONBL
	-orientation_search 				1
	-initial_rotation 				0.0 0.0 0.0
	-initial_translation 				0.0 0.0 0.0
	-delta_rotation 				0.30 0.30 0.50
	-delta_translation 				0.30 0.30 0.50
	-resolution_scaling 				10 4 1
	-delta_scaling				1.00 0.10 0.10
	-weight_gobj_overlap 			0.20 0.00 0.00
	-tolerance 					0.00001
	-min_corr_coef 				0.5
	-min_method 				1
	-object_func 				1
	-max_nr_iter 				1000
	

! Now, transfer all the grayscale image and masks of the baseline to the middle domain. Same for the solid block
/transform
  -in                        in1
  -out                       int1
  -Tmat_file_name            "TRANSFORMATIONBL
  -img_interpol_option       2						
  -el_size_mm_out            -1.000 -1.000 -1.000
  

/write_v020
  -name                      int1
  -filename                  "GST1
  -compress_type             none
  -version_020               false	


! Transform the mask of the baseline image to the middle domain to create a gobj for the future registrations
/transform
  -in                        msk1
  -out                       mskt1
  -Tmat_file_name            "TRANSFORMATIONBL
  -img_interpol_option       0
  -el_size_mm_out            -1.000 -1.000 -1.000


/togobj_from_aim
  -input                     mskt1
  -gobj_filename             "GOBJT1
  -min_elements              0
  -max_elements              0
  -curvature_smooth          1



/transform
  -in                        tbmsk1
  -out                       tbmskt1
  -Tmat_file_name            "TRANSFORMATIONBL
  -img_interpol_option       0
  -el_size_mm_out            -1.000 -1.000 -1.000


/togobj_from_aim
  -input                     tbmskt1
  -gobj_filename             "TBGOBJT1
  -min_elements              0
  -max_elements              0
  -curvature_smooth          1



/transform
  -in                        ctmsk1
  -out                       ctmskt1
  -Tmat_file_name            "TRANSFORMATIONBL
  -img_interpol_option       0
  -el_size_mm_out            -1.000 -1.000 -1.000


/togobj_from_aim
  -input                     ctmskt1
  -gobj_filename             "CTGOBJT1
  -min_elements              0
  -max_elements              0
  -curvature_smooth          1	  
	



! Now create a solid block from each aim file to find the common region
/norm_max
	-input					in1
	-output					blk1
	-max					127
	-type_out				char
	
/delete
	-name					in1


/set_value
	-input					blk1
	-value_object			5
	-value_background		0

/transform
  -in                        blk1
  -out                       blkt1
  -Tmat_file_name            "TRANSFORMATIONBL
  -img_interpol_option       0
  -el_size_mm_out            -1.000 -1.000 -1.000


	




! Now, register the followup image to the middle domain to incorporate the 500 slice offset
/register
	-in1 					in1_ref
	-gobj_filename_in1 				"GOBJT_REF
	-in2 					in2
	-gobj_filename_in2 				gobj2_file
	-Tmat_file_name 				"TRANSFORMATIONFU
	-orientation_search 				1
	-initial_rotation 				0.0 0.0 0.0
	-initial_translation 				0.0 0.0 0.0
	-delta_rotation 				0.30 0.30 0.50
	-delta_translation 				0.30 0.30 0.50
	-resolution_scaling 				10 4 1
	-delta_scaling				1.00 0.10 0.10
	-weight_gobj_overlap 			0.20 0.00 0.00
	-tolerance 					0.00001
	-min_corr_coef 				0.5
	-min_method 				1
	-object_func 				1
	-max_nr_iter 				1000


/delete
	-name					in1_ref	

! Now, transfer all the grayscale image and masks of the followup to the middle domain. Same for the solid block
/transform
  -in                        in2
  -out                       int2
  -Tmat_file_name            "TRANSFORMATIONFU
  -img_interpol_option       2						
  -el_size_mm_out            -1.000 -1.000 -1.000
  

/write_v020
  -name                      int2
  -filename                  "GST2
  -compress_type             none
  -version_020               false	


! Transform the mask of the followup image to the middle domain to create a gobj for the future registrations
/transform
  -in                        msk2
  -out                       mskt2
  -Tmat_file_name            "TRANSFORMATIONFU
  -img_interpol_option       0
  -el_size_mm_out            -1.000 -1.000 -1.000


/togobj_from_aim
  -input                     mskt2
  -gobj_filename             "GOBJT2
  -min_elements              0
  -max_elements              0
  -curvature_smooth          1



/transform
  -in                        tbmsk2
  -out                       tbmskt2
  -Tmat_file_name            "TRANSFORMATIONFU
  -img_interpol_option       0
  -el_size_mm_out            -1.000 -1.000 -1.000


/togobj_from_aim
  -input                     tbmskt2
  -gobj_filename             "TBGOBJT2
  -min_elements              0
  -max_elements              0
  -curvature_smooth          1



/transform
  -in                        ctmsk2
  -out                       ctmskt2
  -Tmat_file_name            "TRANSFORMATIONFU
  -img_interpol_option       0
  -el_size_mm_out            -1.000 -1.000 -1.000


/togobj_from_aim
  -input                     ctmskt2
  -gobj_filename             "CTGOBJT2
  -min_elements              0
  -max_elements              0
  -curvature_smooth          1	  
	



! Now create a solid block from each aim file to find the common region
/norm_max
	-input					in2
	-output					blk2
	-max					127
	-type_out				char
	
/delete
	-name					in2

/set_value
	-input					blk2
	-value_object			10
	-value_background		10

/transform
  -in                        blk2
  -out                       blkt2
  -Tmat_file_name            "TRANSFORMATIONFU
  -img_interpol_option       0
  -el_size_mm_out            -1.000 -1.000 -1.000



	
! =============================================================================	
! Now, the solid blocks are added and thresholded to isolate the common region

! Add the segmented images; common regions will get a value of 3
/concat
	-input1 					blkt1
	-input2 					blkt2
	-output 					maskcomb
	-common_region_only 			false
	-add_not_overlay 				true
	-make_edge 				false
	-shift_ofin2 				0 0 0
	-turnangle 					0.0
	-turnpoint_global 				-1 -1


/delete
	-name					blkt1
	
/delete
	-name					blkt2
		
/threshold
	-input					maskcomb
	-output					comm
	-lower_in_perm_aut_al				15
	-upper_in_perm_aut_al				15
	-value_in_range				127
	-unit					5

/delete
	-name					maskcomb
	
	
/set_value
	-input					comm
	-value_object				1
	-value_background				0

! Crop the mask using the common mask  
/multiply_volumes
  -input1                    ctmskt1
  -input2                    comm
  -output                    ctmskb1
  -common_region_only        false
  -multiply_zero_pixels      true

/set_value
  -input                     ctmskb1
  -value_object              1
  -value_background          0

/rename
  -old_name                  ctmskb1
  -new_name                  ctmskt1 
  
/multiply_volumes
  -input1                    tbmskt1
  -input2                    comm
  -output                    tbmskb1
  -common_region_only        false
  -multiply_zero_pixels      true

/set_value
  -input                     tbmskb1
  -value_object              1
  -value_background          0

/rename
  -old_name                  tbmskb1
  -new_name                  tbmskt1  

/multiply_volumes
  -input1                    ctmskt2
  -input2                    comm
  -output                    ctmskb2
  -common_region_only        false
  -multiply_zero_pixels      true

/set_value
  -input                     ctmskb2
  -value_object              1
  -value_background          0

/rename
  -old_name                  ctmskb2
  -new_name                  ctmskt2 
  
/multiply_volumes
  -input1                    tbmskt2
  -input2                    comm
  -output                    tbmskb2
  -common_region_only        false
  -multiply_zero_pixels      true

/delete
	-name					comm


/set_value
  -input                     tbmskb2
  -value_object              1
  -value_background          0

/rename
  -old_name                  tbmskb2
  -new_name                  tbmskt2  
  
! Segmentation

/offset_add
  -input                     int1
  -add_offset                1 1 1


/fill_offset_duplicate
  -input                     int1


/sup_divide
  -input                     int1
  -supdim_numbers            4 4 1
  -testoff_pixels            1 1 1
  -suppos_pixels_local       -1 -1 -1
  -subdim_pixels             -1 -1 -1
  
/copy
  -in                        int1
  -out                       int12

/seg_gauss
  -input                     int1
  -output                    int1_gauss
  -sigma                     0.800000
  -support                   1
  -lower_in_perm_aut_al      320.000000
  -upper_in_perm_aut_al      3000.000000
  -value_in_range            126
  -unit                      2  

/delete
	-name					int1

/seg_gauss
  -input                     int12
  -output                    int12_gauss
  -sigma                     0.800000
  -support                   1
  -lower_in_perm_aut_al      450.000000
  -upper_in_perm_aut_al      3000.000000
  -value_in_range            127
  -unit                      2  

/delete
	-name					int12
	
/multiply_volumes
  -input1                    int1_gauss
  -input2                    tbmskt1
  -output                    tb_input1
  -common_region_only        false
  -multiply_zero_pixels      true

/delete
	-name					int1_gauss
  
/multiply_volumes
  -input1                    int12_gauss
  -input2                    ctmskt1
  -output                    ct_input1
  -common_region_only        false
  -multiply_zero_pixels      true

/delete
	-name					int12_gauss

/add_aims
  -input1                    tb_input1
  -input2                    ct_input1
  -output                    comb1

/delete
	-name					tb_input1
	
/delete
	-name					ct_input1
  
! SECOND SCAN ==============================================================================  

/offset_add
  -input                     int2
  -add_offset                1 1 1


/fill_offset_duplicate
  -input                     int2


/sup_divide
  -input                     int2
  -supdim_numbers            4 4 1
  -testoff_pixels            1 1 1
  -suppos_pixels_local       -1 -1 -1
  -subdim_pixels             -1 -1 -1
  
/copy
  -in                        int2
  -out                       int22

/seg_gauss
  -input                     int2
  -output                    int2_gauss
  -sigma                     0.800000
  -support                   1
  -lower_in_perm_aut_al      320.000000
  -upper_in_perm_aut_al      3000.000000
  -value_in_range            126
  -unit                      2  

/delete
	-name					int2

/seg_gauss
  -input                     int22
  -output                    int22_gauss
  -sigma                     0.800000
  -support                   1
  -lower_in_perm_aut_al      450.000000
  -upper_in_perm_aut_al      3000.000000
  -value_in_range            127
  -unit                      2  

/delete
	-name					int22
	
/multiply_volumes
  -input1                    int2_gauss
  -input2                    tbmskt2
  -output                    tb_input2
  -common_region_only        false
  -multiply_zero_pixels      true

/delete
	-name					int2_gauss
  
/multiply_volumes
  -input1                    int22_gauss
  -input2                    ctmskt2
  -output                    ct_input2
  -common_region_only        false
  -multiply_zero_pixels      true

/delete
	-name					int22_gauss

/add_aims
  -input1                    tb_input2
  -input2                    ct_input2
  -output                    comb2

/delete
	-name					tb_input2
	
/delete
	-name					ct_input2

! ====== Now make sure both scans have the same region selected ========
/bounding_box_cut
  -input                     comb1
  -output                    comb1trim
  -z_only                    false
  -border                    0 0 0

/rename
  -old_name                  comb1trim
  -new_name                  comb1c
  
/set_value
  -input                     comb1c
  -value_object              1
  -value_background          1

/bounding_box_cut
  -input                     comb2
  -output                    comb2trim
  -z_only                    false
  -border                    0 0 0
  
/rename
  -old_name                  comb2trim
  -new_name                  comb2c

/set_value
  -input                     comb2c
  -value_object              1
  -value_background          1
  
  
! Create common coordinates by multiplying the images
/multiply_volumes
  -input1                    comb1c			
  -input2                    comb2c
  -output                    combmsk
  -common_region_only        true
  -multiply_zero_pixels      true

/delete
	-name					comb1c
	
/delete
	-name					comb2c
	
/set_value
  -input                     combmsk
  -value_object              1
  -value_background          1

! Multiply the mask in the common coordinates to match the cootdinates for the images    
/multiply_volumes
  -input1                    comb1			
  -input2                    combmsk
  -output                    comb1u
  -common_region_only        false
  -multiply_zero_pixels      true
  
/multiply_volumes
  -input1                    comb2			
  -input2                    combmsk
  -output                    comb2u
  -common_region_only        false
  -multiply_zero_pixels      true

!/delete
!	-name					combmsk

/delete
	-name					comb1
	
/delete
	-name					comb2
	
/bounding_box_cut
  -input                     comb1u
  -output                    comb1ut
  -z_only                    false
  -border                    0 0 0

/rename
  -old_name                  comb1ut
  -new_name                  comb1u 

/write_v020
  -name                      comb1u
  -filename                  "SEG1
  -compress_type             none              ! Non-binarized since may be used out of IPL
  -version_020               false 
  
! Flatten the top and bottom surfaces    
/sub_get
  -input                     comb1u
  -output                    sub_input1
  -pos                       -1 -1 "top
  -dim                       -1 -1 "dim
  -global_pos_flag           false
	
!/delete
!	-name					comb1u

/bounding_box_cut
  -input                     sub_input1
  -output                    sub_input1t
  -z_only                    false
  -border                    0 0 0

/delete
	-name					sub_input1
	
/norm_max
	-input					sub_input1t
	-output					sub_input_char1
	-max					127
	-type_out				char

/delete
	-name					sub_input1t

! component labeling to remove the unconnected regions
/cl_ow_rank_extract
	-input_output			sub_input_char1
	-first_rank				1
	-last_rank				1
	-connect_boundary		false
	-value_in_range			127
	
! Writes the CL´ed input for the solver
/copy
  -in                        sub_input_char1
  -out                       fe_input1

/delete
	-name					sub_input_char1
	
/write_v020
  -name                      fe_input1
  -filename                  "FE_FLAT_1
  -compress_type             none              ! Non-binarized since may be used out of IPL
  -version_020               false 
  
/bounding_box_cut
  -input                     comb2u
  -output                    comb2ut
  -z_only                    false
  -border                    0 0 0

/rename
  -old_name                  comb2ut
  -new_name                  comb2u 

/write_v020
  -name                      comb2u
  -filename                  "SEG2
  -compress_type             none              ! Non-binarized since may be used out of IPL
  -version_020               false 
  
/sub_get
  -input                     comb2u
  -output                    sub_input2
  -pos                       -1 -1 "top
  -dim                       -1 -1 "dim
  -global_pos_flag           false

! Create an overlay of the segmented images registered using matched angle method, to check the registration success/failure
/set_value
  -input                     comb1u
  -value_object              1
  -value_background          0
  
/set_value
  -input                     comb2u
  -value_object              2
  -value_background          0
  
/concat
	-input1 					comb1u
	-input2 					comb2u
	-output 					overlay12
	-common_region_only 			false
	-add_not_overlay 				true
	-make_edge 				false
	-shift_ofin2 				0 0 0
	-turnangle 					0.0
	-turnpoint_global 				-1 -1

/delete
	-name					comb1u

/delete
	-name					comb2u
	

/bounding_box_cut
  -input                     sub_input2
  -output                    sub_input2t
  -z_only                    false
  -border                    0 0 0

/delete
	-name					sub_input2
 
/norm_max
	-input					sub_input2t
	-output					sub_input_char2
	-max					127
	-type_out				char

/delete
	-name					sub_input2t

! component labeling to remove the unconnected regions
/cl_ow_rank_extract
	-input_output			sub_input_char2
	-first_rank				1
	-last_rank				1
	-connect_boundary		false
	-value_in_range			127
	
! Writes the CL´ed input for the solver
/copy
  -in                        sub_input_char2
  -out                       fe_input2

/delete
	-name					sub_input_char2
	
/write_v020
  -name                      fe_input2
  -filename                  "FE_FLAT_2 
  -compress_type             none          ! Non-binarized since may be used out of IPL
  -version_020               false



	
! Now, crop all the masks using the common region (combmsk12 has the same size for the first and second scans)
! Whole bone mask

/multiply_volumes
  -input1                    mskt1			
  -input2                    combmsk
  -output                    mskt1u
  -common_region_only        false
  -multiply_zero_pixels      true

/delete
	-name					mskt1

/bounding_box_cut
  -input                     mskt1
  -output                    mskt1t
  -z_only                    false
  -border                    0 0 0

/rename mskt1t mskt1
	
/norm_max
	-input					mskt1u
	-output					mskt1u_char
	-max					127
	-type_out				char

/delete
	-name					mskt1u
	
/write_v020
  -name                      mskt1u_char
  -filename                  "MSKT1
  -compress_type             none
  -version_020               false 

/delete
	-name					mskt1u_char
	
  
/multiply_volumes
  -input1                    mskt2			
  -input2                    combmsk
  -output                    mskt2u
  -common_region_only        false
  -multiply_zero_pixels      true

/delete
	-name					mskt2

/bounding_box_cut
  -input                     mskt2
  -output                    mskt2t
  -z_only                    false
  -border                    0 0 0

/rename mskt2t mskt2
	
/norm_max
	-input					mskt2u
	-output					mskt2u_char
	-max					127
	-type_out				char

/delete
	-name					mskt2u

/write_v020
  -name                      mskt2u_char
  -filename                  "MSKT2
  -compress_type             none
  -version_020               false 

/delete
	-name					mskt2u_char

	
! Cortical mask	

/multiply_volumes
  -input1                    ctmskt1			
  -input2                    combmsk
  -output                    ctmskt1u
  -common_region_only        false
  -multiply_zero_pixels      true

/delete
	-name					ctmskt1

/bounding_box_cut
  -input                     ctmskt1u
  -output                    ctmskt1ut
  -z_only                    false
  -border                    0 0 0

/rename ctmskt1ut ctmskt1u
	
/norm_max
	-input					ctmskt1u
	-output					ctmskt1u_char
	-max					127
	-type_out				char

/delete
	-name					ctmskt1u
	
/write_v020
  -name                      ctmskt1u_char
  -filename                  "CTMSKT1
  -compress_type             none
  -version_020               false 

/delete
	-name					ctmskt1u_char
	
  
/multiply_volumes
  -input1                    ctmskt2			
  -input2                    combmsk
  -output                    ctmskt2u
  -common_region_only        false
  -multiply_zero_pixels      true

/delete
	-name					ctmskt2

/bounding_box_cut
  -input                     ctmskt2u
  -output                    ctmskt2ut
  -z_only                    false
  -border                    0 0 0

/rename ctmskt2ut ctmskt2u
	
/norm_max
	-input					ctmskt2u
	-output					ctmskt2u_char
	-max					127
	-type_out				char

/delete
	-name					ctmskt2u

/write_v020
  -name                      ctmskt2u_char
  -filename                  "CTMSKT2
  -compress_type             none
  -version_020               false 

/delete
	-name					ctmskt2u_char


! Trabecular mask	

/multiply_volumes
  -input1                    tbmskt1			
  -input2                    combmsk
  -output                    tbmskt1u
  -common_region_only        false
  -multiply_zero_pixels      true

/delete
	-name					tbmskt1

/bounding_box_cut
  -input                     tbmskt1u
  -output                    tbmskt1ut
  -z_only                    false
  -border                    0 0 0

/rename tbmskt1ut tbmskt1u
	
/norm_max
	-input					tbmskt1u
	-output					tbmskt1u_char
	-max					127
	-type_out				char

/delete
	-name					tbmskt1u
	
/write_v020
  -name                      tbmskt1u_char
  -filename                  "TBMSKT1
  -compress_type             none
  -version_020               false 

/delete
	-name					tbmskt1u_char
	
  
/multiply_volumes
  -input1                    tbmskt2			
  -input2                    combmsk
  -output                    tbmskt2u
  -common_region_only        false
  -multiply_zero_pixels      true

/delete
	-name					tbmskt2

/bounding_box_cut
  -input                     tbmskt2u
  -output                    tbmskt2ut
  -z_only                    false
  -border                    0 0 0

/rename tbmskt2ut tbmskt2u	

/norm_max
	-input					tbmskt2u
	-output					tbmskt2u_char
	-max					127
	-type_out				char

/delete
	-name					tbmskt2u

/write_v020
  -name                      tbmskt2u_char
  -filename                  "TBMSKT2
  -compress_type             none
  -version_020               false 

/delete
	-name					tbmskt2u_char
..
$ EXIT
