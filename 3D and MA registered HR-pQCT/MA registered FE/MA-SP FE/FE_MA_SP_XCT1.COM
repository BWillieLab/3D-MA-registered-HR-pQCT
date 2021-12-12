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
$! The current script performs an alternative version of MA-registered microFE, where instead of masking each image using its
$! own mask, both images are masked using a dilated mask, effectively eliminating the influence of contouring error and interpolation error of the masks. This approach is called MA-SP (same periosteal).
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
$! @FE_MA_SP_XCT1.com    c0000001.aim    c0000002.aim		TOP		BOTTOM
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
$
$ FE_FLAT_1   		:== 'DEV''DIR_READ''FILENAME'_FE_MA_SP_FLAT.AIM
$ FE_FLAT_2   		:== 'DEV2''DIR_READ2''FILENAME2'_FE_MA_SP_FLAT.AIM
$
$ GST_REF 		:== 'DEV''DIR_READ''FILENAME'_MA_REF.AIM
$ GOBJT_REF 		:== 'DEV''DIR_READ''FILENAME'_MA_REF.GOBJ
$
$ 	define gobj1_file 'DIR1'".GOBJ"
$ 	define gobj2_file 'DIR2'".GOBJ"
$
$ 	define msk1_file 'DIR1'"_BLCK.AIM" 	! This is different between XCT1 (_BLCK) & XCT2 (_CORT_MASK)
$ 	define msk2_file 'DIR2'"_BLCK.AIM"
$
$ 	define tbmsk1_file 'DIR1'"_TRAB_BLCK.AIM"	! This is different between XCT1 (_TRAB_BLCK) & XCT2 (_TRAB_MASK)
$ 	define tbmsk2_file 'DIR2'"_TRAB_BLCK.AIM"
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
	-name 					msk1
	-filename					msk1_file
	

/read
	-name 					msk2
	-filename 					msk2_file
	

/read
	-name 					tbmsk1
	-filename					tbmsk1_file
	

/read
	-name 					tbmsk2
	-filename 					tbmsk2_file
	
	
/subtract_aims
	-input1					msk1
	-input2					tbmsk1
	-output					ctmsk1
	
/subtract_aims
	-input1					msk2
	-input2					tbmsk2
	-output					ctmsk2
	
	
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
  



! Transform the mask of the baseline image to the middle domain to create a gobj for the future registrations
/transform
  -in                        msk1
  -out                       mskt1
  -Tmat_file_name            "TRANSFORMATIONBL
  -img_interpol_option       0
  -el_size_mm_out            -1.000 -1.000 -1.000


! Dilate the mask of the 1st image to be used for both images, hence no effect of the periosteal mask
/offset_add
	-input					mskt1
	-add_offset					6 6 6

/dilation
  -input                     mskt1
  -output                    dil_mskt1
  -dilate_distance           5
  -continuous_at_boundary    0 0 0
  -use_previous_margin       false
  -metric                    101

/rename
  -old_name                  dil_mskt1
  -new_name                  mskt1  
  
  
/set_value
	-input					mskt1
	-value_object				127
	-value_background				0
	

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
  



! Transform the mask of the followup image to the middle domain to create a gobj for the future registrations
/transform
  -in                        msk2
  -out                       mskt2
  -Tmat_file_name            "TRANSFORMATIONFU
  -img_interpol_option       0
  -el_size_mm_out            -1.000 -1.000 -1.000


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

  
  
! Segmented images

/fft_laplace_hamming
  -input                     int1
  -output                    lh
  -redim_pow2                0 0 0
  -laplace_eps               0.500000
  -lp_cut_off_freq           0.400000
  -hamming_amp               1.000000

/delete
  -name                      int1
  
/norm_max
  -input                     lh
  -output                    norm
  -max                       200000.0
  -type_out                  short	

/delete
  -name                      lh
  

/threshold
  -input                     norm
  -output                    th1
  -lower_in_perm_aut_al      400.000000
  -upper_in_perm_aut_al      1000.000000
  -value_in_range            127
  -unit                      0
  
/delete
  -name                      norm

/cl_nr_extract
  -input                     th1
  -output                    cl1
  -min_number                25
  -max_number                0
  -value_in_range            127

/rename cl1 th1
  
/multiply_volumes
  -input1                    mskt1
  -input2                    comm
  -output                    mskb1
  -common_region_only        false
  -multiply_zero_pixels      true


  
/set_value
  -input                     mskb1
  -value_object              1
  -value_background          0
  
/multiply_volumes
  -input1                    th1
  -input2                    mskb1
  -output                    comb1
  -common_region_only        false
  -multiply_zero_pixels      true

/delete
	-name					th1
	
/delete
	-name					mskb1
  
! SECOND SCAN ==============================================================================  

/fft_laplace_hamming
  -input                     int2
  -output                    2h
  -redim_pow2                0 0 0
  -laplace_eps               0.500000
  -lp_cut_off_freq           0.400000
  -hamming_amp               1.000000

/delete
	-name					int2
	
/norm_max
  -input                     2h
  -output                    norm
  -max                       200000.0
  -type_out                  short	

/delete
  -name                      2h
  

/threshold
  -input                     norm
  -output                    th2
  -lower_in_perm_aut_al      400.000000
  -upper_in_perm_aut_al      1000.000000
  -value_in_range            127
  -unit                      0
  
/delete
  -name                      norm

/cl_nr_extract
  -input                     th2
  -output                    cl2
  -min_number                25
  -max_number                0
  -value_in_range            127

/rename cl2 th2

/multiply_volumes
  -input1                    mskt1
  -input2                    comm
  -output                    mskb2
  -common_region_only        false
  -multiply_zero_pixels      true

	
/set_value
  -input                     mskb2
  -value_object              1
  -value_background          0
  
/multiply_volumes
  -input1                    th2
  -input2                    mskb2
  -output                    comb2
  -common_region_only        false
  -multiply_zero_pixels      true

/delete
	-name					th2
	
/delete
	-name					mskb2

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

  
! Flatten the top and bottom surfaces    
/sub_get
  -input                     comb1u
  -output                    sub_input1
  -pos                       -1 -1 "top
  -dim                       -1 -1 "dim
  -global_pos_flag           false
	

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

  
/sub_get
  -input                     comb2u
  -output                    sub_input2
  -pos                       -1 -1 "top
  -dim                       -1 -1 "dim
  -global_pos_flag           false




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


..
$ EXIT
