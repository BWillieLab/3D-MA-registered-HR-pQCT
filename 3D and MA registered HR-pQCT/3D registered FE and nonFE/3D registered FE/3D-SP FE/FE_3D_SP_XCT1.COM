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
$! The current script performs an alternative version of 3D-registered microFE, where instead of masking each image using its
$! own mask, both images are masked using a dilated mask, effectively eliminating the influence of contouring error and interpolation error of the masks. This approach is called 3D-SP (same periosteal).
$!
$!========================================================================================
$! This script:
$!				1. Aligns the grayscale images of two scans from XCT device
$!              2. Segmentes the images in the common domain
$!              3. Masks both images using the dilated mask of the 1st image to remove the effect of masking
$!              4. Flattens the segmented images at the top and bottom surfaces for applying FE boundary conditions
$!              5. Maps both segmented images in the common coordinates
$
$! INPUTS:
$!              a. The grayscale ".AIM" file of the 1st image
$!              b. The grayscale ".AIM" file of the 2nd image
$!              c. The top slice number (indexed from 1) to be flattened (the slice number at which the flat region starts)
$!              d. The top slice number (indexed from 1) to be flattened (the slice number at which the flat region ends)
$
$! OUTPUTS:
$!              A. The segmented image of the 1st image cropped by the common region
$!              B. The segmented image of the 2nd image cropped by the common region
$!              C. The segmented image of the 1st image cropped by the common region and flattened
$!              D. The segmented image of the 2nd image cropped by the common region and flattened
$
$! Version history:
$
$! call as:
$! @FE_3D_SP_XCT1.com    c0000001.aim    c0000002.aim		TOP		BOTTOM
$!==========================================================================
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
$ top1=F$INTEGER(P3)
$ bottom1=F$INTEGER(P4)
$
$ define comm_file 'DIR1'"_2COMMON1_BLK.AIM"
$ define tmat_file 'DEV2''DIR_READ2''FILENAME'_'FILENAME2'_TMAT.DAT
$
$ 	define msk1_file 'DIR1'"_BLCK.AIM" 	! This is different between XCT1 (_BLCK) & XCT2 (_CORT_MASK)
$ 	define msk2_file 'DIR2'"_BLCK.AIM"
$
$  JOBFILE := 'DIR'_SEG_3Dreg_CONTOUR_FE
$ FE_FLAT_1   		:== 'DIR1'_FE_3D_SP_FLAT.AIM
$ FE_FLAT_2   		:== 'DIR2'_FE_3D_SP_FLAT.AIM
$ FE_COMMON_1   		:== 'DIR1'_FE_3D_SP.AIM
$ FE_COMMON_2   		:== 'DIR2'_FE_3D_SP.AIM
$
$
$ top=top1-1
$ bottom=bottom1-1
$ dim= bottom - top 
$
$ IPL_BATCH

! STEP 1:
! Read the grayscale and masks of the 1st image
/read
  -name                      aim1
  -filename                  "P1

/read
	-name 					msk1
	-filename				msk1_file

! Remove slices at the top and bottom to compensate for dilation increase in size
/offset_add
	-input					msk1
	-add_offset					6 6 6

! Dilate the mask of the 1st image to be used for both images
/dilation
  -input                     msk1
  -output                    dil_msk1
  -dilate_distance           5
  -continuous_at_boundary    0 0 0
  -use_previous_margin       false
  -metric                    101

  
/rename
  -old_name                  dil_msk1
  -new_name                  msk1  
  
/set_value
  -input                     msk1
  -value_object              1
  -value_background          0
  
! ------ read common aim, clean it and create a cortex mask file ---------
/read
  -name                      comm
  -filename                  comm_file

/set_value
  -input                     comm
  -value_object              1
  -value_background          0
  
  
! ------ combine images: 126=trab 127=cort bone ----------

! Segmentation starts with Laplace-Hamming filter
/fft_laplace_hamming
  -input                     aim1
  -output                    lh
  -redim_pow2                0 0 0
  -laplace_eps               0.500000
  -lp_cut_off_freq           0.400000
  -hamming_amp               1.000000


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


! Crop the dilated mask using the common mask  
/multiply_volumes
  -input1                    msk1
  -input2                    comm
  -output                    mskb1
  -common_region_only        false
  -multiply_zero_pixels      true

/set_value
  -input                     mskb1
  -value_object              1
  -value_background          0


! crop the segmented image using the common mask  
/multiply_volumes
  -input1                    th1
  -input2                    mskb1
  -output                    comb1
  -common_region_only        false
  -multiply_zero_pixels      true

/delete
	-name					th1
  
! SECOND SCAN ==============================================================================  
/read
  -name                      aim2
  -filename                  "P2

! ------ read common aim, clean it and create a cortex mask file ---------
/read
  -name                      comm
  -filename                  comm_file

/set_value
  -input                     comm
  -value_object              1
  -value_background          0

! Transform the components of the 2nd image to the domain of the 1st image  
/transform
	-in 					aim2
	-out 					aimt2
	-Tmat_file_name 				tmat_file
	-img_interpol_option 				2		          ! Cubic interpolation for the grayscale images
	-el_size_mm_out 				-1 -1 -1
	-Tmat_invert				false
  
/delete
	-name					aim2

/fft_laplace_hamming
  -input                     aimt2
  -output                    2h
  -redim_pow2                0 0 0
  -laplace_eps               0.500000
  -lp_cut_off_freq           0.400000
  -hamming_amp               1.000000


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


! For the 2nd image also, we use the dilated mask of the 1st image
/multiply_volumes
  -input1                    msk1
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

! ====== Now make sure both scans have the same region selected ========

! Trim the images
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

/delete
	-name					combmsk

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
  -filename                  "FE_COMMON_1
  -compress_type             none              ! Non-binarized since may be used out of IPL
  -version_020               false 
  
! Flatten the top and bottom surfaces  
/sub_get
  -input                     comb1u
  -output                    sub_input1
  -pos                       -1 -1 "top
  -dim                       -1 -1 "dim
  -global_pos_flag           false
	
/delete
	-name					comb1u

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
  -filename                  "FE_COMMON_2
  -compress_type             none              ! Non-binarized since may be used out of IPL
  -version_020               false 
  
/sub_get
  -input                     comb2u
  -output                    sub_input2
  -pos                       -1 -1 "top
  -dim                       -1 -1 "dim
  -global_pos_flag           false
	
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
..
$ EXIT
