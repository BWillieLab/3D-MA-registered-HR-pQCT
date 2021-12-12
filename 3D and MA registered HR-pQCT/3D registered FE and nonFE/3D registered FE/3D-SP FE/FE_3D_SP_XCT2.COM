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
$!				1. Aligns the grayscale images of two scans from XCT2 device
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
$! @FE_3D_SP_XCT2.com    c0000001.aim    c0000002.aim		TOP		BOTTOM
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
$  JOBFILE := 'DIR'_SEG_3Dreg_CONTOUR_FE
$ FE_FLAT_1   		:== 'DIR1'_FE_3D_SP_FLAT.AIM
$ FE_FLAT_2   		:== 'DIR2'_FE_3D_SP_FLAT.AIM
$
$ 	define ctmask1_file 'DIR1'"_CORT_MASK.AIM" 	! This is different between XCT1 (_BLCK) & XCT2 (_CORT_MASK)
$ 	define ctmask2_file 'DIR2'"_CORT_MASK.AIM"
$
$ 	define tbmask1_file 'DIR1'"_TRAB_MASK.AIM"	! This is different between XCT1 (_TRAB_BLCK) & XCT2 (_TRAB_MASK)
$ 	define tbmask2_file 'DIR2'"_TRAB_MASK.AIM"
$
$
$ top=top1-1
$ bottom=bottom1-1
$ dim= bottom - top 
$
$
$ IPL_BATCH

! STEP 1:
! Creating a component-labeled mask, cropped to 104 slices to match standard UPAT_LIST default output
/read
  -name                      aimt1
  -filename                  "P1

! ------ read common aim, clean it and create a cortex mask file ---------
/read
  -name                      comm
  -filename                  comm_file
  
/read
  -name                      ctmask1
  -filename                  ctmask1_file

/set_value
  -input                     ctmask1
  -value_object              1
  -value_background          0

  
/read
  -name                      tbmask1
  -filename                  tbmask1_file
  
/set_value
  -input                     tbmask1
  -value_object              1
  -value_background          0


! Crop the mask using the common mask  
/multiply_volumes
  -input1                    ctmask1
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
  -new_name                  ctmask1 
  
/multiply_volumes
  -input1                    tbmask1
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
  -new_name                  tbmask1 

! Create and dilate the periosteal mask  
/concat
	-input1 					ctmask1
	-input2 					tbmask1
	-output 					msk1
	-common_region_only 			false
	-add_not_overlay 				true
	-make_edge 				false
	-shift_ofin2 				0 0 0
	-turnangle 					0.0
	-turnpoint_global 				-1 -1 

/offset_add
	-input					msk1
	-add_offset					6 6 6


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
  
/subtract_aims
  -input1                    msk1
  -input2                    tbmask1
  -output                    ctmask1
  
/set_value
  -input                     ctmask1
  -value_object              1
  -value_background          0
  	
! ------ combine images: 126=trab 127=cort bone ----------
/offset_add
  -input                     aimt1
  -add_offset                1 1 1


/fill_offset_duplicate
  -input                     aimt1


/sup_divide
  -input                     aimt1
  -supdim_numbers            4 4 1
  -testoff_pixels            1 1 1
  -suppos_pixels_local       -1 -1 -1
  -subdim_pixels             -1 -1 -1
  
/copy
  -in                        aimt1
  -out                       aimt12

/seg_gauss
  -input                     aimt1
  -output                    aimt1_gauss
  -sigma                     0.800000
  -support                   1
  -lower_in_perm_aut_al      320.000000
  -upper_in_perm_aut_al      3000.000000
  -value_in_range            126
  -unit                      2  

/delete
	-name					aimt1

/seg_gauss
  -input                     aimt12
  -output                    aimt12_gauss
  -sigma                     0.800000
  -support                   1
  -lower_in_perm_aut_al      450.000000
  -upper_in_perm_aut_al      3000.000000
  -value_in_range            127
  -unit                      2  

/delete
	-name					aimt12
	
/multiply_volumes
  -input1                    aimt1_gauss
  -input2                    tbmask1
  -output                    tb_input1
  -common_region_only        false
  -multiply_zero_pixels      true

/delete
	-name					aimt1_gauss
  
/multiply_volumes
  -input1                    aimt12_gauss
  -input2                    ctmask1
  -output                    ct_input1
  -common_region_only        false
  -multiply_zero_pixels      true

/delete
	-name					aimt12_gauss

/add_aims
  -input1                    tb_input1
  -input2                    ct_input1
  -output                    comb1

/delete
	-name					tb_input1
	
/delete
	-name					ct_input1
  
! SECOND SCAN ==============================================================================  
/read
  -name                      aim2
  -filename                  "P2

! ------ read common aim, clean it and create a cortex mask file ---------
/read
  -name                      comm
  -filename                  comm_file

  
/read
  -name                      tbmask2
  -filename                  tbmask2_file

/transform
	-in 					tbmask2
	-out 					tbmskt2
	-Tmat_file_name 				tmat_file
	-img_interpol_option 				0
	-el_size_mm_out 				-1 -1 -1
	-Tmat_invert				false
	
/delete
	-name					tbmask2
	
/set_value
  -input                     tbmskt2
  -value_object              1
  -value_background          0
 

/multiply_volumes
  -input1                    tbmskt2
  -input2                    comm
  -output                    tbmskb2
  -common_region_only        false
  -multiply_zero_pixels      true

/set_value
  -input                     tbmskb2
  -value_object              1
  -value_background          0


/rename
  -old_name                  tbmskb2
  -new_name                  tbmskt2 
  
/subtract_aims
  -input1                    msk1        ! Use the dilated mask of the 1st image also for the 2nd image
  -input2                    tbmskt2
  -output                    ctmskt2


/set_value
  -input                     ctmskt2
  -value_object              1
  -value_background          0

  

/delete
	-name					comm
	
! ------ combine images: 126=trab 127=cort bone ----------
/transform
	-in 					aim2
	-out 					aimt2
	-Tmat_file_name 				tmat_file
	-img_interpol_option 				2
	-el_size_mm_out 				-1 -1 -1
	-Tmat_invert				false
	
/delete
	-name					aim2

/offset_add
  -input                     aimt2
  -add_offset                1 1 1


/fill_offset_duplicate
  -input                     aimt2


/sup_divide
  -input                     aimt2
  -supdim_numbers            4 4 1
  -testoff_pixels            1 1 1
  -suppos_pixels_local       -1 -1 -1
  -subdim_pixels             -1 -1 -1
  
/copy
  -in                        aimt2
  -out                       aimt22

/seg_gauss
  -input                     aimt2
  -output                    aimt2_gauss
  -sigma                     0.800000
  -support                   1
  -lower_in_perm_aut_al      320.000000
  -upper_in_perm_aut_al      3000.000000
  -value_in_range            126
  -unit                      2  


/delete
	-name					aimt2

/seg_gauss
  -input                     aimt22
  -output                    aimt22_gauss
  -sigma                     0.800000
  -support                   1
  -lower_in_perm_aut_al      450.000000
  -upper_in_perm_aut_al      3000.000000
  -value_in_range            127
  -unit                      2  


/delete
	-name					aimt22

	
/multiply_volumes
  -input1                    aimt2_gauss
  -input2                    tbmskt2
  -output                    tb_input2
  -common_region_only        false
  -multiply_zero_pixels      true

/delete
	-name					aimt2_gauss
	
/multiply_volumes
  -input1                    aimt22_gauss
  -input2                    ctmskt2
  -output                    ct_input2
  -common_region_only        false
  -multiply_zero_pixels      true

/delete
	-name					aimt22_gauss
	
/delete
	-name					ctmskt2

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

/copy
  -in                        comb1trim
  -out                       comb1c
  
/set_value
  -input                     comb1c
  -value_object              1
  -value_background          1

/bounding_box_cut
  -input                     comb2
  -output                    comb2trim
  -z_only                    false
  -border                    0 0 0
  
/copy
  -in                        comb2trim
  -out                       comb2c

/set_value
  -input                     comb2c
  -value_object              1
  -value_background          1
  
! Now create the mask of the common region
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
/cl_rank_extract
	-input			sub_input_char1
	-output					cl
	-first_rank				1
	-last_rank				1
	-connect_boundary		false
	-value_in_range			1

/multiply_volumes
  -input1                    sub_input_char1
  -input2                    cl
  -output                    fe_input1
  -common_region_only        false
  -multiply_zero_pixels      true

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
/cl_rank_extract
	-input			sub_input_char2
	-output					cl
	-first_rank				1
	-last_rank				1
	-connect_boundary		false
	-value_in_range			1

/multiply_volumes
  -input1                    sub_input_char2
  -input2                    cl
  -output                    fe_input2
  -common_region_only        false
  -multiply_zero_pixels      true

/delete
	-name					sub_input_char2
	
/write_v020
  -name                      fe_input2
  -filename                  "FE_FLAT_2 
  -compress_type             none          ! Non-binarized since may be used out of IPL
  -version_020               false
..
$ EXIT
