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
$! 3D registration process to align two images and to find the common volume for density and microstructural outcomes for XCT2 scans.
$! The first step was to find the transformation matrix that would align the moving image onto the reference image.
$! In the second step, solid blocks were created and aligned to find the common volume. Finally, during the third step,
$! the common volume was used to modify the original contour lines (i.e., Scanco GOBJ files).
$!========================================================================================
$! This script:
$!				1. Finds the transformation matrix that aligns the two images
$!              2. Identifies the common region in the domain of the 1st image and transforms it to the 2nd image domain
$!              3. Modifies the periosteal contours to only include the common region
$!              4. Create the common cortical and trabecular masks in the 1st image domain
$
$! INPUTS:
$!              a. The grayscale ".AIM" file of the 1st image
$!              b. The grayscale ".AIM" file of the 2nd image
$
$! OUTPUTS:
$!              A. The transformation matrix aligning the 2nd image onto the 1st image (.DAT file)
$!              B. The common block in the domain of the 1st image
$!              C. The common whole bone masks in the domain of each image
$!              D. The common trabecular and cortical bone masks in the domain of the 1st image
$!              E. The grayscale images in the 1st image domain
$!              F. The 3D-registered periosteal contours
$
$! Version history:
$
$! call as:
$! @3DREGISTER_2IMAGES_XCT2.com    c0000001.aim    c0000002.aim
$!==========================================================================
$ 	define in1_file 'p1'
$ 	define in2_file 'p2'
$
$!============================================================================
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
$ TRANSFORMATION    	:== 'DEV2''DIR_READ2''FILENAME'_'FILENAME2'_TMAT.DAT
$ COMM1   		:== 'DEV''DIR_READ''FILENAME'_2COMMON1.AIM
$ COMMCT1   		:== 'DEV''DIR_READ''FILENAME'_2COMMONCORT1.AIM
$ COMMTB1   		:== 'DEV''DIR_READ''FILENAME'_2COMMONTRAB1.AIM
$ COMM2 		:== 'DEV2''DIR_READ2''FILENAME2'_2COMMON2.AIM
$ COMMCT2 		:== 'DEV2''DIR_READ2''FILENAME2'_2COMMONCORT2.AIM
$ COMMTB2 		:== 'DEV2''DIR_READ2''FILENAME2'_2COMMONTRAB2.AIM
$ GST1 		:== 'DEV''DIR_READ''FILENAME'_TL_GRAY.AIM
$ GST2 		:== 'DEV2''DIR_READ2''FILENAME2'_TL_GRAY.AIM
$
$ COMMGOBJ1           :== 'DEV''DIR_READ''FILENAME'_3D.GOBJ
$ COMMGOBJ2           :== 'DEV2''DIR_READ2''FILENAME2'_3D.GOBJ
$
$ COMM1BLK   		:== 'DEV''DIR_READ''FILENAME'_2COMMON1_BLK.AIM
$
$ 	define gobj1_file 'DIR1'".GOBJ"
$ 	define gobj2_file 'DIR2'".GOBJ"
$
!$ 	define ins1_file 'DIR1'"_SEG.AIM"
!$ 	define ins2_file 'DIR2'"_SEG.AIM"
!$ 	define ins3_file 'DIR3'"_SEG.AIM"
$
$ 	define ctmsk1_file 'DIR1'"_CORT_MASK.AIM" 	! This is different between XCT1 (_BLCK) & XCT2 (_CORT_MASK)
$ 	define ctmsk2_file 'DIR2'"_CORT_MASK.AIM"
$
$ 	define tbmsk1_file 'DIR1'"_TRAB_MASK.AIM"	! This is different between XCT1 (_TRAB_BLCK) & XCT2 (_TRAB_MASK)
$ 	define tbmsk2_file 'DIR2'"_TRAB_MASK.AIM"
$
$!============================================================================
$!
$ 	ipl_scanco_prog := $um:ipl_scanco_m.exe
$!
$ ON ERROR THEN EXIT
$!
$ 	ipl_scanco_prog

! Read the grayscale 'AIM' file for the 1st scan
/read
	-name 					in1
	-filename					in1_file

! Read the grayscale 'AIM' file for the 2nd scan
/read
	-name 					in2
	-filename 					in2_file
	
! Read the cortical masks
/read
	-name 					ctmsk1
	-filename					ctmsk1_file
	
/delete
	-name					ctmsk1_file


/read
	-name 					ctmsk2
	-filename 					ctmsk2_file
	
/delete
	-name					ctmsk2_file


! Read the trabecular masks
/read
	-name 					tbmsk1
	-filename					tbmsk1_file
	
/delete
	-name					tbmsk1_file


/read
	-name 					tbmsk2
	-filename 					tbmsk2_file
	
/delete
	-name					tbmsk2_file
	

	
! Find the transformation matrix that aligns the 2nd image onto the 1st image	
/register
	-in1 					in1
	-gobj_filename_in1 		gobj1_file
	-in2 					in2
	-gobj_filename_in2 		gobj2_file
	-Tmat_file_name 		"TRANSFORMATION
	-orientation_search 	1
	-initial_rotation 		0.0 0.0 0.0
	-initial_translation 	0.0 0.0 0.0
	-delta_rotation 		0.30 0.30 0.50
	-delta_translation 		0.30 0.30 0.50
	-resolution_scaling 	10 4 1
    -delta_scaling          1.00 0.10 0.10
    -tolerance              0.00001
    -min_corr_coef          0.5
    -min_method             1
    -object_func            1
    -max_nr_iter            1000
	

! Transform the grayscale of tge 2nd image to the domain of the 1st image	
/transform
	-in 					in2
	-out 					int2
	-Tmat_file_name 				"TRANSFORMATION
	-img_interpol_option 				2
	-el_size_mm_out 				-1 -1 -1
	-Tmat_invert				false
	
/write_v020
	-name                      				int2
	-filename                  				"GST2
	-compress_type             			none
	-version_020               				false
	
	
/delete
	-name					int2


/write_v020
	-name                      				in1
	-filename                  				"GST1
	-compress_type             			none
	-version_020               				false
	

! Now create a solid block from each aim file to find the common region without losing the periosteal surface voxels
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
	-value_background		5


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
	-in 					blk2
	-out 					blkt2
	-Tmat_file_name 		"TRANSFORMATION
	-img_interpol_option 	0
	-el_size_mm_out 		-1 -1 -1
	-Tmat_invert			false
	

! Add the solid blocks of the two images, and threshold using a value of 15 (5 + 10 = regions common between images)
/concat
	-input1 					blk1
	-input2 					blkt2
	-output 					added_blocks
	-common_region_only 			false
	-add_not_overlay 				true
	-make_edge 				false
	-shift_ofin2 				0 0 0
	-turnangle 					0.0
	-turnpoint_global 				-1 -1

		
/threshold
	-input					added_blocks
	-output					common_block
	-lower_in_perm_aut_al				15
	-upper_in_perm_aut_al				15
	-value_in_range				127
	-unit					5

/delete
	-name					added_blocks
	
	
/set_value
	-input					common_block
	-value_object				127
	-value_background				0



! Transform the common_block to the domain of the first image
/transform
	-in 					common_block
	-out 					common_blockt
	-Tmat_file_name 		"TRANSFORMATION
	-img_interpol_option 	0
	-el_size_mm_out 		-1 -1 -1
	-Tmat_invert			true

!Smooth the transformed mask to remove the effect of rotation (e.g., holes)(dilation of the mask at this level is ok)
/median_filter           
    -input                  common_blockt
	-output                 common_blockt_med
	-support                1	

/rename
  -old_name                  common_blockt_med
  -new_name                  common_blockt

/copy
  -in                        common_block
  -out                       common_block_intact
  
! Crop the common_block by the gobjs to create the common gobj	
/gobj_maskaimpeel_ow		
	-input_output			common_block
	-gobj_filename			gobj1_file
	-peel_iter				0

	
/bounding_box_cut
	-input					common_block
	-output					common_blocktrim
	-z_only					false
	-border					0 0 0

/togobj_from_aim
	-input						common_blocktrim
	-gobj_filename					"COMMGOBJ1
	-min_elements					0
	-max_elements					0
	-curvature_smooth				1
	
/write_v020
	-name                   common_blocktrim
	-filename               "COMM1
	-compress_type          none
	-version_020            false



/gobj_maskaimpeel_ow		
	-input_output			common_blockt
	-gobj_filename			gobj2_file
	-peel_iter				0
	
	
/bounding_box_cut
	-input					common_blockt
	-output					common_blockttrim
	-z_only					false
	-border					0 0 0
	
/delete
	-name					common_blockt


/togobj_from_aim
	-input						common_blockttrim
	-gobj_filename					"COMMGOBJ2
	-min_elements					0
	-max_elements					0
	-curvature_smooth				1


/write_v020
	-name                      				common_blockttrim
	-filename                  				"COMM2
	-compress_type             			none
	-version_020               				false



/set_value
	-input					common_block_intact
	-value_object				127
	-value_background				0

/bounding_box_cut
	-input					common_block_intact
	-output					common_blockttrim
	-z_only					false
	-border					0 0 0

/write_v020
	-name                      				common_blockttrim
	-filename                  				"COMM1BLK
	-compress_type             			none
	-version_020               				false
	

! ================ Finding common region of cortical and trabecular masks =====================
!!! CORTICAL BONE
/transform
	-in 					ctmsk2
	-out 					ctmskt2
	-Tmat_file_name 				"TRANSFORMATION
	-img_interpol_option 				0
	-el_size_mm_out 				-1 -1 -1
	-Tmat_invert				false
	
/delete
	-name					ctmsk2

	
! Keeping the common region from each scan
/set_value
	-input					common_block
	-value_object				1
	-value_background				0

/set_value
	-input					ctmsk1
	-value_object				1
	-value_background				0

/set_value
	-input					ctmskt2
	-value_object				1
	-value_background				0 
	

/concat
	-input1 					common_block
	-input2 					ctmsk1
	-output 					ct1common
	-common_region_only 			false
	-add_not_overlay 				true
	-make_edge 				false
	-shift_ofin2 				0 0 0
	-turnangle 					0.0
	-turnpoint_global 				-1 -1

/threshold
	-input					ct1common
	-output					ct1commonm
	-lower_in_perm_aut_al				2
	-upper_in_perm_aut_al				2
	-value_in_range				127
	-unit					5
	
/delete
	-name					ct1common
	
/set_value
	-input					ct1commonm
	-value_object				127
	-value_background				0
	
	
	
/concat
	-input1 					common_block
	-input2 					ctmskt2
	-output 					ct2common
	-common_region_only 			false
	-add_not_overlay 				true
	-make_edge 				false
	-shift_ofin2 				0 0 0
	-turnangle 					0.0
	-turnpoint_global 				-1 -1

/threshold
	-input					ct2common
	-output					ct2commonm
	-lower_in_perm_aut_al				2
	-upper_in_perm_aut_al				2
	-value_in_range				127
	-unit					5

/delete
	-name					ct2common
	
/set_value
	-input					ct2commonm
	-value_object				127
	-value_background				0
	
/bounding_box_cut
	-input					ct1commonm
	-output					ct1commonmtrim
	-z_only					false
	-border					0 0 0
	
/delete
	-name					ct1commonm
	
/write_v020
	-name                      				ct1commonmtrim
	-filename                  				"COMMCT1
	-compress_type             			none
	-version_020               				false


/bounding_box_cut
	-input					ct2commonm
	-output					ct2commonmtrim
	-z_only					false
	-border					0 0 0
	
/delete
	-name					ct2commonm
	
/write_v020
	-name                      				ct2commonmtrim
	-filename                  				"COMMCT2
	-compress_type             			none
	-version_020               				false
	
/delete
	-name					ct1commonmtrim
	
/delete
	-name					ct2commonmtrim
	

!!! TRABECULAR BONE
/transform
	-in 					tbmsk2
	-out 					tbmskt2
	-Tmat_file_name 				"TRANSFORMATION
	-img_interpol_option 				0
	-el_size_mm_out 				-1 -1 -1
	-Tmat_invert				false
	
/delete
	-name					tbmsk2
	


! Keeping the common region from each scan
/set_value
	-input					common_block
	-value_object				1
	-value_background				0

/set_value
	-input					tbmsk1
	-value_object				1
	-value_background				0

/set_value
	-input					tbmskt2
	-value_object				1
	-value_background				0 


/concat
	-input1 					common_block
	-input2 					tbmsk1
	-output 					tb1common
	-common_region_only 			false
	-add_not_overlay 				true
	-make_edge 				false
	-shift_ofin2 				0 0 0
	-turnangle 					0.0
	-turnpoint_global 				-1 -1

/threshold
	-input					tb1common
	-output					tb1commonm
	-lower_in_perm_aut_al				2
	-upper_in_perm_aut_al				2
	-value_in_range				127
	-unit					5
	
/delete
	-name					tb1common
	
/set_value
	-input					tb1commonm
	-value_object				127
	-value_background				0
	
	
	
/concat
	-input1 					common_block
	-input2 					tbmskt2
	-output 					tb2common
	-common_region_only 			false
	-add_not_overlay 				true
	-make_edge 				false
	-shift_ofin2 				0 0 0
	-turnangle 					0.0
	-turnpoint_global 				-1 -1


/delete
	-name					common_block
	
	
/threshold
	-input					tb2common
	-output					tb2commonm
	-lower_in_perm_aut_al				2
	-upper_in_perm_aut_al				2
	-value_in_range				127
	-unit					5
	
/delete
	-name					tb2common
	
/set_value
	-input					tb2commonm
	-value_object				127
	-value_background				0
	
/bounding_box_cut
	-input					tb1commonm
	-output					tb1commonmtrim
	-z_only					false
	-border					0 0 0
	
/delete
	-name					tb1commonm

	
/write_v020
	-name                      				tb1commonmtrim
	-filename                  				"COMMTB1
	-compress_type             			none
	-version_020               				false


/bounding_box_cut
	-input					tb2commonm
	-output					tb2commonmtrim
	-z_only					false
	-border					0 0 0
	
/delete
	-name					tb2commonm
	
/write_v020
	-name                      				tb2commonmtrim
	-filename                  				"COMMTB2
	-compress_type             			none
	-version_020               				false
	
	
/delete
	-name					tb1commonmtrim
	
/delete
	-name					tb2commonmtrim
	
..
$ exit
