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
$! Two types of overlaps can be calculated for 3D and MA registrations. First, the percentage of overlap for bone masks before and after registration
$! to quantify the degree of repositioning error. Second, the overlap between bone voxels after alignment during 3D or MA registration, which provides an estimate of registration accuracy.
$! Please note that in case of motion, the latter overlap tend to underestimate registration accuracy. This is because motion can distort bone shape
$! which cannot be aligned using rigid registration.
$!========================================================================================
$! This script:
$!				1. Calculates the overlap % of the 3D- and MA-registered and unregistered masks
$!              2. Calculates the overlap % of the flattened 3D- and MA-registered and unregistered masks
$!              3. Calculates the overlap % of the flattened 3D- and MA-registered common voxels relative to all voxels in segmented image 1
$
$! INPUTS:
$!              a. The grayscale ".AIM" file of the 1st image
$!              b. The grayscale ".AIM" file of the 2nd image
$!              c. The top slice number (indexed from 1) of the flattened 3D-registered region
$!              d. The top slice number (indexed from 1) of the flattened 3D-registered region
$!              e. The top slice number (indexed from 0) of the flattened MA-registered region
$!              f. The top slice number (indexed from 0) of the flattened MA-registered region
$
$! OUTPUTS:
$!				A. Calculated overlap%s written to a log file
$
$! call as:
$! @3D_MA_REG_OVERLAPS.com    c0000001.aim    c0000002.aim		TOP_3D 	   BOTTOM_3D      TOP_MA 	   BOTTOM_MA
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
$
$ COMMGOBJ1           :== 'DEV''DIR_READ''FILENAME'_3D.GOBJ
$ MSKT1   		:== 'DEV''DIR_READ''FILENAME'_MA_MASKT1.AIM
$ 	define gobj1_file 'DIR1'".GOBJ"
$
$ FE_COMMON_1   		:== 'DIR1'_FE_3D_DP_FLAT.AIM
$ FE_COMMON_2   		:== 'DIR2'_FE_3D_DP_FLAT.AIM
$ MA_COMMON_1   		:== 'DIR1'_FE_MA_DP_FLAT.AIM
$ MA_COMMON_2   		:== 'DIR2'_FE_MA_DP_FLAT.AIM
$
$ top1=F$INTEGER(P3)
$ bottom1=F$INTEGER(P4)
$ top=top1-1
$ bottom=bottom1-1
$ dim= bottom - top
$
$ top1=F$INTEGER(P5)
$ bottom1=F$INTEGER(P6)
$ top_ma=top1
$ bottom_ma=bottom1
$ dim_ma= bottom - top
$
$ IPL_BATCH

! PART 1: Overlap of the registered and unregistered masks showing the amount of misalignment
! 3D-registration
! calculate the original volume
/gobj_to_aim
  -gobj_filename             gobj1_file
  -output                    unreg_mask
  -peel_iter                 0
  
/vox_scanco_param
  -input                     unreg_mask
  -region_number             0
  
  
! calculate the 3d_registered volume
/gobj_to_aim
  -gobj_filename             "COMMGOBJ1
  -output                    reg_mask
  -peel_iter                 0
  
/vox_scanco_param
  -input                     reg_mask
  -region_number             0
  
! calculate the flattened 3d_registered volume
/bounding_box_cut
  -input                     reg_mask
  -output                    reg_maskt
  -z_only                    false
  -border                    0 0 0
  
/sub_get
  -input                     reg_maskt
  -output                    reg_mask_flat
  -pos                       -1 -1 "top
  -dim                       -1 -1 "dim
  -global_pos_flag           false

	
/vox_scanco_param
  -input                     reg_mask_flat
  -region_number             0


! PART 2: The ratio of the matching voxels in the overlapped registered scans to total voxels of the 1st image, showing 3D-registeration quality
/read
  -name                      reg_seg1
  -filename                  "FE_COMMON_1
  
/set_value
  -input                     reg_seg1
  -value_object              1
  -value_background          0
  
/read
  -name                      reg_seg2
  -filename                  "FE_COMMON_2

/set_value
  -input                     reg_seg2
  -value_object              2
  -value_background          0


! Remove the non-matching voxels from the overlapped images
/concat
	-input1 					reg_seg1
	-input2 					reg_seg2
	-output 					overlay
	-common_region_only 			false
	-add_not_overlay 				true
	-make_edge 				false
	-shift_ofin2 				0 0 0
	-turnangle 					0.0
	-turnpoint_global 				-1 -1
	
/threshold
	-input					overlay
	-output					overlay_filtered
	-lower_in_perm_aut_al				3
	-upper_in_perm_aut_al				3
	-value_in_range				127
	-unit					5
	

! calculate the voxels in the 3D-registered segmented 1st image
/vox_scanco_param
  -input                     reg_seg1
  -region_number             0
  
  
! calculate the mathcing voxels of the two images
/vox_scanco_param
  -input                     overlay_filtered
  -region_number             0
  
  
  

! MA-registration
! calculate the MA_registered volume

/read
  -name                      mskt1
  -filename                  "MSKT1
  
/bounding_box_cut
  -input                     mskt1
  -output                    reg_mask
  -z_only                    false
  -border                    0 0 0
  
/vox_scanco_param
  -input                     reg_mask
  -region_number             0
  
! calculate the flattened 3d_registered volume
/bounding_box_cut
  -input                     reg_mask
  -output                    reg_maskt
  -z_only                    false
  -border                    0 0 0
  
/sub_get
  -input                     reg_maskt
  -output                    reg_mask_flat
  -pos                       -1 -1 "top_ma
  -dim                       -1 -1 "dim_ma
  -global_pos_flag           false

	
/vox_scanco_param
  -input                     reg_mask_flat
  -region_number             0


! PART 2: The ratio of the matching voxels in the overlapped registered scans to total voxels of the 1st image, showing 3D-registeration quality
/read
  -name                      reg_seg1
  -filename                  "MA_COMMON_1
  
/set_value
  -input                     reg_seg1
  -value_object              1
  -value_background          0
  
/read
  -name                      reg_seg2
  -filename                  "MA_COMMON_2

/set_value
  -input                     reg_seg2
  -value_object              2
  -value_background          0


! Remove the non-matching voxels from the overlapped images
/concat
	-input1 					reg_seg1
	-input2 					reg_seg2
	-output 					overlay
	-common_region_only 			false
	-add_not_overlay 				true
	-make_edge 				false
	-shift_ofin2 				0 0 0
	-turnangle 					0.0
	-turnpoint_global 				-1 -1
	
/threshold
	-input					overlay
	-output					overlay_filtered
	-lower_in_perm_aut_al				3
	-upper_in_perm_aut_al				3
	-value_in_range				127
	-unit					5
	

! calculate the voxels in the 3D-registered segmented 1st image
/vox_scanco_param
  -input                     reg_seg1
  -region_number             0
  
  
! calculate the mathcing voxels of the two images
/vox_scanco_param
  -input                     overlay_filtered
  -region_number             0
   
..
$
$ DIRECT :== 'DEV''DIR_READ'
$ SET DEFAULT 'DIRECT'
$ CREATE CHECK.TXT
$ EXIT