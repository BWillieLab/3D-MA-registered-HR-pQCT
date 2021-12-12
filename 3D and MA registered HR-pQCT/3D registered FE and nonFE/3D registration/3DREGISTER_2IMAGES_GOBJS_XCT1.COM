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
$!              1. Modifies the endosteal, and cortical contours to only include the common region
$
$! INPUTS:
$!              a. The grayscale ".AIM" file of the 1st image
$!              b. The grayscale ".AIM" file of the 2nd image
$
$! OUTPUTS:
$!              A. The 3D-registered AUTO periosteal contours
$!              B. The 3D-registered cortical contours
$!              C. The 3D-registered posority contours
$!              D. The 3D-registered endosteal contours
$
$! Version history:
$
$! call as:
$! @3DREGISTER_2IMAGES_XCT1.com    c0000001.aim    c0000002.aim
$!==========================================================================
$ 	define in1_file 'p1'
$ 	define in2_file 'p2'
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
$ COMM2 		:== 'DEV2''DIR_READ2''FILENAME2'_2COMMON2.AIM
$ BASEauto           :== 'DEV''DIR_READ''FILENAME'_AUTO2_3D.GOBJ					
$ FOLLOW06auto           :== 'DEV2''DIR_READ2''FILENAME2'_AUTO2_3D.GOBJ				
$ BASEcrtx           :== 'DEV''DIR_READ''FILENAME'_CRTX_3D.GOBJ					
$ FOLLOW06crtx           :== 'DEV2''DIR_READ2''FILENAME2'_CRTX_3D.GOBJ				
$ BASEcrtxpo          :== 'DEV''DIR_READ''FILENAME'_CRTXPO_3D.GOBJ					
$ FOLLOW06crtxpo           :== 'DEV2''DIR_READ2''FILENAME2'_CRTXPO_3D.GOBJ			
$ BASEtrab          :== 'DEV''DIR_READ''FILENAME'_TRAB_3D.GOBJ					
$ FOLLOW06trab           :== 'DEV2''DIR_READ2''FILENAME2'_TRAB_3D.GOBJ					
$
$ 	define auto_gobj1_file 'DIR1'"_AUTO2.GOBJ"
$ 	define auto_gobj2_file 'DIR2'"_AUTO2.GOBJ"
$
$ 	define crtx_gobj1_file 'DIR1'"_CRTX.GOBJ"
$ 	define crtx_gobj2_file 'DIR2'"_CRTX.GOBJ"
$
$ 	define crtxpo_gobj1_file 'DIR1'"_CRTXPO.GOBJ"
$ 	define crtxpo_gobj2_file 'DIR2'"_CRTXPO.GOBJ"
$
$ 	define trab_gobj1_file 'DIR1'"_TRAB.GOBJ"
$ 	define trab_gobj2_file 'DIR2'"_TRAB.GOBJ"
$
$!============================================================================
$!
$!
$ 	ipl_scanco_prog := $um:ipl_scanco_m.exe
$!
$ ON ERROR THEN EXIT
$!
$ 	ipl_scanco_prog
	

! Read the common region 'AIM' file
/read
	-name 					common
	-filename					"COMM1

/bounding_box_cut
	-input					common
	-output					commontrim
	-z_only					false
	-border					0 0 0

/delete
	-name					common
	
! Now cut the commmon region in each domain with gobj
/gobj_maskaimpeel_ow		
	-input_output			commontrim
	-gobj_filename			auto_gobj1_file
	-peel_iter				0
	
	
/bounding_box_cut
	-input					commontrim
	-output					commontrimt
	-z_only					false
	-border					0 0 0
	
/delete
	-name					commontrim
	
/togobj_from_aim
	-input						commontrimt
	-gobj_filename					"BASEauto
	-min_elements					0
	-max_elements					0
	-curvature_smooth				1
	
! =========================================================
! Now cut the commmon region in each domain with gobj
/gobj_maskaimpeel_ow		
	-input_output			commontrim
	-gobj_filename			crtx_gobj1_file
	-peel_iter				0

/bounding_box_cut
	-input					commontrim
	-output					commontrimt
	-z_only					false
	-border					0 0 0
	
/togobj_from_aim
	-input						commontrimt
	-gobj_filename					"BASEcrtx
	-min_elements					0
	-max_elements					0
	-curvature_smooth				1
	
! =========================================================
! Now cut the commmon region in each domain with gobj
/gobj_maskaimpeel_ow		
	-input_output			commontrim
	-gobj_filename			crtxpo_gobj1_file
	-peel_iter				0

/bounding_box_cut
	-input					commontrim
	-output					commontrimt
	-z_only					false
	-border					0 0 0	
	
/togobj_from_aim
	-input						commontrimt
	-gobj_filename					"BASEcrtxpo
	-min_elements					0
	-max_elements					0
	-curvature_smooth				1
	
/delete
	-name					commontrimt
	
! =========================================================		
! Now cut the commmon region in each domain with gobj
/gobj_maskaimpeel_ow		
	-input_output			commontrim
	-gobj_filename			trab_gobj1_file
	-peel_iter				0

/bounding_box_cut
	-input					commontrim
	-output					commontrimt
	-z_only					false
	-border					0 0 0
	
/togobj_from_aim
	-input						commontrimt
	-gobj_filename					"BASEtrab
	-min_elements					0
	-max_elements					0
	-curvature_smooth				1
	
/delete
	-name					commontrimt
	

! =================== 2nd image ===========================	

/read
	-name 					common06
	-filename					"COMM2
	
/set_value
	-input					common06
	-value_object			127
	-value_background		0
	
/bounding_box_cut
	-input					common06
	-output					common06trim
	-z_only					false
	-border					0 0 0

/delete
	-name					common06
	
! Now cut the commmon region in each domain with gobj
/gobj_maskaimpeel_ow		
	-input_output			common06trim
	-gobj_filename			auto_gobj2_file
	-peel_iter				0

/bounding_box_cut
	-input					common06trim
	-output					common06trimt
	-z_only					false
	-border					0 0 0
	
/delete
	-name					common06trim
	
/togobj_from_aim
	-input						common06trimt
	-gobj_filename					"FOLLOW06auto
	-min_elements					0
	-max_elements					0
	-curvature_smooth				1
	
! ==============================================
! Now cut the commmon region in each domain with gobj
/gobj_maskaimpeel_ow		
	-input_output			common06trim
	-gobj_filename			crtx_gobj2_file
	-peel_iter				0

/bounding_box_cut
	-input					common06trim
	-output					common06trimt
	-z_only					false
	-border					0 0 0
	
/togobj_from_aim
	-input						common06trimt
	-gobj_filename					"FOLLOW06crtx
	-min_elements					0
	-max_elements					0
	-curvature_smooth				1
	
! =========================================================
	
! Now cut the commmon region in each domain with gobj
/gobj_maskaimpeel_ow		
	-input_output			common06trim
	-gobj_filename			crtxpo_gobj2_file
	-peel_iter				0

/bounding_box_cut
	-input					common06trim
	-output					common06trimt
	-z_only					false
	-border					0 0 0
	
/togobj_from_aim
	-input						common06trimt
	-gobj_filename					"FOLLOW06crtxpo
	-min_elements					0
	-max_elements					0
	-curvature_smooth				1	
! =========================================================	
! Now cut the commmon region in each domain with gobj
/gobj_maskaimpeel_ow		
	-input_output			common06trim
	-gobj_filename			trab_gobj2_file
	-peel_iter				0

/bounding_box_cut
	-input					common06trim
	-output					common06trimt
	-z_only					false
	-border					0 0 0
	
/togobj_from_aim
	-input						common06trimt
	-gobj_filename					"FOLLOW06trab
	-min_elements					0
	-max_elements					0
	-curvature_smooth				1

/delete
	-name					common06trimt	

..
$ exit
