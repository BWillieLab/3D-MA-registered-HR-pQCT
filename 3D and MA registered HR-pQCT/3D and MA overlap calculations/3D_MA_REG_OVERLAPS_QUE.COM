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
$!				1. Calls the "3D_MA_REG_OVERLAPS.COM" script, then opens the log file containing the overlap%s, and write the results to a CSV file
$!  * Please note that the output numbers will be integers since floating numbers cannot be handled using IPL
$!
$! INPUTS:
$!              a. The grayscale ".AIM" file of the 1st image
$!              b. The grayscale ".AIM" file of the 2nd image
$!              c. The top slice number (indexed from 1) of the flattened 3D-registered region
$!              d. The top slice number (indexed from 1) of the flattened 3D-registered region
$!              e. The top slice number (indexed from 1) of the flattened MA-registered region
$!              f. The top slice number (indexed from 1) of the flattened MA-registered region
$
$! OUTPUTS: Calculated numbers in a csv file
$!				A. The overlap % of the 3D-registered and unregistered masks
$!              B. The overlap % of the flattened 3D-registered and unregistered masks
$!              C. The overlap % of the 3D-registered common voxels relative to all voxels in segmented image 1
$
$! call as:
$! @3D_MA_REG_OVERLAPS_QUE.com    c0000001.aim    c0000002.aim		TOP_3D 	   BOTTOM_3D      TOP_MA 	   BOTTOM_MA
$!==========================================================================
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
$
$ PATNO = F$EXTRACT(15,8,DIR_READ)
$
$ CHECK_FILE :== 'DEV''DIR_READ'CHECK.TXT
$
$ QUEUE_NAME :== SYS$SLOW
$ LOG_VOX :== 'DEV''DIR_READ''FILENAME'_3DREGISTER_2SCANS_OVERLAP.LOG
$ SCRIPT_VOX = "DISK1:[RESEARCH]3D_REG_OVERLAPS.COM"
$ INPUT1 = P1
$ INPUT2 = P2
$ INPUT3 = P3
$ INPUT4 = P4
$ INPUT5 = P5
$ INPUT6 = P6
$ SUBMIT/QUEUE='QUEUE_NAME'/PARAM=('INPUT1','INPUT2','INPUT3','INPUT4','INPUT5','INPUT6')/LOG=('LOG_VOX') 'SCRIPT_VOX'
$
$ WRITE_FILE :== "DISK1:[RESEARCH.RESULTS]3D_MA_OVERLAPS.CSV"   ! the text file in which the overlap of 2 3D-registered scans is appended
$
$ DIRECT :== "DISK1:[RESEARCH]"
$ SET DEFAULT 'DIRECT'
$
$! First, check if the log file is created
$CHECK2:
$   IF F$SEARCH(CHECK_FILE) .NES. "" THEN GOTO EXITCHECK2
$   GOTO CHECK2
$EXITCHECK2:
$
$ DIRECT :== 'DEV''DIR_READ'
$ SET DEFAULT 'DIRECT'
$ DELETE CHECK.TXT;*
$ DIRECT :== "DISK1:[RESEARCH]"
$ SET DEFAULT 'DIRECT'
$! Adding some idle time to make sure that the log file is closed by the program, otherwise, we will get an error
$ i = 0
$ WRITE SYS$OUTPUT ""
$ WRITE SYS$OUTPUT "========= NOW, OPENING THE LOG FILE TO READ THE OVERLAPS ========="
$IDLE:
$   IF i .EQ. 1000000 THEN GOTO EXITIDLE
$   i = i+1
$   GOTO IDLE
$EXITIDLE:
$
$! Now, opening the log file to read the number of voxels
$ CLOSE/NOLOG FOO
$ OPEN/READ FOO 'LOG_VOX'
$! Iterating through the lines of the log file to find the line that contains the number of voxels
$LOOP:
$	READ FOO LINE
$   IF F$LOCATE("obj_nr",LINE) .EQ. 3 THEN GOTO EXITLOOP 
$   GOTO LOOP
$EXITLOOP:
$
$ ! Now, extracting the number of voxels and converting it to an integer
$ NUM = F$EXTRACT(15,30,LINE)
$ UNREG_MASK = F$INTEGER(NUM)
$
$! Continue reading the file for the number fo voxels of the second AIM file
$! 3D-registration
$LOOP1:
$	READ FOO LINE
$   IF F$LOCATE("obj_nr",LINE) .EQ. 3 THEN GOTO EXITLOOP1 
$   GOTO LOOP1
$EXITLOOP1:
$
$ NUM = F$EXTRACT(15,30,LINE)
$ REG_MASK_TD = F$INTEGER(NUM)
$
$LOOP2:
$	READ FOO LINE
$   IF F$LOCATE("obj_nr",LINE) .EQ. 3 THEN GOTO EXITLOOP2
$   GOTO LOOP2
$EXITLOOP2:
$
$ NUM = F$EXTRACT(15,30,LINE)
$ REG_MASK_FLAT_TD = F$INTEGER(NUM)
$
$LOOP3:
$	READ FOO LINE
$   IF F$LOCATE("obj_nr",LINE) .EQ. 3 THEN GOTO EXITLOOP3
$   GOTO LOOP3
$EXITLOOP3:
$
$ NUM = F$EXTRACT(15,30,LINE)
$ REG_SEG1 = F$INTEGER(NUM)
$
$LOOP4:
$	READ FOO LINE
$   IF F$LOCATE("obj_nr",LINE) .EQ. 3 THEN GOTO EXITLOOP4
$   GOTO LOOP4
$EXITLOOP4:
$
$ NUM = F$EXTRACT(15,30,LINE)
$ REG_SEG_OVERLAY_TD = F$INTEGER(NUM)
$
$! MA-registration
$LOOP5:
$	READ FOO LINE
$   IF F$LOCATE("obj_nr",LINE) .EQ. 3 THEN GOTO EXITLOOP5 
$   GOTO LOOP5
$EXITLOOP5:
$
$LOOP6:
$	READ FOO LINE
$   IF F$LOCATE("obj_nr",LINE) .EQ. 3 THEN GOTO EXITLOOP6 
$   GOTO LOOP6
$EXITLOOP6:
$
$ NUM = F$EXTRACT(15,30,LINE)
$ REG_MASK_MA = F$INTEGER(NUM)
$
$LOOP7:
$	READ FOO LINE
$   IF F$LOCATE("obj_nr",LINE) .EQ. 3 THEN GOTO EXITLOOP7
$   GOTO LOOP7
$EXITLOOP7:
$
$LOOP8:
$	READ FOO LINE
$   IF F$LOCATE("obj_nr",LINE) .EQ. 3 THEN GOTO EXITLOOP8
$   GOTO LOOP8
$EXITLOOP8:
$
$LOOP9:
$	READ FOO LINE
$   IF F$LOCATE("obj_nr",LINE) .EQ. 3 THEN GOTO EXITLOOP9
$   GOTO LOOP9
$EXITLOOP9:
$
$ NUM = F$EXTRACT(15,30,LINE)
$ REG_MASK_FLAT_MA = F$INTEGER(NUM)
$
$LOOP10:
$	READ FOO LINE
$   IF F$LOCATE("obj_nr",LINE) .EQ. 3 THEN GOTO EXITLOOP10
$   GOTO LOOP10
$EXITLOOP10:
$
$LOOP11:
$	READ FOO LINE
$   IF F$LOCATE("obj_nr",LINE) .EQ. 3 THEN GOTO EXITLOOP11
$   GOTO LOOP11
$EXITLOOP11:
$
$LOOP12:
$	READ FOO LINE
$   IF F$LOCATE("obj_nr",LINE) .EQ. 3 THEN GOTO EXITLOOP12
$   GOTO LOOP12
$EXITLOOP12:
$
$LOOP13:
$	READ FOO LINE
$   IF F$LOCATE("obj_nr",LINE) .EQ. 3 THEN GOTO EXITLOOP13
$   GOTO LOOP13
$EXITLOOP13:
$
$LOOP14:
$	READ FOO LINE
$   IF F$LOCATE("obj_nr",LINE) .EQ. 3 THEN GOTO EXITLOOP14
$   GOTO LOOP14
$EXITLOOP14:
$
$
$ NUM = F$EXTRACT(15,30,LINE)
$ REG_SEG1_MA = F$INTEGER(NUM)
$
$LOOP15:
$	READ FOO LINE
$   IF F$LOCATE("obj_nr",LINE) .EQ. 3 THEN GOTO EXITLOOP15
$   GOTO LOOP15
$EXITLOOP15:
$
$
$ NUM = F$EXTRACT(15,30,LINE)
$ REG_SEG_OVERLAY_MA = F$INTEGER(NUM)
$
$
$! Calculating the volume fraction of the common region relative to the mask of the baseline scan
$! MA-registration
$ MSKDIV = UNREG_MASK/100
$ OVERLAP_TD = REG_MASK_TD/MSKDIV
$ LAPSTR_TD = F$STRING(OVERLAP_TD)
$
$ OVERLAP_FLAT_TD = REG_MASK_FLAT_TD/MSKDIV
$ LAPSTR_FLAT_TD = F$STRING(OVERLAP_FLAT_TD)
$
$ SEGDIV = REG_SEG1/100
$ OVERLAY_TD = REG_SEG_OVERLAY_TD/SEGDIV
$ LAPSTR_OVERLAY_TD = F$STRING(OVERLAY_TD)
$
$! MA-registration
$ OVERLAP_MA = REG_MASK_MA/MSKDIV
$ LAPSTR_MA = F$STRING(OVERLAP_MA)
$
$ OVERLAP_FLAT_MA = REG_MASK_FLAT_MA/MSKDIV
$ LAPSTR_FLAT_MA = F$STRING(OVERLAP_FLAT_MA)
$
$ SEGDIV_MA = REG_SEG1_MA/100
$ OVERLAY_MA = REG_SEG_OVERLAY_MA/SEGDIV_MA
$ LAPSTR_OVERLAY_MA = F$STRING(OVERLAY_MA)
$
$ TIME = F$TIME()
$
$! Concatenating the patient information and overlap % to be appended into the text file
$ CAT1 = PATNO + "," + FILENAME + "_" + FILENAME2 + "," + TIME + "," + LAPSTR_TD + "," + LAPSTR_FLAT_TD + "," + LAPSTR_OVERLAY_TD
$ CAT = CAT1 + "," + LAPSTR_MA + "," + LAPSTR_FLAT_MA + "," + LAPSTR_OVERLAY_MA
$
$! Open the output text file
$ CLOSE/NOLOG OUT
$ OPEN/APPEND OUT 'WRITE_FILE'
$ WRITE OUT CAT
$ CLOSE OUT
$
$ DIRECT :== "DISK1:[RESEARCH]"
$ SET DEFAULT 'DIRECT'
$ EXIT