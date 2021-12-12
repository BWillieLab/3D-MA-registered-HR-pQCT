# -*- coding: utf-8 -*-
"""
$!========================================================================================
$! Developed by Seyedmahdi(Mahdi) Hosseinitabatabaei (mahdi.tabatabaei@mail.mcgill.ca; https://www.researchgate.net/profile/Seyedmahdi-Hosseinitabatabaei)
$! Dr. B.Willie's lab, Shriners Hospital for Children - McGill University, Montreal, Canada
$! Version P1 (07-DEC-2021)
$!========================================================================================
$! For the 3D-TB model, bscan specific boundary condition files (.BCD) are needed to align the loading vectors with each scan
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
$!              a. The list of scans with their XCT number, measurements number, rotation angles
$!              b. Template BCD file
$
$! OUTPUTS:
$!              A. Scan specific BCD files
$
$! Version history:
$!				N/A
$!========================================================================================
"""


import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import csv
from decimal import Decimal

# First, reading csv file containing the list of scans with their rotation angles 
scan_list = "C:\\folder\\sub folder\\Scans_rotation_angles.csv"

# initializing the titles and rows list 
fields = [] 
rows = [] 

with open(scan_list, 'r') as csvfile: 
    # creating a csv reader object 
    csvreader = csv.reader(csvfile) 
      
    # extracting field names through first row 
    fields = next(csvreader) 
  
    # extracting each data row one by one 
    for row in csvreader: 
        rows.append(row) 


# Extract the scan number and rotation angles
for find in rows:
    xct_num = str(find[0])
    meas_num = str(find[1])
    vox = float(find[4])
    stat = str(find[5])
    if vox == 82:
        total_len = 110*82/1000
    else:
        total_len = 168*60.7/1000
    
    if stat == 'ref':
        rot_x = 0
        rot_y = 0
        rot_z = 0
    elif stat == 'mov':
        rot_x = float(find[6])
        rot_y = float(find[7])
        rot_z = float(find[8])

    app_strain = -0.01
    
    
    # Transform the boundary conditions using the rotation angles
    rot_x = np.radians(rot_x)
    rot_y = np.radians(rot_y)
    rot_z = np.radians(rot_z)
    
    cx, sx = np.cos(rot_x), np.sin(rot_x)
    Rx = np.array(((1,0,0),(0, cx, -sx), (0,sx, cx)))

    cy, sy = np.cos(rot_y), np.sin(rot_y)
    Ry = np.array(((cy,0, sy), (0,1,0),(-sy,0, cy)))

    cz, sz = np.cos(rot_z), np.sin(rot_z)
    Rz = np.array(((cz, -sz,0), (sz, cz,0),(0,0,1)))

    # Inverse transformation
    Rx1 = np.linalg.inv(Rx)
    Ry1 = np.linalg.inv(Ry)
    Rz1 = np.linalg.inv(Rz)
    
    R = np.dot(Rx1,np.dot(Ry1,Rz1))      # Transformation matrix (inverse)
    
    v0 = np.array([0,0,app_strain*total_len])
    
    
    v = np.dot(R,v0)       # Transform the loading vector

    # Create the final loading vector
    u1 = v[0]    # x-component
    u2 = v[1]    # y-component
    u3 = v[2]    # z-component
    
    # Identify the writing starting index for each vector component
    if u1 >=0:
        u1_index = 43
        gap1 = 10
    elif u1 <0:
        u1_index = 42
        gap1 = 11
        
    if u2 >=0:
        u2_index = 67
        gap2 = 10
    elif u2 <0:
        u2_index = 66
        gap2 = 11
        
    if u3 >=0:
        u3_index = 91
        gap3 = 10
    elif u3 <0:
        u3_index = 90
        gap3 = 11

    # convert floating numbers to exponential notation with 4 digits right from the decimal point
    u1_str ='%.4E' % Decimal(u1)
    u2_str ='%.4E' % Decimal(u2)
    u3_str ='%.4E' % Decimal(u3)

    # Convert the E notation to d to be compatible with the BCD template
    u1_str = u1_str.replace("E","d")
    u2_str = u2_str.replace("E","d")
    u3_str = u3_str.replace("E","d")
    
    # Open the template BCD file
    bcd_temp = open("C:\\folder\\sub folder\\Template_BDC.BCD", 'r')
    temp_lines = bcd_temp.readlines()

    # Replace the BCD file components using the calculated vector
    
    n = 1
    replacement = 'C'
    # Replace character at nth index position
    temp_lines[7] = temp_lines[7][0:u1_index] + u1_str + temp_lines[7][u1_index+gap1: ]   #u1
    temp_lines[8] = temp_lines[8][0:u2_index] + u2_str + temp_lines[8][u2_index+gap2: ]   #u2
    temp_lines[9] = temp_lines[9][0:u3_index] + u3_str + temp_lines[9][u3_index+gap3: ]   #u3

    # writing to file
    bcd_filename = str(xct_num) + "_" + str(meas_num) + "_FE_3D_TB.BCD"
    write_bcd = open("C:\\folder\\sub folder\\BCDs\\" + bcd_filename, 'w')
    write_bcd.writelines(temp_lines)
    write_bcd.close()
