# -*- coding: utf-8 -*-
"""
$!========================================================================================
$! Developed by Seyedmahdi(Mahdi) Hosseinitabatabaei (mahdi.tabatabaei@mail.mcgill.ca; https://www.researchgate.net/profile/Seyedmahdi-Hosseinitabatabaei)
$! Dr. B.Willie's lab, Shriners Hospital for Children - McGill University, Montreal, Canada
$! Version P1 (07-DEC-2021)
$!========================================================================================
$! For the 3D-TB model, during post-processing in python, reaction forces need to be projected onto the unit vector representing the direction of applied compressive load and aggregated to obtain the total reaction force.
$! Finally, stiffness and failure load were calculated.
$!========================================================================================
$! Important Parameters:
$!
$!========================================================================================
$! This script:
$!				1. Parses required data from te LISTING and POSTLIST files to calculate microFE stiffness, and failure load
$
$! INPUTS:
$!              a. The list of scans with their XCT number, measurements number, rotation angles
$!              b. Complete LISTING files
$!              c. Complete POSTLIST files
$
$! OUTPUTS:
$!              A. MicroFE outcomes written to a CSV file
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
import time
import datetime


start = time.time()
filename = "C:\\folder\\sub folder\\Scans_rotation_angles.csv"
 
# initializing the titles and rows list 
fields = [] 
rows = [] 
  
# reading csv file 
with open(filename, 'r') as csvfile: 
    # creating a csv reader object 
    csvreader = csv.reader(csvfile) 
      
    # extracting field names through first row 
    fields = next(csvreader) 
  
    # extracting each data row one by one 
    for row in csvreader: 
        rows.append(row) 
  
#print(rows)
for find in rows:
#    find = [row for row in rows if (str(xct) in row[0]) and (str(meas) in row[1])]

    xct_num = str(find[0])
    meas_num = str(find[1])
    vox = float(find[4])
    stat = str(find[5])
    
    if stat == 'ref':
        rot_x = 0
        rot_y = 0
        rot_z = 0
    elif stat == 'mov':
        rot_x = float(find[6])
        rot_y = float(find[7])
        rot_z = float(find[8])

    listing_file_name = xct_num + "_" + meas_num + "_FE_3D_TB.LISTING"
    postlist_file_name = xct_num + "_" + meas_num + "_FE_3D_TB.POSTLIST"


    listing_file = "C:\\folder\\sub folder\\listing_and_postlist\\" + listing_file_name
    postlist_file = "C:\\folder\\sub folder\\listing_and_postlist\\" + postlist_file_name

    
    # Find the number of nodes and slices in the model
    fp = open(listing_file)
    for i, line in enumerate(fp):
        if i == 32:
            n_slices = int(line[-4:-1])
            if n_slices == 110:
                vox = 82
                E_mod = 6829
            elif n_slices == 168:
                vox = 60.7
                E_mod = 8748
            total_len = vox*n_slices
        elif i == 60:
            n_elems = int(line[-10:-1])
        elif i == 61:
            n_nodes = int(line[-10:-1])
            break
    fp.close()
    
    
    app_strain = -0.01
    
    print("identifying the loading vector")
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
    
    v0 = np.array([0,0,app_strain*total_len/1000])
    
    
    v = np.dot(R,v0)       # Transform the loading vector
    
    # Create the final loading vector
    u1 = v[0]    # x-component
    u2 = v[1]    # y-component
    u3 = v[2]    # z-component
    
    
    u_mag = math.sqrt(u1**2+u2**2+u3**2)    # loading vector magnitude
       
    # Loading unit vector
    u1_norm = u1/u_mag
    u2_norm = u2/u_mag
    u3_norm = u3/u_mag
    
    print("Calculating the total reaction force alogn the loading vector")
    # Read the reaction force components of the loading surface
    RF_tot = 0
    fp = open(listing_file)
    for i, line in enumerate(fp):
        if (line[12:15] == str(n_slices)) and (line[98:101].isdigit()):
            rf_x = float(line[60:75])
            rf_y = float(line[77:91])
            rf_z = float(line[92:107])
            # Project the force components to the loading direction
            rf_x *= u1_norm
            rf_y *= u2_norm
            rf_z *= u3_norm
            
            # Calculate the magnitude of the projected reaction force
            RF_sum = rf_x + rf_y + rf_z
            # Add the force to the total force count
            RF_tot = RF_tot + RF_sum
    fp.close()
    
    print("reading elemental energy equivalent strains")
    
    # Now, read the postlist file to obtain energy equivalent strain values
    edens = []
    fp = open(postlist_file)
    for i, line in enumerate(fp):
        if (line[5:6] == "2") and (line[28:31].isdigit()):
            eseden = float(line[23:38])   # strain energy density
            edens.append(eseden)
    fp.close()
    
    edens = np.array(edens)
    
    print("Final calculations...")
    
    # Calculate stiffness
    Stiffness = -RF_tot/u_mag
    
    
    n_elems = len(edens)
    
    v_crit = 0.02      # Pistoia's citerior: critical volume of 2%
      
    # sort the data in ascending order 
    ordered = np.sort(edens) 
      
    # get the cdf values of energy equivalent strain 
    cdf = np.arange(n_elems) / float(n_elems) 
    
    u_thresh=ordered[abs((cdf-(1-v_crit))) <= 0.00005]
    
       
    u_thresh = u_thresh[0]
    
    EES = u_thresh
    EES_crit = 0.007    # Pistoia's citerior: critical value of 7000 ue
    factor = EES_crit/EES
    
    failure_load = factor*RF_tot
    
    
    to_write = str(xct_num) +','+ str(meas_num) +','+stat +','+ str(vox) +','+str(find[6])+\
    ','+str(find[7])+','+str(find[8])+','+str(E_mod)+','+str(u1)+','+str(u2)+','+str(u3)+','+str(u_mag)+','+str(RF_tot)+','+\
    ','+str(u_thresh)+','+str(EES)+','+str(EES_crit)+','+str(v_crit)+','+str(factor)+','+str(app_strain)+','\
    +str(Stiffness)+','+str(failure_load)+','+','+str(datetime.datetime.now())+'\n'

    file = open('"C:\\folder\\sub folder\\listing_and_postlist\\ipl_3D_tb_results.csv','a')
    file.write(to_write)
    file.close()

    
    print('DONE!')
    print(datetime.datetime.now())


