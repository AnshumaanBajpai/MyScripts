#!/afs/crc.nd.edu/x86_64_linux/python/2.7.4/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 11:47:30 2013

@author: abajpai1
"""

##################################################################################
##################################################################################
##               This code was written by Anshumaan Bajpai                      ##
##    during the Christmas week of 2013 when rest of the world was celebrating  ##
##   and should be remembered for the devotion and commitment of the student.   ##
##    This code would have been impossible without Mike Penninger whose         ##
############   encouragement and technical help made this possible.       ########
##################################################################################
##################################################################################
#############   The code reads DOSCAR files for following systems:   #############
#############                    LORBIT = 10/11                      #############
#############          spin polarized and non spin polarized.        #############
##################################################################################
##################################################################################
#############   It takes input from the user about the atoms and     #############
#############   orbitals needed to be analyzed and produces output   #############
#############         files and plots with the required data         #############
##################################################################################
##################################################################################

import os
import  string
import numpy as np
import  matplotlib.pyplot as plt
#from matplotlib.pyplot import *
#from matplotlib.patches import *
#from scipy import *
#import scipy.optimize as optimize
#import pylab as pyplt
import shutil
from shutil import rmtree


###################################################################################################################################################################################################
################################################################### Define Function to read POSCAR lines ##############################################
###################################################################################################################################################################################################



def posinfo():								# Gathers information form POSCAR
    contcar=open('CONTCAR', 'r')						# Next few lines reads all the initial data of POSCAR. Title/Scale factor/lattice parameters
    name=contcar.readline()						#
    scalefac=contcar.readline()						#
    lata=contcar.readline()						#
    latb=contcar.readline()						#
    latc=contcar.readline()						#
    temp=contcar.readline()
    #print temp
    if not temp.split()[0].isdigit:					# Why will it not be digit!!! Somehing fishy					
        #print "is digit"
        atomnum =temp
        atomname= []
    else:
        #print "is not digit"
        atomname=temp
        atomnum=contcar.readline()
    temp=contcar.readline().lower()	
    print atomname
    print atomnum
    contcar.close()



###################################################################################################################################################################################################
################################################################### This section of the program splits the DOSCAR and gives file for individual atoms##############################################
###################################################################################################################################################################################################
dos_to_split = open("DOSCAR","r")
count = 1
while count<6:
    temp = dos_to_split.readline()
    count +=1

details = dos_to_split.readline().split()
en_lev = int(details[2])
fermi_level = float(details[3])
#print en_lev

count = 1
while count <=en_lev:
    temp = dos_to_split.readline()
    #print temp
    count +=1

cwd = os.getcwd()
nwd = cwd + "/split_DOSCAR"
if os.path.isdir(nwd):
    #os.chdir(nwd)
    shutil.rmtree(nwd)
    #os.rmdir(nwd)
    
os.mkdir(nwd)
os.chdir(nwd)
   
atom_no = 1   
while True:
    temp = dos_to_split.readline()
    if temp == "":
        break
    elif len(temp.split()) == 5 and int(temp.split()[2]) == en_lev and float(temp.split()[3]) == fermi_level:
        file1 = open("Atom_"+str(atom_no),"w")
        count = 1
        while count <= en_lev:
            temp = dos_to_split.readline()
            file1.write(temp)
            count +=1
            
        file1.close()
        atom_no +=1
        
dos_to_split.close()
            
            
os.chdir(cwd)

#print os.getcwd() 
##################################################################################################################################################################################################           
###############################################################################  Atom splitting complete  ########################################################################################            
##################################################################################################################################################################################################    


##################################################################################################################################################################################################
############################################################################### function to do the plotting in case of a single set of outputs####################################################
##################################################################################################################################################################################################
def plot_single(A):
    fol_name= raw_input("Enter a name for the folder for saving the plot(s) for this calculation.\n")
    print "************************************************************************************************************"
    os.mkdir(cwd+"/"+fol_name)
    os.chdir(cwd+"/"+fol_name)
    plt.figure(fol_name)
    plt.plot(Atom_matrix[1][0],A,"r-")
    plt.title("Density of states vs Fermi normalized energy")
    plt.xlabel("Energy (eV)")
    plt.ylabel("No. of states (a.u.)")
    plt.xlim()
    plt.ylim()
    plt.axhline(y=0, color = 'k')
    plt.axvline(x=0, color = 'k')
    plt.savefig("total_DOS_des"+".pdf")
                                
    os.chdir(cwd)    
##################################################################################################################################################################################################
############################################################################### function to do the plotting in case of a single set of outputs complete ##########################################
##################################################################################################################################################################################################


##################################################################################################################################################################################################
############################################################################### function to do the plotting in case of a multiple outputs  #######################################################
##################################################################################################################################################################################################
def plot_multiple(total_DOS_s,total_DOS_s_up,total_DOS_s_down,total_DOS_p,total_DOS_p_up,total_DOS_p_down,total_DOS_d,total_DOS_d_up,total_DOS_d_down,total_DOS_f,total_DOS_f_up,total_DOS_f_down,total_DOS,total_DOS_up,total_DOS_down):
    plot = [total_DOS_s,total_DOS_s_up,total_DOS_s_down,total_DOS_p,total_DOS_p_up,total_DOS_p_down,total_DOS_d,total_DOS_d_up,total_DOS_d_down,total_DOS_f,total_DOS_f_up,total_DOS_f_down,total_DOS,total_DOS_up,total_DOS_down]
    filenames = ["s","p","d","f","total"]
    fol_name= raw_input("Enter a name for the folder for saving the plot(s) for this calculation.\n")
    print "************************************************************************************************************"
    os.mkdir(cwd+"/"+fol_name)
    os.chdir(cwd+"/"+fol_name)
    for pl_cv in range((len(plot))/3):
        plt.figure(pl_cv)
        plt.plot(matrix[1][0],plot[3*pl_cv+1],"g-",matrix[1][0],plot[3*pl_cv+2],"b-")
        plt.title("Density of states vs Fermi normalized energy")
        plt.xlabel("Energy (eV)")
        plt.ylabel("No. of states (a.u.)")
        plt.xlim()
        plt.ylim()
        plt.axhline(y=0, color = 'k')
        plt.axvline(x=0, color = 'k')
        plt.legend((filenames[pl_cv]+"_up",filenames[pl_cv]+"_down"),loc="best")
        plt.savefig(filenames[pl_cv]+"up_down.pdf")
        
        
        plt.figure(pl_cv+50)
        plt.plot(matrix[1][0],plot[3*pl_cv])
        plt.title("Density of states vs Fermi normalized energy")
        plt.xlabel("Energy (eV)")
        plt.ylabel("No. of states (a.u.)")
        plt.xlim()
        plt.ylim()
        plt.axhline(y=0, color = 'k')
        plt.axvline(x=0, color = 'k')
        #plt.legend((filenames[pl_cv],filenames[pl_cv]+"_up",filenames[pl_cv]+"_down"),loc=4)
        plt.savefig(filenames[pl_cv]+".pdf")
        
        pl_cv +=1
        
    os.chdir(cwd)    
##################################################################################################################################################################################################
############################################################################### function to do the plotting in case of a multiple outputs complete#######################################################
##################################################################################################################################################################################################


##################################################################################################################################################################################################
############################################################################### orbital assignment for 101/102#######################################################
##################################################################################################################################################################################################

def orbital_assignment_10():
    print "orbital =\ts_up\ts_down\tp_up\tp_down\td_up\td_down\tf_up\tf_down"
    print "Ref.no. =\t1\t2\t3\t4\t5\t6\t7\t8" 
    
##################################################################################################################################################################################################
############################################################################### orbital assignment for 101/102 complete#######################################################
##################################################################################################################################################################################################
#################################################################################################################################################################################################
############################################################################### orbital assignment for 111/112#######################################################
##################################################################################################################################################################################################
def orbital_assignment_11():
    print "s: 1 to 2\np: 3 to 8\nd: 9 to 18\nf:19 to 32\nFor detailed description refer to the README file."
   
##################################################################################################################################################################################################
############################################################################### orbital assignment for 111/112 complete#######################################################
##################################################################################################################################################################################################
   
   
   
print "LORBIT = LB, SPIN POLARIZED = SP"      #Introduce user to the terminology for the next step
#if terminology == "C":
system = raw_input("Identify your system? \nif LB = 11 and SP = 1, enter 111 \nif LB = 11 and SP = 2, enter 112 \nif LB = 10 and SP = 1, enter 101 \nif LB = 10 and SP = 2, enter 102\nAny other number/letter to exit\n")        #Gathers information about the type of system to be analyzed
print "************************************************************************************************************"


if system == "111" or system == "112":
    print "Your system is being saved as a matrix, this may take a moment...or 2"
    print "####################################################################################################\n####################################################################################################"
    def atom_111_112():     # This function creates a matrix atom =[0,Atom 1, Atom 2,.......], where each atom is a matrix with [0, sup, sdown,p_zup.....]
        doscar = open("DOSCAR","r")
            
            
        i = 1       #Next 4 lines are to skip the general inforamtion and arrive at the line that contains parameters 
        while i<6:
            temp = doscar.readline()
                #print temp
            i += 1
            
                
        parameters  = doscar.readline().split()     #Read the line that contains paremters
        E_max       = float(parameters[0])          #Maximum energy in the system
        E_min       = float(parameters[1])          #Minimum energy in the system
        #E_scaled     = E_min - 0.1*(E_max - E_min)
        nl          = float(parameters[2])          #No. of energy levels in the system
        E_fermi     = float(parameters[3])          #Fermi energy of the system
            #print E_max, E_min, nl, E_fermi
        
        
        NOA = 0     #Next 8 lines count the number of atoms in the system
        while True:
            temp = doscar.readline()
            if temp == "":
                break
            elif len(temp.split()) == 5 and float(temp.split()[2]) == nl and float(temp.split()[3]) == E_fermi:
                    NOA = NOA +1
            else:
                continue
            #print NOA
        doscar.close()
    

        doscar = open("DOSCAR","r")     #In the next 26 lines, we create an array for energy and Fermi level normalized energy. We also figure out the no of lines below fermilevel 
        i = 1
        E = []      #Energy list
        E_f = []    #Fermi Level normalized energy list
        while i < 7:
            temp = doscar.readline()
            #print temp
            i += 1
        i = 1
        while i <= nl:
            temp = doscar.readline()
            temp_E = [float(temp.split()[0])]
            temp_E_f = [float(temp.split()[0])-E_fermi]
            E += temp_E
            E_f += temp_E_f
            i += 1
        #print E
        
            
        atom = [E]      #In the next 90 lines, we try to create the atom matrix with first element as a list of energy,next elements are Atom1, Atom2 and so on. Each Atom is a list with elements as list E_f, sup, sdown, and so on.
            #print atom
        i = 1
        while i <= NOA:
            
            A = [E_f]
            temp = doscar.readline()
            sup           = []
            sdown         = []
            p_yup         = []
            p_ydown       = []
            p_zup         = []
            p_zdown       = []
            p_xup         = []
            p_xdown       = []
            d_xyup        = []
            d_xydown      = []
            d_yzup        = []
            d_yzdown      = []
            d_z2up        = []
            d_z2down      = []
            d_xzup        = []
            d_xzdown      = []
            d_x2y2up      = []
            d_x2y2down    = []
            f_y3x2y2up    = []
            f_y3x2y2down  = []
            f_xyzup       = []
            f_xyzdown     = []
            f_yz2up       = []
            f_yz2down     = []
            f_z3up        = []
            f_z3down      = []
            f_xz2up       = []
            f_xz2down     = []
            f_zx2y2up     = []
            f_zx2y2down   = []
            f_xx23y2up    = []
            f_xx23y2down  = []
                
            j = 1
            if system == "111":
                while j <= nl:
                    temp = doscar.readline().split()
                    a = 17-len(temp)
                    for k in range(a):
                        temp = temp + ["0"]
                            
                    sup         = sup + [0.5*float(temp[1])]
                    sdown       = sdown + [-0.5*float(temp[1])]
                    p_yup       = p_yup + [0.5*float(temp[2])]
                    p_ydown     = p_ydown + [-0.5*float(temp[2])]
                    p_zup       = p_zup + [0.5*float(temp[3])]
                    p_zdown     = p_zdown + [-0.5*float(temp[3])]
                    p_xup       = p_xup + [0.5*float(temp[4])]
                    p_xdown     = p_xdown + [-0.5*float(temp[4])]
                    d_xyup      = d_xyup + [0.5*float(temp[5])]
                    d_xydown    = d_xydown + [-0.5*float(temp[5])]
                    d_yzup      = d_yzup + [0.5*float(temp[6])]
                    d_yzdown    = d_yzdown + [-0.5*float(temp[6])]
                    d_z2up      = d_z2up + [0.5*float(temp[7])]
                    d_z2down    = d_z2down + [-0.5*float(temp[7])]
                    d_xzup      = d_xzup + [0.5*float(temp[8])]
                    d_xzdown    = d_xzdown + [-0.5*float(temp[8])]
                    d_x2y2up    = d_x2y2up + [0.5*float(temp[9])]
                    d_x2y2down  = d_x2y2down + [-0.5*float(temp[9])]
                    f_y3x2y2up    = f_y3x2y2up + [0.5*float(temp[10])]
                    f_y3x2y2down  = f_y3x2y2down + [-0.5*float(temp[10])]
                    f_xyzup       = f_xyzup +[0.5*float(temp[11])]
                    f_xyzdown     = f_xyzdown +[-0.5*float(temp[11])]
                    f_yz2up       = f_yz2up +[0.5*float(temp[12])]
                    f_yz2down     = f_yz2down +[-0.5*float(temp[12])]
                    f_z3up        = f_z3up +[0.5*float(temp[13])]
                    f_z3down      = f_z3down +[-0.5*float(temp[13])]
                    f_xz2up       = f_xz2up +[0.5*float(temp[14])]
                    f_xz2down     = f_xz2down +[-0.5*float(temp[14])]
                    f_zx2y2up     = f_zx2y2up +[0.5*float(temp[15])]
                    f_zx2y2down   = f_zx2y2down +[-0.5*float(temp[15])]
                    f_xx23y2up    = f_xx23y2up +[0.5*float(temp[16])]
                    f_xx23y2down  = f_xx23y2down +[-0.5*float(temp[16])]                  
                    j += 1
                    
            else:
                while j <= nl:
                    temp = doscar.readline().split()
                    a = 33-len(temp)
                    for k in range(a):
                        temp = temp + ["0"]
                        
                    sup         = sup + [0.5*float(temp[1])]
                    sdown       = sdown + [-0.5*float(temp[2])]
                    p_yup       = p_yup + [0.5*float(temp[3])]
                    p_ydown     = p_ydown + [-0.5*float(temp[4])]
                    p_zup       = p_zup + [0.5*float(temp[5])]
                    p_zdown     = p_zdown + [-0.5*float(temp[6])]
                    p_xup       = p_xup + [0.5*float(temp[7])]
                    p_xdown     = p_xdown + [-0.5*float(temp[8])]
                    d_xyup      = d_xyup + [0.5*float(temp[9])]
                    d_xydown    = d_xydown + [-0.5*float(temp[10])]
                    d_yzup      = d_yzup + [0.5*float(temp[11])]
                    d_yzdown    = d_yzdown + [-0.5*float(temp[12])]
                    d_z2up      = d_z2up + [0.5*float(temp[13])]
                    d_z2down    = d_z2down + [-0.5*float(temp[14])]
                    d_xzup      = d_xzup + [0.5*float(temp[15])]
                    d_xzdown    = d_xzdown + [-0.5*float(temp[16])]
                    d_x2y2up    = d_x2y2up + [0.5*float(temp[17])]
                    d_x2y2down  = d_x2y2down + [-0.5*float(temp[18])]
                    f_y3x2y2up    = f_y3x2y2up + [0.5*float(temp[19])]
                    f_y3x2y2down  = f_y3x2y2down + [-0.5*float(temp[20])]
                    f_xyzup       = f_xyzup +[0.5*float(temp[21])]
                    f_xyzdown     = f_xyzdown +[-0.5*float(temp[22])]
                    f_yz2up       = f_yz2up +[0.5*float(temp[23])]
                    f_yz2down     = f_yz2down +[-0.5*float(temp[24])]
                    f_z3up        = f_z3up +[0.5*float(temp[25])]
                    f_z3down      = f_z3down +[-0.5*float(temp[26])]
                    f_xz2up       = f_xz2up +[0.5*float(temp[27])]
                    f_xz2down     = f_xz2down +[-0.5*float(temp[28])]
                    f_zx2y2up     = f_zx2y2up +[0.5*float(temp[29])]
                    f_zx2y2down   = f_zx2y2down +[-0.5*float(temp[30])]
                    f_xx23y2up    = f_xx23y2up +[0.5*float(temp[31])]
                    f_xx23y2down  = f_xx23y2down +[-0.5*float(temp[32])]
                    j +=1

                      
              
            A = A + [sup] +[sdown] + [p_yup]+ [p_ydown] + [p_zup] + [p_zdown] + [p_xup] + [p_xdown]  + [d_xyup] + [d_xydown] + [d_yzup] +[d_yzdown] + [d_z2up] + [d_z2down] + [d_xzup] + [d_xzdown] + [d_x2y2up] + [d_x2y2down] + [f_y3x2y2up] + [f_y3x2y2down] + [f_xyzup] + [f_xyzdown] + [f_yz2up] + [f_yz2down] + [f_z3up] + [f_z3down] + [f_xz2up] + [f_xz2down] + [f_zx2y2up] + [f_zx2y2down] + [f_xx23y2up] + [f_xx23y2down]
            atom = atom + [A]
            i += 1
            
        doscar.close()
        return atom
        
        
        
        
    def nl_fermi_11():
        nl_f = 0
        i = 0
        while i < len(Atom_matrix[1][0]):
            if Atom_matrix[1][0][i]<=0:
                nl_f += 1
                i += 1
                  #print nl
            else:
                break
        return nl_f
            
              
    Atom_matrix = atom_111_112()
    nl_f = nl_fermi_11()


    print "Your system has been saved as a matrix. Let's do some analysis."
    print "The DOSCAR has been split into individual atoms and saved in the folder named split_DOSCAR in your current working directory.  You can refer to split_DOSCAR for any additional analysis.\n"
    print "####################################################################################################\n####################################################################################################"
        #print nl_f
        #print Atom_matrix[1][6]
        #print Atom_matrix[1][8]
        #print Atom_matrix[2][6]
        #print Atom_matrix[2][8]
        #print Atom_matrix[3][6]
        #print Atom_matrix[3][8]
        #print Atom_matrix[4][6]
        #print Atom_matrix[4][8]
        #print Atom_matrix[5][6]
        #print Atom_matrix[5][8]
        
        
        #cwd = os.getcwd()
        #nwd_name = raw_input("What would you like to name the new directory where the data will be saved?\n")
        #nwd = cwd + "/" + "Atom_data"
        #os.mkdir(nwd)
        #os.chdir(nwd)
    while True:
        query1 = raw_input("Do you wish to study properties for all the atoms in your system or only specific atoms in a particluar layer/adorbent etc.?\nEnter A for all atoms;\nEnter R for a range of atoms\nEnter S for specific atoms;\nEnter D if you are done with your analysis.\n")
        print "************************************************************************************************************"
        if query1 == "A" or query1 == "R":
            if query1 == "A":
                matrix = Atom_matrix
            else:
                while True:
                    lowlim = raw_input("Please enter the lower limit of the range.\n")
                    highlim = raw_input("Please enter the upper limit of the range.\n")
                    print "************************************************************************************************************"                            
                    if lowlim.isdigit() and highlim.isdigit():
                        low = int(lowlim)
                        high = int(highlim)
                        if 0<low<=high<len(Atom_matrix):
                            #print Atom_matrix[0]
                            matrix = [Atom_matrix[0]] + Atom_matrix[low:high+1]
                            break
                        else:
                            print "Make sure that the lower limit is less than the upper limit and both of them are less than or equal to the number of atoms in your system."
                    else:
                        print "Both the limits need to be integers greater than 0 and less than or equal to the number of atoms in the system and lower limit should be less than upper limit."
               #print matrix[1][1]
            while True:
                query2 = raw_input("Do you need total DOS for all orbitals?\nEnter A for all;\nEnter S for specific;\nEnter D if you are done.\n")
                print "************************************************************************************************************"
                if query2 == "A":
                    total_DOS_s_up = [0]*len(matrix[0])
                    total_DOS_p_up = [0]*len(matrix[0])
                    total_DOS_d_up = [0]*len(matrix[0])
                    total_DOS_f_up = [0]*len(matrix[0])
                    total_DOS_s_down = [0]*len(matrix[0])
                    total_DOS_p_down = [0]*len(matrix[0])
                    total_DOS_d_down = [0]*len(matrix[0])
                    total_DOS_f_down = [0]*len(matrix[0])
                        
                    for j in range(1,len(matrix)):
                        for k in range(1,3,2):
                            total_DOS_s_up = np.array(np.array(total_DOS_s_up) + np.array(matrix[j][k])).tolist()
                        for k in range(2,4,2):
                            total_DOS_s_down = np.array(np.array(total_DOS_s_down) + np.array(matrix[j][k])).tolist()
                        for k in range(3,9,2):
                            total_DOS_p_up = np.array(np.array(total_DOS_p_up) + np.array(matrix[j][k])).tolist()
                        for k in range(4,10,2):
                            total_DOS_p_down = np.array(np.array(total_DOS_p_down) + np.array(matrix[j][k])).tolist()
                        for k in range(9,19,2):
                            total_DOS_d_up = np.array(np.array(total_DOS_d_up) + np.array(matrix[j][k])).tolist()
                        for k in range(10,20,2):
                            total_DOS_d_down = np.array(np.array(total_DOS_d_down) + np.array(matrix[j][k])).tolist()
                        for k in range(19,33,2):
                            total_DOS_f_up = np.array(np.array(total_DOS_f_up) + np.array(matrix[j][k])).tolist()
                        for k in range(20,34,2):
                            total_DOS_f_down = np.array(np.array(total_DOS_f_down) + np.array(matrix[j][k])).tolist()
                       
                        
                    total_DOS_up = np.array(np.array(total_DOS_s_up) + np.array(total_DOS_p_up) + np.array(total_DOS_d_up) + np.array(total_DOS_f_up)).tolist()
                    total_DOS_down = np.array(np.array(total_DOS_s_down) + np.array(total_DOS_p_down) + np.array(total_DOS_d_down) + np.array(total_DOS_f_down)).tolist()
                    total_DOS = np.array(np.array(total_DOS_up) - np.array(total_DOS_down)).tolist()
                    total_DOS_s = np.array(np.array(total_DOS_s_up) - np.array(total_DOS_s_down)).tolist()
                    total_DOS_p = np.array(np.array(total_DOS_p_up) - np.array(total_DOS_p_down)).tolist()
                    total_DOS_d = np.array(np.array(total_DOS_d_up) - np.array(total_DOS_d_down)).tolist()
                    total_DOS_f = np.array(np.array(total_DOS_f_up) - np.array(total_DOS_f_down)).tolist()
                        
                    #print total_DOS_s_up, total_DOS_f_up, total_DOS_p_up
                    int_total_DOS_up = float(np.trapz(total_DOS_up,matrix[1][0]))
                    int_total_DOS_up_fm = float(np.trapz(total_DOS_up[:nl_f],matrix[1][0][:nl_f]))
                    int_total_DOS_down = float(np.trapz(total_DOS_down,matrix[1][0]))
                    int_total_DOS_down_fm = float(np.trapz(total_DOS_down[:nl_f],matrix[1][0][:nl_f]))
                    int_total_DOS = float(np.trapz(total_DOS,matrix[1][0]))
                    int_total_DOS_fm = float(np.trapz(total_DOS[:nl_f],matrix[1][0][:nl_f]))
                    int_total_DOS_s_up = float(np.trapz(total_DOS_s_up,matrix[1][0]))
                    int_total_DOS_s_up_fm = float(np.trapz(total_DOS_s_up[:nl_f],matrix[1][0][:nl_f]))
                    int_total_DOS_s_down = float(np.trapz(total_DOS_s_down,matrix[1][0]))
                    int_total_DOS_s_down_fm = float(np.trapz(total_DOS_s_down[:nl_f],matrix[1][0][:nl_f]))
                    int_total_DOS_p_up = float(np.trapz(total_DOS_p_up,matrix[1][0]))
                    int_total_DOS_p_up_fm = float(np.trapz(total_DOS_p_up[:nl_f],matrix[1][0][:nl_f]))
                    int_total_DOS_p_down = float(np.trapz(total_DOS_p_down,matrix[1][0]))
                    int_total_DOS_p_down_fm = float(np.trapz(total_DOS_p_down[:nl_f],matrix[1][0][:nl_f]))
                    int_total_DOS_d_up = float(np.trapz(total_DOS_d_up,matrix[1][0]))
                    int_total_DOS_d_up_fm = float(np.trapz(total_DOS_d_up[:nl_f],matrix[1][0][:nl_f]))
                    int_total_DOS_d_down = float(np.trapz(total_DOS_d_down,matrix[1][0]))
                    int_total_DOS_d_down_fm = float(np.trapz(total_DOS_d_down[:nl_f],matrix[1][0][:nl_f]))
                    int_total_DOS_f_up = float(np.trapz(total_DOS_f_up,matrix[1][0]))
                    int_total_DOS_f_up_fm = float(np.trapz(total_DOS_f_up[:nl_f],matrix[1][0][:nl_f]))
                    int_total_DOS_f_down = float(np.trapz(total_DOS_f_down,matrix[1][0]))
                    int_total_DOS_f_down_fm = float(np.trapz(total_DOS_f_down[:nl_f],matrix[1][0][:nl_f]))
                    int_total_DOS_s = float(np.trapz(total_DOS_s,matrix[1][0]))
                    int_total_DOS_s_fm = float(np.trapz(total_DOS_s[:nl_f],matrix[1][0][:nl_f]))
                    int_total_DOS_p = float(np.trapz(total_DOS_p,matrix[1][0]))
                    int_total_DOS_p_fm = float(np.trapz(total_DOS_p[:nl_f],matrix[1][0][:nl_f]))
                    int_total_DOS_d = float(np.trapz(total_DOS_d,matrix[1][0]))
                    int_total_DOS_d_fm = float(np.trapz(total_DOS_d[:nl_f],matrix[1][0][:nl_f]))
                    int_total_DOS_f = float(np.trapz(total_DOS_f,matrix[1][0]))
                    int_total_DOS_f_fm = float(np.trapz(total_DOS_f[:nl_f],matrix[1][0][:nl_f]))
                        
                    fnq = raw_input("Enter a name for the output file.\n")
                    print "************************************************************************************************************"
                    file1 = fnq +"_non_spin_polarized"
                    file2 = fnq +"_spin_polarized"
                    f1 = open(file1,"w")
                    f1.write("E_f_Norm" "\t" "Total_s" "\t\t" "Total_p" "\t\t" "Total_d" "\t\t" "Total_f" "\t\t" "Total DOS" "\n")
                    for line in zip(np.array(matrix[1][0]),np.array(total_DOS_s),np.array(total_DOS_p),np.array(total_DOS_d),np.array(total_DOS_f),np.array(total_DOS)):
                        f1.write("\t".join(str("%.8f" %x) for x in line) +"\n")
                    f1.write("Total integral s = "+ str(int_total_DOS_s) + "\n")
                    f1.write("Total integral p = "+ str(int_total_DOS_p) + "\n")
                    f1.write("Total integral d = "+ str(int_total_DOS_d) + "\n")
                    f1.write("Total integral f = "+ str(int_total_DOS_f) + "\n")
                    f1.write("Total integral DOS = "+ str(int_total_DOS) + "\n")
                    f1.write("Total integral s below fermi level = "+ str(int_total_DOS_s_fm) + "\n")
                    f1.write("Total integral p below fermi level = "+ str(int_total_DOS_p_fm) + "\n")
                    f1.write("Total integral d below fermi level = "+ str(int_total_DOS_d_fm) + "\n")
                    f1.write("Total integral f below fermi level = "+ str(int_total_DOS_f_fm) + "\n")
                    f1.write("Total integral DOS below fermi level = "+ str(int_total_DOS_fm) + "\n")
                    f1.close()
                        
                        
                    f2 = open(file2,"w")
                    f2.write("E_f_Norm" "\t" "Total_s_up" "\t" "Total_s_down" "\t" "Total_p_up" "\t" "Total_p_down" "\t" "Total_d_up" "\t" "Total_d_down" "\t" "Total_f_up" "\t" "Total_f_down" "\t" "Total_up" "\t" "Total_down" "\n" )
                    for line in zip(np.array(matrix[1][0]),np.array(total_DOS_s_up),np.array(total_DOS_s_down),np.array(total_DOS_p_up),np.array(total_DOS_p_down),np.array(total_DOS_d_up),np.array(total_DOS_d_down),np.array(total_DOS_f_up),np.array(total_DOS_f_down),np.array(total_DOS_up),np.array(total_DOS_down)):
                        f2.write("\t".join(str("%.8f" %x) for x in line) +"\n")
                    f2.write("Total integral s_up = "+ str(int_total_DOS_s_up) + "\n")
                    f2.write("Total integral s_down = "+ str(int_total_DOS_s_down) + "\n")
                    f2.write("Total integral p_up = "+ str(int_total_DOS_p_up) + "\n")
                    f2.write("Total integral p_down = "+ str(int_total_DOS_p_down) + "\n")
                    f2.write("Total integral d_up = "+ str(int_total_DOS_d_up) + "\n")
                    f2.write("Total integral d_down = "+ str(int_total_DOS_d_down) + "\n")
                    f2.write("Total integral f_up = "+ str(int_total_DOS_f_up) + "\n")
                    f2.write("Total integral f_down = "+ str(int_total_DOS_f_down) + "\n")
                    f2.write("Total integral DOS_up = "+ str(int_total_DOS_up) + "\n")
                    f2.write("Total integral DOS_down = "+ str(int_total_DOS_down) + "\n")
                    f2.write("Total integral s_up below fermi level = "+ str(int_total_DOS_s_up_fm) + "\n")
                    f2.write("Total integral s_down below fermi level = "+ str(int_total_DOS_s_down_fm) + "\n")
                    f2.write("Total integral p_up below fermi level = "+ str(int_total_DOS_p_up_fm) + "\n")
                    f2.write("Total integral p_down below fermi level = "+ str(int_total_DOS_p_down_fm) + "\n")
                    f2.write("Total integral d_up below fermi level = "+ str(int_total_DOS_d_up_fm) + "\n")
                    f2.write("Total integral d_down below fermi level = "+ str(int_total_DOS_d_down_fm) + "\n")
                    f2.write("Total integral f_up below fermi level = "+ str(int_total_DOS_f_up_fm) + "\n")
                    f2.write("Total integral f_down below fermi level = "+ str(int_total_DOS_f_down_fm) + "\n")
                    f2.write("Total integral DOS_up below fermi level = "+ str(int_total_DOS_up_fm) + "\n")
                    f2.write("Total integral DOS_down below fermi level = "+ str(int_total_DOS_down_fm) + "\n")
                    f2.close()
                        
                        
                    plot_multiple(total_DOS_s,total_DOS_s_up,total_DOS_s_down,total_DOS_p,total_DOS_p_up,total_DOS_p_down,total_DOS_d,total_DOS_d_up,total_DOS_d_down,total_DOS_f,total_DOS_f_up,total_DOS_f_down,total_DOS,total_DOS_up,total_DOS_down)
                    #filenames = ["s","p","d","f","total"]
                    #fol_name= raw_input("Enter a name for saving the plots for this calculation\n")
                    #os.mkdir(cwd+"/"+fol_name)
                    #os.chdir(cwd+"/"+fol_name)
                    #for pl_cv in range((len(plot))/3):
                    #    plt.figure(pl_cv)
                    #    plt.plot(matrix[1][0],plot[3*pl_cv],"r-",matrix[1][0],plot[3*pl_cv+1],"g-",matrix[1][0],plot[3*pl_cv+2],"b-")
                    #    plt.title("Density of states vs Fermi normalized energy")
                        #    plt.xlabel("Energy (eV)")
                        #    plt.ylabel("No. of states (a.u.)")
                        #    plt.xlim(-1,1)
                        #    plt.ylim(-0.2,0.2)
                        #    plt.axhline(y=0, color = 'k')
                        #    plt.legend((filenames[pl_cv],filenames[pl_cv]+"_up",filenames[pl_cv]+"_down"))
                        #    plt.savefig(filenames[pl_cv]+".pdf")
                        #    pl_cv +=1
                        #    
                        #os.chdir(cwd)    
                        
                        
                        
                    #print "You still have the same selection of atoms as you had in the previous step."



                
                        #int_total_DOS_up_fermi = float(np.trapz(total_DOS_up_fermi,energy_fermi_fermi))
                        #int_total_DOS_down_fermi = float(np.trapz(total_DOS_down_fermi,energy_fermi_fermi))
                  #print "Total integral spin up DOS ", int_total_DOS_up
                  #print "Total integral spin up DOS below fermi level ",int_total_DOS_up_fermi
                  ##print "Total integral spin down DOS ",int_total_DOS_down
                  #print "Total integral spin down DOS below fermi level  ",int_total_DOS_down_fermi

                elif query2 == "S":
                    orbital_matrix = []
                    while True:
                        print "************************************************************************************************************"
                        orbital = raw_input("Please enter the orbitals that you wish to add for all atoms;\nYou can refer to the README file for orbital no.s or type orbital--h to see the assignment\nYou need to enter all the orbitals that you wish to add in a single line each followed by a space.\nFor eg., if you need to add s_up, s_down and d_up: You need to enter:1 2 5\nEnter the orbital no.s;\n")
                        #print 'if you would like to see the orbital assignments please type orbital'
                        print "************************************************************************************************************"
                        good = True            
                        if orbital == "orbital--h":
                            orbital_assignment_11()
                            good = False
                            
                        else:
                            try:
                                orbital_matrix = [int(orbital.rstrip().split()[i]) for i in range(len(orbital.strip().split(' ')))] 
                            
                            except:
                                ValueError
                                good = False
                                print "You need to input digits for the orbitals. You can type orbital--h at the next prompt to see the orbital assignment."
                        
                            for val in orbital_matrix:
                                if val > 32 or val < 1:
                                    good = False
                                    print 'You have specified an orbital that does not exist. Please re-enter your selection'
                        #print orbital_matrix
                        if good:    
                            break
                                                    
                        #if orbital.isdigit() and int(orbital) == 0:
                         #   break
#                            elif orbital.isdigit() and int(orbital) <33:
#                                orbital_matrix += [int(orbital)]
#                                print orbital_matrix
#                            else:
#                                print "The orbital no. needs to be a number from 1 to 32."
                    print orbital_matrix            
                    total_DOS_des = [0]*len(matrix[0])
                    for j in range(1,len(matrix)):
                        for k in orbital_matrix:
                            if k%2 == 1:
                                total_DOS_des = np.array(np.array(total_DOS_des) + np.array(matrix[j][k])).tolist()
                            else:
                                total_DOS_des = np.array(np.array(total_DOS_des) - np.array(matrix[j][k])).tolist()
                        
                    int_total_DOS_des = float(np.trapz(total_DOS_des,matrix[1][0]))
                    int_total_DOS_des_fm = float(np.trapz(total_DOS_des[:nl_f],matrix[1][0][:nl_f]))

                    fnq2 = raw_input("Enter a name for the output file.\n")
                    print "************************************************************************************************************"
                        
                    f3 = open(fnq2,"w")
                    f3.write("E_f_Norm" "\t\t" "Total_DOS_desired" "\n")
                    for line in zip(np.array(matrix[1][0]),np.array(total_DOS_des)):
                        f3.write("\t\t".join(str(x) for x in line) +"\n")
                    f3.write("Total_integral DOS_desired = " + str(int_total_DOS_des) +"\n")
                    f3.write("Total_integral DOS_desired below fermi level = " + str(int_total_DOS_des_fm) +"\n")
                    f3.close()
                        
                    plot_single(total_DOS_des)
                        
                    #print "You still have the same selection of atoms as you had in the previous step."
                        
                        
                        
                elif query2 == "D":
                    break
                else:
                    print "Please chose the correct letter"
            
        elif  query1 == "S":
            total_DOS_des = [0]*len(Atom_matrix[0])
            while True:                    
                    #print matrix
                posinfo()
                atom_no = raw_input("Please enter the atom no. or 0 if you are done adding the atoms.\n")
                print "************************************************************************************************************"
                if atom_no.isdigit():
                    atom_no = int(atom_no)
                    if atom_no == 0:
                        break
                    elif atom_no >len(Atom_matrix) -1:
                        print "This atom doesn't exit in your system."
                    else:
                        query3 = raw_input("Do you wish to add all orbitals of this atom?\nEnter A for all;\nEnter S for specific;\nAny other letter/number if you are done with this atom.\n")
                        print "************************************************************************************************************"
                        if query3 == "A":
                            for i in range(1,33):
                                if i%2 == 1:
                                    total_DOS_des = np.array(np.array(total_DOS_des) + np.array(Atom_matrix[atom_no][i])).tolist()
                                else:
                                    total_DOS_des = np.array(np.array(total_DOS_des) - np.array(Atom_matrix[atom_no][i])).tolist()
                        elif query3 == "S":
                            orbital_matrix = []
                            while True:
                                print "************************************************************************************************************"
                                orbital = raw_input("Please enter the orbitals that you wish to add for all atoms;\nYou can refer to the README file for orbital no.s or type orbital--h to see the assignment\nYou need to enter all the orbitals that you wish to add in a single line each followed by a space.\nFor eg., if you need to add s_up, s_down and d_up: You need to enter:1 2 5\nEnter the orbital no.s;\n")
                                #print 'if you would like to see the orbital assignments please type orbital'
                                print "************************************************************************************************************"
                                good = True
                                if orbital == "orbital--h":
                                    orbital_assignment_11()
                                    good = False
                                    
                                else:
                                    try:
                                        orbital_matrix = [int(orbital.rstrip().split()[i]) for i in range(len(orbital.strip().split(' ')))]
                                        
                                    except:
                                        ValueError
                                        good = False
                                        print "You need to input digits for the orbitals. You can type orbital--h at the next prompt to see the orbital assignment."
                                                                                    
                                    for val in orbital_matrix:
                                        if val > 32 or val < 1:
                                            good = False
                                            print 'You have specified an orbital that does not exist. Please re-enter your selection'
                                                
                                                
                                                #print orbital_matrix
                                if good:
                                    break
                                                           
                            for i in orbital_matrix:
                                if i%2 == 1:
                                    total_DOS_des = np.array(np.array(total_DOS_des) + np.array(Atom_matrix[atom_no][i])).tolist()
                                else:
                                    total_DOS_des = np.array(np.array(total_DOS_des) - np.array(Atom_matrix[atom_no][i])).tolist()
                                                    
                        else:
                            print "Please choose the correct letter."
                                
                else:
                    print "The atom no. has to be a number less then or equal to the number of atoms in your system."
                
            int_total_DOS_des = float(np.trapz(total_DOS_des,Atom_matrix[1][0]))                 
            int_total_DOS_des_fm = float(np.trapz(total_DOS_des[:nl_f],Atom_matrix[1][0][:nl_f]))
            fnq3 = raw_input("Enter a name for the output file.\n")
            print "************************************************************************************************************"
            f4 = open(fnq3,"w")
            f4.write("E_f_Norm" "\t\t" "Total_DOS_desired" "\n")
            for line in zip(np.array(Atom_matrix[1][0]),np.array(total_DOS_des)):
                f4.write("\t\t".join(str(x) for x in line) +"\n")
            f4.write("Total_integral DOS_desired = " + str(int_total_DOS_des) +"\n")
            f4.write("Total_integral DOS_desired below fermi level = " + str(int_total_DOS_des_fm) +"\n")
            f4.close()
                
            plot_single(total_DOS_des)
                
                
                #fol_name= raw_input("Enter a name for the folder for saving the plots for this calculation\n")
                #os.mkdir(cwd+"/"+fol_name)
                #os.chdir(cwd+"/"+fol_name)
                #plt.figure(fol_name)
                #plt.plot(matrix[1][0],total_DOS_des,"r-")
                #plt.title("Density of states vs Fermi normalized energy")
                #plt.xlabel("Energy (eV)")
                #plt.ylabel("No. of states (a.u.)")
                #plt.xlim(-1,1)
                #plt.ylim(-0.2,0.2)
                #plt.axhline(y=0, color = 'k')
                #plt.savefig("total_DOS_des"+".pdf")
                                
                #os.chdir(cwd)    






                    
        elif query1 == "D":
            print "Thank you for using DOS_analysis.\nThe program will now exit."
            break
            
        else:
            print "Please choose a correct letter."
                
                
                
elif system == "101" or system == "102":
    print "Your system is being saved as a matrix, this may take a moment...or 2"
    print "####################################################################################################\n####################################################################################################"

    def atom_101_102():           # This function creates a matrix atom =[0,Atom 1, Atom 2,.......], where each atom is a matrix with [0, sup, sdown,p_zup.....]
        doscar = open("DOSCAR","r")
        i = 1       #Next 4 lines are to skip the general inforamtion and arrive at the line that contains parameters 
        while i<6:
            temp = doscar.readline()
            #print temp
            i += 1
            
                
        parameters  = doscar.readline().split()     #Read the line that contains paremters
        E_max       = float(parameters[0])          #Maximum energy in the system
        E_min       = float(parameters[1])          #Minimum energy in the system
        nl          = float(parameters[2])          #No. of energy levels in the system
        E_fermi     = float(parameters[3])          #Fermi energy of the system
            #print E_max, E_min, nl, E_fermi
        
        
        NOA = 0     #Next 8 lines count the number of atoms in the system
        while True:
            temp = doscar.readline()
            if temp == "":
                break
            elif len(temp.split()) == 5 and float(temp.split()[2]) == nl and float(temp.split()[3]) == E_fermi:
                NOA = NOA +1
            else:
                continue
            #print NOA
        doscar.close()
        

        doscar = open("DOSCAR","r")     #In the next 26 lines, we create an array for energy and Fermi level normalized energy. We also figure out the no of lines below fermilevel 
        i = 1
        E = []      #Energy list
        E_f = []    #Fermi Level normalized energy list
        while i < 7:
            temp = doscar.readline()
            #print temp
            i += 1
        i = 1
        while i <= nl:
            temp = doscar.readline()
            temp_E = [float(temp.split()[0])]
            temp_E_f = [float(temp.split()[0])-E_fermi]
            E += temp_E
            E_f += temp_E_f
            i += 1
            #print E
            
            
        atom = [E]      #In the next 90 lines, we try to create the atom matrix with first element as a list of energy,next elements are Atom1, Atom2 and so on. Each Atom is a list with elements as list E_f, sup, sdown, and so on.
        #print atom
        i = 1
        while i <= NOA:
            A = [E_f]
            temp = doscar.readline()
            sup           = []
            sdown         = []
            pup           = []
            pdown         = []
            dup           = []
            ddown         = []
            fup           = []
            fdown         = []
                                
            j = 1
            if system == "101":
                while j <= nl:
                    temp = doscar.readline().split()
                    a = 5-len(temp)
                    for k in range(a):
                        temp = temp + ["0"]
                            
                    sup         = sup + [0.5*float(temp[1])]
                    sdown       = sdown + [-0.5*float(temp[1])]
                    pup         = pup + [0.5*float(temp[2])]
                    pdown       = pdown + [-0.5*float(temp[2])]
                    dup         = dup + [0.5*float(temp[3])]
                    ddown       = ddown + [-0.5*float(temp[3])]
                    fup         = fup + [0.5*float(temp[4])]
                    fdown       = fdown + [-0.5*float(temp[4])]
                    j += 1
                   
            else:
                while j <= nl:
                    temp = doscar.readline().split()
                    a = 9-len(temp)
                    for k in range(a):
                        temp = temp + ["0"]
                            
                    sup         = sup + [0.5*float(temp[1])]
                    sdown       = sdown + [-0.5*float(temp[2])]
                    pup         = pup + [0.5*float(temp[3])]
                    pdown       = pdown + [-0.5*float(temp[4])]
                    dup         = dup + [0.5*float(temp[5])]
                    ddown       = ddown + [-0.5*float(temp[6])]
                    fup         = fup + [0.5*float(temp[7])]
                    fdown       = fdown + [-0.5*float(temp[8])]
                    j +=1

                      
              
            A = A + [sup] +[sdown] + [pup]+ [pdown] + [dup] + [ddown] + [fup] + [fdown]
            atom = atom + [A]
            i += 1
            
        doscar.close()
        return atom

        
        
        
        
    def nl_fermi_10():
        nl_f = 0
        i = 0
        while i < len(Atom_matrix[1][0]):
            if Atom_matrix[1][0][i]<=0:
              nl_f += 1
              i += 1
                  #print nl_f
            else:
                break
        return nl_f
            
              
    Atom_matrix = atom_101_102()
    nl_f = nl_fermi_10()


    print "Your system has been saved as a matrix. Let's do some analysis."
    print "The DOSCAR has been split into individual atoms and saved in the folder named split_DOSCAR in your current working directory.  You can refer to split_DOSCAR for any additional analysis.\n"
    print "####################################################################################################\n####################################################################################################"
        #print Atom_matrix
        #print nl_f
        #print Atom_matrix[0]
        
        
        #cwd = os.getcwd()
        #nwd_name = raw_input("What would you like to name the new directory where the data will be saved?\n")
        #nwd = cwd + "/" + nwd_name
        #os.mkdir(nwd)
        #os.chdir(nwd)
    while True:
        query1 = raw_input("Do you wish to study properties for all the atoms in your system or only specific atoms in a particluar layer/adorbent etc.?\nEnter A for all atoms;\nEnter R for a range of atoms;\nEnter S for specific atoms;\nEnter D if you are done with your analysis.\n")
        print "************************************************************************************************************"
        if query1 == "A" or query1 == "R":
            if query1 == "A":
                matrix = Atom_matrix
            else:
                while True:
                    lowlim = raw_input("Please enter the lower limit of the range.\n")
                    highlim = raw_input("Please enter the upper limit of the range.\n")
                    print "************************************************************************************************************"
                        
                    if lowlim.isdigit() and highlim.isdigit():
                        low = int(lowlim)
                        high = int(highlim)
                        if 0<low<=high<len(Atom_matrix):
                                #print Atom_matrix[0]
                            matrix = [Atom_matrix[0]] + Atom_matrix[low:high+1]
                            break
                        else:
                            print "Make sure that the lower limit is less than the upper limit and both of them are less than or equal to the number of atoms in your system."
                                
                    else:
                        print "Both the limits need to be integers less than the number of atoms in the system with lower limit less than or equal to the upper limit."
                    #print matrix[1][1]
            while True:
                query2 = raw_input("Do you need total DOS for all orbitals?\nEnter A for all;\nEnter S for specific;\nEnter D if you are done.\n")
                print "************************************************************************************************************"
                if query2 == "A":
                    total_DOS_s_up = [0]*len(matrix[0])
                    total_DOS_p_up = [0]*len(matrix[0])
                    total_DOS_d_up = [0]*len(matrix[0])
                    total_DOS_f_up = [0]*len(matrix[0])
                    total_DOS_s_down = [0]*len(matrix[0])
                    total_DOS_p_down = [0]*len(matrix[0])
                    total_DOS_d_down = [0]*len(matrix[0])
                    total_DOS_f_down = [0]*len(matrix[0])
                        
                    for j in range(1,len(matrix)):
                        total_DOS_s_up = np.array(np.array(total_DOS_s_up) + np.array(matrix[j][1])).tolist()
                        total_DOS_s_down = np.array(np.array(total_DOS_s_down) + np.array(matrix[j][2])).tolist()
                        total_DOS_p_up = np.array(np.array(total_DOS_p_up) + np.array(matrix[j][3])).tolist()
                        total_DOS_p_down = np.array(np.array(total_DOS_p_down) + np.array(matrix[j][4])).tolist()
                        total_DOS_d_up = np.array(np.array(total_DOS_d_up) + np.array(matrix[j][5])).tolist()
                        total_DOS_d_down = np.array(np.array(total_DOS_d_down) + np.array(matrix[j][6])).tolist()
                        total_DOS_f_up = np.array(np.array(total_DOS_f_up) + np.array(matrix[j][7])).tolist()
                        total_DOS_f_down = np.array(np.array(total_DOS_f_down) + np.array(matrix[j][8])).tolist()
                       
                        
                    total_DOS_up = np.array(np.array(total_DOS_s_up) + np.array(total_DOS_p_up) + np.array(total_DOS_d_up) + np.array(total_DOS_f_up)).tolist()
                    total_DOS_down = np.array(np.array(total_DOS_s_down) + np.array(total_DOS_p_down) + np.array(total_DOS_d_down) + np.array(total_DOS_f_down)).tolist()
                    total_DOS = np.array(np.array(total_DOS_up) - np.array(total_DOS_down)).tolist()
                    total_DOS_s = np.array(np.array(total_DOS_s_up) - np.array(total_DOS_s_down)).tolist()
                    total_DOS_p = np.array(np.array(total_DOS_p_up) - np.array(total_DOS_p_down)).tolist()
                    total_DOS_d = np.array(np.array(total_DOS_d_up) - np.array(total_DOS_d_down)).tolist()
                    total_DOS_f = np.array(np.array(total_DOS_f_up) - np.array(total_DOS_f_down)).tolist()
                        
                        #print total_DOS_s_up, total_DOS_f_up, total_DOS_p_up
                    int_total_DOS_up = float(np.trapz(total_DOS_up,matrix[1][0]))
                    int_total_DOS_up_fm = float(np.trapz(total_DOS_up[:nl_f],matrix[1][0][:nl_f]))
                    int_total_DOS_down = float(np.trapz(total_DOS_down,matrix[1][0]))
                    int_total_DOS_down_fm = float(np.trapz(total_DOS_down[:nl_f],matrix[1][0][:nl_f]))
                    int_total_DOS = float(np.trapz(total_DOS,matrix[1][0]))
                    int_total_DOS_fm = float(np.trapz(total_DOS[:nl_f],matrix[1][0][:nl_f]))
                    int_total_DOS_s_up = float(np.trapz(total_DOS_s_up,matrix[1][0]))
                    int_total_DOS_s_up_fm = float(np.trapz(total_DOS_s_up[:nl_f],matrix[1][0][:nl_f]))
                    int_total_DOS_s_down = float(np.trapz(total_DOS_s_down,matrix[1][0]))
                    int_total_DOS_s_down_fm = float(np.trapz(total_DOS_s_down[:nl_f],matrix[1][0][:nl_f]))
                    int_total_DOS_p_up = float(np.trapz(total_DOS_p_up,matrix[1][0]))
                    int_total_DOS_p_up_fm = float(np.trapz(total_DOS_p_up[:nl_f],matrix[1][0][:nl_f]))
                    int_total_DOS_p_down = float(np.trapz(total_DOS_p_down,matrix[1][0]))
                    int_total_DOS_p_down_fm = float(np.trapz(total_DOS_p_down[:nl_f],matrix[1][0][:nl_f]))
                    int_total_DOS_d_up = float(np.trapz(total_DOS_d_up,matrix[1][0]))
                    int_total_DOS_d_up_fm = float(np.trapz(total_DOS_d_up[:nl_f],matrix[1][0][:nl_f]))
                    int_total_DOS_d_down = float(np.trapz(total_DOS_d_down,matrix[1][0]))
                    int_total_DOS_d_down_fm = float(np.trapz(total_DOS_d_down[:nl_f],matrix[1][0][:nl_f]))
                    int_total_DOS_f_up = float(np.trapz(total_DOS_f_up,matrix[1][0]))
                    int_total_DOS_f_up_fm = float(np.trapz(total_DOS_f_up[:nl_f],matrix[1][0][:nl_f]))
                    int_total_DOS_f_down = float(np.trapz(total_DOS_f_down,matrix[1][0]))
                    int_total_DOS_f_down_fm = float(np.trapz(total_DOS_f_down[:nl_f],matrix[1][0][:nl_f]))
                    int_total_DOS_s = float(np.trapz(total_DOS_s,matrix[1][0]))
                    int_total_DOS_s_fm = float(np.trapz(total_DOS_s[:nl_f],matrix[1][0][:nl_f]))
                    int_total_DOS_p = float(np.trapz(total_DOS_p,matrix[1][0]))
                    int_total_DOS_p_fm = float(np.trapz(total_DOS_p[:nl_f],matrix[1][0][:nl_f]))
                    int_total_DOS_d = float(np.trapz(total_DOS_d,matrix[1][0]))
                    int_total_DOS_d_fm = float(np.trapz(total_DOS_d[:nl_f],matrix[1][0][:nl_f]))
                    int_total_DOS_f = float(np.trapz(total_DOS_f,matrix[1][0]))
                    int_total_DOS_f_fm = float(np.trapz(total_DOS_f[:nl_f],matrix[1][0][:nl_f]))
                        
                    fnq = raw_input("Enter a name for the output file.\n")
                    print "************************************************************************************************************"
                    file1 = fnq +"_non_spinpolarized"
                    file2 = fnq +"_spinpolarized"
                    f1 = open(file1,"w")
                    f1.write("E_f_Norm" "\t" "Total_s" "\t\t" "Total_p" "\t\t" "Total_d" "\t\t" "Total_f" "\t\t" "Total DOS" "\n")
                    for line in zip(np.array(matrix[1][0]),np.array(total_DOS_s),np.array(total_DOS_p),np.array(total_DOS_d),np.array(total_DOS_f),np.array(total_DOS)):
                        f1.write("\t".join(str("%.8f" %x) for x in line) +"\n")
                    f1.write("Total integral s = "+ str(int_total_DOS_s) + "\n")
                    f1.write("Total integral p = "+ str(int_total_DOS_p) + "\n")
                    f1.write("Total integral d = "+ str(int_total_DOS_d) + "\n")
                    f1.write("Total integral f = "+ str(int_total_DOS_f) + "\n")
                    f1.write("Total integral DOS = "+ str(int_total_DOS) + "\n")
                    f1.write("Total integral s below fermi level = "+ str(int_total_DOS_s_fm) + "\n")
                    f1.write("Total integral p below fermi level = "+ str(int_total_DOS_p_fm) + "\n")
                    f1.write("Total integral d below fermi level = "+ str(int_total_DOS_d_fm) + "\n")
                    f1.write("Total integral f below fermi level = "+ str(int_total_DOS_f_fm) + "\n")
                    f1.write("Total integral DOS below fermi level = "+ str(int_total_DOS_fm) + "\n")
                    f1.close()
                        
                        
                    f2 = open(file2,"w")
                    f2.write("E_f_Norm" "\t" "Total_s_up" "\t" "Total_s_down" "\t" "Total_p_up" "\t" "Total_p_down" "\t" "Total_d_up" "\t" "Total_d_down" "\t" "Total_f_up" "\t" "Total_f_down" "\t" "Total_up" "\t" "Total_down" "\n" )
                    for line in zip(np.array(matrix[1][0]),np.array(total_DOS_s_up),np.array(total_DOS_s_down),np.array(total_DOS_p_up),np.array(total_DOS_p_down),np.array(total_DOS_d_up),np.array(total_DOS_d_down),np.array(total_DOS_f_up),np.array(total_DOS_f_down),np.array(total_DOS_up),np.array(total_DOS_down)):
                        f2.write("\t".join(str("%.8f" %x) for x in line) +"\n")
                    f2.write("Total integral s_up = "+ str(int_total_DOS_s_up) + "\n")
                    f2.write("Total integral s_down = "+ str(int_total_DOS_s_down) + "\n")
                    f2.write("Total integral p_up = "+ str(int_total_DOS_p_up) + "\n")
                    f2.write("Total integral p_down = "+ str(int_total_DOS_p_down) + "\n")
                    f2.write("Total integral d_up = "+ str(int_total_DOS_d_up) + "\n")
                    f2.write("Total integral d_down = "+ str(int_total_DOS_d_down) + "\n")
                    f2.write("Total integral f_up = "+ str(int_total_DOS_f_up) + "\n")
                    f2.write("Total integral f_down = "+ str(int_total_DOS_f_down) + "\n")
                    f2.write("Total integral DOS_up = "+ str(int_total_DOS_up) + "\n")
                    f2.write("Total integral DOS_down = "+ str(int_total_DOS_down) + "\n")
                    f2.write("Total integral s_up below fermi level = "+ str(int_total_DOS_s_up_fm) + "\n")
                    f2.write("Total integral s_down below fermi level = "+ str(int_total_DOS_s_down_fm) + "\n")
                    f2.write("Total integral p_up below fermi level = "+ str(int_total_DOS_p_up_fm) + "\n")
                    f2.write("Total integral p_down below fermi level = "+ str(int_total_DOS_p_down_fm) + "\n")
                    f2.write("Total integral d_up below fermi level = "+ str(int_total_DOS_d_up_fm) + "\n")
                    f2.write("Total integral d_down below fermi level = "+ str(int_total_DOS_d_down_fm) + "\n")
                    f2.write("Total integral f_up below fermi level = "+ str(int_total_DOS_f_up_fm) + "\n")
                    f2.write("Total integral f_down below fermi level = "+ str(int_total_DOS_f_down_fm) + "\n")
                    f2.write("Total integral DOS_up below fermi level = "+ str(int_total_DOS_up_fm) + "\n")
                    f2.write("Total integral DOS_down below fermi level = "+ str(int_total_DOS_down_fm) + "\n")
                    f2.close()
                        
                    plot_multiple(total_DOS_s,total_DOS_s_up,total_DOS_s_down,total_DOS_p,total_DOS_p_up,total_DOS_p_down,total_DOS_d,total_DOS_d_up,total_DOS_d_down,total_DOS_f,total_DOS_f_up,total_DOS_f_down,total_DOS,total_DOS_up,total_DOS_down)
                        
                    #print "You still have the same selection of atoms as you had in previous step."




                        #int_total_DOS_up_fermi = float(np.trapz(total_DOS_up_fermi,energy_fermi_fermi))
                        #int_total_DOS_down_fermi = float(np.trapz(total_DOS_down_fermi,energy_fermi_fermi))
                  #print "Total integral spin up DOS ", int_total_DOS_up
                  #print "Total integral spin up DOS below fermi level ",int_total_DOS_up_fermi
                  ##print "Total integral spin down DOS ",int_total_DOS_down
                  #print "Total integral spin down DOS below fermi level  ",int_total_DOS_down_fermi

                elif query2 == "S":
                    orbital_matrix = []
                    while True:
                        print "************************************************************************************************************"
                        orbital = raw_input("Please enter the orbitals that you wish to add for all atoms;\nYou can refer to the README file for orbital no.s or type orbital--h to see the assignment\nYou need to enter all the orbitals that you wish to add in a single line each followed by a space.\nFor eg., if you need to add s_up, s_down and d_up: You need to enter:1 2 5\nEnter the orbital no.s;\n")
                        #print 'if you would like to see the orbital assignments please type orbital'
                        print "************************************************************************************************************"
                        good = True            
                        if orbital == "orbital--h":
                            orbital_assignment_10()
                            good = False
                            
                        else:
                            try:
                                orbital_matrix = [int(orbital.rstrip().split()[i]) for i in range(len(orbital.strip().split(' ')))] 
                            
                            except:
                                ValueError
                                good = False
                                print "You need to input digits for the orbitals. You can type orbital--h at the next prompt to see the orbital assignment."
                        
                            for val in orbital_matrix:
                                if val > 8 or val < 1:
                                    good = False
                                    print 'You have specified an orbital that does not exist. Please re-enter your selection'
                        #print orbital_matrix
                        if good:    
                            break
                                
                    total_DOS_des = [0]*len(matrix[0])
                    for j in range(1,len(matrix)):
                        for k in orbital_matrix:
                            if k%2 == 1:
                                total_DOS_des = np.array(np.array(total_DOS_des) + np.array(matrix[j][k])).tolist()
                            else:
                                total_DOS_des = np.array(np.array(total_DOS_des) - np.array(matrix[j][k])).tolist()
                        
                    int_total_DOS_des = float(np.trapz(total_DOS_des,matrix[1][0]))
                    int_total_DOS_des_fm = float(np.trapz(total_DOS_des[:nl_f],matrix[1][0][:nl_f]))
                        
                    fnq2 = raw_input("Enter a name for the output file.\n")
                    print "************************************************************************************************************"
                    f3 = open(fnq2,"w")
                    f3.write("E_f_Norm" "\t\t" "Total_DOS_desired" "\n")
                    for line in zip(np.array(matrix[1][0]),np.array(total_DOS_des)):
                        f3.write("\t\t".join(str(x) for x in line) +"\n")
                    f3.write("Total_integral DOS_desired = " + str(int_total_DOS_des) +"\n")
                    f3.write("Total_integral DOS_desired below fermi level = " + str(int_total_DOS_des_fm) +"\n")
                    f3.close()
                        
                    plot_single(total_DOS_des)
                        
                    #print "You still have the same selection of atoms as you had in previous step."
                elif query2 == "D":
                    break
                else:
                    print "Please choose a correct letter"
            
        elif  query1 == "S":
            total_DOS_des = [0]*len(Atom_matrix[0])
            while True:                    
                    #print matrix
                posinfo()
                atom_no = raw_input("Please enter the atom no. or 0 if you are done adding the atoms\n")
                print "************************************************************************************************************"
                if atom_no.isdigit():
                    atom_no = int(atom_no)
                    if atom_no ==0:
                        break
                    elif atom_no > len(Atom_matrix)-1:
                        print "This atom doesn't exist in your system."
                            
                    else:
                        query3 = raw_input("Do you wish to add all orbitals of this atom?\nEnter A for all;\nEnter S for specific.\n")
                        print "************************************************************************************************************"
                        if query3 == "A":
                            for i in range(1,9):
                                if i%2 == 1:
                                    total_DOS_des = np.array(np.array(total_DOS_des) + np.array(Atom_matrix[atom_no][i])).tolist()
                                else:
                                    total_DOS_des = np.array(np.array(total_DOS_des) - np.array(Atom_matrix[atom_no][i])).tolist()
                        elif query3 == "S":
                            orbital_matrix = []
                            while True:
                                print "************************************************************************************************************"
                                orbital = raw_input("Please enter the orbitals that you wish to add for all atoms;\nYou can refer to the README file for orbital no.s or type orbital--h to see the assignment\nYou need to enter all the orbitals that you wish to add in a single line each followed by a space.\nFor eg., if you need to add s_up, s_down and d_up: You need to enter:1 2 5\nEnter the orbital no.s;\n")
                                #print 'if you would like to see the orbital assignments please type orbital'
                                print "************************************************************************************************************"
                                good = True
                                if orbital == "orbital--h":
                                    orbital_assignment_10()
                                    good = False
                                    
                                else:
                                    try:
                                        orbital_matrix = [int(orbital.rstrip().split()[i]) for i in range(len(orbital.strip().split(' ')))]
                                        
                                    except:
                                        ValueError
                                        good = False
                                        print "You need to input digits for the orbitals. You can type orbital--h at the next prompt to see the orbital assignment."
                                                                                    
                                    for val in orbital_matrix:
                                        if val > 8 or val < 1:
                                            good = False
                                            print 'You have specified an orbital that does not exist. Please re-enter your selection'
                                                
                                                
                                                #print orbital_matrix
                                if good:
                                    break
                                                           
                            for i in orbital_matrix:
                                if i%2 == 1:
                                    total_DOS_des = np.array(np.array(total_DOS_des) + np.array(Atom_matrix[atom_no][i])).tolist()
                                else:
                                    total_DOS_des = np.array(np.array(total_DOS_des) - np.array(Atom_matrix[atom_no][i])).tolist()

                                
                            
                        else:
                            print "Please choose the correct letter."
                                
                else:
                    print "The atom no. has to be a number less than or equal to the number of atoms in your system."
                
            int_total_DOS_des = float(np.trapz(total_DOS_des,Atom_matrix[1][0]))                 
            int_total_DOS_des_fm = float(np.trapz(total_DOS_des[:nl_f],Atom_matrix[1][0][:nl_f]))
            fnq3 = raw_input("Enter a name for the output file.\n")
            print "************************************************************************************************************"
            f4 = open(fnq3,"w")
            f4.write("E_f_Norm" "\t\t" "Total_DOS_desired" "\n")
            for line in zip(np.array(Atom_matrix[1][0]),np.array(total_DOS_des)):
                f4.write("\t\t".join(str(x) for x in line) +"\n")
            f4.write("Total_integral DOS_desired = " + str(int_total_DOS_des) +"\n")
            f4.write("Total_integral DOS_desired below fermi level = " + str(int_total_DOS_des_fm) +"\n")
            f4.close()
                
            plot_single(total_DOS_des)





                    
        elif query1 == "D":
            print "Thank you for using DOS_analysis.\nThe program will now exit."
            break
            
        else:
            print "Please choose a correct letter."
                
else:
    print "Thank you for using DOS_analysis.\nThe program will now exit."
        
#else:
 #   print "Thank you for using DOS_analysis.\nThe program will now exit."
    

        
                        
