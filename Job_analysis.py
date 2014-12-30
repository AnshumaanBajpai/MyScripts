#!/afs/crc.nd.edu/x86_64_linux/python/2.7.4/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 15:04:08 2014

@author: abajpai1
"""
##Electronic energy plot, ##magnetic moment#
############################################################################################################################################################################
########################################################## This code is to be used for analyzing a vasp calculation ######################################################## 
############################################################### and figuring out issues with your job, if any.##############################################################
################################################### In case you find this script useful for your work, feel free to take ###################################################
#####################################################  Anshumaan Bajpai and Michael Penninger out for lunch or something. ##################################################
############################################################################################################################################################################

import string, numpy as np, os, matplotlib.pyplot as plt, re
import shutil
from shutil import rmtree
from StringIO import StringIO
from pylab import *
import math
import argparse, linecache
from matplotlib import rcParams

plt.close('all')


############################################################################################################################################################################
######################################################### This section gathers information from the CONTCAR ################################################################
############################################### We are using the CONTCAR to find out which atoms/dimensions have been flagged TRUE #########################################
################################################## as these will be considered when we analyze the forces on atoms #########################################################
############################################################################################################################################################################

def continfo():
    contcar = open("CONTCAR",'r')
    #total_cont = contcar.readlines()
    #print total_cont
    info = []                                   # will contain the list of information about the coordinates fo each atom and its flags for fixed or relaxed 
    i =1
    atom_list = []                              # will contain the list of atoms including the number of their instances
    title = contcar.readline()
    sf = contcar.readline().split()   #scale factor
    #print sf
    lattice = [[float(i) for i in contcar.readline().split()]] + [[float(i) for i in contcar.readline().split()]] + [[float(i) for i in contcar.readline().split()]] # stores all the three unit vectors as a single string
    #print lattice
    atom_type = contcar.readline().split()
    nea = [int(i) for i in contcar.readline().split()]
    for i in range(len(atom_type)):
        atom_list = atom_list + [atom_type[i]]*nea[i]
    #print atom_list
    #print atom
            
    if contcar.readline().split()[0]=="Selective":
        data = np.genfromtxt("CONTCAR", dtype = [('X','f8'),('Y','f8'),('Z','f8'),('X_rel','S1'),('Y_rel','S1'),('Z_rel','S1')], skip_header=9, skip_footer=len(atom_list))
        info = [atom_list] + [np.array(data['X']).tolist()] + [np.array(data['Y']).tolist()] + [np.array(data['Z']).tolist()] + [np.array(data['X_rel']).tolist()] + [np.array(data['Y_rel']).tolist()] + [np.array(data['Z_rel']).tolist()]
    else:
        data = np.genfromtxt("CONTCAR", dtype = [('X','f8'),('Y','f8'),('Z','f8')], skip_header=8, skip_footer=len(atom_list))
        info = [atom_list] + [np.array(data['X']).tolist()] + [np.array(data['Y']).tolist()] + [np.array(data['Z']).tolist()] + [['T']*len(atom_list)] + [['T']*len(atom_list)] + [['T']*len(atom_list)]
    
    contcar.close()
    return [info] + [lattice] +[sf]         # returns a list with three elements: first is hte position and the flags for atoms, second is lattice vectors and third is scale factor
    
atom_des = continfo()[0]
#print atom_des
lattice = continfo()[1]
#print lattice
sf = float(continfo()[2][0])
#print hihi
#print sf
#print lattice
#print atom_des
############################################################################################################################################################################
#################################################### Data from CONTCAR has been gathered ###################################################################################
############################################################################################################################################################################





############################################################################################################################################################################
############################################################ FINDING ISPIN from INCAR ######################################################################################
############################################################################################################################################################################

outcar = open("OUTCAR",'r')
outcar_s = outcar.read()
searchISPIN = re.search(r'ISPIN  = (.* ) (.*) (.*)',outcar_s, re.M|re.I)            # searches through the OUTCAR for ISPIN tag and assigns its value to the variable ISPIN
ISPIN  = int(searchISPIN.group(1))
#print ISPIN
outcar.close()


############################################################################################################################################################################
###################################################################### ISPIN found #########################################################################################
############################################################################################################################################################################

############################################################################################################################################################################
############################################################################### function to do the plotting ################################################################
############################################################################################################################################################################

def plot_variable(A,B,C,D,E,F):
    if not os.path.isdir(cwd+"/job_analysis_plots"):
        os.mkdir(cwd+"/job_analysis_plots")
    os.chdir(cwd+"/job_analysis_plots")
    if os.path.isfile("E_t_F_vs_Geo.pdf"):
        os.remove("E_t_F_vs_Geo.pdf")
    
    f = figure(1)
    ax = f.add_axes([0,0,1,1])
    subplots_adjust(hspace = .2)
    
    ax1 = subplot(311)
    ax1.plot(A,B,'ro--')
    #ax1.plot(A,B,'g--')
    ylabel('Energy (eV)',fontsize = 10, color = 'r')
    ax1.tick_params(axis='y', colors='r')
    #ax1.yaxis.set_label_position("right")
    
    
    ax2 = subplot(312, sharex=ax1)
    ax2.plot(A,C,'ro--')
    ylabel("time/(elec.step)  (s)",fontsize = 10, color = 'r')
    ax2.tick_params(axis='y', colors='r')
    #ax2.yaxis.label.set_color('r')
    #ax2.yaxis.set_label_position("right")
    
    ax3 = subplot(313, sharex=ax1)
    ax3.plot(A,D,'ro--')
    ylabel("Max. force (eV/Angst)", fontsize = 10)
    ylim([0,3])
    xlabel("# of Geometric steps")
    #ax3.yaxis.set_label_position("right")
    ax4 = ax2.twinx()
    ax4.plot(A,E,'g+')
    ylabel("Total time (h)", fontsize = 10, color = 'g')
    ax4.tick_params(axis = 'y', colors = 'g')
    #ax4.yaxis.set_label_position
    
    ax5 = ax1.twinx()
    ax5.plot(A,F,'g*--')
    ylabel("Delta E (eV)", fontsize = 10, color = 'g')
    ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax5.tick_params(axis = 'y', colors = 'g') 
    
    xticklabels = ax1.get_xticklabels() + ax2.get_xticklabels()
    setp(xticklabels, visible = False)
    text(0.5, 0.92, "Parameters vs # of Geometric step",verticalalignment = 'center', horizontalalignment = 'center', transform = ax.transAxes)
    
    savefig("E_t_F_vs_Geo.pdf")
    show()
    os.chdir(cwd)    

############################################################################################################################################################################
################################################################# Sifting OSZICAR for energies: Code obtained from Mike#####################################################
############################################################################################################################################################################

itera=[]  #this sets up a blank array that will allow you write to it or append later
geoenergy=[] #list of energies after each geometric step
magmom=[]    # list of magnetic moments from OSZICAR file: Will  have information only if ISPIN = 1
coleleciter=[]  # this is a temporary list for storing no. of electronic iterations each geometric step. Adds data to eleciter after each geometric step
colelecenergy=[]	
colelecdelta=[]
eleciter=[]
elecenergy=[]
elecdelta=[]
elecdelta_f=[]

####Pull electronic and geometry step energies
####Supoorts both ALGO=Fast and Normal 
####Fast should have 5 DAV steps and the rest RMM
####Normal should have all DAV steps

for line in open('OSZICAR'):
    if "DAV:" in line or "RMM:" in line:
        values = line.split()
        #print values[1]
        coleleciter.append(float(values[1]))  ##append the values to open blank array
        colelecenergy.append(float(values[2]))
        colelecdelta.append(float(values[3]))
        
    elif "F=" in line:  ### look for geometry step
        line=line.split()
        itera.append(float(line[0]))
        geoenergy.append(float(line[4]))
        if ISPIN == 2:
            magmom.append(float(line[-2]))
        else:
            magmom.append(float(0))

        eleciter.append(coleleciter)
        elecenergy.append(colelecenergy)
        elecdelta.append(colelecdelta)
        elecdelta_f.append(colelecdelta[-1:][0])
        coleleciter=[]
        colelecenergy=[]	
        colelecdelta=[]
           
############################################################################################################################################################################
################################################################# Sifting OSZICAR for energies complete#####################################################################
############################################################################################################################################################################

############################################################################################################################################################################
#################################################### Sifting OUTCAR for forces and electronic step time ####################################################################
############################################################################################################################################################################

outcar = open('OUTCAR','r')
iteration = 0
cwd = os.getcwd()
nwd = cwd + "/forces"
if os.path.isdir(nwd):
    shutil.rmtree(nwd)
            
os.mkdir(nwd)
os.chdir(nwd)
iter_list =[]
max_force = []
t_estep = []
while True:
    temp = outcar.readline()
    if temp =='':
        break
    elif len(temp.split()) > 0 and temp.split()[0] == "LOOP:":
        t_estep.append(float(temp.split()[-1]))
    else:
        if len(temp.split())>0 and temp.split()[0] == "POSITION":
            iteration +=1
            iter_list.append(iteration)
            #print iter_list
            skip = outcar.readline()
            force = open("iter_"+str(iteration),'w')
            fti = []                                           # force list this iteration
            for i in range(len(atom_des[0])):
                f = outcar.readline()
                force.write(f)
                for j in (4,5,6):
                    if atom_des[j][i] == "T":
                        fti.append(abs(float(f.split()[j-1])))
                        
            force.close()
            max_force.append(max(fti))
            #max_force = max_force + [max(fti)] 

os.chdir(cwd)
outcar.close()    

############################################################################################################################################################################
#################################################### Sifting OUTCAR for forces and electronic step time complete ###########################################################
############################################################################################################################################################################

############################################################################################################################################################################
#################################################### Calculating average time for each electronic step #####################################################################
############################################################################################################################################################################

t_eleciter_avg = []
t_geom_ex = [0]
#print t_geom_ex[0]
#print hihi
for t in range(len(iter_list)):
    temp_t_iter = t_estep[:len(elecenergy[t])]
    temp_t_geom = sum(temp_t_iter)
    t_geom_ex.append(temp_t_geom + t_geom_ex[t])
    t_avg = sum(temp_t_iter)/len(elecenergy[t])
    t_eleciter_avg.append(t_avg)
    t_estep = t_estep[len(elecenergy[t]):]

t_geom = t_geom_ex[1:]
t_geom = [i/3600 for i in t_geom]
#print len(t_geom)
#print t_geom
#print hihi
############################################################################################################################################################################
#################################################### Calculating average time for each electronic step complete ############################################################
############################################################################################################################################################################
geoenergy = [float("%.3f" %i) for i in geoenergy]
#print len(iter_list)
#print len(geoenergy)
#print len(t_eleciter_avg)
#print len(max_force)

plot_variable(iter_list,geoenergy,t_eleciter_avg,max_force,t_geom,elecdelta_f)

#print len(iter_list)
#print hihi

############################################################################################################################################################################
#################################################### Gathers data from XDATCAR #############################################################################################
############################################################################################################################################################################

def xdatinfo():
    xdatcar = open("XDATCAR",'r')
    #total_cont = contcar.readlines()
    #print total_cont
    x_i = []
    y_i = []
    z_i = []
    x_f = []
    y_f = []
    z_f = []
    while True:        
        temp = xdatcar.readline()
        if temp == "":
            break
        elif temp.split()[0] == "Direct" and temp.split()[-1] == "1":
            for i in range(len(atom_des[0])):
                atom_pos = xdatcar.readline().split()
                x_i.append(float(atom_pos[0]))
                y_i.append(float(atom_pos[1]))
                z_i.append(float(atom_pos[2]))
                
        elif temp.split()[0] == "Direct" and temp.split()[-1] == str(len(iter_list)):
            for i in range(len(atom_des[0])):
                atom_pos = xdatcar.readline().split()
                x_f.append(float(atom_pos[0]))
                y_f.append(float(atom_pos[1]))
                z_f.append(float(atom_pos[2]))
                
        else:
            continue
#    print x_i
#    print x_f
    x_disp = np.array(np.array(x_f)-np.array(x_i)).tolist()
#    x_disp = [abs(i) for i in x_disp]
    y_disp = np.array(np.array(y_f)-np.array(y_i)).tolist()
#    y_disp = [abs(i) for i in y_disp]
    #print y_disp
    z_disp = np.array(np.array(z_f)-np.array(z_i)).tolist()
#    z_disp = [abs(i) for i in z_disp]
    x_mov = []
    y_mov = []
    z_mov = []
    counter = 0
    for j in [x_disp, y_disp, z_disp]:
        counter +=1
        #print 'this is J'
        #print j
        #exit
        #print hihi
        for i in range(len(atom_des[0])):
            if -0.5 <= j[i] <= 0.5:
                if counter ==1:
                    x_mov.append(j[i])
                elif counter == 2:                    
                    y_mov.append(j[i])
                else:
                    z_mov.append(j[i])
                    
            elif j[i] > 0.5:
                if counter == 1:
                    x_mov.append(j[i]-1)
                elif counter == 2:
                    y_mov.append(j[i]-1)
                else:
                    z_mov.append(j[i]-1)
            else:
                if counter == 1:
                    x_mov.append(j[i]+1)
                elif counter == 2:
                    y_mov.append(j[i]+1)
                else:
                    z_mov.append(j[i]+1)
                
#    print counter
#    print counter2        
    xdatcar.close()
#    print y_mov
    return [x_mov]+[y_mov]+[z_mov]

xdata = xdatinfo()

############################################################################################################################################################################
#################################################### Data gathered from XDATCAR ############################################################################################
############################################################################################################################################################################

############################################################################################################################################################################
#################################################### calculating distances along the vectors ###############################################################################
############################################################################################################################################################################
#print lattice
#print type(sf)
lat_x = sf*math.sqrt(math.pow(lattice[0][0], 2) + math.pow(lattice[0][1], 2) + math.pow(lattice[0][2], 2))
#print lat_x
lat_y = sf*math.sqrt(math.pow(lattice[1][0], 2) + math.pow(lattice[1][1], 2) + math.pow(lattice[1][2], 2))
#print lat_y
lat_z = sf*math.sqrt(math.pow(lattice[2][0], 2) + math.pow(lattice[2][1], 2) + math.pow(lattice[2][2], 2))
#print lat_z

#print xdata[1][0]
x_dist = [abs(i*lat_x) for i in xdata[0]]
#print x_dist
 
y_dist = [abs(i*lat_y) for i in xdata[1]]
z_dist = [abs(i*lat_z) for i in xdata[2]]
#print "ads"
############################################################################################################################################################################
#################################################### distances along the vectors ###########################################################################################
############################################################################################################################################################################

############################################################################################################################################################################
#################################################### Calculating absolute distances ########################################################################################
############################################################################################################################################################################


x_abs = np.array(xdata[0])*lattice[0][0]*sf + np.array(xdata[1])*lattice[1][0]*sf + np.array(xdata[2])*lattice[2][0]*sf
y_abs = np.array(xdata[0])*lattice[0][1]*sf + np.array(xdata[1])*lattice[1][1]*sf + np.array(xdata[2])*lattice[2][1]*sf
z_abs = np.array(xdata[0])*lattice[0][2]*sf + np.array(xdata[1])*lattice[1][2]*sf + np.array(xdata[2])*lattice[2][2]*sf
#print x_abs
atom_disp_sq = np.square(x_abs) + np.square(y_abs) + np.square(z_abs)
atom_disp = np.array(np.sqrt(atom_disp_sq)).tolist()
#print atom_disp

############################################################################################################################################################################
#################################################### absolute distances calcualted #########################################################################################
############################################################################################################################################################################

############################################################################################################################################################################
#################################################### Plotting bar charts for distances #####################################################################################
############################################################################################################################################################################
atoms = tuple(atom_des[0])

os.chdir(cwd+"/job_analysis_plots")
if os.path.isfile("displacements.pdf"):
    os.remove("displacements.pdf")
f2 = figure(2)
ax = f2.add_axes([0,0,1,1])
subplots_adjust(hspace = .3)
atom_pos = np.arange(len(atoms))


ax1 = subplot(411)
ax1.bar(atom_pos,x_dist, width = 0.5, color = 'r')
ylabel('LV1_dis (A)',fontsize = 10)
ax1.tick_params(axis='y', labelsize = 10)
legend("X", loc = "best")
    
    
ax2 = subplot(412, sharex=ax1)
ax2.bar(atom_pos,y_dist, width = 0.5, color = 'r')
ylabel("LV2_dis (A)",fontsize = 10)
ax2.tick_params(axis='y', labelsize = 10)
legend("Y",loc = "best")
#ax2.yaxis.set_label_position("right")
    
ax3 = subplot(413, sharex=ax1)
ax3.bar(atom_pos,z_dist, width = 0.5, color = 'r')
ylabel("LV3_dis (A)", fontsize = 10)
ax3.tick_params(axis='y', labelsize = 10)
legend("Z", loc = "best")

ax4 = subplot(414, sharex=ax1)
ax4.bar(atom_pos, atom_disp, width = 0.5, color = 'r')
ylabel("Displacement (A)", fontsize = 10)
xlabel("Atoms")
#ax3.yaxis.set_label_position("right")
    
#xticklabels = ax1.get_xticklabels() + ax2.get_xticklabels()
#setp(xticklabels, visible = False)
xticks(atom_pos,atoms)
text(0.5, 0.92, "Atomic displacements from intial to final",verticalalignment = 'center', horizontalalignment = 'center', transform = ax.transAxes)
text(0.1, 0.01, "LV = Lattice Vector", transform = ax.transAxes)

savefig("displacements.pdf")    
show()
os.chdir(cwd)    
############################################################################################################################################################################
#################################################### bar charts for distances complete #####################################################################################
############################################################################################################################################################################
#print magmom

if ISPIN == 2:
    os.chdir(cwd+"/job_analysis_plots")
    if os.path.isfile("magnetic_moment.pdf"):
        os.remove("magnetic_moment.pdf")

    f3 = figure(3)
    ax = f3.add_axes([0,0,1,1])
    ax1 = subplot(111)
    ax1.plot(iter_list,magmom)
    xlabel("# of Geometric steps")
    ylabel("Magnetic moment")
    text(0.5,0.92, "Magnetic moment vs geometric step",verticalalignment = 'center', horizontalalignment = 'center', transform = ax.transAxes)
    plt.savefig("magnetic_moment.pdf")
    show()
    os.chdir(cwd)
    
