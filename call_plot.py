import numpy as np
import matplotlib.pyplot as plt
from os import listdir
import sys

data_path="./data/"
plot_path="./plots/"
allfiles=listdir(data_path)
allfiles_sorted=sorted(allfiles)

nx=int(sys.argv[1])
ny=int(sys.argv[2])
plotvar=int(sys.argv[3])
ICtype=sys.argv[4]
if plotvar<0 or plotvar>7:
    print("Error as wrong variable type selected!")
    print("Default variable of density taken instead.")
    plotvar=0
    # sys.exit()
print("Begin plotting on Python script: nx="+str(nx)+" ny="+str(ny))

def choose_plotvar(x):
    return {
            0:'Density',
            1:'U_vel',
            2:'V_vel',
            3:'W_vel',
            4:'Pressure',
            5:'Bx',
            6:'By',
            7:'Bz'
            }.get(x,'Unknown')

pvar=choose_plotvar(plotvar)

frame_no=0
for fil_name in allfiles_sorted:
    fname=data_path+fil_name
    data=np.loadtxt(fname)
    Z=np.reshape(data[:,plotvar+2],(ny,nx))
    X=np.reshape(data[:,0],(ny,nx))
    Y=np.reshape(data[:,1],(ny,nx))
    plt.title(pvar,loc='left')
    plt.title(fil_name)
    plt.xlabel('x')
    plt.ylabel('y')
    CS=plt.contour(X,Y,Z,30,cmap=plt.cm.gray)
    #CS=plt.contour(X,Y,Z,10,colors='k')
    #plt.clabel(CS,inline=True,fontsize=10)
    plt.colorbar()
    plt.savefig(plot_path + ICtype + '_' + pvar + str(frame_no).zfill(3)+'.png')
    plt.clf()
    frame_no=frame_no+1
    if (frame_no % 10==0):
        sys.stdout.write('\rPlotting file at '+str(frame_no)+'/'+str( len(allfiles)))
        sys.stdout.flush()

print("")
print("Completed plots!")

