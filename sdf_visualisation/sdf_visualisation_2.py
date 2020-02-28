import sdf
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

#################################################################
# USER DATA
files_dir = "../rtfiles/pau-search/L1/"
prefix = "AdS5xS5_L1_"
var_name = "gb_tx"
time_level = "2" # default: time level n (1 -> nm1, 2 -> n, 3 -> np1)
num_procs = 6
xmin = 0.8
xmax = 1.0

#################################################################
# Loading dataset
files = np.empty(num_procs, dtype=object)
for i in range(0,num_procs):
    files[i] = files_dir + prefix + var_name + "_tl" + time_level + "_" + str(i) + ".sdf"
    
data2d = sdf.loadSDFFiles(files)
data2d_time=sdf.groupByTime(data2d)
    
#################################################################
# Make the plot
tt = 2 #time slice
print("Plotting time = %f" % data2d_time[0][tt][0].time)

x1 = list(range(0,len(data2d_time[0][tt])))
x2 = list(range(0,len(data2d_time[0][tt])))

#plt.ion()
fig = plt.figure()
ax = Axes3D(fig)
#ax = fig.add_subplot(111, projection='3d')

#i = 1
#print("Showing proc = %d" % i)

#x1[i] = np.linspace(data2d_time[0][tt][i].bbox[0], data2d_time[0][tt][i].bbox[1],
#                    num=data2d_time[0][tt][i].shape[0])
#x2[i] = np.linspace(data2d_time[0][tt][i].bbox[2], data2d_time[0][tt][i].bbox[3],
#                    num=data2d_time[0][tt][i].shape[1])
#X,Y = np.meshgrid(x1[i],x2[i])
#ax.plot_surface(X,Y,data2d_time[0][tt][i]._data.T,
#                rstride=1,cstride=1,color='c', alpha=0.6)
                
#print(data2d_time[0][tt][i]._data.T)

for i in range(0,len(data2d_time[0][tt])):
    x1[i] = np.linspace(data2d_time[0][tt][i].bbox[0], data2d_time[0][tt][i].bbox[1],
                        num=data2d_time[0][tt][i].shape[0])
    x2[i] = np.linspace(data2d_time[0][tt][i].bbox[2], data2d_time[0][tt][i].bbox[3],
                        num=data2d_time[0][tt][i].shape[1])
    X,Y = np.meshgrid(x1[i],x2[i])
    ax.plot_surface(X,Y,data2d_time[0][tt][i]._data.T,
                    rstride=1,cstride=1, cmap=cm.coolwarm, alpha=0.6)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel(var_name)
ax.set_xlim3d(xmin, xmax)
#plt.savefig('test.png')
plt.show()
