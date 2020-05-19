from mpl_toolkits.mplot3d import Axes3D 
import sys
import getopt
import matplotlib.pyplot as plt
import numpy as np

def treat(f_name):
    with open(f_name, "r") as f:
        l = f.readlines()
    dims = l[0].split(" ")[:-1]
    dims = list(map(int, dims))
    data = l[1].split(" ")[:-1]
    data= list(map(int, data))
    
    arr = np.array(data)
    arr = np.reshape(arr, tuple(dims))
    return arr

def plot_viz(ff, output_name):
    fig = plt.figure(figsize=(20,14))
    ax = fig.gca(projection='3d')
    ff = ff[:,:,:]
    tod = ff >0
    colors = np.empty(tod.shape, dtype=object)
    colors[ff==2] = 'blue'
    colors[ff==1] = 'blue'
    #colors[ff==0] = 'yellow'
    ax.voxels(tod,facecolors=colors, edgecolor='k')
    plt.savefig(output_name)

if __name__=="__main__":
    argv= sys.argv[1:]
    inputfile = ''
    outputfile = ''
    try:
       opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
       print( 'test.py -i <inputfile> -o <outputfile>')
       sys.exit(2)
    for opt, arg in opts:
       if opt == '-h':
          print( 'test.py -i <inputfile> -o <outputfile>')
          sys.exit()
       elif opt in ("-i", "--ifile"):
          inputfile = arg
       elif opt in ("-o", "--ofile"):
          outputfile = arg
    
    lattice = treat(inputfile)
    plot_viz(lattice, outputfile)

