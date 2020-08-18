import numpy as np 
import tqdm 
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.pyplot as plt
import os
import io
from PIL import Image

dir_n = "C:/Users/samue/Documents/annealing_simulation/3d/output_files/thierno/"
prefix = "iter_"


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

def plot_viz(ff):
    fig = plt.figure(figsize=(20,14))
    ax = fig.gca(projection='3d')
    ff = ff[:,:,:]
    tod = ff >0
    colors = np.empty(tod.shape, dtype=object)
    colors[ff==2] = 'green'
    colors[ff==1] = 'blue'
    #colors[ff==0] = 'yellow'
    ax.voxels(tod,facecolors=colors, edgecolor='k')
    return fig

def get_im(dir_n, prefix, n):
    f_name = os.path.join(dir_n,prefix+str(n)+".txt")
    arr = treat(f_name)
    p = plot_viz(arr)
    
    buf = io.BytesIO()
    p.savefig(buf)
    buf.seek(0)
    
    img = Image.open(buf)
    return img

def get_height(dir_n,prefix,n):
    f_name = os.path.join(dir_n,prefix+str(n)+".txt")
    arr = treat(f_name)
    
    mx,my,mz=arr.shape
    depth = np.zeros((mx,my))
    for x in range(mx):
        for y in range(my):
            
            for z in range(mz-1,-1,-1):
                p = arr[x,y,z]
                if p>0:
                    depth[x,y] = z
                    break
    
    return depth

def get_himg(dir_n,prefix,n):
    tt = get_height(dir_n,prefix,n)
    f = plt.figure()
    plt.imshow(tt)
    buf= io.BytesIO()
    f.savefig(buf)
    plt.close(f)
    buf.seek(0)

    img = Image.open(buf)
    return img

l_img=[]
for n in tqdm.tqdm(range(11)):    
    img = get_himg(dir_n,prefix,n)
    l_img.append(img)

l_img[0].save(r'C:\Users\samue\Documents\annealing_simulation\3d\output_files\thierno\or_en4.gif', save_all=True, append_images=l_img[1:], duration=100, loop=0)    