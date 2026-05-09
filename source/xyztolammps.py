import random
import numpy as np
# import matplotlib.pyplot as plt
import knot_xyz as xyz
import os

def read_xyz(file):
    with open(file) as f:
        while True:
            try:
                n_atoms = int(f.readline())
            except ValueError:
                break
            if not n_atoms:
                break
            f.readline()
            
            # 将坐标保存为narray
            coords = []
            for i in range(n_atoms):
                line = f.readline().split()
                coords.append([float(x) for x in line[1:]])
            yield coords
def write_lammps(data,filename,type="close",Lx=200,Ly=200,Lz=200):
    N=data.shape[0]
    with open(filename,'w') as f:
        f.write("#LAMMPS input file\n")
        f.write('{} atoms\n'.format(N))
        # 写入bond数目
        if (type=="open"):
            f.write('{} bonds\n'.format(N-1))
        elif(type=="close"):
            f.write('{} bonds\n'.format(N))
        # 写入angle数目
        if (type=="open"):
            f.write('{} angles\n'.format(N-2))
        elif(type=="close"):
            f.write('{} angles\n'.format(N))

        # 写入原子类型数目
        f.write('\n1 atom types\n')
        # 写入bond类型数目
        f.write('1 bond types\n')
        # 写入angle类型数目
        f.write('1 angle types\n')

        # 写入box的大小
        # min_x=np.min(data[:,0])
        # data[:,0] = data[:,0] - min_x
        # min_y=np.min(data[:,1])
        # data[:,1] = data[:,1] - min_y
        # min_z=np.min(data[:,2])
        # data[:,2] = data[:,2] - min_z
        # f.write('\n{} {} xlo xhi\n'.format(-300,300))
        # f.write('{} {} ylo yhi\n'.format(-300,300))
        # f.write('{} {} zlo zhi\n'.format(-300,300))
        print(np.min(data[:,0]),np.min(data[:,1]),np.min(data[:,2]))
        print(np.max(data[:,0]),np.max(data[:,1]),np.max(data[:,2]))

        # 写入质量
        f.write('\nMasses\n\n1 1.0\n')

        # 写入原子坐标
        f.write('\nAtoms\n\n')
        for i in range(N):
            f.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(i+1,1,1,data[i,0],data[i,1],data[i,2]))
        # 写入bond信息
        f.write('\nBonds\n\n')
        if (type=="open"):
            for i in range(N-1):
                f.write('{}\t{}\t{}\t{}\n'.format(i+1,1,i+1,i+2))
        elif(type=="close"):
            for i in range(N-1):
                f.write('{}\t{}\t{}\t{}\n'.format(i+1,1,i+1,i+2))
            f.write('{}\t{}\t{}\t{}\n'.format(N,1,N,1))
        # 写入angle信息
        f.write('\nAngles\n\n')
        if (type=="open"):
            for i in range(N-2):
                f.write('{}\t{}\t{}\t{}\t{}\n'.format(i+1,1,i+1,i+2,i+3))
        elif(type=="close"):
            for i in range(N-2):
                f.write('{}\t{}\t{}\t{}\t{}\n'.format(i+1,1,i+1,i+2,i+3))
            f.write('{}\t{}\t{}\t{}\t{}\n'.format(N-1,1,N-1,N,1))
            f.write('{}\t{}\t{}\t{}\t{}\n'.format(N,1,N,1,2))

# 定义计算两个向量的距离的函数
def distance(x1:np.ndarray, x2):
    return np.sqrt(np.sum((x1-x2)**2))

# 定义计算两个向量的角度的函数
def angle(x1:np.ndarray, x2:np.ndarray):
    return np.arccos(np.dot(x1,x2)/(np.linalg.norm(x1)*np.linalg.norm(x2)))

index=0
c=0
core=10
for data in read_xyz("tail.xyz"):
    L=len(data)
    # 定义总体范围
    population_range = range(0, L)
    # 设置采样次数
    num_samples = L


    population_list = list(population_range)
    data=np.array(data)

    # 用于存储每个segment的距离和角度
    dis_seg=[]
    ang_seg=[]
    
    flag = False
# 进行无放回采样
    # for out in range(0,L):
    #     if not population_list:
    #         break
    #     # 检查总体是否还有足够的元素可供采样
    #     # 从总体中采样一个元素
    #     sampled_index = random.randint(0, len(population_list) - 1)
    #     i = population_list[sampled_index]
    #     for j in range(0,L-1):
    #         if abs(i-j)>=core and abs(i-j)<=L-core:
    #             # 计算距离
    #             dis=distance(data[i],data[j])
    #             # 计算角度

    #             if dis<=2:
    #                 ang=angle(data[i-1]-data[i],data[j-1]-data[j])/3.14*180
                    
    #                 if ang>80 and ang<100:
    #                 # if dis>=0 and dis<=7 :
    #                     print(index,i,j)
    #                     print(data[i])
    #                     print(data[j])
    #                     data[[i, j]] = data[[j, i]]
    #                     data[[i-1, j-1]] = data[[j-1, i-1]]
    #                     data[[(i+1)%L, (j+1)%L]] = data[[(j+1)%L, (i+1)%L]]

    #                     # save_pos("top_{}.pos".format(L),data,'cenline from 41 knot Lknot{} f{}'.format(L,'data'),dir="changes/")
    #                     write_lammps(data,"lammps"+str(L)+".data",type="close",Lx=L,Ly=L,Lz=L)
    #                     # xyz.save_xyz("top_{}.xyz".format(L),data,'cenline from 41 knot Lknot{} f{}'.format(L,'data'),dir="changes/")
    #                     c+=1
    #                     flag = True
    #                     break        
    #     if flag:
    #         break
    #     del population_list[sampled_index]
    # index+=1


write_lammps(data,"lammps_{}_{}.data".format(type,L),type="close",Lx=L,Ly=L,Lz=L)