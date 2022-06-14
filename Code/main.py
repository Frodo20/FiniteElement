from FiniteEle.rod import *
from FiniteEle.triangle import *
import os
import numpy as np

'''
输入数据包括（以杆的二维问题为例）
    ND = 2 单元节点个数内置为2

    Input.txt

    node_dof = data[0][0]  节点自由度数(一般为二维或三维)
    node_num = data[1][0]  节点个数
    uni_num = data[2][0]  单元个数
    kind_ma = data[3][0]  单元的种类数
    X = data[4]  节点x坐标
    Y = data[5]  节点y坐标
    ndoe_uni = [data[6], data[7]]  每个单元对应节点编号
    constraints = data[8]  受约束的位移分量编号
    material[i] = [data[9][i], data[10][i]]  特征类
    uni_ma = data[11]  每个单元对应的特征数编号
    force = data[12]  载荷数据
'''

'''
输出数据包括（以杆的二维问题为例）

    Result.txt

    节点位移  node_dof*node_num
    单元节点位移(局部坐标) 2*uni_num
    单元节点位移（整体坐标） 4*uni_num
    两者区别在于, 局部坐标系中仅有一个坐标x, 总体坐标中有坐标xy
    单元内力 2*uni_num
    单元应力 uni_num
    单元节点力 4*uni_num
    结构节点力 node_num
    约束反力 node_num

'''
rootP = './Input/'

def readInput(file_name):  #获取Input
    data = []
    f = open(rootP+file_name, 'r')  
    data0 = f.readlines()  
    for row in data0:
        tmp = row.split('\t')  
        tmp[-1] = tmp[-1].replace('\n', '')  #删去换行符
        for col in range(len(tmp)):
            tmp[col] = float(tmp[col]) #将字符串强制转换为浮点数，包括坐标、弹性模量、面积、力
        data.append(tmp)  
    
    return data

def dataoper(data, type=1): #对获取数据进行处理，将某些值化成整数
    for row in range(len(data)):
        for col in range(len(data[row])):
            if type==1:
                if row!=4 and row!=8 and row!=9 and row!=11:
                    data[row][col] = int(data[row][col])
            if type==2:
                if row!=4 and row!=5 and row!=9 and row!=10 and row!=12:
                    data[row][col] = int(data[row][col])
            if type==3:
                if row!=4 and row!=5 and row!=6 and row!=10 and row!=11 and row!=13:
                    data[row][col] = int(data[row][col])
            if type==4:
                if row!=4 and row!=5 and row!=10 and row!=11 and row!=12 and row!=14:
                    data[row][col] = int(data[row][col])
    return data

def out_result(data, file_name):  #输出output
    f = open('./Output/'+file_name, 'w')  
    for row in range(len(data)):
        for col in range(len(data[row])):
            s = str(data[row][col])
            f.write(s)
            f.write('\t')
        f.write('\n')
    f.close()

if __name__ == '__main__':
    ''' node_coord = np.array([[0, 0, 1, 1],
                          [-1, 0, 0, -1]])
    node_unit = np.array([[1, 2, 3, 1, 2, 1],
                        [2, 3, 4, 4, 4, 3]])
    force = [0, 0, 0, 0, 0, -1, 0, 0]
    constraints = [1, 2, 3, 4]
    materials = [Material(1,1)]
    uni_ma = [1, 1, 1, 1, 1, 1]
    
    rod1 = Rod(node_dof=2, node_num=4, node_coord=node_coord, uni_num=6, node_unit=node_unit,
               force = force, constraints = constraints, materials=materials, uni_ma=uni_ma, uni_type='rod')
    
    
    x=rod1.cholesky()
    print(x) '''
    for root, dirs, files in os.walk(rootP):
        #print(files)
        for file in files:
            #print(file)
            if (file=='RodInput.txt'):
                data = readInput(file)
                #print(data)
                if data[0][0]==1: #一维问题
                    data = dataoper(data,1)
                    node_coord = np.array([data[4]])
                    node_unit = np.array([data[5],
                                        data[6]])
                    force = data[11]
                    constraints = data[7]
                    kind_ma = data[3][0]
                    materials = []
                    for i  in range(kind_ma):
                        materials.append(Material(data[8][i],data[9][i]))
                    uni_ma = data[10]

                if data[0][0]==2: #二维问题
                    data = dataoper(data,2)
                    node_coord = np.array([data[4],
                                        data[5]])
                    node_unit = np.array([data[6],
                                        data[7]])
                    force = data[12]
                    constraints = data[8]
                    kind_ma = data[3][0]
                    materials = []
                    for i  in range(kind_ma):
                        materials.append(Material(data[9][i],data[10][i]))
                    uni_ma = data[11]
                    
                if data[0][0]==3: #三维问题
                    data = dataoper(data,3)
                    node_coord = np.array([data[4],
                                        data[5],
                                        data[6]])
                    node_unit = np.array([data[7],
                                        data[8]])
                    force = data[13]
                    constraints = data[9]
                    kind_ma = data[3][0]
                    materials = []
                    for i  in range(kind_ma):
                        materials.append(Material(data[10][i],data[11][i]))
                    uni_ma = data[12]
                    
                rod = Rod(node_dof=data[0][0], node_num=data[1][0], node_coord=node_coord, uni_num=data[2][0], node_unit=node_unit,
                force = force, constraints = constraints, materials=materials, uni_ma=uni_ma, uni_type='rod')

                x=rod.cholesky()

                data0 = rod.Output()
                out_result(data0, 'RodResult.txt')

            if (file=='TriangleInput.txt'):
                print(111)
                data1 = readInput(file)

                if data1[0][0] ==2: #三角形单元默认二维问题
                    data1 = dataoper(data1,4)
                    node_coord1 = np.array([data1[4],
                                        data1[5]])
                    node_unit1 = np.array([data1[6],
                                        data1[7],
                                        data1[8]])
                    force1 = data1[14]
                    constraints1 = data1[9]
                    kind_ma1 = data1[3][0]
                    materials1 = []
                    for i  in range(kind_ma1):
                        materials1.append(Material(data1[10][i],data1[11][i],data1[12][i]))
                    uni_ma1 = data1[13]

                    triangle = Triangle(node_dof=data1[0][0], node_num=data1[1][0], node_coord=node_coord1, 
                    uni_num=data1[2][0], node_unit=node_unit1, force=force1, constraints=constraints1, materials=materials1, uni_ma=uni_ma1)

                    x=triangle.cholesky()

                    data2 = triangle.Output()
                    out_result(data2, 'TriangleResult.txt')
        
    #print(x)