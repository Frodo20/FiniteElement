# -*- coding: UTF-8 -*-
import numpy as np

#材料类，用于获取杆件的材料系数
class Material:
    def __init__(self, e=0, a=0, mu=0):
        self.E = e #弹性模量
        self.A = a #杆的面积
        self.mu = mu

#单元类，用于获取单元的物理常数，包括杆长、刚度、坐标变换等
class Unit:
    def __init__(self,  unit_id, material, node_dof, node_i, node_j, u_type = 'rod', node_k=0, node_l=0):
        if u_type == 'triangle': #三角形类
            self.node_k = node_k #节点3
        if u_type == 'quadrangle': #四边形类
            self.node_k = node_k
            self.node_l = node_l
        self.id = unit_id #单元编号
        self.node_i = node_i #节点1
        self.node_j = node_j #节点2
        self.type = u_type #单元类型
        self.material = material #单元材料常数
        self.dof = node_dof #节点自由度
        
    def Length(self, node_i, node_j): #获取单元长度
        return np.sqrt((node_i.x-node_j.x)**2 + (node_i.y-node_j.y)**2 + (node_i.z-node_j.z)**2)
    
    def triangle_area(self): #获取三角形的面积
        a = self.Length(self.node_i, self.node_j)
        b = self.Length(self.node_j, self.node_k)
        c = self.Length(self.node_i, self.node_k)
        p = (a+b+c)/2
        return np.sqrt( p*(p-a)*(p-b)*(p-c))

    def triangle_B(self):
        bi = self.node_j.y - self.node_k.y; ci = self.node_k.x - self.node_j.x
        bj = self.node_k.y - self.node_i.y; cj = self.node_i.x - self.node_k.x
        bk = self.node_i.y - self.node_j.y; cm = self.node_j.x - self.node_i.x
        A = self.triangle_area()
        B =1/2/A*np.array(
            [
                [bi, 0, bj, 0, bk, 0],
                [0, ci, 0, cj, 0, cm],
                [ci, bi, cj, bj, cm, bk]
            ]
        )
        return B

    def D(self, d_type=0):
        mu = self.material.mu
        e = self.material.E
        if d_type == 0:  # 平面应变
            d = np.array(
                [
                    [1-mu, mu, 0],
                    [mu, 1-mu, 0],
                    [0, 0, (1 - 2*mu)/2]
                ]
            )
            return d*e/(1+mu)/(1-2*mu)
        elif d_type == 1:  # 平面应力
            d = np.array(
                [
                    [1, mu, 0],
                    [mu, 1, 0],
                    [0, 0, (1 - mu)/2]
                ]
            )
            return d*e/(1-mu**2)


    def Stiffness_u(self): #获取单元刚度阵（局部坐标）
        if self.type == 'rod':
            #局部坐标系中： i -> j
            e = self.material.E
            a = self.material.A
            l = self.Length(self.node_i,self.node_j)
            #print(l)
    
            k_u = np.array([[e*a/l, -e*a/l], [-e*a/l, e*a/l]])

        
    def TransformMatrix(self): #计算坐标变换矩阵
        #局部坐标系中： j-i
        l = self.Length(self.node_i,self.node_j)
        if self.type == 'rod':
            if self.dof == 3:
                T = np.zeros((2,6))
                T[0,0] = (self.node_j.x - self.node_i.x) / l
                T[0,1] = (self.node_j.y - self.node_i.y) / l
                T[0,2] = (self.node_j.z - self.node_i.z) / l
                T[1,3] = (self.node_j.x - self.node_i.x) / l
                T[1,4] = (self.node_j.y - self.node_i.y) / l
                T[1,5] = (self.node_j.z - self.node_i.z) / l
                
            if self.dof == 2:
                T  = np.zeros((2,4))
                T[0,0] = (self.node_j.x - self.node_i.x) / l
                T[0,1] = (self.node_j.y - self.node_i.y) / l
                T[1,2] = (self.node_j.x - self.node_i.x) / l
                T[1,3] = (self.node_j.y - self.node_i.y) / l
        #print(T)        
        return T
        
    def TransformMatrixT(self): #计算变换矩阵的转置矩阵，实际上就是求转置
        
        T = self.TransformMatrix()
        
        if self.type == 'rod':
            
            TT = np.zeros((2*self.dof,2))
            for i in range(2):
                for j in range(2*self.dof):
                    TT[j][i] = T[i][j]  
        #print(TT)        
        return TT
        
    def Stiffness_g(self): #计算单元刚度阵（总体坐标系下）利用坐标变换 K_g=TT*k_u*T
        #print(11111111)
        k_g=0
        if self.type == 'rod':
            
            k_g = np.zeros((2*self.dof,2*self.dof))
            S = np.zeros((2*self.dof,2))
            k_u = self.Stiffness_u()
            T = self.TransformMatrix()
            TT = self.TransformMatrixT()
            
            for i in range(2*self.dof):
                for j in range(2):
                    for k in range(2):
                        S[i,j]+=TT[i,k]*k_u[k,j]
                    
            for i in range(2*self.dof):
                for j in range(2*self.dof):
                    for k in range(2):
                        k_g[i,j]+=S[i,k]*T[k,j] 

        if self.type == 'triangle':
            t = 1 #平面三角形单元 厚度默认取1
            B = self.triangle_B()
            A = self.triangle_area()
            D = self.D()
            k_g = np.dot(B.T, np.dot(D, B))*A*t

        #print(k_g) 
        return k_g