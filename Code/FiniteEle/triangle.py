# -*- coding: UTF-8 -*-
import numpy as np
from FiniteEle.node import *
from FiniteEle.unit import *

#平面三角形问题
#内置ND=3
class Triangle:
    def __init__(self, node_dof, node_num, node_coord, uni_num, node_unit, force, constraints, materials, uni_ma, uni_type='rod'):
        ''' 
        force 外力载荷
        constraints 受约束的位移编号
        node_dof 节点自由度
        uni_num 单元个数
        node_num 节点个数
        node_unit 记录单元与节点的对应关系 3*uni_num
        uni_ma 记录单元对应的材料常数类别 [node_num] 
        '''
        self.force = force
        self.constraints = constraints
        self.node_dof = node_dof
        self.uni_num = uni_num
        self.node_num = node_num
        self.node_unit = node_unit
        self.x = 0
        
        nodes = []
        units = []
        
        for i in range(node_num):
            nodes.append(Node(i+1, node_dof, node_coord[:, i]))
        for i in range(uni_num):
            units.append(Unit(i+1, materials[uni_ma[i]-1], node_dof, nodes[node_unit[0,i]-1], nodes[node_unit[1,i]-1], 'triangle', nodes[node_unit[2,i]-1]))
            
        self.nodes = nodes
        self.units = units
        
        
    def s_change(self,i): #实现单元中节点位移的id与整体中节点位移id的转换
        if self.node_dof==3: #自由度为3
            s_change=np.zeros(6)
            s_change[0] = (self.node_unit[0,i]-1)*self.node_dof+1 
            s_change[1] = (self.node_unit[0,i]-1)*self.node_dof+2
            s_change[2] = (self.node_unit[0,i]-1)*self.node_dof+3
            
            s_change[3] = (self.node_unit[1,i]-1)*self.node_dof+1
            s_change[4] = (self.node_unit[1,i]-1)*self.node_dof+2
            s_change[5] = (self.node_unit[1,i]-1)*self.node_dof+3

            s_change[6] = (self.node_unit[2,i]-1)*self.node_dof+1
            s_change[7] = (self.node_unit[2,i]-1)*self.node_dof+2
            s_change[8] = (self.node_unit[2,i]-1)*self.node_dof+3
        
        if self.node_dof==2: #自由度为2
            s_change=np.zeros(6)
            s_change[0] = (self.node_unit[0,i]-1)*self.node_dof+1 
            s_change[1] = (self.node_unit[0,i]-1)*self.node_dof+2

            s_change[2] = (self.node_unit[1,i]-1)*self.node_dof+1
            s_change[3] = (self.node_unit[1,i]-1)*self.node_dof+2

            s_change[4] = (self.node_unit[2,i]-1)*self.node_dof+1
            s_change[5] = (self.node_unit[2,i]-1)*self.node_dof+2
            
        return s_change
        
    def ID_diag(self): #存储主角元在一维数组中的编号
        #计算并保存相关的最小节点号
        ID_diag = np.zeros(self.node_num*self.node_dof)
        ID_diag[0]=1
        
        #按节点循环
        for node in self.nodes:
        #for k in range(self.node_num):
            temp = 100 #中间变量记录最小节点号
            for uni in self.units:
            #for i in range(self.uni_num):
                for j in range(3): #一个单元三个节点
                    if self.node_unit[j,int(uni.id-1)] == node.id: #判断单元中是否有k点
                        for l in range(3):
                            if self.node_unit[l,int(uni.id-1)] < temp:
                                temp = self.node_unit[l,int(uni.id-1)] #寻找与k点相关的最小行号
            
            #确定k点所在行号的ID_diag数组
            for  i in range(self.node_dof):
                J = self.node_dof*(node.id-1)+i+1
                if J != 1: #递归获取ID_diag数组
                    ID_diag[J-1] = ID_diag[J-2]+self.node_dof*(node.id-temp)+i+1
        #print(len(ID_diag))        
        return ID_diag
    
    def AK_g(self): #组集合成结构刚度阵
        ID_diag = self.ID_diag() #获取ID_diag矩阵
        #print(ID_diag)
        AK = np.zeros(int(ID_diag[-1])) #ID_diag[-1]即为一维半带宽存储的刚度矩阵的容量
        for uni in self.units:
        #for m in range(self.uni_num):
            ''' node_i = unit.node_i.id
            node_j = unit.node_j.id '''
            s_change = self.s_change(uni.id-1)
            K = uni.Stiffness_g() #单元刚度阵（总体坐标下）
            for i in range(self.node_dof*3):
                for j in range(self.node_dof*3):
                    ''' row_i = (node_i-1) * self.node_dof + m
                    row_j = (node_j-1) * self.node_dof + m
                    AK[ID_diag[row_i]-m:ID_diag[row_i]+1] += K[m,0:m+1]
                    AK[ID_diag[row_j]-m:ID_diag[row_j]+1] += K[m+self.node_dof, self.node_dof:m+1+self.node_dof]
                    if node_i > node_j:
                        start_id = ID_diag[row_i] - m -self.node_dof * (node_i - node_j)
                    else:
                        start_id = ID_diag[row_j] - m -self.node_dof * (node_j - node_i)
                        
                    AK[start] '''
                    s_changeS = s_change[i]-s_change[j]
                    if s_changeS>=0: #将单元刚度阵的元素逐一加到总体结构刚度阵中
                        id = ID_diag[int(s_change[i]-1)]-(s_change[i]-s_change[j])
                        #print(id)
                        AK[int(id-1)] += K[i,j]
        #print(AK)
                    
        return AK
    
    def constraints_handling(self): #约束处理-置大数法
        ID_diag = self.ID_diag()
        AK = self.AK_g()
        for i in self.constraints:
            #print(i)
            row = i-1
            AK[int(ID_diag[row]-1)] = 1e30
        return ID_diag, AK
    
    def cholesky(self):  #乔累斯基法解矩阵方程
        ID_diag, AK = self.constraints_handling()
        #print(AK)
        F = self.force.copy()
        
        for i in range(self.node_dof*self.node_num): #按行循环
            if i != 0: #第一项不需要分解
                mi = i+1-(ID_diag[i]-ID_diag[i-1])+1 #第一个非零元列号
                if mi != i+1: #第i+1行非对角元情况
                    for j in range(int(mi-1),i): #在带宽内，按列循环，不包括对角元
                        if j != mi-1:
                            mj = j+1+1-(ID_diag[j]-ID_diag[j-1]) #第一个非零元列号mj
                            igp = ID_diag[i] - (i-j) #当前被分解元素Kij的一维编号
                            mij = mi if mj < mi else mj #取出较大的mij
                            jgp = ID_diag[j]-j-1
                            if mij <= j:
                                for k in range(int(mij-1), j):
                                    ik = ID_diag[i] - i - 1 + k + 1 #Lik的一维地址
                                    jk = jgp + k + 1 #Ljk的一维地址
                                    kk = ID_diag[k] #dkk的一维地址
                                    AK[int(igp-1)]-=AK[int(ik-1)]*AK[int(kk-1)]*AK[int(jk-1)] #计算分解的中间值
                        else: #j=mi-1
                            igp = ID_diag[i-1]+1 #最左端非零元的一维地址
                        ii = ID_diag[j] #对角元djj的一维地址
                        AK[int(igp-1)] = AK[int(igp-1)]/AK[int(ii-1)] #计算Lij
                        F[i]-=AK[int(igp-1)]*AK[int(ii-1)]*F[j] #计算分解载荷项
                    ij = ID_diag[i]
                    for k in range(int(mi-1),i):
                        ii = ID_diag[i] - i + k #Lik的一维地址
                        jj = ID_diag[k] #dkk的一维地址
                        AK[int(ij-1)] -= AK[int(ii-1)]**2 * AK[int(jj-1)] #计算dii
            ij = ID_diag[i] #对角元dii的一维地址
            F[i] = F[i]/AK[int(ij-1)] #完成载荷分解
        for m in range(1, self.node_dof*self.node_num): #回代求解
            i = self.node_num*self.node_dof-m
            mi = i+1-(ID_diag[i]-ID_diag[i-1])+1  #第i行最左端非零元素列号
            if (mi-1) != i: 
                iig = ID_diag[i]-i-1
                for k in range(int(mi-1), i):
                    ij = iig+k+1  
                    F[k] -= AK[int(ij-1)]*F[i] #求解位移
                
        self.x = F
        return F
    
    def Output(self):
        
        deg = np.zeros((self.uni_num, 3*self.node_dof)) #单元节点位移列阵
        #dee = np.zeros((self.uni_num, 2)) #单元节点位移（局部坐标系）
        #fee = np.zeros((self.uni_num, 2)) #单元节点力（局部坐标系）
        #sigma = np.zeros(self.uni_num) #单元应力
        f = np.zeros((self.uni_num, 3*self.node_dof)) #节点力（总体坐标）
        P = np.zeros(self.node_dof*self.node_num) #结构节点力阵
        R = np.zeros(self.node_dof*self.node_num) #约束反力
        
        for uni in self.units:
        #for i in range(self.uni_num):
            
            TK = uni.Stiffness_g()
            s_change = self.s_change(uni.id-1)
            
            for j in range(self.node_dof*3):
                deg[int(uni.id-1),j] = self.x[int(s_change[j]-1)]
            
                
            for j in range(self.node_dof*3):
                for k in range(2):
                    f[int(uni.id-1),j] += TK[j,k] * deg[int(uni.id-1),k]
                    
            for j in range(self.node_dof*2):
                P[int(s_change[j]-1)] += f[int(uni.id-1),j]
                
        for i in range(self.node_dof*self.node_num):
            R[i] = P[i] - self.force[i]
        
        #输出部分
        # 保留小数（四位）
        x = np.zeros(self.node_dof*self.node_num)
        for i in range(self.node_dof*self.node_num):  
            #print(1111)
            x[i] = round(self.x[i], 4)
            if x[i] == 0:
                x[i] = 0 

        for uni in self.units:

        #for i in range(self.uni_num):
            s_change = self.s_change(uni.id-1)
            for j in range(self.node_dof*3):  
                deg[int(uni.id-1),j] = round(deg[int(uni.id-1),j], 4)
                if deg[int(uni.id-1),j]==0:
                    deg[int(uni.id-1),j] = 0 
            ''' for j in range(2):  
                dee[int(uni.id-1),j] = round(dee[int(uni.id-1),j], 4)
                if dee[int(uni.id-1),j]==0:
                    dee[int(uni.id-1),j]=0 '''
            ''' for j in range(2):  # 单元内力
                fee[int(uni.id-1),j] = round(fee[int(uni.id-1),j], 4)
                if fee[int(uni.id-1),j]==0:
                    fee[int(uni.id-1),j]=0 '''
            ''' sigma[int(uni.id-1)] = round(sigma[int(uni.id-1)], 3)  # 单元应力
            if sigma[int(uni.id-1)]==0:
                sigma[int(uni.id-1)]=0 '''
            for j in range(3*self.node_dof):  
                f[int(uni.id-1),j] = round(f[int(uni.id-1),j], 4)
                if f[int(uni.id-1),j]==0:
                    f[int(uni.id-1),j]=0

        for i in range(self.node_dof*self.node_num):
            P[i] = round(P[i], 4)
            if P[i]==0:
                P[i]=0
        for i in range(self.node_dof*self.node_num):  
            R[i] = round(R[i], 4)
            if R[i]==0:
                R[i]=0


        # 输出
        data0 = []
        data0.append(['节点位移'])
        data0.append(x)
        data0.append('\n')

        ''' data0.append(['单元节点位移(局部坐标系)'])
        for i in range(self.uni_num):
            data0.append(dee[i])
        data0.append('\n') '''

        data0.append(['单元节点位移'])
        for i in range(self.uni_num):
            data0.append(deg[i])
        data0.append('\n')
        
        ''' data0.append(['单元内力'])
        for i in range(self.uni_num):
            data0.append(fee[i])
        data0.append('\n') '''

        ''' data0.append(['单元应力'])
        data0.append(sigma)
        data0.append('\n') '''

        data0.append(['单元节点力'])
        for i in range(self.uni_num):
            data0.append(f[i])
        data0.append('\n')

        data0.append(['结构节点力'])
        data0.append(P)
        data0.append('\n')

        data0.append(['约束反力'])
        data0.append(R)
        
        return data0

        
                    
                        
                
            
    
                    
                
                
            
                    
                        
                
        
        
        