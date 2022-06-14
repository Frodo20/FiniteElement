# -*- coding: UTF-8 -*-
import numpy as np

#节点类
class Node:
    def __init__(self, node_id, dof, node_coord):
        self.id = node_id #节点id
        self.dof =dof #节点自由度
        #节点坐标
        self.x = node_coord[0]
        self.y = node_coord[1] if dof > 1 else 0
        if dof>2:
            self.z = node_coord[2]
        else:
            self.z = 0