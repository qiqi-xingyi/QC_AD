# --*-- conding:utf-8 --*--
# @Time : 3/18/25 1:43 AM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : active_space_selector.py

import numpy as np
from pyscf import scf
import re

class ActiveSpaceSelector:
    def __init__(self, threshold=0.2):
        self.threshold = threshold  # ratio threshold for MO selection

    def run_scf(self, mol):
        if mol.spin == 0:
            mf = scf.RHF(mol)
        else:
            mf = scf.ROHF(mol)
        mf.kernel()
        return mf

    def select_active_space(self, mol, mf, residue_list, ligand_info, pdb_path):
        """
        1. 标记关键原子
        2. 投影 ratio >= threshold
        3. 返回 (active_e, mo_count, mo_start)
        """
        # => 同你现有的 find_active_space_by_projection() 类似
        # 这里可以写成一个私有方法，也可以直接写在此处

        # ...
        # 伪代码：
        # active_e, active_o, mo_start = your_projection_function(...)
        # return active_e, active_o, mo_start
