# --*-- conding:utf-8 --*--
# @Time : 2/14/25 9:08 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : Quantum_dock.py
from numpy.random import geometric
from pyscf import gto, scf, mcscf
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import ParityMapper
from qiskit_nature.second_q.transformers import ActiveSpaceTransformer

def read_xyz_file(xyz_path):
    """
    读取简单的.xyz文件:
      第1行=原子数
      第2行=注释
      从第3行起= symbol  x  y  z
    返回 geometry = [(symbol,(x,y,z)), ...]
    """
    geometry = []
    with open(xyz_path, 'r') as f:
        lines = f.readlines()

    atom_lines = lines[2:]
    for line in atom_lines:
        parts = line.split()
        if len(parts) < 4:
            continue
        symbol = parts[0]
        x = float(parts[1])
        y = float(parts[2])
        z = float(parts[3])
        geometry.append( (symbol,(x,y,z)) )
    return geometry



if __name__=="__main__":

    geometric = read_xyz_file('./subsystems_xyz/sub_A_267_PHE/res_plus_ligand.xyz')

    print(type(geometric))
