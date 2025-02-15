# --*-- conding:utf-8 --*--
# @Time : 2/14/25 9:08 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : Quantum_dock.py

import numpy as np
from pyscf import gto, scf, mcscf, tools
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

    atom_str = '\n'.join([f"{atom} {x:.3f} {y:.3f} {z:.3f}" for atom, (x, y, z) in geometric])

    mol = gto.M(
        atom=atom_str,
        basis='sto-3g',
        charge=0,
        spin=0
    )

    print("总电子数:", mol.nelectron)
    print("总轨道数:", mol.nao)

    # 计算 Hartree-Fock 电子结构
    mf = scf.RHF(mol)
    mf.kernel()

    # # 自动选择活性空间
    # mo_occ = mf.mo_occ  # 轨道占据数
    # homo_idx = np.where(mo_occ > 0)[0][-1]  # 找到最后一个占据轨道
    # lumo_idx = homo_idx + 1  # LUMO 轨道索引

    n_active_orbitals = 40  # 选取 HOMO-3 ~ LUMO+3 共 8 个轨道
    n_active_electrons = 42  # 选取 10 个电子

    # print(f"自动选择 HOMO: {homo_idx}, LUMO: {lumo_idx}")
    print(f"选取 {n_active_electrons} 个电子，{n_active_orbitals} 个轨道")

    # Qiskit 计算哈密顿量
    driver = PySCFDriver(atom=mol.atom, basis="sto-3g", charge=0, spin=0)
    problem = driver.run()

    # 应用活性空间
    active_space = ActiveSpaceTransformer(
        num_electrons=n_active_electrons,
        num_spatial_orbitals=n_active_orbitals
    )
    problem = active_space.transform(problem)


    # 量子比特映射
    qubit_mapper = ParityMapper()
    qubit_op = qubit_mapper.map(problem.hamiltonian.second_q_op())

    print("Qubit Hamiltonian:", qubit_op)



