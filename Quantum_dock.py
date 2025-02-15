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

    # 格式化原子数据
    atom_str = '\n'.join([f"{atom} {x:.3f} {y:.3f} {z:.3f}" for atom, (x, y, z) in geometric])

    # 生成 PySCF 分子对象
    mol = gto.M(
        atom=atom_str,
        basis='sto-3g',  # 可更改为 '6-31G(d)'
        charge=0,
        spin=0
    )

    # 打印分子信息
    print("总电子数:", mol.nelectron)
    print("总轨道数:", mol.nao)

    # 计算 Hartree-Fock 电子结构
    mf = scf.RHF(mol)
    mf.kernel()

    # CASSCF 预计算
    casscf_guess_orbitals = 20  # 先用 20 轨道计算密度矩阵
    mc = mcscf.CASSCF(mf, casscf_guess_orbitals, mol.nelectron)
    mc.kernel()

    # 计算 Mulliken 电子密度，并导出 NOONs
    dm = mc.make_rdm1()
    noons = mc.mo_occ  # 获取自然轨道占据数 (NOONs)

    # 选择 NOONs 在 0.02 - 1.98 之间的轨道
    active_orbitals = [i for i, occ in enumerate(noons) if 0.02 < occ < 1.98]
    n_active_orbitals = len(active_orbitals)

    # 计算活性电子数
    n_active_electrons = int(sum(np.round([occ for occ in noons if 0.02 < occ < 1.98])))
    # **保证活性电子数 <= 2 * 活性轨道数**
    n_active_electrons = min(n_active_electrons, 2 * n_active_orbitals)

    print(f"自动选择的活性轨道数: {n_active_orbitals}")
    print(f"自动选择的活性电子数: {n_active_electrons}")

    # 重新运行 CASSCF 计算（优化活性空间）
    mc = mcscf.CASSCF(mf, n_active_orbitals, n_active_electrons)
    mc.kernel()

    # 生成 Qiskit Nature PySCFDriver
    driver = PySCFDriver(molecule=mol, basis="sto-3g", charge=0, spin=0)
    problem = driver.run()

    # 限定活性空间
    active_space = ActiveSpaceTransformer(
        num_electrons=n_active_electrons,
        num_spatial_orbitals=n_active_orbitals
    )
    problem = active_space.transform(problem)

    qubit_mapper = ParityMapper()
    qubit_op = qubit_mapper.map(problem.hamiltonian.second_q_op())

    print("Qubit Hamiltonian:", qubit_op)



