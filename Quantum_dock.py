# --*-- conding:utf-8 --*--
# @Time : 2/14/25 9:08 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : Quantum_dock.py

from pyscf import gto, scf, mcscf
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import QubitConverter, ParityMapper
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

    # 选择活性空间
    n_active_orbitals = min(8, mol.nao)  # 不能超过总轨道数
    n_active_electrons = min(10, mol.nelectron)  # 不能超过总电子数

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

    # 量子比特转换
    qubit_converter = QubitConverter(ParityMapper())
    qubit_op = qubit_converter.convert(problem.hamiltonian)

    print("Qubit Hamiltonian:", qubit_op)



