# --*-- conding:utf-8 --*--
# @Time : 2/12/25 10:13 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : Molecular_Docking.py

import numpy as np

# ================
# 1) 导入Qiskit Nature相关模块
# ================
from qiskit_nature.drivers.second_quantization import ElectronicStructureMoleculeDriver
from qiskit_nature.drivers import Molecule
from qiskit_nature.units import DistanceUnit
from qiskit_nature.problems.second_quantization import ElectronicStructureProblem
from qiskit_nature.transformers.second_quantization.electronic import FreezeCoreTransformer
from qiskit_nature.converters.second_quantization import QubitConverter
from qiskit_nature.mappers.second_quantization import JordanWignerMapper

# VQE相关
from qiskit.algorithms import VQE
from qiskit.algorithms.optimizers import COBYLA
from qiskit_nature.circuit.library import UCCSD

# Qiskit 量子模拟器
from qiskit import Aer
from qiskit.utils import QuantumInstance


# ================
# 2) 定义一个帮助函数：给定Molecule对象，返回基态能量
# ================
def compute_ground_state_energy(molecule: Molecule, basis: str = "sto3g") -> float:
    """
    给定一个分子Molecule对象和基组，返回使用VQE计算出的基态能量（单位：Hartree）。
    """
    # 1) 构造驱动
    driver = ElectronicStructureMoleculeDriver(
        molecule=molecule,
        basis=basis,
        driver_type="PYSCF"
    )
    # 2) 构建电子问题
    problem = ElectronicStructureProblem(driver)

    # 3) (可选) 冻结核（FreezeCore）等近似，减小体系规模
    problem = FreezeCoreTransformer().transform(problem)

    # 4) 获取费米子哈密顿量
    second_q_ops = problem.second_q_ops()
    main_op = second_q_ops[0]

    # 5) 映射到量子比特 (Jordan-Wigner)
    num_particles = problem.molecule_data_transformed.num_alpha, problem.molecule_data_transformed.num_beta
    num_spin_orbitals = problem.molecule_data_transformed.num_spin_orbitals

    qubit_converter = QubitConverter(
        mapper=JordanWignerMapper(),
        two_qubit_reduction=False
    )

    # 将费米子哈密顿量转为量子比特哈密顿量
    qubit_hamiltonian = qubit_converter.convert(main_op, num_particles=num_particles)

    # 6) 定义VQE, 选用 UCCSD ansatz
    ansatz = UCCSD(
        qubit_converter=qubit_converter,
        num_particles=num_particles,
        num_spin_orbitals=num_spin_orbitals,
    )

    optimizer = COBYLA(maxiter=200)

    vqe_solver = VQE(
        ansatz=ansatz,
        optimizer=optimizer
    )

    # 7) 选用本地模拟器：Statevector 或者就 Qasm 都行，这里用无噪声 statevector
    quantum_instance = QuantumInstance(
        backend=Aer.get_backend("statevector_simulator"),
        shots=1
    )

    vqe_solver.quantum_instance = quantum_instance

    # 8) 运行VQE，得到基态能量
    result = vqe_solver.compute_minimum_eigenvalue(qubit_hamiltonian)
    energy_hartree = result.eigenvalue.real

    return energy_hartree

if __name__ == '__main__':

    # ================
    # 3) 准备简化的“蛋白-配体”三种结构
    # ================
    # 演示用：将H2O当作“receptor”，H2当作“ligand”，再把H2O和H2放在一起做“complex”
    # 坐标单位：Å (Angstrom)；坐标可能随意放置以演示
    # ================

    # -- Receptor (H2O) --
    receptor_molecule = Molecule(
        geometry=[
            ("O", (0.0000, 0.0000, 0.0000)),
            ("H", (0.7572, 0.5860, 0.0000)),
            ("H", (-0.7572, 0.5860, 0.0000))
        ],
        multiplicity=1,    # 单重态
        charge=0,
        units=DistanceUnit.ANGSTROM
    )

    # -- Ligand (H2) --
    ligand_molecule = Molecule(
        geometry=[
            ("H", (0.0000, 0.0000, 0.0000)),
            ("H", (0.0000, 0.0000, 0.74))  # ~0.74 Å
        ],
        multiplicity=1,
        charge=0,
        units=DistanceUnit.ANGSTROM
    )

    # -- Complex (H2O + H2)，把两个分子稍微分开一点 --
    complex_molecule = Molecule(
        geometry=[
            # H2O
            ("O", (0.0000, 0.0000, 0.0000)),
            ("H", (0.7572, 0.5860, 0.0000)),
            ("H", (-0.7572, 0.5860, 0.0000)),
            # H2 (放在X轴+2.0 Å处，使其与水分子保持一定距离)
            ("H", (2.0000, 0.0000, 0.0000)),
            ("H", (2.0000, 0.0000, 0.74))
        ],
        multiplicity=1,
        charge=0,
        units=DistanceUnit.ANGSTROM
    )


    # ================
    # 4) 分别计算三者的基态能量
    # ================
    E_receptor = compute_ground_state_energy(receptor_molecule, basis="sto3g")
    E_ligand   = compute_ground_state_energy(ligand_molecule,   basis="sto3g")
    E_complex  = compute_ground_state_energy(complex_molecule,  basis="sto3g")

    print("Receptor (H2O) Energy:  ", E_receptor, "Hartree")
    print("Ligand   (H2)  Energy:  ", E_ligand,   "Hartree")
    print("Complex (H2O+H2) Energy:", E_complex,  "Hartree")

    # 计算“结合能”
    E_binding = E_complex - (E_receptor + E_ligand)
    print("Binding Energy (Hartree) =", E_binding)
