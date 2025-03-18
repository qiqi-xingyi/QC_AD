# --*-- conding:utf-8 --*--
# @Time : 3/14/25 2:03 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : Bing_E_Quantum.py

# Test group: 1c5z

import os
import csv
import json
import numpy as np
from typing import List, Tuple

# Qiskit & Qiskit Nature
from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import ParityMapper
from qiskit_nature.second_q.circuit.library import HartreeFock, UCCSD
from qiskit_nature.second_q.problems import ElectronicStructureProblem
from qiskit_nature.second_q.transformers import ActiveSpaceTransformer

# IBM Quantum runtime
from qiskit_ibm_runtime import QiskitRuntimeService, Session
try:
    from qiskit_ibm_runtime import Estimator as EstimatorV2
except ImportError:
    from qiskit.primitives import Estimator as EstimatorV2

# PySCF for initial SCF
from pyscf import gto, scf
from scipy.optimize import minimize
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager


# -----------------------------
# 1) 读取配置 (TOKEN, INSTANCE)
# -----------------------------
def read_config(file_path):
    config = {}
    with open(file_path, "r") as f:
        for line in f:
            if "=" in line:
                key, value = line.strip().split("=")
                config[key.strip()] = value.strip()
    return config


# -----------------------------
# 2) 将完整 PDB 转为 XYZ (保留所有原子)
# -----------------------------
def pdb_to_xyz(pdb_path: str) -> str:
    """
    读取整个PDB（蛋白+配体），将所有 ATOM/HETATM 行转换为
    PySCFDriver 能识别的 xyz 字符串: "C x y z; O x y z; ..."
    不做任何物理截断。
    """
    xyz_list = []
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # 识别元素
                atom_element = line[76:78].strip()
                if not atom_element:
                    atom_element = line[12:16].strip(" 1234567890")

                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                xyz_list.append(f"{atom_element} {x:.3f} {y:.3f} {z:.3f}")
    return "; ".join(xyz_list)


# -----------------------------
# 3) 解析 PLIP 结果以获取“关键残基 + 配体”信息
#    （仅用于估算活性电子数，不做物理截断）
# -----------------------------
def parse_plip_interaction_file(plip_path: str):
    """
    在这里根据你的PLIP分析文件(1c5z_interaction.txt)，
    找到与配体相互作用的重要残基。

    示例：文件中显示SER B 190, VAL B 213, GLY B 219, ...
    以及配体“MOL A 1”。

    你可以用更复杂的正则或表格解析，这里演示写死或简单匹配。
    返回两个对象:
      - residue_list = [("SER","B","190"), ("VAL","B","213"), ...]
      - ligand_info = ("MOL", "A", "1")
    """
    # 此处直接写死示例:
    residue_list = [
        ("SER", "B", "190"),
        ("VAL", "B", "213"),
        ("GLY", "B", "219")
    ]
    ligand_info = ("MOL", "A", "1")

    return residue_list, ligand_info


# -----------------------------
# 4) 根据 PLIP 结果，统计关键残基+配体的“价电子总数”
#    然后用以估计“活性电子数”和“活性轨道数”。
# -----------------------------
VALENCE_DICT = {
    "H": 1, "He": 2,
    "Li": 1, "Be": 2, "B": 3, "C": 4, "N": 5, "O": 6, "F": 7, "Ne": 8,
    "Na":1, "Mg":2, "Al":3, "Si":4, "P":5, "S":6, "Cl":7, "Ar":8,
}

def estimate_active_space_from_plip(pdb_path: str,
                                    residue_list: List[Tuple[str,str,str]],
                                    ligand_info: Tuple[str,str,str],
                                    additional_orbitals=4):
    """
    1) 读取 PDB 全部行，但只在“关键残基 + 配体”行上累加其原子价电子；
    2) 将这些电子数视作活性电子数 active_e；
    3) 初步设定 active_o = active_e // 2 + additional_orbitals (可自由调整)。

    注意：这是非常粗略的做法。也可以更复杂：先做 Mulliken 分析或局部轨道分析。
    """
    # residue_list = [("SER","B","190"), ...]
    # ligand_info  = ("MOL","A","1")

    # 转成可检索结构
    residue_set = set(residue_list)
    lig_resname, lig_chain, lig_resnum = ligand_info

    active_e = 0

    with open(pdb_path, "r") as f:
        for line in f:
            if (line.startswith("ATOM") or line.startswith("HETATM")):
                res_name = line[17:20].strip()
                chain_id = line[21].strip()
                res_num  = line[22:26].strip()

                # 判断此原子是否属“关键残基”or“配体”
                is_residue = (res_name, chain_id, res_num) in residue_set
                is_ligand  = (res_name == lig_resname and
                              chain_id == lig_chain and
                              res_num  == lig_resnum)
                if not (is_residue or is_ligand):
                    continue

                # 若是关键残基/配体，则加上其价电子
                atom_element = line[76:78].strip()
                if not atom_element:
                    atom_element = line[12:16].strip(" 1234567890")

                val_e = VALENCE_DICT.get(atom_element, 0)
                active_e += val_e

    # 粗略设定活性轨道数
    # 例如: active_orbitals = active_e//2 + additional_orbitals
    active_o = active_e//2 + additional_orbitals
    return active_e, active_o


# -----------------------------
# 5) 自定义 VQE 类
# -----------------------------
class CustomVQE:
    def __init__(self, service, qubit_op, ansatz,
                 shots=50, min_qubit_num=10, maxiter=20, optimization_level=3):
        self.service = service
        self.qubit_op = qubit_op
        self.ansatz = ansatz

        self.shots = shots
        self.min_qubit_num = min_qubit_num
        self.maxiter = maxiter
        self.optimization_level = optimization_level

        self.backend = self._select_backend(self.min_qubit_num)
        self.energy_list = []
        self.iter_count = 0

    def _select_backend(self, min_qubits):
        backend = self.service.least_busy(
            simulator=False,
            operational=True,
            min_num_qubits=min_qubits
        )
        return backend

    def _generate_pass_manager(self):
        pm = generate_preset_pass_manager(
            target=self.backend.target,
            optimization_level=self.optimization_level
        )
        return pm

    def cost_func(self, params, ansatz_isa, hamiltonian_isa, estimator):
        pub = (ansatz_isa, [hamiltonian_isa], [params])
        result = estimator.run(pubs=[pub]).result()
        energy = result[0].data.evs[0]

        self.iter_count += 1
        self.energy_list.append(energy)
        print(f"Iteration {self.iter_count}, Energy = {energy}")
        return energy

    def run_vqe(self):
        pm = self._generate_pass_manager()
        ansatz_isa = pm.run(self.ansatz)
        hamiltonian_isa = self.qubit_op.apply_layout(layout=ansatz_isa.layout)

        x0 = np.random.random(ansatz_isa.num_parameters)

        with Session(backend=self.backend) as session:
            estimator = EstimatorV2(mode=session)
            estimator.options.default_shots = self.shots

            res = minimize(
                self.cost_func,
                x0,
                args=(ansatz_isa, hamiltonian_isa, estimator),
                method="cobyla",
                options={'maxiter': self.maxiter}
            )
        return self.energy_list, res.x


# -----------------------------
# 6) 主流程：保留完整PDB，但用 PLIP 结果自动估计活性空间
# -----------------------------
def main():
    # (A) 读取配置
    config = read_config("config.txt")
    service = QiskitRuntimeService(
        channel='ibm_quantum',
        instance=config["INSTANCE"],
        token=config["TOKEN"]
    )

    # (B) PDB/PLIP 文件路径
    pdb_path  = "./data_set/data/2_benchmark_binidng_sites/1c5z/1c5z_Binding_mode.pdb"
    plip_path = "./data_set/data/2_benchmark_binidng_sites/1c5z/1c5z_interaction.txt"

    # (C) 将整个体系转XYZ
    xyz_str = pdb_to_xyz(pdb_path)

    # (D) 根据 PLIP 结果获取关键残基+配体
    residue_list, ligand_info = parse_plip_interaction_file(plip_path)

    # (E) 估算活性电子数 & 轨道数 (active_e, active_o)
    #     (仅用于后续 ActiveSpaceTransformer)
    active_e, active_o = estimate_active_space_from_plip(
        pdb_path,
        residue_list,
        ligand_info,
        additional_orbitals=4   # 你想多留几条轨道
    )
    print(f"[PLIP-based] Est. active_e = {active_e}, active_o = {active_o}")

    # (F) 对完整体系做 PySCFDriver -> ElectronicStructureProblem
    #     (这里你可指定真实总电荷和自旋；示例中写0)
    total_charge = 0
    spin         = 0

    driver = PySCFDriver(
        atom=xyz_str,
        basis="sto3g",
        charge=total_charge,
        spin=spin,
        unit=DistanceUnit.ANGSTROM
    )
    es_problem = driver.run()
    print("Created ElectronicStructureProblem for full system.")

    # (G) ActiveSpaceTransformer: 用 (active_e, active_o)
    #    让 Qiskit 在完整 SCF 波函数基础上，只相关处理这部分电子/轨道
    ast = ActiveSpaceTransformer(
        num_electrons=active_e,
        num_spatial_orbitals=active_o
    )
    red_problem = ast.transform(es_problem)

    # (H) 映射到量子比特
    from qiskit_nature.second_q.mappers import ParityMapper
    op = red_problem.hamiltonian.second_q_op()
    mapper = ParityMapper()
    qubit_op = mapper.map(op)

    # (I) 构造 UCCSD ansatz
    num_spatial_orbs = red_problem.num_spatial_orbitals
    num_alpha = red_problem.num_alpha
    num_beta  = red_problem.num_beta
    hf_init = HartreeFock(num_spatial_orbs, (num_alpha, num_beta), mapper)
    ansatz = UCCSD(
        num_spatial_orbitals=num_spatial_orbs,
        num_particles=(num_alpha, num_beta),
        mapper=mapper,
        initial_state=hf_init
    )

    # (J) 自定义 VQE
    vqe = CustomVQE(
        service=service,
        qubit_op=qubit_op,
        ansatz=ansatz,
        shots=100,
        min_qubit_num=30,  # 你可根据设备规模调大
        maxiter=15,
        optimization_level=3
    )
    energy_list, best_params = vqe.run_vqe()
    final_energy = energy_list[-1]
    print(f"\nFinal Energy (ActiveSpace) = {final_energy} Hartree")

    # (K) 保存结果
    out_dir = "results_plip_active"
    os.makedirs(out_dir, exist_ok=True)

    with open(os.path.join(out_dir,"energy_iterations.csv"), "w", newline="") as cf:
        writer = csv.writer(cf)
        writer.writerow(["Iter","Energy"])
        for i,e in enumerate(energy_list):
            writer.writerow([i+1, e])

    with open(os.path.join(out_dir,"optimized_params.json"), "w") as jf:
        json.dump({"best_params": best_params.tolist()}, jf, indent=4)

    print(f"Done. Results saved to {os.path.abspath(out_dir)}")


if __name__=="__main__":
    main()


