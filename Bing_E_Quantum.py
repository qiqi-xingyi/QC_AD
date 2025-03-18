# --*-- conding:utf-8 --*--
# @Time : 3/14/25 2:03 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : Bing_E_Quantum.py

# Test group: 1c5z

import os
import re
import csv
import json
import numpy as np
from typing import List, Tuple

# Qiskit & Qiskit Nature
from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.problems import ElectronicStructureProblem
from qiskit_nature.second_q.mappers import ParityMapper
from qiskit_nature.second_q.circuit.library import HartreeFock, UCCSD
from qiskit_nature.second_q.transformers import ActiveSpaceTransformer

# IBM Quantum runtime
from qiskit_ibm_runtime import QiskitRuntimeService, Session

try:
    from qiskit_ibm_runtime import Estimator as EstimatorV2
except ImportError:
    from qiskit.primitives import Estimator as EstimatorV2

# PySCF
from pyscf import gto, scf
from scipy.optimize import minimize
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager


# --------------------------------------------------
# 1) 读取 config (TOKEN, INSTANCE)
# --------------------------------------------------
def read_config(file_path):
    config = {}
    with open(file_path, "r") as f:
        for line in f:
            if "=" in line:
                key, val = line.strip().split("=")
                config[key.strip()] = val.strip()
    return config


# --------------------------------------------------
# 2) 解析 PLIP 文件, 获取 (resName, chainID, resNum) + (ligResName, ligChain, ligNum)
# --------------------------------------------------
def parse_plip_interaction_file(plip_path: str):
    """
    根据 1c5z_interaction.txt 的文本格式，示例解析：
    - 找到**蛋白残基**行，如 "| 190   | SER     | B | ... "
    - 找到**配体**信息，如 "MOL:A:1 (MOL) - SMALLMOLECULE"

    返回:
      residue_list: List[(resName, chainID, resNum)]
      ligand_info:  (ligResName, ligChainID, ligResNum)
    """
    residue_list = []
    ligand_info = None

    with open(plip_path, "r") as f:
        content = f.read()

    # 1) 从 "MOL:A:1 (MOL) - SMALLMOLECULE" 抓取 (MOL, A, 1)
    #   你可以用更精确正则. 此处示例：
    lig_pattern = re.compile(r"(MOL):([A-Z]):(\d+)\s+\(MOL\)")
    match_lig = lig_pattern.search(content)
    if match_lig:
        ligName = match_lig.group(1)  # "MOL"
        ligChain = match_lig.group(2)  # "A"
        ligNum = match_lig.group(3)  # "1"
        ligand_info = (ligName, ligChain, ligNum)

    # 2) 找蛋白-配体相互作用行(如 "| 190   | SER     | B        | 1 ")
    #   简单正则示例, 注意可能匹配多行
    #   这里可根据 "RESNR" "RESTYPE" "RESCHAIN" 这几列解析
    #   假设行类似: "| 190   | SER     | B        | 1   | MOL | A ...."
    residue_pattern = re.compile(
        r"\|\s+(\d+)\s+\|\s+([A-Z]{3})\s+\|\s+([A-Z])\s+\|\s+1\s+\|\s+MOL\s+\|\s+A.*"
    )
    # 说明：这个正则假设 "1" 是 RESNR_LIG, "MOL" 是 RESTYPE_LIG, "A" 是 RESCHAIN_LIG,
    # 你需根据真实文件表格列对齐.
    all_res = residue_pattern.findall(content)
    # each match is like ("190","SER","B") ...
    for (resnum, restype, chainid) in all_res:
        residue_list.append((restype, chainid, resnum))

    # 如果上面没匹配到，你可以再写一个单独的解析氢键表格/疏水表格
    # 这里做示例.

    return residue_list, ligand_info


# --------------------------------------------------
# 3) 读取 PDB 并存为 PySCF 的 Mole (不截断任何原子)
# --------------------------------------------------
def build_full_mol_from_pdb(pdb_path: str, charge=0, spin=0, basis="sto3g"):
    """
    读取 PDB 全原子 (ATOM/HETATM), 构造 PySCF Mole 对象
    """
    coords = []
    at_symbols = []

    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_element = line[76:78].strip()
                if not atom_element:
                    atom_element = line[12:16].strip(" 1234567890")
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                at_symbols.append(atom_element)
                coords.append((x, y, z))

    mol = gto.Mole()
    mol.build(
        atom=[(at_symbols[i], coords[i]) for i in range(len(coords))],
        basis=basis,
        charge=charge,
        spin=spin
    )
    return mol


# --------------------------------------------------
# 4) 在 PySCF Mole 中标记“关键区”原子, 并对MO做投影, 决定活性空间
# --------------------------------------------------
def find_active_space_by_projection(mol, mf,
                                    residue_list: List[Tuple[str, str, str]],
                                    ligand_info: Tuple[str, str, str],
                                    pdb_path: str,
                                    threshold=0.2):
    """
    1. residue_list, ligand_info: 来自 parse_plip_interaction_file()
       用于判断哪些(原子index)属于关键残基或配体
    2. 读取 pdb_path 再次，将每个atom与PySCF中mol.atom_index对应, 标记是否关键区
    3. 对于每个 MO, 计算其在关键区上的贡献:
         ratio = sum( |C_{i,MO}|^2 for i in AO_in_key_region )
         / sum( |C_{all,MO}|^2 )
       其中 C_{i,MO} 是 MO系数, i遍历 AO
    4. 如果 ratio >= threshold, 认为该MO主要分布在关键区.
       - 若是占据轨道(occ>0.1), 就选入活性空间 => active_electrons += 2 (或1)
       - 若是虚轨道(occ<0.1), 就选入活性空间 => active_orbitals_count += 1
    最后返回 ( mo_occ_list, mo_index_list ) 给 Qiskit Nature's ActiveSpaceTransformer.
    """

    # (A) 标记PDB中的所有原子 => 是否关键区
    #     需要将 PDB 读取成同顺序(与 PySCF )?
    #     但实际 PySCF 没有保存"resName, chainID" 等信息.
    #     我们只能比对 (x,y,z) 近似方式, 或
    #     先存下 mol.atom_coords() 与 pdb coords, 做最近邻匹配.
    #     这里做一个简单近似匹配(可能不够严格).
    pdb_coords = []
    pdb_key_atom = []
    residue_set = set(residue_list)  # {(SER, B, 190), ...}
    lig_resname, lig_chain, lig_num = ligand_info

    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                resName = line[17:20].strip()
                chainID = line[21].strip()
                resNum = line[22:26].strip()

                # 是否关键
                is_residue = (resName, chainID, resNum) in residue_set
                is_ligand = (resName == lig_resname and chainID == lig_chain and resNum == lig_num)

                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                pdb_coords.append((x, y, z))
                pdb_key_atom.append(is_residue or is_ligand)

    # PySCF Mole 原子顺序 => mol.atom_coord(i)
    # 可能顺序不一定与 PDB 一致.
    # 这里做一个"坐标最近邻"匹配(有时对于大体系可能1~2位小数的差).
    coords_mf = mf.mol.atom_coords()

    # 建立 atom_index => boolean (是否关键)
    key_atom_idx = set()
    used_pdb_idx = set()

    for i_mol in range(mf.mol.natm):
        cx, cy, cz = coords_mf[i_mol]
        # 找到 pdb_coords 中最近点
        min_dist = 9999
        min_j = -1
        for j, (px, py, pz) in enumerate(pdb_coords):
            dx = px - cx
            dy = py - cy
            dz = pz - cz
            dist = dx * dx + dy * dy + dz * dz
            if dist < min_dist:
                min_dist = dist
                min_j = j
        # if that min_j is a "key" atom => i_mol is key
        if pdb_key_atom[min_j]:
            key_atom_idx.add(i_mol)
        # 标记这个pdb原子被用了
        used_pdb_idx.add(min_j)

    nmo = mf.mo_coeff.shape[1]
    mo_occ = mf.mo_occ
    mo_coeff = mf.mo_coeff

    # (B) 构造 AO->原子index 映射
    # PySCF 的 basis 函数是按原子分块的
    ao_loc = mf.mol.ao_loc_nr()  # array of length natm+1
    # ao_loc[i] 到 ao_loc[i+1] 对应 第i个原子的所有AO
    # i 即原子index

    # (C) 计算 MO 在关键区上的贡献
    chosen_mo_idx = []
    chosen_occ = []  # 占据度

    for mo_i in range(nmo):
        # mo_coeff[:, mo_i] => AO基函数到 MO_i 的系数
        # sum_of_squares_all = sum( c^2 )
        # sum_of_squares_key = sum( c^2 for c in AO in key atoms )
        cvec = mo_coeff[:, mo_i]
        sum_sq_all = 0.0
        sum_sq_key = 0.0
        for at_i in range(mf.mol.natm):
            # AO范围
            ao_start = ao_loc[at_i]
            ao_end = ao_loc[at_i + 1]
            c_sub = cvec[ao_start:ao_end]
            # sum of squares
            sum_sub = np.sum(c_sub * c_sub)
            sum_sq_all += sum_sub
            if at_i in key_atom_idx:
                sum_sq_key += sum_sub

        ratio = sum_sq_key / sum_sq_all if sum_sq_all > 1e-15 else 0.0
        # 如果 ratio >= threshold => 认为此 MO 主要分布在关键区
        if ratio >= threshold:
            chosen_mo_idx.append(mo_i)
            chosen_occ.append(mo_occ[mo_i])

    # (D) 计算活性电子数 & 轨道数
    #   对 chosen_mo_idx 里那些 occ>~0.1 => 累计电子
    active_e = 0
    for occ_val in chosen_occ:
        if occ_val > 1.9:
            active_e += 2
        elif occ_val > 0.1:
            # 半占(ROHF?)
            active_e += 1

    # "活性轨道数" 就是 len(chosen_mo_idx)
    active_o = len(chosen_mo_idx)

    print(f"Projection-based MO selection: found {active_o} MOs with ratio>={threshold}, total active_e={active_e}")

    # (E) 返回  (list_of_mo_indices, active_e)
    #  供 Qiskit's ActiveSpaceTransformer 使用:  set_active_space()
    #  注意: ActiveSpaceTransformer 默认用 (num_electrons, num_spatial_orbitals),
    #  也支持 .transform( problem, active_orbitals=..., active_window=(start,end)... )
    #  但不能直接给它不连续的MO indices.
    #  如果 chosen_mo_idx 不是连续区间, 需要再实现 "fragmented active space" 方式.
    #  目前 Qiskit Nature 仅支持**连续**的MO子区间.
    #  => 你可能需要取 min_idx~max_idx 的区间.
    #  这就是ActiveSpaceTransformer的局限.
    #  如果你需要更灵活, 可能要自己构造二次量子算符.

    min_idx = min(chosen_mo_idx) if chosen_mo_idx else 0
    max_idx = max(chosen_mo_idx) if chosen_mo_idx else 1
    num_selected = max_idx - min_idx + 1
    # 这可能比 chosen_mo_idx 多, 如果中间有一些MO ratio<threshold

    # 这里仅演示.
    # 最简单: 取 [min_idx, max_idx+1) 作为活性空间 =>  contiguous range
    return active_e, num_selected, min_idx


# --------------------------------------------------
# 5) 自定义 VQE 类
# --------------------------------------------------
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


# --------------------------------------------------
# 6) 主流程
# --------------------------------------------------
def main():
    # 1) 读取config(TOKEN, INSTANCE)
    config = read_config("config.txt")
    service = QiskitRuntimeService(
        channel='ibm_quantum',
        instance=config["INSTANCE"],
        token=config["TOKEN"]
    )

    # 2) parse PLIP => residue_list, ligand_info
    plip_file = "data_set/data/2_benchmark_binidng_sites/1c5z/1c5z_interaction.txt"
    residue_list, ligand_info = parse_plip_interaction_file(plip_file)
    print("Residues from PLIP:", residue_list)
    print("Ligand from PLIP:", ligand_info)

    # 3) 构造 PySCF Mole(完整体系)
    pdb_path = "data_set/data/2_benchmark_binidng_sites/1c5z/1c5z_Binding_mode.pdb"
    total_charge = 0
    spin = 1
    mol = build_full_mol_from_pdb(pdb_path, charge=total_charge, spin=spin, basis="sto3g")

    # 4) 全体系SCF
    mf = scf.RHF(mol) if spin == 0 else scf.ROHF(mol)
    mf.kernel()

    # 5) 基于投影 => 选出 "关键区" MOs
    #    目前 QiskitNature 只能处理 contiguous range
    active_e, n_orb_contig, mo_start = find_active_space_by_projection(
        mol, mf, residue_list, ligand_info, pdb_path,
        threshold=0.2
    )

    mo_end = mo_start + n_orb_contig
    print(f"ActiveSpace => electrons={active_e}, mo_range=[{mo_start}, {mo_end})")

    # 6) 构造 QiskitNature Problem
    driver = PySCFDriver(
        atom=mol.atom,
        basis=mol.basis,
        charge=mol.charge,
        spin=mol.spin,
        unit=DistanceUnit.ANGSTROM
    )
    es_problem = driver.run()

    # 7) 用 ActiveSpaceTransformer, 但要用 "active_window" 形式
    ast = ActiveSpaceTransformer(
        num_electrons=active_e,
        num_spatial_orbitals=n_orb_contig,
        # 仅支持 contiguous window.
        # QiskitNature 0.6+ : param active_window=(mo_start, mo_end)
        active_window=(mo_start, mo_end)
    )
    red_problem = ast.transform(es_problem)

    # 8) 构造量子比特算符
    op = red_problem.hamiltonian.second_q_op()
    mapper = ParityMapper()
    qubit_op = mapper.map(op)

    # 9) UCCSD ansatz
    n_so = red_problem.num_spatial_orbitals
    alpha = red_problem.num_alpha
    beta = red_problem.num_beta
    hf_init = HartreeFock(n_so, (alpha, beta), mapper)
    ansatz = UCCSD(
        num_spatial_orbitals=n_so,
        num_particles=(alpha, beta),
        mapper=mapper,
        initial_state=hf_init
    )

    # 10) VQE
    vqe = CustomVQE(
        service=service,
        qubit_op=qubit_op,
        ansatz=ansatz,
        shots=100,
        min_qubit_num=30,
        maxiter=15,
        optimization_level=3
    )
    energies, opt_params = vqe.run_vqe()
    final_e = energies[-1]
    print("Final active-space energy:", final_e)

    # 11) 保存
    os.makedirs("results_projection", exist_ok=True)
    with open("results_projection/energy.csv", "w", newline="") as cf:
        writer = csv.writer(cf)
        writer.writerow(["Iter", "Energy"])
        for i, e in enumerate(energies):
            writer.writerow([i + 1, e])
    with open("results_projection/params.json", "w") as jf:
        json.dump({"opt_params": opt_params.tolist()}, jf, indent=4)

    print("\nAll done. Results in 'results_projection/'.")


if __name__ == "__main__":
    main()



