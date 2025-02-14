# --*-- conding:utf-8 --*--
# @Time : 2/12/25 10:13 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : Molecular_Docking.py

import os
from collections import defaultdict

from qiskit_nature.drivers.second_quantization import ElectronicStructureMoleculeDriver
from qiskit_nature.drivers import Molecule
from qiskit_nature.units import DistanceUnit
from qiskit_nature.problems.second_quantization import ElectronicStructureProblem
from qiskit_nature.transformers.second_quantization.electronic import FreezeCoreTransformer
from qiskit_nature.converters.second_quantization import QubitConverter
from qiskit_nature.mappers.second_quantization import JordanWignerMapper

from qiskit.algorithms import VQE
from qiskit.algorithms.optimizers import COBYLA
from qiskit_nature.circuit.library import UCCSD
from qiskit.utils import QuantumInstance
from qiskit import Aer

# -------------------------------
# 1) 先把之前的PDBStructure类黏贴或import (简化版)
# -------------------------------
class PDBStructure:
    def __init__(self, pdb_file=None):
        self.atoms = []
        self.conect_info = {}
        self.ter_records = []
        self.other_lines = []
        if pdb_file:
            self.read_pdb(pdb_file)

    def read_pdb(self, filepath):
        with open(filepath, 'r') as f:
            for line in f:
                record_type = line[0:6].strip()
                if record_type in ("ATOM","HETATM"):
                    atom_data = self._parse_atom_line(line)
                    self.atoms.append(atom_data)
                elif record_type=="CONECT":
                    self._parse_conect_line(line)
                elif record_type=="TER":
                    self.ter_records.append(line.strip())
                else:
                    self.other_lines.append(line.rstrip("\n"))

    def _parse_atom_line(self, line):
        atom_serial = int(line[6:11].strip())
        atom_name   = line[12:16].strip()
        res_name    = line[17:20].strip()
        chain_id    = line[21].strip()
        res_seq     = line[22:26].strip()
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        element = line[76:78].strip() if len(line)>=78 else ""

        return {
            "record_type": line[0:6].strip(),  # "ATOM"/"HETATM"
            "serial": atom_serial,
            "name": atom_name,
            "resName": res_name,
            "chainID": chain_id,
            "resSeq": res_seq,
            "x": x,
            "y": y,
            "z": z,
            "element": element
        }

    def _parse_conect_line(self, line):
        fields = line.split()
        if len(fields)<2:
            return
        atom_main = int(fields[1])
        bonded_list = []
        for f in fields[2:]:
            try:
                bonded_list.append(int(f))
            except:
                pass

        if atom_main not in self.conect_info:
            self.conect_info[atom_main] = set()
        for b in bonded_list:
            self.conect_info[atom_main].add(b)
            # 可选:对称
            if b not in self.conect_info:
                self.conect_info[b] = set()
            self.conect_info[b].add(atom_main)

    def group_by_residue(self):
        residue_dict = defaultdict(list)
        for atm in self.atoms:
            key = (atm["chainID"], atm["resSeq"], atm["resName"])
            residue_dict[key].append(atm)
        return residue_dict

# -------------------------------
# 2) 定义：把 (symbol,(x,y,z)) 形式对接到 Qiskit 并算能量
# -------------------------------
def compute_energy_vqe(geometry, basis="sto3g"):
    """
    给定 geometry=[(symbol,(x,y,z)), ...], 用Qiskit Nature (PySCF) + VQE 计算基态能量(Hartree).
    """
    molecule = Molecule(
        geometry=geometry,
        charge=0,
        multiplicity=1,
        units=DistanceUnit.ANGSTROM
    )
    driver = ElectronicStructureMoleculeDriver(molecule=molecule, basis=basis, driver_type="PYSCF")
    problem = ElectronicStructureProblem(driver)

    # 可选: 冻结核
    transf = FreezeCoreTransformer()
    problem = transf.transform(problem)

    second_q_ops = problem.second_q_ops()
    main_op = second_q_ops[0]

    num_particles = (problem.molecule_data_transformed.num_alpha,
                     problem.molecule_data_transformed.num_beta)
    num_spin_orbitals = problem.molecule_data_transformed.num_spin_orbitals

    converter = QubitConverter(mapper=JordanWignerMapper(), two_qubit_reduction=False)
    qubit_hamiltonian = converter.convert(main_op, num_particles=num_particles)

    ansatz = UCCSD(
        qubit_converter=converter,
        num_particles=num_particles,
        num_spin_orbitals=num_spin_orbitals
    )
    optimizer = COBYLA(maxiter=100)
    vqe_solver = VQE(ansatz=ansatz, optimizer=optimizer)
    quantum_instance = QuantumInstance(backend=Aer.get_backend("statevector_simulator"), shots=1)
    vqe_solver.quantum_instance = quantum_instance

    result = vqe_solver.compute_minimum_eigenvalue(qubit_hamiltonian)
    energy = result.eigenvalue.real
    return energy

# -------------------------------
# 3) 主流程: 读取pdb -> 选第一个片段 + 配体 -> 计算能量
# -------------------------------
def main():
    PDB_FILE_PATH = "example.pdb"  # 修改为你的实际PDB路径
    pdb_struct = PDBStructure(PDB_FILE_PATH)

    residues = pdb_struct.group_by_residue()
    # 过滤出蛋白氨基酸(Chain A...) vs 配体(2RV ?)
    # 这里简化: 第一个片段(蛋白) => residue_keys[0]
    residue_keys = list(residues.keys())
    residue_keys.sort()  # sort by chain/resSeq/resName

    print("All residue keys:", residue_keys)

    # 取第一个(蛋白)片段
    first_res_key = residue_keys[0]
    first_res_atoms = residues[first_res_key]  # list of dict
    print(f"Selected first residue: {first_res_key}, #atoms={len(first_res_atoms)}")

    # 取配体: 假定resName='2RV' (根据你pdb中的真实配体名称)
    ligand_atoms = []
    for (chain,resSeq,resName), atoms_list in residues.items():
        if resName == "2RV":  # or another check
            ligand_atoms = atoms_list
            break

    if not ligand_atoms:
        print("Warning: Did not find '2RV' in the PDB. Perhaps the PDB has different ligand name?")
        return

    print(f"Selected ligand 2RV: #atoms={len(ligand_atoms)}")

    # 4) 拼接 geometry
    #   Qiskit geometry = [(symbol,(x,y,z)), ...]
    #   from each atom, symbol=atom["element"] or guess from atom["name"]
    def atom_to_xyz(atom):
        symbol = atom["element"] if atom["element"] else atom["name"][0]  # fallback
        return (symbol, (atom["x"], atom["y"], atom["z"]))

    fragment_geometry = [atom_to_xyz(a) for a in first_res_atoms]
    ligand_geometry = [atom_to_xyz(a) for a in ligand_atoms]

    combined_geometry = fragment_geometry + ligand_geometry
    print(f"Combined #atoms = {len(combined_geometry)}")

    # 5) 用 VQE 计算能量
    energy_hartree = compute_energy_vqe(combined_geometry, basis="sto3g")
    print(f"Energy for (first protein fragment + ligand) subsystem = {energy_hartree:.6f} Hartree")

if __name__ == "__main__":
    main()

