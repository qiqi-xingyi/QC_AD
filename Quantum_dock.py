# --*-- conding:utf-8 --*--
# @Time : 2/14/25 9:08 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : Quantum_dock.py

from pyscf import gto, scf
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import ParityMapper
from qiskit_nature.second_q.transformers import ActiveSpaceTransformer, FreezeCoreTransformer


def xyz_to_string(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    atom_lines = lines[2:]
    result = "; ".join(line.strip() for line in atom_lines)

    return result

if __name__=="__main__":

    file_path = "subsystems_xyz/sub_A_267_PHE/res_plus_ligand.xyz"
    converted_string = xyz_to_string(file_path)
    print(converted_string)

    driver = PySCFDriver(atom=converted_string, charge=0, spin=0, basis="sto3g")
    problem = driver.run()
    print('result:', problem)

    # transformer = FreezeCoreTransformer()
    # problem = transformer.transform(problem)

    transformer = ActiveSpaceTransformer(num_electrons=40, num_spatial_orbitals=40)
    reduced_problem = transformer.transform(problem)

    print('reduced_problem:', reduced_problem)

    op = reduced_problem.hamiltonian.second_q_op()

    mapper = ParityMapper()
    qubit_op = mapper.map(op)

    print("Qubit Hamiltonian Terms:", len(qubit_op))
    print("Qubit Num:", qubit_op.num_qubits)




