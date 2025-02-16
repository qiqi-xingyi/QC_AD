# --*-- conding:utf-8 --*--
# @Time : 2/14/25 9:08 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : Quantum_dock.py

from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import ParityMapper
from qiskit_nature.second_q.transformers import ActiveSpaceTransformer
# from qiskit_nature.second_q.circuit.library import UCCSD
from vqe import VQE
from qiskit_ibm_runtime import QiskitRuntimeService
import json
import csv


def xyz_to_string(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    atom_lines = lines[2:]
    result = "; ".join(line.strip() for line in atom_lines)

    return result

def read_config(file_path):

    config = {}
    try:
        with open(file_path, "r") as file:
            for line in file:
                key, value = line.strip().split("=")
                config[key.strip()] = value.strip()
    except Exception as e:
        print(f"Fail to read token file: {e}")
        return None
    return config

def save_results(energy_list, optimized_params, energy_filename="energy_results.csv", params_filename="optimized_params.json"):
    """
    Saves the VQE energy list and optimized parameters.

    Parameters:
    - energy_list: List of computed energy values during optimization.
    - optimized_params: Optimized parameters from VQE.
    - energy_filename: Filename for saving energy values (CSV).
    - params_filename: Filename for saving parameters (JSON).
    """
    # Save energy list as CSV
    with open(energy_filename, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Iteration", "Energy"])
        for i, energy in enumerate(energy_list):
            writer.writerow([i + 1, energy])
    print(f"Saved energy results to {energy_filename}")

    # Save optimized parameters as JSON
    params_dict = {"optimized_params": optimized_params.tolist()}
    with open(params_filename, "w") as jsonfile:
        json.dump(params_dict, jsonfile, indent=4)
    print(f"Saved optimized parameters to {params_filename}")

if __name__=="__main__":

    file_path = "2RV_subsyetem.xyz"
    converted_string = xyz_to_string(file_path)
    print(converted_string)

    driver = PySCFDriver(atom=converted_string, charge=0, spin=0, basis="sto3g")
    problem = driver.run()
    print('result:', problem)

    transformer = ActiveSpaceTransformer(num_electrons=10, num_spatial_orbitals=8)
    reduced_problem = transformer.transform(problem)

    print('reduced_problem:', reduced_problem)

    op = reduced_problem.hamiltonian.second_q_op()

    mapper = ParityMapper()
    qubit_op = mapper.map(op)

    print("Qubit Hamiltonian Terms:", len(qubit_op))
    print("Qubit Num:", qubit_op.num_qubits)

############

    config_path = "config.txt"

    config = read_config(config_path)

    service = QiskitRuntimeService(
        channel='ibm_quantum',
        instance=config["INSTANCE"],
        token=config["TOKEN"]
    )

    qubits = qubit_op.num_qubits + 3

    print("Qubit Num Chosen:", qubits)

    vqe = VQE(service=service, hamiltonian=qubit_op, min_qubit_num=qubits, shots=50, maxiter=30)

    ene_list, ground_state = vqe.run_vqe()








