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
import os


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

def save_results(energy_list, optimized_params, save_dir="results",
                 energy_filename="energy_results.csv", params_filename="optimized_params.json"):
    """
    Saves the VQE energy list and optimized parameters with improvements:
    1. Allows specifying the save directory.
    2. Finds the lowest energy and writes it at the beginning.
    3. Creates the directory if it does not exist.

    Parameters:
    - energy_list: List of computed energy values during optimization.
    - optimized_params: Optimized parameters from VQE.
    - save_dir: Directory to save the files.
    - energy_filename: Filename for saving energy values (CSV).
    - params_filename: Filename for saving parameters (JSON).
    """

    # Ensure the save directory exists
    os.makedirs(save_dir, exist_ok=True)

    # Get the minimum energy
    min_energy = min(energy_list)
    min_index = energy_list.index(min_energy) + 1  # Iteration starts from 1

    # Save energy list as CSV (with lowest energy at the top)
    energy_path = os.path.join(save_dir, energy_filename)
    with open(energy_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Iteration", "Energy"])
        writer.writerow(["Lowest Energy", min_energy])  # Write the lowest energy first
        writer.writerow([])  # Empty row for clarity
        for i, energy in enumerate(energy_list):
            writer.writerow([i + 1, energy])
    print(f"Saved energy results to {energy_path}")

    # Save optimized parameters as JSON
    params_dict = {"optimized_params": optimized_params.tolist()}
    params_path = os.path.join(save_dir, params_filename)
    with open(params_path, "w") as jsonfile:
        json.dump(params_dict, jsonfile, indent=4)
    print(f"Saved optimized parameters to {params_path}")

if __name__=="__main__":

    file_path = "subsystem/2RV_reactive_fragments.xyz"
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

    vqe = VQE(service=service, hamiltonian=qubit_op, min_qubit_num=qubits, shots=50, maxiter=20)

    ene_list, ground_state = vqe.run_vqe()

    save_results(energy_list=ene_list, optimized_params=ground_state, save_dir="vqe_results", energy_filename="energy_results_1.csv", params_filename="optimized_params_1.json")











