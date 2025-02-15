# --*-- conding:utf-8 --*--
# @Time : 2/14/25 9:08 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : Quantum_dock.py
from networkx.algorithms.tournament import hamiltonian_path
from pyscf import gto, scf
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import ParityMapper
from qiskit_nature.second_q.transformers import ActiveSpaceTransformer


def xyz_to_string(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    atom_lines = lines[2:]

    result = "; ".join(line.strip() for line in atom_lines)

    return result

if __name__=="__main__":

    # file_path = "subsystems_xyz/sub_A_267_PHE/res.xyz"  # 请替换为你的xyz文件路径
    # converted_string = xyz_to_string(file_path)
    # print(converted_string)


    # 定义苯丙氨酸的坐标
    phe_geometry = """N  22.705  13.849  13.767
                    C  23.996  14.373  14.234
                    C  24.354  15.649  13.450
                    O  25.018  15.565  12.415
                    C  24.929  13.195  13.947
                    C  24.387  12.551  12.678
                    C  22.909  12.937  12.641
                    H  24.048  14.576  15.297
                    H  25.979  13.472  13.858
                    H  24.851  12.490  14.774
                    H  24.888  12.964  11.805
                    H  24.540  11.472  12.673
                    H  22.652  13.424  11.700
                    H  22.291  12.045  12.748
                    C  17.293  10.433  14.307
                    C  9.912  10.241  19.781
                    C  9.055  11.287  19.048
                    N  16.708  11.169  15.351
                    C  16.047  10.677  16.417
                    N  15.617  11.640  17.243
                    C  14.530  11.527  18.220
                    C  14.631  10.395  19.280
                    C  13.441  10.329  20.273
                    N  12.918  11.684  20.569
                    C  11.595  11.919  20.797
                    C  10.627  10.733  21.072
                    C  9.610  11.063  22.192
                    O  11.146  13.065  20.892
                    C  13.831  12.818  20.343
                    C  14.107  12.886  18.828
                    O  15.853  9.485  16.629
                    C  16.781  9.216  13.837
                    C  17.447  8.503  12.843
                    C  18.631  9.003  12.301
                    O  19.303  8.274  11.349
                    C  20.269  7.359  11.863
                    C  19.113  10.245  12.721
                    C  18.440  10.962  13.713
                    F  20.247  10.745  12.186
                    F  20.775  6.697  10.826
                    F  19.720  6.464  12.689
                    F  21.271  7.992  12.477
                    H  16.884  12.166  15.345
                    H  16.082  12.537  17.202
                    H  9.279  9.389  20.032
                    H  10.641  9.842  19.075
                    H  8.574  10.848  18.176
                    H  8.269  11.681  19.690
                    H  13.678  11.205  17.616
                    H  9.657  12.127  18.702
                    H  11.208  9.897  21.463
                    H  15.867  8.809  14.236
                    H  17.054  7.556  12.507
                    H  18.831  11.910  14.047
                    H  15.531  10.559  19.867
                    H  14.754  9.410  18.834
                    H  13.705  9.802  21.191
                    H  12.681  9.762  19.737
                    H  8.974  10.206  22.416
                    H  10.120  11.334  23.118
                    H  8.957  11.896  21.931
                    H  14.749  12.585  20.882
                    H  13.495  13.753  20.790
                    H  13.226  13.280  18.317
                    H  14.901  13.617  18.692"""

    # 运行 PySCF 计算
    driver = PySCFDriver(atom=phe_geometry, charge=0, spin=0, basis="sto3g")
    problem = driver.run()
    print('result:', problem)

    transformer = ActiveSpaceTransformer(num_electrons=12, num_spatial_orbitals=12)
    reduced_problem = transformer.transform(problem)

    print('reduced_problem:', reduced_problem)

    hamiltonian = reduced_problem.hamiltonian
    print('hamiltonian:', hamiltonian)

    op = hamiltonian.second_q_op()

    print('op')

    mapper = ParityMapper()

    print('mapper:', mapper)

    qubit_op = mapper.map(op)

    print("Qubit Hamiltonian:", qubit_op)
    print('Qubit Hamiltonian 计算完成')

    print("Qubit Num:", qubit_op.num_qubits)



