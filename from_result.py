# --*-- conding:utf-8 --*--
# @Time : 2/15/25 10:47 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : from_result.py

import datetime

def hartree_to_units(energy_hartree):
    """
    Converts energy from Hartree to eV, kcal/mol, and kJ/mol.

    Parameters:
    - energy_hartree: Energy in Hartree

    Returns:
    - Dictionary with converted energy values
    """
    hartree_to_ev = 27.2114  # 1 Hartree = 27.2114 eV
    hartree_to_kcalmol = 627.509  # 1 Hartree = 627.509 kcal/mol
    hartree_to_kjmol = 2625.5  # 1 Hartree = 2625.5 kJ/mol

    return {
        "eV": energy_hartree * hartree_to_ev,
        "kcal/mol": energy_hartree * hartree_to_kcalmol,
        "kJ/mol": energy_hartree * hartree_to_kjmol
    }

def compute_binding_energy(E_complex, E_protein, E_ligand):
    """
    计算结合能 (Binding Energy).

    参数:
    - E_complex: 复合体系的最低能量
    - E_protein: 游离蛋白的最低能量
    - E_ligand: 游离配体的最低能量

    返回:
    - E_binding: 结合能 (负值表示有利结合)
    """
    E_binding = E_complex - (E_protein + E_ligand)
    return E_binding


def save_results(results, filename=None):
    """
    将计算结果保存到文本文件

    参数:
    - results: 包含所有计算结果的字典
    - filename: 自定义文件名 (可选)
    """
    # 生成带时间戳的默认文件名
    if not filename:
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"binding_energy_results_{timestamp}.txt"

    # 构建文件内容
    content = f"""\
================================
 Binding Energy Calculation Report
================================

Calculation Time: {results['timestamp']}

Energies from VQE(Hartree):
- Complex System: {results['E_complex']:.4f}
- Protein Alone: {results['E_protein']:.4f}
- Ligand Alone:  {results['E_ligand']:.4f}

Results:
[Binding Energy]
Hartree:    {results['E_binding']:.4f}
eV:        {results['converted']['eV']:.4f}
kcal/mol:  {results['converted']['kcal/mol']:.4f}
kJ/mol:    {results['converted']['kJ/mol']:.4f}

"""

    # 写入文件
    with open(filename, 'w') as f:
        f.write(content)

    return filename

if __name__ == '__main__':

    E_complex = -5.6799  # The minimum energy of a composite system
    E_protein = -1.1267  # Minimum energy of free protein
    E_ligand = -4.4925  # Minimum energy of free ligand

    # Calculate binding energy
    E_binding = compute_binding_energy(E_complex, E_protein, E_ligand)
    converted = hartree_to_units(E_binding)

    print(f"Binding Energy: {E_binding:.4f} Hartree")

    energy_hartree = E_binding
    converted_energy = hartree_to_units(energy_hartree)

    print(f"Energy in eV: {converted_energy['eV']:.4f} eV")
    print(f"Energy in kcal/mol: {converted_energy['kcal/mol']:.4f} kcal/mol")
    print(f"Energy in kJ/mol: {converted_energy['kJ/mol']:.4f} kJ/mol")

    results = {
        'timestamp': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'E_complex': E_complex,
        'E_protein': E_protein,
        'E_ligand': E_ligand,
        'E_binding': E_binding,
        'converted': converted
    }

    saved_file = save_results(results)
    print(f"\nResults saved to: {saved_file}")