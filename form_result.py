# --*-- conding:utf-8 --*--
# @Time : 2/15/25 10:47 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : form_result.py

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

if __name__ == '__main__':


    # 假设我们已经分别计算出这三者的最低能量：
    E_complex = -75.34  # 复合体系的最低能量
    E_protein = -50.12  # 游离蛋白的最低能量
    E_ligand = -25.01  # 游离配体的最低能量

    # 计算结合能
    E_binding = compute_binding_energy(E_complex, E_protein, E_ligand)

    print(f"结合能 (Binding Energy): {E_binding:.4f} Hartree")

    # 示例：假设 VQE 计算出的最低能量是 -75.34 Hartree
    energy_hartree = -75.34
    converted_energy = hartree_to_units(energy_hartree)

    print(f"Energy in eV: {converted_energy['eV']:.4f} eV")
    print(f"Energy in kcal/mol: {converted_energy['kcal/mol']:.4f} kcal/mol")
    print(f"Energy in kJ/mol: {converted_energy['kJ/mol']:.4f} kJ/mol")
