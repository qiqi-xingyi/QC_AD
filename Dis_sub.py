# --*-- conding:utf-8 --*--
# @Time : 2/15/25 4:25 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : Dis_sub.py

import Bio.PDB
import os

def extract_reactive_fragments(input_pdb, output_pdb, ligand_resname, distance_cutoff=2.4):
    parser = Bio.PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_pdb)
    io = Bio.PDB.PDBIO()
    io.set_structure(structure)

    ligand_atoms = []
    selected_atoms = set()

    # 找到配体原子
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() == ligand_resname:
                    ligand_atoms.extend(residue.get_atoms())

    # 仅保留配体相互作用的氨基酸原子（而不是整个氨基酸）
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    for ligand_atom in ligand_atoms:
                        if atom - ligand_atom < distance_cutoff:
                            # 仅保留有作用的原子
                            if atom.get_name() not in ["C", "O", "CA", "N"]:  # 去掉主链
                                selected_atoms.add(atom)
                            break

    # 保存新 PDB 文件，仅包含相互作用的原子
    class SelectAtoms(Bio.PDB.Select):
        def accept_atom(self, atom):
            return atom in selected_atoms

    io.save(output_pdb, select=SelectAtoms())
    print(f"Extracted reactive fragments saved to {output_pdb}")

    # 转换为 XYZ 文件
    output_xyz = os.path.splitext(output_pdb)[0] + ".xyz"
    convert_pdb_to_xyz(output_pdb, output_xyz)
    print(f"Converted PDB to XYZ format: {output_xyz}")

def convert_pdb_to_xyz(input_pdb, output_xyz):
    with open(input_pdb, "r") as pdb_file, open(output_xyz, "w") as xyz_file:
        atoms = []
        for line in pdb_file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                element = line[76:78].strip()
                if not element:
                    element = line[12:14].strip().upper()
                if element in {"C", "H", "O", "N", "S", "P", "F", "Cl", "Br", "I"}:
                    x, y, z = line[30:38].strip(), line[38:46].strip(), line[46:54].strip()
                    atoms.append(f"{element} {x} {y} {z}")
        xyz_file.write(f"{len(atoms)}\n\n")
        xyz_file.write("\n".join(atoms) + "\n")

if __name__ == '__main__':
    input_pdb = "data_set/EC5026_5Apart.pdb"
    output_pdb = "2RV_reactive_fragments.pdb"
    ligand_resname = "2RV"
    extract_reactive_fragments(input_pdb, output_pdb, ligand_resname)


