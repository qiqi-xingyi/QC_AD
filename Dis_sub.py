# --*-- conding:utf-8 --*--
# @Time : 2/15/25 4:25 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : Dis_sub.py

import Bio.PDB


def extract_nearby_residues(input_pdb, output_pdb, ligand_resname, distance_cutoff=2.3):
    parser = Bio.PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_pdb)
    io = Bio.PDB.PDBIO()
    io.set_structure(structure)

    ligand_atoms = []
    protein_atoms = []
    selected_residues = set()

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() == ligand_resname:
                    ligand_atoms.extend(residue.get_atoms())
                else:
                    protein_atoms.extend(residue.get_atoms())

    for protein_atom in protein_atoms:
        for ligand_atom in ligand_atoms:
            if protein_atom - ligand_atom < distance_cutoff:
                selected_residues.add(protein_atom.get_parent())

    class SelectResidues(Bio.PDB.Select):
        def accept_residue(self, residue):
            return residue in selected_residues

    io.save(output_pdb, select=SelectResidues())
    print(f"Extracted system saved to {output_pdb}")


if __name__ == '__main__':

    input_pdb = "data_set/EC5026_5Apart.pdb"  # 你的 PDB 文件路径
    output_pdb = "2RV_subsystem.pdb"  # 提取的 PDB 文件
    ligand_resname = "2RV"  # 配体的名字（根据你的 PDB 文件调整）
    extract_nearby_residues(input_pdb, output_pdb, ligand_resname)
