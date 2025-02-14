# --*-- conding:utf-8 --*--
# @Time : 2/14/25 4:43 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : Main.py

from utils import PDBStructure

if __name__ == "__main__":
    # 示例使用
    pdb_path = "data_set/EC5026_5Apart.pdb"  # 放你的pdb路径
    pdb_struct = PDBStructure(pdb_path)

    print(f"总共读取到 {len(pdb_struct.atoms)} 个原子记录.")
    print(f"连接信息(CONECT) keys: {list(pdb_struct.conect_info.keys())[:10]} ...")

    # 按残基分组
    by_res = pdb_struct.group_by_residue()
    print("Residues found:", list(by_res.keys())[:5], "...")

    # 找某个配体resName, 例如 '2RV'
    ligand_atoms = pdb_struct.get_atoms(record_type="HETATM", resName="2RV")
    print(f"配体2RV共 {len(ligand_atoms)} 个原子.")
