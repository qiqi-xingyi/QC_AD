# --*-- conding:utf-8 --*--
# @Time : 2/15/25 4:25 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : Dis_sub.py

# !/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PDB子系统提取工具
用法：python extract_subsystem.py input.pdb [options]
"""

# --*-- coding:utf-8 --*--
from Bio.PDB import PDBParser, Select, PDBIO

# ========== 调试参数设置区 ==========
INPUT_PDB = "./data_set/EC5026_5Apart.pdb"  # 输入文件路径
LIGAND_ID = "B:603"  # 配体位置 (格式: 链:残基号)
CUTOFF = 4.0  # 选择半径 (Å)
OUTPUT_PDB = "subsystem.pdb"  # 输出文件名


# ===================================

class SubsystemSelector(Select):
    """ 基于距离的子系统选择器 """

    def __init__(self):
        self.ligand_chain, self.ligand_resnum = LIGAND_ID.split(":")
        self.ligand_resnum = int(self.ligand_resnum)
        self.cutoff_sq = CUTOFF ** 2
        self.ligand_atoms = []
        self.selected_residues = set()

    def process_ligand(self, residue):
        """ 识别并记录配体原子坐标 """
        if (residue.parent.id == self.ligand_chain
                and residue.id[1] == self.ligand_resnum):
            self.ligand_atoms.extend(atom.get_coord() for atom in residue)
            return True
        return False

    def accept_residue(self, residue):
        # 自动包含配体本身
        if self.process_ligand(residue):
            return True

        # 检查蛋白质残基
        for atom in residue:
            for lig_coord in self.ligand_atoms:
                distance_sq = sum(
                    (a - b) ** 2
                    for a, b in zip(atom.get_coord(), lig_coord)
                )
                if distance_sq < self.cutoff_sq:
                    self.selected_residues.add(residue)
                    return True
        return False

    def accept_atom(self, atom):
        return atom.get_parent() in self.selected_residues


if __name__ == '__main__':
    try:
        # 初始化选择器
        selector = SubsystemSelector()

        # 解析结构
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('input', INPUT_PDB)

        # 验证配体存在
        if not selector.ligand_atoms:
            raise ValueError(
                f"❌ 配体 {LIGAND_ID} 未找到！请检查：\n"
                f"1. 输入文件是否包含链 {LIGAND_ID.split(':')[0]}\n"
                f"2. 是否存在残基号 {LIGAND_ID.split(':')[1]}"
            )

        # 保存子系统
        io = PDBIO()
        io.set_structure(structure)
        io.save(OUTPUT_PDB, selector)

        # 打印调试信息
        print("✅ 调试结果")
        print(f"输入文件: {INPUT_PDB}")
        print(f"找到配体: {LIGAND_ID}")
        print(f"包含邻近残基数: {len(selector.selected_residues)}")
        print(f"输出文件: {OUTPUT_PDB}")

        # 可视化提示
        print("\n🔍 建议使用PyMOL验证结果:")
        print(f"cmd.load('{OUTPUT_PDB}')")
        print("cmd.show('sticks', 'all')")
        print("cmd.zoom()")

    except FileNotFoundError:
        print(f"🔥 文件未找到: {INPUT_PDB}\n"
              f"请检查路径是否正确，当前工作目录: {os.getcwd()}")
    except Exception as e:
        print(f"❌ 发生错误: {str(e)}\n"
              "💡 调试建议:")
        print("1. 检查配体标识符格式 (例如 B:603)")
        print("2. 用文本编辑器打开PDB文件验证配体是否存在")
        print("3. 尝试减小CUTOFF值")