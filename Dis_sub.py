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

import argparse
import math
from Bio.PDB import PDBParser, Select, PDBIO


class SubsystemSelector(Select):
    """ 基于距离的子系统选择器 """

    def __init__(self, ligand_id, cutoff=4.0):
        self.ligand_atoms = []
        self.cutoff_sq = cutoff ** 2
        self.ligand_chain, self.ligand_resnum = ligand_id.split(":")
        self.ligand_resnum = int(self.ligand_resnum)
        self.selected_residues = set()

    def process_ligand(self, residue):
        """ 识别并记录配体原子坐标 """
        if (residue.parent.id == self.ligand_chain and
                residue.id[1] == self.ligand_resnum):
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

                distance_sq = sum((a - b) ** 2 for a, b in zip(atom.get_coord(), lig_coord)

                if distance_sq < self.cutoff_sq:

                    self.selected_residues.add(residue)
                    return True

        return False

    def accept_atom(self, atom):
        # 保留所有选中残基的原子
        return atom.get_parent() in self.selected_residues


def main():
    parser = argparse.ArgumentParser(description='PDB子系统提取工具')
    parser.add_argument('input', default='data_set/EC5026_5Apart.pdb', help='输入PDB文件路径')
    parser.add_argument('-l', '--ligand', required=True, default='B:603',
                        help='配体位置 (格式: 链:残基号 如 B:603)')
    parser.add_argument('-c', '--cutoff', type=float, default=3.0,
                        help='选择半径 (Å，默认4.0)')
    parser.add_argument('-o', '--output', default='subsystem.pdb',
                        help='输出文件名 (默认subsystem.pdb)')

    args = parser.parse_args()

    try:
        # 创建选择器
        selector = SubsystemSelector(args.ligand, args.cutoff)

        # 解析结构
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('input', args.input)

        # 验证配体存在
        if not selector.ligand_atoms:
            raise ValueError(f"未找到配体 {args.ligand}，请检查输入")

        # 保存子系统
        io = PDBIO()
        io.set_structure(structure)
        io.save(args.output, selector)

        # 统计信息
        print(f"成功生成子系统文件: {args.output}")
        print(f"包含配体 {args.ligand} 及其 {len(selector.selected_residues)} 个邻近残基")
        print(f"选择半径: {args.cutoff} Å")

    except Exception as e:
        print(f"错误: {str(e)}")
        print("建议检查：")
        print("1. 配体标识符格式是否正确 (例如 B:603)")
        print("2. 输入文件是否包含指定配体")
        print("3. 文件路径是否正确")


if __name__ == '__main__':
    main()