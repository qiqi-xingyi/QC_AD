# --*-- conding:utf-8 --*--
# @Time : 2/15/25 4:12 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : Amber.py

"""
基于MM/PBSA分解能量自动选择关键残基
用法：python select_subsystem.py decomp.dat [-c CUTOFF] [-o OUTPUT]
"""

import argparse
import pandas as pd
from Bio.PDB import PDBParser, PDBIO, Select


class EnergyBasedSelector(Select):
    """ 基于能量贡献的残基选择器 """

    def __init__(self, energy_data, cutoff=-1.0):
        self.key_residues = set(
            (row['Residue'], row['Chain'])
            for _, row in energy_data.iterrows()
            if row['Total'] < cutoff
        )

    def accept_residue(self, residue):
        res_id = residue.get_id()[1]
        chain = residue.get_parent().id
        return (str(res_id), chain) in self.key_residues


def parse_decomp_file(filename):
    """ 解析AMBER分解能量文件 """
    try:
        # 自动检测文件格式
        with open(filename) as f:
            first_line = f.readline().strip()

        if first_line.startswith('Residue'):
            # 新版MMPBSA格式
            df = pd.read_csv(filename, delim_whitespace=True, skiprows=1)
            df.columns = [c.strip() for c in df.columns]
        else:
            # 旧版固定宽度格式
            df = pd.read_fwf(
                filename,
                colspecs=[(0, 7), (8, 13), (14, 20), (21, 28), (29, 36), (37, 44)],
                names=['Residue', 'Chain', 'Van der Waals', 'Electrostatic', 'Polar', 'Non-Polar']
            )
            df['Total'] = df.sum(axis=1)

        return df
    except Exception as e:
        raise ValueError(f"文件解析失败: {str(e)}")


def main():
    parser = argparse.ArgumentParser(description='基于能量贡献选择子体系')
    parser.add_argument('input', help='MMPBSA分解结果文件')
    parser.add_argument('-c', '--cutoff', type=float, default=-1.0,
                        help='能量阈值(kcal/mol)，默认-1.0')
    parser.add_argument('-o', '--output', default='./subsystem/subsystem.pdb',
                        help='输出PDB文件名')
    parser.add_argument('-s', '--structure', default='./data_set/EC5026_5Apart.pdb',
                        help='原始PDB结构文件')
    args = parser.parse_args()

    # 步骤1：解析能量分解文件
    try:
        energy_data = parse_decomp_file(args.input)
        print(f"成功解析能量数据，共{len(energy_data)}个残基")
    except Exception as e:
        print(f"错误: {str(e)}")
        return

    # 步骤2：选择关键残基
    selector = EnergyBasedSelector(energy_data, args.cutoff)
    selected_res = len(selector.key_residues)
    print(f"选择{selected_res}个关键残基 (ΔG < {args.cutoff} kcal/mol)")

    # 步骤3：从原始结构提取子体系
    try:
        parser = PDBParser()
        structure = parser.get_structure('original', args.structure)

        io = PDBIO()
        io.set_structure(structure)
        io.save(args.output, selector)

        print(f"子体系已保存至 {args.output}")
        print("\nPyMOL查看命令:")
        print(f"load {args.output}")
        print(f"load {args.structure}, original")
        print("align subsystem, original")
    except Exception as e:
        print(f"结构处理错误: {str(e)}")


if __name__ == '__main__':
    main()