# --*-- conding:utf-8 --*--
# @Time : 2/13/25 8:17 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : PDBStructure.py

import re
from collections import defaultdict


class PDBStructure:
    """
    用于读取并存储PDB文件中的信息，
    包括:
      - ATOM / HETATM 行 (原子坐标, 元素, 残基名, 链ID, 等)
      - CONECT 行 (配体/其他HETATM原子间的键信息)

    可以后续为量子计算的片段化/构造分子坐标等做准备。
    """

    def __init__(self, pdb_file=None):
        self.atoms = []  # 存所有原子(包括ATOM/HETATM)
        self.conect_info = {}  # {atom_serial: [bonded_atom_serials]}
        self.ter_records = []  # 存TER行出现的位置/信息
        self.title = []  # 可存TITLE/HEADER等行 (可选)
        self.other_lines = []  # 存储不处理的其他行

        if pdb_file:
            self.read_pdb(pdb_file)

    def read_pdb(self, filepath):
        """
        读取PDB文件，对ATOM/HETATM/CONECT/TER等感兴趣的记录进行解析。
        """
        with open(filepath, 'r') as f:
            for line in f:
                record_type = line[0:6].strip()  # 例如"ATOM","HETATM","CONECT","TER"

                if record_type in ("ATOM", "HETATM"):
                    atom_data = self._parse_atom_line(line)
                    self.atoms.append(atom_data)

                elif record_type == "CONECT":
                    self._parse_conect_line(line)

                elif record_type == "TER":
                    self.ter_records.append(line.strip())

                else:
                    # 也许保留其他记录(HEADER,TITLE,REMARK,END等)以备后用
                    self.other_lines.append(line.rstrip("\n"))

    def _parse_atom_line(self, line):
        """
        解析ATOM/HETATM行的固定宽度字段, 返回一个dict或自定义的Atom对象.
        参考PDB format对各列字段的定义.
        """

        # 以下针对固定宽度的列进行截取/strip
        # 下标(基于0)    长度
        #  0-5   ->  record type (ATOM/HETATM)
        #  6-11  ->  serial (atom序号)
        # 12     ->  空格
        # 12-16  ->  atom name
        # 17     ->  altLoc
        # 17-20  ->  resName
        # 21     ->  chainID
        # 22-25  ->  resSeq
        # 26     ->  iCode
        # 30-37  ->  x
        # 38-45  ->  y
        # 46-53  ->  z
        # 54-59  ->  occupancy
        # 60-65  ->  tempFactor
        # 76-77  ->  element
        # 78-79  ->  charge
        # 不同PDB版本略有差异,这里做一个常见处理.

        atom_serial = int(line[6:11].strip())  # 原子序号
        atom_name = line[12:16].strip()  # 原子名
        alt_loc = line[16].strip()  # 可选
        res_name = line[17:20].strip()  # 残基名
        chain_id = line[21].strip()  # 链ID
        res_seq = line[22:26].strip()  # 残基序号(注意可能包含文字)
        i_code = line[26].strip()  # 插入代码
        x = float(line[30:38].strip())  # x坐标
        y = float(line[38:46].strip())  # y坐标
        z = float(line[46:54].strip())  # z坐标
        occupancy = line[54:60].strip()  # 占有率
        temp_factor = line[60:66].strip()  # B因子
        element = line[76:78].strip() if len(line) >= 78 else ""
        charge = line[78:80].strip() if len(line) >= 80 else ""

        # 构造成dict(也可用自定义类)
        atom_dict = {
            "record_type": line[0:6].strip(),  # "ATOM"/"HETATM"
            "serial": atom_serial,
            "name": atom_name,
            "altLoc": alt_loc,
            "resName": res_name,
            "chainID": chain_id,
            "resSeq": res_seq,
            "iCode": i_code,
            "x": x,
            "y": y,
            "z": z,
            "occupancy": occupancy,
            "tempFactor": temp_factor,
            "element": element,
            "charge": charge
        }
        return atom_dict

    def _parse_conect_line(self, line):
        """
        解析CONECT行, 将配体/异质原子之间的连接关系记录下来.
        一般格式:
          "CONECT" + atom_serial + bonded_serial(s) ...
        例如: CONECT  455  458  472  478
        """

        # 取固定宽度,可以也只取[6:]按空格split
        fields = line.split()
        # fields[0] 应该是 "CONECT"
        if len(fields) < 2:
            return

        atom_serial_main = int(fields[1])
        bonded_list = []
        for f in fields[2:]:
            try:
                bonded_list.append(int(f))
            except ValueError:
                pass

        if atom_serial_main not in self.conect_info:
            self.conect_info[atom_serial_main] = set()
        # 把bonded_list加入
        for b in bonded_list:
            self.conect_info[atom_serial_main].add(b)

        # 如果需要，我们也可以做对称:
        for b in bonded_list:
            if b not in self.conect_info:
                self.conect_info[b] = set()
            self.conect_info[b].add(atom_serial_main)

    # ========== 一些查询辅助函数 ==========

    def get_atoms(self, record_type=None, chainID=None, resName=None):
        """
        简易查询函数: 根据record_type("ATOM"/"HETATM"), chainID, resName等条件
        返回匹配的原子列表.
        """
        results = []
        for atom in self.atoms:
            if record_type and atom["record_type"] != record_type:
                continue
            if chainID and atom["chainID"] != chainID:
                continue
            if resName and atom["resName"] != resName:
                continue
            results.append(atom)
        return results

    def group_by_residue(self):
        """
        把ATOM/HETATM分组: (chainID, resSeq, resName) -> list of atoms
        方便后续识别蛋白片段与配体, 以及做坐标封端等处理.
        """
        residue_dict = defaultdict(list)
        for atom in self.atoms:
            key = (atom["chainID"], atom["resSeq"], atom["resName"])
            residue_dict[key].append(atom)
        return residue_dict

    def is_ligand(self, resName):
        """
        简易判断: 如果resName是 'HOH'或蛋白标准三字母(ALA, PHE, GLY...),
        则认为不是配体. 否则可能是配体(例如 '2RV').

        实际可维护一个标准残基列表 / ligand列表 / water等.
        """
        # 粗略处理:
        protein_residues = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
            'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
            # etc. ...
        }
        if resName in protein_residues or resName == 'HOH':
            return False
        return True

    # ...更多功能可依需求扩展...


if __name__ == "__main__":
    # 示例使用
    pdb_path = "example.pdb"  # 放你的pdb路径
    pdb_struct = PDBStructure(pdb_path)

    print(f"总共读取到 {len(pdb_struct.atoms)} 个原子记录.")
    print(f"连接信息(CONECT) keys: {list(pdb_struct.conect_info.keys())[:10]} ...")

    # 按残基分组
    by_res = pdb_struct.group_by_residue()
    print("Residues found:", list(by_res.keys())[:5], "...")

    # 找某个配体resName, 例如 '2RV'
    ligand_atoms = pdb_struct.get_atoms(record_type="HETATM", resName="2RV")
    print(f"配体2RV共 {len(ligand_atoms)} 个原子.")
