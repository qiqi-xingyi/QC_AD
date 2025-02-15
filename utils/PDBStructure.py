# --*-- conding:utf-8 --*--
# @Time : 2/13/25 8:17 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : PDBStructure.py

import re
from collections import defaultdict
import os

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

    def get_xyz_geometry(self, atoms_list):
        """
        将给定 atoms_list(里面是 dict {x,y,z,element,...})
        转成量子流程常用的 geometry=[(symbol,(x,y,z)),...].
        如果 element为空,可回退 atom_name首字母.
        """
        out_geometry = []
        for atm in atoms_list:
            symbol = atm["element"] if atm["element"] else atm["name"][0]
            coords = (atm["x"], atm["y"], atm["z"])
            out_geometry.append((symbol, coords))
        return out_geometry

    def export_subsystems_for_quantum(self, output_dir, ligand_resname="2RV"):
        """
        读取自身解析得到的atoms信息，将全PDB系统分解成若干子系统，
        并将每个子系统的原子坐标以XYZ文件形式输出到指定文件夹(output_dir).

        示例逻辑:
        1) 查找所有蛋白残基(默认: 标准氨基酸, chain==A等),
        2) 查找配体(resName == ligand_resname),
        3) 对每个蛋白残基 + 配体 => 生成xyz文件,
        4) 也输出该残基单独xyz, 以及配体单独xyz(只一次).

        可根据需要作微调.
        """
        # 0) 如果文件夹不存在就创建
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # 1) 将原子按 (chain,resSeq,resName) 分组
        residue_dict = defaultdict(list)
        for atm in self.atoms:
            key = (atm["chainID"], atm["resSeq"], atm["resName"])
            residue_dict[key].append(atm)

        # 2) 识别配体碎片(假设resName==ligand_resname)
        #    也可能有多个配体, 这里简化只取第一个
        ligand_key = None
        ligand_atoms = []
        for (chain, rseq, rname), at_list in residue_dict.items():
            if rname == ligand_resname:
                ligand_key = (chain, rseq, rname)
                ligand_atoms = at_list
                break

        if not ligand_atoms:
            print(f"[Warning] No ligand with resName={ligand_resname} found in PDB.")
        else:
            # 先输出一下配体单独的XYZ(如果需要)
            ligand_filename = os.path.join(output_dir,
                                           f"ligand_{ligand_key[0]}_{ligand_key[1]}_{ligand_key[2]}.xyz")
            self._write_xyz(ligand_atoms, ligand_filename,
                            comment=f"Ligand {ligand_key}")

        # 3) 输出每个蛋白残基(当子系统)
        #    定义: 标准氨基酸 or chain=='A' etc. 这里举例只排除 "HOH" 和 ligand_resname
        protein_keys = []
        for key in residue_dict.keys():
            chain, rseq, rname = key
            # 简单判断: 如果 rname != 'HOH' and != ligand_resname => 视作蛋白
            if rname not in ("HOH", ligand_resname):
                protein_keys.append(key)

        # 4) 针对每个蛋白残基 => 输出单独XYZ, 以及与配体组合XYZ
        for key in protein_keys:
            chain, rseq, rname = key
            res_atoms = residue_dict[key]

            # 4.1 输出蛋白残基单独xyz
            res_filename = os.path.join(output_dir, f"res_{chain}_{rseq}_{rname}.xyz")
            self._write_xyz(res_atoms, res_filename,
                            comment=f"Single residue: {chain} {rseq} {rname}")

            # 4.2 如果有配体 => 输出蛋白残基 + 配体
            if ligand_atoms:
                combined = res_atoms + ligand_atoms
                comb_filename = os.path.join(output_dir, f"res_{chain}_{rseq}_{rname}_plus_ligand.xyz")
                self._write_xyz(combined, comb_filename,
                                comment=f"Residue {chain} {rseq} {rname} + ligand {ligand_key}")

        print(f"[Done] Exported subsystem XYZ files into '{output_dir}'.")

    # ---------------------------------
    # 以下是一个辅助写xyz的方法
    # ---------------------------------
    def _write_xyz(self, atoms_list, xyz_path, comment=""):
        """
        将给定 atoms_list(原子dict的列表) 写成 .xyz 文件.
        coords => x,y,z
        element => element or fallback.

        .xyz格式:
          第一行: 原子数
          第二行: 注释(可选)
          第三行起: symbol  x  y  z
        """
        # 准备 geometry data
        geometry = []
        for atm in atoms_list:
            symbol = atm["element"] if atm["element"] else atm["name"][0]
            x, y, z = atm["x"], atm["y"], atm["z"]
            geometry.append((symbol, x, y, z))

        with open(xyz_path, "w") as f:
            f.write(f"{len(geometry)}\n")
            f.write(comment + "\n")
            for (sym, x, y, z) in geometry:
                f.write(f"{sym}  {x:.3f}  {y:.3f}  {z:.3f}\n")

        print(f"[Exported] {xyz_path}  (#atoms={len(geometry)})")

    def export_subsystems_for_quantum_in_subfolders(self, output_dir, ligand_resname="2RV"):
        """
        将蛋白-配体体系分解为 [每个氨基酸(残基) + 配体] 的子系统，
        并在 output_dir 下为每个子系统创建一个专门的文件夹，里面存:
          - res.xyz (该残基单独)
          - ligand.xyz (配体单独)
          - res_plus_ligand.xyz (二者组合)

        这样可方便二体展开或结合能计算的后续处理。
        """
        # 如果output_dir不存在就创建
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # 1) 将原子按 (chain,resSeq,resName) 分组
        residue_dict = defaultdict(list)
        for atm in self.atoms:
            key = (atm["chainID"], atm["resSeq"], atm["resName"])
            residue_dict[key].append(atm)

        # 2) 找到配体(假设resName=ligand_resname)
        ligand_key = None
        ligand_atoms = []
        for (chain, rseq, rname), at_list in residue_dict.items():
            if rname == ligand_resname:
                ligand_key = (chain, rseq, rname)
                ligand_atoms = at_list
                break

        if not ligand_atoms:
            print(f"[Warning] No ligand with resName='{ligand_resname}' found in PDB.")
        else:
            print(f"Found ligand {ligand_key}, #atoms={len(ligand_atoms)}")

        # 3) 准备输出：对每个蛋白残基 => 创建子文件夹 => 写 xyz
        #    这里简单过滤: 排除 HOH(水) 与 ligand_resname
        protein_keys = []
        for (chain, rseq, rname) in residue_dict.keys():
            if rname not in ("HOH", ligand_resname):
                protein_keys.append((chain, rseq, rname))
        protein_keys.sort()

        # 4) 逐个处理蛋白残基
        for key in protein_keys:
            chain, rseq, rname = key
            res_atoms = residue_dict[key]
            # 子系统文件夹名：如 sub_A_267_PHE
            subfolder_name = f"sub_{chain}_{rseq}_{rname}"
            subfolder_path = os.path.join(output_dir, subfolder_name)
            if not os.path.exists(subfolder_path):
                os.makedirs(subfolder_path)

            # 4.1 输出res.xyz (该残基单独)
            res_filename = os.path.join(subfolder_path, "res.xyz")
            self._write_xyz(res_atoms, res_filename,
                            comment=f"Residue {chain} {rseq} {rname}")

            # 4.2 如果有配体 => 写 ligand.xyz (同一个子文件夹)
            if ligand_atoms:
                ligand_filename = os.path.join(subfolder_path, "ligand.xyz")
                # 这里可以把配体也单独写一次
                self._write_xyz(ligand_atoms, ligand_filename,
                                comment=f"Ligand {ligand_key}")

                # 4.3 写res_plus_ligand.xyz
                combined_atoms = res_atoms + ligand_atoms
                comb_filename = os.path.join(subfolder_path, "res_plus_ligand.xyz")
                self._write_xyz(combined_atoms, comb_filename,
                                comment=f"Residue+Ligand => {key}+{ligand_key}")

        print(f"[Done] Exported all subsystems into subfolders under '{output_dir}'.")





