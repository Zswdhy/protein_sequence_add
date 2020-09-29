import time
import pandas as pd
import mygene
import requests
import re
import random


def get_ensembl_id(path):
    res = pd.read_excel(path)
    res = pd.DataFrame(res)
    # ensembl ID
    Gene_ID = res["Gene_ID"].tolist()
    return res, Gene_ID


def emsemblid_change_uniprotid(data, ensembl_id_list):
    ensembl_id_list = list(set(ensembl_id_list))
    print("维度", data.shape, "data_len", len(data), type(data))

    uniprot_id_list = []
    # 在ensembl数据库中检索不到的ID
    not_ensembl = []
    # 以ensemble ID为键，uniprot ID为值
    ensembl_mapping_uniprot = {}

    mg = mygene.MyGeneInfo()
    xli = ensembl_id_list
    out = mg.querymany(xli, scopes="ensembl.gene", fields="uniprot", species="human")

    for item in out:
        ensembl_mapping_uniprot[item["query"]] = []

        # 此时需要考虑将数据删除，在uniprot中检索不到的ensembl
        if "uniprot" not in item:
            not_ensembl.append(item['query'])
        # 符合要求的ensembl
        else:
            # print("Swiss-Prot", item["query"], item["uniprot"]["Swiss-Prot"], type(item["uniprot"]["Swiss-Prot"]))
            if isinstance(item["uniprot"]["Swiss-Prot"], list):
                ensembl_mapping_uniprot[item["query"]] = item["uniprot"]["Swiss-Prot"]
                uniprot_id_list += item["uniprot"]["Swiss-Prot"]
            elif isinstance(item["uniprot"]["Swiss-Prot"], str):
                ensembl_mapping_uniprot[item["query"]].append(item["uniprot"]["Swiss-Prot"])
                uniprot_id_list.append(item["uniprot"]["Swiss-Prot"])

            if "TrEMBL" in item["uniprot"]:
                # print("TrEMBL", item["query"], item["uniprot"]["TrEMBL"], type(item["uniprot"]["TrEMBL"]))
                if isinstance(item["uniprot"]["TrEMBL"], list):
                    ensembl_mapping_uniprot[item["query"]] += item["uniprot"]["TrEMBL"]
                    uniprot_id_list += item["uniprot"]["TrEMBL"]
                elif isinstance(item["uniprot"]["TrEMBL"], str):
                    ensembl_mapping_uniprot[item["query"]].append(item["uniprot"]["TrEMBL"])
                    uniprot_id_list.append(item["uniprot"]["TrEMBL"])

    return ensembl_mapping_uniprot, uniprot_id_list


# 请求的时候，利用set集合，存储为字典，减少检索的次数
def get_protein_seq(ensembl_mapping_uniprot):
    protein_seq = {}
    # https://www.uniprot.org/uniprot/P04629.fasta
    url = " https://www.uniprot.org/uniprot/"
    for key, values in ensembl_mapping_uniprot.items():
        protein_seq[key] = []
        for val in values:
            request_url = url + val + ".fasta"
            print("request_url", request_url)
            response = requests.get(request_url).text
            time.sleep(random.random())
            res = re.match(r"^>.*\s((\w+|\s)+)", response).group(1).replace("\n", "")
            res = "".join(res.split("\n"))
            protein_seq[key].append(res)

    return protein_seq


def add_protein():
    for indexs in original_data.index:
        # 逐行查看
        # print(remove_data.loc[indexs].values[0:-1])
        add_len = 25
        # Norm_peptide
        Norm_peptide = original_data.loc[indexs].values[1]
        # ensembl ID
        ensembl = original_data.loc[indexs].values[5]
        # 突变序列
        acid_change = original_data.loc[indexs].values[7]
        left, right = acid_change.split("/")
        # 原始蛋白质序列，突变位置
        acid_change_position = original_data.loc[indexs].values[10]

        if ensembl in protein_seq:
            for item in protein_seq[ensembl]:
                # 正常序列起始位置
                pos = item.find(Norm_peptide)
                if pos == -1:
                    continue
                # 病变位置
                pro_pop = pos + acid_change_position

                # 正常情况
                l = item[pro_pop - add_len - 1:pro_pop - 1]
                r = item[pro_pop + len(left) - 1:pro_pop + add_len]
                # 左边值不够
                if pos + acid_change_position - add_len < 0:
                    l = "-" * (pos - 1) + item[:pos + acid_change_position - 1]
                # 右边值不够
                if pos + acid_change_position + add_len > len(item):
                    r = item[pos + acid_change_position + len(left) - 1:] + "-" * (
                            add_len - (len(item) - pos) + 2)

                print(item)
                print(Norm_peptide, "pos:", pos, "pro_pop", pro_pop)
                print("left:", l, len(l), "acid_change:", acid_change, "right:", r, len(r))
                print(ensembl)
                norm_protein = l + left + r
                mult_protein = l + right + r
                print("正常", norm_protein, len(norm_protein))
                print("突变", mult_protein, len(mult_protein))
                print("-" * 200)
                # df.loc[indexs, "姓名"] = "david"
                original_data.loc[indexs, "Norm_peptide"] = norm_protein
                original_data.loc[indexs, "Mut_peptide"] = mult_protein
                break
                
    original_data.to_excel(f"{rename}_change.xlsx")
    return original_data


if __name__ == '__main__':
    old_time = time.time()
    path = r"Y088N_WES_SNP_VQSR_SELECTED.xlsx"

    # 截取文件名称，命名时使用
    rename = path.split(".")[0]

    original_data, ensembl_id_list = get_ensembl_id(path)
    ensembl_mapping_uniprot, _ = emsemblid_change_uniprotid(original_data, ensembl_id_list)

    protein_seq = get_protein_seq(ensembl_mapping_uniprot)
    add_protein()
    new_time = time.time()
    print("run time", new_time - old_time)
    print("process over")
