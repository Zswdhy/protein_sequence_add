import time
import pandas as pd
import mygene
import requests
import re


def get_ensembl_id(path):
    res = pd.read_excel(path)
    res = pd.DataFrame(res)
    # 正常的protein seq
    Norm_peptide = res["Norm_peptide"].tolist()
    # 突变的protein seq
    Mut_peptide = res["Mut_peptide"].tolist()
    # ensembl ID
    Gene_ID = res["Gene_ID"].tolist()
    # 突变序列
    Amino_Acid_Change = res["Amino_Acid_Change"].tolist()
    # 突变位置
    peptide_position = res["peptide_position"].tolist()
    return res, Gene_ID


def emsemblid_change_uniprotid(data, ensembl_id_list):
    """

    :param data:原始数据
    :param ensembl_id_list:Gene_id列表
    :return:
    """
    print("维度", data.shape, "data_len", len(data), type(data))
    # uniprot ID
    uniprot_id_list = []
    # 不存在的ensembl ID
    not_ensembl = []
    # 映射 以ensembl ID 为键，uniprot ID 为值
    ensembl_mapping_uniprot = {}

    mg = mygene.MyGeneInfo()
    xli = ensembl_id_list
    out = mg.querymany(xli, scopes="ensembl.gene", fields="uniprot", species="human")

    for item in out:
        """ item 
        {'query': 'ENSG00000214562', '_id': '728130', '_score': 22.992172, 
        'uniprot': {'Swiss-Prot': ['Q5VT03', 'Q8IVF1'], 'TrEMBL': ['A0A075B6P9', 'U3KPT3']}}
        """
        # 此时需要考虑将数据删除，在uniprot中检索不到的ensembl
        if "uniprot" not in item:
            not_ensembl.append(item['query'])
        # 符合要求的ensembl
        else:
            # 一个ensembl对应多个uniprotID
            if isinstance(item["uniprot"]["Swiss-Prot"], list):
                ensembl_mapping_uniprot[item["query"]] = item["uniprot"]["Swiss-Prot"][0]
                uniprot_id_list.append(item["uniprot"]["Swiss-Prot"][0])
            else:
                ensembl_mapping_uniprot[item["query"]] = item["uniprot"]["Swiss-Prot"]
                uniprot_id_list.append(item["uniprot"]["Swiss-Prot"])
    print("待删除的ID", not_ensembl)
    # 将不存在的ensembl ID 写入 txt文本
    with open(f"not_ensembl_{rename}.txt", "a") as f_w:
        for item in set(not_ensembl):
            f_w.write(item + "\n")
        f_w.close()

    print("不符合要求的ensemblID", len(not_ensembl))
    # 删除
    data = data[~data.Gene_ID.isin(not_ensembl)]
    print("修改后的维度", data.shape)
    return data, ensembl_mapping_uniprot, uniprot_id_list


# 请求的时候，利用set集合，存储为字典，减少检索的次数
def get_protein_seq(uniprot_id_list, ensembl_mapping_uniprot):
    # print("uniprot_id_list",uniprot_id_list)
    uniprot_set = list(set(uniprot_id_list))
    # print(uniprot_set)
    print("去重后的URL链接个数", len(uniprot_set))
    # print(remove_data)
    print("处理后的数据", type(remove_data), "长度", len(remove_data))
    # 扩增长度
    add_len = 25

    protein_seq = {}
    #  https://www.uniprot.org/uniprot/P04629.fasta
    url = " https://www.uniprot.org/uniprot/"
    for uniprot in uniprot_set:
        request_url = url + uniprot + ".fasta"
        print("request_url", request_url)
        response = requests.get(request_url).text
        # print("response",response)
        res = re.match(r"^>.*\s((\w+|\s)+)", response).group(1).replace("\n", "")
        protein_seq[list(ensembl_mapping_uniprot.keys())[list(ensembl_mapping_uniprot.values()).index(uniprot)]] = res

    for indexs in remove_data.index:
        # 逐行查看
        # print(remove_data.loc[indexs].values[0:-1])

        # Norm_peptide
        Norm_peptide = remove_data.loc[indexs].values[1]
        # ensembl ID
        ensembl = remove_data.loc[indexs].values[5]
        # 突变序列
        acid_change = remove_data.loc[indexs].values[7]
        left, right = acid_change.split("/")
        # 原始蛋白质序列，突变位置
        acid_change_position = remove_data.loc[indexs].values[10]
        # print(Norm_peptide, ensembl, acid_change, left, right, acid_change_position)
        # print("-" * 200)
        if ensembl in protein_seq:
            # 正常序列起始位置
            pos = protein_seq[ensembl].find(Norm_peptide)
            # 病变位置
            pro_pop = pos + acid_change_position - 1
            l = protein_seq[ensembl][pro_pop - add_len:pro_pop]
            r = protein_seq[ensembl][pro_pop + 1:pro_pop + add_len + 1]

            norm_protein = l + left + r
            if right == "-":
                mult_protein = l + right + right + r
            else:
                mult_protein = l + right + r
            print(print(remove_data.loc[indexs].values[0:-1]))
            print("正常", norm_protein)
            print("突变", mult_protein)
            print("-" * 200)
            # df.loc[indexs, "姓名"] = "david"
            remove_data.loc[indexs, "Norm_peptide"] = norm_protein
            remove_data.loc[indexs, "Mut_peptide"] = mult_protein

    remove_data.to_excel(f"{rename}_change.xlsx")
    return remove_data


if __name__ == '__main__':
    old_time = time.time()
    # path = "../data/example.xls"
    path = r"Y088N_WES_SNP_VQSR_SELECTED.xlsx"

    # 截取文件名称，命名时使用
    rename = path.split(".")[0]

    old_data, ensembl_id_list = get_ensembl_id(path)
    remove_data, ensembl_mapping_uniprot, uniprot_id_list = emsemblid_change_uniprotid(old_data, ensembl_id_list)
    get_protein_seq(uniprot_id_list, ensembl_mapping_uniprot)
    new_time = time.time()
    print("run time", new_time - old_time)
    print("process over")
