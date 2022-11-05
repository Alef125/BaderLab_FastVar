"""
Nov 4, 2022
This module, uses S2G and G2CT dictionaries to build the final annotation file
    to be used in the FastVar
"""
import pickle as pkl
import pandas as pd
# from pandas_plink import read_plink, read_plink1_bin
from pandas_plink._read import _read_bim


class CellTypeSNPAnnotation:
    def __init__(self,
                 s2g_dict_filepath: str,
                 g2ct_dict_filepath: str):
        with open(s2g_dict_filepath, 'rb') as file:
            s2g_dict = pkl.load(file)
        with open(g2ct_dict_filepath, 'rb') as file:
            g2ct_dict = pkl.load(file)
        self.s2ct_dict = {}
        for _snp in s2g_dict.keys():
            _gene = s2g_dict[_snp]
            if _gene in g2ct_dict.keys():
                self.s2ct_dict[_snp] = g2ct_dict[_gene]
            else:
                # print("Bad gene: " + _gene)
                pass
        # Final dict length: 52122
        self.all_genes_snps = list(s2g_dict.keys())
        self.snps_for_cell_types = {}

    def find_snps_for_cell_types(self):
        for _snp, _cell_type in self.s2ct_dict.items():
            if _cell_type in self.snps_for_cell_types.keys():
                self.snps_for_cell_types[_cell_type].append(_snp)
            else:
                self.snps_for_cell_types[_cell_type] = [_snp]

    def make_annotation_file(self,
                             all_snps_filepath: str,
                             filepath_to_save: str) -> None:
        all_snps = _read_bim(all_snps_filepath)['snp'].values.tolist()
        self.find_snps_for_cell_types()
        all_cell_types = list(self.snps_for_cell_types.keys())  # 169 cell_types
        annotations = pd.DataFrame(index=all_snps)
        annotations['GCTA'] = 1
        annotations['All genes'] = 0
        annotations['All genes'].loc[self.all_genes_snps] = 1
        for cell_type in all_cell_types[:50]:
            annotations[cell_type] = 0
            annotations[cell_type].loc[self.snps_for_cell_types[cell_type]] = 1
        annotations.to_csv(filepath_to_save, sep=' ', index=False)
