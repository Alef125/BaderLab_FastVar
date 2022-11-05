"""
Oct 15, 2022
This module, handles cell-type markers.

'./Cell-type marker/PanglaoDB_markers_27_Mar_2020.tsv':
8286 rows,
14 columns:
['species', 'official gene symbol', 'cell type', 'nicknames',
       'ubiquitousness index', 'product description', 'gene type',
       'canonical marker', 'germ layer', 'organ', 'sensitivity_human',
       'sensitivity_mouse', 'specificity_human', 'specificity_mouse']
"""
import pandas as pd
import pickle as pkl


class CellTypeMarker:
    """
    This class, handles all queries about cell-type markers
    """
    def __init__(self, path_to_tsv: str):
        self.cell_marker_info = pd.read_csv(path_to_tsv, sep='\t')
        # print(self.cell_marker_info[['official gene symbol', 'cell type', 'nicknames']].iloc[1:100])
        self.gene_to_cell_type_dict = pd.Series(self.cell_marker_info['cell type'].values,
                                                index=self.cell_marker_info['official gene symbol']).to_dict()

    def save_gene_to_cell_type_dict(self, filepath_to_save: str):
        with open(filepath_to_save, 'wb') as file:
            pkl.dump(self.gene_to_cell_type_dict, file, protocol=pkl.HIGHEST_PROTOCOL)
