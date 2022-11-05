"""
Created: 15 Oct 2022
Project: Cell-type marker heritability
"""
from CellMarkerHandler import CellTypeMarker
from S2G import S2GTranslator
from CellMarkerAnnotationMaker import CellTypeSNPAnnotation
import pickle as pkl


def main():
    s2g_dict_filepath = './SNP to Gene/S2G_dict.pkl'
    g2ct_dict_filepath = './Cell-type marker/G2CT_dict.pkl'
    # ################ Making SNP-to-Gene dictionary ####################
    # s2g = S2GTranslator(path_to_tsv='./SNP to Gene/gwas_gold_standards.191108.tsv')
    # s2g = S2GTranslator(path_to_bim='./SNP to Gene/unrelated_bed.bim',
    #                     folder_of_chromosomes='./SNP to Gene/Chromosomes - Detailed')
    # s2g.save_snp_to_gene_dict(filepath_to_save=s2g_dict_filepath)
    # ################ Making Gene-to-Cell type dictionary ####################
    # cell_type_marker = CellTypeMarker(path_to_tsv='./Cell-type marker/PanglaoDB_markers_27_Mar_2020.tsv')
    # cell_type_marker.save_gene_to_cell_type_dict(filepath_to_save=g2ct_dict_filepath)
    # ################ Merging => Making SNP-to-Cell type dictionary ###############
    annotation_maker = CellTypeSNPAnnotation(s2g_dict_filepath=s2g_dict_filepath,
                                             g2ct_dict_filepath=g2ct_dict_filepath)
    annotation_maker.make_annotation_file(all_snps_filepath='./SNP to Gene/unrelated_bed.bim',
                                          filepath_to_save='./CT_Annotations.txt')


if __name__ == "__main__":
    main()
