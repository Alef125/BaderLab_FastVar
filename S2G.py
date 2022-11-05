"""
Oct 24, 2022
This module, handles SNP to Gene module.

'./SNP to Gene/gwas_gold_standards.191108.tsv':
2435 rows x 31 columns
Columns with example:
[
association_info.ancestry                                                                 EUR
association_info.doi                                                                      NaN
association_info.gwas_catalog_id                                                   GCST000386
association_info.neg_log_pval                                                         323.301
association_info.otg_id                                                            GCST000386
association_info.pubmed_id                                                         19414484.0
association_info.url                                                                      NaN
gold_standard_info.evidence.class                                              expert curated
gold_standard_info.evidence.confidence                                                   High
gold_standard_info.evidence.curated_by                                            Eric Fauman
gold_standard_info.evidence.description     UGT1A1 encodes an enzyme with bilirubin glucur...
gold_standard_info.evidence.pubmed_id                                                01898728
gold_standard_info.evidence.source                                                        NaN
gold_standard_info.gene_id                                                    ENSG00000242366
gold_standard_info.highest_confidence                                                    High
metadata.comments                                                                         NaN
metadata.date_added                                                                2019-05-17
metadata.reviewed_by                                                              Ed Mountjoy
metadata.set_label                                                                     ProGeM
metadata.submitted_by                                                             Eric Fauman
metadata.tags                                                                 metabolite|mQTL
sentinel_variant.alleles.alternative                                                        T
sentinel_variant.alleles.reference                                                          G
sentinel_variant.locus_GRCh37.chromosome                                                    2
sentinel_variant.locus_GRCh37.position                                              234672639
sentinel_variant.locus_GRCh38.chromosome                                                    2
sentinel_variant.locus_GRCh38.position                                              233763993
sentinel_variant.rsid                                                               rs6742078
trait_info.ontology                                                               HMDB0000054
trait_info.reported_trait_name                                               Bilirubin levels
trait_info.standard_trait_name                                                      Bilirubin
]


'./SNP to Gene/Chromosomes - Detailed/cS2G.2.SGscore'  (replace 2 with chromosome number):
SNP	GENE	cS2G	INFO
2:10514_C_T	FAM110C	1	|ABC=1
2:10515_G_A	FAM110C	1	|ABC=1
2:10533_A_T	FAM110C	1	|ABC=1
2:10554_G_A	FAM110C	1	|ABC=1
2:10560_G_A	FAM110C	1	|ABC=1
2:10566_G_A	FAM110C	1	|ABC=1
2:10572_AC_A	FAM110C	1	|ABC=1
2:10574_C_G	FAM110C	1	|ABC=1
...

'./SNP to Gene/unrelated_bed.bim':
1	rs3115860	0	753405	C	A
1	rs3131970	0	753425	T	C
1	rs12184325	0	754105	T	C
1	rs3131969	0	754182	A	G
...

'./SNP to Gene/allsnps.txt'
rs201106462 1 10352
rs534229142 1 10511
1:10616_CCGCCGTTGCAAAGGCGCGCCG_C 1 10616
rs575272151 1 11008
rs544419019 1 11012
rs540538026 1 13110
rs62635286 1 13116
rs62028691 1 13118
rs531730856 1 13273
1:13289_CCT_C 1 13289
rs568927457 1 13453
...
"""
import pandas as pd
# from warnings import warn
import pickle as pkl
import os


BIM_HEADERS = ['Chrom', 'SNP', 'Pos', 'bp_coord', 'Minor_allele', 'Major_allele']


class S2GTranslator:
    """
    This class, handles all queries about cell-type markers
    """
    def __init__(self,
                 path_to_tsv: str = None,
                 folder_of_chromosomes: str = None,
                 path_to_bim: str = None,
                 path_to_all_snps: str = None):
        if path_to_tsv:
            self.cell_marker_info = pd.read_csv(path_to_tsv, sep='\t')
            # print(self.cell_marker_info.iloc[8])
        elif folder_of_chromosomes:
            if path_to_bim:
                self.bim_df = pd.read_csv(path_to_bim, sep='\t', header=None, names=BIM_HEADERS)
            elif path_to_all_snps:
                raise NotImplementedError("path_to_all_snps not implemented")
            else:
                raise AttributeError("At least one of the path_to_bim or path_to_all_snps should be assigned")
            # ###################################################################################
            self.chromosomes_snp_dict = {}
            for _chrom_id in range(1, 23):
                _chromosome_snp_dict = {}
                _chrom_file_name = 'cS2G.' + str(_chrom_id) + '.SGscore'
                _chrom_filepath = os.path.join(folder_of_chromosomes, _chrom_file_name)
                _chrom_df = pd.read_csv(_chrom_filepath, sep='\t')
                for _, _row in _chrom_df.iterrows():
                    _snp_info = _row['SNP'].split(':')
                    _gene = _row['GENE']
                    if _snp_info[0] != str(_chrom_id):  # 'rs...', direct mapping
                        _snp_key = _snp_info[0]
                        _chromosome_snp_dict[_snp_key] = _gene
                        # warn("Chromosomal index of this SNP is not consistent")
                    else:
                        _snp_details = _snp_info[1].split('_')  # SNP_id, original_base, variant_base
                        _snp_key = _snp_details[0]  # _row['SNP']
                        # Note: _snp_key's duplication is not checked
                        _chromosome_snp_dict[_snp_key] = _gene
                self.chromosomes_snp_dict[_chrom_id] = _chromosome_snp_dict
            self.snp_to_gene_dict = {}
            self.make_snp_to_gene_dict()
        else:
            raise AttributeError("At least one of the path_to_tsv or folder_of_chromosomes should be assigned")

    def make_snp_to_gene_dict(self):
        for _, _row in self.bim_df.iterrows():
            _chrom_dict = self.chromosomes_snp_dict[_row['Chrom']]
            _snp = _row['SNP']
            if _snp in _chrom_dict.keys():
                self.snp_to_gene_dict[_snp] = _chrom_dict[_snp]
            elif _row['bp_coord'] in _chrom_dict:
                self.snp_to_gene_dict[_snp] = _chrom_dict[_row['bp_coord']]
            else:
                print("Bad item " + str(_row['Chrom']) + ":" + _row['SNP'])
        print(len(self.snp_to_gene_dict.keys()))

    def save_snp_to_gene_dict(self, filepath_to_save: str):
        with open(filepath_to_save, 'wb') as file:
            pkl.dump(self.snp_to_gene_dict, file, protocol=pkl.HIGHEST_PROTOCOL)
