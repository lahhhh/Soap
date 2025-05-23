#pragma once

#include <QApplication>

#define SOAP_RESOURCE_DIR qApp->applicationDirPath() + "/Resources/"

#define FILE_SINGLE_CELL_RNA_PIPELINE SOAP_RESOURCE_DIR "single_cell_rna_pipeline.txt"

#define FILE_FAQ SOAP_RESOURCE_DIR "SOAP_GUIDE" 

#define FILE_HUMAN_CELLMARKER SOAP_RESOURCE_DIR "CellMarker/CELLMARKER_CELLMARKER_HUMAN"
#define FILE_MOUSE_CELLMARKER SOAP_RESOURCE_DIR "CellMarker/CELLMARKER_CELLMARKER_MOUSE"
#define FILE_HUMAN_PANGLAODB SOAP_RESOURCE_DIR "CellMarker/CELLMARKER_PANGLAODB_HUMAN"
#define FILE_MOUSE_PANGLAODB SOAP_RESOURCE_DIR "CellMarker/CELLMARKER_PANGLAODB_MOUSE"

#define FILE_GO_PATH2NAME SOAP_RESOURCE_DIR "Enrichment/GO_PATH2NAME"

#define FILE_KEGG_PATH2NAME_HUMAN SOAP_RESOURCE_DIR "Enrichment/KEGG_PATH2NAME_HUMAN"
#define FILE_KEGG_PATH2NAME_MOUSE SOAP_RESOURCE_DIR "Enrichment/KEGG_PATH2NAME_MOUSE"

#define FILE_GO_PATH2SYMBOL_HUMAN SOAP_RESOURCE_DIR "Enrichment/GO_PATH2SYMBOL_HUMAN"
#define FILE_KEGG_PATH2SYMBOL_HUMAN SOAP_RESOURCE_DIR "Enrichment/KEGG_PATH2SYMBOL_HUMAN"

#define FILE_GO_SYMBOL2PATH_HUMAN SOAP_RESOURCE_DIR "Enrichment/GO_SYMBOL2PATH_HUMAN"
#define FILE_KEGG_SYMBOL2PATH_HUMAN SOAP_RESOURCE_DIR "Enrichment/KEGG_SYMBOL2PATH_HUMAN"

#define FILE_GO_PATH2SYMBOL_MOUSE SOAP_RESOURCE_DIR "Enrichment/GO_PATH2SYMBOL_MOUSE"
#define FILE_KEGG_PATH2SYMBOL_MOUSE SOAP_RESOURCE_DIR "Enrichment/KEGG_PATH2SYMBOL_MOUSE"

#define FILE_GO_SYMBOL2PATH_MOUSE SOAP_RESOURCE_DIR "Enrichment/GO_SYMBOL2PATH_MOUSE"
#define FILE_KEGG_SYMBOL2PATH_MOUSE SOAP_RESOURCE_DIR "Enrichment/KEGG_SYMBOL2PATH_MOUSE"

#define FILE_CNV_GENE_LOCATION_HUMAN SOAP_RESOURCE_DIR "CNV/CNV_GENE_LOCATION_HUMAN"
#define FILE_CNV_GENE_LOCATION_MOUSE SOAP_RESOURCE_DIR "CNV/CNV_GENE_LOCATION_MOUSE"

#define FILE_HUMAN_GRCH38_BLACKLIST_BED SOAP_RESOURCE_DIR "Genome/Human/GRCh38_unified_blacklist.bed"
#define FILE_HUMAN_GRCH38_2BIT SOAP_RESOURCE_DIR "Genome/Human/single_sequences.2bit"
#define FILE_HUMAN_HG38_GTF SOAP_RESOURCE_DIR "Genome/Human/hg38.gtf"
#define FILE_MOUSE_MM10_GTF SOAP_RESOURCE_DIR "Genome/Mouse/mm10.gtf"
#define FILE_HUMAN_HG38_MASK_GTF SOAP_RESOURCE_DIR "Genome/Human/hg38_rmsk.gtf"
#define FILE_MOUSE_MM10_MASK_GTF SOAP_RESOURCE_DIR "Genome/Mouse/mm10_rmsk.gtf"
#define FILE_HUMAN_GENE_LOCATION SOAP_RESOURCE_DIR "Genome/Human/hg38_gene_location"
#define FILE_HUMAN_TRANSCRIPT_LENGTH SOAP_RESOURCE_DIR "Genome/Human/HG38_REFSEQ_TRANSCRIPT_LENGTH"

#define FILE_SCENT_HUMAN_PPI SOAP_RESOURCE_DIR "SCENT/SCENT_PPI_HUMAN"

#define FILE_HUMAN_GENOME_GENOMIC_RANGE_SIF SOAP_RESOURCE_DIR "Genome/Human/HUMAN_GENOME_GENOMIC_RANGE.sif"

#define FILE_HUMAN_GSEA_CURATED SOAP_RESOURCE_DIR "GSEA/CURATED_HUMAN_202401"
#define FILE_MOUSE_GSEA_CURATED SOAP_RESOURCE_DIR "GSEA/CURATED_MOUSE_202401"
#define FILE_HUMAN_GSEA_ONTOLOGY SOAP_RESOURCE_DIR "GSEA/ONTOLOGY_HUMAN_202401"
#define FILE_MOUSE_GSEA_ONTOLOGY SOAP_RESOURCE_DIR "GSEA/ONTOLOGY_MOUSE_202401"

#define FILE_JASPAR2024_HUMAN SOAP_RESOURCE_DIR "Motif/JASPAR_HUMAN_2024"
#define FILE_PANDO_MOTIF SOAP_RESOURCE_DIR "Motif/PANDO_MOTIF"

#define FILE_EYE_ICON_PNG SOAP_RESOURCE_DIR "Image/eyeIcon.png"
#define FILE_DELETE_ICON_PNG SOAP_RESOURCE_DIR "Image/deleteIcon.png"
#define FILE_SOAP_ICON_JPG SOAP_RESOURCE_DIR "Image/mainwindowicon.jpg" 

#define FILE_CHARACTER_A_PNG SOAP_RESOURCE_DIR "Image/character_a.png"
#define FILE_CHARACTER_C_PNG SOAP_RESOURCE_DIR "Image/character_c.png"
#define FILE_CHARACTER_G_PNG SOAP_RESOURCE_DIR "Image/character_g.png"
#define FILE_CHARACTER_T_PNG SOAP_RESOURCE_DIR "Image/character_t.png"

#define FILE_HPCA_CELL_TYPE SOAP_RESOURCE_DIR "Annotation/hpca.sif"

#define FILE_CICERO_HG38_CHROMSIZE SOAP_RESOURCE_DIR "Cicero/hg38_chrom_sizes.txt"

#define FILE_CELLCHAT_HUMAN_INTERACTION SOAP_RESOURCE_DIR "CellChat/CELLCHAT_HUMAN_INTERACTION"
#define FILE_CELLCHAT_HUMAN_COMPLEX SOAP_RESOURCE_DIR "CellChat/CELLCHAT_HUMAN_COMPLEX"
#define FILE_CELLCHAT_HUMAN_GENEINFO SOAP_RESOURCE_DIR "CellChat/CELLCHAT_HUMAN_GENEINFO"
#define FILE_CELLCHAT_HUMAN_COFACTOR SOAP_RESOURCE_DIR "CellChat/CELLCHAT_HUMAN_COFACTOR"
#define FILE_CELLCHAT_MOUSE_INTERACTION SOAP_RESOURCE_DIR "CellChat/CELLCHAT_MOUSE_INTERACTION"
#define FILE_CELLCHAT_MOUSE_COMPLEX SOAP_RESOURCE_DIR "CellChat/CELLCHAT_MOUSE_COMPLEX"
#define FILE_CELLCHAT_MOUSE_GENEINFO SOAP_RESOURCE_DIR "CellChat/CELLCHAT_MOUSE_GENEINFO"
#define FILE_CELLCHAT_MOUSE_COFACTOR SOAP_RESOURCE_DIR "CellChat/CELLCHAT_MOUSE_COFACTOR"