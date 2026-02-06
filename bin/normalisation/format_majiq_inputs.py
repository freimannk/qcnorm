#!/usr/bin/env python3

import polars as pl
import argparse



def filter_phenotype( df: pl.DataFrame,
    phenotype_meta: pl.DataFrame, KEEP_XY=False):

    VALID_CHRS =  ["1","10","11","12","13","14","15","16","17","18","19",
                        "2","20","21","22","3","4","5","6","7","8","9"]
    VALID_GENE_TYPES = ["lincRNA", "lncRNA","protein_coding","IG_C_gene","IG_D_gene","IG_J_gene",
                       "IG_V_gene", "TR_C_gene","TR_D_gene","TR_J_gene", "TR_V_gene",
                       "3prime_overlapping_ncrna","known_ncrna", "processed_transcript",
                       "antisense","sense_intronic","sense_overlapping", "leafcutter"]
    if(KEEP_XY):
        VALID_CHRS = VALID_CHRS + ["X","Y"]

    if "chromosome" not in phenotype_meta.columns:
            raise ValueError("phenotype_meta must contain 'chromosome'")
    phenotype_meta = phenotype_meta.filter(pl.col("chromosome").is_in(VALID_CHRS))

    if "gene_type" not in phenotype_meta.columns:
            raise ValueError("phenotype_meta must contain 'gene_type'")
    phenotype_meta = phenotype_meta.filter(pl.col("gene_type").is_in(VALID_GENE_TYPES))

    # keep only matching phenotype rows in df
    df = df.join(
        phenotype_meta.select("phenotype_id"),
        on="phenotype_id",
        how="inner",
    )
    return df

def sample_qc(df: pl.DataFrame,
    sample_metadata: pl.DataFrame,
    filter_rna_qc=False,
    filter_genotype_qc=False):

    if filter_rna_qc:
        if "rna_qc_passed" not in sample_metadata.columns:
            raise ValueError("sample_metadata missing rna_qc_passed")
        sample_metadata = sample_metadata.filter(pl.col("rna_qc_passed"))

    if filter_genotype_qc:
        if "genotype_qc_passed" not in sample_metadata.columns:
            raise ValueError("sample_metadata missing genotype_qc_passed")
        sample_metadata = sample_metadata.filter(pl.col("genotype_qc_passed"))

    keep_samples = sample_metadata["sample_id"].to_list()

    keep_cols = ["phenotype_id"] + keep_samples
    df = df.select([c for c in keep_cols if c in df.columns])
    return df



def main(phenotype_metadata,psi_matrix,sample_metadata,QC,KEEP_XY):

    print("Loading phenotype metadata...")
    ph_meta = (
        pl.read_csv(phenotype_metadata, separator="\t",schema_overrides={"chromosome": pl.Utf8})
        .drop(["phenotype_id", "quant_id"], strict=False)
    )

    print("Loading phenotype matrix...")
    pheno = pl.read_csv(psi_matrix, separator="\t")
    pheno = (
        pheno
        .with_columns(
            pl.col("phenotype_id").str.replace("^gene:", "")
        )
        .with_columns(
            parts=pl.col("phenotype_id").str.split(":")
        )
        .with_columns([
            pl.col("parts").list.get(0).alias("gene_id"),
            (pl.col("parts").list.get(1) + ":" + pl.col("parts").list.get(2)).alias("lsv_coord"),
            (pl.col("parts").list.get(0) + ":" +
             pl.col("parts").list.get(1) + ":" +
             pl.col("parts").list.get(2)).alias("lsv_id"),
        ])
        .drop("parts")
    )

    gene_counts = (
    pheno
    .group_by("lsv_coord", maintain_order=True) 
    .agg([
        pl.col("gene_id").n_unique().alias("unique_genes")
    ])
    )
    # keep only lsv_coord with exactly 1 unique gene_id
    ok_coords = gene_counts.filter(pl.col("unique_genes") == 1).select("lsv_coord")

    # filter pheno to keep only these lsv_coord rows
    pheno = pheno.join(ok_coords, on="lsv_coord", how="inner")

    # -------------------------
    # Build phenotype metadata
    # -------------------------
    pheno_meta = (
        pheno
        .select(["gene_id", "phenotype_id", "lsv_id"])
        .join(ph_meta, on="gene_id", how="left")
        .with_columns([
            pl.col("lsv_id").alias("group_id"),
            pl.col("lsv_id").alias("quant_id"),
        ])
        .drop("lsv_id")
        .filter(pl.col("chromosome").is_not_null())
    )

    cols = pheno_meta.columns
    pheno_meta = pheno_meta.select(
        ["phenotype_id", "quant_id", "group_id", "gene_id"] + cols[3:13]
    )

    # Filter phenotype matrix
    pheno = pheno.join(
        pheno_meta.select("phenotype_id"),
        on="phenotype_id",
        how="inner"
    )

    print(f"Traits: {pheno.height}")
    print(f"LSVs: {pheno_meta.select(pl.col('group_id').n_unique()).item()}")

    pheno = pheno.drop(["gene_id", "lsv_id", "lsv_coord"])

    # -------------------------
    # QC filtering
    # -------------------------
    sm = pl.read_csv(sample_metadata, separator="\t")

    if(KEEP_XY):
        pheno = filter_phenotype(pheno, pheno_meta,True)
    else:
        pheno = filter_phenotype(pheno, pheno_meta)

    if(QC):
        pheno = sample_qc(pheno, sm, True, True)
    else:
        pheno = sample_qc(pheno,sm)

    #pheno = pheno.to_pandas()
    pheno.write_csv(f"{args.output_prefix}_pheno_matrix.tsv.gz", separator="\t",compression="gzip")
    pheno_meta.write_csv(f"{args.output_prefix}_majiq_metadata.tsv.gz", separator="\t",compression="gzip")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Format majiq input files.")
    parser.add_argument('-p', '--phenotype_metadata', required=True,type=str, help="phenotype_metadata file.")
    parser.add_argument('-m', '--psi_matrix', required=True,type=str, help="MAJIQ psi_matrix")
    parser.add_argument('-s', '--sample_metadata', required=True, type=str, help="Sample_metadata file.")
    parser.add_argument('-o', '--output_prefix', required=True, type=str, help="Output file prefix")
    parser.add_argument('-q', '--qc', type=str, help="")
    parser.add_argument('-x', '--keep_xy', type=str, help="")

    args = parser.parse_args()
    main(args.phenotype_metadata,args.psi_matrix,args.sample_metadata,args.qc,args.keep_xy)
