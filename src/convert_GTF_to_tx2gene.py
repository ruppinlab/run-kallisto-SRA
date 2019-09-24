from gtfparse import read_gtf

df = read_gtf(snakemake.input[0])
df = df[["transcript_id", "gene_name"]].drop_duplicates()
df = df.loc[~(df["transcript_id"] == '')]
df.to_csv(snakemake.output[0], index=False)
