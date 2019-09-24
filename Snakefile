from os.path import join
import pandas as pd

# URL for downloading transcript sequences file
GENCODE_URL="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/GRCh37_mapping/gencode.v31lift37.transcripts.fa.gz"
GENCODE_GTF_URL="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/GRCh37_mapping/gencode.v31lift37.annotation.gtf.gz"

# Directories
OUTPUT_DIR="output"
DATA_DIR="data"
FASTQ_DIR=join(OUTPUT_DIR, "FASTQ")
SAMPLE_OUTPUT_DIR=join(OUTPUT_DIR, "{sample}")

# Input Files
SRA_FILE=join("..", "raw_data", "PRJEB23709", "sra", "{sample}.sra")
SAMPLES_FILE=join(DATA_DIR, "SRR_Acc_list.txt")

# Intermediate Files
TRANSCRIPT_SEQ_FILE=join(DATA_DIR, "ref_transcriptome.fq.gz")
GENCODE_GTF_GZ=join(DATA_DIR, "gencode_comprehensive-annotation-GRCh37.gtf.gz")
GENCODE_GTF=join(DATA_DIR, "gencode_comprehensive-annotation-GRCh37.gtf")
TX2GENE_FILE=join(DATA_DIR, "tx2gene-GRCh37.csv")
INDEX_FILE=join(DATA_DIR, "GRCh37.idx")
RAW_FASTQ1_FILE=join(FASTQ_DIR, "{sample}_1.fastq")
RAW_FASTQ2_FILE=join(FASTQ_DIR, "{sample}_2.fastq")

# Output Files
H5_ABUNDANCE_FILE = join(SAMPLE_OUTPUT_DIR, "abundance.h5")
TPM_FILE=join(OUTPUT_DIR, "TPM.tsv")
df = pd.read_csv(SAMPLES_FILE, sep="\t")

rule all:
    input:
        TPM_FILE
        # expand(H5_ABUNDANCE_FILE, sample=df["sample"])

# rule for combining kallisto
rule get_TPM_kallisto:
    conda:
        "envs/process-mRNA.yml"
    input:
        TX2GENE_FILE,
        SAMPLES_FILE
        # expand(H5_ABUNDANCE_FILE, sample=df["sample"])
    output:
        TPM_FILE
    script:
        "src/combine_kallisto_output.R"

rule extract_tx2gene:
    conda:
        "envs/process-mRNA.yml"
    input:
        GENCODE_GTF
    output:
        TX2GENE_FILE
    script:
        "src/convert_GTF_to_tx2gene.py"

rule decompress_GTF_file:
    input:
        GENCODE_GTF_GZ
    output:
        GENCODE_GTF
    shell:
        "gzip -d {input}"

rule download_GTF_file:
    params:
        GENCODE_GTF_URL
    output:
        GENCODE_GTF_GZ
    shell:
        "wget {params} -O {output}"

# rules for running kallisto
rule process_fastq:
    input:
        fq1=RAW_FASTQ1_FILE,
        fq2=RAW_FASTQ2_FILE,
        i=INDEX_FILE
    output:
        H5_ABUNDANCE_FILE
    params:
        SAMPLE_OUTPUT_DIR
    shell:
        "module load kallisto; kallisto quant -i {input.i} -o {params} -t 8 -b 100 {input.fq1} {input.fq2}"

# rule for downloading FASTQ files from SRA
# rule download_SRA_files:
#     output:
#         SRA_FILE
#     shell:
#         "module load sratoolkit; prefetch {sample}"

rule extract_FASTQ_from_SRA:
    params:
        FASTQ_DIR
    output:
        RAW_FASTQ1_FILE,
        RAW_FASTQ2_FILE
    shell:
        "module load sratoolkit; fastq-dump -O {params} -I --split-files {sample}"

# rules for preparing kallisto index
rule download_gencode_GRCh38p12_transcript_sequences:
    params:
        GENCODE_URL
    output:
        TRANSCRIPT_SEQ_FILE
    shell:
        "wget {params} -O {output}"


rule create_kallisto_index:
    input:
        TRANSCRIPT_SEQ_FILE
    output:
        INDEX_FILE
    shell:
        "module load kallisto; kallisto index -i {output} {input}"
