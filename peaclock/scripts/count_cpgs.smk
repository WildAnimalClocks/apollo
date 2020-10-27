import os
from Bio import SeqIO
import peaclockfunks as qcfunk


rule all:
    input:
        os.path.join(config["outdir"], "reports","cpg_counts.csv"),
        os.path.join(config["outdir"], "reports","cpg_wide.csv"),
        expand(os.path.join(config["outdir"],"gathered_reads","{barcode}.fastq"), barcode=config["barcodes"])

rule gather_demuxed_reads:
    input:
    params:
        barcode_path = os.path.join(config["outdir"],"demultiplexed_reads","{barcode}"),
        barcode = "{barcode}"
    output:
        reads = os.path.join(config["outdir"],"gathered_reads","{barcode}.fastq")
    run:
        with open(output.reads, "w") as fw:
            reads = []
            for r,d,f in os.walk(params.barcode_path):
                for fn in f:
                    if fn.endswith(".fastq") or fn.endswith(".fq"):
                        for record in SeqIO.parse(os.path.join(r, fn),"fastq"):
                            if len(record) > config["min_length"] and len(record) < config["max_length"]:
                                reads.append(record)
            SeqIO.write(reads,fw,"fastq")
            print(qcfunk.green(f"Barcode {params.barcode}: ") + f"{len(reads)}")

rule paramether:
    input:
        reads = rules.gather_demuxed_reads.output.reads,
        refs = config["genes"],
        cpg_sites = config["cpg_sites"],
        matrix_file = config["matrix_file"],
    params:
        sample = "{barcode}"
    output:
        counts_long = os.path.join(config["outdir"],"counts","{barcode}.cpg_counts.csv"),
        counts_wide = os.path.join(config["outdir"],"counts","{barcode}.cpg_wide.csv")
    shell:
        """
        paramether.py \
            --reads {input.reads:q} \
            --reference {input.genes:q} \
            --cpg_csv {input.cpg_sites:q} \
            --substitution_matrix {input.matrix:q} \
            --sample {params.sample} \
            --report {output.counts_wide:q} \
            --counts {output.counts_long:q}
        """

rule gather_reports:
    input:
        expand(os.path.join(config["outdir"],"counts","{barcode}.cpg_wide.csv"), barcode = config["barcodes"])
    output:
        os.path.join(config["outdir"],"reports","cpg_wide.csv")
    run:

        header = qcfunk.make_cpg_header(config['cpg_csv'])
        with open(output[0],"w") as fw:
            fw.write(f"{header}\n")
            # fw.write("sample,gm7_103,gm7_112,gm7_133,gm7_144,gm7_148,gm7_159,gm7_64,gm7_67,gm7_76,gm7_79,hsp4_100,hsp4_105,hsp4_120,hsp4_126,hsp4_131,hsp4_144,hsp4_70,kcns1_100,kcns1_103,kcns1_111,kcns1_116,kcns1_119,kcns1_121,kcns1_125,kcns1_134,kcns1_139,kcns1_144,kcns1_146,kcns1_162,kcns1_165,kcns1_168,kcns1_178,kcns1_30,kcns1_52,kcns1_54,kcns1_57,kcns1_65,kcns1_80,prima1_127,prima1_141,prima1_54\n")
            for fn in input:
                with open(str(fn),"r") as f:
                    for l in f:
                        fw.write(l.rstrip() + '\n')
        # shell("cat temp.txt {input} > {output} && rm temp.txt")

rule gather_counts:
    input:
        expand(os.path.join(config["outdir"],"counts","{barcode}.cpg_counts.csv"), barcode = config["barcodes"])
    output:
        os.path.join(config["outdir"],"reports","cpg_counts.csv")
    run:
        with open(output[0],"w") as fw:
            fw.write("sample,cpg_site,total,c,t,a,g,-\n")
            for fn in input:
                with open(str(fn),"r") as f:
                    for l in f:
                        fw.write(l.rstrip() + '\n')

            # shell("cat temp.txt '{input}' > '{output}' && rm temp.txt")


