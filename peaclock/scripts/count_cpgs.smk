rule ParaMethR:
    input:
        refs = config["reference_fasta"],
        reads = os.path.join(config["output_path"],"demultiplexed_reads","{barcode}.fastq"),
        cpg = config["cpg_sites"],
        matrix = config["matrix"],
    params:
        sample = "{barcode}"
    output:
        counts_long = os.path.join(config["output_path"],"counts","{barcode}.cpg_counts.csv"),
        counts_wide = os.path.join(config["output_path"],"counts","{barcode}.cpg_wide.csv")
    shell:
        """
        paramether.py \
            --reads {input.reads:q} \
            --reference {input.refs:q} \
            --cpg_csv {input.cpg:q} \
            --substitution_matrix {input.matrix:q} \
            --sample {params.sample} \
            --report {output.counts_wide:q} \
            --counts {output.counts_long:q}
        """

rule gather_reports:
    input:
        expand(os.path.join(config["output_path"],"counts","{barcode}.cpg_wide.csv"), barcode = config["barcodes"])
    output:
        os.path.join(config["output_path"],"reports","cpg_wide.csv")
    run:
        with open(output[0],"w") as fw:
            fw.write("sample,gm7_103,gm7_112,gm7_133,gm7_144,gm7_148,gm7_159,gm7_64,gm7_67,gm7_76,gm7_79,hsp4_100,hsp4_105,hsp4_120,hsp4_126,hsp4_131,hsp4_144,hsp4_70,kcns1_100,kcns1_103,kcns1_111,kcns1_116,kcns1_119,kcns1_121,kcns1_125,kcns1_134,kcns1_139,kcns1_144,kcns1_146,kcns1_162,kcns1_165,kcns1_168,kcns1_178,kcns1_30,kcns1_52,kcns1_54,kcns1_57,kcns1_65,kcns1_80,prima1_127,prima1_141,prima1_54\n")
            for fn in input:
                with open(str(fn),"r") as f:
                    for l in f:
                        fw.write(l.rstrip() + '\n')
        # shell("cat temp.txt {input} > {output} && rm temp.txt")

rule gather_counts:
    input:
        expand(os.path.join(config["output_path"],"counts","{barcode}.cpg_counts.csv"), barcode = config["barcodes"])
    output:
        os.path.join(config["output_path"],"reports","cpg_counts.csv")
    run:
        with open(output[0],"w") as fw:
            fw.write("sample,cpg_site,total,c,t,a,g,-\n")
            for fn in input:
                with open(str(fn),"r") as f:
                    for l in f:
                        fw.write(l.rstrip() + '\n')

            # shell("cat temp.txt '{input}' > '{output}' && rm temp.txt")
