import argparse
import os
import sys
import csv
from collections import Counter
import collections
import parasail
from Bio import SeqIO
from Bio import Seq

def parse_args():
    parser = argparse.ArgumentParser(description='ParaMethR')

    parser.add_argument("--reads", action="store", type=str, dest="reads")
    parser.add_argument("--references", action="store", type=str, dest="references")
    parser.add_argument("--cpg_csv", action="store", type=str, dest="cpg_csv")
    parser.add_argument("--substitution_matrix", action="store", type=str, dest="substitution_matrix")
    parser.add_argument("--report", action="store", type=str, dest="report")
    parser.add_argument("--sample", action="store", type=str, dest="sample")
    parser.add_argument("--counts", action="store", type=str, dest="counts")

    return parser.parse_args()

def get_best_reference(query, ref_dict, matrix):
    best_reference_alignment = {
        "reference": "None_NA",
        "identity": 0,
        "coverage": 0
        }
    for ref in ref_dict:
        
        new_reference_alignment = align_read(query, ref, ref_dict[ref], matrix)
        if new_reference_alignment["coverage"] > 0.7:
            if best_reference_alignment["identity"] < new_reference_alignment["identity"]: 
                best_reference_alignment = new_reference_alignment
            else:
                pass
    return best_reference_alignment

def align_read(query, ref_id, reference, matrix, gap_open=3, gap_extension=2):
    result_trace = None
    traceback = None
    result_trace = parasail.sw_trace_striped_sat(query, reference, gap_open, gap_extension, matrix)
    traceback = result_trace.get_traceback('|', '.', ' ')
    # print(traceback.ref)
    # print(traceback.comp)
    # print(traceback.query)
    query_start = result_trace.cigar.beg_query
    reference_start = result_trace.cigar.beg_ref
    # cigar = result.cigar.decode.decode("UTF-8")
    result_stats = None
    result_stats = parasail.sw_stats_striped_sat(query, reference, gap_open, gap_extension, matrix)
    alignment_covers = int(result_stats.length) / len(reference)
    # print(ref_id, len(reference), alignment_covers, result_stats.matches, len(reference))
    return {
            "reference":ref_id,
            "query_start": query_start,
            "reference_start": reference_start,
            "matches": result_stats.matches,
            "coverage": alignment_covers,
            "aln_len": result_stats.length,
            "len": result_stats.len_ref,
            "identity": result_stats.matches / result_stats.len_ref,
            "ref": traceback.ref,
            "comp": traceback.comp,
            "query": traceback.query
        }

def process_file(reads,references,cpg_dict,sample,cpg_counter,nuc_matrix):

    
    counts = Counter()
    # error = collections.defaultdict(list)

    for record in SeqIO.parse(reads, "fastq"):

        # print("*****")
        # print("the record id is",record.id)
        stats = get_best_reference(str(record.seq), references, nuc_matrix)
        best_ref,direction = stats["reference"].split("_")

        background_error_rate = get_background_error_rate(stats)

        # output_file.write('{},{},{},{},{},{},{}\n'.format(sample,best_ref,direction,stats["identity"],stats["query_start"],stats["reference_start"],background_error_rate))
        alignment_covers = int(stats["aln_len"]) / int(stats["len"]) # doesn't account for gaps so can be > 1

        if stats["identity"] > 0.75:

            read_seq = ''
            if direction == "forward":
                read_seq = str(record.seq)
            else:
                read_seq = str(record.seq.reverse_complement())
            ref_seq = references[best_ref + "_forward"]
            
            alignment = align_read(read_seq, best_ref, ref_seq, nuc_matrix)
            
            # print('\n', alignment["reference"], '\n', alignment["query"],'\n', alignment["comp"],'\n', alignment["ref"])
            # print(best_ref,'\t', direction, '\t',alignment_covers, '\t', alignment["len"], '\t',alignment["aln_len"], '\t',alignment["identity"],'\t', alignment["query_start"],'\t', alignment["reference_start"])

            for site in cpg_dict[best_ref]:

                read_variant = get_site(site[1],alignment)
                cpg_counter[site[0]][read_variant]+=1


                # if read_variant not in ["C","T","-"]:
                #     print(site)
                #     print('\n', alignment["reference"], '\n', alignment["query"],'\n', alignment["comp"],'\n', alignment["ref"])
                

            counts[best_ref]+=1

        else:
            counts["None"]+=1


    return counts, cpg_counter


def get_background_error_rate(stats):

    T_count = 0
    C_count = 0
    for i in range(int(stats["aln_len"])):
        if stats["ref"][i] == "T":
            T_count +=1
            if stats["query"][i] == 'C':
                C_count +=1

    return C_count/T_count 

def get_site(cpg_index, stats):
    adjusted_index = cpg_index - stats["reference_start"]
    current_index = 0

    variant = ""
    for i in range(int(stats["aln_len"])):
        if adjusted_index == current_index:
            # print("CPG index:", cpg_index, stats["ref"][i],stats["query"][i])
            variant = stats["query"][i]

        if stats["ref"][i] != '-':
            current_index +=1

    return variant

def load_cpg_dict(cpg_csv):
    cpg_dict= collections.defaultdict(list)
    with open(str(cpg_csv),"r") as f:
        cpg_file = csv.DictReader(f)
        for row in cpg_file:
            position = int(row["position(1-based)"]) - 1
            cpg_dict[row["gene"].lower()].append((row["gene"].lower()+ "_" + row["position(1-based)"],position))
    return cpg_dict

def make_cpg_counter(cpg_csv):
    cpg_counts= {}
    with open(str(cpg_csv),"r") as f:
        cpg_file = csv.DictReader(f)
        for row in cpg_file:
            cpg_counts[row["gene"].lower()+ "_" + row["position(1-based)"]] = Counter()
    return cpg_counts

def load_reference_dict(ref_file):

    references = {}
    for record in SeqIO.parse(ref_file, "fasta"):
         references[record.id + "_forward"] = str(record.seq)
         references[record.id + "_reverse"] = str(record.seq.reverse_complement())
    return references

if __name__ == '__main__':

    args = parse_args()
    references = load_reference_dict(args.references)
    cpg_dict = load_cpg_dict(args.cpg_csv)
    cpg_counter = make_cpg_counter(args.cpg_csv)
    fw = open(str(args.report),"w")
    fw2 = open(str(args.counts),"w")
    nuc_matrix = parasail.Matrix(str(args.substitution_matrix))

    counts, cpg_counts = process_file(str(args.reads), references, cpg_dict, args.sample, cpg_counter,nuc_matrix)

    count_str = str(args.sample) + ","

    for i in sorted(cpg_counts):


        c_and_t = cpg_counts[i]["C"] + cpg_counts[i]["T"]
        prop = "NA"
        if c_and_t > 50:
            prop = round(cpg_counts[i]["C"] / c_and_t, 3)
        count_str += f"{prop},"
        
        

        x = [j for j in cpg_counts[i]]
        y = [cpg_counts[i][j] for j in cpg_counts[i]]

        print_string = i + '\t'
        for index in range(len(x)):
            print_string += f"{x[index]}\t{y[index]}\t"

        total = sum(cpg_counts[i].values())

        c = cpg_counts[i]["C"]
        t = cpg_counts[i]["T"]
        a = cpg_counts[i]["A"]
        g= cpg_counts[i]["G"]
        gap = cpg_counts[i]["-"]
        fw2.write(f"{args.sample},{i},{total},{c},{t},{a},{g},{gap}\n")

    count_str = count_str.rstrip(',')
    fw.write(count_str+'\n')

    fw.close()
    fw2.close()