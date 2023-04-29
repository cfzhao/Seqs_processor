import argparse
import subprocess
from Bio import SeqIO
import os
import re


def prepare_ref_for_hybpiper(ref_dir):
    res = []
    for i in os.listdir(ref_dir):
        if i.endswith(("fasta", "fa")):
            og_id = re.findall("OG[0-9]+", i)[0]
            og_fasta = SeqIO.to_dict(SeqIO.parse(ref_dir+"/"+i, "fasta"))
            for v in og_fasta.values():
                v.id = v.id+"-"+og_id
                v.description = ""
                res.append(v)
    commond1 = "mkdir ogs_assembly/ref_for_hybpiper"
    subprocess.run(commond1, shell=True)
    SeqIO.write(res, "ogs_assembly/ref_for_hybpiper/ref_for_hyb.fasta", "fasta")


def conbine_filtered_reads(species_name):
    commond4 = "cat ./ogs_assembly/filtered_reads/*fasta > ./ogs_assembly/filtered_reads/all_reads.fasta"
    commond5 = "sed 's/>@/>/g' -i ./ogs_assembly/filtered_reads/all_reads.fasta"
    commond6 = "seqtk seq -1 ./ogs_assembly/filtered_reads/all_reads.fasta > ogs_assembly/reads_for_hybpiper/" + species_name + "_reads1.fasta"
    commond7 = "seqtk seq -2 ./ogs_assembly/filtered_reads/all_reads.fasta > ogs_assembly/reads_for_hybpiper/" + species_name + "_reads2.fasta"
    commond8 = "seqtk seq -F '#' ogs_assembly/reads_for_hybpiper/" + species_name + "_reads1.fasta > ogs_assembly/reads_for_hybpiper/" + species_name + "_reads1.fastq"
    commond9 = "seqtk seq -F '#' ogs_assembly/reads_for_hybpiper/" + species_name + "_reads2.fasta > ogs_assembly/reads_for_hybpiper/" + species_name + "_reads2.fastq"
    commond10 = "mkdir ogs_assembly/reads_for_hybpiper"
    subprocess.run(commond10, shell=True)
    subprocess.run(commond4, shell=True)
    subprocess.run(commond5, shell=True)
    subprocess.run(commond6, shell=True)
    subprocess.run(commond7, shell=True)
    subprocess.run(commond8, shell=True)
    subprocess.run(commond9, shell=True)

def main():
    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, usage="%(prog)s [options]",
                                   description="Ogs_assembler @Biofeng.20230409")
    pars.add_argument("-r1", dest="fq_file_1", type=str,
                      help="Input file(s) with forward paired-end reads (*.fq/.gz/.tar.gz).")
    pars.add_argument("-r2", dest="fq_file_2", type=str,
                      help="Input file(s) with reverse paired-end reads (*.fq/.gz/.tar.gz).")
    pars.add_argument("-ref", dest="reference", type=str, help="Input a file(directory) with references.")
    pars.add_argument("-k1", dest="filter_kmer", type=int, help="Kmer setting for filtering reads. Default:29")
    pars.add_argument("-k2", dest="assemble_kmer", type=int, help="Kmer setting for assembling reads. Default:21, 33, 55")
    pars.add_argument("-t", dest="threads", type=int, help="Threads setting for reads enrichment")
    pars.add_argument("-n", dest="species_name", type=str, help="Species name for file naming")
    pars.add_argument("-o", dest="output_dir", type=str, help="Output dir")
    args = pars.parse_args()
    
    #commond0 = "mkdir easy353_res"
    #subprocess.run(commond0, shell=True)
    if args.filter_kmer and not args.threads:
        easy353_commond = "easy353.py -1 " + str(args.fq_file_1) + " -2 " + str(args.fq_file_2) \
                          + " -r " + str(args.reference) + " -k1 " + str(args.filter_kmer) + " -o ./ogs_assembly_res -f 1"
    elif args.threads and not args.filter_kmer:
        easy353_commond = "easy353.py -1 " + str(args.fq_file_1) + " -2 " + str(args.fq_file_2) \
                          + " -r " + str(args.reference) + " -t1 " + str(args.threads) + " -o ./ogs_assembly_res -f 1"
    elif args.filter_kmer and args.threads:
        easy353_commond = "easy353.py -1 " + str(args.fq_file_1) + " -2 " + str(args.fq_file_2) \
                          + " -r " + str(args.reference) + " -k1 " + str(args.filter_kmer) + " -t1 " + str(args.threads) + " -o ./ogs_assembly_res -f 1"
    else:
        easy353_commond = "easy353.py -1 " + str(args.fq_file_1) + " -2 " + str(args.fq_file_2) \
                          + " -r " + str(args.reference) + " -o ./ogs_assembly_res -f 1"
    subprocess.run(easy353_commond, shell=True)

    prepare_ref_for_hybpiper(str(args.reference))
    conbine_filtered_reads(str(args.species_name))

    if args.assemble_kmer:
        hybpiper_commond = "hybpiper assemble -r " + str(args.species_name) + "_reads1.fastq " \
                           + str(args.species_name)\
                           + "_reads2.fastq -t_dna ogs_assembly/ref_for_hybpiper/ref_for_hyb.fasta --hybpiper_output ogs_assembly/"\
                           + str(args.output_dir) + " --kvals " + str(args.assemble_kmer)
    else:
        hybpiper_commond = "hybpiper assemble -r " + str(args.species_name) + "_reads1.fastq " \
                           + str(args.species_name) \
                           + "_reads2.fastq -t_dna ogs_assembly/ref_for_hybpiper/ref_for_hyb.fasta --hybpiper_output ogs_assembly/" \
                           + str(args.output_dir)
    subprocess.run(hybpiper_commond, shell=True)


if __name__ == "__main__":
    main()
