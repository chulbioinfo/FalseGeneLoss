# Python 3.6.9
# Copyright 2020 Juwan Kim

"""
[FalseGeneLoss.py]
Version 1.0
    (1) OS: Ubuntu 18.04.3 LTS (GNU/Linux 5.4.0-42-generic x86_64)
    (2) Necessary input files:
        - Assembly A from NCBI (.fasta, should be softmasked by windowmasker)
        - Assembly A's annotation from NCBI (.gff)
        - Assembly B (.fasta, should be softmasked by windowmasker)
        - Assembly B's annotation from CAT (.gff, Comparative annotation toolkit; https://github.com/ComparativeGenomicsToolkit/Comparative-Annotation-Toolkit)
        - Hal file from the cactus alignment between Assembly A & B (.hal)
        - Bam file of assembly A with illumina genomic reads (.bam, should be sorted and indexed)
        - Bam file of assembly B with illumina genomic reads (.bam, should be sorted and indexed)
        - NCBI gap coordinate file of assembly A
        - NCBI gap coordinate file of assembly B
        - CAT working directory and out directory
        - VGP02All_merged.merged_by_bedtools.bed --> Manual BED file for the duplicated chr 29 segments of zebra finch

        [optional]
        - (if want to perform false indel/SNP analyses) false indel/SNP bed file based on assembly A
        - (if want to compare NCBI remap result with cactus) NCBI remap excel file

    (3) Related program version
         - paftools version: included in minimap2 2.17-r974-dirty
         - minimap2 version: 2.17-r974-dirty
         - bedtools version: v2.27.0
         - samtools version: 1.9-177-g796cf22
         - halLiftover version: halLiftover v2.1
         - newcpgreport version: EMBOSS:6.6.0.0
         - makeblastdb & blastn version: 2.6.0+
         - blat version: v. 36x2

    (4) variable_suffix:
        [variable name]_s: string
        [variable name]_chr: character (=length 1 string)
        [variable name]_l: list
        [variable name]_nl: nested list
        [variable name]_i: integer
        [variable name]_bool: boolean
        [variable name]_d: dictonary
        [variable name]_f: file
        [variable name]_fl: float
        [variable name]_rec: gffutils record

    (5) Genome assembly naming: [GenomeVersion]_[SpeciesName]
    Genome version: vgp for VGP genome
                    pre for Previous genome
    Species name: taegut for zebra finch
                  calann for Anna's hummingbird
                  ornana for Platypus
                  anates for Climbing perch
"""

# used modules
import gffutils, os, sqlite3, re, argparse, sys, pickle, subprocess, copy, openpyxl
from Bio.SeqIO import FastaIO
from multiprocessing import Pool
from itertools import islice

# global commands
paftools_command_s = "PATH_TO_PAFTOOLS"  # version: included in minimap2 2.17-r974-dirty
minimap2_command_s = "PATH_TO_MINIMAP2"  # version: 2.17-r974-dirty
bedtools_command_s = "PATH_TO_BEDTOOLS"  # version: v2.27.0
samtools_command_s = "PATH_TO_SAMTOOLS"  # version: 1.9-177-g796cf22
halLiftover_commnd_s = "PATH_TO_HALLIFTOVER"  # version: halLiftover v2.1:
newcpgreport_command_s = "PATH_TO_NEWCPGREPORT"  # version: EMBOSS:6.6.0.0
makeblastdb_command_s = "PATH_TO_MAKEBLASTDB"  # version: 2.6.0+
blastn_command_s = "PATH_TO_BLASTN"  # version: 2.6.0+
blat_command_s = "PATH_TO_BLAT"  # version: v. 36x2

home_dir_s = os.path.join("PATH", "TO", "THE", "RESULT_DIR")

# global function

def db_maker(gff_s, name_s, working_dir_s):
    """
    function to return pre-built annotation DB or newly generate it if not exist
    :param gff_s: string, path to gff file
    :param name_s: string, name of annotation DB
    :param working_dir_s: string, path to the directory a the DB would be made
    :return: gffutils DB
    """
    annot_db_s = os.path.join(working_dir_s, "%s.db" % name_s)
    # (0) return pre-bulit annotation DB
    if file_check(annot_db_s) is True:
        db = gffutils.FeatureDB(dbfn=annot_db_s)
    # (1) if there's no pre-bulit annotation DB:
    else:
        assert file_check(annot_db_s) is False
        dir_check(working_dir_s)
        # (1-1) try to generate new DB
        try:
            db = gffutils.create_db(data=gff_s, dbfn=annot_db_s, force=True, keep_order=True,
                                    merge_strategy='create_unique', sort_attribute_values=True)
        # (1-2) or generate Error
        except FileNotFoundError:
            assert file_check(gff_s) is False
            print("[FileNotFoudError] The path %s seems not found. Please check the path again. Exit..." % gff_s)
            sys.exit()
    return db


def dir_check(path_s):
    """
    function to check if given directory exists or generate it
    :param path_s: string, path to directory
    :return: none
    """
    if os.path.isdir(path_s) is False:
        os.makedirs(path_s)
        print("[PATH] %s is made" % path_s)


def file_check(path_s):
    """
    function to check whether given file exists
    :param path_s: string, path to file
    :return: none
    """
    if os.path.isfile(path_s) is False or os.stat(path_s).st_size == 0:
        return False
    else:
        return True


def concatenate_fasta(input_path_s, output_path_s, header_s):
    """
    function to concatenate multifasta file into single sequence fasta file
    :param input_path_s: string, path to input fasta file
    :param output_path_s: string, path to output concatenated fasta file
    :return: X
    """
    concatenated_seq_s = ""
    with open(input_path_s) as input_path_f:
        with open(output_path_s, "w") as output_path_f:
            count_i = 0
            for values in FastaIO.SimpleFastaParser(input_path_f):
                seq_header_s, seq_s = values[0], values[1]
                concatenated_seq_s += seq_s
                count_i += 1
                if divmod(count_i, 10000)[1] == 0:
                    print(len(concatenated_seq_s))
            idx_i = 0
            chunk_size_i = 100000000
            chunk_seq_l = [concatenated_seq_s[y - chunk_size_i:y] for y in range(chunk_size_i, len(concatenated_seq_s) + chunk_size_i, chunk_size_i)]
            for chunk_seq_s in chunk_seq_l:
                idx_i += 1
                output_path_f.write(">%s_%d\n%s\n" % (header_s, idx_i, chunk_seq_s))
    genome_bedfile_s = os.path.join(result_dir, "basic_inputs", "%s.genome.bed" % header_s)
    os.system("%s faidx %s" % (samtools_command_s, output_path_s))
    os.system("awk -F'[\"\\t\"]' '{OFS=\"\\t\"} {print $1,0,$2,$1,0,\".\"}' %s > %s" % (output_path_s + ".fai", genome_bedfile_s))


def window_level_gc_n_repeat(window_size_i, genome_bed_file_s, softmasked_bed_file_s, mode_s, fasta_file_s, gc_repeat_window_file_s):
    """
    function to calculate GC and repeat content with given window size (non-overlapping)
    :param window_size_i: integer, size of window
    :param genome_bed_file_s: string, path to genome bed file
    :param softmasked_bed_file_s: string, path to softmasking bed file
    :param mode_s: string, mode (ref or tgt or missing or aligned)
    :param fasta_file_s: string, path to fasta file
    :param gc_repeat_window_file_s: string, path to gc and repeat window file
    :return: X
    """
    os.system("%s makewindows -b %s -w %d |" % (bedtools_command_s, genome_bed_file_s, window_size_i) +  # generate windows
              "awk -F'[\"\\t\"]' '{OFS=\"\\t\"} {if ($3-$2 == %d) {print $0}}'|" % window_size_i +  # remove remnants
              "%s sort -i stdin |" % bedtools_command_s +
              "%s intersect -wao -a stdin -b %s |" % (bedtools_command_s, softmasked_bed_file_s) +  # calculate repeat ratio
              "%s sort -i stdin |" % bedtools_command_s +
              "%s groupby -i stdin -g 1,2,3 -c 7 -o sum | " % bedtools_command_s +
              "awk -F'[\"\\t\"]' '{OFS=\"\\t\"} {print $1, $2, $3, \"%s|\"$4*100/%d, 0, \".\"}'|" % (mode_s, window_size_i) +
              "%s sort -i stdin |" % bedtools_command_s +
              "%s nuc -fi %s -bed stdin | grep -v \"#1_usercol\" |" % (bedtools_command_s, fasta_file_s) +  # calculate gc ratio
              "awk -F'[\"\\t\",\"|\"]' '{OFS=\"\\t\"} {if (($10+$11+$12+$13) >= %d) {print $1, $4, $5, ($11+$12)/%d*100, ($10+$11+$12+$13)}}' > %s" %
              (window_size_i, window_size_i, gc_repeat_window_file_s))


def aln_id_2_gene_id_dictionary():
    """
    function to return dict which links alingment ID (transcript ID) to the name of gene in the projected annotation DB
    :return: dict, aln_id_2_gene_id_d, key: aln_id (e.g. XM_..), value: gene name (gene-XXX)
             dict. gene_id_2_aln_id_d, key: gene name (gene-XXX), value: aln_id (e.g. XM_..)
    """
    aln_2_gene_d = {}
    gene_2_aln_d = {}
    for trans_rec in tgt_all_trans_recs_l:
        parent_gene_rec = [gene_rec for gene_rec in tgtDB.parents(id=trans_rec.id, featuretype="gene")][0]
        parent_gene_name_s = get_new_name(gene_rec=parent_gene_rec, mode_s="tgt")
        aln_id_s = trans_rec.attributes['alignment_id'][0]
        aln_2_gene_d[aln_id_s] = parent_gene_name_s
        gene_2_aln_d[parent_gene_name_s] = aln_id_s
    return aln_2_gene_d, gene_2_aln_d


def gene_name_2_gene_id():
    """
    function to return dict which links gene ID to the name of gene in annotation DBs
    :return: dict (nested), gene_name_2_gene_id_d,
                            key: mode (ref: VGP, tgt: OLD), value: dict with key:gene name, value:gene id
             dict (nested), gene_id_2_gene_name_d,
                            key: mode (ref: VGP, tgt: OLD), value: dict with key:gene id, value:gene name
    """
    # pickle files to save the dict (gene_name -- gene_id)
    gene_name_2_gene_id_file_s = os.path.join(result_dir, "gene_name_2_gene_id_d.pickle")
    gene_id_2_gene_name_file_s = os.path.join(result_dir, "gene_id_2_gene_name_d.pickle")
    try:
        with open(gene_name_2_gene_id_file_s, "rb") as gene_name_2_gene_id_file_f,\
                open(gene_id_2_gene_name_file_s, "rb") as gene_id_2_gene_name_file_f:
            gene_name_2_id_d = pickle.load(gene_name_2_gene_id_file_f)
            gene_id_2_name_d = pickle.load(gene_id_2_gene_name_file_f)
    except FileNotFoundError:
        gene_name_2_id_d = {"ref": {}, "tgt": {}}
        gene_id_2_name_d = {"ref": {}, "tgt": {}}
        # gene record files (mode: ref (VGP) or tgt (Prior))
        for gene_rec in ref_all_gene_recs_l:
            gene_id_s = gene_rec.id
            gene_name_s = get_new_name(gene_rec=gene_rec, mode_s="ref")
            gene_name_2_id_d["ref"][gene_name_s] = gene_id_s
            gene_id_2_name_d["ref"][gene_id_s] = gene_name_s
        for gene_rec in tgt_all_gene_recs_l:
            gene_id_s = gene_rec.id
            gene_name_s = get_new_name(gene_rec=gene_rec, mode_s="tgt")
            gene_name_2_id_d["tgt"][gene_name_s] = gene_id_s
            gene_id_2_name_d["tgt"][gene_id_s] = gene_name_s
        with open(gene_name_2_gene_id_file_s, "wb") as gene_name_2_gene_id_file_f,\
                open(gene_id_2_gene_name_file_s, "wb") as gene_id_2_gene_name_file_f:
            pickle.dump(gene_name_2_id_d, gene_name_2_gene_id_file_f, pickle.HIGHEST_PROTOCOL)
            pickle.dump(gene_id_2_name_d, gene_id_2_gene_name_file_f, pickle.HIGHEST_PROTOCOL)
    return gene_name_2_id_d, gene_id_2_name_d


def scaffold_2_size_dictionary(mode_s):
    """
    function to return dict with scaffold size from fasta index file
    :param mode_s: string, "ref" for VGP, "tgt" for prior assembly
    :return: dict, chrom_2_size_d, key: scaffold name (string), value: scaffold size (int)
    """
    chrom_2_size_d = {}
    if mode_s == "ref":
        faidx_file_s = ref_fasta_file_s + ".fai"
    else:
        assert mode_s == "tgt"
        faidx_file_s = tgt_fasta_file_s + ".fai"
    if file_check(faidx_file_s) is False:
        os.system("%s faidx %s" % (samtools_command_s, faidx_file_s))
    with open(faidx_file_s) as faidx_file_f:
        for line_s in faidx_file_f:
            if len(line_s) < 1:
                continue
            line_l = line_s.rstrip("\n").split("\t")
            chrom_s, size_i = line_l[0], int(line_l[1])
            chrom_2_size_d[chrom_s] = size_i
    return chrom_2_size_d


def blast_2_bed(gene_name_s):
    """
    function to read blast output & reformat into bed-style string
    :param gene_name_s: string, name of gene
    :return: out_bedline_l: list, list of bed-style string of blast hits
    """
    blast_file_s = os.path.join(result_dir, "missingGenes", "BLAST", gene_name_s, "blast_%s.tabs" % gene_name_s)
    out_bedline_l = []
    hit_count_i = 0
    with open(blast_file_s, "r") as blast_file_f:
        for line_s in blast_file_f:
            hit_count_i += 1
            line_l = line_s.rstrip("\n").split("\t")
            query_id_s = "%s^%s^hit_%d" % (gene_name_s, line_l[0], hit_count_i)
            blast_chrom_s = line_l[1]
            blast_start_i, blast_end_i = sorted([int(line_l[8]) - 1, int(line_l[9])])
            blast_bedline_s = "\t".join(map(str, [blast_chrom_s, blast_start_i, blast_end_i, query_id_s, "0", "."]))
            out_bedline_l.append(blast_bedline_s)
    return out_bedline_l


def duplication_gene_list(duplication_gene_file_s):
    """
    return duplicated gene list
    :param duplication_gene_file_s: path to duplicated gene list file
    :return: list of duplicated genes
    """
    dup_gene_ids_l = []
    with open(duplication_gene_file_s, "r") as duplication_gene_file_f:
        for line_s in duplication_gene_file_f:
            dup_gene_id_s = line_s.rstrip("\n")
            dup_gene_ids_l.append(dup_gene_id_s)
    return dup_gene_ids_l


def gene_and_transcript_records(mode_s, dup_gene_ids_l):
    """
    function to return all gene and transcript records from annotation DB
    :param mode_s: string, ref for VGP, tgt for Projected
    :return: list, gene_recs_l (gene records), transcript_recs_l (transcirpt records)
    """
    gene_recs_l_file_s = os.path.join(result_dir, "%s.gene_recs_l.pickle" % mode_s)
    transcript_recs_l_file_s = os.path.join(result_dir, "%s.transcript_recs_l.pickle" % mode_s)
    try:  # try to get pre_stored pickles
        with open(gene_recs_l_file_s, "rb") as gene_recs_l_file_f:
            gene_recs_l = pickle.load(gene_recs_l_file_f)
        with open(transcript_recs_l_file_s, "rb") as transcript_recs_l_file_f:
            transcript_recs_l = pickle.load(transcript_recs_l_file_f)
    except FileNotFoundError:  # generate pickles + return list
        gene_recs_l, transcript_recs_l = [], []
        transcript_ids_l = []
        duplicated_transcript_ids_l = []
        gene_rec_2_transcript_id_d = {}
        # (1-1-1) from input files
        for dup_gene_id_s in dup_gene_ids_l:
            for dup_mRNA_rec in refDB.children(id=dup_gene_id_s, featuretype="mRNA"):
                dup_trans_id_s = dup_mRNA_rec.attributes["transcript_id"][0]
                print(dup_trans_id_s)
                duplicated_transcript_ids_l.append(dup_trans_id_s)
                transcript_ids_l.append(dup_trans_id_s)
        # (1) VGP
        if mode_s == "ref":
            # (1-1) collect duplicated transcript ID
            # (1-1-2) from symbols
            for gene_rec in refDB.all_features(featuretype="gene"):
                for mRNA_rec in refDB.children(id=gene_rec.id, featuretype="mRNA"):
                    trans_id = mRNA_rec.attributes["transcript_id"][0]
                    gene_rec_2_transcript_id_d[gene_rec] = trans_id
                    # check duplicated transcript
                    if trans_id in transcript_ids_l:
                        print(str(gene_rec))
                        duplicated_transcript_ids_l.append(trans_id)
                    transcript_ids_l.append(trans_id)
            # (1-2) exclude duplicated transcript ID
            for gene_rec in refDB.all_features(featuretype="gene"):
                trans_id = gene_rec_2_transcript_id_d[gene_rec]
                if trans_id not in duplicated_transcript_ids_l:
                    gene_recs_l.append(gene_rec)
            for gene_rec in gene_recs_l:
                transcript_recs_l.append([mRNA_rec for mRNA_rec in refDB.children(id=gene_rec.id, featuretype="mRNA")][0])
        # (2) Old
        else:
            assert mode_s == "tgt"
            # (2-1) collect duplicated transcript ID
            for gene_rec in tgtDB.all_features(featuretype="gene"):
                for mRNA_rec in tgtDB.children(id=gene_rec.id, featuretype="transcript"):
                    trans_id = mRNA_rec.attributes["source_transcript_name"][0]
                    gene_rec_2_transcript_id_d[gene_rec] = trans_id
                    # check duplicated transcript
                    if trans_id in transcript_ids_l:
                        print(str(gene_rec))
                        duplicated_transcript_ids_l.append(trans_id)
                    transcript_ids_l.append(trans_id)
            # (2-2) collect duplicated transcript ID
            for gene_rec in tgtDB.all_features(featuretype="gene"):
                trans_id = gene_rec_2_transcript_id_d[gene_rec]
                if trans_id not in duplicated_transcript_ids_l:
                    gene_recs_l.append(gene_rec)
            for gene_rec in gene_recs_l:
                transcript_recs_l.append([mRNA_rec for mRNA_rec in tgtDB.children(id=gene_rec.id, featuretype="transcript")][0])
        with open(gene_recs_l_file_s, "wb") as gene_recs_l_file_f:
            pickle.dump(gene_recs_l, gene_recs_l_file_f, pickle.HIGHEST_PROTOCOL)
        with open(transcript_recs_l_file_s, "wb") as transcript_recs_l_file_f:
            pickle.dump(transcript_recs_l, transcript_recs_l_file_f, pickle.HIGHEST_PROTOCOL)
    return gene_recs_l, transcript_recs_l


def make_blast_db(fasta_s, database_s):
    """
    function to make blast database
    :param fasta_s: string, fasta filepath to generate blastDB
    :param database_s: string, name of blastDB
    :return: none
    """
    os.system("%s -dbtype nucl -in %s -out %s" % (makeblastdb_command_s, fasta_s, database_s))


def exon_blast(info_l):
    """
    function to do blast
    :param info_l: [geneName_s,tgtDB]
    :return:
    """
    gene_rec, blast_db_path_s = info_l
    gene_name_s = get_new_name(gene_rec=gene_rec, mode_s="ref")
    dir_check(os.path.join(result_dir, "missingGenes", "BLAST", gene_name_s))
    query_fasta_file_s = os.path.join(result_dir, "missingGenes", "BLAST", gene_name_s, "%s.fasta" % gene_name_s)
    exon_fasta_file_s = os.path.join(result_dir, "CDS", "ref", gene_name_s, "EXON_%s.fasta" % gene_name_s)
    blast_result_file_s = os.path.join(result_dir, "missingGenes", "BLAST", gene_name_s, "blast_%s.tabs" % gene_name_s)
    with open(exon_fasta_file_s,"r") as exon_fasta_file_f, open(query_fasta_file_s,"w") as query_fasta_file_f:
        for values in FastaIO.SimpleFastaParser(exon_fasta_file_f):
            new_header_s, seq_s = values[0].split("_")[0].replace("exon","EXON_"), values[1]
            query_fasta_file_f.write(">%s\n%s\n" % (new_header_s, seq_s))
    os.system("%s -db %s -query %s -task blastn -perc_identity 90 -qcov_hsp_perc 50 -dust no -outfmt '6 std qcovus' -out %s" %
              (blastn_command_s, blast_db_path_s, query_fasta_file_s, blast_result_file_s))


def blat_preprocessing(trans_rec):
    """
    preprocessing blat inputs: db and query for each transcript
    :param trans_rec: record, transcript record (from gffutils)
    :return: X, just generate db and query file of each transcript
    """
    parent_gene_rec = [gene_rec for gene_rec in tgtDB.parents(id=trans_rec.id, featuretype="gene")][0]
    aln_id_s = trans_rec.attributes['alignment_id'][0]
    ref_trans_id_s = trans_rec.attributes['source_transcript'][0]
    parent_gene_name_s = get_new_name(gene_rec=parent_gene_rec, mode_s="tgt")

    dir_check(os.path.join(result_dir, "blat_dna", "db", parent_gene_name_s))
    dir_check(os.path.join(result_dir, "blat_dna", "query", parent_gene_name_s))
    blat_db_s = os.path.join(result_dir, "blat_dna", "db", parent_gene_name_s, "%s.fasta" % parent_gene_name_s)
    blat_query_s = os.path.join(result_dir, "blat_dna", "query", parent_gene_name_s, "%s.fasta" % parent_gene_name_s)
    ref_cdsfile_s = os.path.join(result_dir, "CDS", "ref", parent_gene_name_s, "merged_%s.fasta" % parent_gene_name_s)
    tgt_cdsfile_s = os.path.join(result_dir, "CDS", "tgt", parent_gene_name_s, "merged_%s.fasta" % parent_gene_name_s)

    is_ref_multiple_of_three_bool = True
    # VGP CDS --> DB, Prev CDS --> Query
    with open(ref_cdsfile_s) as ref_cdsfile_f, open(tgt_cdsfile_s) as tgt_cdsfile_f:
        with open(blat_db_s, "w") as blat_db_f, open(blat_query_s, "w") as blat_query_f:
            for values in FastaIO.SimpleFastaParser(ref_cdsfile_f):
                header_s, ref_cds_seq_s = ref_trans_id_s, values[1]
                try:
                    assert divmod(len(ref_cds_seq_s), 3)[1] == 0
                # if VGP CDS sequence is not multiple of three ==> Excluded from our analysis
                except AssertionError:
                    print(parent_gene_name_s)
                    is_ref_multiple_of_three_bool = False
                blat_db_f.write(">%s\n%s" % (header_s, ref_cds_seq_s))
            for values in FastaIO.SimpleFastaParser(tgt_cdsfile_f):
                header_s, tgt_cds_seq_s = aln_id_s, values[1]
                blat_query_f.write(">%s\n%s" % (header_s, tgt_cds_seq_s))
    if is_ref_multiple_of_three_bool is False:
        os.system("rm %s" % blat_db_s)
        os.system("rm %s" % blat_query_s)
    return parent_gene_name_s, is_ref_multiple_of_three_bool


def blat(trans_rec):
    """
    function to do blat with tgtCDS as query and refCDS as db for each gene, DNA level
    except dna level, it's the same option with blat in CAT
    :param trans_rec: record, transcript record (from gffutils)
    :return: X
    """
    parent_gene_rec = [gene_rec for gene_rec in tgtDB.parents(id=trans_rec.id, featuretype="gene")][0]
    parent_gene_name_s = get_new_name(gene_rec=parent_gene_rec, mode_s="tgt")
    blat_db_s = os.path.join(result_dir, "blat_dna", "db", parent_gene_name_s, "%s.fasta" % parent_gene_name_s)
    blat_query_s = os.path.join(result_dir, "blat_dna", "query", parent_gene_name_s, "%s.fasta" % parent_gene_name_s)
    blat_dna_pslfile_s = os.path.join(result_dir, "blat_dna", "out", "%s.psl" % parent_gene_name_s)
    os.system("%s -noHead -minIdentity=0 %s %s %s" % (blat_command_s, blat_db_s, blat_query_s, blat_dna_pslfile_s))


def ignore_1_or_2bp_introns(raw_rec_l):
    """
    function to remove the false 1,2bp introns because of assembly errors
    :param raw_rec_l: list, element: gffutils record
    :return: new_starts_l, new_ends_l: list, element: integer, coordinates of records starts & ends
    """
    new_starts_l, new_ends_l = [], []
    raw_rec_d = {0: ["default", "default"]}
    # (1) find less than 10bp intron
    for i in range(0, len(raw_rec_l) - 1):
        if i not in raw_rec_d:
            raw_rec_d[i] = ["default", "default"]
        if i + 1 not in raw_rec_d:
            raw_rec_d[i + 1] = ["default", "default"]
        # ...---[ REC i ]--------------[ REC i+1 ]---...
        #               ^              ^
        #        raw_rec[i].end   raw_rec[i+1].start
        #        [d,merge_start]     [merge_end,d]
        # **** if interval is less that 10bp: start merging ****
        if raw_rec_l[i].end + 11 >= raw_rec_l[i + 1].start:  # if length of interval between CDSs is less than 10bp
            raw_rec_d[i][1] = "merge_start"
            raw_rec_d[i + 1][0] = "merge_end"
    # (2) skip the false introns and annotate new recs
    # in this way, we can merge more than three falsely-splited records as well
    for i in range(0, len(raw_rec_l)):
        if raw_rec_d[i][0] == "default":
            start_i = raw_rec_l[i].start
            new_starts_l.append(start_i)
        if raw_rec_d[i][1] == "default":
            end_i = raw_rec_l[i].end
            new_ends_l.append(end_i)
    assert len(new_starts_l) == len(new_ends_l)
    assert (new_starts_l == sorted(new_starts_l)) and (new_ends_l == sorted(new_ends_l))
    assert len(new_starts_l) > 0
    return new_starts_l, new_ends_l


def cpg_island_prediction():
    """
    function to perform CpG island prediction (bed format) by newcpgreport in EMBOSS & format it to BED
    :return: X
    """
    cgi_embl_file_s = os.path.join(species_fasta_dir_s, "%s.newcpgreport" % ref_s)
    out_bed_file_s = os.path.join(species_fasta_dir_s, "%s.newcpgreport.bed" % ref_s)
    # (1) Performing CGI prediction by newcpgreport
    if file_check(cgi_embl_file_s) is False:
        os.system("%s -sequence %s -outfile %s" % (newcpgreport_command_s, ref_fasta_file_s, cgi_embl_file_s))
    # (2) Reformat EMBL file to BED manually
    with open(cgi_embl_file_s) as cgi_embl_file_f, open(out_bed_file_s, "w") as out_bed_file_f:
        embl_contents_l = cgi_embl_file_f.read().split("\n//\n")
        for cpgreport_s in embl_contents_l:
            cpgreport_l = cpgreport_s.split("\nXX\n")
            try:
                assert cpgreport_l[0].startswith("ID")
                chrom_s = cpgreport_l[0].split(" ")[3]
                print(chrom_s)
            except AssertionError:
                assert len(cpgreport_l) == 1
                continue
            assert cpgreport_l[1].startswith("DE")
            assert cpgreport_l[2].startswith("CC")
            assert cpgreport_l[3].startswith("FH")
            assert len(cpgreport_l) == 4
            cpg_infos_l = cpgreport_l[3].split("\nFT   CpG island       ")[1:]
            count_i = 0
            for cpg_info_s in cpg_infos_l:
                count_i += 1
                cpg_info_line_l = cpg_info_s.split("\nFT                    /")
                print(cpg_info_line_l)
                assert ".." in cpg_info_line_l[0]
                start_i, end_i = [int(x) for x in cpg_info_line_l[0].split("..")]
                assert start_i <= end_i
                assert "size=" in cpg_info_line_l[1]
                size_i = int(cpg_info_line_l[1].split("=")[-1])
                assert size_i == end_i-start_i+1
                assert "Sum C+G=" in cpg_info_line_l[2]
                sum_c_g_i = int(cpg_info_line_l[2].split("=")[-1])
                assert "Percent CG=" in cpg_info_line_l[3]
                percent_gc_fl = float(cpg_info_line_l[3].split("=")[-1])
                assert "ObsExp=" in cpg_info_line_l[4]
                # for last line
                if "numislands" in cpg_info_line_l[4]:
                    obs_exp_s = cpg_info_line_l[4].split("\nFT   numislands       ")[0]
                    num_islands_i = int(cpg_info_line_l[4].split("\nFT   numislands       ")[-1])
                    assert len(cpg_infos_l) == num_islands_i
                else:
                    obs_exp_s = cpg_info_line_l[4]
                obs_exp_fl = float(obs_exp_s.split("=")[-1])
                out_bed_file_f.write("\t".join(map(str,[chrom_s, start_i-1, end_i, "CpGIsland_%s_%i|size=%i|Sum C+G=%i|Percent CG=%f|ObsExp=%s" %
                                                        (chrom_s, count_i, size_i, sum_c_g_i, percent_gc_fl, obs_exp_fl), 0, "+"]))+"\n")


def ref_trans_id_2_gene_id():
    """
    function to return dict linking transcript ID to gene ID
    :return:
    """
    ref_transcript_id_2_gene_id_d = {}
    for gene_rec in ref_all_gene_recs_l:
        for mRNA_rec in refDB.children(id=gene_rec.id, featuretype="mRNA"):
            trans_id = mRNA_rec.attributes["transcript_id"][0]
            ref_transcript_id_2_gene_id_d[trans_id] = gene_rec.id
    return ref_transcript_id_2_gene_id_d


def rreplace(input_s, exceptions_l):
    """
    function to replace the exceptional characters to "replaced" in a given string
    :param input_s: string, input string
    :param exceptions_l: list, exceptional characters
    :return: return_s: string, replaced string
    """
    return_s = ""
    for each_chr in input_s:
        if each_chr in exceptions_l:
            return_s += "replaced"
        else:
            return_s += each_chr
    return return_s


def get_new_name(gene_rec, mode_s):
    """
    function to name each gene record since some gene include "/" within gene symbol
    :param gene_rec: record, from annotation DB
    :param mode_s: string, ref for VGP and tgt for Old
    :return: string, newly named gene (e.g. gene-XXX=XM_XXX)
    """
    pass_symbol_l = ["/"]  # symbols to replace (since gene symbol would be used as a directory name)
    trans_id_s = "default"
    if mode_s == "ref":
        gene_symbol_s = gene_rec.id
        new_symbol_s = rreplace(input_s=gene_symbol_s, exceptions_l=pass_symbol_l)
        for mRNA_rec in refDB.children(id=gene_rec.id, featuretype="mRNA"):
            trans_id_s = mRNA_rec.attributes["transcript_id"][0]
        # one mRNA_rec for one gene_rec
        return new_symbol_s + "=" + trans_id_s
    else:
        assert mode_s == "tgt"
        for trans_rec in tgtDB.children(id=gene_rec.id, featuretype="transcript"):
            trans_id_s = trans_rec.attributes["source_transcript"][0]
        gene_symbol_s = ref_trans_id_2_gene_id_d[trans_id_s]
        new_symbol_s = rreplace(input_s=gene_symbol_s, exceptions_l=pass_symbol_l)
        # one trans_rec for one gene_rec
        return new_symbol_s + "=" + trans_id_s


def get_softmasked(work_dir, fasta_file_s, mode_s):
    """
    function to get the coordinates of windowmasked (lower cases) region from assemblies
    :param work_dir: string, path of species working directory
    :param fasta_file_s: string, path of input fasta file
    :param mode_s: string, ref for VGP, tgt for Old assembly
    :return:
    """
    dir_check(os.path.join(work_dir, "basic_inputs"))
    lowercase_bed_file_s = os.path.join(work_dir, "basic_inputs", "%s_softmasked.bed" % mode_s)
    with open(fasta_file_s) as fasta_file_f, open(lowercase_bed_file_s, "w") as lowercase_bed_file_f:
        for values in FastaIO.SimpleFastaParser(fasta_file_f):
            header_s = values[0]
            seq_s = values[1]
            lower_cases_rex = re.compile("[atcgn]+")
            coordinates_nl = list(x.span() for x in re.finditer(pattern=lower_cases_rex, string=seq_s))
            for coordinate_l in coordinates_nl:
                lowercase_bed_file_f.write("\t".join(map(str, [header_s, coordinate_l[0], coordinate_l[1]]))+"\n")


def compare_with_remap():
    """
    function to compare totally missing and exon deletion of our analysis with NCBI remap
    :return: nothing, just generates the summary file
    """
    working_dir_s = os.path.join(result_dir, "remap")
    # (1) get chromosome name (Super_scaffold) and its NCBI RefSeq accession (NC_XX or NW_XX, without ".1" suffix)
    chr_name_2_accession_d = dict()
    with open(ref_assembly_report_file_s, "r") as assembly_report_file_f:
        for line_s in assembly_report_file_f:
            if line_s.startswith("#"):
                continue
            line_l = line_s.rstrip("\n").split("\t")
            seqname_s, refseq_id_s = line_l[0], line_l[6]
            if "." in refseq_id_s:
                refseq_id_s = refseq_id_s.split(".")[0]
            chr_name_2_accession_d[seqname_s] = refseq_id_s
    # (2) reading remap excel file to get the "nohit" regions
    remap_excel_file_wb = openpyxl.load_workbook(remap_excel_file_s)
    excel_sheet_names_l = remap_excel_file_wb.sheetnames
    chr_names_l = [sheet_name_s for sheet_name_s in excel_sheet_names_l if "chr" in sheet_name_s]
    chr_names_l.remove("chr summary")
    remap_discrepancy_summary_d = {}  # {chr.1_1 = ["nohit", "NW_XXX", 13, 200, ...]}
    # placed scaffolds: chr 1, chr 2, ... , unplaced scaffolds: assigned together with chr Un
    for chr_name_s in chr_names_l:  # for each chromosomes, save all the discrepancy
        chr_sheet_ws = remap_excel_file_wb.get_sheet_by_name(chr_name_s)
        discrepancy_count_i = 0
        rows = chr_sheet_ws.rows
        next(rows) # to exclude header line
        for row in rows:
            discrepancy_count_i += 1
            discrepancy_id_s = chr_name_s.replace(" ",".") + "_" + str(discrepancy_count_i)
            discrepancy_type_s, seqname_s = row[0].value, row[2].value
            query_start_i, query_stop_i = int(row[4].value), int(row[5].value)
            target_name_s, target_seq_s = row[11].value, row[12].value
            try:
                target_start_i, target_end_i = int(row[13].value), int(row[14].value)
            except TypeError:
                target_start_i, target_end_i = None, None
            discrepancy_summary_l = [discrepancy_type_s, seqname_s, query_start_i, query_stop_i, target_name_s, target_seq_s, target_start_i, target_end_i]
            remap_discrepancy_summary_d[discrepancy_id_s] = discrepancy_summary_l

    # (3) subset only NoHit result (not aligned by remap)
    nohit_bedfile_s = os.path.join(working_dir_s, "nohit_summary.bed")
    nohit_bedlines_l = list()
    with open(nohit_bedfile_s, "w") as nohit_bedfile_f:
        discrepancy_summaries_nl = list(remap_discrepancy_summary_d.values())
        nohit_summaries_nl = list(filter(lambda x: x[0] == "NoHit", discrepancy_summaries_nl))
        for nohit_summary_l in nohit_summaries_nl:
            seqname_s, ref_start_i, ref_end_i = nohit_summary_l[1], nohit_summary_l[2], nohit_summary_l[3]
            nohit_bedline_s = "\t".join(map(str,[chr_name_2_accession_d[seqname_s], ref_start_i-1, ref_end_i]))
            nohit_bedlines_l.append(nohit_bedline_s)
        nohit_bedfile_f.write("\n".join(nohit_bedlines_l))


def cut_neighbor_sequences(seq_s, flanking_i):
    """
    cut the flanking sequences
    :param seq_s: string, seq
    :param flanking_i: size of flanking seq
    :return: strings, cut (start), cut (the rest), cut (last)
    """
    assert type(seq_s) is str
    return seq_s[0:flanking_i], seq_s[flanking_i:-flanking_i], seq_s[-flanking_i:]


def identity(seq1_s, seq2_s):
    """
    function to calculate identity between string; only considering substitution
    :param seq1_s: string, seq1
    :param seq2_s: string, seq2
    :return: integer, identity
    """
    ident_i = 0
    assert len(seq1_s) == len(seq2_s)
    seq1_s = seq1_s.upper()
    seq2_s = seq2_s.upper()
    for pair_t in zip(seq1_s, seq2_s):
        if pair_t[0] == pair_t[1]:
            ident_i += 1
    return ident_i


def convert_header_of_assmebly_statistics(assembly_report_file_s):
    """
    function to convert genbank id of each scaffold to refseq id
    :param assembly_report_file_s: string, path to assembly report file (basic input)
    :return: genbank_2_refseq_d: {genbank_id: refseq_id}
    """
    assert file_check(assembly_report_file_s) is True
    genbank_2_refseq_d = {}
    with open(assembly_report_file_s) as assembly_report_file_f:
        for line_s in assembly_report_file_f:
            line_l = line_s.rstrip("\n").split("\t")
            # header line or mitochondrial scaffold
            if line_s.startswith("#") or line_l[2] == "MT":
                continue
            genbank_id_s, refseq_id_s = line_l[4], line_l[6].split(".")[0]
            if genbank_id_s != "na" and refseq_id_s != "na":
                pass
            elif genbank_id_s == "na" and refseq_id_s != "na":  # if genbank id is not available --> use refseq id
                genbank_id_s = refseq_id_s
            elif genbank_id_s != "na" and refseq_id_s == "na":  # if refseq id is not available --> use genbank id
                refseq_id_s = genbank_id_s.split(".")[0]
            genbank_2_refseq_d[genbank_id_s] = refseq_id_s
    return genbank_2_refseq_d


def convert_genomic_gaps_to_bed(gap_info_file_s, gap_info_bedfile_s, genbank_2_refseq_d):
    """
    funcion to converte genomic gap file (from NCBI) to bed file
    :param gap_info_file_s: string, path to genomic gap file from NCBI
    :param gap_info_bedfile_s: string, path to output bed file
    :param genbank_2_refseq_d: dict, from convert_header_of_assmebly_statistics
    :return:
    """
    assert file_check(gap_info_file_s) is True
    with open(gap_info_file_s) as gap_info_file_f, open(gap_info_bedfile_s, "w") as gap_info_bedfile_f:
        count_i = 0
        gap_info_file_f.readline()
        for line_s in gap_info_file_f:
            count_i += 1
            line_l = line_s.rstrip("\n").split("\t")
            genbank_id_s, start_s, end_s, gap_length_i, gap_type_s, linkage_evidence_s = line_l
            id_s = "gap_%d_%s_%s" % (count_i, gap_type_s, linkage_evidence_s)
            refseq_id_s = genbank_2_refseq_d.get(genbank_id_s, genbank_id_s.split(".")[0])
            gap_info_bedfile_f.write("\t".join(map(str, [refseq_id_s, int(start_s) - 1, int(end_s), id_s, 0, "+"])) + "\n")


def gc_and_repeat_content_of_genes(seq_and_gene_name_l):
    """
    function to return GC and repeat content of gene
    :param header_n_seq_l:
    :return:
    """
    rawseq_s, gene_name_s = seq_and_gene_name_l
    gene_2_info_d = {}
    gene_2_info_d["gene"] = gene_name_s
    # GC content
    upper_seq_s = rawseq_s.upper()
    upper_wo_n_i = len(upper_seq_s.replace("N", ""))
    gc_count_i = upper_seq_s.count("G") + upper_seq_s.count("C")
    try:
        gc_content_fl = float(gc_count_i / upper_wo_n_i * 100)
    except ZeroDivisionError:
        gc_content_fl = 0
    gene_2_info_d["gc"] = gc_content_fl
    # repeat content
    upper_i = sum(1 for c in rawseq_s if c.isupper())
    lower_i = sum(1 for c in rawseq_s if c.islower())
    assert lower_i + upper_i == len(rawseq_s)
    gene_2_info_d["repeat"] = float(lower_i / len(rawseq_s) * 100)
    return gene_2_info_d


def check_scaffold_type(scaffold_id_s):
    """
    classify scaffold type (chromosome, unlocalized scaffold, unplaced scaffold)
    :param scaffold_id_s: string, scaffold id
    :return: string, type of scaffold
    """
    if scaffold_id_s in scaff_2_chrom_d.keys():
        if scaffold_id_s in unlocalized_d.keys():
            scaffold_type_s = "unlocalized_scaffold"
        else:
            scaffold_type_s = "chromosome"
    else:
        scaffold_type_s = "unplaced_scaffold"
    return scaffold_type_s


def classify_variant_id_error(ref_bool, tgt_bool, fgl_s, variant_id_error_type_s, variant_id_s):
    """
    function to classify (no_errors, previous_error (=False gene loss), vgp_error, filter_out)
    :param ref_bool:
    :param tgt_bool: boolean, True: no assembly error in tgt (prior) assembly, False: assembly error in ref (VGP) assembly,
    :param fgl_s: string, type of false gene loss: one of frameshift, prematurestopcodon, or intronexonjunctiondisruption
    :param variant_id_error_type_s:
    :param variant_id_s:
    :return:
    """
    class_s = ""
    fgl_ids_l, vgp_error_ids_l = [], []
    if fgl_s == "frameshift":
        if ref_bool is True and tgt_bool is True:
            class_s = "no_errors"
        elif ref_bool is True and tgt_bool is False:
            if "INS" in variant_id_error_type_s or "DEL" in variant_id_error_type_s:
                class_s = "previous_error"
                fgl_ids_l.append(variant_id_s)
            else:
                class_s = "filter_out"
        elif ref_bool is False and tgt_bool is True:
            if "INS" in variant_id_error_type_s or "DEL" in variant_id_error_type_s:
                class_s = "vgp_error"
                vgp_error_ids_l.append(variant_id_s)
            else:
                class_s = "filter_out"
        elif ref_bool is False and tgt_bool is False:
            class_s = "filter_out"
        else:
            assert "NA" in [ref_bool, tgt_bool]
            class_s = "filter_out"

    else:
        assert fgl_s in ["prematurestopcodon", "intronexonjunctiondisruption"]
        if ref_bool is True and tgt_bool is True:
            class_s = "no_errors"
        elif ref_bool is True and tgt_bool is False:
            class_s = "previous_error"
            fgl_ids_l.append(variant_id_s)
        elif ref_bool is False and tgt_bool is True:
            class_s = "vgp_error"
            vgp_error_ids_l.append(variant_id_s)
        elif ref_bool is False and tgt_bool is False:
            class_s = "filter_out"
    return class_s


class FalseGeneLoss():
    #######################
    # Basic parsing codes #
    #######################
    def gff2cds(self, info_l):
        """
        function to
        :param info_l: list, [gene_rec: gene record from gffutils DB, mode_s: string, "ref(VGP)" or "tgt(Prior)"]
        :return:
        """
        # (0) basic variables
        gene_rec, mode_s = info_l
        if mode_s == "ref":
            fasta_s = ref_fasta_file_s
            DB = refDB
        else:
            assert mode_s == "tgt"
            fasta_s = tgt_fasta_file_s
            DB = tgtDB
        gene_id_s = gene_rec.id
        gene_name_s = get_new_name(gene_rec=gene_rec, mode_s=mode_s)
        gene_strand = gene_rec.strand
        gene_chrom_s = gene_rec.chrom
        working_dir_s = os.path.join(result_dir, "CDS", mode_s, gene_name_s)
        dir_check(working_dir_s)

        # (1) CDS records: for detecting sequence-level false gene loss in prior assembly
        cds_recs_l = sorted([cds_rec for cds_rec in DB.children(id=gene_id_s, featuretype="CDS")], key=lambda cds_rec: cds_rec.start)
        # (1-1) cds gff file: cds record of gff file, automatically correct errors (--> for reference: VGP)
        cds_gfffile_s = os.path.join(working_dir_s, "%s.gff3" % gene_name_s)
        with open(cds_gfffile_s, "w") as cds_gfffile_f:
            cds_gff_lines_l = []
            for cds_rec in cds_recs_l:
                cds_gff_lines_l.append(str(cds_rec))
            cds_gfffile_f.write("\n".join(cds_gff_lines_l))
        # (1-2) cds bed file: cds record of bed file
        #                     additional check to include 1 or 2bp introns in CDS which was due to frameshift error in assembly (--> for target: Prior)
        cds_bedfile_s = os.path.join(working_dir_s, "cds_%s.bed" % gene_name_s)
        new_starts_l, new_ends_l = ignore_1_or_2bp_introns(raw_rec_l=cds_recs_l)
        original_cds_count_i = len(new_starts_l)
        cds_order_i = 0
        revised_cds_bed_lines_l = []
        for (start_i, end_i) in zip(new_starts_l, new_ends_l):
            cds_order_i += 1
            cds_start_i, cds_end_i = start_i - 1, end_i
            cds_position_type_s = self.check_position(order_i=cds_order_i, count_i=original_cds_count_i, strand_c=gene_strand, mode_s="CDS")
            cds_level_id_s = "|".join(map(str, [gene_name_s, cds_order_i, original_cds_count_i, cds_position_type_s]))
            cds_bed_line_s = "\t".join(map(str, [gene_chrom_s, cds_start_i, cds_end_i, cds_level_id_s, 0, gene_strand]))
            revised_cds_bed_lines_l.append(cds_bed_line_s)
        with open(cds_bedfile_s, "w") as cds_bedfile_f:
            cds_bedfile_f.write("\n".join(revised_cds_bed_lines_l)+"\n")
        # (1-3) convert bed/gff file to fasta
        cds_fastafile_s = os.path.join(working_dir_s, "%s.fasta" % gene_name_s)
        merged_cds_fastafile_s = os.path.join(working_dir_s, "merged_%s.fasta" % gene_name_s)
        if mode_s == "ref":
            os.system("%s getfasta -fi %s -bed %s -s -name+ -fullHeader > %s" % (bedtools_command_s, fasta_s, cds_gfffile_s, cds_fastafile_s))
        else:
            assert mode_s == "tgt"
            os.system("%s getfasta -fi %s -bed %s -s -name+ -fullHeader > %s" % (bedtools_command_s, fasta_s, cds_bedfile_s, cds_fastafile_s))
        # (1-4) generate merged fasta file
        with open(cds_fastafile_s, "r") as cds_fastafile_f, open(merged_cds_fastafile_s, "w") as merged_cds_fastafile_f:
            header_2_seq_nl = []
            for values in FastaIO.SimpleFastaParser(cds_fastafile_f):
                header_s, seq_s = values[0], values[1]
                header_2_seq_nl.append([header_s, seq_s])
            # sort by CDS start coordinate: int(header_2_seq_l[0].split(":")[-1].split("-")[0])
            if gene_strand == "+":
                sorted_header_2_seq_nl = sorted(header_2_seq_nl, key=lambda header_2_seq_l: int(header_2_seq_l[0].split(":")[-1].split("-")[0]))
            else:
                assert gene_strand == "-"
                sorted_header_2_seq_nl = sorted(header_2_seq_nl, key=lambda header_2_seq_l: int(header_2_seq_l[0].split(":")[-1].split("-")[0]), reverse=True)
            sorted_seqs_l = [header_2_seq_l[1] for header_2_seq_l in sorted_header_2_seq_nl]
            merged_cds_seq_s = "".join(sorted_seqs_l)
            merged_cds_fastafile_f.write(">%s_CDS\n%s" % (gene_name_s, merged_cds_seq_s))

        # (2) Exon records: for blast
        exon_recs_l = sorted([exon_rec for exon_rec in DB.children(id=gene_id_s, featuretype="exon")], key=lambda exon_rec: exon_rec.start)
        exon_bedfile_s = os.path.join(working_dir_s, "exon_%s.bed" % gene_name_s)
        # (2-1) exon bed file: exon record of bed file
        #                     additional check to include 1 or 2bp introns in exon which was due to frameshift error in assembly (--> for both assemblies)
        exon_bed_lines_l = []
        new_starts_l, new_ends_l = ignore_1_or_2bp_introns(raw_rec_l=exon_recs_l)
        original_exon_count_i = len(new_starts_l)
        exon_order_i = 0
        for (start_i, end_i) in zip(new_starts_l, new_ends_l):
            exon_order_i += 1
            exon_start_i, exon_end_i = start_i - 1, end_i
            exon_position_type_s = self.check_position(order_i=exon_order_i, count_i=original_exon_count_i, strand_c=gene_strand, mode_s="Exon")
            exon_level_id_s = "|".join(map(str, [gene_name_s, exon_order_i, original_exon_count_i, exon_position_type_s]))
            exon_bed_line_s = "\t".join(map(str, [gene_chrom_s, exon_start_i, exon_end_i, exon_level_id_s, 0, gene_strand]))
            exon_bed_lines_l.append(exon_bed_line_s)
        with open(exon_bedfile_s, "w") as exon_bedfile_f:
            exon_bedfile_f.write("\n".join(exon_bed_lines_l)+"\n")
        # (2-2) generate exon fasta file whose records were named after the order of exons
        exon_fastafile_s = os.path.join(working_dir_s, "exon_%s.fasta" % gene_name_s)
        exon_by_order_fasta_file_s = os.path.join(working_dir_s, "EXON_%s.fasta" % gene_name_s)
        os.system("%s getfasta -fi %s -bed %s -s -fullHeader > %s" % (bedtools_command_s, fasta_s, exon_bedfile_s, exon_fastafile_s))
        with open(exon_fastafile_s) as exon_fastafile_f, open(exon_by_order_fasta_file_s,'w') as rename_exon_fastafile_f:
            header_2_seq_d = {}
            header_2_count_d = {}
            count_i = 0
            for values in FastaIO.SimpleFastaParser(exon_fastafile_f):
                count_i += 1
                header_s = values[0]
                header_2_seq_d[header_s] = values[1] # exon sequence
                header_2_count_d[header_s] = count_i
            for key_s in header_2_seq_d.keys():
                exon_seq_s = header_2_seq_d[key_s]
                if gene_strand == "+":
                    exon_id_by_order_s = "exon" + str(header_2_count_d[key_s])
                else:
                    assert gene_strand == "-"
                    exon_id_by_order_s = "exon" + str(count_i+1-header_2_count_d[key_s])
                rename_exon_fastafile_f.write(">%s_%s\n%s\n" % (exon_id_by_order_s, gene_name_s, exon_seq_s))
        return ">%s_CDS\n%s" % (gene_name_s, merged_cds_seq_s)

    def mp_gff2cds(self, num_processors_i):
        """
        Multiprocess version of GFF2CDS
        :param num_processors_i: number of multiprocess cores
        :return: list, out: list of strings, but actually wouldn't use it
        """
        dir_check(os.path.join(result_dir, "CDS"))
        modes_l = ['tgt', 'ref']
        for mode_s in modes_l:
            dir_check(os.path.join(result_dir, "CDS", mode_s))
            consensus_cds_file_s = os.path.join(result_dir, "CDS", mode_s, "%s_consensus.cds.fasta" % mode_s)
            if mode_s == "tgt":
                gene_recs = tgt_all_gene_recs_l
            else:
                assert mode_s == "ref"
                gene_recs = ref_all_gene_recs_l
            infos_l = [[gene_rec, mode_s] for gene_rec in gene_recs]
            with Pool(num_processors_i) as p:
                out = p.map(self.gff2cds, infos_l)
            p.close()
            p.join()
            with open(consensus_cds_file_s, "w") as consensus_cds_file_f:
                consensus_cds_file_f.write("\n".join(out))

    ##############################
    # 1. Missing region analysis #
    ##############################
    def missing_chrom(self, num_processors_i):
        """
        function to perform missing region analysis
        :param num_processors_i: int, number of processors to use
        :return: X
        """
        # (1) make genome bed files from faidx files
        working_dir_s = os.path.join(result_dir, "missing_chrom")
        ref_faidx_file_s = ref_fasta_file_s + ".fai"
        tgt_faidx_file_s = tgt_fasta_file_s + ".fai"
        ref_genome_bed_file_s = os.path.join(result_dir, "basic_inputs", "ref.genome.bed")
        tgt_genome_bed_file_s = os.path.join(result_dir, "basic_inputs", "tgt.genome.bed")
        if (file_check(ref_faidx_file_s) is False) or (file_check(tgt_faidx_file_s) is False):
            os.system("%s faidx %s" % (samtools_command_s, ref_fasta_file_s))
            os.system("%s faidx %s" % (samtools_command_s, tgt_fasta_file_s))
        os.system("awk -F '\\t' -v OFS='\\t' '{print $1, 0, $2, $1, 0, \".\"}' %s > %s" % (ref_faidx_file_s, ref_genome_bed_file_s))
        os.system("awk -F '\\t' -v OFS='\\t' '{print $1, 0, $2, $1, 0, \".\"}' %s > %s" % (tgt_faidx_file_s, tgt_genome_bed_file_s))

        # (2) convert genomic gap coordinates into bed format
        ref_gap_info_bedfile_s = os.path.join(result_dir, "basic_inputs", "ref.gaps.bed")
        tgt_gap_info_bedfile_s = os.path.join(result_dir, "basic_inputs", "tgt.gaps.bed")
        ref_genbank_2_refseq_d = convert_header_of_assmebly_statistics(assembly_report_file_s=ref_assembly_report_file_s)
        tgt_genbank_2_refseq_d = convert_header_of_assmebly_statistics(assembly_report_file_s=tgt_assembly_report_file_s)
        convert_genomic_gaps_to_bed(gap_info_file_s=ref_gap_info_file_s, gap_info_bedfile_s=ref_gap_info_bedfile_s, genbank_2_refseq_d=ref_genbank_2_refseq_d)
        convert_genomic_gaps_to_bed(gap_info_file_s=tgt_gap_info_file_s, gap_info_bedfile_s=tgt_gap_info_bedfile_s, genbank_2_refseq_d=tgt_genbank_2_refseq_d)

        # (3) Cactus summary: halliftover and get the not aligned region not overlapping with N
        dir_check(os.path.join(working_dir_s, "1_halliftover"))
        ref_cactus_aligned_bed_file_s = os.path.join(working_dir_s, "1_halliftover", "ref.cactus.aligned.bed")
        ref_cactus_not_aligned_bed_file_s = os.path.join(working_dir_s, "1_halliftover", "ref.cactus.not_aligned.bed")
        if file_check(ref_cactus_not_aligned_bed_file_s) is False:
            os.system("%s subtract -a %s -b %s > %s" % (bedtools_command_s, tgt_genome_bed_file_s, tgt_gap_info_bedfile_s,
                                                        tgt_genome_bed_file_s + ".remove_gap.bed"))
            os.system("%s --inMemory %s %s %s %s %s" %
                      (halLiftover_commnd_s, hal_path_s, tgt_s, tgt_genome_bed_file_s + ".remove_gap.bed",
                       ref_s, ref_cactus_aligned_bed_file_s+".not_merged.bed"))
            os.system("%s sort -i %s | %s merge -i stdin | %s subtract -a stdin -b %s > %s" %
                      (bedtools_command_s, ref_cactus_aligned_bed_file_s + ".not_merged.bed",
                       bedtools_command_s,
                       bedtools_command_s, ref_gap_info_bedfile_s, ref_cactus_aligned_bed_file_s))
            os.system("%s subtract -a %s -b %s | %s subtract -a stdin -b %s > %s " %
                      (bedtools_command_s, ref_genome_bed_file_s, ref_cactus_aligned_bed_file_s,
                       bedtools_command_s, ref_gap_info_bedfile_s, ref_cactus_not_aligned_bed_file_s))

        # (4) Minimap2 summary:
        dir_check(os.path.join(working_dir_s, "2_minimap2"))
        ref_minimap2_aligned_paf_file_s = os.path.join(working_dir_s, "2_minimap2", "ref.minimap2.aligned.paf")
        ref_minimap2_aligned_bed_file_s = os.path.join(working_dir_s, "2_minimap2", "ref.minimap2.aligned.bed")
        ref_minimap2_aligned_bed3_file_s = os.path.join(working_dir_s, "2_minimap2", "ref.minimap2.aligned.three_col.bed")
        ref_minimap2_not_aligned_bed_file_s = os.path.join(working_dir_s, "2_minimap2", "ref.minimap2.not_aligned.bed")
        if file_check(ref_minimap2_not_aligned_bed_file_s) is False:
            os.system("%s -x asm5 -r 50 -c -g 1000 -t %d --no-long-join %s %s > %s" %
                      (minimap2_command_s, num_processors_i, ref_fasta_file_s, tgt_fasta_file_s, ref_minimap2_aligned_paf_file_s))
            os.system("%s %s > %s ; cut -f 1,2,3 %s | %s sort -i stdin | %s merge -i stdin > %s" %
                      (paftools_command_s, ref_minimap2_aligned_paf_file_s, ref_minimap2_aligned_bed_file_s,
                       ref_minimap2_aligned_bed_file_s, bedtools_command_s, bedtools_command_s, ref_minimap2_aligned_bed3_file_s))
            os.system("%s subtract -a %s -b %s | %s subtract -a stdin -b %s > %s" %
                      (bedtools_command_s, ref_genome_bed_file_s, ref_minimap2_aligned_bed3_file_s,
                       bedtools_command_s, ref_gap_info_bedfile_s, ref_minimap2_not_aligned_bed_file_s))

        # (5) Overlapping missing regions: supported from both aligner
        dir_check(os.path.join(working_dir_s, "3_consensus"))
        consensus_missing_bed_file_s = os.path.join(working_dir_s, "3_consensus", "consensus.missing.bed")
        consensus_aln_bed_file_s = os.path.join(working_dir_s, "3_consensus", "consensus.aligned.bed")
        no_align_scaffold_level_s = os.path.join(working_dir_s, "3_consensus", "all_scaffolds.missing_ratio.tsv")
        no_align_chromosome_level_s = os.path.join(working_dir_s, "3_consensus", "chromosomes.missing_ratio.tsv")
        missing_summary_file_s = os.path.join(working_dir_s, "3_consensus", "genome.summary.tsv")

        # (5-1) get the coordinates of missing & aligned regions
        os.system("%s intersect -a %s -b %s | %s subtract -a stdin -b %s > %s" %
                  (bedtools_command_s, ref_cactus_not_aligned_bed_file_s, ref_minimap2_not_aligned_bed_file_s,
                   bedtools_command_s, ref_gap_info_bedfile_s, consensus_missing_bed_file_s))
        os.system("%s subtract -a %s -b %s | %s subtract -a stdin -b %s > %s" %
                  (bedtools_command_s, ref_genome_bed_file_s, ref_gap_info_bedfile_s,
                   bedtools_command_s, consensus_missing_bed_file_s, consensus_aln_bed_file_s))
        ref_scaffold_2_type_d = {scaffold_id_s: {"missing":0, "aligned":0, "gap":0} for scaffold_id_s in ref_chrom_2_size_d}
        with open(consensus_missing_bed_file_s,"r") as consensus_missing_bed_file_f, \
                open(consensus_aln_bed_file_s, "r") as consensus_aln_bed_file_f, \
                open(ref_gap_info_bedfile_s, "r") as ref_gap_info_bedfile_f, \
                open(no_align_scaffold_level_s, "w") as no_align_scaffold_level_f,\
                open(no_align_chromosome_level_s, "w") as no_align_chromosome_level_f,\
                open(missing_summary_file_s, "w") as missing_summary_file_f:
            for line_s in consensus_missing_bed_file_f:
                line_l = line_s.rstrip("\n").split("\t")
                scaffold_id_s, start_i, end_i = line_l[0], int(line_l[1]), int(line_l[2])
                not_aligned_length_i = end_i - start_i
                ref_scaffold_2_type_d[scaffold_id_s]["missing"] += not_aligned_length_i
            for line_s in consensus_aln_bed_file_f:
                line_l = line_s.rstrip("\n").split("\t")
                scaffold_id_s, start_i, end_i = line_l[0], int(line_l[1]), int(line_l[2])
                aligned_length_i = end_i - start_i
                ref_scaffold_2_type_d[scaffold_id_s]["aligned"] += aligned_length_i
            for line_s in ref_gap_info_bedfile_f:
                line_l = line_s.rstrip("\n").split("\t")
                scaffold_id_s, start_i, end_i = line_l[0], int(line_l[1]), int(line_l[2])
                gap_length_i = end_i - start_i
                ref_scaffold_2_type_d[scaffold_id_s]["gap"] += gap_length_i
            header_l = ["chrom", "size", "missing_length", "aligned_length", "gap_length", "missing_ratio", "aligned_ratio", "gap_ratio"]
            no_align_scaffold_level_f.write("\t".join(header_l) + "\n")
            no_align_chromosome_level_f.write("\t".join(header_l) + "\n")
            scaffold_ids_l = sorted(list(ref_chrom_2_size_d.keys()))
            genome_missing_length_i = 0
            genome_aligned_length_i = 0
            genome_gap_length_i = 0
            genome_size_i = 0
            for scaffold_id_s in scaffold_ids_l:
                size_i = ref_chrom_2_size_d[scaffold_id_s]
                missing_length_i = ref_scaffold_2_type_d[scaffold_id_s]["missing"]
                aligned_length_i = ref_scaffold_2_type_d[scaffold_id_s]["aligned"]
                gap_length_i = ref_scaffold_2_type_d[scaffold_id_s]["gap"]
                genome_size_i += size_i
                genome_missing_length_i += missing_length_i
                genome_aligned_length_i += aligned_length_i
                genome_gap_length_i += gap_length_i
                missing_ratio_fl = round(missing_length_i/size_i*100, 3)
                aligned_ratio_fl = round(aligned_length_i/size_i*100, 3)
                gap_ratio_fl = round(gap_length_i/size_i*100, 3)
                no_align_scaffold_level_f.write("\t".join(map(str,[scaff_2_chrom_d.get(scaffold_id_s, scaffold_id_s),
                                                                   size_i, missing_length_i, aligned_length_i, gap_length_i,
                                                                   missing_ratio_fl, aligned_ratio_fl, gap_ratio_fl]))+'\n')
            missing_summary_file_f.write("\t".join(["type", "length", "ratio"]) + "\n")
            missing_summary_file_f.write("\t".join(["missing", str(genome_missing_length_i),
                                                    "%f" % round(genome_missing_length_i / genome_size_i*100, 1) + "%"]) + "\n")
            missing_summary_file_f.write("\t".join(["aligned", str(genome_aligned_length_i),
                                                    "%f" % round(genome_aligned_length_i / genome_size_i*100, 1) + "%"]) + "\n")
            missing_summary_file_f.write("\t".join(["gap", str(genome_gap_length_i),
                                                    "%f" % round(genome_gap_length_i / genome_size_i*100, 1) + "%"]) + "\n")

    ############################
    # 2. Missing gene analysis #
    ############################
    def main_missing_genes_analysis(self, num_processors_i):
        """
        function to calculate missing ratio of genes and exons
        :param num_processors_i: num of multiprocessing cores
        :return: none, generating summaries
        """
        working_dir_s = os.path.join(result_dir, "missingGenes")
        basic_input_dir_s = os.path.join(result_dir, "basic_inputs")
        dir_check(working_dir_s)
        dir_check(os.path.join(working_dir_s, "BLAST"))
        dir_check(basic_input_dir_s)
        ref_all_cdsfile_s = os.path.join(result_dir, "CDS", "ref", "ref_consensus.cds.fasta")
        tgt_all_cdsfile_s = os.path.join(result_dir, "CDS", "tgt", "tgt_consensus.cds.fasta")
        if file_check(ref_all_cdsfile_s) is False or file_check(tgt_all_cdsfile_s) is False:
            self.mp_gff2cds(num_processors_i=process_num_i)
        # (1) Based on genome-wide alignment, calculate unaligned ratio of genes, intergenic regions, exons, cds, and introns
        # (1-1) input1: coordinates of each type of genetic features
        ref_gene_bed_file_s = os.path.join(basic_input_dir_s, "ref.genes_before_flt.bed")
        ref_intergenic_bed_file_s = os.path.join(basic_input_dir_s, "ref.intergenic.bed")
        ref_exon_bed_file_s = os.path.join(basic_input_dir_s, "ref.exon.bed")
        ref_cds_bed_file_s = os.path.join(basic_input_dir_s, "ref.cds.bed")
        ref_intron_bed_file_s = os.path.join(basic_input_dir_s, "ref.intron.bed")
        ref_non_coding_bed_file_s = os.path.join(basic_input_dir_s, "ref.noncoding_genes.bed")
        if True:# file_check(ref_gene_bed_file_s) is False:
            self.gene_2_bed(mode_s="ref")
        if True:# file_check(ref_intergenic_bed_file_s) is False:
            self.intergenic_2_bed(mode_s="ref")
        if True:# file_check(ref_exon_bed_file_s) is False:
            self.exon_2_bed(mode_s="ref")
        if True:# file_check(ref_cds_bed_file_s) is False:
            self.cds_2_bed(mode_s="ref")
        if True:# file_check(ref_intron_bed_file_s) is False:
            self.intron_2_bed(mode_s="ref")
        # (1-2) input2: consensus missing coordinates
        missing_bed_file_s = os.path.join(result_dir, "missing_chrom", "3_consensus", "consensus.missing.bed")
        if file_check(missing_bed_file_s) is False:
            self.missing_chrom(num_processors_i=process_num_i)
        # (1-3) calculate the intersection
        compare_missing_file_s = os.path.join(working_dir_s, "compare.missing.tsv")
        with open(compare_missing_file_s,"w") as compare_missing_file_f:
            compare_missing_file_f.write("\t".join(["type", "missing", "all", "ratio"])+"\n")
            for position_s in ["genes_before_flt", "intergenic", "exon", "intron", "cds", "noncoding_genes"]:
                all_rec_bed_file_s = os.path.join(result_dir, "basic_inputs", "ref.%s.bed" % position_s)
                missing_rec_bed_file_s = os.path.join(working_dir_s, "ref.%s.not_aligned.bed" % position_s)
                os.system("%s intersect -a %s -b %s -wao | " % (bedtools_command_s, all_rec_bed_file_s, missing_bed_file_s) +
                          "%s groupby -i stdin -g 1,2,3,4 -c 13 -o sum |" % bedtools_command_s +
                          "awk -F'[\"\\t\"]' '{OFS=\"\\t\"} {print $1, $2, $3, $4, $5, $5/($3-$2)*100}' > %s" % missing_rec_bed_file_s)
                all_cmd = "%s sort -i %s | %s merge -i stdin | awk -F'[\"\\t\"]' '{OFS=\"\\t\"} {sum += $3-$2} END {print sum}'" % (bedtools_command_s,
                                                                                                                                    all_rec_bed_file_s,
                                                                                                                                    bedtools_command_s)
                all_ps = subprocess.Popen(all_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                all_output = int(all_ps.communicate()[0])
                missing_cmd = "%s sort -i %s| %s merge -i stdin |" % (bedtools_command_s, all_rec_bed_file_s, bedtools_command_s)+\
                              "%s intersect -a stdin -b %s |" % (bedtools_command_s, missing_bed_file_s) + \
                              "awk -F'[\"\\t\"]' '{OFS=\"\\t\"} {sum += $3-$2} END {print sum}'"
                missing_ps = subprocess.Popen(missing_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                missing_output = int(missing_ps.communicate()[0])
                compare_missing_file_f.write("\t".join(map(str,
                                                           [position_s, missing_output, all_output, missing_output/all_output*100])) +"\n")
        # (2) Perform exon blast & summary the result with the genome-wide alignment evidences
        ref_gene_missing_file_s = os.path.join(working_dir_s, "ref.genes.not_aligned.bed")
        ref_exon_missing_file_s = os.path.join(working_dir_s, "ref.exon.not_aligned.bed")
        exon_lost_bed_file_s = os.path.join(working_dir_s, "lost_exons.bed")
        gene_2_lost_exon_pickle_s = os.path.join(working_dir_s, "gene_2_lost_exon_d.pickle")
        missing_summary_file_s = os.path.join(working_dir_s, "totally_missing_and_exon_deletion.tabs")
        # (2-1) Exon Blast
        if file_check(missing_summary_file_s) is False:
            # Check if blast db exists
            dir_check(os.path.join(result_dir,"blastDB"))
            ref_blast_db_s = os.path.join(result_dir, "blastDB", ref_s)
            if file_check(ref_blast_db_s + ".nhr") is False:
                make_blast_db(fasta_s=ref_fasta_file_s, database_s=ref_blast_db_s)
            tgt_blast_db_s = os.path.join(result_dir, "blastDB", tgt_s)
            if file_check(tgt_blast_db_s + ".nhr") is False:
                make_blast_db(fasta_s=tgt_fasta_file_s, database_s=tgt_blast_db_s)
            # Perform exon-wide blast
            infos_l = [[gene_rec, tgt_blast_db_s] for gene_rec in ref_all_gene_recs_l]
            with Pool(num_processors_i) as p:
                out = p.map(exon_blast, infos_l)
            p.close()
            p.join()
        # (2-2) Summary blast & genome-wide alignment results
        try:
            with open(gene_2_lost_exon_pickle_s,"rb") as gene_2_lost_exon_pickle_f:
                gene_2_lost_exon_d = pickle.load(gene_2_lost_exon_pickle_f)
        except FileNotFoundError:
            with open(ref_gene_missing_file_s) as ref_gene_missing_file_f,\
                    open(ref_exon_missing_file_s) as ref_exon_missing_file_f,\
                    open(missing_summary_file_s,"w") as missing_summary_file_f,\
                    open(exon_lost_bed_file_s,"w") as exon_lost_bed_file_f,\
                    open(gene_2_lost_exon_pickle_s, "wb") as gene_2_lost_exon_pickle_f:
                missing_summary_file_f.write("\t".join(map(str, ["gene", "num_of_CDS", "cactus_aligned", "cactus_not_aligned",
                                                                 "blast_hits", "blast_no_hits", "type"])) + "\n")
                # a.Genome-wide alignment summary
                gene_id_2_missing_d = {}
                # a-1.gene level
                for line_s in ref_gene_missing_file_f:
                    line_l = line_s.rstrip("\n").split("\t")
                    gene_name_s, missing_ratio_gene_fl = line_l[3], float(line_l[5])
                    if gene_name_s not in gene_name_2_gene_id_d["ref"]:
                        continue
                    assert gene_name_s not in gene_id_2_missing_d
                    gene_id_2_missing_d[gene_name_s] = {"gene_aln": missing_ratio_gene_fl, "exon": {}}
                # a-2.exon level
                for line_s in ref_exon_missing_file_f:
                    line_l = line_s.rstrip("\n").split("\t")
                    gene_name_s, missing_ratio_exon_fl = line_l[3].split("|")[0], float(line_l[5])
                    if gene_name_s not in gene_name_2_gene_id_d["ref"]:
                        continue
                    gene_id_s = gene_name_2_gene_id_d["ref"][gene_name_s]
                    gene_strand_s = refDB[gene_id_s].strand
                    exon_order_i, exon_count_i = int(line_l[3].split("|")[1]), int(line_l[3].split("|")[2])
                    exon_id_s = "EXON_%d" % exon_order_i if gene_strand_s == "+" else "EXON_%d" % (exon_count_i + 1 - exon_order_i)
                    assert gene_name_s in gene_id_2_missing_d
                    assert exon_id_s not in gene_id_2_missing_d[gene_name_s]["exon"]
                    gene_id_2_missing_d[gene_name_s]["exon"][exon_id_s] = {}
                    gene_id_2_missing_d[gene_name_s]["exon"][exon_id_s]["missing"] = missing_ratio_exon_fl

                # b.blast summary
                gene_2_lost_exon_d = {}
                for gene_rec in ref_all_gene_recs_l:
                    gene_name_s = gene_id_2_gene_name_d["ref"][gene_rec.id]
                    ref_exon_bedfile_s = os.path.join(result_dir,"CDS","ref", gene_name_s, "exon_%s.bed" % gene_name_s)
                    exon_id_2_bed_d = {}
                    sorted_exon_bedlines_nl = []
                    # b-1. sort exons by their coordinates
                    with open(ref_exon_bedfile_s) as ref_exon_bedfile_f:
                        exon_bedlines_nl = [line_s.rstrip("\n").split("\t") for line_s in ref_exon_bedfile_f.readlines() if len(line_s) > 1]
                        for bed_line_l in exon_bedlines_nl:
                            bed_line_l[1], bed_line_l[2] = int(bed_line_l[1]), int(bed_line_l[2])
                            sorted_exon_bedlines_nl.append(bed_line_l)
                        if gene_rec.strand == "+":
                            sorted_exon_bedlines_nl = sorted(sorted_exon_bedlines_nl, key=lambda bed_l: bed_l[1])
                        else:
                            assert gene_rec.strand == "-"
                            sorted_exon_bedlines_nl = sorted(sorted_exon_bedlines_nl, key=lambda bed_l: bed_l[1], reverse=True)
                    # b-2. name the exons
                    num_of_exons_i = len(sorted_exon_bedlines_nl)
                    for exon_order_i in range(0,num_of_exons_i):
                        exon_id_2_bed_d["EXON_%d" % (exon_order_i + 1)] = sorted_exon_bedlines_nl[exon_order_i]

                    exon_ids_l = list(exon_id_2_bed_d.keys())
                    # b-3. sort exons by their coordinates
                    exon_file_s = os.path.join(result_dir, "missingGenes", "BLAST", gene_name_s, "%s.fasta" % gene_name_s)
                    blast_result_file_s = os.path.join(result_dir, "missingGenes", "BLAST", gene_name_s, "blast_%s.tabs" % gene_name_s)
                    with open(blast_result_file_s) as blast_result_file_f, open(exon_file_s) as exon_file_f:
                        exon_id_2_seq_d = {}
                        for values in FastaIO.SimpleFastaParser(exon_file_f):
                            header_s, seq_s = values[0], values[1]
                            exon_id_2_seq_d[header_s] = seq_s
                        # *** exon under 15bp: excluded from our analysis
                        for exon_id_s in exon_ids_l:
                            exon_exception_s = "default" if len(exon_id_2_seq_d[exon_id_s]) >= 15 else "short_than_15bp"
                            gene_id_2_missing_d[gene_name_s]["exon"][exon_id_s]["blast_hit"] = exon_exception_s
                    # b-4. summary blast result file: only save the hits with more than 90% of qcovus
                        blast_lines_l = list(filter(lambda x: x.startswith("#") is False, blast_result_file_f.readlines()))
                        if len(blast_lines_l) > 0:  # == if there is any hit
                            for blast_line_s in blast_lines_l:
                                blast_line_l = blast_line_s.rstrip("\n").split("\t")
                                query_s = blast_line_l[0]
                                identity_fl, qcoverage_fl = float(blast_line_l[2]), float(blast_line_l[-1])
                                if qcoverage_fl >= 90:
                                    assert query_s in gene_id_2_missing_d[gene_name_s]["exon"]
                                    gene_id_2_missing_d[gene_name_s]["exon"][query_s]["blast_hit"] = "hit"
                    # b-5. compare genome-wide alignment vs. blast result
                    gene_missing_fl = gene_id_2_missing_d[gene_name_s]["gene_aln"]
                    # missing: from genome-wide alignment
                    aln_exons_l = [exon_id_s for exon_id_s in gene_id_2_missing_d[gene_name_s]["exon"]
                                   if gene_id_2_missing_d[gene_name_s]["exon"][exon_id_s]["missing"] < 100]
                    not_aln_exons_l = [exon_id_s for exon_id_s in gene_id_2_missing_d[gene_name_s]["exon"]
                                       if gene_id_2_missing_d[gene_name_s]["exon"][exon_id_s]["missing"] == 100]
                    # blast_hit: from exon-wide blast
                    blast_hit_exons_l = [exon_id_s for exon_id_s in gene_id_2_missing_d[gene_name_s]["exon"]
                                         if gene_id_2_missing_d[gene_name_s]["exon"][exon_id_s]["blast_hit"] == "hit"]
                    blast_lost_exons_l = [exon_id_s for exon_id_s in gene_id_2_missing_d[gene_name_s]["exon"]
                                          if gene_id_2_missing_d[gene_name_s]["exon"][exon_id_s]["blast_hit"] == "default"]
                    assert len(aln_exons_l) + len(not_aln_exons_l) == num_of_exons_i
                    # Totally missing: 100% missing in genome-wide alignment + no exon blast hit
                    if gene_missing_fl == 100 and len(blast_lost_exons_l) == num_of_exons_i:
                        type_s = "totally_missing"
                    # Exon deletion: not in totally missing, but including missing exon supported by both genome-wide alignment + blast
                    elif len(list(set(not_aln_exons_l) & set(blast_lost_exons_l))) >= 1:
                        type_s = "exon_deletion"
                        lost_exons_l = list(set(not_aln_exons_l) & set(blast_lost_exons_l))
                        gene_nickname_s = gene_name_s.split("=")[0]
                        assert gene_nickname_s not in gene_2_lost_exon_d
                        gene_2_lost_exon_d[gene_nickname_s] = []
                        for lost_exon_id_s in lost_exons_l:
                            lost_exon_bedline_l = exon_id_2_bed_d[lost_exon_id_s]
                            lost_exon_bedline_l[3] = gene_name_s + "|" + lost_exon_bedline_l[3]
                            exon_lost_bed_file_f.write("\t".join(map(str, lost_exon_bedline_l)) + "\n")
                            gene_2_lost_exon_d[gene_nickname_s].append(lost_exon_bedline_l[1:3])
                    else:
                        type_s = "NA"
                    missing_summary_file_f.write("\t".join(map(str,[gene_name_s, num_of_exons_i, ",".join(aln_exons_l), ",".join(not_aln_exons_l),
                                                                    ",".join(blast_hit_exons_l), ",".join(blast_lost_exons_l), type_s])) + "\n")
                pickle.dump(gene_2_lost_exon_d, gene_2_lost_exon_pickle_f, pickle.HIGHEST_PROTOCOL)
        # (2-3) Classify exons by their position: 5'UTR, 1st coding exon, internal coding exon, last coding exon, 3'UTR
        ref_utr_n_cds_bedfile_s = os.path.join(working_dir_s, "lost_utr_n_cds.bed")
        ref_filtered_bed_file_s = os.path.join(result_dir, "gene_structure", "ref", "ref.filter.all_exons.bed")
        softmasked_bed_file_s = os.path.join(basic_input_dir_s, "ref_softmasked.bed")
        if file_check(ref_filtered_bed_file_s) is False:
            print("Warning: gene structure module is needed...gene structure analysis begins.")
            self.gene_structure(mode_s="ref")
        if file_check(softmasked_bed_file_s) is False:
            get_softmasked(fasta_file_s=ref_fasta_file_s, mode_s="ref", work_dir=result_dir)
        position_2_summary_d = {"five-prime-UTR": {"count": 0, "missing": 0}, "FirstCDS": {"count": 0, "missing": 0},
                                "InternalCDS": {"count": 0, "missing": 0}, "LastCDS": {"count": 0, "missing": 0},
                                "three-prime-UTR": {"count": 0, "missing": 0}}
        # a. filter UTR or coding exon records only
        with open(ref_utr_n_cds_bedfile_s, "w") as ref_utr_n_cds_bedfile_f,\
            open(ref_filtered_bed_file_s) as ref_filtered_bed_file_f:
            for line_s in ref_filtered_bed_file_f:
                line_l = line_s.rstrip("\n").split("\t")
                chrom_s, start_i, end_i = line_l[0], int(line_l[1]), int(line_l[2])
                id_s = line_l[3]
                missing_bool = False
                try:
                    type_s,gene_name_s, coding_exon_type_s = id_s.split("|")[0], id_s.split("|")[1], id_s.split("|")[6]
                except IndexError:
                    type_s = id_s.split("|")[0]
                    assert type_s not in ["cds", "three-prime-UTR", "five-prime-UTR"]
                    continue
                if type_s in ["cds", "three-prime-UTR", "five-prime-UTR"]:
                    position_s = id_s.split("|")[5]
                    if coding_exon_type_s in ["more_than_four_coding_exon", "ThreeCodingExon"]:
                        position_2_summary_d[position_s]["count"] += 1
                    # check if exon is missing
                    if gene_name_s in gene_2_lost_exon_d:
                        lost_exon_ranges_l = gene_2_lost_exon_d[gene_name_s]
                        for range_l in lost_exon_ranges_l:
                            range_start_i, range_end_i = range_l
                            if range_start_i <= start_i <= end_i <= range_end_i:
                                missing_bool = True
                        # if there is a missing exon, check position
                        if missing_bool is True:
                            if coding_exon_type_s in ["more_than_four_coding_exon", "ThreeCodingExon"]:
                                position_2_summary_d[position_s]["missing"] += 1
                    missing_info_s = "missing" if missing_bool is True else "not missing"
                    new_id_s = id_s + "|%s" % missing_info_s
                    new_line_l = copy.deepcopy(line_l)
                    new_line_l[3] = new_id_s
                    ref_utr_n_cds_bedfile_f.write("\t".join(new_line_l) + "\n")
        # b. calculate GC and repeat content
        ref_utr_n_cds_tsvfile_s = os.path.join(working_dir_s, "lost_utr_n_cds.tsv")
        os.system("%s intersect -wao -a %s -b %s|" % (bedtools_command_s, ref_utr_n_cds_bedfile_s, softmasked_bed_file_s) +
                  "%s groupby -i stdin -g 1,2,3,4,6 -c 10 -o sum |" % bedtools_command_s +
                  "awk -F'[\"\t\"]' '{OFS=\"\t\"} {print $1,$2,$3,$4\"|\"$6/($3-$2),0,$5}' |" +
                  "%s nuc -fi %s -s -bed stdin |" % (bedtools_command_s, ref_fasta_file_s) +
                  "cut -f4,10,11,15 | grep -v \"4_usercol\" |" +
                  "awk -F'[\"\\t\"\"|\"]' '{OFS=\"\\t\"} " +
                  "{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,($14+$15)/$16}' > %s" % (ref_utr_n_cds_tsvfile_s))
        # c. summary the nummber and missing ratio of exons on each position (UTRs or coding exons)
        ref_lost_ratio_summary_file_s = os.path.join(working_dir_s, "lost_depending_on_exon_position.tsv")
        with open(ref_lost_ratio_summary_file_s,"w") as ref_lost_ratio_summary_file_f:
            ref_lost_ratio_summary_file_f.write("\t".join(["position", "count", "missing", "missing_ratio"]) + "\n")
            for position_s in position_2_summary_d:
                count_i = position_2_summary_d[position_s]["count"]
                missing_i = position_2_summary_d[position_s]["missing"]
                missing_ratio_fl = missing_i/count_i * 100
                ref_lost_ratio_summary_file_f.write("\t".join(map(str,[position_s, count_i, missing_i, missing_ratio_fl])) + "\n")

    def gene_2_bed(self, mode_s):
        """
        function to get coordinate of gene records
        :param mode_s: string, ref or tgt
        :return: X
        """
        if mode_s == "ref":
            gene_recs_l = ref_all_gene_recs_l
        else:
            assert mode_s == "tgt"
            gene_recs_l = tgt_all_gene_recs_l
        working_dir_s = os.path.join(result_dir, "basic_inputs")
        gene_bedfile_s = os.path.join(working_dir_s, "%s.genes.bed" % mode_s)
        gene_bedlines_l = []
        with open(gene_bedfile_s, "w") as gene_bedfile_f:
            for gene_rec in gene_recs_l:
                gene_name_s = gene_id_2_gene_name_d[mode_s][gene_rec.id]
                gene_bed_line_s = "\t".join(map(str, [gene_rec.chrom, gene_rec.start-1, gene_rec.end,
                                                      gene_name_s, 0, gene_rec.strand]))
                gene_bedlines_l.append(gene_bed_line_s)
            gene_bedfile_f.write("\n".join(gene_bedlines_l))
        # for intergenic bed file
        original_gff_s = args.originalGFF
        genes_before_flt_bedfile_s = os.path.join(working_dir_s, "%s.genes_before_flt.bed" % mode_s)
        # for noncoding bed file
        noncoding_genes_bedfile_s = os.path.join(working_dir_s, "%s.noncoding_genes.bed" % mode_s)
        original_db = db_maker(gff_s=original_gff_s, name_s="original.%s" % mode_s,
                               working_dir_s=working_dir_s)
        gene_before_flt_bedlines_l = []
        noncoding_genes_bedlines_l = []
        with open(genes_before_flt_bedfile_s, "w") as genes_before_flt_bedfile_f, \
                open(noncoding_genes_bedfile_s, "w") as noncoding_genes_bedfile_f:
            for gene_rec in original_db.all_features(featuretype="gene"):
                gene_bed_line_s = "\t".join(map(str, [gene_rec.chrom.split(".")[0], gene_rec.start-1, gene_rec.end,
                                                      gene_rec.id, 0, gene_rec.strand]))
                gene_before_flt_bedlines_l.append(gene_bed_line_s)
                if gene_rec.attributes["gene_biotype"][0] in ["lncRNA", "rRNA", "tRNA", "snRNA", "snoRNA"]:
                    noncoding_genes_bedlines_l.append(gene_bed_line_s)
            genes_before_flt_bedfile_f.write("\n".join(gene_before_flt_bedlines_l))
            noncoding_genes_bedfile_f.write("\n".join(noncoding_genes_bedlines_l))


    def intergenic_2_bed(self, mode_s):
        """
        function to get coordinate of gene records
        :param mode_s: string, ref or tgt
        :return: X
        """
        genes_before_flt_bedfile_s = os.path.join(result_dir, "basic_inputs", "%s.genes_before_flt.bed" % mode_s)
        genome_bedfile_s = os.path.join(result_dir, "basic_inputs", "%s.genome.bed" % mode_s)
        intergenic_bedfile_s = os.path.join(result_dir, "basic_inputs", "%s.intergenic.bed" % mode_s)
        os.system("%s sort -i %s | %s merge -i stdin | %s subtract -a %s -b stdin > %s" %
                  (bedtools_command_s, genes_before_flt_bedfile_s,
                   bedtools_command_s,
                   bedtools_command_s, genome_bedfile_s, intergenic_bedfile_s))

    def cds_2_bed(self, mode_s):
        """
        function to get (revised) CDS coordinates
        :param mode_s: string, ref or tgt
        :return: X
        """
        cds_dir_s = os.path.join(result_dir, "CDS", mode_s)
        cds_bedfile_s = os.path.join(result_dir, "basic_inputs", "%s.cds.bed" % mode_s)
        os.system("rm %s" % cds_bedfile_s)
        os.system('find %s -path "*/cds_*.bed" | while read F; do cat ${F} >> %s ; done' % (cds_dir_s, cds_bedfile_s))

    def exon_2_bed(self, mode_s):
        """
        function to get (revised) CDS coordinates
        :param mode_s: string, ref or tgt
        :return: X
        """
        exon_dir_s = os.path.join(result_dir, "CDS", mode_s)
        exon_bedfile_s = os.path.join(result_dir, "basic_inputs", "%s.exon.bed" % mode_s)
        os.system("rm %s" % exon_bedfile_s)
        os.system('find %s -path "*/exon_*.bed" | while read F; do cat ${F} >> %s ; done' % (exon_dir_s, exon_bedfile_s))

    def intron_2_bed(self, mode_s):
        """
        function to get (revised) intron coordinates
        :param mode_s: string, ref or tgt
        :return: X
        """
        ref_gene_name_2_intron_d, ref_intron_recs_l = self.generate_intron_records_from_ref_assembly()
        intron_bedfile_s = os.path.join(result_dir, "basic_inputs", "%s.intron.bed" % mode_s)
        ref_gene_ids_with_intron_l = ref_gene_name_2_intron_d.keys()
        with open(intron_bedfile_s, "w") as intron_bedfile_f:
            for gene_name_s in ref_gene_ids_with_intron_l:
                gene_id_s = gene_name_2_gene_id_d["ref"][gene_name_s]
                gene_chrom_s = refDB[gene_id_s].chrom
                gene_strand_chr = refDB[gene_id_s].strand
                if gene_strand_chr == "+":
                    sorted_intron_recs_l = sorted(ref_gene_name_2_intron_d[gene_name_s],
                                                  key=lambda intron_rec: intron_rec.start)
                else:
                    assert gene_strand_chr == "-"
                    sorted_intron_recs_l = sorted(ref_gene_name_2_intron_d[gene_name_s],
                                                  key=lambda intron_rec: intron_rec.start, reverse=True)
                # for each intron, check intron length + give ID for each intron
                sorted_intron_recs_l = [intron_rec for intron_rec in sorted_intron_recs_l if intron_rec.end-intron_rec.start+1 >= 20]
                intron_order_i = 0
                intron_count_i = len(sorted_intron_recs_l)
                for ref_intron_rec in sorted_intron_recs_l:
                    intron_order_i += 1
                    intron_id_s = "|".join(map(str, [gene_name_s, intron_order_i, intron_count_i]))
                    bed_file_l = [gene_chrom_s, ref_intron_rec.start, ref_intron_rec.end, intron_id_s, 0, gene_strand_chr]
                    intron_bedfile_f.write("\t".join(map(str,bed_file_l))+"\n")

    #############################
    # 3. Splitted gene analysis #
    #############################
    def main_split_genes_analysis(self, fgl_s):
        """
        function to detect split genes (fragmented or intrascaffold split) based on two independent data sources
            (1) CAT result: from function "split_candidates" ==> genes with "possible_split_gene_locations" ID (string, cat_id_s)
                            ==> ID (cat_id_s): (1) gene regions (gene_locus_id_s) & (2) split regions (split_locus_id_s)
            (2) BLAST result: from Exon-wide blast
            ==> THe genes with the blast hits intersecting with both gene and split regions ==> split genes
        :param fgl_s: string, fragmented 0r intrascaffold split
        :return: X
        """
        # basic input variables
        working_dir_s = os.path.join(result_dir, fgl_s)
        dir_check(working_dir_s)
        split_bed_file_s = os.path.join(working_dir_s, "cat_report_%s_genes.bed" % fgl_s)
        blast_bed_file_s = os.path.join(working_dir_s, "blast_hits.bed")
        intersection_bed_file_s = os.path.join(working_dir_s, "intersection.cat_vs_blast.bed")
        summary_of_blast_hit_s = os.path.join(working_dir_s, "summary.cat_vs_blast.txt")
        final_gene_names_file_s = os.path.join(working_dir_s, "%s_genes.txt" % fgl_s)

        # (1) summary of CAT split genes
        # gene_2_split_d: dict, key=gene name(string), value=dictionary {gene_locus_id_s(string):[], split_locus_id_s(string):[]}
        split_gene_names_l, gene_2_split_d = self.split_candidates(fgl_s=fgl_s)

        # (2) summary of blast hit of the split genes
        all_blast_bedline_l = []
        for split_gene_name_s in split_gene_names_l:
            blast_bedline_l = blast_2_bed(gene_name_s=split_gene_name_s)
            all_blast_bedline_l += blast_bedline_l
        with open(blast_bed_file_s, "w") as blast_bed_file_f:
            blast_bed_file_f.write("\n".join(all_blast_bedline_l))

        # (3) intersect CAT and blast hits
        os.system("%s intersect -a %s -b %s -wo > %s" % (bedtools_command_s, split_bed_file_s, blast_bed_file_s, intersection_bed_file_s))

        # (4) check if the intersection is within a single gene
        with open(intersection_bed_file_s, "r") as intersection_bed_file_f:
            for line_s in intersection_bed_file_f:
                line_l = line_s.rstrip("\n").split("\t")
                # cat_id_s = ID of genic region (gene_locus_id_s) or splited region (split_locus_id_s)
                cat_id_s, blast_id_s = line_l[3], line_l[9]
                cat_gene_name_s = cat_id_s.split("^")[0]
                blast_gene_name_s, blast_exon_s = blast_id_s.split("^")[0:2]
                if cat_gene_name_s == blast_gene_name_s:
                    gene_2_split_d[cat_gene_name_s][cat_id_s].append(blast_exon_s)

        # (5) check if a blast hit was found within both original gene & split coordiantes
        final_gene_names_l = []
        with open(summary_of_blast_hit_s,"w") as summary_of_blast_hit_f:
            for gene_name_s in gene_2_split_d.keys():
                summary_line_l = []
                cat_id_with_blast_hit_l = []
                gene_locus_id_s = "%s^gene^1" % gene_name_s  # ID of genic region
                gene_or_split_ids_l = list(gene_2_split_d[gene_name_s].keys())
                summary_line_l.append(gene_name_s)
                for cat_id_s in gene_or_split_ids_l:
                    blast_hits_in_cat_l = gene_2_split_d[gene_name_s][cat_id_s]
                    summary_line_l.append("%s:%s" % (cat_id_s, ",".join(blast_hits_in_cat_l)))
                    # check if splited region includes blast hit
                    is_blast_hit_in_cat_bool = (len(blast_hits_in_cat_l) >= 1)
                    if is_blast_hit_in_cat_bool is True:
                        cat_id_with_blast_hit_l.append(cat_id_s)
                summary_of_blast_hit_f.write("\t".join(summary_line_l)+"\n")
                # check if both gene and split region (from CAT) include one or more blast hit
                split_ids_with_blast_hit_l = copy.deepcopy(cat_id_with_blast_hit_l)
                split_ids_with_blast_hit_l.remove(gene_locus_id_s)
                did_gene_id_include_blast_hit_bool = (gene_locus_id_s in cat_id_with_blast_hit_l)
                did_split_id_include_blast_hit_bool = (len(split_ids_with_blast_hit_l) >= 1)
                if did_split_id_include_blast_hit_bool is True and did_gene_id_include_blast_hit_bool is True:
                    final_gene_names_l.append(gene_name_s)

        # (6) save the cat coordinates & gene names
        final_gene_names_l = list(set(final_gene_names_l))
        with open(final_gene_names_file_s, 'w') as final_gene_names_file_f:
            final_gene_names_file_f.write("\n".join(final_gene_names_l))

    def split_candidates(self, fgl_s):
        """
        function to report cat-based fragmented/intrascaffold split genes
        :param fgl_s: string, fragmented or intrascaffoldsplit
        :return: candidates_genes_l: list, gene names of fragmented or intrascaffoldsplit candidates
                 gene_2_split_d: nested dict, key: gene name, value: dict with key: split it (genename^split_type^order)
        """
        # basic input variables
        working_dir_s = os.path.join(result_dir, fgl_s)
        dir_check(working_dir_s)
        candidates_file_s = os.path.join(working_dir_s, "cat_report_%s_genes.txt" % fgl_s)
        split_bed_file_s = os.path.join(working_dir_s, "cat_report_%s_genes.bed" % fgl_s)
        gene_2_split_d = {}
        # (1) get "possible_split_gene_locations" information reported from cat
        #     + classify them into fragmented or intra-scaffold split
        candidates_genes_l = []
        split_bed_lines_l = []
        for trans_rec in tgt_all_trans_recs_l:
            parent_gene_rec = [gene_rec for gene_rec in tgtDB.parents(id=trans_rec.id, featuretype="gene")][0]
            parent_gene_name_s = get_new_name(gene_rec=parent_gene_rec, mode_s="tgt")
            if "possible_split_gene_locations" in trans_rec.attributes.keys():
                chrom_s = trans_rec.chrom
                split_count_i = 0
                gene_locus_id_s = "%s^gene^1" % parent_gene_name_s
                gene_2_split_d[parent_gene_name_s] = {}
                gene_2_split_d[parent_gene_name_s][gene_locus_id_s] = []
                gene_bedline_s = "\t".join(map(str,[chrom_s, parent_gene_rec.start-1, parent_gene_rec.end, gene_locus_id_s, "0", "."]))
                for split_loci_s in trans_rec.attributes["possible_split_gene_locations"]:
                    type_s = "default"
                    split_chrom_s = split_loci_s.split(":")[0]
                    split_start_s, split_end_s = split_loci_s.split(":")[1].split("-")
                    if (split_chrom_s == chrom_s) and (fgl_s=="intrascaffoldsplit"): # == intrascaffold split
                        type_s = "intrascaffoldsplit"
                    elif (split_chrom_s != chrom_s) and (fgl_s == "fragmented"):  # == fragmented
                        type_s = "fragmented"
                    split_count_i += 1
                    split_id_s = "%s^%s^%d" % (parent_gene_name_s, type_s, split_count_i)
                    split_bedline_s = "\t".join([split_chrom_s, split_start_s, split_end_s, split_id_s, "0", "."])
                    candidates_genes_l.append(parent_gene_name_s)
                    split_bed_lines_l.append(split_bedline_s)
                    gene_2_split_d[parent_gene_name_s][split_id_s] = []
                if split_count_i >= 1:
                    split_bed_lines_l.append(gene_bedline_s)
            else:
                pass
        candidates_genes_l = sorted(list(set(candidates_genes_l)))
        # candidate genes and bed lines
        with open(candidates_file_s, "w") as candidates_file_f, open(split_bed_file_s, "w") as split_bed_file_f:
            candidates_file_f.write("\n".join(candidates_genes_l))
            split_bed_file_f.write("\n".join(split_bed_lines_l))
        return candidates_genes_l, gene_2_split_d

    ###############################
    # 4. (sequence level) mpileup #
    ###############################
    def mpileup(self, info_l):
        """
        function to perform mpileup of samtools
        :param info_l: list, bed_l (6 columns) + mode_s (string, ref for vgp, tgt for old),
                                                fgl_s (string, error type),
                                                flanking_i (integer, flanking sequence length)
        :return: X
        """
        # 1st. get the necessary info for samtools' mpileup: bed, bam, fasta, chr_2_seq dict.
        bed_l, mode_s, fgl_s, flanking_i = info_l[:-3], info_l[-3], info_l[-2], info_l[-1]
        if mode_s == "ref":
            bam_s = ref_bamfile_s
            genome_fasta_s = ref_fasta_file_s
        else:
            assert mode_s == "tgt"
            bam_s = tgt_bamfile_s
            genome_fasta_s = tgt_fasta_file_s
        chrom_s = bed_l[0]
        start_i, end_i = int(bed_l[1]) + 1 - flanking_i, int(bed_l[2]) + flanking_i
        variant_id_s = bed_l[3]
        assert (start_i <= end_i)
        mpileup_block_s = "%s:%d-%d" % (chrom_s, start_i, end_i)
        mpileup_file_s = os.path.join(result_dir, fgl_s, "mpileup", mode_s, "%s.mpileup" % variant_id_s)
        os.system("%s mpileup %s -r %s --reference %s -Bx -s -aa --min-BQ 0 -o %s" %
                  (samtools_command_s, bam_s, mpileup_block_s, genome_fasta_s, mpileup_file_s))

    def summary_mpileup(self, info_l):
        """
        function to compare CAT FS sites with the mpileup result
        :param info_l: nested list, ref_bedline_l, tgt_bedline_l, fgl_s (error type), flanking_i
        :return: return_l, list
        """
        ref_bedline_l, tgt_bedline_l, fgl_s, flanking_i = info_l
        ref_var_id_s, tgt_var_id_s = ref_bedline_l[3], tgt_bedline_l[3]
        ref_mpileup_file_s = os.path.join(result_dir, fgl_s, "mpileup", "ref", "%s.mpileup" % ref_var_id_s)
        tgt_mpileup_file_s = os.path.join(result_dir, fgl_s, "mpileup", "tgt", "%s.mpileup" % tgt_var_id_s)
        return_l = []

        # (1) check each raw mpileup file to count the number of indel/mismatch
        for mode_s, var_id_s, mpileup_file_s in zip(["ref","tgt"], [ref_var_id_s, tgt_var_id_s], [ref_mpileup_file_s, tgt_mpileup_file_s]):
            assert file_check(mpileup_file_s) is True
            with open(mpileup_file_s) as mpileup_file_f:
                mpileups_nl = [line_s.rstrip("\n").split("\t") for line_s in mpileup_file_f.readlines()]
                # (1-1) average_mapped_readnum_i: average number of mapped reads
                # when less than 10 reads were mapped: continue
                average_mapped_readnum_i = sum([int(mpileup_l[3]) for mpileup_l in mpileups_nl])/len(mpileups_nl)
                if average_mapped_readnum_i < 10:
                    for mpileup_l in mpileups_nl:
                        assert len(mpileup_l) == 7
                        chrom_s, locus_s, nucleotide_chr, read_count_s, reads_s = mpileup_l[0:5]
                        return_l.append("\t".join(map(str, [mode_s, var_id_s, chrom_s, locus_s, nucleotide_chr,
                                                            read_count_s, "LESS_THAN_TEN_READS", 0, 0, 0, 0, 0])) + "\n")
                        continue
                # (1-2) block_size_i: to check no missing sites of mpileup result
                if mode_s == "ref":
                    block_size_i = int(ref_bedline_l[2]) - int(ref_bedline_l[1]) + 2 * flanking_i
                else:
                    assert mode_s == "tgt"
                    block_size_i =  int(tgt_bedline_l[2]) - int(tgt_bedline_l[1]) + 2 * flanking_i
                assert len(mpileups_nl) == block_size_i
                # (1-3) "Read" mpileup result
                # for each line, classify whether there are homozygous / heterozygous variants
                for mpileup_l in mpileups_nl:
                    assert len(mpileup_l) == 7
                    chrom_s, locus_s, nucleotide_chr, read_count_s, reads_s = mpileup_l[0:5]
                    read_count_i = int(read_count_s)
                    # if no mapped reads: go to the next site
                    if read_count_i == 0:
                        return_l.append("\t".join(map(str, [mode_s, var_id_s, chrom_s, locus_s, nucleotide_chr,
                                                            read_count_s, "NO_MAPPING_READ", 0, 0, 0, 0, 0])) + "\n")
                        continue
                    number_of_reads_i, number_of_insertion_i, number_of_deletion_i,\
                    deletion_on_site_i, number_of_match_i, number_of_mismatch_i,\
                    error_type_2_count_d = self.read_mpileup(rawreads_s=reads_s, read_count_i=read_count_i)

                    insertion_pct_fl = float(number_of_insertion_i / read_count_i * 100)
                    deletion_pct_fl = float(number_of_deletion_i / read_count_i * 100)
                    deletion_on_site_pct_fl = float(deletion_on_site_i / read_count_i * 100)
                    match_pct_fl = float(number_of_match_i / read_count_i * 100)
                    mismatch_pct_fl = float(number_of_mismatch_i / read_count_i * 100)
                    read_info_2_variant_d = {"insertion": [insertion_pct_fl, "HOMOINS_START"],
                                             "deletion": [deletion_pct_fl, "HOMODEL_START"],
                                             "deletion_on_site": [deletion_on_site_pct_fl, "HOMODEL_MID"],
                                             "match": [match_pct_fl, "HOMOMATCH"],
                                             "mismatch": [mismatch_pct_fl, "HOMOMISMAT"]}
                    variant_type_l = []
                    read_infos_l = ["insertion", "deletion", "deletion_on_site", "match", "mismatch"]
                    for read_info_s in read_infos_l:
                        pct_fl = read_info_2_variant_d[read_info_s][0]
                        if pct_fl >= 80:
                            variant_type_l.append(read_info_2_variant_d[read_info_s][1])
                    variant_type_s = "_".join(variant_type_l)
                    return_l.append("\t".join(map(str,[mode_s, var_id_s, chrom_s, locus_s, nucleotide_chr, read_count_s,
                                                       variant_type_s, insertion_pct_fl, deletion_pct_fl, deletion_on_site_pct_fl,
                                                       match_pct_fl, mismatch_pct_fl,
                                                       list(error_type_2_count_d["insertion"].items()),
                                                       list(error_type_2_count_d["deletion"].items()),
                                                       list(error_type_2_count_d["mismatch"].items())])) + "\n")
        return return_l

    def read_mpileup(self, rawreads_s, read_count_i):
        """
        function to extract how many insertions, deletions, matches, mismatches are in the mpileup results
        from http://samtools.sourceforge.net/pileup.shtml & http://www.htslib.org/doc/samtools-mpileup.html...
            (1) "A pattern `\+[0-9]+[ACGTNacgtn]+' indicates there is an insertion between this reference position and the next reference position"
            (2) "`ACGTN' for a mismatch on the forward strand and `acgtn' for a mismatch on the reverse strand."
            (3) "A symbol `^' marks the start of a read segment which is a contiguous subsequence on the read separated by `N/S/H' CIGAR operations."
            (4) "The ASCII of the character following `^' minus 33 gives the mapping quality."
            (5) "A symbol `$' marks the end of a read segment."
        :param rawreads_s: string, reads string (e.g. ,c$..,,,,,.......,,.,,,.,,.....^T,) from mpileup
        :param read_count_i: integer, reads number (e.g. 31)
        :return:read info list [sumPlaceHolders_i, sumIns_i, sumDel_i, deletionCount_i, matchCount_i, mismatchCount_i]
        """
        # (0) basic inputs
        reads_s = copy.deepcopy(rawreads_s)
        assert (type(reads_s) == str)
        assert (type(read_count_i) == int)
        lenString_i = len(reads_s)
        # if no read mapped ==> return 0 value
        if lenString_i == 0:
            return [0, 0, 0, 0, 0, 0]
        error_type_2_count_d = {"mismatch": {}, "insertion": {}, "deletion": {}}

        # (1) insertion
        # (1-1) get every types of insertion on the loci
        insertion_rexs_l = []
        insertion_count_rex = re.compile("\+[0-9]+")
        insertion_types_l = list(set(re.findall(pattern=insertion_count_rex, string=reads_s)))
        # for each type of insertions (# e.g. +11, +2, ...),
        for insertion_type_s in insertion_types_l:
            insertion_size_i = int(insertion_type_s.replace("+", ""))  # e.g. +11 ==> 11
            # trick: to catch occurence of specific insertion size regardless of the digits of number, using regex
            #        (1) if insertion size =1  ==> insertion_catch_s: [1]
            #        (2) if insertion size =11  ==> insertion_catch_s: [1][1]
            #       in this way, we can handle +1T or +11AAAAAAAAAAA as different types
            insertion_catch_s = "[" + "][".join([x for x in str(insertion_size_i)]) + "]"
            insertion_rex_s = "\+%s[ACGTNacgtn]{%d}" % (insertion_catch_s, insertion_size_i)  # e.g. 11 ==> "+[1][1][ACGTNacgtn]{11}" (+11AAAAAAAAAAA)
            insertion_rex = re.compile(insertion_rex_s)
            insertion_rexs_l.append(insertion_rex)
        # (2) deletion
        # in mpileup, there are two diffrent ways to represent the deletion
        # (2-1) start of deletion (similar with insertion)
        #      : get every types of deletion on the loci
        deletion_rexs_l = []
        deletion_count_rex = re.compile("-[0-9]+")
        deletion_types_l = list(set(re.findall(pattern=deletion_count_rex, string=reads_s)))
        for deletion_type_s in deletion_types_l:
            deletion_size_i = int(deletion_type_s.replace("-", ""))
            deletion_catch_s = "[" + "][".join([x for x in str(deletion_size_i)]) + "]"
            deletion_rex_s = "-%s[ACGTNacgtn]{%d}" % (deletion_catch_s, deletion_size_i)
            deletion_rex = re.compile(deletion_rex_s)
            deletion_rexs_l.append(deletion_rex)

        # "\^[ -~]": regex to catch all printable ascii characters, http://samtools.sourceforge.net/pileup.shtml
        read_start_rex = re.compile("\^[ -~]")
        # trick: in order to catch mismatch, all insertion/deletions should be removed since they include [ACTGNactgn]
        reads_s = str(re.sub(pattern=read_start_rex, repl="", count=1000000000000, string=reads_s))
        number_of_insertion_i = 0
        number_of_deletion_i = 0
        for insertion_rex in insertion_rexs_l:
            insertion_types_l = re.findall(pattern=insertion_rex, string=reads_s)
            number_of_insertion_i += len(insertion_types_l)
            insertion_dedup_l = list(set(insertion_types_l))
            for insertion_type_s in insertion_dedup_l:
                try:
                    error_type_2_count_d["insertion"][insertion_type_s.upper()] += insertion_types_l.count(insertion_type_s)
                except KeyError:
                    error_type_2_count_d["insertion"][insertion_type_s.upper()] = insertion_types_l.count(insertion_type_s)
            reads_s = str(re.sub(pattern=insertion_rex, repl="", count=1000000000000, string=reads_s))
        for deletion_rex in deletion_rexs_l:
            deletion_types_l = re.findall(pattern=deletion_rex, string=reads_s)
            number_of_deletion_i += len(deletion_types_l)
            deletion_dedup_l = list(set(deletion_types_l))
            for deletion_type_s in deletion_dedup_l:
                try:
                    error_type_2_count_d["deletion"][deletion_type_s.upper()] += deletion_types_l.count(deletion_type_s)
                except KeyError:
                    error_type_2_count_d["deletion"][deletion_type_s.upper()] = deletion_types_l.count(deletion_type_s)

            reads_s = str(re.sub(pattern=deletion_rex, repl="", count=1000000000000, string=reads_s))

        number_of_match_i = reads_s.count(".") + reads_s.count(",")
        # since we already removed all the ascii codes following "^", just use simple regex to catch mismatch
        mismatch_rex = re.compile("[ACGTNacgtn]")
        mismatch_types_l = re.findall(pattern=mismatch_rex, string=reads_s)
        number_of_mismatch_i = len(mismatch_types_l)
        mismatch_dedup_l = list(set(mismatch_types_l))
        for mismatch_type_s in mismatch_dedup_l:
            try:
                error_type_2_count_d["mismatch"][mismatch_type_s.upper()] += mismatch_types_l.count(mismatch_type_s)
            except KeyError:
                error_type_2_count_d["mismatch"][mismatch_type_s.upper()] = mismatch_types_l.count(mismatch_type_s)

        # *: deletion on the site
        # (2-2) deletion on the site: marked as *
        deletion_on_site_i = reads_s.count("*")
        number_of_reads_i = number_of_match_i + number_of_mismatch_i + deletion_on_site_i
        assert (read_count_i == number_of_reads_i)
        return [number_of_reads_i, number_of_insertion_i, number_of_deletion_i,
                deletion_on_site_i, number_of_match_i, number_of_mismatch_i,
                error_type_2_count_d]

    def bed_for_mpileup(self,fgl_s):
        """
        function to prepare bed lines from the ref & tgt bed file
        :param fgl_s: types of fgl
        :return: list of bedlines
        """
        ref_bedlines_nl, tgt_bedlines_nl = [], []
        ref_fgl_ids_l, tgt_fgl_ids_l = [], []
        if fgl_s == "frameshift":
            ref_bedfile_s = os.path.join(result_dir, "frameshift", "psl_2_frameshift_candidates", "ref_less_than_10bp.bed")
            tgt_bedfile_s = os.path.join(result_dir, "frameshift", "psl_2_frameshift_candidates", "tgt_less_than_10bp.bed")
            with open(ref_bedfile_s, "r") as ref_bedfile_f, open(tgt_bedfile_s,"r") as tgt_bedfile_f:
                for line_s in ref_bedfile_f:
                    line_l = line_s.rstrip("\n").split("\t")
                    ref_bedlines_nl.append(line_l)
                for line_s in tgt_bedfile_f:
                    line_l = line_s.rstrip("\n").split("\t")
                    tgt_bedlines_nl.append(line_l)
            # (2-1) sort with the FS ID (refCoding~ or tgtCoding~)
            ref_bedlines_nl = sorted(ref_bedlines_nl, key=lambda x: "gene-" + x[3].split("gene-")[-1])
            tgt_bedlines_nl = sorted(tgt_bedlines_nl, key=lambda x: "gene-" + x[3].split("gene-")[-1])
            ref_fgl_ids_l = ["gene-" + x[3].split("gene-")[-1] for x in ref_bedlines_nl]
            tgt_fgl_ids_l = ["gene-" + x[3].split("gene-")[-1] for x in tgt_bedlines_nl]

        if fgl_s == "prematurestopcodon":
            ref_bedlines_nl, tgt_bedlines_nl = [], []
            ref_bedfile_s = os.path.join(result_dir, "prematurestopcodon", "ref_premature_stopcodon_candidates.bed")
            tgt_bedfile_s = os.path.join(result_dir, "prematurestopcodon", "tgt_premature_stopcodon_candidates.bed")
            with open(ref_bedfile_s,"r") as ref_bedfile_f, open(tgt_bedfile_s,"r") as tgt_bedfile_f:
                for line_s in ref_bedfile_f:
                    line_l = line_s.rstrip("\n").split("\t")
                    ref_bedlines_nl.append(line_l)
                for line_s in tgt_bedfile_f:
                    line_l = line_s.rstrip("\n").split("\t")
                    tgt_bedlines_nl.append(line_l)
            # (2-1) sort with the FS ID (refCoding~ or tgtCoding~)
            ref_bedlines_nl = sorted(ref_bedlines_nl, key=lambda x: x[3])  # x[3] ==> premature stop codon id
            tgt_bedlines_nl = sorted(tgt_bedlines_nl, key=lambda x: x[3])
            ref_fgl_ids_l = [x[3] for x in ref_bedlines_nl]
            tgt_fgl_ids_l = [x[3] for x in tgt_bedlines_nl]

        if "intronexonjunctiondisruption" in fgl_s:
            ref_bedlines_nl, tgt_bedlines_nl = [], []
            ref_bedfile_s = os.path.join(result_dir, "intronexonjunctiondisruption", "ref.splicing_junction_disruption_candidates.bed")
            tgt_bedfile_s = os.path.join(result_dir, "intronexonjunctiondisruption", "tgt.splicing_junction_disruption_candidates.bed")
            with open(ref_bedfile_s,"r") as ref_bedfile_f, open(tgt_bedfile_s,"r") as tgt_bedfile_f:
                for line_s in ref_bedfile_f:
                    line_l = line_s.rstrip("\n").split("\t")
                    ref_bedlines_nl.append(line_l)
                for line_s in tgt_bedfile_f:
                    line_l = line_s.rstrip("\n").split("\t")
                    tgt_bedlines_nl.append(line_l)
            # (2-1) sort with the FS ID (refCoding~ or tgtCoding~)
            ref_bedlines_nl = sorted(ref_bedlines_nl, key=lambda x: x[3])  # x[3] ==> junction disruption id
            tgt_bedlines_nl = sorted(tgt_bedlines_nl, key=lambda x: x[3])
            ref_fgl_ids_l = ["_".join(x[3].split("_")[3:]) for x in ref_bedlines_nl]
            tgt_fgl_ids_l = ["_".join(x[3].split("_")[3:]) for x in tgt_bedlines_nl]
        assert (ref_fgl_ids_l == tgt_fgl_ids_l)
        return ref_bedlines_nl, tgt_bedlines_nl

    def mp_main_mpileup_of_sequence_level_fgl(self, fgl_s, num_processors_i=15, flanking_size_i=2):
        """
        main function for multiprocessing mpileup; perform mpileup, read and summary the variants, and filter candidates
        :param fgl_s: string, type of false gene loss: one of "frameshift", "prematurestopcodon", "intronexonjunctiondisruption"
        :param num_processors_i: integer, number of processors
        :param flanking_size_i: integer, length of flanking sequences
        :return:
        """
        # (1) make list of bed-formated lines
        ref_bedlines_nl, tgt_bedlines_nl = self.bed_for_mpileup(fgl_s)

        # (2-1) (ref:VGP) perform mpileup with false gene loss loci +- neighbor size
        ref_mpileup_dir_s = os.path.join(result_dir, fgl_s, "mpileup", "ref")
        dir_check(ref_mpileup_dir_s)
        ref_mpileups_l = os.listdir(ref_mpileup_dir_s)
        if len(ref_mpileups_l) == 0:
            ref_mpileup_nl = [ref_bedline_l + ["ref", fgl_s, flanking_size_i] for ref_bedline_l in ref_bedlines_nl]
            with Pool(num_processors_i) as p:
                out = p.map(self.mpileup, ref_mpileup_nl)

        # (2-2) (tgt:Prior) perform mpileup with false gene loss loci +- neighbor size
        tgt_mpileup_dir_s = os.path.join(result_dir, fgl_s, "mpileup", "tgt")
        dir_check(tgt_mpileup_dir_s)
        tgt_mpileups_l = os.listdir(tgt_mpileup_dir_s)
        if len(tgt_mpileups_l) == 0:
            tgt_mpileup_nl = [tgt_bedline_l + ["tgt",fgl_s, flanking_size_i] for tgt_bedline_l in tgt_bedlines_nl]
            with Pool(num_processors_i) as p:
                out = p.map(self.mpileup, tgt_mpileup_nl)

        # (2-3) (locus level) count number of insertion, deletion, and mismatch
        mpileup_summary_header_l = ["mode", "variant_id", "chrom", "variant_loci", "ref_nt", "read_count",
                                    "mpileup_type", "insertion(%)", "deletion(%)", "deletion_mark(%)",
                                    "match(%)", "mismatch(%)", "insertion_types", "deletion_types", "mismatch_types"]
        mpileup_summaryfile_s = os.path.join(result_dir, fgl_s, "mpileup_summary.tabs")
        if file_check(mpileup_summaryfile_s) is False:
            assert len(ref_bedlines_nl) == len(tgt_bedlines_nl)
            ref_n_tgt_bedlines_nl = [list(ref_n_tgt_bedlines_t) + [fgl_s, flanking_size_i]
                                     for ref_n_tgt_bedlines_t in zip(ref_bedlines_nl, tgt_bedlines_nl)]
            with Pool(num_processors_i) as p:
                out = p.map(self.summary_mpileup, ref_n_tgt_bedlines_nl)
            with open(mpileup_summaryfile_s, "w") as mpileup_summaryfile_f:
                mpileup_summaryfile_f.write("\t".join(mpileup_summary_header_l) + "\n")
                for summary_mpileup_lines_l in out:
                    if len(summary_mpileup_lines_l) >= 1:
                        for summary_mpileup_line_s in summary_mpileup_lines_l:
                            mpileup_summaryfile_f.write(summary_mpileup_line_s)

        # (3) (variant level) summarize the number of insertion, deletion, and mismatch
        var_id_2_mpileup_line_d = {}
        locus_summaryfile_s = os.path.join(result_dir, fgl_s, "locus_wise.mpileup_summary.tabs")
        with open(mpileup_summaryfile_s) as mpileup_summaryfile_f, open(locus_summaryfile_s, "w") as locus_summaryfile_f:
            # (3-1) var_id_2_mpileup_line_d: key -> variant id, value -> dict (key -> "ref" or "tgt", value -> mpileup summary line)
            mpileup_summaryfile_f.readline()
            for line_s in mpileup_summaryfile_f:
                line_l = line_s.rstrip("\n").split("\t")
                mode_s, variant_id_s, mpileup_type_s = line_l[0], "gene-" + line_l[1].split("gene-")[-1], line_l[6]
                try:
                    var_id_2_mpileup_line_d[variant_id_s][mode_s].append(mpileup_type_s)
                except KeyError:
                    var_id_2_mpileup_line_d[variant_id_s] = {"ref": [], "tgt": []}
                    var_id_2_mpileup_line_d[variant_id_s][mode_s].append(mpileup_type_s)
            # (3-2) for variant id: summarize error type
            classes_l, locus_summary_l, vgp_error_ids_l, fgl_ids_l = [], [], [], []
            variant_ids_l = sorted(list(var_id_2_mpileup_line_d.keys()))
            for variant_id_s in variant_ids_l:
                ref_mpileup_types_l = var_id_2_mpileup_line_d[variant_id_s].get("ref", [])
                tgt_mpileup_types_l = var_id_2_mpileup_line_d[variant_id_s].get("tgt", [])
                ref_all_types_s = ",".join(list(set(ref_mpileup_types_l)))
                tgt_all_types_s = ",".join(list(set(tgt_mpileup_types_l)))
                # (3-3) if both VGP and prior assemblies have mapped reads:
                #       check if either assembly include any of errors ("HOMOINS", "HOMODEL", "HOMOMISMAT", "NO_MAPPING_READ", "LESS_THAN_TEN_READS")
                if len(ref_mpileup_types_l) >= 1 and len(tgt_mpileup_types_l) >= 1:
                    # (3-4) check any potential error flag was found in any of assemblies
                    assembly_error_types_l = []
                    ref_bool, tgt_bool = True, True
                    error_types_l = ["HOMOINS", "HOMODEL", "HOMOMISMAT", "NO_MAPPING_READ", "LESS_THAN_TEN_READS"]
                    for error_s in error_types_l:
                        if error_s in ref_all_types_s:
                            assembly_error_types_l.append("ref;" + error_s)
                            ref_bool = False
                        if error_s in tgt_all_types_s:
                            assembly_error_types_l.append("tgt;" + error_s)
                            tgt_bool = False
                    # (3-4) based on the errors in assemblies, classify error type of each locus
                    variant_id_error_type_s = ",".join(assembly_error_types_l)
                    class_s = classify_variant_id_error(ref_bool=ref_bool, tgt_bool=tgt_bool, fgl_s=fgl_s,
                                                        variant_id_error_type_s=variant_id_error_type_s, variant_id_s=variant_id_s)
                    classes_l.append(class_s)
                    if class_s == "vgp_error":
                        vgp_error_ids_l.append(variant_id_s)
                    if class_s == "previous_error":
                        fgl_ids_l.append(variant_id_s)
                    locus_summary_l.append("\t".join(map(str, [variant_id_s, variant_id_error_type_s, class_s,
                                                               ref_all_types_s, tgt_all_types_s, ",".join(ref_mpileup_types_l),",".join(tgt_mpileup_types_l)])))
            class_stats_l = [class_s + ":" + str(classes_l.count(class_s)) for class_s in ["vgp_error", "previous_error", "no_errors", "filter_out"]]
            locus_summaryfile_f.write("#" + ",".join(class_stats_l) + "\n")
            locus_summaryfile_f.write("#" +"\t".join(map(str,["id","errors","class","ref_all","tgt_all","ref_mpileup","tgt_mpileup"]))+"\n")
            locus_summaryfile_f.write("\n".join(locus_summary_l))

        # (4) final summary of filtered bed files
        ref_flt_bedfile_s = os.path.join(result_dir, fgl_s, "ref_filtered_by_mpileup.bed")
        tgt_flt_bedfile_s = os.path.join(result_dir, fgl_s, "tgt_filtered_by_mpileup.bed")
        ref_bedlines_nl, tgt_bedlines_nl = self.bed_for_mpileup(fgl_s)
        assert len(ref_bedlines_nl) == len(tgt_bedlines_nl)
        with open(ref_flt_bedfile_s, "w") as ref_flt_bedfile_f, open(tgt_flt_bedfile_s, "w") as tgt_flt_bedfile_f:
            ref_flt_bedlines_l, tgt_flt_bedlines_l = [], []
            for i in range(0, len(ref_bedlines_nl)):
                ref_bedline_l, tgt_bedline_l = ref_bedlines_nl[i], tgt_bedlines_nl[i]
                ref_var_id_s = "gene-" + ref_bedline_l[3].split("gene-")[-1]
                tgt_var_id_s = "gene-" + tgt_bedline_l[3].split("gene-")[-1]
                assert ref_var_id_s == tgt_var_id_s
                if ref_var_id_s in fgl_ids_l:
                    ref_flt_bedlines_l.append("\t".join(ref_bedline_l))
                    tgt_flt_bedlines_l.append("\t".join(tgt_bedline_l))
            ref_flt_bedfile_f.write("\n".join(ref_flt_bedlines_l))
            tgt_flt_bedfile_f.write("\n".join(tgt_flt_bedlines_l))

    def mp_summary_mpileup_whole_genome(self, num_processors_i):
        """
        perform multiprocess-search of summarizing mpileup from the whole genome sequences
        :return: X
        """
        working_dir_s = os.path.join(result_dir, "indel")
        dir_check(working_dir_s)
        for mode_s in [ref_s, tgt_s]:
            mpileup_summary_file_s = os.path.join(working_dir_s, "%s.mplieup.indel.bed" % mode_s)
            # (1) summary raw mpileup file
            if file_check(mpileup_summary_file_s) is False:
                mpileup_file_s = "PATH_TO_MPILEUP_FILE"
                if file_check(mpileup_file_s) is False:
                    os.system("%s mpileup %s --reference %s -Bx -s -aa --min-BQ 0 -o %s" %
                              (samtools_command_s, ref_bamfile_s, ref_fasta_file_s, mpileup_file_s))
                assert file_check(mpileup_file_s) is True
                # chunk size =  10,000,000 lines
                with open(mpileup_file_s) as mpileup_file_f, open(mpileup_summary_file_s,"w") as mpileup_indel_file_f:
                    idx_i = 0
                    chunk_size_i = 10000000
                    while True:
                        next_n_lines = list(islice(mpileup_file_f, chunk_size_i))
                        idx_i += 1
                        print("[Read mpileup] ...... %d0Mb" % idx_i)
                        if not next_n_lines:
                            break
                        with Pool(num_processors_i) as p:
                            out = p.map(self.summary_mpileup_for_indel, next_n_lines)
                        for out_s in out:
                            if len(out_s) > 0:
                                mpileup_indel_file_f.write(out_s)
            # (1-1) select loci with more than 1 seq.difference supported by ten or more reads
            #       column 4: read count, column 5~8: num. of reads with insertion, deletion, deletion with mark (*), and mismatch
            loci_with_seq_diff_file_s = os.path.join(working_dir_s, "%s.1_seq_difference.10reads.txt" % mode_s)
            if file_check(loci_with_seq_diff_file_s) is False:
                os.system("awk -F'[\"\\t\"]' '{OFS=\"\\t\"}" +
                          "{if ($4 >= 10 && ($5+$6+$7+$8) >= 1) {print $0}}' %s > %s" % (mpileup_summary_file_s, loci_with_seq_diff_file_s))
            # (2) collect the sites with more than 20% of sequence difference in reads + more than 10 reads (heterozygous allele or assembly error)
            if mode_s == ref_s:
                to_be_excluded_file_s = os.path.join(working_dir_s, "%s.mplieup.exclude.bed" % mode_s)
                false_indel_file_s = os.path.join(working_dir_s, "false_indel.bed")
                false_snp_file_s = os.path.join(working_dir_s, "false_snp.bed")
                if file_check(false_indel_file_s) is False or file_check(false_snp_file_s) is False:
                    # column 2: start coordinates, column 4: read count
                    # column 9~12: % of reads with insertion, deletion, deletion with mark (*), and mismatch
                    os.system("awk -F'[\"\\t\"]' '{OFS=\"\\t\"}" +
                              "{if ( ($2>2) && ($4>=10) && ( ($9 >= 20) || ($10 >= 20) || ($11 >= 20) || ($12 >= 20) ) ) {print $1, $2-2, $3+2}}' %s > %s" %
                              (mpileup_summary_file_s, to_be_excluded_file_s))
                    # exclude loci with more than 10 reads & 20% or more of sequence differences (flanking size=+-2bp)
                    if file_check(false_indel_file_s) is False:
                        os.system("%s intersect -v -a %s -b %s > %s" %
                                  (bedtools_command_s, indel_from_vg_file_s, to_be_excluded_file_s, false_indel_file_s))
                    if file_check(false_snp_file_s) is False:
                        os.system("%s intersect -v -a %s -b %s > %s" %
                                  (bedtools_command_s, snp_from_vg_file_s, to_be_excluded_file_s, false_snp_file_s))

    def summary_mpileup_for_indel(self, line_s):
        """
        function to summary mpileup result against whole genome and count number of ins, del, and mismatch
        :param line_s: string, mpileup result
        :return:
        """
        return_s = ""
        mpileup_l = line_s.rstrip("\n").split("\t")
        assert len(mpileup_l) == 7
        chrom_s, locus_s, nucleotide_chr, read_count_s, reads_s = mpileup_l[0:5]
        read_count_i = int(read_count_s)
        # if no mapped reads on the site: go to the next site
        if read_count_i == 0:
            pass
        else:
            # ins_pct_fl: % of insertion (marked like "+3T" in mpileup)
            # del_pct_fl: % of deletion (marked like "-1A" in mpileup)
            # del_mark_pct_fl: % of deletion on the site (marked as * in mpileup)
            # mismatch_pct_fl: % of mismatch on the site (marked as other nucleotide in mpileup)
            num_reads_i, num_ins_i, num_del_i, num_del_mark_i, \
            num_match_i, num_mismatch_i, error_type_2_count_d = self.read_mpileup(rawreads_s=reads_s, read_count_i=read_count_i)
            ins_pct_fl = float(num_ins_i / read_count_i * 100)
            del_pct_fl = float(num_del_i / read_count_i * 100)
            del_mark_pct_fl = float(num_del_mark_i / read_count_i * 100)
            mismatch_pct_fl = float(num_mismatch_i / read_count_i * 100)
            return_s = "\t".join(map(str, [chrom_s, locus_s, int(locus_s) + 1,
                                           read_count_i,num_ins_i, num_del_i, num_del_mark_i, num_mismatch_i,
                                           ins_pct_fl, del_pct_fl, del_mark_pct_fl, mismatch_pct_fl])) + "\n"
        return return_s

    ###############################
    # 5. Frameshift gene analysis #
    ###############################
    def mp_blat(self, num_processors_i):
        """
        multiprocess function of blat
        :param num_processors_i: integer, number of multiprocess cores
        :return: X
        """
        dir_check(os.path.join(result_dir, "blat_dna", "db"))
        dir_check(os.path.join(result_dir, "blat_dna", "query"))
        dir_check(os.path.join(result_dir, "blat_dna", "out"))
        trans_recs_l = tgt_all_trans_recs_l
        # 1. prepare DB (refCDSs) and query (tgtCDSs) sequences
        summary_of_preprocessing_file_s = os.path.join(os.path.join(result_dir, "blat_dna", "preprocessing.summary.tsv"))
        with Pool(num_processors_i) as p:
            out = p.map(blat_preprocessing, trans_recs_l)
            with open(summary_of_preprocessing_file_s, "w") as summary_of_preprocessing_file_f:
                for parent_gene_name_s, blat_exception_bool in out:
                    summary_of_preprocessing_file_f.write("%s\t%s\n" % (parent_gene_name_s, str(blat_exception_bool)))
        p.close()
        p.join()
        # 2. perform blat
        with Pool(num_processors_i) as p:
            out = p.map(blat, trans_recs_l)
        p.close()
        p.join()
        # 3. summary the blat results
        self.blat_2_summary_psl()

    def blat_2_summary_psl(self):
        """
        function to summary the blat result as the way similar with CAT
        :param mode_s: string, cds or dna (according to blat mode, currenty only dna is available)
        :return: X
        """
        blat_pslfiles_l = os.listdir(os.path.join(result_dir, "blat_dna", "out"))
        blat_pslfiles_l = list(filter(lambda x: x.endswith("psl"), blat_pslfiles_l))
        blat_pslfiles_l = [os.path.join(result_dir, "blat_dna", "out", x) for x in blat_pslfiles_l]
        blat_summary_pslfile_s = os.path.join(result_dir, "blat_dna", "pslConsensus.psl")
        with open(blat_summary_pslfile_s, "w") as blat_summary_pslfile_f:
            for blat_pslfile_s in blat_pslfiles_l:
                # when there is 0-size blat result file --> exclude
                if os.stat(blat_pslfile_s).st_size == 0:
                    continue
                # if multiple blat alignment
                #    ==> pick the alignment with the highest coverage (sum of match / querysize)
                else:
                    with open(blat_pslfile_s) as blat_pslfile_f:
                        lines_nl = []
                        for line_s in blat_pslfile_f:
                            line_l = line_s.rstrip("\n").split("\t")
                            match_i, mismatch_i, repmatch_i = int(line_l[0]), int(line_l[1]), int(line_l[2])
                            sum_of_match_i = match_i + mismatch_i + repmatch_i
                            query_size_i = int(line_l[10])
                            coverage_fl = round(float(sum_of_match_i / query_size_i), 5)
                            line_l.append(coverage_fl)
                            lines_nl.append(line_l)
                        # highest coverage hit would be regared as "longest_psl_result_s"
                        longest_psl_result_s = sorted(lines_nl, key=lambda x: -x[-1])[0]
                        blat_summary_pslfile_f.write("\t".join(longest_psl_result_s[:-1]) + "\n")  # exclude the coverage (coverage_fl)

    def psl_2_frameshift(self):
        """
        detect frameshift based on blat result file, psl
        :return: dict, aln_id_2_psl_summary_d ==> key: aln_id, value:{"tgt": {"start": tgtFsStarts_l, "end": tgtFsEnds_l},
                                                                      "ref": {"start": refFsStarts_l, "end": refFsEnds_l}}
        """
        # (1) check if psl file from dna-level blat result is available, or perform multiprocess blat
        aln_ids_l = sorted(list(gene_id_2_aln_id_d.values()))
        psl_summary_file_s = os.path.join(result_dir, "blat_dna", "pslConsensus.psl")
        if os.path.isfile(psl_summary_file_s) is False:
            self.mp_blat(num_processors_i=process_num_i)
        # (2) based on psl from (1), extract FS inducing indels ==> not multiple of 3, and there should be no gap on other sequence
        aln_id_2_psl_summary_d = {}
        with open(psl_summary_file_s) as psl_summary_file_f:
            for line_s in psl_summary_file_f.readlines():
                line_l = line_s.rstrip("\n").split("\t")
                aln_id_s = line_l[9]
                assert aln_id_s in aln_ids_l
                block_sizes_l = [int(x) for x in line_l[18].split(",")[:-1]]
                tgt_block_starts_l = [int(x) for x in line_l[19].split(",")[:-1]] # 0-based
                tgt_block_ends_l = [sum(x) for x in zip(block_sizes_l, tgt_block_starts_l)] # 0-based
                ref_block_starts_l = [int(x) for x in line_l[20].split(",")[:-1]] # 0-based
                ref_block_ends_l = [sum(x) for x in zip(block_sizes_l, ref_block_starts_l)] # 0-based
                assert (len(block_sizes_l) == len(tgt_block_starts_l) == len(ref_block_starts_l))
                tgt_frameshift_starts_l, tgt_frameshift_ends_l = [], []
                ref_frameshift_starts_l, ref_frameshift_ends_l = [], []
                # (2-1). check: if more than one block? --> indels
                if len(block_sizes_l) > 1:
                    for i in range(0, len(block_sizes_l) - 1):
                        tgt_indel_size_i = tgt_block_starts_l[i + 1] - tgt_block_ends_l[i]
                        ref_indel_size_i = ref_block_starts_l[i + 1] - ref_block_ends_l[i]
                        assert not (tgt_indel_size_i == ref_indel_size_i == 0)
                        assert (tgt_indel_size_i >= 0 and ref_indel_size_i >= 0)
                        # (2-2). check: if there is an indel whose size is not multiple of three --> indel that can lead to frameshift
                        tgt_frameshift_insertion_bool = (tgt_indel_size_i % 3 != 0 and ref_indel_size_i == 0)
                        ref_frameshift_insertion_bool = (tgt_indel_size_i == 0 and ref_indel_size_i % 3 != 0)
                        if tgt_frameshift_insertion_bool is True or ref_frameshift_insertion_bool is True:
                            tgt_frameshift_starts_l.append(tgt_block_ends_l[i])
                            tgt_frameshift_ends_l.append(tgt_block_starts_l[i + 1])
                            ref_frameshift_starts_l.append(ref_block_ends_l[i])
                            ref_frameshift_ends_l.append(ref_block_starts_l[i + 1])
                # (2-3). if there is at least one frameshift indel --> Save
                if (len(tgt_frameshift_starts_l) + len(tgt_frameshift_ends_l) + len(ref_frameshift_starts_l) + len(ref_frameshift_ends_l)) >= 1:
                    assert (len(tgt_frameshift_starts_l) == len(tgt_frameshift_ends_l) == len(ref_frameshift_starts_l) == len(ref_frameshift_ends_l))
                    aln_id_2_psl_summary_d[aln_id_s] = {"tgt": {"start": tgt_frameshift_starts_l, "end": tgt_frameshift_ends_l},
                                                        "ref": {"start": ref_frameshift_starts_l, "end": ref_frameshift_ends_l}}
        return aln_id_2_psl_summary_d

    def frameshift_2_genomic(self, varInfo_l):
        """
        funciton to get the genome coordinate of frameshift-inducing indels
        :param varInfo_l: [alignmentID, {"tgt": {"start": tgtFsStarts_l, "end": tgtFsEnds_l},"ref": {"start": refFsStarts_l, "end": refFsEnds_l}}]
        :return:
        """
        aln_id_s = varInfo_l[0]
        gene_name_s = aln_id_2_gene_id_d[aln_id_s]
        assembly_2_indels_d = varInfo_l[1]
        frameshift_bedline_d = {"ref": [], "tgt": []}
        for mode_s in ["ref", "tgt"]:
            count_i = 0
            cds_frameshift_starts_l = assembly_2_indels_d[mode_s]["start"]
            cds_frameshift_ends_l = assembly_2_indels_d[mode_s]["end"]
            assert len(cds_frameshift_starts_l) == len(cds_frameshift_ends_l)
            # (1) frameshift inducing indels ==> get +- 1bp sequences
            for cds_frameshift_loci_t in zip(cds_frameshift_starts_l, cds_frameshift_ends_l):
                count_i += 1
                cds_frameshift_start_i = cds_frameshift_loci_t[0]  # +1 for 0-based --> 1-based, -1 for including 1 left base
                cds_frameshift_end_i = cds_frameshift_loci_t[1] + 1  # +0 for 0-based --> 1-based, +1 for including 1 right base
                # (2) convert coding sequence coordinates to genome coordinates
                # variant_id_s = ["ref" or "tgt"]^"name of a gene"^["CodingInsertion" or "CodingDeletion"]_"count_i"
                variant_id_s = "^".join(["frameshift", gene_name_s, str(count_i)])
                dir_check(os.path.join(result_dir, "frameshift", "psl_2_frameshift_candidates", "gene_wise", variant_id_s))
                cds_fasta_file_s = os.path.join(result_dir, "CDS", mode_s, gene_name_s, "%s.fasta" % gene_name_s)
                var_bedfile_s = os.path.join(result_dir, "frameshift", "psl_2_frameshift_candidates", "gene_wise", variant_id_s, "%s.bed" % variant_id_s)
                assert os.path.isfile(cds_fasta_file_s) is True
                with open(cds_fasta_file_s) as each_exon_fasta_file_f, open(var_bedfile_s,"w") as var_bedfile_f:
                    # a. cds_lengths_l: length of sorted coding exons (5' --> 3')
                    headers_l, cds_seqs_l = [], []
                    for values in FastaIO.SimpleFastaParser(each_exon_fasta_file_f):
                        headers_l.append(values[0])
                        cds_seqs_l.append(values[1])
                    cds_loci_l = [x.split(":")[-1].split("(")[0] for x in headers_l]
                    cds_lengths_l = [abs(int(cds_locus_s.split("-")[1]) - int(cds_locus_s.split("-")[0])) for cds_locus_s in cds_loci_l]
                    assert cds_lengths_l == [len(cds_seq_s) for cds_seq_s in cds_seqs_l]
                    chr_s, sum_i, strand_str = headers_l[0].split(":")[-2], sum(cds_lengths_l), headers_l[0][-3:].strip("(").strip(")")
                    if strand_str == "-":
                        cds_lengths_l.reverse()
                        cds_loci_l.reverse()
                    # b. cumulative_cds_lengths_l: cumulative length of sorted coding exons (5' --> 3')
                    cumulative_cds_lengths_l = []
                    if len(cds_lengths_l) == 1:  # single coding exon:
                        cumulative_cds_lengths_l, start_cds_order_i, end_cds_order_i = cds_lengths_l, 0, 0
                    else:  # multiple coding exon:
                        for i in range(0, len(cds_lengths_l)):
                            cumulative_cds_lengths_l.append(sum(cds_lengths_l[0:i + 1]))
                        # insert (0,0): for multiple exons
                        cumulative_cds_lengths_l.insert(0, 0)
                        # check the order of coding exon including cds_frameshift_start_i
                        start_cumulative_cds_lengths_l = copy.deepcopy(cumulative_cds_lengths_l)
                        start_cumulative_cds_lengths_l.append(cds_frameshift_start_i)
                        start_cumulative_cds_lengths_l = sorted(start_cumulative_cds_lengths_l)
                        start_cds_order_i = start_cumulative_cds_lengths_l.index(cds_frameshift_start_i) - 1
                        # check the order of coding exon including cds_frameshift_end_i
                        end_cumulative_cds_lengths_l = copy.deepcopy(cumulative_cds_lengths_l)
                        end_cumulative_cds_lengths_l.append(cds_frameshift_end_i)
                        end_cumulative_cds_lengths_l = sorted(end_cumulative_cds_lengths_l)
                        end_cds_order_i = end_cumulative_cds_lengths_l.index(cds_frameshift_end_i) - 1

                    # (4) genomic coordinates
                    if strand_str == "+":
                        exon_with_frameshift_start_i, exon_with_frameshift_end_i = sorted([int(cds_loci_l[start_cds_order_i].split("-")[0]),int(cds_loci_l[end_cds_order_i].split("-")[0])])
                        if len(cds_lengths_l) > 1:
                            # ==> since adding (0, 0), it returns the summed length of coding exons before the exon with frameshift
                            exon_with_frameshift_start_i -= cumulative_cds_lengths_l[start_cds_order_i]
                            exon_with_frameshift_end_i -= cumulative_cds_lengths_l[end_cds_order_i]
                        start_i = exon_with_frameshift_start_i + cds_frameshift_start_i
                        end_i = exon_with_frameshift_end_i + cds_frameshift_end_i
                    # if direction is (-) ==> we should add 1 for start_i and end_i since...
                    else:
                        assert strand_str == "-"
                        exon_with_frameshift_start_i, exon_with_frameshift_end_i = sorted([int(cds_loci_l[start_cds_order_i].split("-")[1]),int(cds_loci_l[end_cds_order_i].split("-")[1])])
                        if len(cds_lengths_l) > 1:
                            # ==> since add of 0, it would return the sum of the length of exons before the FS exon
                            exon_with_frameshift_start_i += cumulative_cds_lengths_l[end_cds_order_i]
                            exon_with_frameshift_end_i += cumulative_cds_lengths_l[start_cds_order_i]
                        start_i = exon_with_frameshift_start_i - cds_frameshift_end_i + 1
                        end_i = exon_with_frameshift_end_i - cds_frameshift_start_i + 1
                    # start-1: due to 1-based to 0-based bed file again
                    frameshift_bedline_s = "\t".join(map(str, [chr_s, start_i - 1, end_i, variant_id_s, strand_str]))
                    var_bedfile_f.write(frameshift_bedline_s)
                    frameshift_bedline_d[mode_s].append(frameshift_bedline_s)
        return frameshift_bedline_d

    def mp_frameshift_2_genome(self, num_processors_i):
        """
        function for .
        :param num_processors_i: int, number of processors to use
        :return: X
        """
        # (1) working directory
        dir_check(os.path.join(result_dir, "frameshift", "psl_2_frameshift_candidates", "gene_wise"))
        frameshift_size_file_s = os.path.join(result_dir, "frameshift", "psl_2_frameshift_candidates", "frameshift_size.tsv")
        ref_less_than_10bp_bedfile_s = os.path.join(result_dir, "frameshift", "psl_2_frameshift_candidates", "ref_less_than_10bp.bed")
        tgt_less_than_10bp_bedfile_s = os.path.join(result_dir, "frameshift", "psl_2_frameshift_candidates", "tgt_less_than_10bp.bed")
        varInfos_l = []
        # (2) PSL --> frameshift coordinates (in protein)
        # tmp_aln_id_gene_psl_summary_file: gene name, aln_id, FS coordinates summary
        aln_id_2_psl_summary_d = self.psl_2_frameshift()
        tmp_aln_id_gene_psl_summary_file_s = os.path.join(result_dir, "frameshift", "remove.gene.aln_id.psl.txt")
        with open(tmp_aln_id_gene_psl_summary_file_s, "w") as tmp_aln_id_gene_psl_summary_file_f:
            for aln_id in aln_id_2_psl_summary_d:
                gene_name_s = aln_id_2_gene_id_d[aln_id]
                ref_size_l, tgt_size_l = [], []
                ref_starts_l = aln_id_2_psl_summary_d[aln_id]["ref"]["start"]
                ref_ends_l = aln_id_2_psl_summary_d[aln_id]["ref"]["end"]
                for start_n_end_t in zip(ref_starts_l, ref_ends_l):
                    ref_size_l.append(abs(start_n_end_t[1]-start_n_end_t[0]))
                tgt_starts_l = aln_id_2_psl_summary_d[aln_id]["tgt"]["start"]
                tgt_ends_l = aln_id_2_psl_summary_d[aln_id]["tgt"]["end"]
                for start_n_end_t in zip(tgt_starts_l, tgt_ends_l):
                    tgt_size_l.append(abs(start_n_end_t[1]-start_n_end_t[0]))
                tmp_aln_id_gene_psl_summary_file_f.write("\t".join(map(str,[gene_name_s, aln_id,
                                                                           tgt_starts_l, tgt_ends_l, ref_starts_l, ref_ends_l,
                                                                           tgt_size_l, ref_size_l]))+"\n")

        # (3) protein --> genome of FS sites from frameshift_2_genomic
        # (3-1) some trick to locate "indel" in genomic scale ==> +-1 bp
        for value_l in aln_id_2_psl_summary_d.items():
            varInfos_l.append(value_l)
        with Pool(num_processors_i) as p:
            frameshift_bedlines_nl = p.map(self.frameshift_2_genomic, varInfos_l)
        # (3-2) prot2genome_bedfile_s: bed formmated frameshift candidates from protein alignment
        ref_prot2genome_bedfile_s = os.path.join(result_dir, "frameshift", "psl_2_frameshift_candidates", "ref.raw_candidates.bed")
        tgt_prot2genome_bedfile_s = os.path.join(result_dir, "frameshift", "psl_2_frameshift_candidates", "tgt.raw_candidates.bed")
        with open(ref_prot2genome_bedfile_s, "w") as ref_prot2genome_bedfile_f, \
                open(tgt_prot2genome_bedfile_s, "w") as tgt_prot2genome_bedfile_f:
            ref_bedlines_l = []
            tgt_bedlines_l = []
            for frameshift_bedline_d in frameshift_bedlines_nl:
                ref_bedlines_l += frameshift_bedline_d["ref"]
                tgt_bedlines_l += frameshift_bedline_d["tgt"]
            ref_bedlines_l = sorted(ref_bedlines_l, key=lambda x: x.split("\t")[3])
            tgt_bedlines_l = sorted(tgt_bedlines_l, key=lambda x: x.split("\t")[3])
            ref_prot2genome_bedfile_f.write("\n".join(ref_bedlines_l))
            tgt_prot2genome_bedfile_f.write("\n".join(tgt_bedlines_l))
        # (3-3) dictionary to save frameshift id & its frameshift size for both ref and tgt
        less_than_10bp_d = {}
        with open(ref_prot2genome_bedfile_s) as ref_prot2genome_bedfile_f,\
            open(tgt_prot2genome_bedfile_s)  as tgt_prot2genome_bedfile_f,\
            open(frameshift_size_file_s, "w") as frameshift_size_file_f,\
            open(ref_less_than_10bp_bedfile_s, "w") as ref_less_than_10bp_bedfile_f,\
            open(tgt_less_than_10bp_bedfile_s, "w") as tgt_less_than_10bp_bedfile_f:
            for line_s in ref_prot2genome_bedfile_f:
                line_l = line_s.rstrip("\n").split("\t")
                frameshift_size_i = int(line_l[2]) - int(line_l[1]) - 2
                variant_id_s = line_l[3]
                assert frameshift_size_i >= 0
                if variant_id_s not in less_than_10bp_d:
                    less_than_10bp_d[variant_id_s] = {"ref":[line_s, frameshift_size_i], "tgt":[]}
            for line_s in tgt_prot2genome_bedfile_f:
                line_l = line_s.rstrip("\n").split("\t")
                frameshift_size_i = int(line_l[2]) - int(line_l[1]) - 2
                variant_id_s = line_l[3]
                assert frameshift_size_i >= 0
                less_than_10bp_d[variant_id_s]["tgt"] = [line_s, frameshift_size_i]
            frameshift_size_file_f.write("\t".join(["frameshift_id", "ref_size", "tgt_size", "less_than_10bp"])+"\n")
            print("------------------")
            print(list(less_than_10bp_d.items())[0:10])
            for indel_id_s in less_than_10bp_d:
                print(indel_id_s)
                ref_indel_size_i = less_than_10bp_d[indel_id_s]["ref"][1]
                tgt_indel_size_i = less_than_10bp_d[indel_id_s]["tgt"][1]
                # (3-4) only FS indels < 10bp on both assemblies would pass
                less_than_10bp_bool = ref_indel_size_i < 10 and tgt_indel_size_i < 10
                frameshift_size_file_f.write("\t".join(map(str, [indel_id_s, ref_indel_size_i, tgt_indel_size_i, less_than_10bp_bool])) + "\n")
                if less_than_10bp_bool is True:
                    ref_indel_bedline_s = less_than_10bp_d[indel_id_s]["ref"][0]
                    tgt_indel_bedline_s = less_than_10bp_d[indel_id_s]["tgt"][0]
                    ref_less_than_10bp_bedfile_f.write(ref_indel_bedline_s)
                    tgt_less_than_10bp_bedfile_f.write(tgt_indel_bedline_s)

    #########################################
    # 6. premature stop codon gene analysis #
    #########################################

    def variant_matrix_2_premature_stop_codon_dict(self):
        """
        from the variant matrix files ==> convert them to dict
        :return: aln_id_2_premature_stop_codon_variant_d: dict, key: alignment ID (string), value: variant bed info (nested list)
        """
        working_dir_s = os.path.join(result_dir, "variant_from_sqlDB")
        # (0) generate variant matrix (SQL DB --> TSV)
        dir_check(working_dir_s) # CDS=alignment mode of CAT
        # (1) use sql DB to extract the coordinates of premature stop codon candidates reported by CAT
        trans_mode_2_index_d = {'augTMR': 2, 'augTM': 1,'transMap': 0}  # currently: only transMap mode
        variant_matrix_header_l = ["AlignmentId", "chromosome", "start", "stop", "name", "score", "strand",
                                   "thickStart", "thickStop", "rgb", "blockCount", "blockSizes", "blockStarts"]
        gene_name_2_variant_bed_d = {}
        conn = sqlite3.connect(os.path.join(db_dir, "%s.db" % tgt_s))
        for trans_rec in tgt_all_trans_recs_l:
            cur = conn.cursor()
            parent_gene_rec = [gene_rec for gene_rec in tgtDB.parents(id=trans_rec.id, featuretype="gene")][0]
            parent_gene_name_s = get_new_name(gene_rec=parent_gene_rec, mode_s="tgt")
            aln_id_s = trans_rec.attributes["alignment_id"][0]
            trans_mode_info_nl = [[trans_mode_s, trans_mode_2_index_d[trans_mode_s]] for trans_mode_s in trans_rec.attributes["transcript_modes"]]
            # if multiple transcript prediction mode is available ==> select the first mode
            trans_mode_s = sorted(trans_mode_info_nl, key=lambda trans_mode_info_l: trans_mode_info_l[1])[-1][0]
            cur.execute('select * from CDS_%s_Evaluation where AlignmentID="%s"' % (trans_mode_s, aln_id_s))
            rows = cur.fetchall()
            for row in rows:
                rowline_s = "\t".join(str(x) for x in row)
                try:
                    gene_name_2_variant_bed_d[parent_gene_name_s].append(rowline_s)
                except:
                    gene_name_2_variant_bed_d[parent_gene_name_s] = []
                    gene_name_2_variant_bed_d[parent_gene_name_s].append(rowline_s)
        conn.close()
        variant_types_l = ["CodingInsertion", "CodingDeletion", "CodingMult3Insertion", "CodingMult3Deletion",
                           "NonCodingInsertion", "NonCodingDeletion", "InFrameStop"]
        # (2) save each variant record (AlignmentID~blockStarts) according to each variant type (CodingInsertion~InFrameStop)
        for variant_type_s in variant_types_l:
            variant_matrix_file_s = os.path.join(working_dir_s, "%s.tsv" % variant_type_s)
            variant_genes_file_s = os.path.join(working_dir_s, "%s.geneIDs.tsv" % variant_type_s)
            with open(variant_matrix_file_s, "w") as variant_matrix_file_f, \
                    open(variant_genes_file_s, "w") as variant_genes_file_f:
                variant_matrix_file_f.write("\t".join(variant_matrix_header_l) + "\n")
                for gene_name_s in gene_name_2_variant_bed_d.keys():
                    # NonCodingInsertion or NonCodingDeletion ==> exclude (not in use)
                    gene_variant_records_l = list(filter(lambda x: (variant_type_s in x) and ("Non" + x not in variant_type_s),
                                                         gene_name_2_variant_bed_d[gene_name_s]))
                    if len(gene_variant_records_l) > 0:
                        variant_matrix_file_f.write("\n".join(gene_variant_records_l) + "\n")
                        for gene_variant_record_s in gene_variant_records_l:
                            aln_id_s =  gene_variant_record_s.split("\t")[0]
                            variant_genes_file_f.write(aln_id_2_gene_id_d[aln_id_s] + "\n")

        # (3)
        aln_id_2_premature_stop_codon_variant_d = {}
        variant_matrix_file_s = os.path.join(working_dir_s, "InFrameStop.tsv")
        with open(variant_matrix_file_s) as variant_matrix_file_f:
            # First line ==> header
            variant_matrix_file_f.readline()
            for line_s in variant_matrix_file_f:
                line_l = line_s.rstrip("\n").split("\t")
                aln_id_s = line_l[0]
                variant_bed_line_l = line_l[1:]
                gene_name_s = aln_id_2_gene_id_d[aln_id_s]
                variant_bed_line_l.insert(0,gene_name_s)
                try:
                    aln_id_2_premature_stop_codon_variant_d[aln_id_s].append(variant_bed_line_l)
                except KeyError:
                    aln_id_2_premature_stop_codon_variant_d[aln_id_s] = [variant_bed_line_l]
        return aln_id_2_premature_stop_codon_variant_d

    def matching_tgt_stop_codon_to_ref_assembly(self, bed_l):
        """
        function to (1) convert variant BED string into BED file  & halLiftover it to ref assembly
        :param bed_l: list, bed-like: (from bed-like string in InFrameStop.tsv)
        [geneName,chromosome,start,stop,name,score,strand,thickStart,thickStop,rgb,blockCount,blockSizes,neighbor_size]
        :return: X, make BED file
        """
        # (1) basic variables
        gene_name_s, chrom_s, flanking_i, strand_chr = bed_l[0], bed_l[1], int(bed_l[-1]), bed_l[6]
        var_start_i, var_end_i = sorted([int(bed_l[2]), int(bed_l[3])])
        assert var_start_i + 1 < var_end_i
        start_i = var_start_i - flanking_i
        end_i = var_end_i + flanking_i
        # (2) convert variant name and generate uniqueID: "prematurestopcodon^[gene name]^[start]"
        variant_id_s = "^".join(map(str, ["prematurestopcodon", gene_name_s, start_i]))
        # (3) set the working directory and tmp files
        working_dir_s = os.path.join(result_dir, "prematurestopcodon", "flanking_%d_sequences" % flanking_i, gene_name_s)
        dir_check(working_dir_s)
        tgt_bedfile_s = os.path.join(working_dir_s, "0_%s.tgt.bed" % variant_id_s)
        tgt_fastafile_s = os.path.join(working_dir_s, "1_%s.tgt.fasta" % variant_id_s)
        ref_bedfile_s = os.path.join(working_dir_s,"2_%s.ref.bed" % variant_id_s)
        ref_flt_bedfile_s = os.path.join(working_dir_s, "3_flt_%s.ref.bed" % variant_id_s)
        ref_fastafile_s = os.path.join(working_dir_s,"4_%s.ref.fasta" % variant_id_s)
        tgt_n_ref_fasta_file_s = os.path.join(working_dir_s, "5_%s.merged.fasta" % variant_id_s)
        tgt_wo_flanking_bed_file_s = os.path.join(working_dir_s, "6_%s.original.tgt.bed" % variant_id_s)
        ref_wo_flanking_bed_file_s = os.path.join(working_dir_s, "7_%s.original.ref.bed" % variant_id_s)
        # (4) [prior --> VGP] liftover + check if the coordinates are uniquely within the same gene in both assemblies (with flanking sequences)
        # (4-1) (prior) bed file --> fasta file --> halLiftOver to VGP assembly
        with open(tgt_bedfile_s,"w") as tgt_variant_bed_file_f:
            tgt_variant_bed_file_f.write("\t".join(str(x) for x in [chrom_s, start_i, end_i, variant_id_s, 0, strand_chr]))
        os.system("%s getfasta -fi %s -bed %s -s -name -fullHeader > %s" % (bedtools_command_s, tgt_fasta_file_s, tgt_bedfile_s, tgt_fastafile_s))
        os.system("%s %s %s %s %s %s" % (halLiftover_commnd_s, hal_path_s, tgt_s, tgt_bedfile_s, ref_s, ref_bedfile_s))
        # (4-2) check uniquely liftover on the same gene
        bedlines_within_ref_gene_l = self.check_unique_lifted_bed(bed_file_s=ref_bedfile_s, var_id_s=variant_id_s)
        if len(bedlines_within_ref_gene_l) == 1:
            # (4-3) convert bed to fasta (would be used for comparing flanking sequences)
            os.system("%s getfasta -fi %s -bed %s -s -name -fullHeader > %s" % (bedtools_command_s, ref_fasta_file_s, ref_flt_bedfile_s, ref_fastafile_s))
            with open(ref_fastafile_s, "r") as ref_fastafile_f, open(tgt_fastafile_s, "r") as tgt_fastafile_f,\
                    open(tgt_n_ref_fasta_file_s, "w") as tgt_n_ref_fasta_file_f:
                tgt_n_ref_fasta_file_f.write(">TGT_%s>REF_%s" % (tgt_fastafile_f.read().lstrip(">"), ref_fastafile_f.read().lstrip(">")))

            # (5) [prior --> VGP] liftover + check if the coordinates are uniquely within the same gene in both assemblies (w/o flanking sequences)
            # (5-1) (prior) bed file --> fasta file --> halLiftOver to VGP assembly
            tgt_var_locus_s = "\t".join(str(x) for x in [chrom_s, var_start_i, var_end_i, variant_id_s, 0, strand_chr])
            with open(tgt_wo_flanking_bed_file_s, "w") as tgt_wo_flanking_bed_file_f:
                tgt_wo_flanking_bed_file_f.write(tgt_var_locus_s)
            os.system("%s %s %s %s %s %s" % (halLiftover_commnd_s, hal_path_s, tgt_s, tgt_wo_flanking_bed_file_s, ref_s, ref_wo_flanking_bed_file_s))

            # (5-2) [ref assembly, flanking X] check if liftover is uniquely done within the same gene between old and vgp assembly
            bedlines_original_ref_gene_l = self.check_unique_lifted_bed(bed_file_s=ref_wo_flanking_bed_file_s, var_id_s=variant_id_s)
            try:
                assert len(bedlines_original_ref_gene_l) == 1
                ref_var_locus_s = bedlines_original_ref_gene_l[0]
                return [tgt_var_locus_s, ref_var_locus_s]
            except AssertionError:  # when there are no actual sequences (insertion in prior or deletion in VGP)
                assert len(bedlines_original_ref_gene_l) == 0
                return [tgt_var_locus_s, variant_id_s + "_not_enough_alignment"]

    def check_unique_lifted_bed(self, bed_file_s, var_id_s):
        """
        check if a VGP coordinate which was from liftover is within the same gene of prior assembly.
        :param bed_file_s:
        :param var_id_s:
        :return:
        """
        bedlines_within_ref_gene_l = []
        gene_name_s = var_id_s.split("^")[1]
        gene_id_s = gene_name_2_gene_id_d["ref"][gene_name_s]
        gene_chrom_s, gene_start_i, gene_end_i = refDB[gene_id_s].chrom, refDB[gene_id_s].start, refDB[gene_id_s].end
        if "original" in bed_file_s:  # w/o flanking sequences
            log_file_s = bed_file_s.replace("/7_", "/7.5_").replace(".ref.bed", "log_after_lifting.ref.txt")
            tmp_bedfile_s = bed_file_s.replace("/7_", "/7.5_tmp_")
            flt_bedfile_s = bed_file_s.replace("/7_", "/8_flt_")
        else:  # with flanking sequences
            log_file_s = bed_file_s.replace("/2_", "/2.5_").replace(".ref.bed", "log_after_lifting.ref.txt")
            tmp_bedfile_s = bed_file_s.replace("/2_", "/2.5_tmp_")
            flt_bedfile_s = bed_file_s.replace("/2_", "/3_flt_")
        # (1) check if the liftover bed is within the same gene
        with open(bed_file_s) as bed_file_f, open(log_file_s, "w") as log_file_f, open(tmp_bedfile_s, "w") as tmp_bedfile_f:
            if file_check(bed_file_s) is True:
                bed_lines_l = list(filter(lambda x: x.startswith("#") is False and len(x) > 1, bed_file_f.read().split("\n")))
                bed_lines_l = sorted(bed_lines_l, key=lambda bedline_s: (bedline_s.split("\t")[0], int(bedline_s.split("\t")[1])))
                for bed_line_s in bed_lines_l:
                    bed_line_l = bed_line_s.split("\t")
                    chrom_s, start_i, end_i, variant_id_s, strand_chr = bed_line_l[0], int(bed_line_l[1]), int(bed_line_l[2]), bed_line_l[3], bed_line_l[5]
                    if (chrom_s == gene_chrom_s) and (max(end_i, gene_end_i) == gene_end_i) and (min(start_i+1, gene_start_i) == gene_start_i):
                        tmp_bedfile_f.write(bed_line_s+"\n")
                    else:
                        log_file_f.write("<<< FILTERED OUT...\n %s BECAUSE NO ACTUAL ALIGNMENTS WERE IN HAL FILE >>>\n" % bed_line_s)
            else:
                log_file_f.write("NO BED INFORMATION: %s" % bed_file_s)
                return []
        # (2) merge liftover bed blocks if they are within 5 bp
        os.system("%s merge -s -d 5 -c 4,5,6 -o distinct,distinct,distinct -i %s > %s" % (bedtools_command_s, tmp_bedfile_s, flt_bedfile_s))
        with open(flt_bedfile_s) as flt_bedfile_f:
            for bedline_s in flt_bedfile_f:
                bedlines_within_ref_gene_l.append(bedline_s.rstrip("\n"))
        return bedlines_within_ref_gene_l

    def mp_premature_stop_codon(self, num_processors_i, aln_id_2_premature_stop_codon_variant_d, flanking_i=5):
        """
        function to (1) make tgt BED & (2) halliftover tgtBED to ref assembly & (3) compare the sequences
        :param num_processors_i: int, number of multiple process
        :param aln_id_2_premature_stop_codon_variant_d: dict, key: aln_id (string),  value: BED-like list (list)
        :param flanking_i: int, size of flanking sequences
        :return: X
        """
        working_dir_s = os.path.join(result_dir, "prematurestopcodon")
        dir_check(os.path.join(working_dir_s, "flanking_%d_sequences" % flanking_i))
        # (1) liftover PMS of old assembly to VGP assembly
        variant_beds_nl = []
        for value_l in aln_id_2_premature_stop_codon_variant_d.values():
            for variant_bed_l in value_l:
                variant_bed_l.append(flanking_i)
                variant_beds_nl.append(variant_bed_l)
        with Pool(num_processors_i) as p:
            tgt_n_ref_locus_nl = p.map(self.matching_tgt_stop_codon_to_ref_assembly, variant_beds_nl)
            print("length of locus of premature stop codon (before filtering no ref record cases) : %d" % len(tgt_n_ref_locus_nl))
            tgt_n_ref_locus_nl = [x for x in tgt_n_ref_locus_nl if x is not None]
            tgt_n_ref_locus_nl = [x for x in tgt_n_ref_locus_nl if x[-1].endswith("_not_enough_alignment") is False]
            print("length of locus of premature stop codon (after filtering no ref record cases) : %d" % len(tgt_n_ref_locus_nl))
            ref_var_loci_bedfile_s = os.path.join(working_dir_s, "ref_loci_with_unique_lifted_loci.bed")
            tgt_var_loci_bedfile_s = os.path.join(working_dir_s, "tgt_loci_with_unique_lifted_loci.bed")
            with open(ref_var_loci_bedfile_s,"w") as ref_var_loci_bedfile_f, open(tgt_var_loci_bedfile_s,"w") as tgt_var_loci_bedfile_f:
                for tgt_var_locus_s, ref_var_locus_s in tgt_n_ref_locus_nl:
                    tgt_var_loci_bedfile_f.write(tgt_var_locus_s + "\n")
                    ref_var_loci_bedfile_f.write(ref_var_locus_s + "\n")
            p.close()
            p.join()
        # (2) compare flanking sequences between vgp and old assembly
        neigbor_compare_file_s = os.path.join(working_dir_s, "genome_compare.tsv")
        tmp_compare_file_s = os.path.join(working_dir_s, "summary.genome_compare.tmp")
        ref_candidates_file_s = os.path.join(working_dir_s, "ref_premature_stopcodon_candidates.bed")
        tgt_candidates_file_s = os.path.join(working_dir_s, "tgt_premature_stopcodon_candidates.bed")
        ref_all_bedfile_s = os.path.join(working_dir_s, "ref_loci_with_unique_lifted_loci.bed")
        tgt_all_bedfile_s = os.path.join(working_dir_s, "tgt_loci_with_unique_lifted_loci.bed")
        with open(neigbor_compare_file_s,"w") as neigbor_compare_file_f, \
                open(ref_all_bedfile_s) as ref_all_bedfile_f, open(tgt_all_bedfile_s) as tgt_all_bedfile_f,\
                open(ref_candidates_file_s,"w") as ref_candidates_file_f, open(tgt_candidates_file_s,"w") as tgt_candidates_file_f:
            # (3-1) read unique bed lines and convert it to dict
            bed_id_2_bedline_d = {"ref": {}, "tgt": {}}
            for mode_s, file_f in [["ref",ref_all_bedfile_f], ["tgt", tgt_all_bedfile_f]]:
                for line_s in file_f:
                    line_l = line_s.rstrip("\n").split("\t")
                    if len(line_l) == 6:
                        bed_id_s = line_l[3]
                        bed_id_2_bedline_d[mode_s][bed_id_s] = line_s
            # (3-2) for each PMS record, check flanking fasta
            for variant_bed_l in variant_beds_nl:
                gene_name_s, chrom_s, flanking_i, variant_s, strand_chr = variant_bed_l[0], variant_bed_l[1], variant_bed_l[-1],variant_bed_l[4], variant_bed_l[6]
                start_i = min(int(variant_bed_l[2]), int(variant_bed_l[3])) - flanking_i
                assert variant_s == "InFrameStop"
                variant_id_s = "^".join(map(str, ["prematurestopcodon", gene_name_s, start_i]))
                size_i = variant_bed_l[9]
                tmp_working_dir_s = os.path.join(working_dir_s, "flanking_%d_sequences" % flanking_i, gene_name_s)
                tgt_fastafile_s = os.path.join(tmp_working_dir_s, "1_%s.tgt.fasta" % variant_id_s)
                ref_fastafile_s = os.path.join(tmp_working_dir_s, "4_%s.ref.fasta" % variant_id_s)
                tgt_n_ref_fasta_file_s = os.path.join(tmp_working_dir_s, "5_%s.merged.fasta" % variant_id_s)
                if file_check(tgt_fastafile_s) and file_check(ref_fastafile_s) and file_check(tgt_n_ref_fasta_file_s) is True:
                    # (3-3) if all liftover is successfully performed, generate dict with sequences and headers
                    with open(tgt_n_ref_fasta_file_s) as tgt_n_ref_fasta_file_f:
                        mode_2_seq_with_flanking_d = {}
                        for values in FastaIO.SimpleFastaParser(tgt_n_ref_fasta_file_f):
                            header_s = values[0]
                            locus_info_s = header_s.split("::")[-1]
                            if header_s.startswith("REF_"):
                                ref_locus_info_s = locus_info_s
                                mode_2_seq_with_flanking_d[ref_s] = values[1]
                            else:
                                assert header_s.startswith("TGT_")
                                tgt_locus_info_s = locus_info_s
                                mode_2_seq_with_flanking_d[tgt_s] = values[1]
                        # (3-4) if both old and vgp sites have actual sequences (longer thant 2*flanking_i), compare flanking sequences
                        ref_genomic_seq_s, tgt_genomic_seq_s = mode_2_seq_with_flanking_d[ref_s].upper(), mode_2_seq_with_flanking_d[tgt_s].upper()
                        if min(len(ref_genomic_seq_s),len(tgt_genomic_seq_s)) >= 2*flanking_i:
                            ref_left_seq_s, ref_mid_seq_s, ref_right_seq_s = cut_neighbor_sequences(seq_s=ref_genomic_seq_s, flanking_i=flanking_i)
                            tgt_left_seq_s, tgt_mid_seq_s, tgt_right_seq_s = cut_neighbor_sequences(seq_s=tgt_genomic_seq_s, flanking_i=flanking_i)
                            ident_left_i, ident_right_i = identity(ref_left_seq_s, tgt_left_seq_s), identity(ref_right_seq_s, tgt_right_seq_s)
                            return_s = "\t".join(map(str, [variant_id_s, ref_genomic_seq_s, tgt_genomic_seq_s,
                                                           ref_left_seq_s, tgt_left_seq_s, ident_left_i,
                                                           ref_mid_seq_s, tgt_mid_seq_s, abs(len(ref_mid_seq_s)-len(tgt_mid_seq_s)),
                                                           ref_right_seq_s, tgt_right_seq_s, ident_right_i,
                                                           str(size_i), tgt_locus_info_s, ref_locus_info_s]))+"\n"
                            neigbor_compare_file_f.write(return_s)
                            if ident_left_i + ident_right_i == 10 and ref_mid_seq_s != tgt_mid_seq_s:
                                ref_candidates_file_f.write(bed_id_2_bedline_d["ref"].get(variant_id_s,""))
                                tgt_candidates_file_f.write(bed_id_2_bedline_d["tgt"].get(variant_id_s,""))
        os.system("cut -f1,2,4,5,6,9,10,11,12,13 %s > %s" % (neigbor_compare_file_s, tmp_compare_file_s))
    ############################################
    # 7. splicing junction dirsuption analysis #
    ############################################
    def splicing_junction(self, num_processors_i):
        """
        function to detect splicing junction disruption
        :param num_processors_i: int, number of processors
        :return:
        """
        working_dir_s = os.path.join(result_dir, "intronexonjunctiondisruption")
        dir_check(working_dir_s)
        ref_disrupted_junction_candidate_bedfile_s = os.path.join(working_dir_s, "ref.splicing_junction_disruption_candidates.bed")
        tgt_disrupted_junction_candidate_bedfile_s = os.path.join(working_dir_s, "tgt.splicing_junction_disruption_candidates.bed")
        # (1) In VGP assembly, collect every introns, check if its length is > 30bp (biological intron), save as a bed file
        ref_gene_with_intron_l = ref_gene_name_2_intron_d.keys()
        ref_donor_and_acceptor_bed_l = self.convert_intron_2_bed_file(genes_with_intron_l=ref_gene_with_intron_l, mode_s="ref", minimum_intron_i=30)
        ref_intron_junction_bed_file_s = os.path.join(working_dir_s, "0.ref_junctions_over_30bp_introns.bed")
        with open(ref_intron_junction_bed_file_s, "w") as ref_intron_junction_bed_file_f:
            ref_intron_junction_bed_file_f.write("\n".join(ref_donor_and_acceptor_bed_l))
        # (2) In prior assembly, collect every introns, check if its length is > 30bp (biological intron), and save as a bed file
        tgt_gene_with_intron_l = []
        for gene_rec in tgt_all_gene_recs_l:
            gene_name_s = get_new_name(gene_rec=gene_rec, mode_s="tgt")
            intron_recs_l = [intron_rec for intron_rec in tgtDB.children(id=gene_rec.id, featuretype="intron")]
            if len(intron_recs_l) > 0 and gene_name_s in gene_name_2_gene_id_d["ref"]:
                tgt_gene_with_intron_l.append(gene_name_s)
        tgt_donor_and_acceptor_bed_l = self.convert_intron_2_bed_file(genes_with_intron_l=tgt_gene_with_intron_l, mode_s="tgt", minimum_intron_i=30)
        tgt_intron_junction_bed_file_s = os.path.join(working_dir_s, "2.tgt_junctions_over_30bp_introns.bed")
        with open(tgt_intron_junction_bed_file_s, "w") as tgt_intron_junction_bed_file_f:
            tgt_intron_junction_bed_file_f.write("\n".join(tgt_donor_and_acceptor_bed_l))
        # (3) VGP splicing junctions --> liftover to prior assembly
        #     check if the projected intron junctions are also annotated as splicing junction in the projected annotation by CAT
        tgt_projected_bed_file_s = os.path.join(working_dir_s, "1.tgt_projected.bed")
        tgt_projected_junction_tmp_file_s = os.path.join(working_dir_s, "3.tgt_junctions_with_projected.tmp")
        tgt_projected_junction_bed_file_s = os.path.join(working_dir_s, "3.tgt_junctions_with_projected.bed")
        os.system("%s %s %s %s %s %s" % (halLiftover_commnd_s, hal_path_s, ref_s, ref_intron_junction_bed_file_s, tgt_s, tgt_projected_bed_file_s))
        # (3-1) intersect projected junctions to annotated junction with 100% coverage
        os.system("%s intersect -f 1.00 -wo -a %s -b %s > %s"
                  % (bedtools_command_s, tgt_projected_bed_file_s, tgt_intron_junction_bed_file_s, tgt_projected_junction_tmp_file_s))
        with open(tgt_projected_junction_tmp_file_s, "r") as tgt_projected_junction_tmp_file_f, \
                open(tgt_projected_junction_bed_file_s, "w") as tgt_projected_junction_bed_file_f:
            projected_nl, junctions_nl = [], []
            projected_d = {}  # key: coord of tgt projected junctions, value: dict with "count (how many intersections?)", "index (Where is it?)"
            index_i = 0
            for line_s in tgt_projected_junction_tmp_file_f:
                line_l = line_s.rstrip("\n").split("\t")
                projected_l, junction_l = line_l[:6], line_l[6:]
                projected_s = "\t".join(projected_l)
                projected_nl.append(projected_s)
                if projected_s not in projected_d:
                    projected_d[projected_s] = {"count": 0, "index": []}
                projected_d[projected_s]["count"] += 1
                projected_d[projected_s]["index"].append(index_i)
                junctions_nl.append("\t".join(junction_l))
                index_i += 1
            input_nl = [[projected_s, projected_d, junctions_nl] for projected_s in projected_nl]
            with Pool(num_processors_i) as p:
                out = p.map(self.check_whether_projected_and_tgt_junctions_are_identical, input_nl)
                out_not_sorted_beds_l = [return_value for return_value in out if return_value is not None]
                out_not_sorted_beds_l = list(set(out_not_sorted_beds_l)) # projected to the same gene
                out_bed_l = sorted(out_not_sorted_beds_l, key=lambda bedline_s:[bedline_s.split("\t")[0],int(bedline_s.split("\t")[1])])
            tgt_projected_junction_bed_file_f.write("\n".join(out_bed_l))

        # (4) check if junctions are canonical
        # (4-1) ref and projected junction: convert bed to fasta
        ref_intron_junction_fasta_file_s = os.path.join(working_dir_s, "0.ref_junctions_over_30bp_introns.fasta")
        tgt_projected_junction_fasta_file_s = os.path.join(working_dir_s, "3.tgt_junctions_with_projected.fasta")
        os.system("%s getfasta -fi %s -bed %s -s -name -fullHeader > %s" %
                  (bedtools_command_s, ref_fasta_file_s, ref_intron_junction_bed_file_s, ref_intron_junction_fasta_file_s))
        os.system("%s getfasta -fi %s -bed %s -s -name -fullHeader > %s" %
                  (bedtools_command_s, tgt_fasta_file_s, tgt_projected_junction_bed_file_s, tgt_projected_junction_fasta_file_s))
        # (4-2) save junction sequences into nested dictionary (gene_2_junction_d)
        gene_2_junction_d = {} # nested dict, key1 = gene_name, key2 = intron_id, key3 = "donor" or "acceptor", key4: "ref" or "tgt", value = junction seq
        with open(ref_intron_junction_fasta_file_s) as ref_intron_junction_fasta_file_f,\
                open(tgt_projected_junction_fasta_file_s) as tgt_projected_junction_fasta_file_f:
            for values in FastaIO.SimpleFastaParser(ref_intron_junction_fasta_file_f):
                header_s, seq_s = values[0], values[1].upper()
                assert header_s.count("^") == 2
                header_l = header_s.split("^")
                gene_name_s, intron_id_s, junction_type_s = header_l[1], header_l[2].split("=")[0], header_l[2].split("=")[-1].split("(")[0]
                if gene_name_s not in gene_2_junction_d:
                    gene_2_junction_d[gene_name_s] = {}
                if intron_id_s not in gene_2_junction_d[gene_name_s]:
                    gene_2_junction_d[gene_name_s][intron_id_s] = {"donor": {"ref": "", "tgt": ""}, "acceptor": {"ref": "", "tgt": ""}}
                assert junction_type_s in ["donor", "acceptor"]
                gene_2_junction_d[gene_name_s][intron_id_s][junction_type_s]["ref"] = seq_s
            for values in FastaIO.SimpleFastaParser(tgt_projected_junction_fasta_file_f):
                header_s, seq_s = values[0], values[1].upper()
                assert header_s.count("^") == 2
                header_l = header_s.split("^")
                gene_name_s, intron_id_s, junction_type_s = header_l[1], header_l[2].split("=")[0], header_l[2].split("=")[-1].split("(")[0]
                assert gene_name_s in gene_2_junction_d
                assert intron_id_s in gene_2_junction_d[gene_name_s]
                assert junction_type_s in ["donor", "acceptor"]
                gene_2_junction_d[gene_name_s][intron_id_s][junction_type_s]["tgt"] = seq_s
        # (4-3) compare VGP(ref) and prior(tgt) junction sequences --> collect intron id with splicing junction disruptions
        canonical_splicing_junctions_l = ["GTAG", "ATAC", "GCAG"]
        candidate_junction_ids_l = []
        summary_lines_l = []
        for gene_name_s in gene_2_junction_d:
            for intron_id_s in gene_2_junction_d[gene_name_s]:
                type_s = "default"
                ref_donor_seq_s = gene_2_junction_d[gene_name_s][intron_id_s]["donor"]["ref"]
                ref_acceptor_seq_s = gene_2_junction_d[gene_name_s][intron_id_s]["acceptor"]["ref"]
                tgt_donor_seq_s = gene_2_junction_d[gene_name_s][intron_id_s]["donor"]["tgt"]
                tgt_acceptor_seq_s = gene_2_junction_d[gene_name_s][intron_id_s]["acceptor"]["tgt"]
                is_ref_canonical_bool = (ref_donor_seq_s + ref_acceptor_seq_s) in canonical_splicing_junctions_l
                is_tgt_canonical_bool = (tgt_donor_seq_s + tgt_acceptor_seq_s) in canonical_splicing_junctions_l
                is_donor_same_bool = ref_donor_seq_s == tgt_donor_seq_s
                is_acceptor_same_bool = ref_acceptor_seq_s == tgt_acceptor_seq_s
                if (is_ref_canonical_bool is True) & (is_donor_same_bool is False) \
                        & (len(ref_donor_seq_s) == 2) & (len(tgt_donor_seq_s) == 2) & ("N" not in tgt_donor_seq_s):
                    candidate_junction_ids_l.append("SplicingJunctionDisruption^%s^%s=donor" % (gene_name_s, intron_id_s))
                    type_s = "error"
                if (is_ref_canonical_bool is True) & (is_acceptor_same_bool is False) &\
                        (len(ref_acceptor_seq_s) == 2) & (len(tgt_acceptor_seq_s) == 2) & ("N" not in tgt_acceptor_seq_s):
                    candidate_junction_ids_l.append("SplicingJunctionDisruption^%s^%s=acceptor" % (gene_name_s, intron_id_s))
                    type_s = "error"
                summary_line_s = "\t".join(map(str,[gene_name_s, intron_id_s, ref_donor_seq_s, ref_acceptor_seq_s, tgt_donor_seq_s, tgt_acceptor_seq_s,
                                                    is_ref_canonical_bool, is_tgt_canonical_bool, is_donor_same_bool, is_acceptor_same_bool, type_s]))
                summary_lines_l.append(summary_line_s)
        junction_disruption_summary_file_s = os.path.join(working_dir_s, "summary_on_splicing_junction_disruption.summary.tabs")
        with open(junction_disruption_summary_file_s, "w") as junction_disruption_summary_file_f:
            junction_disruption_summary_file_f.write("\t".join(["gene_id", "intron_id", "ref_donor", "ref_acceptor", "tgt_donor", "tgt_acceptor",
                                                                "is_ref_canonical", "is_tgt_canonical", "is_donor_same", "is_acceptor_same", "type"])+"\n")
            junction_disruption_summary_file_f.write("\n".join(summary_lines_l))

        # (5) collect splicing junction disruptions --> save as a bed file
        with open(ref_disrupted_junction_candidate_bedfile_s, "w") as ref_disrupted_junction_candidate_bedfile_f,\
                open(tgt_disrupted_junction_candidate_bedfile_s, "w") as tgt_disrupted_junction_candidate_bedfile_f,\
                open(ref_intron_junction_bed_file_s, "r") as ref_intron_junction_bed_file_f,\
                open(tgt_projected_junction_bed_file_s, "r") as tgt_projected_bed_file_f:
            ref_candidates_i = 0
            tgt_candidates_i = 0
            ref_gene_and_type_ids_l = []
            ref_bedlines_l = []
            tgt_index_and_lines_nl = []
            for line_s in ref_intron_junction_bed_file_f:
                line_l = line_s.rstrip("\n").split('\t')
                junction_id_s = line_l[3]
                if junction_id_s in candidate_junction_ids_l:
                    ref_gene_name_s, ref_type_s = junction_id_s.split("^")[1], junction_id_s.split("^")[-1]
                    ref_gene_and_type_id_s = "^".join([ref_gene_name_s, ref_type_s])
                    ref_gene_and_type_ids_l.append(ref_gene_and_type_id_s)
                    ref_bedlines_l.append(line_s)
                    ref_candidates_i += 1
            ref_disrupted_junction_candidate_bedfile_f.write("".join(ref_bedlines_l))
            for line_s in tgt_projected_bed_file_f:
                line_l = line_s.rstrip("\n").split('\t')
                junction_id_s = line_l[3]
                if junction_id_s in candidate_junction_ids_l:
                    tgt_gene_name_s, tgt_type_s = junction_id_s.split("^")[1], junction_id_s.split("^")[-1]
                    tgt_gene_and_type_id_s = "^".join([tgt_gene_name_s, tgt_type_s])
                    assert tgt_gene_and_type_id_s in ref_gene_and_type_ids_l
                    tgt_index_i = ref_gene_and_type_ids_l.index(tgt_gene_and_type_id_s)
                    tgt_index_and_lines_nl.append([tgt_index_i, line_s])
                    tgt_candidates_i += 1
            sorted_tgt_index_and_lines_l = sorted([tgt_index_and_line for tgt_index_and_line in tgt_index_and_lines_nl],
                                                  key=lambda tgt_index_and_line:tgt_index_and_line[0])
            sorted_tgt_lines_l = [tgt_index_and_line[1] for tgt_index_and_line in sorted_tgt_index_and_lines_l]
            tgt_disrupted_junction_candidate_bedfile_f.write("".join(sorted_tgt_lines_l))
            assert ref_candidates_i == tgt_candidates_i

    def convert_intron_2_bed_file(self, genes_with_intron_l, mode_s, minimum_intron_i=30):
        """
        function to convert intron records into bed file
        (1) collect every introns, (2) check if its length is > 30bp (biological intron), and (3) liftover to prior assembly
        :param gene_ids_with_intron_l: list, gene id with one or more intron records
        :param mode_s: string, "ref" for VGP, "tgt" for Prior
        :return: donor_and_acceptor_bed_l: list, with bed-style string
        """
        assert mode_s in ["ref", "tgt"]
        db = refDB if mode_s == "ref" else tgtDB
        donor_and_acceptor_bed_l = []
        for gene_name_s in genes_with_intron_l:
            gene_id_s = gene_name_2_gene_id_d[mode_s][gene_name_s]
            gene_chrom_s = db[gene_id_s].chrom
            gene_strand_chr = db[gene_id_s].strand
            if gene_strand_chr == "+":
                reverse_bool = False
            else:
                assert gene_strand_chr == "-"
                reverse_bool = True
            sorted_intron_recs_l = sorted([intron_rec for intron_rec in db.children(id=gene_id_s, featuretype="intron")],
                                          key=lambda intron_rec: intron_rec.start, reverse=reverse_bool)
            # (1-2) check if its length is > 30bp (biological intron) & make ID for each donor and acceptor
            intron_order_i = 0
            for intron_rec in sorted_intron_recs_l:
                intron_order_i += 1
                intron_length_i = intron_rec.end-intron_rec.start+1
                if intron_length_i < minimum_intron_i: # non-biological intron
                    continue
                assert intron_length_i >= minimum_intron_i
                if gene_strand_chr == "+":
                    donor_start_i, donor_end_i = intron_rec.start-1, (intron_rec.start-1)+2
                    acceptor_start_i, acceptor_end_i = intron_rec.end-2, intron_rec.end
                else:
                    assert gene_strand_chr == "-"
                    donor_start_i, donor_end_i = intron_rec.end-2, intron_rec.end
                    acceptor_start_i, acceptor_end_i = intron_rec.start-1, (intron_rec.start-1)+2
                donor_intron_id_s = "SplicingJunctionDisruption^%s^intron_%d=donor" % (gene_name_s, intron_order_i)
                acceptor_intron_id_s = "SplicingJunctionDisruption^%s^intron_%d=acceptor" % (gene_name_s, intron_order_i)
                donor_bed_s = "\t".join(map(str, [gene_chrom_s, donor_start_i, donor_end_i, donor_intron_id_s, 0, gene_strand_chr]))
                acceptor_bed_s = "\t".join(map(str, [gene_chrom_s, acceptor_start_i, acceptor_end_i, acceptor_intron_id_s, 0, gene_strand_chr]))
                donor_and_acceptor_bed_l.append(donor_bed_s)
                donor_and_acceptor_bed_l.append(acceptor_bed_s)
        return donor_and_acceptor_bed_l

    def check_whether_projected_and_tgt_junctions_are_identical(self, input_l):
        """
        function to check if projected loci & junctions are perfectly overlapping (Same gene, same type (donor or acceptor?))
        :param input_l: nested list, [projected_s, projected_d, junctions_nl]
        :return:
        """
        projected_s, projected_d, junctions_nl = input_l
        projected_l = projected_s.split("\t")
        projected_gene_name_s, projected_type_s = projected_l[3].split("^")[1], projected_l[3].split("=")[-1].split("(")[0]
        projected_gene_and_type_id_s = "^".join([projected_gene_name_s, projected_type_s])
        all_indices_l = projected_d[projected_s]["index"]
        junction_gene_and_type_ids_l = []
        for index_i in all_indices_l:
            junction_s = junctions_nl[index_i]
            junction_l = junction_s.split("\t")
            junction_gene_name_s, junction_type_s = junction_l[3].split("^")[1], junction_l[3].split("=")[-1].split("(")[0]
            junction_gene_and_type_id_s = "^".join([junction_gene_name_s, junction_type_s])
            junction_gene_and_type_ids_l.append(junction_gene_and_type_id_s)
        if junction_gene_and_type_ids_l.count(projected_gene_and_type_id_s) == 1:
            return projected_s
        else:
            return None

    def generate_intron_records_from_ref_assembly(self):
        """
        function to generate intron information of reference assembly
        :return: gene_name_2_intron_d: {geneID_s:[rec1,rec2,rec3],...}
                 intron_recs_l: nested list, [[rec1, rec1.chrom, rec1.start, rec1.end, rec1.strand],...]
        """
        working_dir_s = os.path.join(result_dir, "intronexonjunctiondisruption")
        dir_check(working_dir_s)
        gene_name_2_intron_d = {}
        intron_recs_l = []
        valid_gene_id_l = [gene_rec.id for gene_rec in ref_all_gene_recs_l]
        ref_gene_name_2_intron_pickle_s = os.path.join(working_dir_s, "ref_gene_name_2_intron_d.pickle")
        ref_intron_rec_pickle_s = os.path.join(working_dir_s, "intron_recs_l.pickle")
        if os.path.isfile(ref_gene_name_2_intron_pickle_s) is False:
            with open(ref_gene_name_2_intron_pickle_s, 'wb') as ref_gene_name_2_intron_pickle_f,\
                    open(ref_intron_rec_pickle_s, 'wb') as ref_intron_rec_pickle_f:
                # generate intron from ref gff by using create_introns
                try:
                    introns = [intron_rec for intron_rec in refDB.all_features(featuretype="intron")]
                    assert len(introns) > 0
                except AssertionError:
                    introns = list(refDB.create_introns(exon_featuretype="exon", grandparent_featuretype="gene"))
                    refDB.update(introns, merge_strategy="create_unique")
                for ref_intron_rec in introns:
                    parent_gene_rec = [gene_rec for gene_rec in refDB.parents(ref_intron_rec.attributes["Parent"][0], featuretype="gene")][0]
                    gene_id_s = parent_gene_rec.id
                    if gene_id_s not in valid_gene_id_l:
                        continue
                    else:
                        gene_name_s = gene_id_2_gene_name_d["ref"][gene_id_s]
                        try:
                            gene_name_2_intron_d[gene_name_s].append(ref_intron_rec)
                        except:
                            gene_name_2_intron_d[gene_name_s] = [ref_intron_rec]
                        intron_recs_l.append(ref_intron_rec)
                # save intron dict and recoreds to pickle file
                pickle.dump(gene_name_2_intron_d, ref_gene_name_2_intron_pickle_f, pickle.HIGHEST_PROTOCOL)
                pickle.dump(intron_recs_l, ref_intron_rec_pickle_f, pickle.HIGHEST_PROTOCOL)
        else:
            with open(ref_gene_name_2_intron_pickle_s, 'rb') as ref_gene_name_2_intron_pickle_f,\
                    open(ref_intron_rec_pickle_s, 'rb') as ref_intron_rec_pickle_f:
                gene_name_2_intron_d = pickle.load(ref_gene_name_2_intron_pickle_f)
                intron_recs_l = pickle.load(ref_intron_rec_pickle_f)
        return gene_name_2_intron_d, intron_recs_l

    ##################################
    # 8. N in coding region analysis #
    ##################################
    def main_n_in_coding_region(self):
        """
        function to find the projected genes with N
        :return:
        """
        working_dir_s = os.path.join(result_dir, "nincodingregion")
        dir_check(working_dir_s)
        n_in_coding_region_summary_file_s = os.path.join(working_dir_s, "nincodingregion.summary.tabs")

        # 1. prepare necessary files ==> CDS and blat result (nt-wise)
        ref_all_cdsfile_s = os.path.join(result_dir, "CDS", "ref", "ref_consensus.cds.fasta")
        tgt_all_cdsfile_s = os.path.join(result_dir, "CDS", "tgt", "tgt_consensus.cds.fasta")
        if file_check(ref_all_cdsfile_s) is False or file_check(tgt_all_cdsfile_s) is False:
            self.mp_gff2cds(num_processors_i=process_num_i)

        # 2. get CDS sequence and make them dictionary (header 2 seq)
        ref_header_2_cds_d = {}
        tgt_header_2_cds_d = {}
        with open(ref_all_cdsfile_s, "r") as ref_all_cdsfile_f, open(tgt_all_cdsfile_s, "r") as tgt_all_cdsfile_f:
            for values in FastaIO.SimpleFastaParser(ref_all_cdsfile_f):
                ref_header_2_cds_d[values[0]] = values[1]
            for values in FastaIO.SimpleFastaParser(tgt_all_cdsfile_f):
                tgt_header_2_cds_d[values[0]] = values[1]

        # 3. check if N in coding regions
        with open(n_in_coding_region_summary_file_s,"w") as n_in_coding_region_summary_file_f:
            info_nl = [[tgt_header_2_cds_d[key].upper(),ref_header_2_cds_d[key].upper(),key.split("_CDS")[0]] for key in tgt_header_2_cds_d]
            info_nl = [x for x in list(filter(lambda x: ("N" in x[0]) & ("N" not in x[1]), info_nl))]
            n_in_coding_region_genes_l = list(set([x[2] for x in info_nl]))
            n_in_coding_region_genes_l = sorted(list(set(n_in_coding_region_genes_l)))
            n_in_coding_region_summary_file_f.write("\n".join(n_in_coding_region_genes_l))

    ######################
    # 9. summary results #
    ######################
    def read_mapping_result_summary(self):
        """
        function to summary the number of "vgp_error","previous_error","no_errors","filter_out" based on read mapping result
        :return: X, save the output file with the number of each error type of each species
        """
        outfile_s = os.path.join(result_dir, "supplementary.count_variants.tabs")
        species_l = ["Zebra finch", "Anna's hummingbird", "Platypus", "Climbing perch"]
        with open(outfile_s, "w") as outfile_f:
            outfile_f.write("\t".join(["FGL_type", "error_type", "species", "count"])+"\n")
            for fgl_s in ["frameshift", "prematurestopcodon", "intronexonjunctiondisruption"]:
                # read each species' mpileup result for each type of false gene loss
                finch_summary_file_s = os.path.join(home_dir_s, "zebra_finch", "tm", "CAT2MISSING_pre_taegut", fgl_s, "locus_wise.mpileup_summary.tabs")
                hummingbird_summary_file_s = os.path.join(home_dir_s, "anna", "tm", "CAT2MISSING_pre_calann", fgl_s, "locus_wise.mpileup_summary.tabs")
                platypus_summary_file_s = os.path.join(home_dir_s, "platypus", "tm", "CAT2MISSING_pre_ornana", fgl_s, "locus_wise.mpileup_summary.tabs")
                perch_summary_file_s = os.path.join(home_dir_s, "climbing_perch", "tm", "CAT2MISSING_pre_anates", fgl_s, "locus_wise.mpileup_summary.tabs")
                # convert them to dictionary
                finch_error_d = self.mpileup_summary_to_dict(finch_summary_file_s)
                hummingbird_error_d = self.mpileup_summary_to_dict(hummingbird_summary_file_s)
                platypus_error_d = self.mpileup_summary_to_dict(platypus_summary_file_s)
                perch_error_d = self.mpileup_summary_to_dict(perch_summary_file_s)
                for error_type_s in ["vgp_error","previous_error","no_errors","filter_out"]:
                    for species_s, species_d in zip(species_l, [finch_error_d, hummingbird_error_d, platypus_error_d, perch_error_d]):
                        outfile_f.write("\t".join([fgl_s, error_type_s, species_s, str(species_d[error_type_s])]) + "\n")

    def mpileup_summary_to_dict(self, file_s):
        """
        function to summary the number of errors ("vgp_error","previous_error","no_errors","filter_out") to dict
        :param file_s: string, path to locus wise mpileup summary file
        :return: dict, key: error type (string), value: the number of the error (integer)
        """
        with open(file_s) as file_f:
            first_line_s = file_f.readline()
        error_type_2_count_d = {x.split(":")[0]:int(x.split(":")[1]) for x in first_line_s.lstrip("#").split(",")}
        return error_type_2_count_d

    def ratio_of_fgl(self):
        """
        function to calculate ratio of false gene loss (structure, sequence, structure & sequence level)
        :return: X
        """
        fgl_genes_d = self.count_fgl_number()
        structure_fgl_types_l = ["totallymissing","exonDeletion","fragmented","intrascaffoldsplit"]
        sequence_fgl_types_l = ["frameshift","prematurestopcodon","intronexonjunctiondisruption","nincodingregion"]
        fgl_types_l = structure_fgl_types_l + sequence_fgl_types_l
        fgl_gene_names_l, structure_gene_name_l, sequence_gene_name_l = [], [], []
        for fgl_s in fgl_types_l:
            fgl_gene_names_l += fgl_genes_d[fgl_s]
            if fgl_s in structure_fgl_types_l:
                structure_gene_name_l += fgl_genes_d[fgl_s]
            if fgl_s in sequence_fgl_types_l:
                sequence_gene_name_l += fgl_genes_d[fgl_s]
        fgl_gene_names_l = sorted(list(set(fgl_gene_names_l)))
        structure_gene_name_l = sorted(list(set(structure_gene_name_l)))
        sequence_gene_name_l = sorted(list(set(sequence_gene_name_l)))
        both_gene_name_l = sorted(list(set(structure_gene_name_l) & set(sequence_gene_name_l)))
        structure_only_l = sorted(list(set(structure_gene_name_l) - set(both_gene_name_l)))
        sequence_only_l = sorted(list(set(sequence_gene_name_l) - set(both_gene_name_l)))
        ref_all_genes_i = len(ref_all_gene_recs_l)
        normal_i = ref_all_genes_i - len(fgl_gene_names_l)

        fgl_fl = float(len(fgl_gene_names_l)/ref_all_genes_i)*100
        structure_fl = float(len(structure_only_l)/ref_all_genes_i)*100
        sequence_fl = float(len(sequence_only_l)/ref_all_genes_i)*100
        both_fl = float(len(both_gene_name_l)/ref_all_genes_i)*100
        normal_fl = float(normal_i/ref_all_genes_i)*100
        percent_summary_file_s = os.path.join(result_dir, "total.summary.tsv")
        with open(percent_summary_file_s, "w") as percent_summary_file_f:
            percent_summary_file_f.write("\t".join(["type", "count", "percent", "ratio_of_fgl", "sum_of_fgl"])+"\n")
            percent_summary_file_f.write("normal\t%d\t%f\t%f\t%d\n" % (normal_i, normal_fl, fgl_fl, len(fgl_gene_names_l)))
            percent_summary_file_f.write("structure\t%d\t%f\t%f\t%d\n" % (len(structure_only_l), structure_fl, fgl_fl, len(fgl_gene_names_l)))
            percent_summary_file_f.write("sequence\t%d\t%f\t%f\t%d\n" % (len(sequence_only_l), sequence_fl, fgl_fl, len(fgl_gene_names_l)))
            percent_summary_file_f.write("structure_and_sequence\t%d\t%f\t%f\t%d\n" % (len(both_gene_name_l), both_fl, fgl_fl, len(fgl_gene_names_l)))

    def count_fgl_number(self):
        """
        function to count and summary the false gene loss found in the prior assembly
        :return:
        """
        # summary dictionary
        fgl_2_genes_d = {"totallymissing": [],"exonDeletion": [],"fragmented" : [],"intrascaffoldsplit": [],
                         "frameshift": [], "prematurestopcodon": [],"intronexonjunctiondisruption": [],"nincodingregion": []}
        # (1) totally missing and exon deletion
        totally_missing_summary_file_s = os.path.join(result_dir, "missingGenes", "totally_missing_and_exon_deletion.tabs")
        totally_missing_genes_l = []
        exon_deletion_genes_l = []
        with open(totally_missing_summary_file_s) as totally_missing_summary_file_f:
            totally_missing_summary_file_f.readline()
            for line_s in totally_missing_summary_file_f:
                line_l = line_s.rstrip("\n").split("\t")
                gene_id_s, type_s = line_l[0], line_l[-1]
                if type_s == "totally_missing":
                    totally_missing_genes_l.append(gene_id_s)
                if type_s == "exon_deletion":
                    exon_deletion_genes_l.append(gene_id_s)
        fgl_2_genes_d["totallymissing"] = sorted(list(set(totally_missing_genes_l)))
        fgl_2_genes_d["exonDeletion"] = sorted(list(set(exon_deletion_genes_l)))

        # (2) fragmented & intrascaffold_split
        fragmented_summary_file_s = os.path.join(result_dir, "fragmented", "fragmented_genes.txt")
        with open(fragmented_summary_file_s) as fragmented_summary_file_f:
            fragmented_genes_l = [x.rstrip("\n").split("\t")[0] for x in fragmented_summary_file_f.readlines() if len(x)>= 1]
        fgl_2_genes_d["fragmented"] = sorted(list(set(fragmented_genes_l)))
        intra_split_summary_file_s = os.path.join(result_dir, "intrascaffoldsplit", "intrascaffoldsplit_genes.txt")
        with open(intra_split_summary_file_s) as intra_split_summary_file_f:
            intra_split_genes_l = [x.rstrip("\n").split("\t")[0] for x in intra_split_summary_file_f.readlines() if len(x)>= 1]
        fgl_2_genes_d["intrascaffoldsplit"] = sorted(list(set(intra_split_genes_l)))

        # (3) count N in Coding Region
        n_in_coding_region_summary_file_s = os.path.join(result_dir, "nincodingregion", "nincodingregion.summary.tabs")
        n_in_coding_region_genes = []
        with open(n_in_coding_region_summary_file_s) as n_in_coding_region_summary_file_f:
            for line_s in n_in_coding_region_summary_file_f:
                gene_id_s = line_s.rstrip("\n")
                n_in_coding_region_genes.append(gene_id_s)
        fgl_2_genes_d["nincodingregion"] = sorted(list(set(n_in_coding_region_genes)))

        # (4) frameshift, premature stop codon, splicing junction disruption
        for fgl_s in ["frameshift", "prematurestopcodon", "intronexonjunctiondisruption"]:
            fgl_genes_l = []
            mpileup_file_s = os.path.join(result_dir, fgl_s, "ref_filtered_by_mpileup.bed")
            with open(mpileup_file_s) as mpileup_file_f:
                for line_s in mpileup_file_f:
                    if line_s.startswith("#"):
                        continue
                    line_l = line_s.rstrip("\n").split("\t")
                    gene_id_s = line_l[3].split("^")[1]
                    fgl_genes_l.append(gene_id_s)
                fgl_genes_l = sorted(list(set(fgl_genes_l)))
                fgl_2_genes_d[fgl_s] = fgl_genes_l
        # summary the number and list of genes with errors
        for fgl_s in fgl_2_genes_d:
            for gene_name_s in fgl_2_genes_d[fgl_s]:
                if gene_name_s not in gene_name_2_gene_id_d["ref"]:
                    fgl_2_genes_d[fgl_s].remove(gene_name_s)
        count_file_s, genes_file_s = os.path.join(result_dir,"total.count.tabs"), os.path.join(result_dir, "total.genes.tabs")
        with open(count_file_s, "w") as count_file_f, open(genes_file_s, "w") as genes_file_f:
            for fgl_s in fgl_2_genes_d.keys():
                count_file_f.write("%s\t%d\n" % (fgl_s, len(fgl_2_genes_d[fgl_s])))
                genes_file_f.write("%s\t%s\n" % (fgl_s, ",".join(fgl_2_genes_d[fgl_s])))
        return fgl_2_genes_d

    def mp_gc_n_repeat_of_genes(self,num_processors_i):
        """
        function to calculate GC content of genes (in VGP assembly)
        :param num_processors_i: integer, number of processors
        :return: X
        """
        working_dir_s = os.path.join(result_dir, "basic_inputs")
        cds_dir_s = os.path.join(result_dir, "CDS")
        gene_gc_n_repeat_s = os.path.join(working_dir_s, "gc_and_repeat_of_genes.stats")
        dir_check(working_dir_s)
        fgl_2_genes_d = self.count_fgl_number()
        # (1) Get GC & repeat content (gene-level)
        seq_and_gene_name_nl = []
        gene_bed_file_s = os.path.join(working_dir_s, "ref.genes.bed")
        gene_fasta_file_s = os.path.join(working_dir_s, "ref.genes.fasta")
        if file_check(gene_fasta_file_s) is False:
            os.system("%s getfasta -fi %s -bed %s -name -s > %s" % (bedtools_command_s, ref_fasta_file_s, gene_bed_file_s, gene_fasta_file_s))
        with open(gene_fasta_file_s, "r") as gene_fasta_file_f:
            for values in FastaIO.SimpleFastaParser(gene_fasta_file_f):
                gene_name_s = values[0].split("(")[0]
                seq_s = values[1]
                seq_and_gene_name_nl.append([seq_s, gene_name_s])
        # (2) summary it to dict (gene-level)
        with Pool(num_processors_i) as p:
            gene_2_gc_n_repeat_d = {}
            gene_2_infos_l = p.map(gc_and_repeat_content_of_genes, seq_and_gene_name_nl)
            for gene_2_info_d in gene_2_infos_l:
                gene_name_s = gene_2_info_d["gene"]
                gene_2_gc_n_repeat_d[gene_name_s] = {"gc": gene_2_info_d["gc"], "repeat": gene_2_info_d["repeat"]}
        # (3) save
        with open(gene_gc_n_repeat_s, "w") as gene_gc_n_repeat_f:
            gene_gc_n_repeat_f.write("\t".join(["gene", "GCcontents", "repeat_percent", "type"])+"\n")
            missing_gene_names_l = fgl_2_genes_d["totallymissing"]
            all_gene_names_l = sorted(list(gene_id_2_gene_name_d["ref"][gene_rec.id] for gene_rec in ref_all_gene_recs_l))
            for gene_name_s in all_gene_names_l:
                type_s = "fgl" if gene_name_s in missing_gene_names_l else "control"
                gene_gc_n_repeat_f.write("\t".join(map(str, [gene_name_s, gene_2_gc_n_repeat_d[gene_name_s]["gc"],
                                                             gene_2_gc_n_repeat_d[gene_name_s]["repeat"], type_s]))+"\n")


    def genome_gc_n_repeat_summary(self,window_size_i=10000):
        """
        calculate GC and repeat content of VGP and prior assemblies (with a given window size)
        :param window_size_i: integer, size of windows (non-overlapping)
        :return: X
        """
        working_dir_s = os.path.join(result_dir, "window_size_gc_and_repeat")
        dir_check(working_dir_s)
        # Aligned vs. missing regions
        consensus_missing_bedfile_s = os.path.join(result_dir, "missing_chrom", "3_consensus", "consensus.missing.bed")
        consensus_aligned_bedfile_s = os.path.join(result_dir, "missing_chrom", "3_consensus", "consensus.aligned.bed")
        consensus_missing_fastafile_s = os.path.join(result_dir, "missing_chrom", "3_consensus", "consensus.missing.fasta")
        consensus_aligned_fastafile_s = os.path.join(result_dir, "missing_chrom", "3_consensus", "consensus.aligned.fasta")
        # manual exclusion of chr29 of zebra finch
        if ref_s == "vgp_taegut":
            chr_29_dup_bedfile_s = "PATH TO VGP02All_merged.merged_by_bedtools.bed"
            consensus_missing_rm_dup_bedfile_s = os.path.join(result_dir, "missing_chrom", "3_consensus", "consensus.missing.rm_chr29_dup.bed")
            consensus_aligned_rm_dup_bedfile_s = os.path.join(result_dir, "missing_chrom", "3_consensus", "consensus.aligned.rm_chr29_dup.bed")
            os.system("%s subtract -a %s -b %s > %s" % (bedtools_command_s, consensus_missing_bedfile_s, chr_29_dup_bedfile_s, consensus_missing_rm_dup_bedfile_s))
            os.system("%s subtract -a %s -b %s > %s" % (bedtools_command_s, consensus_aligned_bedfile_s, chr_29_dup_bedfile_s, consensus_aligned_rm_dup_bedfile_s))
            consensus_missing_bedfile_s = consensus_missing_rm_dup_bedfile_s
            consensus_aligned_bedfile_s = consensus_aligned_rm_dup_bedfile_s
        os.system("%s getfasta -fi %s -bed %s > %s" % (bedtools_command_s, ref_fasta_file_s, consensus_missing_bedfile_s, consensus_missing_fastafile_s))
        os.system("%s getfasta -fi %s -bed %s > %s" % (bedtools_command_s, ref_fasta_file_s, consensus_aligned_bedfile_s, consensus_aligned_fastafile_s))
        concat_missing_fastafile_s = os.path.join(result_dir, "missing_chrom", "3_consensus", "concat.missing.fasta")
        concat_aligned_fastafile_s = os.path.join(result_dir, "missing_chrom", "3_consensus", "concat.aligned.fasta")
        concatenate_fasta(input_path_s=consensus_missing_fastafile_s, output_path_s=concat_missing_fastafile_s, header_s="missing")
        concatenate_fasta(input_path_s=consensus_aligned_fastafile_s, output_path_s=concat_aligned_fastafile_s, header_s="aligned")
        # reference (VGP) and target (Prior)
        for mode_s, fasta_file_s in zip(["ref", "tgt", "missing", "aligned"],
                                        [ref_fasta_file_s, tgt_fasta_file_s, concat_missing_fastafile_s, concat_aligned_fastafile_s]):
            genome_bed_file_s = os.path.join(result_dir, "basic_inputs", "%s.genome.bed" % mode_s)
            softmasked_bed_file_s = os.path.join(result_dir, "basic_inputs", "%s_softmasked.bed" % mode_s)
            if file_check(softmasked_bed_file_s) is False:
                get_softmasked(work_dir=result_dir, fasta_file_s=fasta_file_s, mode_s=mode_s)
            gc_repeat_window_file_s = os.path.join(result_dir, "window_size_gc_and_repeat", "%s.summary.tsv" % mode_s)
            window_level_gc_n_repeat(window_size_i=window_size_i, genome_bed_file_s=genome_bed_file_s, softmasked_bed_file_s=softmasked_bed_file_s,
                                     mode_s=mode_s, fasta_file_s=fasta_file_s, gc_repeat_window_file_s=gc_repeat_window_file_s)

    def circos_maker(self, scaff_2_chrom_d, window_size_i=10000):
        """
        function to generate input data for circos plot
        :param scaff_2_chrom_d: dict, scaffold id 2 chromosome
        :param window_size_i: integer, window size to calculate GC, repeat content and missing ratio
        :return:
        """
        out_dir_circos_s = os.path.join(result_dir, "chromosomeMap","out_circos")
        dir_check(out_dir_circos_s)
        # (1) chromosome to scaffold id
        chrom_2_scaff_d = {}
        for scaff_id_s in scaff_2_chrom_d:
            chrom_id_s = scaff_2_chrom_d[scaff_id_s]
            chrom_2_scaff_d[chrom_id_s] = scaff_id_s
        # (2) generate dict to store scaffold information (chromosome/unlocalized/unplaced scaffold)
        scaff_2_info_d = {}
        ref_scaffold_ids_l = sorted(list(ref_chrom_2_size_d.keys()))
        for scaffold_id_s in ref_scaffold_ids_l:
            size_i = ref_chrom_2_size_d[scaffold_id_s]
            scaffold_type_s = check_scaffold_type(scaffold_id_s=scaffold_id_s)
            scaff_2_info_d[scaffold_id_s] = {"all_genes": [], "missing_genes": [], "scaffold_gc": 0, "scaffold_repeat": 0,
                                             "size": size_i, "type": scaffold_type_s, "missing_ratio": 0, "aligned_ratio": 0, "gap_ratio": 0,
                                             "gap_count": 0, "gap_length": 0,"partial_genes": []}
        # (3) Circos input 1: chromosome size file
        scaffold_size_file_s = os.path.join(out_dir_circos_s, "chromosomes.size.csv")
        with open(scaffold_size_file_s, "w") as scaffold_size_file_f:
            scaffold_size_file_f.write(",".join(["chr","start","end","value1","value2"])+"\n")
            for scaffold_id_s in ref_chrom_2_size_d:
                size_i = ref_chrom_2_size_d[scaffold_id_s]
                scaffold_size_file_f.write(",".join(map(str,[scaffold_id_s, 1, size_i,"NA","NA"]))+"\n")
        # (4) Circos input 2:missing, repeat, and gc
        ref_genome_bed_file_s = os.path.join(result_dir, "basic_inputs", "ref.genome.bed")
        missing_bed_file_s = os.path.join(result_dir, "missing_chrom", "3_consensus", "consensus.missing.bed")
        windowmasked_bed_file_s = os.path.join(result_dir, "basic_inputs","ref_softmasked.bed")
        scaffold_nuc_file_s = os.path.join(out_dir_circos_s, "chromosomes.nuc")
        scaffold_gc_file_s = os.path.join(out_dir_circos_s, "chromosomes.gc.csv")
        scaffold_repeat_file_s = os.path.join(out_dir_circos_s, "chromosomes.repeat.csv")
        scaffold_aln_file_s = os.path.join(out_dir_circos_s, "chromosomes.aln.csv")
        os.system("%s makewindows -b %s -w %d |" % (bedtools_command_s, ref_genome_bed_file_s, window_size_i) +
                  "%s intersect -wao -a stdin -b %s |" % (bedtools_command_s, windowmasked_bed_file_s) +
                  "%s sort -i stdin |" % bedtools_command_s +
                  "%s groupby -i stdin -g 1,2,3,4 -c 7 -o sum |" % bedtools_command_s +
                  "awk -F'[\"\t\"]' '{OFS=\"\t\"} {print $1,$2,$3,\"window|\"$5/($3-$2),0,\".\"}'|" +
                  "%s intersect -wao -a stdin -b %s |" % (bedtools_command_s, missing_bed_file_s) +
                  "%s sort -i stdin |" % bedtools_command_s +
                  "%s groupby -i stdin -g 1,2,3,4 -c 13 -o sum |" % bedtools_command_s +
                  "awk -F'[\"\t\"]' '{OFS=\"\t\"} {print $1,$2,$3,$4\"|\"$5/($3-$2),0,\".\"}'|" +
                  "%s nuc -fi %s -bed stdin > %s" % (bedtools_command_s, ref_fasta_file_s, scaffold_nuc_file_s))
        with open(scaffold_nuc_file_s, "r") as scaffold_nuc_file_f, \
                open(scaffold_repeat_file_s, "w") as scaffold_repeat_file_f, \
                open(scaffold_gc_file_s, "w") as scaffold_gc_file_f,\
                open(scaffold_aln_file_s, "w") as scaffold_aln_file_f:
            scaffold_nuc_file_f.readline()
            header_s = ",".join(["chr","start","value"])+"\n"
            scaffold_gc_file_f.write(header_s)
            scaffold_repeat_file_f.write(header_s)
            scaffold_aln_file_f.write(header_s)
            for line_s in scaffold_nuc_file_f:
                line_l = line_s.rstrip("\n").split('\t')
                chrom_s, start_i, id_s = line_l[0], round((int(line_l[1]) + int(line_l[2]))/2, 0), line_l[3]
                repeat_fl, missing_fl = id_s.split("|")[1:]
                numA, numC, numG, numT = line_l[8], line_l[9], line_l[10], line_l[11]
                try:
                    gc_fl = (int(numC) + int(numG)) / sum([int(x) for x in [numA, numC, numG, numT]])
                except ZeroDivisionError:
                    continue
                scaffold_repeat_file_f.write(",".join(map(str,[chrom_s, start_i, repeat_fl]))+"\n")
                scaffold_gc_file_f.write(",".join(map(str,[chrom_s, start_i, gc_fl]))+"\n")
                scaffold_aln_file_f.write(",".join(map(str,[chrom_s, start_i, 1-float(missing_fl)]))+"\n")
        # (5) Circos input 3: gap
        gap_bedfile_s = os.path.join(out_dir_circos_s, "vgp_gap_plot.csv")
        ref_gap_info_bedfile_s = os.path.join(result_dir, "basic_inputs", "ref.gaps.bed")
        with open(ref_gap_info_bedfile_s, "r") as ref_gap_info_bedfile_f, \
                open(gap_bedfile_s, "w") as gap_bedfile_f:
            gap_bedfile_f.write(",".join(["chr", "start", "end"]) + "\n")
            for line_s in ref_gap_info_bedfile_f:
                line_l = line_s.rstrip("\n").split("\t")
                new_line_s = ",".join(line_l[:3]) + "\n"
                gap_bedfile_f.write(new_line_s)
        # (6) Circos input 4: gene density
        gene_density_csv_file_s = os.path.join(out_dir_circos_s, "vgp_genes_plot.csv")
        gene_density_bed_file_s = os.path.join(out_dir_circos_s, "gene_density.200kb.bed")
        ref_gene_bed_file_s = os.path.join(result_dir, "basic_inputs", "ref.genes.bed")
        os.system("%s makewindows -b %s -w 200000 | %s intersect -a stdin -b %s -c > %s" %
                  (bedtools_command_s, ref_genome_bed_file_s,
                   bedtools_command_s, ref_gene_bed_file_s, gene_density_bed_file_s))
        with open(gene_density_bed_file_s, "r") as gene_density_bed_file_f, \
                open(gene_density_csv_file_s, "w") as gene_density_csv_file_f:
            gene_density_csv_file_f.write(",".join(["chr", "start", "end"]) + "\n")
            for line_s in gene_density_bed_file_f:
                line_l = line_s.rstrip("\n").split("\t")
                scaffold_id_s, start_i, end_i, gene_density_s = line_l[0], int(line_l[1]), int(line_l[2]), line_l[3]
                new_line_s = ",".join([scaffold_id_s, str(round((start_i + end_i) / 2)), gene_density_s]) + "\n"
                gene_density_csv_file_f.write(new_line_s)
        # (7) Circos input 5: missing genes
        fgl_2_gene_d = self.count_fgl_number()
        missing_gene_names_l = fgl_2_gene_d["totallymissing"]
        scaffold_missing_file_s = os.path.join(out_dir_circos_s, "chromosomes.missing_loci.csv")
        scaffold_label_file_s = os.path.join(out_dir_circos_s, "chromosomes.missing_label.csv")
        with open(scaffold_missing_file_s, "w") as scaffold_missing_file_f:
            scaffold_missing_file_f.write(",".join(["chr", "start", "end", "label"]) + "\n")
            for gene_rec in ref_all_gene_recs_l:
                scaffold_id_s = gene_rec.chrom
                start_i, end_i = gene_rec.start, gene_rec.end
                gene_name_s = gene_id_2_gene_name_d["ref"][gene_rec.id]
                scaff_2_info_d[scaffold_id_s]["all_genes"].append(gene_name_s)
                if gene_name_s in missing_gene_names_l:
                    scaffold_missing_file_f.write(",".join(map(str,[scaffold_id_s, start_i, end_i, gene_name_s])) + "\n")
                    scaff_2_info_d[scaffold_id_s]["missing_genes"].append(gene_name_s)
        with open(scaffold_label_file_s, "w") as scaffold_label_file_f:
            scaffold_label_file_f.write(",".join(["chr", "start", "label"]) + "\n")
            for scaffold_id_s in ref_chrom_2_size_d:
                missing_gene_num_i = len(scaff_2_info_d[scaffold_id_s]["missing_genes"])
                all_gene_num_i = len(scaff_2_info_d[scaffold_id_s]["all_genes"])
                size_i = ref_chrom_2_size_d[scaffold_id_s]
                missing_gene_ratio_s = str(round(missing_gene_num_i/all_gene_num_i*100,1))+"%" if all_gene_num_i > 0 else "NA"
                scaffold_label_file_f.write(",".join(map(str, [scaffold_id_s, round(size_i/2,0), missing_gene_ratio_s]))+"\n")
        # (8) Circos input 6: chromosome name
        chromosome_name_file_s = os.path.join(out_dir_circos_s, "chromosomes.name.csv")
        with open(chromosome_name_file_s, "w") as chromosome_name_file_f:
            header_s = ",".join(["chr","locus","value"]) + "\n"
            chromosome_name_file_f.write(header_s)
            scaffold_name_lines_nl = list()
            for scaffold_id_s in ref_scaffold_ids_l:
                half_size_i = round(scaff_2_info_d[scaffold_id_s]["size"]/2)
                try:
                    scaffold_name_s = scaff_2_chrom_d[scaffold_id_s].replace("chr.","").replace(",u",":u")
                except KeyError:
                    scaffold_name_s = scaffold_id_s
                name_line_l = ",".join(map(str,[scaffold_id_s, half_size_i, scaffold_name_s]))
                scaffold_name_lines_nl.append(name_line_l)
            chromosome_name_file_f.write("\n".join(scaffold_name_lines_nl))

        # (9) summary statistics of each scaffold
        out_dir_correl_s = os.path.join(result_dir, "chromosomeMap", "out_correlation")
        dir_check(out_dir_correl_s)
        # (9-1) partial gene flags
        for gene_rec in ref_all_gene_recs_l:
            if gene_rec.attributes["gene_biotype"][0] == "protein_coding":
                scaffold_id_s = gene_rec.chrom.split(".")[0]
                if scaffold_id_s in ref_scaffold_ids_l:
                    if gene_rec.attributes.get("partial", "NA")[0] == "true":
                        scaff_2_info_d[scaffold_id_s]["partial_genes"].append(gene_rec.id)

        # (9-2) GC and repeat for each scaffold
        ref_genome_nuc_file_s = os.path.join(out_dir_correl_s, "ref.genome.nuc")
        os.system("%s nuc -fi %s -bed %s > %s" % (bedtools_command_s, ref_fasta_file_s, ref_genome_bed_file_s, ref_genome_nuc_file_s))
        scaff_2_gc_or_repeat_d = {}
        for scaffold_id_s in ref_chrom_2_size_d:
            scaff_2_gc_or_repeat_d[scaffold_id_s] = {"GC": 0, "Repeat": 0}
        with open(windowmasked_bed_file_s) as windowmasked_bed_file_f:
            for line_s in windowmasked_bed_file_f:
                line_l = line_s.rstrip("\n").split('\t')
                chrom_s, start_i, end_i = line_l[0], int(line_l[1]), int(line_l[2])
                size_i = end_i - start_i
                scaff_2_gc_or_repeat_d[chrom_s]["Repeat"] += size_i
        with open(ref_genome_nuc_file_s) as ref_genome_nuc_file_f:
            ref_genome_nuc_file_f.readline()
            for line_s in ref_genome_nuc_file_f:
                line_l = line_s.rstrip("\n").split('\t')
                chrom_s, gc_fl = line_l[0], (int(line_l[9]) + int(line_l[10])) / (int(line_l[14]) - int(line_l[12]))
                scaff_2_gc_or_repeat_d[chrom_s]["GC"] += gc_fl
        for scaffold_id_s in ref_chrom_2_size_d:
            repeat_fl = scaff_2_gc_or_repeat_d[scaffold_id_s]["Repeat"]
            gc_fl = scaff_2_gc_or_repeat_d[scaffold_id_s]["GC"]
            repeat_content_fl = repeat_fl / scaff_2_info_d[scaffold_id_s]["size"] * 100
            gc_content_fl = gc_fl * 100
            scaff_2_info_d[scaffold_id_s]["scaffold_gc"] = gc_content_fl
            scaff_2_info_d[scaffold_id_s]["scaffold_repeat"] = repeat_content_fl

        # (9-3) Count of gap
        ref_gap_info_bedfile_s = os.path.join(out_dir_correl_s, "ref_genomic_gaps.bed")
        ref_genbank_2_refseq_d = convert_header_of_assmebly_statistics(assembly_report_file_s=ref_assembly_report_file_s)
        convert_genomic_gaps_to_bed(gap_info_file_s=ref_gap_info_file_s, gap_info_bedfile_s=ref_gap_info_bedfile_s, genbank_2_refseq_d=ref_genbank_2_refseq_d)
        with open(ref_gap_info_bedfile_s, "r") as ref_gap_info_bedfile_f:
            for line_s in ref_gap_info_bedfile_f:
                line_l = line_s.rstrip("\n").split("\t")
                scaffold_id_s = line_l[0]
                gap_size_i = abs(int(line_l[2]) - int(line_l[1]))
                scaff_2_info_d[scaffold_id_s]["gap_count"] += 1
                scaff_2_info_d[scaffold_id_s]["gap_length"] += gap_size_i
        missing_ratio_file_s = os.path.join(result_dir, "missing_chrom", "3_consensus", "all_scaffolds.missing_ratio.tsv")
        with open(missing_ratio_file_s) as missing_ratio_file_f:
            missing_ratio_file_f.readline()
            for line_s in missing_ratio_file_f:
                line_l = line_s.rstrip("\n").split("\t")
                scaffold_id_s = chrom_2_scaff_d.get(line_l[0],line_l[0])
                missing_ratio_fl = float(line_l[5])
                aligned_ratio_fl = float(line_l[6])
                gap_ratio_fl = float(line_l[7])
                scaff_2_info_d[scaffold_id_s]["missing_ratio"] = missing_ratio_fl
                scaff_2_info_d[scaffold_id_s]["aligned_ratio"] = aligned_ratio_fl
                scaff_2_info_d[scaffold_id_s]["gap_ratio"] = gap_ratio_fl

        # (9-4) final summary file
        summary_chromosomes_file_s = os.path.join(out_dir_correl_s, "%s.stat_summary.tabs" % ref_s)
        with open(summary_chromosomes_file_s, "w") as summary_chromosomes_file_f:
            headers_l = ["scaffold", "type", "size", "scaffold_gc", "scaffold_repeat",
                         "num_of_genes", "num_of_missing_genes", "missing_gene_pct", "num_of_partial_genes",
                         "partial_gene_pct", "missing_ratio", "aligned_ratio", "gap_ratio", "gap_count_freq", "gap_length_freq"]
            summary_chromosomes_file_f.write("\t".join(headers_l)+"\n")
            for scaffold_id_s in ref_scaffold_ids_l:
                chrom_id_s = scaff_2_chrom_d.get(scaffold_id_s, scaffold_id_s)
                type_s = scaff_2_info_d[scaffold_id_s]["type"]
                size_i = scaff_2_info_d[scaffold_id_s]["size"]
                scaffold_gc_fl = scaff_2_info_d[scaffold_id_s]["scaffold_gc"]
                scaffold_repeat_fl = scaff_2_info_d[scaffold_id_s]["scaffold_repeat"]
                num_of_genes_i = len(list(set(scaff_2_info_d[scaffold_id_s]["all_genes"])))
                num_of_missing_genes_i = len(list(set(scaff_2_info_d[scaffold_id_s]["missing_genes"])))
                try:
                    missing_gene_pct_fl = round(num_of_missing_genes_i/num_of_genes_i*100,5)
                except ZeroDivisionError:
                    missing_gene_pct_fl = "NA"
                num_of_partial_genes = len(list(set(scaff_2_info_d[scaffold_id_s]["partial_genes"])))
                try:
                    partial_gene_pct_fl = round(num_of_partial_genes/num_of_genes_i*100,5)
                except ZeroDivisionError:
                    partial_gene_pct_fl = "NA"
                missing_ratio_fl = scaff_2_info_d[scaffold_id_s]["missing_ratio"]
                aligned_ratio_fl = scaff_2_info_d[scaffold_id_s]["aligned_ratio"]
                gap_ratio_fl = scaff_2_info_d[scaffold_id_s]["gap_ratio"]
                gap_count_freq_fl = scaff_2_info_d[scaffold_id_s]["gap_count"]/size_i
                gap_length_freq_fl = scaff_2_info_d[scaffold_id_s]["gap_length"] / size_i
                summary_chromosomes_file_f.write("\t".join(map(str,[chrom_id_s, type_s, size_i, scaffold_gc_fl, scaffold_repeat_fl,
                                                                    num_of_genes_i, num_of_missing_genes_i, missing_gene_pct_fl,
                                                                    num_of_partial_genes, partial_gene_pct_fl, missing_ratio_fl, aligned_ratio_fl, gap_ratio_fl, gap_count_freq_fl, gap_length_freq_fl]))+"\n")
    ##############################
    # 10. gene structure analysis #
    ##############################
    def add_UTR_manually(self, DB, mode_s):
        """
        function to manually add UTR features on the GFF file by using GffUtils DB object
        :param DB: gffutils DB to add manually made UTR records
        :param mode_s: string, type of genome annotation to add UTR records: "ref" (=vgp), "tgt" (=prev, projected by CAT), or "prev" (=prev. annotated by NCBI)
        :return:
        """
        # Only for protein coding genes
        gene_recs_l = [gene_rec for gene_rec in DB.all_features(featuretype="gene") if gene_rec.attributes["gene_biotype"][0] == "protein_coding"]
        gff_lines_l = []
        # For each protein-coding gene, check and annotate UTR if possible
        for gene_rec in gene_recs_l:
            # (1) utr index: to distinguish multiple UTRs in one gene
            utr_index_i = 0
            # (2) basic gene information
            gene_id_s = gene_rec.id
            gene_strand = gene_rec.strand
            gene_chrom_s = gene_rec.chrom
            # (3) if multiple transcripts for one gene:
            # ref or prev: NCBI annotation -->
            for mRNA_rec in DB.children(id=gene_id_s, featuretype="mRNA"):
                mRNA_id_s = mRNA_rec.id
                gene_dbxref_id_s = gene_rec["Dbxref"][0] if mode_s in ("ref", "prev") else mRNA_rec["dbxref"][0]
                if mode_s == "ref":
                    gene_id_s = gene_dbxref_id_s.split(":")[1]
                elif mode_s == "tgt":
                    gene_id_s = gene_rec.attributes["source_gene"][0]
                else:
                    assert mode_s == "prev"
                    gene_id_s = gene_dbxref_id_s.split(" ")[0].split(":")[-1]
                # logic: (1) sort exon and cds based on their start coordinates
                # logic: (2) if no CDS recored: exit
                if gene_strand == "+":
                    left_UTR_s, right_UTR_s = "5", "3"
                else:
                    assert gene_strand == "-"
                    left_UTR_s, right_UTR_s = "3", "5"
                utr_type_2_name_d = {"3": "three_prime_UTR", "5": "five_prime_UTR"}
                cds_recs_l = sorted([cds_rec for cds_rec in DB.children(id=mRNA_id_s, featuretype="CDS")],
                                    key=lambda cds_rec: cds_rec.start)
                exon_recs_l = sorted([exon_rec for exon_rec in DB.children(id=mRNA_id_s, featuretype="exon")],
                                     key=lambda exon_rec: exon_rec.start)
                if len(cds_recs_l) == 0:
                    print(str(gene_rec))
                    continue
                # get the leftmost/rightmost cds and exon
                leftmost_cds_rec = cds_recs_l[0]
                rightmost_cds_rec = cds_recs_l[-1]
                # logic: (3) compare each exon records' start and end with those of leftmost/rightmost CDSs
                leftmost_cds_start_i, leftmost_cds_end_i = leftmost_cds_rec.start, leftmost_cds_rec.end
                rightmost_cds_start_i, rightmost_cds_end_i = rightmost_cds_rec.start, rightmost_cds_rec.end
                for exon_rec in exon_recs_l:
                    left_type1_bool, left_type2_bool, right_type1_bool, right_type2_bool = False, False, False, False
                    exon_start_i, exon_end_i = exon_rec.start, exon_rec.end
                    # (1) leftmost UTR
                    if exon_start_i < leftmost_cds_start_i:
                        left_type1_bool = (leftmost_cds_start_i <= exon_end_i)
                        left_type2_bool = (exon_end_i < leftmost_cds_start_i)
                    # (2) rightmost UTR
                    if exon_end_i > rightmost_cds_end_i:
                        right_type1_bool = (exon_start_i <= rightmost_cds_end_i)
                        right_type2_bool = (rightmost_cds_end_i < exon_start_i)
                        ###########################################################################################
                        #  [Classification of leftmost_UTRs]                                                      #
                        #  For each exon record, if the "start" of exon (^) < the start of 1st coding exon (!):   #
                        #   left_type1: first_cds_start_i(!) <= exon_end_i(*) (^ --> ! --> *)                     #
                        #                   (Left) >>>>[     ]-----[      ]-----[      ]---...                    #
                        #                     EXON:^_________*      ______       ______                           #
                        #                      CDS:    !_____?      ______       ______                           #
                        #                      UTR:_____                                                          #
                        #                         (>>>>[     ]>>>> is possible as well)                           #
                        #                                                                                         #
                        #   left_type2: exon_end_i (*) < first_cds_start_i(!) (^ --> * --> !)                     #
                        #                   (Left) >>>>>>>>>>>-----[      ]-----[      ]---...                    #
                        #                     EXON:^_________*      ______       ______                           #
                        #                      CDS:                 !_____       ______                           #
                        #                      UTR:___________                                                    #
                        #                                                                                         #
                        #  [Classification of rightmost_UTRs]                                                     #
                        #  For each exon record, if the start of 1st coding exon (?) < the "end" of exon (*):     #
                        #  right_type1: exon_start_i (^) <= last_cds_end_i (?)                                    #
                        #                    ...---[      ]-----[      ]---[       ]>>>>> (Right)                 #
                        #                     EXON: ______       ______     ^___________*                         #
                        #                      CDS: ______       ______     !_____?                               #
                        #                      UTR:                                ______                         #
                        #                         (>>>>[     ]>>>> is possible: ? <= * would be proper)           #
                        #                                                                                         #
                        #  right_type2: last_cds_end_i(?) < exon_start_i(^)                                       #
                        #                    ...---[      ]-----[      ]---->>>>>>>>>>>>> (Right)                 #
                        #                     EXON: ______       ______     ^___________*                         #
                        #                      CDS: ______       !____?                                           #
                        #                      UTR:                         _____________                         #
                        #  Note that each exon would be compared:                                                 #
                        #   both type1 and type2 UTRs can exist within one gene                                   #
                        ###########################################################################################
                    if (left_type1_bool | left_type2_bool) is True:
                        if left_type1_bool is True:
                            left_utr_type_s = "utr_%s_1" % left_UTR_s
                            left_utr_start_i = exon_start_i
                            left_utr_end_i = leftmost_cds_start_i - 1
                        else:
                            assert left_type2_bool is True
                            left_utr_type_s = "utr_%s_2" % left_UTR_s
                            left_utr_start_i = exon_start_i
                            left_utr_end_i = exon_end_i
                        utr_index_i += 1
                        info_line_s = ";".join(["ID=%s_%s_%d" % (left_utr_type_s, gene_id_s, utr_index_i),
                                                "Parent=%s" % (mRNA_id_s), "transcript_id=%s" % (mRNA_id_s),
                                                "Dbxref=%s,Genbank:%s" % (gene_dbxref_id_s, mRNA_id_s)])
                        gffline_s = "\t".join(map(str, [gene_chrom_s, "Gnomon", utr_type_2_name_d[left_UTR_s],
                                                        left_utr_start_i, left_utr_end_i, ".", gene_strand, ".",
                                                        info_line_s]))
                        gff_lines_l.append(gffline_s)

                    if (right_type1_bool | right_type2_bool) is True:
                        if right_type1_bool is True:
                            right_utr_type_s = "utr_%s_1" % right_UTR_s
                            right_utr_start_i = rightmost_cds_end_i + 1
                            right_utr_end_i = exon_end_i
                        else:
                            assert right_type2_bool is True
                            right_utr_type_s = "utr_%s_2" % right_UTR_s
                            right_utr_start_i = exon_start_i
                            right_utr_end_i = exon_end_i
                        utr_index_i += 1
                        info_line_s = ";".join(["ID=%s_%s_%d" % (right_utr_type_s, gene_id_s, utr_index_i),
                                                "Parent=%s" % (mRNA_id_s), "transcript_id=%s" % (mRNA_id_s),
                                                "Dbxref=%s,Genbank:%s" % (gene_dbxref_id_s, mRNA_id_s)])
                        gffline_s = "\t".join(map(str, [gene_chrom_s, "Gnomon", utr_type_2_name_d[right_UTR_s],
                                                        right_utr_start_i, right_utr_end_i, ".", gene_strand, ".",
                                                        info_line_s]))
                        gff_lines_l.append(gffline_s)
        return gff_lines_l

    def get_representative(self, DB, input_gff_s, out_gff_s, dup_gene_ids_l):
        """
        fuction to get the longest transcript
        :param DB:
        :param working_dir:
        :param mode_s:
        :return:
        """
        with open(out_gff_s, "w") as out_gff_file_f, open(input_gff_s) as original_gff_file_f:
            for line_s in original_gff_file_f:
                if line_s.startswith("#"):
                    out_gff_file_f.write(line_s)
            # 1. gene level
            for gene_rec in DB.all_features(featuretype="gene"):
                # exclude duplicated genes (start) #
                if gene_rec.id in dup_gene_ids_l:
                    continue
                # exclude duplicated genes ( end ) #
                gene_biotype_s = gene_rec.attributes["gene_biotype"][0]
                if gene_biotype_s not in ["protein_coding", "lncRNA"]:
                    gene_id_s = gene_rec.id
                    gene_chrom_s = gene_rec.chrom.split(".")[0]
                    gene_rec.chrom = gene_chrom_s
                    out_gff_file_f.write(str(gene_rec) + "\n")
                    for all_children_rec in DB.children(id=gene_id_s):
                        all_children_rec.chrom = gene_chrom_s
                        out_gff_file_f.write(str(all_children_rec) + "\n")
                else:
                    assert gene_biotype_s in ["protein_coding", "lncRNA"]
                    # 1-1. mRNA level
                    gene_id_s = gene_rec.id
                    gene_chrom_s = gene_rec.chrom.split(".")[0]
                    gene_rec.chrom = gene_chrom_s
                    RNA_2_len_nl = []
                    if gene_biotype_s == "protein_coding":
                        transcript_recs_l = [mRNA_rec for mRNA_rec in DB.children(id=gene_id_s, featuretype="mRNA")]
                        RNA_type_s = "mRNA"
                    else:
                        assert gene_biotype_s == "lncRNA"
                        transcript_recs_l = [RNA_rec for RNA_rec in DB.children(id=gene_id_s, featuretype="lnc_RNA")]
                        RNA_type_s = "lnc_RNA"
                    if len(transcript_recs_l) == 0:
                        continue
                    assert len(transcript_recs_l) >= 1
                    for RNA_rec in DB.children(id=gene_id_s, featuretype=RNA_type_s):
                        RNA_id_s = RNA_rec.id
                        # 1-2. exons
                        if gene_biotype_s == "protein_coding":
                            cds_l = sorted([cds_rec for cds_rec in DB.children(id=RNA_id_s, featuretype="CDS")],
                                           key=lambda cds_rec: cds_rec.start)
                            RNA_len_i = sum([cds_rec.end - cds_rec.start + 1 for cds_rec in cds_l])
                        else:
                            assert gene_biotype_s == "lncRNA"
                            RNA_len_i = RNA_rec.end - RNA_rec.start + 1
                        RNA_2_len_nl.append([RNA_id_s, RNA_len_i])
                    # longest transcript
                    RNA_2_len_nl = sorted(RNA_2_len_nl, key=lambda RNA_l: RNA_l[1], reverse=True)
                    if len(RNA_2_len_nl) == 0:
                        print(gene_id_s)
                        continue
                    longest_RNA_id_s = RNA_2_len_nl[0][0]
                    longest_RNA_len_i = RNA_2_len_nl[0][1]
                    assert longest_RNA_len_i > 0
                    longest_RNA_rec = DB[longest_RNA_id_s]
                    longest_RNA_rec.chrom = gene_chrom_s
                    longest_exon_recs = sorted([exon_rec for exon_rec in DB.children(id=longest_RNA_id_s, featuretype="exon")],
                                               key=lambda exon_rec: exon_rec.start)
                    for exon_rec in longest_exon_recs:
                        exon_rec.chrom = gene_chrom_s
                    longest_exon_strs = [str(exon_rec) for exon_rec in longest_exon_recs]
                    longest_cds_recs = sorted([cds_rec for cds_rec in DB.children(id=longest_RNA_id_s, featuretype="CDS")],
                                              key=lambda cds_rec: cds_rec.start)
                    for cds_rec in longest_cds_recs:
                        cds_rec.chrom = gene_chrom_s
                    longest_cds_strs = [str(cds_rec) for cds_rec in longest_cds_recs]
                    out_gff_file_f.write(str(gene_rec) + "\n")
                    out_gff_file_f.write(str(longest_RNA_rec) + "\n")
                    assert len(longest_exon_strs) >= 1
                    out_gff_file_f.write("\n".join(longest_exon_strs) + "\n")
                    if len(longest_cds_recs) >= 1:
                        assert gene_biotype_s == "protein_coding"
                        out_gff_file_f.write("\n".join(longest_cds_strs) + "\n")
                    else:
                        assert gene_biotype_s == "lncRNA"

    def CpG_island_analysis(self):
        """
        function to predict CpG islands coordinates + calculate the distance to nearest downstream protein coding genes
        :return: X
        """
        # (0) basic input files
        working_dir_s = os.path.join(result_dir, "CpG_islands")
        dir_check(working_dir_s)
        cpg_island_bed_s = os.path.join(species_fasta_dir_s, "%s.newcpgreport.sorted.bed" % ref_s)
        first_exon_bed_file_s = os.path.join(working_dir_s, "first_exons.bed")
        closest_bed_file_s = os.path.join(working_dir_s, "closest_CpG_islands_of_first_exons.bed")
        gene_2_cgi_info_dict_s = os.path.join(working_dir_s, "gene_2_cgi_info_d.pickle")
        bed_lines_l = []
        # (1) save the first exon records
        for gene_rec in ref_all_gene_recs_l:
            gene_id_s = gene_rec.id
            try:
                gene_name_s = gene_id_2_gene_name_d["ref"][gene_id_s]
            except KeyError:
                continue
            if gene_rec.strand == "+":
                reverse_bool = False
            else:
                assert gene_rec.strand == "-"
                reverse_bool = True
            first_exon_rec = sorted([exon_rec for exon_rec in refDB.children(id=gene_id_s, featuretype="exon")],
                                    key=lambda exon_rec: exon_rec.start, reverse=reverse_bool)[0]
            bedline_s = "\t".join(map(str, [gene_rec.chrom, first_exon_rec.start - 1, first_exon_rec.end,
                                             gene_name_s + "|first_exon", 0, gene_rec.strand]))
            bed_lines_l.append(bedline_s)
        with open(first_exon_bed_file_s, "w") as first_exon_bed_file_f:
            first_exon_bed_file_f.write("\n".join(bed_lines_l))
        os.system("%s sort -i %s | %s closest -id -d -D a -a stdin -b %s > %s" %
                  (bedtools_command_s, first_exon_bed_file_s, bedtools_command_s, cpg_island_bed_s, closest_bed_file_s))
        # (2) save the distance to dict
        gene_2_cgi_info_d = {}
        with open(closest_bed_file_s, "r") as closest_bed_file_f:
            for line_s in closest_bed_file_f:
                line_l = line_s.rstrip("\n").split("\t")
                gene_name_s = line_l[3].split("|")[0]
                overlap_i = abs(int(line_l[-1]))
                if line_l[6] == ".":
                    overlap_i = -1
                gene_2_cgi_info_d[gene_name_s] = [overlap_i, 0 <= overlap_i <= 3000]
        with open(gene_2_cgi_info_dict_s, "wb") as gene_2_cgi_info_dict_f:
            pickle.dump(gene_2_cgi_info_d, gene_2_cgi_info_dict_f, pickle.HIGHEST_PROTOCOL)

    def gene_structure(self, mode_s):
        """
        function to calculate GC content and not-aligned ratio (which would represent missing ratio) nearby genic regions
        :param mode_s: string: ref, prev, or tgt
        :return:
        """
        working_dir_s = os.path.join(result_dir, "gene_structure", mode_s)
        dir_check(working_dir_s)
        utr_tmp_file_s = os.path.join(working_dir_s, "%s_UTR.tmp.gff" % mode_s)
        coordinates_bed_file_s = os.path.join(working_dir_s, "%s.filter.all_exons.bed" % mode_s)
        coordinates_nuc_file_s = os.path.join(working_dir_s, "%s.filter.all_exons.nuc" % mode_s)
        final_summary_file_s = os.path.join(working_dir_s, "%s.gene_structure.tabs" % mode_s)
        gene_2_cgi_info_dict_s = os.path.join(result_dir, "CpG_islands","gene_2_cgi_info_d.pickle")
        # (0) basic input + generate DB
        # before_representative_gff_s: path to gff file of the copy of original annotatin file
        # representative_gff_s: path to gff file with only representative transcript
        # (0-1) ref: VGP assembly, VGP annotation
        if mode_s == "ref":
            original_gff_s = args.originalGFF
            before_representative_gff_s = os.path.join(working_dir_s, "%s.all.gff" % mode_s)
            representative_gff_s = os.path.join(working_dir_s, "%s.gff" % mode_s)
            fasta_file_s = ref_fasta_file_s
            chrom_2_size_d = ref_chrom_2_size_d
            os.system("cp %s %s" % (original_gff_s, before_representative_gff_s))
            before_representative_db = db_maker(gff_s=before_representative_gff_s, name_s="0.copy_before_representative", working_dir_s=working_dir_s)
            self.get_representative(DB=before_representative_db, input_gff_s=before_representative_gff_s, out_gff_s=representative_gff_s, dup_gene_ids_l=dup_gene_ids_l)
        # (0-2) prev: Prior assembly, prior annotation
        elif mode_s == "prev":
            original_gff_s = args.previousGFF
            before_representative_gff_s = os.path.join(working_dir_s, "%s.all.gff" % mode_s)
            representative_gff_s = os.path.join(working_dir_s, "%s.gff" % mode_s)
            fasta_file_s = tgt_fasta_file_s
            chrom_2_size_d = tgt_chrom_2_size_d
            os.system("cp %s %s" % (original_gff_s, before_representative_gff_s))
            before_representative_db = db_maker(gff_s=before_representative_gff_s, name_s="0.copy_before_representative", working_dir_s=working_dir_s)
            self.get_representative(DB=before_representative_db, input_gff_s=before_representative_gff_s, out_gff_s=representative_gff_s, dup_gene_ids_l=dup_gene_ids_l)
        # (0-3) tgt: Prior assembly, projected annotation
        else:
            assert mode_s == "tgt"
            original_gff_s = tarGFF_s
            representative_gff_s = os.path.join(working_dir_s, "%s.gff" % mode_s)
            fasta_file_s = tgt_fasta_file_s
            chrom_2_size_d = tgt_chrom_2_size_d
            # representative_gff_s: convert "transcript" to "mRNA" for convenience (all projected genes: protein coding genes)
            with open(original_gff_s) as original_gff_f, open(representative_gff_s,"w") as copy_gff_f:
                for line_s in original_gff_f:
                    if line_s.startswith("#"):
                        copy_gff_f.write(line_s)
                    else:
                        line_l = line_s.split("\t")
                        if line_l[2] == "transcript":
                            line_l[2] = "mRNA"
                        new_line_s = "\t".join(line_l)
                        copy_gff_f.write(new_line_s)
        copy_DB = db_maker(gff_s=representative_gff_s, name_s="0.copy", working_dir_s=working_dir_s)

        # (1) generating UTRs
        if file_check(utr_tmp_file_s) is False:
            utr_lines_l = self.add_UTR_manually(DB=copy_DB, mode_s=mode_s)
            with open(utr_tmp_file_s, "w") as utr_tmp_file_f:
                utr_tmp_file_f.write("\n".join(utr_lines_l))
            copy_DB.update(data=utr_tmp_file_s, merge_strategy="create_unique")

        # (2) generate intron
        intron_num_i = len([intron_rec for intron_rec in copy_DB.all_features(featuretype="intron")])
        if intron_num_i == 0:
            print("intron count (before): %d" % intron_num_i)
            introns = list(copy_DB.create_introns(exon_featuretype="exon", grandparent_featuretype="gene"))
            copy_DB.update(introns, merge_strategy="create_unique")
            intron_num_i = len([intron_rec for intron_rec in copy_DB.all_features(featuretype="intron")])
            print("intron count (after): %d" % intron_num_i)

        # (3) get CpG islands
        if file_check(gene_2_cgi_info_dict_s) is False:
            self.CpG_island_analysis()
        with open(gene_2_cgi_info_dict_s, "rb") as gene_2_cgi_info_dict_f:
            gene_2_cgi_info_d = pickle.load(gene_2_cgi_info_dict_f)

        # (4) generate bed file
        count_features_d, gene_2_feature_num_d = a.extract_coordinates(db=copy_DB, coordinates_bed_file_s=coordinates_bed_file_s,
                                                                       chrom_2_size_d=chrom_2_size_d, gene_2_cgi_info_d=gene_2_cgi_info_d, mode_s=mode_s)
        # (4-1) summary of total coordinates count
        count_features_file_s = os.path.join(working_dir_s, "count_feature.tsv")
        with open(count_features_file_s, "w") as count_features_file_f:
            gene_biotype_l = ["protein_coding", "lncRNA", "rRNA", "tRNA", "snRNA", "snoRNA"]
            features_l = ["genes", "filtered_genes", "excluded_genes", "partial_genes",
                          "genes_with_less_than_three_exons/introns", "genes_with_not_enough_upstream/downstream_seq(3kb)",
                          "exons", "introns", "cds", "5'utr","3'utr"]
            # header line
            count_features_file_f.write("biotype\t" + "\t".join(features_l) + "\n")
            for gene_biotype_s in gene_biotype_l:
                line_l = [gene_biotype_s]
                for feature_s in features_l:
                    line_l.append(count_features_d[gene_biotype_s][feature_s])
                count_features_file_f.write("\t".join(map(str, line_l))+"\n")
        # (4-2) count coordinates for each gene
        gene_2_feature_num_file_s = os.path.join(working_dir_s, "gene_2_feature_num.tsv")
        with open(gene_2_feature_num_file_s, "w") as gene_2_feature_num_file_f:
            features_l = ["biotype", "five-prime-UTR", "FirstCDS", "InternalCDS", "LastCDS", "three-prime-UTR",
                          "5'UTR-intron", "First-intron", "Internal-intron", "Last-intron", "3'UTR-intron", "FirstCDS_LastCDS_intron",
                          "FirstExon", "InternalExon", "LastExon",
                          "FirstExon_InternalExon_intron", "InternalExon_InternalExon_intron", "InternalExon_LastExon_intron","FirstExon_LastExon_intron"]
            # header line
            gene_2_feature_num_file_f.write("gene_id\t" + "\t".join(features_l) + "\n")
            for gene_id in sorted(list(gene_2_feature_num_d.keys())):
                line_l = [gene_id]
                for feature_s in features_l:
                    line_l.append(gene_2_feature_num_d[gene_id][feature_s])
                gene_2_feature_num_file_f.write("\t".join(map(str, line_l))+"\n")

        # (5) summary the result
        header_l = ["chrom", "start", "end", "general_type", "source_gene", "gene_biotype", "order", "count",
                    "position", "coding_exon_count", "excepted", "excepted_types",
                    "distance_to_upstream_cgi", "cgi_within_3kb", "zero", "strand", "ATcontent", "GCcontent",
                    "numA", "numC", "numG", "numT", "numN", "numOther", "length"]
        general_types_l = ["three-prime-UTR", "five-prime-UTR", "cds","intron","exon",
                           "three-prime-downstream", "five-prime-upstream", "three-inside","five-inside"]
        os.system("%s nuc -fi %s -bed %s > %s" % (bedtools_command_s, fasta_file_s, coordinates_bed_file_s, coordinates_nuc_file_s))
        with open(final_summary_file_s,"w") as final_summary_file_f, open(coordinates_nuc_file_s) as all_exons_nuc_file_f:
            final_summary_file_f.write("\t".join(header_l) + "\n")
            all_exons_nuc_file_f.readline()
            for line_s in all_exons_nuc_file_f:
                new_line_l = line_s.rstrip("\n").replace("|","\t").split("\t")
                general_type_s = new_line_l[3]
                if general_type_s not in general_types_l:
                    continue
                a_i, c_i, g_i, t_i = map(int, new_line_l[18:22])
                acgt_i = sum([a_i, c_i, g_i, t_i])
                at_content_wo_n_fl = float((a_i + t_i) / acgt_i) if acgt_i > 0 else "ALL_filled_N"
                gc_content_wo_n_fl = float((g_i + c_i) / acgt_i) if acgt_i > 0 else "ALL_filled_N"
                new_line_l[16] = at_content_wo_n_fl
                new_line_l[17] = gc_content_wo_n_fl
                new_line_s = "\t".join(map(str,new_line_l)) + "\n"
                final_summary_file_f.write(new_line_s)

    def extract_coordinates(self, db, coordinates_bed_file_s, chrom_2_size_d, gene_2_cgi_info_d, mode_s):
        """
        function to extract exon/intron/CDS/5'UTR/3'UTR/flanking sequences's coordinates
        :param db: database made from gffutils
        :param coordinates_bed_file_s: string, path to output bed file
        :param chrom_2_size_d: dict, key: chromosome name (string), value: chromosome size (int)
        :param gene_2_cgi_info_d: dict, key: gene name (string), value: [dist. to upstream CGI , boolean (dist. to upstream CGI <= 3kb?)] (list)
        :param mode_s: string, mode (ref: VGP or tgt: prior)
        :return:
        """
        # basic input variables
        gene_recs_l = [gene_rec for gene_rec in db.all_features(featuretype="gene")]
        bed_lines_l = []
        gene_2_feature_num_d = {}
        gene_classify_d = {1: "OneCodingExon", 2: "TwoCodingExon", 3: "ThreeCodingExon"}
        gene_types_l = ["protein_coding", "lncRNA", "rRNA", "snoRNA", "snRNA", "tRNA"]
        # dictionary to rename introns: (1) protein coding genes: with UTR and coding exons, (2) other genes: with exons only
        rename_intron_d = {"five-prime-UTR_FirstCDS_intron": "5'UTR-intron", "five-prime-UTR_five-prime-UTR_intron": "5'UTR-intron",
                           "LastCDS_three-prime-UTR_intron": "3'UTR-intron", "FirstCDS_three-prime-UTR_intron": "3'UTR-intron", "three-prime-UTR_three-prime-UTR_intron": "3'UTR-intron", "FirstCDS_InternalCDS_intron": "First-intron",
                           "InternalCDS_InternalCDS_intron": "Internal-intron", "InternalCDS_LastCDS_intron": "Last-intron", "FirstCDS_LastCDS_intron": "FirstCDS_LastCDS_intron",
                           "FirstExon_InternalExon_intron": "FirstExon_InternalExon_intron",
                           "InternalExon_InternalExon_intron": "InternalExon_InternalExon_intron",
                           "InternalExon_LastExon_intron": "InternalExon_LastExon_intron",
                           "FirstExon_LastExon_intron": "FirstExon_LastExon_intron"}
        # for protein coding genes
        coding_gene_positions_l = ["FirstCDS", "InternalCDS", "LastCDS", "5'UTR-intron", "First-intron", "Internal-intron", "Last-intron", "3'UTR-intron", "FirstCDS_LastCDS_intron"]
        # for other type of genes
        gene_positions_l = ["FirstExon", "InternalExon", "LastExon", "FirstExon_InternalExon_intron", "InternalExon_InternalExon_intron", "InternalExon_LastExon_intron", "FirstExon_LastExon_intron"]
        count_features_d = {}
        for gene_type_s in gene_types_l:
            count_features_d[gene_type_s] = {"genes": 0, "filtered_genes": 0, "excluded_genes": 0,
                                             "partial_genes": 0, "genes_with_less_than_three_exons/introns": 0,
                                             "genes_with_not_enough_upstream/downstream_seq(3kb)": 0,
                                             "exons": 0, "introns": 0, "cds":0, "5'utr":0, "3'utr":0}
            for position_s in coding_gene_positions_l:
                count_features_d[gene_type_s][position_s] = 0
            for position_s in gene_positions_l:
                count_features_d[gene_type_s][position_s] = 0
        intron_positions_l = list(rename_intron_d.keys())

        # for loop: each gene record, generate bed lines
        for gene_rec in gene_recs_l:
            # (0) for each gene: basic info
            gene_id_s = gene_rec.id
            assert "|" not in gene_id_s
            gene_biotype_s = db[gene_id_s].attributes["gene_biotype"][0]
            if gene_biotype_s not in gene_types_l:
                continue
            # (0-1) for ref & tgt --> using pre-assigned gene name from gff DB
            if gene_biotype_s == "protein_coding":
                if mode_s == "ref":
                    try:
                        gene_name_s = gene_id_2_gene_name_d["ref"][gene_id_s]
                    except KeyError:  # Exclude duplicated gene records
                        print(gene_id_s)
                        continue
                elif mode_s == "prev":
                    gene_name_s = gene_id_s
                else:
                    assert mode_s == "tgt"
                    try:
                        gene_name_s = gene_id_2_gene_name_d["tgt"][gene_id_s]
                    except KeyError:  # Exclude duplicated gene records
                        print(gene_id_s)
                        continue
            else:
                assert gene_biotype_s in ["lncRNA", "rRNA", "snoRNA", "snRNA", "tRNA"]
                gene_name_s = gene_id_s
            # (0-2) CGI info
            distance_to_cgi_i, within_3kb_bool = gene_2_cgi_info_d.get(gene_name_s, ["NA", "NA"])
            # (0-3) CGI info
            count_features_d[gene_biotype_s]["genes"] += 1
            if gene_biotype_s == "protein_coding":
                transcript_recs_l = [mRNA_rec for mRNA_rec in db.children(id=gene_id_s, featuretype="mRNA")]
            else:
                assert gene_biotype_s in ["lncRNA", "rRNA", "snoRNA", "snRNA", "tRNA"]
                rna_biotype_s = "lnc_RNA" if gene_biotype_s == "lncRNA" else gene_biotype_s
                transcript_recs_l = [mRNA_rec for mRNA_rec in db.children(id=gene_id_s, featuretype=rna_biotype_s)]
            assert len(transcript_recs_l) == 1
            # (1) basic information from gene record:
            gene_strand_s = gene_rec.strand
            gene_chrom_s = gene_rec.chrom
            gene_start_i, gene_end_i = gene_rec.start - 1, gene_rec.end
            gene_2_feature_num_d[gene_id_s] = {"biotype": gene_biotype_s,
                                               "five-prime-UTR":0, "FirstCDS": 0, "InternalCDS": 0, "LastCDS": 0, "three-prime-UTR": 0,
                                               "5'UTR-intron": 0, "First-intron": 0, "Internal-intron": 0, "Last-intron": 0, "3'UTR-intron": 0, "FirstCDS_LastCDS_intron": 0,
                                               "FirstExon": 0, "InternalExon": 0, "LastExon": 0,
                                               "FirstExon_InternalExon_intron": 0, "InternalExon_InternalExon_intron": 0, "InternalExon_LastExon_intron": 0, "FirstExon_LastExon_intron": 0}
            # (2) get sub records:
            # (2-1) exons, UTRs, CDS, introns: sorted based on the start position
            cds_recs_l = sorted([cds_rec for cds_rec in db.children(id=gene_id_s, featuretype="CDS")], key=lambda cds_rec: cds_rec.start)
            exon_recs_l = sorted([exon_rec for exon_rec in db.children(id=gene_id_s, featuretype="exon")], key=lambda exon_rec: exon_rec.start)
            utr5_recs_l = sorted([utr5_rec for utr5_rec in db.children(id=gene_id_s, featuretype="five_prime_UTR")], key=lambda utr5_rec: utr5_rec.start)
            utr3_recs_l = sorted([utr3_rec for utr3_rec in db.children(id=gene_id_s, featuretype="three_prime_UTR")], key=lambda utr3_rec: utr3_rec.start)
            intron_recs_l = sorted([intron_rec for intron_rec in db.children(id=gene_id_s, featuretype="intron")], key=lambda intron_rec: intron_rec.start)
            intron_recs_over_10bp_l = [intron_rec for intron_rec in intron_recs_l if (intron_rec.end - intron_rec.start + 1) > 10]
            intron_recs_with_type_l = [[intron_rec.start, "intron", intron_rec] for intron_rec in intron_recs_over_10bp_l]

            # (3) check exception: number of coding exons/introns, whether located within 3kbp from the end of scaffold, or partial gene?
            gene_exception_types_l = []
            # (3-1) First, internal, last coding exons & introns ==> at least three introns
            if len(intron_recs_over_10bp_l) < 3 and gene_biotype_s == "protein_coding":
                gene_exception_types_l.append("less_than_three_introns")
            if len(intron_recs_over_10bp_l) < 2 and gene_biotype_s == "lncRNA":
                gene_exception_types_l.append("less_than_three_introns")
            # (3-2) Within the 3kbp end of scaffolds
            if (gene_start_i - 3000 < 0) or (gene_end_i + 3000 > chrom_2_size_d[gene_chrom_s]):
                gene_exception_types_l.append("no_3Kb_sequences")
                count_features_d[gene_biotype_s]["genes_with_not_enough_upstream/downstream_seq(3kb)"] += 1
            # (3-3) partial genes
            try:
                partial_flag_s = db[gene_id_s].attributes["partial"][0]
                gene_exception_types_l.append(partial_flag_s)
                count_features_d[gene_biotype_s]["partial_genes"] += 1
            except KeyError:
                partial_flag_s = ""
                pass

            # (0) CDS
            gene_exontype_s = "no_CDS_available"
            cds_recs_with_type_l = [] # would be used for intron records: using GFF-style coordinates
            adj_cds_count_i = 0
            if gene_biotype_s == "protein_coding":
                cds_order_i = 0
                adj_cds_start_l, adj_cds_end_l = ignore_1_or_2bp_introns(raw_rec_l=cds_recs_l)
                adj_cds_count_i = len(adj_cds_start_l)
                if adj_cds_count_i <= 3:
                    gene_exception_types_l.append("less_than_three_CDS")
                gene_exontype_s = gene_classify_d.get(adj_cds_count_i, "more_than_four_coding_exon")
                for (adj_cds_start_i, adj_cds_end_i) in zip(adj_cds_start_l, adj_cds_end_l):
                    cds_exception_types_l = copy.deepcopy(gene_exception_types_l)
                    cds_order_i += 1
                    cds_start_i, cds_end_i = adj_cds_start_i - 1, adj_cds_end_i
                    cds_position_type_s = self.check_position(order_i=cds_order_i, count_i=adj_cds_count_i, strand_c=gene_strand_s, mode_s="CDS")
                    gene_2_feature_num_d[gene_id_s][cds_position_type_s] += 1
                    cds_recs_with_type_l.append([adj_cds_start_i, cds_position_type_s])
                    cds_exception_types_s = ",".join(map(str, cds_exception_types_l))
                    cds_exception_bool_s = "True" if len(cds_exception_types_l) == 0 else "False"
                    cds_level_id_s = "|".join(map(str, ["cds", gene_id_s, gene_biotype_s,
                                                        cds_order_i, adj_cds_count_i, cds_position_type_s,
                                                        gene_exontype_s, cds_exception_bool_s, cds_exception_types_s,
                                                        distance_to_cgi_i, within_3kb_bool]))
                    cds_bed_line_s = "\t".join(map(str, [gene_chrom_s, cds_start_i, cds_end_i, cds_level_id_s, 0, gene_strand_s]))
                    bed_lines_l.append(cds_bed_line_s)

            # (1) Exon
            exon_order_i = 0
            exon_recs_with_type_l = []  # would be used for intron records: using GFF-style coordinates
            adj_exon_start_l, adj_exon_end_l = ignore_1_or_2bp_introns(raw_rec_l=exon_recs_l)
            adj_exon_count_i = len(adj_exon_start_l)
            if adj_exon_count_i < 3 and gene_biotype_s in ["protein_coding", "lncRNA"]:
                gene_exception_types_l.append("less_than_three_exon")
            for (adj_exon_start_i, adj_exon_end_i) in zip(adj_exon_start_l, adj_exon_end_l):
                exon_exception_types_l = copy.deepcopy(gene_exception_types_l)
                exon_order_i += 1
                exon_start_i, exon_end_i = adj_exon_start_i - 1, adj_exon_end_i
                exon_position_type_s = self.check_position(order_i=exon_order_i, count_i=adj_exon_count_i, strand_c=gene_strand_s, mode_s="Exon")
                gene_2_feature_num_d[gene_id_s][exon_position_type_s] += 1
                exon_recs_with_type_l.append([exon_start_i + 1, exon_position_type_s])
                exon_exception_types_s = ",".join(map(str, exon_exception_types_l))
                exon_exception_bool_s = "True" if len(exon_exception_types_l) == 0 else "False"
                exon_level_id_s = "|".join(map(str, ["exon", gene_id_s, gene_biotype_s,
                                                     exon_order_i, adj_exon_count_i, exon_position_type_s,
                                                     gene_exontype_s, exon_exception_bool_s, exon_exception_types_s,
                                                     distance_to_cgi_i, within_3kb_bool]))
                exon_bed_line_s = "\t".join(map(str, [gene_chrom_s, exon_start_i, exon_end_i, exon_level_id_s, 0, gene_strand_s]))
                bed_lines_l.append(exon_bed_line_s)

            # (2) UTR
            utr_recs_with_type_l = []
            for (utr_position_type_s, utr_recs_l) in zip(["five-prime-UTR", "three-prime-UTR"], [utr5_recs_l, utr3_recs_l]):
                utr_order_i = 0
                tmp_utr_recs_with_type_l = []
                utr_count_i = len(utr_recs_l)
                if utr_count_i > 0:
                    adj_utr_start_l, adj_utr_end_l = ignore_1_or_2bp_introns(raw_rec_l=utr_recs_l)
                    adj_utr_count_i = len(adj_utr_start_l)
                    assert utr_count_i == adj_utr_count_i  # Expectation: indels only for CDSs
                    for (adj_utr_start_i, adj_utr_end_i) in zip(adj_utr_start_l, adj_utr_end_l):
                        utr_exception_types_l = copy.deepcopy(gene_exception_types_l)
                        utr_order_i += 1
                        utr_start_i, utr_end_i = adj_utr_start_i - 1, adj_utr_end_i
                        # add exceptions
                        if adj_cds_count_i <= 3:
                            assert ("less_than_three_CDS" in utr_exception_types_l)
                        tmp_utr_recs_with_type_l.append([adj_utr_start_i, utr_position_type_s])
                        gene_2_feature_num_d[gene_id_s][utr_position_type_s] += 1
                        utr_exception_types_s = ",".join(map(str, utr_exception_types_l))
                        utr_exception_bool_s = "True" if len(utr_exception_types_l) == 0 else "False"
                        utr_level_id_s = "|".join(map(str, [utr_position_type_s, gene_id_s, gene_biotype_s,
                                                            utr_order_i, utr_count_i, utr_position_type_s,
                                                            gene_exontype_s, utr_exception_bool_s,utr_exception_types_s,
                                                            distance_to_cgi_i, within_3kb_bool]))
                        utr_bed_line_s = "\t".join(map(str, [gene_chrom_s, utr_start_i, utr_end_i, utr_level_id_s, 0, gene_strand_s]))
                        bed_lines_l.append(utr_bed_line_s)
                    utr_recs_with_type_l += tmp_utr_recs_with_type_l

            # (3) Intron
            intron_order_i = 0
            intron_count_i = len(intron_recs_over_10bp_l)
            if gene_biotype_s == "protein_coding":  # for protein coding genes
                assert len(cds_recs_with_type_l) >= 1
                if gene_strand_s == "+":
                    exon_n_intron_l = sorted(cds_recs_with_type_l + utr_recs_with_type_l + intron_recs_with_type_l, key=lambda info_l: info_l[0])
                else:
                    assert gene_strand_s == "-"
                    exon_n_intron_l = sorted(cds_recs_with_type_l + utr_recs_with_type_l + intron_recs_with_type_l, key=lambda info_l: info_l[0], reverse=True)
            else:  # for non-protein coding genes
                if gene_strand_s == "+":
                    exon_n_intron_l = sorted(exon_recs_with_type_l + intron_recs_with_type_l, key=lambda info_l: info_l[0])
                else:
                    assert gene_strand_s == "-"
                    exon_n_intron_l = sorted(exon_recs_with_type_l + intron_recs_with_type_l, key=lambda info_l: info_l[0], reverse=True)
            rec_index_i = 0
            for rec_l in exon_n_intron_l:
                if rec_l[1] == "intron":
                    intron_exception_types_l = copy.deepcopy(gene_exception_types_l)
                    intron_start_i, intron_end_i = rec_l[-1].start - 1, rec_l[-1].end  # formatted as BED
                    assert (intron_end_i - intron_start_i) > 10
                    intron_order_i += 1
                    # get position information: 5'UTR intron, First intron, middle intron, last intron, 3'UTR intron
                    upstream_index_i = rec_index_i - 1
                    downstream_index_i = rec_index_i + 1
                    assert (upstream_index_i >= 0 and downstream_index_i < len(exon_n_intron_l))
                    upstream_type_s = exon_n_intron_l[upstream_index_i][1]
                    downstream_type_s = exon_n_intron_l[downstream_index_i][1]
                    intron_position_type_s = "%s_%s_intron" % (upstream_type_s, downstream_type_s)
                    try:
                        assert intron_position_type_s in intron_positions_l
                    except:
                        print(exon_n_intron_l)
                        continue
                    gene_2_feature_num_d[gene_id_s][rename_intron_d[intron_position_type_s]] += 1
                    intron_exception_types_s = ",".join(map(str, intron_exception_types_l))
                    intron_exception_bool_s = "True" if len(intron_exception_types_l) == 0 else "False"
                    intron_level_id_s = "|".join(map(str, ["intron", gene_id_s, gene_biotype_s,
                                                           intron_order_i, intron_count_i, intron_position_type_s,
                                                           gene_exontype_s, intron_exception_bool_s, intron_exception_types_s,
                                                           distance_to_cgi_i, within_3kb_bool]))
                    intron_bed_line_s = "\t".join(map(str, [gene_chrom_s, intron_start_i, intron_end_i, intron_level_id_s, 0, gene_strand_s]))
                    bed_lines_l.append(intron_bed_line_s)
                rec_index_i += 1

            # (4) up/down
            # 5' upstream / 3' downstream : 3kbp
            exception_types_l = copy.deepcopy(gene_exception_types_l)
            exception_types_s = ",".join(map(str, exception_types_l))
            exception_bool_s = "True" if len(exception_types_l) == 0 else "False"
            for i in range(0, 30):
                start_i = i * 100
                end_i = (i + 1) * 100
                if gene_strand_s == "+":
                    upstream_start_i = max(db[gene_id_s].start - 1 - end_i, 0)
                    upstream_end_i = max(db[gene_id_s].start - 1 - start_i, 1)
                    downstream_start_i = min(db[gene_id_s].end + start_i, chrom_2_size_d[gene_chrom_s] - 1)
                    downstream_end_i = min(db[gene_id_s].end + end_i, chrom_2_size_d[gene_chrom_s])
                    upinside_start_i = min(db[gene_id_s].start - 1 + start_i, chrom_2_size_d[gene_chrom_s]-1)
                    upinside_end_i = min(db[gene_id_s].start - 1 + end_i, chrom_2_size_d[gene_chrom_s])
                    downinside_start_i = max(db[gene_id_s].end - end_i, 0)
                    downinside_end_i = max(db[gene_id_s].end - start_i, 1)
                else:
                    assert gene_strand_s == "-"
                    upstream_start_i = min(db[gene_id_s].end + start_i, chrom_2_size_d[gene_chrom_s] - 1)
                    upstream_end_i = min(db[gene_id_s].end + end_i, chrom_2_size_d[gene_chrom_s])
                    downstream_start_i = max(db[gene_id_s].start - 1 - end_i, 0)
                    downstream_end_i = max(db[gene_id_s].start - 1 - start_i, 1)
                    upinside_start_i = max(db[gene_id_s].end - end_i, 0)
                    upinside_end_i = max(db[gene_id_s].end - start_i, 1)
                    downinside_start_i = min(db[gene_id_s].start - 1 + start_i, chrom_2_size_d[gene_chrom_s]-1)
                    downinside_end_i = min(db[gene_id_s].start - 1 + end_i, chrom_2_size_d[gene_chrom_s])
                upstream_position_type_s = "five-prime-%sbp" % (str(end_i))
                downstream_position_type_s = "three-prime-%sbp" % (str(end_i))
                upinside_position_type_s = "five-inside-%sbp" % (str(end_i))
                downinside_position_type_s = "three-inside-%sbp" % (str(end_i))
                upstream_level_id_s = "|".join(map(str, ["five-prime-upstream", gene_id_s, gene_biotype_s,
                                                         30 - i, 30, upstream_position_type_s,
                                                         gene_exontype_s, exception_bool_s, exception_types_s,
                                                         distance_to_cgi_i, within_3kb_bool]))
                downstream_level_id_s = "|".join(map(str, ["three-prime-downstream", gene_id_s, gene_biotype_s,
                                                           i, 30, downstream_position_type_s,
                                                           gene_exontype_s, exception_bool_s, exception_types_s,
                                                           distance_to_cgi_i, within_3kb_bool]))
                upinside_level_id_s = "|".join(map(str, ["five-inside", gene_id_s, gene_biotype_s,
                                                         30 - i, 30, upinside_position_type_s,
                                                         gene_exontype_s, exception_bool_s, exception_types_s,
                                                         distance_to_cgi_i, within_3kb_bool]))
                downinside_level_id_s = "|".join(map(str, ["three-inside", gene_id_s, gene_biotype_s,
                                                           i, 30, downinside_position_type_s,
                                                           gene_exontype_s, exception_bool_s, exception_types_s,
                                                           distance_to_cgi_i, within_3kb_bool]))
                upstream_bed_line_s = "\t".join(map(str,[gene_chrom_s, upstream_start_i, upstream_end_i, upstream_level_id_s, 0, gene_strand_s]))
                downstream_bed_line_s = "\t".join(map(str, [gene_chrom_s, downstream_start_i, downstream_end_i, downstream_level_id_s, 0, gene_strand_s]))
                upinside_bed_line_s = "\t".join(map(str,[gene_chrom_s, upinside_start_i, upinside_end_i, upinside_level_id_s, 0, gene_strand_s]))
                downinside_bed_line_s = "\t".join(map(str, [gene_chrom_s, downinside_start_i, downinside_end_i, downinside_level_id_s, 0, gene_strand_s]))
                bed_lines_l.append(upstream_bed_line_s)
                bed_lines_l.append(downstream_bed_line_s)
                bed_lines_l.append(upinside_bed_line_s)
                bed_lines_l.append(downinside_bed_line_s)

            gene_exception_types_s = ",".join(map(str, gene_exception_types_l))
            if ("less_than_three_introns" in gene_exception_types_l) or ("less_than_three_CDS" in gene_exception_types_l):
                count_features_d[gene_biotype_s]["genes_with_less_than_three_exons/introns"] += 1
            gene_exception_bool_s = "True" if len(gene_exception_types_l) == 0 else "False"
            if gene_exception_bool_s == "True":
                count_features_d[gene_biotype_s]["filtered_genes"] += 1
                count_features_d[gene_biotype_s]["cds"] += adj_cds_count_i
                count_features_d[gene_biotype_s]["introns"] += intron_count_i
                count_features_d[gene_biotype_s]["5'utr"] += len(utr5_recs_l)
                count_features_d[gene_biotype_s]["3'utr"] += len(utr3_recs_l)
                count_features_d[gene_biotype_s]["exons"] += adj_exon_count_i
                for position_s in gene_positions_l:
                    count_features_d[gene_biotype_s][position_s] += gene_2_feature_num_d[gene_id_s][position_s]
                for position_s in coding_gene_positions_l:
                    count_features_d[gene_biotype_s][position_s] += gene_2_feature_num_d[gene_id_s][position_s]
            else:
                assert gene_exception_bool_s == "False"
                count_features_d[gene_biotype_s]["excluded_genes"] += 1
            gene_level_id_s = "|".join(map(str, ["gene", gene_id_s, gene_biotype_s,
                                                 1, 1, "default",
                                                 gene_exontype_s, gene_exception_bool_s, gene_exception_types_s,
                                                 distance_to_cgi_i, within_3kb_bool]))
            gene_bed_line_s = "\t".join(map(str, [gene_chrom_s, gene_start_i, gene_end_i, gene_level_id_s, 0, gene_strand_s]))
            bed_lines_l.append(gene_bed_line_s)

        with open(coordinates_bed_file_s, "w") as all_exons_bed_file_f:
            all_exons_bed_file_f.write("\n".join(bed_lines_l) + "\n")
        return count_features_d,  gene_2_feature_num_d

    def check_position(self, order_i, count_i, strand_c, mode_s):
        """
        function to check the order (position, First/Internal/Last) of exon or CDS
        :param order_i: integer, numeric order of the exon or CDS
        :param count_i: integer, total number of exons or CDSs
        :param strand_c: character, strand info
        :param mode_s: string, CDS or Exon
        :return: position_type_s: string, type of position (First/Internal/Last)
        """
        if order_i == 1:
            # considering single coding exon
            if count_i == 1:
                position_type_s = "First%s" % mode_s
            else:
                position_type_s = "First%s" % mode_s if strand_c == "+" else "Last%s" % mode_s
        elif 1 < order_i < count_i:
            position_type_s = "Internal%s" % mode_s
        else:
            assert order_i == count_i
            position_type_s = "Last%s" % mode_s if strand_c == "+" else "First%s" % mode_s
        return position_type_s
    #####################################
    # 1. gene structure analysis (done) #
    #####################################

if __name__ == '__main__':
    if (len(sys.argv) <= 1):
        print("please refer the manual with -h or --help")
        exit()

    work_explanation_s = "choose one of these submodules: allinone, missing, split, frameshift, prematurestopcodon, intronexonjunctiondisruption, nincodingregion, " \
                         "                                summary, circos_maker, gene_structure, false_variant" \
                         "allinone: run all submodules in once" \
                         "missing: analyze missing genes & missing genomic regions" \
                         "split: analyze fragmented & intra-scaffold split genes" \
                         "frameshift: analyze frameshift genes" \
                         "prematurestopcodon: analyze premature stop codon genes" \
                         "intronexonjunctiondisruption: analyze splicing junction disruption genes" \
                         "nincodingregion: analyze N in coding region genes" \
                         "summary: after anlyzing all above modules, start to summarize the results" \
                         "circos_maker: module to generate circos plot" \
                         "gene_structure: module to analyze gene structure-specific characteristics" \
                         "false_variant: module to summary mpileup result (genome-wide)"

    parser = argparse.ArgumentParser(
    description="Usage: [python3] FalseGeneLoss.py")
    parser.add_argument("-w", "--work", action="store", help=work_explanation_s, default="")
    parser.add_argument("-r", "--reference", action="store", help="nickname of VGP assembly. Choose one of these: vgp_taegut, vgp_calann, vgp_ornana, vgp_anates", default=False)
    parser.add_argument("-t", "--target", action="store", help="nickname of prior assembly. Choose one of these: pre_taegut, pre_calann, pre_ornana, pre_anates", default=False)
    parser.add_argument("-p", "--processNumber", action="store", help="number of processors to use", default=15)
    parser.add_argument("-hl", "--hal", action="store", help="path to hal file", default=False)
    parser.add_argument("-g", "--gff3", action="store", help="path to VGP annotation GFF file used in CAT", default=False)
    parser.add_argument("-rb", "--refBam", action="store", help="path to VGP BAM file", default=False)
    parser.add_argument("-tb", "--tgtBam", action="store", help="path to prior BAM file", default=False)
    parser.add_argument("-rm", "--remap", action="store", help="path to NCBI remap file", default=False)
    parser.add_argument("-og", "--originalGFF", action="store", help="path to orignial VGP annotation GFF file", default=False)
    parser.add_argument("-pg", "--previousGFF", action="store", help="path to prior annotation GFF file", default = False)
    parser.add_argument("-id", "--IndelFromVG", action="store", help="path to indel bed file from VG", default=False)
    parser.add_argument("-snp", "--SnpFromVG", action="store", help="path to SNP bed file from VG", default = False)
    parser.add_argument("-fd", "--fastaDir", action="store", help="path to fasta directory", default = False)
    parser.add_argument("-dg", "--duplicatedGenes", action="store", help="path to duplicated gene file", default = False)

    args = parser.parse_args()

    if (args.reference is False) or (args.target is False) or (args.hal is False):
        print("please check the input/reference/target options. Exit")
        sys.exit()
    print("Submodule: %s" % args.work)
    if int(args.processNumber) == int:
        print("Number of process : %s" % str(int(args.processNumber)))

    a = FalseGeneLoss()
    # from input files
    fasta_dir_s = str(args.fastaDir)
    hal_path_s = str(args.hal)
    ref_s = str(args.reference)
    tgt_s = str(args.target)
    work_s = str(args.work)
    process_num_i = int(args.processNumber)
    ref_original_gff_s = str(args.originalGFF)
    ref_bamfile_s = str(args.refBam)
    tgt_bamfile_s = str(args.tgtBam)
    remap_excel_file_s = str(args.remap)
    refGFF_s = str(args.gff3)
    indel_from_vg_file_s = str(args.IndelFromVG)
    snp_from_vg_file_s = str(args.SnpFromVG)
    duplication_gene_file_s = str(args.duplicatedGenes)

    # basic variables (paths and dicts)
    temp_dir = "./"
    db_dir = os.path.join(temp_dir, "out", "databases")
    cons_dir = os.path.join(temp_dir, "out", "consensus_gene_set")
    result_dir = os.path.join(temp_dir, "CAT2MISSING_%s" % tgt_s)
    tarGFF_s = os.path.join(cons_dir, "%s.converted.gff3" % tgt_s)
    species_fasta_dir_s = os.path.join(fasta_dir_s, ref_s.split("_")[1])
    ref_assembly_report_file_s = os.path.join(species_fasta_dir_s, "%s.assembly_report.txt" % ref_s)
    ref_gap_info_file_s = os.path.join(species_fasta_dir_s, "%s.genomic_gaps.txt" % ref_s)
    tgt_assembly_report_file_s = os.path.join(species_fasta_dir_s, "%s.assembly_report.txt" % tgt_s)
    tgt_gap_info_file_s = os.path.join(species_fasta_dir_s, "%s.genomic_gaps.txt" % tgt_s)
    ref_fasta_file_s = os.path.join(species_fasta_dir_s, "%s.fasta" % ref_s)
    tgt_fasta_file_s = os.path.join(species_fasta_dir_s, "%s.fasta" % tgt_s)
    ref_chrom_2_size_d = scaffold_2_size_dictionary(mode_s="ref")
    tgt_chrom_2_size_d = scaffold_2_size_dictionary(mode_s="tgt")
    print("(1) start to generate DB from annotation...")
    refDB = db_maker(gff_s=refGFF_s, name_s="ref", working_dir_s=result_dir)
    tgtDB = db_maker(gff_s=tarGFF_s, name_s="tar", working_dir_s=result_dir)
    print(" --> finished to generate pickle from gffDBs...")
    print("(2) start to generate pickle from gffDBs...")
    dup_gene_ids_l = duplication_gene_list(duplication_gene_file_s=duplication_gene_file_s)
    ref_all_gene_recs_l, ref_all_trans_recs_l = gene_and_transcript_records(mode_s="ref", dup_gene_ids_l=dup_gene_ids_l)
    tgt_all_gene_recs_l, tgt_all_trans_recs_l = gene_and_transcript_records(mode_s="tgt", dup_gene_ids_l=dup_gene_ids_l)
    ref_trans_id_2_gene_id_d = ref_trans_id_2_gene_id()
    gene_name_2_gene_id_d, gene_id_2_gene_name_d = gene_name_2_gene_id()
    aln_id_2_gene_id_d, gene_id_2_aln_id_d = aln_id_2_gene_id_dictionary()
    ref_gene_name_2_intron_d, ref_intron_recs_l = a.generate_intron_records_from_ref_assembly()
    print(" --> finished to generate pickle from gffDBs...")
    if ref_s == "vgp_taegut":
        scaff_2_chrom_d = {"NC_044211": "chr.1", "NC_044212": "chr.1A", "NC_044213": "chr.2", "NC_044214": "chr.3",
                           "NC_044215": "chr.4", "NC_044216": "chr.4A", "NC_044217": "chr.5", "NC_044218": "chr.6",
                           "NC_044219": "chr.7", "NC_044220": "chr.8", "NC_044221": "chr.9", "NC_044222": "chr.10",
                           "NC_044223": "chr.11", "NC_044224": "chr.12", "NC_044225": "chr.13", "NC_044226": "chr.14",
                           "NC_044227": "chr.15", "NC_044228": "chr.16", "NC_044229": "chr.17", "NC_044230": "chr.18",
                           "NC_044231": "chr.19", "NC_044232": "chr.20", "NC_044233": "chr.21", "NC_044234": "chr.22",
                           "NC_044235": "chr.23", "NC_044236": "chr.24", "NC_044237": "chr.25", "NC_044238": "chr.26",
                           "NC_044239": "chr.27", "NC_044240": "chr.28", "NC_044241": "chr.Z", "NC_044242": "chr.29",
                           "NC_044243": "chr.30", "NW_022045287": "chr.5,u1", "NW_022045292": "chr.31", "NW_022045293": "chr.35,u1",
                           "NW_022045299": "chr.11,u1", "NW_022045306": "chr.2,u1", "NW_022045315": "chr.26,u1",
                           "NW_022045319": "chr.29,u1", "NW_022045327": "chr.35,u2", "NW_022045337": "chr.2,u2",
                           "NW_022045339": "chr.1A,u1", "NW_022045341": "chr.32", "NW_022045349": "chr.16,u1",
                           "NW_022045350": "chr.34", "NW_022045351": "chr.33", "NW_022045357": "chr.36", "NW_022045362": "chr.34,u1",
                           "NW_022045374": "chr.35", "NW_022045378": "chr.31,u1", "NW_022045379": "chr.3,u1"}
        unlocalized_d = {"NW_022045287": "5,un1", "NW_022045293": "35,un1", "NW_022045299": "11,un1",
                         "NW_022045306": "2,un1", "NW_022045315": "26,un1", "NW_022045319": "29,un1",
                         "NW_022045327": "35,un2", "NW_022045337": "2,un2", "NW_022045339": "1A,un1",
                         "NW_022045349": "16,un1", "NW_022045362": "34,un1", "NW_022045378": "31,un1",
                         "NW_022045379": "3,un1"}
    elif ref_s == "vgp_calann":
        scaff_2_chrom_d = {"NC_044244": "chr.1","NC_044245": "chr.2","NC_044246": "chr.3","NC_044247": "chr.4","NC_044248": "chr.4A","NC_044249": "chr.4B",
                           "NC_044250": "chr.5","NC_044251": "chr.5A","NC_044252": "chr.6","NC_044253": "chr.7","NC_044254": "chr.8","NC_044255": "chr.9",
                           "NC_044256": "chr.10","NC_044257": "chr.11","NC_044258": "chr.12","NC_044259": "chr.13","NC_044260": "chr.14","NC_044261": "chr.15",
                           "NC_044262": "chr.17","NC_044263": "chr.18","NC_044264": "chr.19","NC_044265": "chr.20","NC_044266": "chr.21","NC_044267": "chr.22",
                           "NC_044268": "chr.23","NC_044269": "chr.24","NC_044270": "chr.25","NC_044271": "chr.26","NC_044272": "chr.27","NC_044273": "chr.28",
                           "NC_044274": "chr.Z","NC_044276": "chr.W","NC_044275": "chr.33"}
        unlocalized_d = {}
    elif ref_s == "vgp_ornana":
        scaff_2_chrom_d = {"NC_041728": "chr.1","NC_041729": "chr.2","NC_041730": "chr.3","NC_041731": "chr.4",
                           "NC_041732": "chr.5","NC_041733": "chr.6","NC_041734": "chr.7","NC_041735": "chr.8",
                           "NC_041736": "chr.9","NC_041737": "chr.10","NC_041738": "chr.11","NC_041739": "chr.12",
                           "NC_041740": "chr.13","NC_041741": "chr.14","NC_041742": "chr.15","NC_041743": "chr.16",
                           "NC_041744": "chr.17","NC_041745": "chr.18","NC_041746": "chr.19","NC_041747": "chr.20",
                           "NC_041748": "chr.21","NC_041749": "chr.X1","NC_041750": "chr.X2","NC_041751": "chr.X3",
                           "NC_041752": "chr.X4","NC_041753": "chr.X5","NW_021638003": "chr.Y1","NW_021638004": "chr.Y2",
                           "NW_021638005": "chr.Y3","NW_021638006": "chr.Y4","NW_021638007": "chr.Y5"}
        unlocalized_d = {}
    else:
        assert ref_s == "vgp_anates"
        scaff_2_chrom_d = {"NC_046610": "chr.1","NC_046611": "chr.2","NC_046612": "chr.3","NC_046613": "chr.4","NC_046614": "chr.5",
                           "NC_046615": "chr.6","NC_046616": "chr.7","NC_046617": "chr.8","NC_046618": "chr.9","NC_046619": "chr.10",
                           "NC_046620": "chr.11","NC_046621": "chr.12","NC_046622": "chr.13","NC_046623": "chr.14","NC_046624": "chr.15",
                           "NC_046625": "chr.16","NC_046626": "chr.17","NC_046627": "chr.18","NC_046628": "chr.19","NC_046629": "chr.21",
                           "NC_046630": "chr.22","NC_046631": "chr.23","NC_046632": "chr.24"}
        unlocalized_d = {}

    try:
        assert (work_s in ["allinone", "missing", "split", "frameshift", "prematurestopcodon", "intronexonjunctiondisruption", "nincodingregion",
                           "summary", "circos_maker", "gene_structure", "false_variant", "missing_chrom"])
    except AssertionError:
        print("Error: %s is not available. please check the help page" % work_s)
        sys.exit()
    print("(3) Submodule [%s] starts..." % work_s)
    if work_s.lower() == "allinone":
        # 1st. CDS information
        a.mp_gff2cds(num_processors_i=process_num_i)
        dir_check(os.path.join(result_dir, "missingGenes"))
        dir_check(os.path.join(result_dir, "missingGenes", "count"))

        # 2nd. totally_missing and exon_deletion
        a.missing_chrom(num_processors_i=process_num_i)
        a.main_missing_genes_analysis(num_processors_i=process_num_i)

        # 3rd. fragmented and intrascaffold split
        dir_check(os.path.join(result_dir, "fragmented"))
        a.main_split_genes_analysis(fgl_s="fragmented")
        dir_check(os.path.join(result_dir, "intrascaffoldsplit"))
        a.main_split_genes_analysis(fgl_s="intrascaffoldsplit")

        # 4th. sequence-level errors

        # frameshift
        dir_check(os.path.join(result_dir, "frameshift"))
        a.mp_frameshift_2_genome(num_processors_i=process_num_i)
        a.mp_main_mpileup_of_sequence_level_fgl(fgl_s="frameshift", num_processors_i=process_num_i, flanking_size_i=2)

        # prematurestopcodon
        dir_check(os.path.join(result_dir, "prematurestopcodon"))
        aln_id_2_premature_stop_codon_variant_d = a.variant_matrix_2_premature_stop_codon_dict()
        a.mp_premature_stop_codon(num_processors_i=process_num_i, aln_id_2_premature_stop_codon_variant_d=aln_id_2_premature_stop_codon_variant_d, flanking_i=5)
        a.mp_main_mpileup_of_sequence_level_fgl(fgl_s="prematurestopcodon", num_processors_i=process_num_i, flanking_size_i=2)

        # splicing junction disruption
        dir_check(os.path.join(result_dir, "intronexonjunctiondisruption"))
        a.splicing_junction(num_processors_i=process_num_i)
        a.mp_main_mpileup_of_sequence_level_fgl(fgl_s="intronexonjunctiondisruption", num_processors_i=process_num_i,flanking_size_i=2)

        # N in coding region
        dir_check(os.path.join(result_dir, "nincodingregion"))
        a.main_n_in_coding_region()

        # 5th. get summary and statistics
        a.count_fgl_number()
        a.ratio_of_fgl()

        # 6th. gene structure analysis
        dir_check(os.path.join(result_dir, "gene_structure"))
        dir_check(os.path.join(result_dir, "gene_structure", "remove"))
        a.gene_structure(mode_s="tgt")
        a.gene_structure(mode_s="ref")
        if ref_s != "vgp_anates":
            a.gene_structure(mode_s="prev")
        cpg_island_prediction()
        a.circos_maker(scaff_2_chrom_d)
        a.mp_summary_mpileup_whole_genome(num_processors_i=process_num_i)

        # 7th. GC and repeat content of genomes
        a.genome_gc_n_repeat_summary(window_size_i=10000)
        a.read_mapping_result_summary()

    if work_s.lower() == "missing_chrom":
        a.missing_chrom(num_processors_i=process_num_i)

    if work_s.lower() == "missing":
        dir_check(os.path.join(result_dir, "missingGenes"))
        a.main_missing_genes_analysis(num_processors_i=process_num_i)
        if ref_s in ["vgp_taegut", "vgp_calann", "vgp_ornana"]:  # climbing perch: No NCBI remap result
            compare_with_remap()
        ref_gene_name_2_intron_d, ref_intron_recs_l = a.generate_intron_records_from_ref_assembly()
        cpg_island_prediction()

    if work_s.lower() == "split":
        dir_check(os.path.join(result_dir, "intrascaffoldsplit"))
        a.main_split_genes_analysis(fgl_s="intrascaffoldsplit")
        dir_check(os.path.join(result_dir, "fragmented"))
        a.main_split_genes_analysis(fgl_s="fragmented")

    if work_s.lower() == "nincodingregion":
        dir_check(os.path.join(result_dir, "nincodingregion"))
        a.main_n_in_coding_region()

    if work_s.lower() == "intronexonjunctiondisruption":
        dir_check(os.path.join(result_dir, "intronexonjunctiondisruption"))
        ref_gene_name_2_intron_d, ref_intron_recs_l = a.generate_intron_records_from_ref_assembly()
        a.splicing_junction(num_processors_i=process_num_i)
        a.mp_main_mpileup_of_sequence_level_fgl(fgl_s="intronexonjunctiondisruption", num_processors_i=process_num_i)

    if work_s.lower() == "frameshift":
        dir_check(os.path.join(result_dir, "frameshift"))
        a.mp_frameshift_2_genome(num_processors_i=process_num_i)
        a.mp_main_mpileup_of_sequence_level_fgl(fgl_s="frameshift", num_processors_i=process_num_i, flanking_size_i=2)

    if work_s.lower() == "prematurestopcodon":
        dir_check(os.path.join(result_dir, "prematurestopcodon"))
        aln_id_2_premature_stop_codon_variant_d = a.variant_matrix_2_premature_stop_codon_dict()
        a.mp_premature_stop_codon(num_processors_i=process_num_i, aln_id_2_premature_stop_codon_variant_d=aln_id_2_premature_stop_codon_variant_d, flanking_i=5)
        a.mp_main_mpileup_of_sequence_level_fgl(fgl_s="prematurestopcodon", num_processors_i=process_num_i)

    if work_s.lower() == "summary":
        a.count_fgl_number()
        a.ratio_of_fgl()
        a.mp_gc_n_repeat_of_genes(num_processors_i=process_num_i)
        a.genome_gc_n_repeat_summary(window_size_i=10000)
        a.circos_maker(scaff_2_chrom_d)
        a.read_mapping_result_summary()

    if work_s.lower() == "circos_maker":
        a.circos_maker(scaff_2_chrom_d)

    if work_s.lower() == "gene_structure":
        dir_check(os.path.join(result_dir, "gene_structure"))
        dir_check(os.path.join(result_dir, "gene_structure", "remove"))
        a.gene_structure(mode_s="ref")
        a.gene_structure(mode_s="tgt")
        if ref_s != "vgp_anates":
            a.gene_structure(mode_s="prev")

    if work_s.lower() == "false_variant":
        a.mp_summary_mpileup_whole_genome(num_processors_i=process_num_i)