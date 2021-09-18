import argparse
import os
import mongoengine
import cmd_code
import utilities
import models
import genome_overview
from flask import Flask
from flask_mongoengine import MongoEngine
from bson.objectid import ObjectId
import getGenomes
import gzip
import refseq_code
import numpy
from collections import defaultdict
import alignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import time
import pandas as pd
import glob

from Bio.Alphabet import generic_nucleotide


parser = argparse.ArgumentParser()

parser.add_argument("-g", "--add_genomes", help="path to list of species")
# parser.add_argument("-i", "--input_file", help="path to list of species")
parser.add_argument("-p", "--add_profiles", help="path to profile folder")

parser.add_argument(
    "-u", "--update_genomes", help="update genomes", action="store_true"
)
parser.add_argument(
    "-f", "--fasta", help="save all regions to fasta files", action="store_true"
)
parser.add_argument(
    "-r",
    "--region_order",
    help="write out order of regions in all genomes ",
    action="store_true",
)

parser.add_argument(
    "-o", "--overview", help="get overview of database", action="store_true"
)
parser.add_argument(
    "-dg", "--delete_genomes", help="delete genomes", action="store_true"
)
parser.add_argument(
    "-dt",
    "--delete_genome_tags",
    help="delete current genome classifications",
    action="store_true",
)
parser.add_argument(
    "-c",
    "--classify",
    help="classify genomes based on their regions ",
    action="store_true",
)
parser.add_argument(
    "-m", "--mlgo", help="write out regions to mlgo format", action="store_true"
)
parser.add_argument("--refseq", help="run full refseq check", action="store_true")
parser.add_argument("--query_db", help="query db", action="store_true")
parser.add_argument("--load_genomes", help="load genomes")

parser.add_argument("--delete_all_profiles", action="store_true")

parser.add_argument

args = parser.parse_args()

if args.add_genomes:
    # Retrieve the genomes
    print("Adding genomes")
    cmd_code.get_genomes(args)

if args.add_profiles:
    # Delete the existing profiles if there is a new one in the Profile folder

    print("Adding profiles")
    cmd_code.delete_profiles(args)

    # Update the profiles
    cmd_code.update_profiles(args)

if args.update_genomes:

    print("Updating genomes")
    cmd_code.update_genomes()

if args.overview:
    print("Overview of database")
    cmd_code.get_overview()

if args.delete_genomes:
    queries = models.GenomeRecords.objects.all().timeout(False).delete()
    print("Deleting all genomes")

if args.delete_genome_tags:
    queries = models.GenomeRecords.objects.all().timeout(False)
    print("Deleting all existing genome classification tags")
    genome_overview.delete_genome_tags(queries)

if args.classify:
    queries = models.GenomeRecords.objects.all().timeout(False)
    print("Classifying the genomes")
    genome_overview.classify_genomes(queries)


if args.fasta:
    profile_names = models.Profile.objects().all()

    for profile in profile_names:
        getGenomes.download_fasta_regions(
            profile.name + "_expanded", filename="cmd", split_strands=False, align=False
        )

if args.region_order:
    genomes = models.GenomeRecords.objects.all().timeout(False)
    getGenomes.write_genome_order(
        genomes, split_strands=False, path="./fasta_folder/genome_order_from_cmd.txt"
    )

if args.mlgo:
    genomes = models.GenomeRecords.objects.all().timeout(False)
    getGenomes.write_mlgo_order(genomes, path="./fasta_folder/mlgo.txt")


if args.query_db:
    genomes = models.GenomeRecords.objects(species="Halomicronema hongdechloris")

    # genomes = models.GenomeRecords.objects(species='Photorhabdus heterorhabditis')

    for genome in genomes:
        print(genome)
        for hit in genome.hits:
            print(hit)


if args.load_genomes:

    if not args.load_genomes[-1] == "/":
        filepath = args.load_genomes + "/"
    else:
        filepath = args.load_genomes

    # filepath = "./files/test_genomes/"
    # filepath = "/Users/gabefoley/Dropbox/PhD/Projects/Phylo_Island/2020/20200817_Checking_full_9000_genomes/genbank" \
    #            "/bacteria/"

    # filepath = "/Users/gabefoley/Dropbox/PhD/Projects/Phylo_Island/2020/20200831_Getting_just_matches/genbank/bacteria/"
    # filepath = "/Users/gabefoley/Dropbox/PhD/Projects/Phylo_Island/2020/20200831_Getting_just_matches/refseq/bacteria/"

    # Just the YP
    # filepath = "/Users/gabefoley/Dropbox/PhD/Projects/Phylo_Island/2020/20201028_Testing_yersinia_pseudotuberculosis" \
    #            "/refseq/bacteria/"
    #
    # # Just the YP
    # filepath = "/Users/gabefoley/Dropbox/PhD/Projects/Phylo_Island/2020/20201119_Testing_Feature_Tables" \
    #            "/refseq/bacteria/"

    # Octopus - eight genomes four each from genbank / refseq with different conditions for testing feature tables
    # filepath = '/Users/gabefoley/Dropbox/PhD/Projects/Phylo_Island/2020/20201215_Checking_all_types_reading_features' \
    #            '/genbank/bacteria/'

    genome_name = [x for x in os.listdir(filepath) if x != ".DS_Store"]

    chunk = numpy.array_split(numpy.array(genome_name), 5)

    print(chunk)

    for genomes in chunk:

        print(genomes)

        for genome in genomes:

            print("NEW GENOME")

            print(filepath)
            print(genome)

            genome_path = glob.glob(filepath + genome + "/*_genomic.fna.gz")[0]

            report_path = glob.glob(filepath + genome + "/*_assembly_report.txt")[0]

            print(genome_path)

            print(report_path)

            # genome_path = filepath + genome + "/" + genome + '_genomic.fna.gz'

            file_from_zip = gzip.open(genome_path, mode="rb")

            outpath = ".".join(genome_path.split(".")[0:-1]) + ".fasta"

            with open(outpath, "w") as query_file:
                for line in file_from_zip:
                    query_file.write(line.decode("utf-8"))

            file_from_zip.close()

            genome_list = []

            genome_dict = refseq_code.get_genome_dict(outpath, report_path)

            # genome_dict = refseq_code.read_genome(outpath, genome)

            if len(genome_dict) > 0:
                print("Found something!")

            print(genome_dict)

            #     for k, v in genome_dict.items():
            #         print (k)
            #         print (v)
            #
            #     print ("Adding genome\n")
            #
            #     print (genome_dict)
            #
            utilities.add_genome(genome_dict)

            # print ("Remove unzipped FASTA file from disk\n")
            # utilities.remove_file(outpath)
            #
            #
            # print ("Remove zip files from disk\n")
            #
            # if os.path.exists(genome_path):
            #
            #     utilities.remove_file(genome_path)
        #
        #
        # print ("Search for profile in new genomes\n")
        #
        # queries = models.GenomeRecords.objects(hits__size=0).timeout(False)
        #
        # profiles = models.Profile.objects(name='TcB_BD_NCBI')
        #
        # for profile in profiles:
        #     cmd_code.get_feature_location_with_profile_cmd(queries, "hmm_outputs", profile.name, "", "", profile.name)
        #
        # del profiles
        # del queries
        #
        # print ("Remove genomes without a hit\n")
        # missing_hit = models.GenomeRecords.objects(hits__size=0)
        #
        # with open(filepath + "removed_genomes.txt", "a") as removed_genomes:
        #     for genome in missing_hit:
        #         removed_genomes.write(genome.description + "\n")
        #
        #
        #
        # models.GenomeRecords.objects(hits__size=0).delete()


if args.delete_all_profiles:
    # Delete the profiles from the database (so we can easily update them with new ones if need be)
    profiles_to_delete = models.Profile.objects()
    profiles_to_delete.delete()
