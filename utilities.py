import models
import os
from flask import flash
import random
import time
import subprocess
import sys
from collections import defaultdict
import json
import regex as re
import shutil
import glob
import tree_code
from Bio.Seq import Seq
from Bio.Alphabet import generic_nucleotide
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, AlignIO, SearchIO
from Bio.Align.Applications import MuscleCommandline
import pickle
import ete3
from configs.auto_classifier import original_classifications
# from genome_overview import models

import getGenomes

region_name_mapper = {
    "A1": "track1",
    "A2": "track2",
    "Chitinase": "track3",
    "TcdA1": "track4",
    "TcB": "track5",
    "TcC": "track6",
    "region1": "track7",
    "A1_expanded": "track1_expanded",
    "A2_expanded": "track2",
    "Chitinase_expanded": "track3_expanded",
    "TcdA1_expanded": "track4_expanded",
    "TcB_expanded": "track5_expanded",
    "TcC_expanded": "track6_expanded",
    "region1_expanded": "track7_expanded",
    "EXISTING:": "track8",
}


def read_fasta(filename):
    """
    Read in a FASTA file
    :param filename:
    :return: Dictionary object containing a SeqRecord
    """
    return SeqIO.to_dict(SeqIO.parse(filename, "fasta"))


def readLinesFromFile(filepath):
    """
    Takes a file and reads each individual line into a set
    :param filepath: Path of the file
    :return: Set containing lines from the file
    """

    content = set()

    with open(filepath, "r") as query_file:
        for line in query_file:
            if len(line) > 1:
                content.add(line.strip())
    return content


def remove_file(*args):
    """
    Remove files in the list from the directory

    :param args: Files to remove
    :return:
    """
    for arg in args:
        os.remove(arg)


def remove_folder(*args):
    for arg in args:
        shutil.rmtree(arg)


def add_genome(genome_results):
    """
    Add a genome into the database
    :param genome_results:
    """
    for record in genome_results:

        current = genome_results[record]
        if type(current) == SeqRecord:
            name = current.id
            species = " ".join(current.annotations.get("organism").split()[0:2])
            plasmid = current.annotations["plasmid"]
            organism = current.annotations["organism"]
            assembly_name = current.annotations["assembly_name"]
            biosample = current.annotations["biosample"]
            bioproject = current.annotations["bioproject"]
            date = current.annotations["date"]
            wgs_project = current.annotations["wgs_project"]
            genome_coverage = current.annotations["genome_coverage"]

            taxid = current.annotations["taxid"]
            assembly_type = current.annotations["assembly_type"]
            release_type = current.annotations["release_type"]
            assembly_level = current.annotations["assembly_level"]
            genome_representation = current.annotations["genome_representation"]
            expected_final_version = current.annotations["expected_final_version"]
            excluded = current.annotations["excluded"]
            genbank_accession_id = current.annotations["genbank_accession_id"]
            refseq_accession_id = current.annotations["refseq_accession_id"]
            r_g_identical = current.annotations["r_g_identical"]
            sequence = str(current.seq)
            description = current.description

            # Check to see if the genome record already exists
            if models.GenomeRecords.objects(name=name):
                print(
                    "The genome record - %s from species - %s already exists in the database"
                    % (name, species)
                )
                continue

            else:
                print(
                    "Adding the genome record - %s from species - %s to the genome database"
                    % (name, species)
                )

                genome = models.GenomeRecords(
                    name=name,
                    species=species,
                    organism=organism,
                    assembly_name=assembly_name,
                    biosample=biosample,
                    bioproject=bioproject,
                    date=date,
                    wgs_project=wgs_project,
                    genome_coverage=genome_coverage,
                    taxid=taxid,
                    assembly_type=assembly_type,
                    release_type=release_type,
                    assembly_level=assembly_level,
                    genome_representation=genome_representation,
                    expected_final_version=expected_final_version,
                    excluded=excluded,
                    genbank_accession_id=genbank_accession_id,
                    refseq_accession_id=refseq_accession_id,
                    r_g_identical=r_g_identical,
                    plasmid=plasmid,
                    description=description,
                    sequence=sequence,
                    tags=[],
                )
                genome.save()


def addSequence(seq_records):
    """
    Add a sequence into the database
    :param seq_records:
    """
    for record in seq_records.values():
        seq_name = record.id
        seq_description = record.description.split(">")[0]
        seq_species = seq_description.split("[")[1].split("]")[0]
        seq_sequence = str(record.seq)

        # Check if the sequence record already exists
        if models.SequenceRecords.objects(name=seq_name):
            print(
                "Sequence with ID - %s from species - %s already exists in the sequence database"
                % (seq_name, seq_species)
                + "\n"
            )
            flash(
                "Sequence with ID - %s from species - %s already exists in the sequence database"
                % (seq_name, seq_species)
                + "\n"
            )

        else:
            print(
                "Adding sequence with ID - %s from species - %s to the sequence database"
                % (seq_name, seq_species)
                + "\n"
            )

            sequence = models.SequenceRecords(
                name=seq_name,
                species=seq_species,
                description=seq_description,
                sequence=seq_sequence,
            )
            sequence.save()


def save_profile(profile, name=None):
    """
    Save a profile into the database
    :param profile: Profile to save
    """
    print("Saving profile")
    if not name:
        name = randstring(5)

    profile_entry = models.Profile(name, profile, {})
    profile_entry.save()


def set_profile_as_reference(profile_ids, region):
    """

    :param ids:
    :param region:
    :return:
    """
    if len(profile_ids) > 1:
        flash("Only select a single record", category="error")
    else:

        # Check for a previous Profile set as this reference
        prev = models.Profile.objects(
            __raw__={"references.%s" % (region): {"$exists": True}}
        )

        if prev:
            prev.update(**{"unset__references__%s" % (region): "*"})

        profile_id = profile_ids[0]

        curr = models.Profile.objects().get(id=profile_id)

        # curr.update(references= {region: "*"})

        curr.update(**{"set__references__%s" % (region): "*"})

        curr.save()

        # Write the new profile to the tmp folder ready to be used
        with open("tmp/" + region + "_profile.hmm", "wb") as profile_path:

            profile_path.write(curr.profile.read())

            # flash("The profile named %s has been set as the reference profile for %s" % (curr.name, region), category='success')


def createAlignment(input, output):

    subprocess.getoutput(f"muscle -in {input} -out {output}")


def search_regions_with_profiles(region_to_search, profile_ids):
    domScore = 1
    region_dict = {}

    regions = models.RegionRecords.objects.get(name=region_to_search)

    profile_folder = f"./tmp/profiles_{regions.name}/"

    os.mkdir(profile_folder)

    fasta_path = f"{profile_folder}{regions.name}.fasta"

    with open(fasta_path, "w+") as fasta_file:
        fasta_file.write(regions.regions.decode())

    while not os.path.exists(fasta_path):
        time.sleep(1)

    if os.path.isfile(fasta_path):

        seqs = read_fasta(fasta_path)

        for name in seqs.keys():
            region_dict[name.replace(".", "***")] = {}

    for profile_id in profile_ids:
        profile = models.Profile.objects().get(id=profile_id)

        with open(profile_folder + profile.name + "_profile.hmm", "wb") as profile_path:
            profile_path.write(profile.profile.read())

        while not os.path.exists(profile_folder + profile.name + "_profile.hmm"):
            time.sleep(1)

        subprocess.getoutput(
            f"hmmsearch -o {profile_folder}{profile.name}.output --domT {domScore} {profile_folder}{profile.name}_profile.hmm {fasta_path}"
        )

        while not os.path.exists(profile_folder + profile.name + ".output"):
            time.sleep(1)

    region_dict, domain_dict = process_hmmer_results(region_dict, profile_folder)

    # Delete the profile folder we used
    remove_folder(profile_folder)

    return region_dict, domain_dict


def process_hmmer_results(region_dict, profile_folder):
    domain_dict = {}
    for infile in glob.glob(profile_folder + "*.output"):

        qresult = SearchIO.read(infile, "hmmer3-text")
        if len(qresult.hits) > 0:

            for hsp in qresult.hsps:
                # hsp = qresult[0][i]

                print(hsp.query.id)
                start = hsp.hit_start
                end = hsp.hit_end
                print("pinto")
                print(qresult.id)
                print(hsp.hit_id)
                print(start)
                print(end)
                if hsp.query.id in region_dict[hsp.hit_id.replace(".", "***")]:
                    region_dict[hsp.hit_id.replace(".", "***")][
                        hsp.query.id + "_multiple_" + randstring(5)
                    ] = (start, end)

                else:
                    region_dict[hsp.hit_id.replace(".", "***")][hsp.query.id] = (
                        start,
                        end,
                    )

                if hsp.query.id not in domain_dict:
                    domain_dict[hsp.query.id] = []
                domain_dict[hsp.query.id].append(hsp.hit_id)

    print(region_dict)

    print(domain_dict)

    return region_dict, domain_dict


def make_alignment_from_regions(aln_name, region_data, tool="MAFFT"):
    # Write out region data to a FASTA file

    fasta_path = "./tmp/" + aln_name + ".fasta"
    aln_path = "./tmp/" + aln_name + ".aln"

    with open(fasta_path, "w+") as fasta_file:
        fasta_file.write(region_data)

    while not os.path.exists(fasta_path):
        time.sleep(1)

    if os.path.isfile(fasta_path):
        if tool == "MAFFT":
            print("Aligning with MAFFT")
            stdoutdata = subprocess.getoutput(
                f"mafft  --reorder {fasta_path} > {aln_path}"
            )

    if os.path.isfile(aln_path):
        return aln_path


def load_list(*args):
    return_list = []
    for x in args:
        with open(x) as f:
            curr = [line.strip() for line in f]
        return_list += curr
    return return_list


def get_sequence_content_dict(region):
    fasta_path = "./tmp/tmp_regions.fasta"

    alns = models.AlignmentRecords.objects().get(name=region)

    # Write out the regions

    with open(fasta_path, "w+") as fasta_file:
        fasta_file.write(alns.alignment.read().decode())

    while not os.path.exists(fasta_path):
        time.sleep(1)

    # Load the regions in as a dictionary

    if os.path.isfile(fasta_path):
        seqs = read_fasta(fasta_path)

    seq_content_dict = {k: v.seq for k, v in seqs.items()}
    print(seq_content_dict)
    return seq_content_dict


def make_tree(alignment_name, alignment, tool):
    aln_path = "./tmp/" + alignment_name + ".aln"
    tree_path = "./tmp/" + alignment_name + ".nwk"

    with open(aln_path, "w+") as fasta_file:
        fasta_file.write(alignment)

    while not os.path.exists(aln_path):
        time.sleep(1)

    if os.path.isfile(aln_path):
        if tool == "FastTree":
            print("Making tree with FastTree")
            stdoutdata = subprocess.getoutput(
                f"fasttree -nosupport {aln_path} > {tree_path}"
            )

    if os.path.isfile(tree_path):
        return tree_path


def get_tree_image(
    tree,
    tree_name,
    tag_dict,
    region_dict,
    region_order_dict,
    sequence_content_dict,
    colour_dict,
    full_names,
    collapse_on_genome_tags,
    display_circular,
    display_circular_180,
):
    tree_path = f"./tmp/{tree_name}.nwk"
    img_path = f"static/img/trees/{tree_name}{'_full' if full_names else ''}{'_collapse' if collapse_on_genome_tags else ''}{'_rd' if region_dict else ''}{'_ro' if region_order_dict else ''}{'_sc' if sequence_content_dict else ''}{'_circ' if display_circular else ''}{'_circ180' if display_circular_180 else ''}.png"
    tag_dict_path = f"./tmp/{tree_name}_tagdict.p"
    region_dict_path = f"./tmp/{tree_name}_regiondict.p"
    region_order_dict_path = f"./tmp/{tree_name}_regionorderdict.p"
    sequence_content_dict_path = f"./tmp/{tree_name}_seqcontentdict.p"
    colour_dict_path = f"./tmp/{tree_name}_colourdict.p"

    print(img_path)

    # Write out tree to file
    with open(tree_path, "w+") as tree_file:
        tree_file.write(tree)

    while not os.path.exists(tree_path):
        time.sleep(1)

    pickle_dict(tag_dict, tag_dict_path)
    pickle_dict(region_dict, region_dict_path)
    pickle_dict(region_order_dict, region_order_dict_path)
    pickle_dict(sequence_content_dict, sequence_content_dict_path)

    print("here is the colour dict")
    print(colour_dict)

    pickle_dict(colour_dict, colour_dict_path)

    if os.path.isfile(tree_path):

        # remove_file(img_path)

        print(tree_path)

        loaded_tree = tree_code.load_tree(tree_path)

        stdoutdata = subprocess.getoutput(
            f'python tree_code.py -t {tree_path} -o {img_path} -td {tag_dict_path} -rd {region_dict_path} -rod {region_order_dict_path} -scd {sequence_content_dict_path} -cd {colour_dict_path} -fn {full_names} -cgt {collapse_on_genome_tags} {" -dc" if display_circular else ""} {" -dco" if display_circular_180 else ""}'
        )

        print(stdoutdata)

        # tree_code.colour_tips(loaded_tree, tag_dict, colour_dict, region_dict, outpath=img_path,
        #                                         custom_layout=False)

        if os.path.isfile(img_path):
            return img_path


def highlight_taxonomy(tree):
    ts = ete3.TreeStyle()
    # disable default PhyloTree Layout
    ts.layout_fn = lambda x: True

    for n in tree.traverse():

        if not n.is_leaf():
            N = ete3.TextFace(
                n.sci_name + " (" + n.rank + ")", fsize=14, fgcolor="black"
            )
            n.add_face(N, 12, position="branch-top")

            nstyle = ete3.NodeStyle()
            nstyle["shape"] = "sphere"
            nstyle["fgcolor"] = "blue"
            nstyle["size"] = 10
            nstyle["hz_line_type"] = 1
            n.set_style(nstyle)

        else:
            pass

    return tree, ts


def get_species_tree_image(tree, tree_name):
    img_path = f"static/img/trees/{tree_name}.png"
    ncbi_img_path = f"static/img/trees/{tree_name}_ncbi.png"

    stree = ete3.PhyloTree(
        tree, sp_naming_function=lambda name: name.split("_taxid_")[1].split("_")[0]
    )

    # Create annotated species tree image
    tax2names, tax2lineages, tax2rank = stree.annotate_ncbi_taxa()
    tcb, ts = highlight_taxonomy(stree)
    stree.render(img_path, tree_style=ts, dpi=300)

    # Create NCBI species tree image
    ncbi = ete3.NCBITaxa()
    taxids = get_taxids(stree)
    ncbi_tree = ncbi.get_topology(taxids)
    ncbi_tree, ts = highlight_taxonomy(ncbi_tree)
    ncbi_tree.render(ncbi_img_path, tree_style=ts, dpi=300)

    if os.path.isfile(img_path):
        return img_path


# Get out a list of the taxonomic IDs stored in the _taxids_ annotation on a tree
def get_taxids(tree):
    taxids = []
    for n in tree.traverse():
        taxids.append(n.taxid)
    return taxids


def get_ml_go_tree_image(tree, name, ancestral_orders, ref_ml_go_dict):
    tree_path = f"./tmp/{name}_ml.nwk"
    img_path = f"static/img/trees/{name}_ml.png"
    ancestral_orders_path = f"./tmp/{name}_ancestralorders.p"
    ref_ml_go_dict_path = f"./tmp/{name}_ref_ml_go_dict.p"

    # Write out tree to file
    with open(tree_path, "w+") as tree_file:
        tree_file.write(tree)

    while not os.path.exists(tree_path):
        time.sleep(1)

    pickle_dict(ancestral_orders, ancestral_orders_path)
    pickle_dict(ref_ml_go_dict, ref_ml_go_dict_path)

    if os.path.isfile(tree_path):

        print("call it")
        print(tree_path)
        print(img_path)
        print(ancestral_orders_path)
        stdoutdata = subprocess.getoutput(
            f"python tree_code.py -t {tree_path} -o {img_path} -ao {ancestral_orders_path} -mlgo {ref_ml_go_dict_path}"
        )

        print(stdoutdata)

        if os.path.isfile(img_path):
            return img_path


def create_pos_dict(regions, profile, trimmed_name, fasta_path):
    profile_path = "./tmp/" + trimmed_name + ".hmm"
    results_outpath = "./tmp/" + trimmed_name + "_results.txt"

    # Write out regions

    with open(fasta_path, "w+") as regions_out:
        regions_out.write(regions.decode())
    # Write out profile

    with open(profile_path, "wb") as profile_out:
        profile_out.write(profile.read())

    # Perform the hmmsearch on the regions file
    os.system(
        "hmmsearch -o"
        + results_outpath
        + " --domT 1 "
        + profile_path
        + " "
        + fasta_path
    )

    seqs = read_fasta(fasta_path)

    pos_dict = get_pos_dict_from_hmm(results_outpath, seqs)

    return pos_dict


def get_pos_dict_from_hmm(path, seqs):
    """
    Read in a hmm output file and extract the positions of the hits for a given set of sequences
    :param path: path of the hmm output file
    :param seqs: SeqRecord of the sequences we want to search for
    :return: A dictionary mapping sequence name -> (start position, end position)
    """
    qresult = SearchIO.read(path, "hmmer3-text")

    pos_dict = {}

    print(len(qresult.hsps))
    print(len(seqs))

    if len(qresult.hsps) > len(seqs):
        print("ERROR: More HSPs than expected")

    for hsp in qresult.hsps:
        pos_dict[hsp.hit.id] = (hsp.hit_start, hsp.hit_end)

    return pos_dict


def trim_to_profile(regions, profile, trimmed_name):
    fasta_path = "./tmp/" + trimmed_name + ".fasta"
    trimmed_path = "./tmp/" + trimmed_name + "_trimmed"

    pos_dict = create_pos_dict(regions, profile, trimmed_name, fasta_path)

    seqs = read_fasta(fasta_path)

    trimmed = []

    failed_seqs = []

    for name, seq in seqs.items():
        if name in pos_dict:
            trimmed_seq = seq.seq.tomutable()
            trimmed_seq = trimmed_seq[int(pos_dict[name][0]) : int(pos_dict[name][1])]
            trimmed.append(SeqRecord(trimmed_seq, seq.name))

            print(trimmed_seq)
        else:
            failed_seqs.append(name)

    # Write the newly trimmed file to disk
    createFasta(trimmed, trimmed_path, align=False)

    # Load and save the newly trimmed file
    with open(trimmed_path + ".fasta", "rb") as trimmed_seqs:
        # Load files into database
        region = models.RegionRecords(name=trimmed_name, regions=trimmed_seqs.read())
        region.save()

    # Return sequences that failed to have a hmmer match
    return failed_seqs


def trim_around_profile(regions, profile1, profile2, pos1, pos2, trimmed_name):
    fasta_path = "./tmp/" + trimmed_name + ".fasta"
    trimmed_path = "./tmp/" + trimmed_name + "_trimmed"

    if profile1:
        pos_dict1 = create_pos_dict(regions, profile1, trimmed_name, fasta_path)
    else:
        pos_dict1 = None

    if profile2:
        pos_dict2 = create_pos_dict(regions, profile2, trimmed_name, fasta_path)
    else:
        pos_dict2 = None

    seqs = read_fasta(fasta_path)
    trimmed = []
    failed_seqs = []

    trimmed, failed_seqs = trim_sequence(seqs, pos_dict1, pos_dict2, pos1, pos2)

    # Write the newly trimmed file to disk
    createFasta(trimmed, trimmed_path, align=False)

    # Load and save the newly trimmed file
    with open(trimmed_path + ".fasta", "rb") as trimmed_seqs:
        # Load files into database
        region = models.RegionRecords(name=trimmed_name, regions=trimmed_seqs.read())
        region.save()

    # Return sequences that failed to have a hmmer match
    return failed_seqs


def trim_sequence(seqs, pos_dict1=None, pos_dict2=None, pos1="start", pos2="start"):
    trimmed = []
    failed_seqs = []

    for name, seq in seqs.items():

        if pos_dict1 == None:
            if pos_dict2 == None:
                raise NameError("Must provide at least one position dictionary")

            else:  # From the start of a sequence to a profile match

                if name in pos_dict2:

                    print("From the start of sequence to profile match")

                    print(pos2)

                    trimmed_seq = seq.seq.tomutable()

                    trimmed_seq = trimmed_seq[
                        : int(pos_dict2[name][0 if pos2 == "start" else 1])
                    ]

                    if trimmed_seq:
                        trimmed.append(SeqRecord(trimmed_seq, seq.name))
                    else:
                        failed_seqs.append(name)

                else:
                    failed_seqs.append(name)

        elif pos_dict2 == None:  # From a profile match to the end of a sequence
            if name in pos_dict1:
                trimmed_seq = seq.seq.tomutable()

                print("From a profile match to the end of the sequence")

                print(pos1)

                trimmed_seq = trimmed_seq[
                    int(pos_dict1[name][0 if pos1 == "start" else 1]) :
                ]
                if trimmed_seq:
                    trimmed.append(SeqRecord(trimmed_seq, seq.name))
                else:
                    failed_seqs.append(name)
            else:
                failed_seqs.append(name)

        else:  # Between two profile matches
            if name in pos_dict1 and name in pos_dict2:

                print("Between two sequences")

                print(pos1)

                print(pos2)
                trimmed_seq = seq.seq.tomutable()

                trimmed_seq = trimmed_seq[
                    int(pos_dict1[name][0 if pos1 == "start" else 1]) : int(
                        pos_dict2[name][0 if pos2 == "start" else 1]
                    )
                ]
                if trimmed_seq:
                    trimmed.append(SeqRecord(trimmed_seq, seq.name))
                else:
                    failed_seqs.append(name)
            else:
                failed_seqs.append(name)

    return trimmed, failed_seqs


def search_for_promoters(mismatch):
    genomes = models.GenomeRecords.objects()

    regions = []

    prom_regex = "(TTGACA.{15,25}TATAAT){s<=" + str(mismatch) + "}"

    for g in genomes:
        # if g.name == "NZ_CP010029.1":
        for hit in g.hits:
            if "expanded" in hit.region:
                print(hit.region)
                print(hit.start)
                print(hit.end)
                if hit.strand == "forward":
                    seq_content = g.sequence[int(hit.start) - 50 : int(hit.start)]

                    print(seq_content)
                    regions.append(seq_content)
                    match = re.findall(prom_regex, seq_content)
                    print(match)
                    if match:
                        hit.promoter = True
                    else:
                        hit.promoter = False

                elif hit.strand == "backward":
                    seq_content = g.sequence[int(hit.end) : int(hit.end) + 50]

                    rev_content = str(Seq(str(seq_content)).reverse_complement())
                    print(rev_content)
                    regions.append(rev_content)
                    match = re.findall(prom_regex, rev_content)
                    print(match)
                    if match:
                        hit.promoter = True
                    else:
                        hit.promoter = False

        g.save()


def clear_all_promoters():
    """
    Clear all the promoters for all the hits
    :return:
    """
    genomes = models.GenomeRecords.objects()

    for g in genomes:
        for hit in g.hits:
            hit.promoter = False
        g.save()


def pickle_dict(dict, outpath):
    with open(outpath, "wb") as handle:
        pickle.dump(dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


def createProfile(align_list):
    SeqIO.write(align_list, "tmp/align.fasta", "fasta")

    hmm_path = "tmp/profile3.hmm"

    createAlignment("tmp/align.fasta", "tmp/align.aln")

    outfile = open(hmm_path, "w")
    result = subprocess.call(
        ["hmmbuild", hmm_path, "tmp/align.aln"], stdout=subprocess.PIPE
    )

    while not os.path.exists(hmm_path):
        time.sleep(1)

    if os.path.isfile(hmm_path):
        file = open(hmm_path, "rb")

        save_profile(file)
        remove_file(hmm_path, "tmp/align.fasta", "tmp/align.aln")


def createFasta(fasta_list, name, align):
    SeqIO.write(fasta_list, name + ".fasta", "fasta")

    while not os.path.exists(name + ".fasta"):
        time.sleep(1)

    if align:
        # print("And now an aln")
        createAlignment(name + ".fasta", name + ".aln")


def check_with_profile(ids, region):
    # Check if a reference profile for this region exists
    profile_reference = models.Profile.objects(
        __raw__={"references.%s" % (region): {"$exists": True}}
    )
    if profile_reference:
        for profile in profile_reference:
            print(
                "Using the %s profile named %s to check for %s regions"
                % (region, profile.name, region)
            )

            eval(
                'checkForFeature.get_feature_location_with_profile(ids, "hmm_outputs'
                + '", "'
                + profile.name
                + '", "'
                + region
                + '", "'
                + region
                + "_loc"
                + '","'
                + region
                + '")'
            )
    else:
        flash(
            "Please set a profile as the %s reference profile first" % (region), "error"
        )


def get_genome_items(
    genome,
    hits="all",
    hidden_type=True,
    show_promoters=False,
    show_stop_codons=False,
    show_existing_features=False,
    checked_regions=None,
):
    """
    Format the items in a genome correctly for a Genome Detail view
    :param self:
    :return:
    """

    items = defaultdict(list)
    stop_codon_glyphs = {}
    promoter_glyphs = {}
    hit_tags = {}
    region_list = []
    genomesize = len(genome.sequence)

    for count, hit in enumerate(genome.hits):

        if (
            (hits == "all")
            or ((hits == "initial") and ("expanded" not in hit.region))
            or ((hits == "expanded") and "expanded" in hit.region)
            or (show_existing_features and "EXISTING" in hit.region)
        ):

            if (
                (checked_regions == None)
                or hit.region in checked_regions
                or hit.region.split("_expanded")[0] in checked_regions
                or (show_existing_features and "EXISTING" in hit.region)
            ):

                if (hidden_type == False) or (
                    hidden_type == True and "hidden" not in hit.tags
                ):

                    # Here is a place to update FASTA ID headers

                    if hit.tags:
                        hit_tags[str(hit.id)] = (
                            str(hit.region)
                            + " ["
                            + str(hit.score)
                            + "] "
                            + str(hit.start)
                            + ":"
                            + str(hit.end)
                            + " "
                            + hit.strand,
                            hit.tags,
                        )

                    if hit.name != None:

                        hit_details = dict()
                        hit_details["id"] = count
                        hit_details["hit_id"] = str(hit.id)
                        hit_details["start"] = hit.start
                        hit_details["end"] = hit.end
                        hit_details["name"] = hit.region
                        hit_details["score"] = hit.score
                        # hit_details['strand'] = 1 if count % 2 == 0 else -1

                        # We set a fake strand here just to force TcC and TcdA1 to a different track
                        if (
                            hit.region == "TcdA1"
                            or hit.region == "TcdA1_expanded"
                            or hit.region == "TcC"
                            or hit.region == "TcC_expanded"
                            or hit.region.startswith("EXISTING:")
                        ):
                            hit_details["strand"] = -1
                        else:
                            hit_details["strand"] = 1

                        # But store the real strand here so we can annotate the hits correctly
                        hit_details["actual_strand"] = hit.strand

                        # Add the promoters:

                        if show_promoters:

                            if hit.promoter:

                                promoter_pos = (
                                    hit.start if hit.strand == "forward" else hit.end
                                )

                                if promoter_pos in promoter_glyphs:

                                    promoter_glyphs[promoter_pos].append(
                                        hit.region + " promoter"
                                    )

                                else:
                                    promoter_glyphs[promoter_pos] = [
                                        hit.region + " promoter"
                                    ]

                        idx1 = 0
                        idx2 = idx1 + 3

                        stop_codons = ["TAG", "TAA", "TGA"]

                        if hit.strand == "backward":
                            sequence = Seq(hit.sequence, generic_nucleotide)

                            hit_sequence = sequence.reverse_complement()
                        else:
                            hit_sequence = hit.sequence
                        #
                        # print ('flipped seq')
                        #
                        # print (hit_sequence)
                        # print ('get the sequence')

                        # print (hit_sequence)
                        #
                        # print (hit.strand)
                        #
                        # print (hit.start)
                        #
                        # print (hit.end)

                        while idx2 <= len(hit.sequence):
                            if hit_sequence[idx1:idx2] in stop_codons:
                                # print('found', idx1)
                                # print (hit_sequence)
                                # print (hit.region)
                                # print (hit_sequence[idx1:idx2 + 20])
                                #
                                # print (hit.start)
                                # print (hit.end)
                                # print (idx1)

                                if hit.strand == "backward":
                                    pos = int(hit.end) - idx1
                                else:
                                    pos = int(hit.start) + idx1

                                # print (pos)

                                if show_stop_codons:

                                    if pos in stop_codons:
                                        stop_codon_glyphs[pos].append(hit.region)
                                    else:
                                        stop_codon_glyphs[pos] = [hit.region]

                            idx1 += 3
                            idx2 += 3

                        region_list.append(hit_details)

                        items[hit.region].append(hit_details)

    # print (items)

    print("stop_codons")
    print(stop_codon_glyphs)
    tracks = build_tracks(
        items, stop_codon_glyphs, promoter_glyphs, show_existing_features
    )

    return tracks, hit_tags, genomesize


# def get_hit_tags((hits, hits='all', hidden_type=True, checked_regions=None):):


def build_tracks(items, stop_codons, promoters, show_existing_features):
    tracks = []

    for region_name in items:

        regions = []

        # NOTE: I'm forcing the strands all to be 1 here to visualise on the same line in the linear genome

        # Here is a place to update FASTA ID headers

        for region in items[region_name]:
            region_dict = {
                "id": region["hit_id"],
                "start": int(region["start"]),
                "end": int(region["end"]),
                "name": region["name"] + " [" + region["score"] + "]",
                "strand": region["strand"],
                "actual_strand": region["actual_strand"],
            }
            regions.append(region_dict)

            # if show_existing_features:
            #     print ('show features')
            #     print (region_name)

            # if show_existing_features and region_name.startswith('EXISTING'):
            #     print('&&&')
            #     print(region_name)
            region_name = region_name.split(" ")[0]

        track = {
            "trackName": region_name_mapper[region_name],
            "trackType": "stranded",
            "visible": "true",
            "inner_radius": 130,
            "outer_radius": 185,
            "trackFeatures": "complex",
            "featureThreshold": 7000000,
            "mouseclick": "linearClick",
            "mouseover_callback": "islandPopup",
            "mouseout_callback": "islandPopupClear",
            "linear_mouseclick": "linearPopup",
            "showLabels": "true",
            "showTooltip": "true",
            "linear_mouseclick": "linearClick",
            "items": regions,
        }

        tracks.append(track)

    if stop_codons:

        stop_codon_regions = []

        count = 0

        for loc, names in stop_codons.items():
            count += 1
            name = " ".join(names)
            stop_codon_dict = {"id": count, "bp": loc, "type": name, "name": name}

            stop_codon_regions.append(stop_codon_dict)

        glyph_track = {
            "trackName": "track1",
            "trackType": "glyph",
            "glyphType": "circle",
            "radius": 155,
            "pixel_spacing": 5,
            "linear_pixel_spacing": 5,
            "glyph_buffer": 5,
            "linear_glyph_buffer": 5,
            "glyphSize": 20,
            "linear_glyphSize": 20,
            "linear_height": 100,
            "showTooltip": "true",
            "items": stop_codon_regions,
        }

        tracks.append(glyph_track)

    if promoters:

        promoter_regions = []

        count = 0

        for loc, names in promoters.items():
            count += 1
            name = names
            promoter_dict = {"id": count, "bp": loc, "type": "promoter", "name": name}

            promoter_regions.append(promoter_dict)

        promoter_track = {
            "trackName": "track3",
            "trackType": "glyph",
            "glyphType": "diamond",
            "radius": 155,
            "pixel_spacing": 5,
            "linear_pixel_spacing": 5,
            "glyph_buffer": 5,
            "linear_glyph_buffer": 5,
            "glyphSize": 40,
            "linear_glyphSize": 40,
            "linear_height": 100,
            "showTooltip": "true",
            "items": promoter_regions,
        }

        tracks.append(promoter_track)

    json_tracks = json.dumps(tracks)

    return json_tracks


def get_associated_dict(genome):
    associated_dict = {}

    print("assoc")

    print(genome.id)

    # assoc_hits = models.AssociatedHits.objects().get(genome_id=str(genome.id))

    aggregate = models.AssociatedHits._get_collection().aggregate(
        [{"$match": {"genome_id": str(genome.id)}}]
    )

    assoc_hits = list(aggregate)

    # for as in aggregate:

    print(assoc_hits)

    for hit in assoc_hits:
        print(hit["_id"])
        print(hit["region1"])
        associated_dict[str(hit["_id"])] = (
            hit["region1"].split("region_")[1]
            + " and "
            + hit["region2"].split("region_")[1]
        )

    return associated_dict


def open_file(filename):
    print("file name is ")
    print(filename)
    file_path = "static/uploads/" + filename

    print("file path is")
    print(file_path)
    while not os.path.exists(file_path):
        time.sleep(1)
    if os.path.isfile(file_path):
        file = open(file_path, "rb")
        return file


def randstring(length=10):
    valid_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    return "".join((random.choice(valid_letters) for i in range(length)))


def sort_func(elem):
    return int(elem.split("_position=_")[1].split("_")[0])


def rename_duplicates(genome_name, old):
    seen = {}
    index = {}
    for pos in range(len(old)):
        check = old[pos]

        if check in seen:
            seen[check] += 1
            if check in index:
                index.pop(check)
        else:
            seen[check] = 1
            index[check] = pos
        old[pos] = check + "_" + str(seen[check])

    print(index)

    # Either add the genome name, or remove the temporarily annotated number from the end

    for pos in range(len(old)):
        if pos in index.values():
            old[pos] = "_".join(old[pos].split("_")[0:-1])
        else:
            old[pos] = (
                "_".join(old[pos].split("_")[0:-1])
                + "_"
                + genome_name
                + "_"
                + old[pos].split("_")[-1]
            )

    return old


def test_auto_classify(queries, skip_tags):


    test_results = ""

    count = 0
    diff_count = 0
    for query in queries:
        if query.name in original_classifications:
            count += 1

            new_classification = models.GenomeTags.objects.get(tag_id=query.name)

            set1 = set(original_classifications[query.name])

            check = list(
                set(new_classification.tags).difference(
                    set(original_classifications[query.name])
                )
            ) + list(
                set(original_classifications[query.name]).difference(
                    set(new_classification.tags)
                )
            )

            check_skip = [x for x in check if x not in skip_tags]

            if check_skip:
                diff_count += 1
                print(
                    "\n\nFound an automatically classified genome that was different\n"
                )
                print(query.name)
                print("Original was ")
                print(
                    [
                        x
                        for x in original_classifications[query.name]
                        if x not in skip_tags
                    ]
                )
                print("Automatic classification was ")
                print(new_classification.tags)

                test_results += "\n\nFound an automatically classified genome that was different\n\n"
                test_results += query.name + "\n"
                test_results += "Original was \n"
                test_results += str(
                    [
                        x
                        for x in original_classifications[query.name]
                        if x not in skip_tags
                    ]
                )
                test_results += "\nAutomatic classification was \n"
                test_results += str(
                    [x for x in new_classification.tags if x not in skip_tags]
                )

                # if original_classifications[query.name] != new_classification.tags[0].split("Auto_")[1]:

    print("\nWrong: " + str(diff_count))
    print("Correct " + str(count - diff_count))
    print("Total " + str(count))
    print()

    test_results = "Total " + str(count) + "\n" + test_results

    test_results = "Correct " + str(count - diff_count) + "\n" + test_results

    test_results = "\nWrong: " + str(diff_count) + "\n" + test_results

    return test_results


def delete_all_tags():

    queries = models.GenomeRecords.objects().all()

    queries.update(tags=[])

    for query in queries:

        # By default we keep the tag 'hidden' as this won't get GenomeTags out of sync

        for hit in query.hits:
            if "hidden" in hit.tags:
                hit.tags = ["hidden"]
            else:
                hit.tags = []

        query.save()

    models.GenomeTags.objects().all().delete()


def get_mlgo_dict(gene_orders):
    mlgo_dict = {}

    lines = gene_orders.split("\n")
    for i in range(0, len(lines)):
        line = lines[i]
        if line.startswith(">"):
            mlgo_dict[line.split(">")[1].strip()] = (
                lines[i + 1].split("$")[0].strip().split(" ")
            )

    print(mlgo_dict)

    return mlgo_dict


def colour_alignment_by_profiles(alignment, profiles):
    colour_dict = {
        "RBD_A": "lightgreen",
        "RBDA": "lightgreen",
        "RBD_C": "blue",
        "RBD_B": "orange",
        "Neuraminidase": "mediumPurple",
        "TcA_RBD": "grey",
        "TcB_BD_seed": "lawnGreen",
        "PF18276_ncbi": "lawnGreen",
        "VRP1_Full": "sandyBrown",
        "Big_1_Full": "lightYellow",
        "TcB_BLAST_500": "orange",
        "TcC_BLAST_500": "blue",
        "Rhs_repeat": "green",
        "Overlap": "pink",
    }

    split = [x for x in alignment.split(">") if len(x) > 0]

    print(split[0:5])

    size = len(split)

    output = {
        split[i]
        .replace(" <unknown description", "")
        .replace(".", "***"): split[i + 1]
        .replace("\n", "")
        if size > i + 1
        else None
        for i in range(0, size, 2)
    }

    for seqname, domains in profiles.region_dict.items():

        print(seqname)
        print(output.keys())

        if seqname not in output:
            print("wowzers")
        # if seqname in profiles.region_dict:

        len_offset = 0
        furtherst_pos = -1

        # newlist = sorted(list_to_be_sorted, key=lambda k: k['name'])

        if seqname in output:
            orig_seq = output[seqname]

            sorted_list = sorted(
                profiles.region_dict[seqname].items(), key=lambda k: k[1]
            )

            # print(sorted_list)

            list_w_overlaps = []

            overlap = False

            for idx, entry in enumerate(sorted_list):
                domain = entry[0]
                pos = entry[1]

                if idx + 1 < len(sorted_list):

                    next_entry = sorted_list[idx + 1]
                    next_domain = next_entry[0]
                    next_pos = next_entry[1]

                    if pos[1] > next_pos[0]:

                        overlap = True
                        print("WARNING: OVERLAP")
                        overlap_pos_1 = next_pos[0]
                        overlap_pos_2 = pos[1]

                        prev_entry = (domain, [pos[0], next_pos[0]])
                        overlap_entry = ("Overlap", [overlap_pos_1, overlap_pos_2])
                        sorted_list[idx + 1] = (
                            next_domain,
                            [overlap_pos_2, next_pos[1]],
                        )

                        list_w_overlaps.append(prev_entry)
                        list_w_overlaps.append(overlap_entry)

                    else:
                        list_w_overlaps.append(entry)

                else:
                    list_w_overlaps.append(entry)

            if overlap:

                print("list with overlaps")
                print(seqname.replace("***", "."))
                print(list_w_overlaps)

            for entry in list_w_overlaps:

                if entry[1][1] < entry[1][0]:
                    print("WARNING: Serious error with overlap")
                    print(seqname.replace("***", "."))
                    print(list_w_overlaps)

            for domain, pos in list_w_overlaps:
                gap_offset = 0
                first_gap_offset = 0
                second_gap_offset = 0

                count = 0
                for aa in orig_seq:
                    if aa == "-":
                        gap_offset += 1
                    else:
                        count += 1
                        if count == pos[0] + 1:
                            first_gap_offset = gap_offset

                        if count == pos[1]:
                            second_gap_offset = gap_offset
                            break
                        else:
                            print("count was " + str(count))

                prev_len = len(output[seqname])

                output[seqname] = (
                    output[seqname][0 : pos[0] + len_offset + first_gap_offset]
                    + '<span style = "background-color:'
                    + colour_dict[domain.split("_multiple_")[0]]
                    + '">'
                    + output[seqname][
                        pos[0]
                        + len_offset
                        + first_gap_offset : pos[1]
                        + len_offset
                        + second_gap_offset
                    ]
                    + "</span>"
                    + output[seqname][pos[1] + len_offset + second_gap_offset :]
                )

                len_offset = len(output[seqname]) - len(orig_seq)

    return output
