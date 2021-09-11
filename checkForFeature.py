import models
import phyloisland
import models
import utilities
import os
import time
from flask import flash
import subprocess
import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import SeqIO, Alphabet
import phyloisland
import glob
import resultread
import genome_overview
import random

def get_feature_location_with_profile(ids, reference, profile_name, recordName, recordLocation, region):
    """
    Annotate a genome sequence with a feature location based on a profile
    :param ids: Genome sequences to annotate
    :param reference: Profile to annotate based on
    :param recordName: Which feature field to update
    :param recordLocation: Which feature location field to update
    :return:
    """

    queries = models.GenomeRecords.objects(id__in=ids).timeout(False)

    for query in queries:

        seq_record = query.sequence

        # Create a path to write the translated genomic sequence to
        # random_id = utilities.randstring(5)

        # Get the nucleotide sequence of the genome
        nuc_seq = Bio.Seq.Seq(str(query.sequence))

        seq_record = SeqRecord(Seq(query.sequence, generic_dna), '', '', '')

        # seq_record = SeqRecord(seq=Seq(nuc_seq, Alphabet()), id='', name='', description='', dbxrefs=[])

        # outpath = reference + "/" + query.name + "/" + query.species + "/" + region + "/"
        outpath = reference + "/" + query.name + "/" + query.species + "/"

        outpath = outpath.replace(" ", "_")

        print (outpath)

        # Check three forward reading frames
        if not os.path.exists(outpath):
            os.makedirs(outpath.replace(" ", "_"))

        for forward in [True, False]:
            for i in list(range(0, 3)):

                strand = "_forward_" + str(i) if forward else "_backward_" + str(i)
                sequence = nuc_seq[i:] if forward else nuc_seq.reverse_complement()[i:]

                cleaned_path = outpath + query.name.replace("/",
                                                               "_") + strand + "_translated_genome.fasta"
                hmmsearch_results = outpath + region + "/" + query.name.replace("/",
                                                                    "_") + strand + "_hmmsearch_results.fasta"
                domScore = 100

                cleaned_path = cleaned_path.replace(" ", "_")
                hmmsearch_results = hmmsearch_results.replace(" ", "_")

                if not os.path.isfile(cleaned_path):

                    print ('Making a new file')

                    # Translate the nucleotide genome sequence to a protein sequence
                    with open(cleaned_path, 'w') as handle:

                        if query.name == "<unknown name>":
                            handle.write(">" + query.description + "\n" + str(
                                sequence.translate(stop_symbol="*")))
                        else:

                            handle.write(">" + query.description + "\n" + str(
                                sequence.translate(stop_symbol="*")))

                    print("Writing the %s sequence with the species %s to %s" % (
                    query.name, query.species, cleaned_path))

                    while not os.path.exists(cleaned_path):
                        time.sleep(1)



                if os.path.isfile(cleaned_path):

                    print (cleaned_path)
                    print (hmmsearch_results)

                    if not os.path.exists(outpath + "/" + region):
                        os.mkdir(outpath + "/" + region)

                    print ('No need to make a new file')
                    stdoutdata = subprocess.getoutput("hmmsearch -o %s --domT %s %s %s" % (
                    hmmsearch_results, domScore, 'tmp/' + region + "_profile.hmm", cleaned_path))

                    print ('woodflock')

                    print (hmmsearch_results)

                    print (domScore)

                    print (cleaned_path)

                    print(stdoutdata)
                    # result = subprocess.call(["hmmsearch -o %s %s %s" % (hmmsearch_results, reference, cleaned_path)])

                    print("The results from the HMM search have been written to %s \n" % hmmsearch_results)
                    # read_hmmer_results(hmmsearch_results)
                    # result = subprocess.call(["hmmsearch", 'files/output.txt', reference, cleaned_path], stdout=subprocess.PIPE)
                    # for x in result:
                    #     print (x)
        # for regions in hmmer_outputs/organism add reg to all_reg
        all_reg = []
        for infile in glob.glob(os.path.join(reference + "/" + query.name + "/" + query.species.replace(" ",
                                                                                                      "_") + '/*')):
            print ('infile')
            print (infile)
            if '.gb' not in infile:
                all_reg.append(infile)
        hmmerout = []
        hmmerout_expanded = []
        print ('here is all reg')
        print (all_reg)
        # add handler to HMMread for output paths
        for reg in all_reg:
            hmmerout.append(resultread.HMMread(reg, query))
            hmmerout_expanded.append(resultread.HMMread(reg, query, expand=True))


            # Update the Genome Overview graphic and save it to the Genome Record

        print ('here are dict regions')
        print (hmmerout)
        print (hmmerout_expanded)

        print ('should we write to gb?')


        # genome_image = genome_overview.writeHMMToImage(hmmerout, reference + "/" + query.name +
        #                                                "/" + query.species.replace(" ", "_"), nuc_seq, query.name, \
        #                                                    query.id, query.species)
        #
        # genome_expanded_image = genome_overview.writeHMMToImage(hmmerout_expanded, reference + "/" + query.name +
        #                                                                                                  "/" +
        #                                                         query.species.replace(" ", "_"), nuc_seq, query.name, query.id, query.species, expand=True)
        #
        # genbank = genome_overview.write_hits_to_gb(hmmerout, reference +"/" + query.name + "/" +
        #                                            query.species.replace(
        #     " ",
        #                                                                                                        "_"),
        #                                            seq_record, query.id, query.species)
        #
        # genbank_expanded = genome_overview.write_hits_to_gb(hmmerout_expanded, reference +"/" + query.name + "/" +
        #                                                     query.species.replace(" ", "_"), seq_record,
        #                                                     query.id, query.species, expand=True)


        # curr = models.GenomeRecords.objects().timeout(False).get(id=query.id)
        #
        # curr.genome_overview.replace(genome_image)
        # curr.genome_expanded_overview.replace(genome_expanded_image)

        #


    del queries



# Function to be called when wanting to generate Genome Diagram and GenBank output
def generate_output(ids, reference, diagram=True, genbank=True, interactive_diagram=False, fasta=True, expand=False):
    query = models.GenomeRecords.query.filter(models.GenomeRecords.uid.in_(ids))
    fasta_region_dict = {}
    for record in query.all():
        seq_record = servers.bio_db.lookup(primary_id=record.name)
        species = record.name.replace(" ", "_")
        species = species.replace(".", "_")

        # for regions in hmmer_outputs/organism add reg to all_reg

        all_reg = []

        for infile in glob.glob(os.path.join(reference + "/" + species + '/*')):
            if '.' not in infile:
                all_reg.append(infile)

        outpath = os.path.join(reference + "/" + species.replace(" ", "_"))

        hmmerout = []

        for reg in all_reg:

            hmmerout.append(resultread.HMMread(reg, record, expand))

            if fasta:
                reg_name = reg.split("/")[-1]

                if reg_name not in fasta_region_dict:
                    fasta_region_dict[reg_name] = []

                # For a given region, create the sequence files
                fasta_output = utilities.createFASTAFromHMMOutput(seq_record, hmmerout, record.name, reg_name)

                # Write the region sequence files out to a genome specific FASTA file
                fasta_outpath = "%s/%s%s.fasta" % (outpath, reg_name, "_expanded" if expand else "")
                utilities.saveFASTA(fasta_output, fasta_outpath)

                # Add the region to a dictionary so we can write out a non-genome specific FASTA file
                fasta_region_dict[reg_name] = fasta_region_dict[reg_name] + fasta_output

                hmmerout = []

        if diagram:
            ToxinGraphicsMain.writeHMMToImage(hmmerout, outpath, seq_record,
                                              species, expand)
            print("Diagram has been written to %s directory" % (reference))

        if genbank:
            # TODO Add better invariant for Shotgun Naming
            ToxinGraphicsMain.writeHmmToSeq(hmmerout, outpath, seq_record, species, expand)

        if interactive_diagram:
            print("output is")
            print(hmmerout)

            tracks = utilities.get_tracks(hmmerout)

    # Write out the non-genome specific FASTA file/s
    if fasta:
        for reg_name, output in fasta_region_dict.items():
            outpath = "%s/hmm_outputs/%s_%s%s.fasta" % (
            os.getcwd(), reg_name, "expanded_" if expand else "", str(len(ids)))
            utilities.saveFASTA(output, outpath)
            print ('porcupine')
            print (reg_name)
            print (outpath)

