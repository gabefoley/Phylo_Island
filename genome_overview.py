from Bio.Graphics import GenomeDiagram
import utilities
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import random
import models


def writeHMMToImage(
    hmm_dict, reference, seq_record, name, query_id, species, expand=False
):
    # Currently taking region to be run simultaneously with hmmer, will need to manipulate some stuff for later
    # Initiation of GenomeDiagram assets"

    # Convert the length of the genome into the length of the translated genome.
    genome_length = round(len(seq_record) / 3)

    # Parameters that can be set

    # The highest number of tracks we'll attempt to add (if a feature overlaps on all tracks it'll just be added to the bottom track)
    # Set max_tracks to None to just add as many as we need to get no overlap
    max_tracks = None

    # The amount we want to allow a feature to overlap with a neighbouring feature
    overlap_amount = 0

    name = reference + "_" + species + "_GenomeDiagram"
    name += "_expanded" if expand else ""

    gd_diagram = GenomeDiagram.Diagram(name)
    max_len = 0
    output_path = "static/user_images/%s_%s%s.png" % (
        query_id,
        species,
        "_expanded" if expand else "",
    )

    output_path = output_path.replace(" ", "_")

    print(output_path)

    print("makea da image")

    start = 0
    # For my work I was considering changing 'region1, 2, and 3' to a3, TcB, and TcC for convenience
    # Up to others though if I fully change that (is just a UI thing tbh)
    region_colours = {
        "A1": "orange",
        "A2": "red",
        "Chitinase": "green",
        "TcdA1": "yellow",
        "TcB": "blue",
        "TcC": "magenta",
        "pore": "grey",
        "region1": "lightblue",
        "region2": "pink",
        "region3": "purple",
        "region4": "black",
    }
    locs = {}
    strand_dict = {}
    strandd = 1
    for result in hmm_dict:
        i = 0
        for reg in result:
            """ Create a dictionary for key = feature type -> value = location """
            if "forward" in reg:
                strandd = 1
                location = reg.split("/")[3] + utilities.randstring(5)
                locs[location] = result[reg].split(":")
                strand_dict[location] = strandd

            elif "backward" in reg:
                strandd = -1
                print("bananas")
                print(reg)
                location = reg.split("/")[3] + utilities.randstring(5)
                locs[location] = result[reg].split(":")
                strand_dict[location] = strandd

            i += 1
            # Prepare for literally the worst code in existence
            """ We have to pull the features from location data in database
            instead of directly from database because of current limitations """
            # Brute force if statements to pull from database
            # TODO - Could likely turn this into a for loop in so

    """ Add the features from dictionary to the seq_record features """
    feature_locations = []
    # If want to add all features to image comment below line out
    seq_record.features = []
    i = 0
    for location in locs:

        print("whoopsy")
        print(location)
        """ Extract start and end values from each location, and add to independent lists """
        """ create and add features based on locations """
        feature = SeqFeature(
            location=FeatureLocation(
                int(locs[location][0]),
                int(locs[location][1]),
                strand=strand_dict[location],
            ),
            type=location[0:-5],
        )
        seq_record.features.append(feature)

    """ Set up the Genome Diagram """
    max_len = max(max_len, len(seq_record))

    gd_track1 = gd_diagram.new_track(
        0, name=name + " Track 1", greytrack=True, start=0, end=len(seq_record)
    )
    gd_feature_set1 = gd_track1.new_set()

    for feature in seq_record.features:
        feature_locations.append((feature.location.start, feature.location.end))

    total_tracks = 1
    # Dictionary to keep track of which locations are at which track
    forward_tracks = {1: []}
    backward_tracks = {1: []}

    for feature in seq_record.features:
        feature_added = False
        current_track = 1

        while current_track <= total_tracks and not feature_added:

            current_strand = feature.location.strand

            if current_strand == 1:
                dict_track = forward_tracks
            elif current_strand == -1:
                dict_track = backward_tracks

            overlap = False

            for loc in dict_track[current_track]:

                if (
                    feature.location.start + overlap_amount in loc
                    or feature.location.end - overlap_amount in loc
                ):
                    overlap = True

            if overlap:

                # If we've reached the highest track, make a new track
                if current_track == total_tracks:
                    if max_tracks is not None and current_track == max_tracks:

                        # We've reached the highest track we want to try and add, so just add it to the bottom track (there will be overlap)
                        gd_feature_set1.add_feature(
                            feature,
                            label=True,
                            name=feature.type,
                            color=region_colours[feature.type],
                            label_position="middle",
                        )
                        dict_track[1].append(feature.location)
                        feature_added = True
                    else:
                        # Add to next highest track
                        current_track += 1
                        total_tracks += 1
                        if current_track not in dict_track.keys():
                            dict_track[current_track] = []
                        if current_track not in forward_tracks.keys():
                            forward_tracks[current_track] = []
                        if current_track not in backward_tracks.keys():
                            backward_tracks[current_track] = []
                        exec(
                            "gd_track"
                            + str(current_track)
                            + "= gd_diagram.new_track(0, name=name + ' Track "
                            + str(current_track)
                            + "', greytrack=True, start=0, end=len(seq_record))"
                        )
                        exec(
                            "gd_feature_set"
                            + str(current_track)
                            + " = gd_track"
                            + str(current_track)
                            + ".new_set()"
                        )

                else:
                    # Check on higher track
                    current_track += 1

            else:

                print("howie")
                print(feature)
                print(feature.type)
                # Add to this track
                exec(
                    "gd_feature_set"
                    + str(current_track)
                    + ".add_feature(feature, label=True, name=feature.type, color=region_colours[feature.type], label_position='middle')"
                )
                if current_track in dict_track:
                    dict_track[current_track].append(feature.location)
                else:
                    dict_track[current_track] = [feature.location]
                feature_added = True

    for track in range(1, total_tracks):
        exec("gd_track" + str(track) + ".add_set(gd_feature_set" + str(track) + ")")
    gd_diagram.draw(
        format="linear", pagesize="A2", fragments=10, start=start, end=len(seq_record)
    )
    gd_diagram.write(output_path, "PNG")
    print("Genome Diagram has been added to file " + output_path)
    fh = open(output_path, "rb")
    return fh


# def add_hit_to_database(query_id, species, region, start, end, expand):
#
#     print ('coco checking ' + str(species) + str(region) + str(start) + str(end))
#
#     models.GenomeRecords.objects(_id= query_id, hits__region=region, hits__start=start, hits__end= end)
#
#     Feed.objects(_id="...", posts__region=region, hits__start).update(set__posts__S__value="updatevalue")
#
#     if models.Hits.objects().get(region=region, start=start, end=end):
#         hit = models.Hits(region=region, score=str(random.randrange(1, 11, 1)), start=str(start), end=str(end),
#                          expand=expand, tags=[])
#
#         hit.save()
#     else:
#         print ('it was already there')
#


def write_hits_to_gb(hmm_dict, reference, seqrecord, query_id, species, expand=False):

    print("writing to genbank")
    print("and reference is ")
    print(reference)
    name = species + "_sequence"
    name += "_expanded" if expand else ""
    output_path = reference + "/" + name + ".gb"

    output_path = output_path.replace(" ", "_")
    print(name)
    print(reference)
    print("and species name is", species)
    seqrecord.name = species[0:9].zfill(9).replace(" ", "_")

    print(seqrecord.name)

    # Write Annotated Sequences to Genbank files to allow easy movement to Artemis
    print("Writing sequences to GenBank File")
    """ Create a dictionary for key = feature type -> value = location """
    locs = {}
    strand_dict = {}
    colour_dict = {
        "A1": "255 165 0",
        "A2": "255 0 0",
        "TcdA1": "255 255 0",
        "TcB": "0 0 255",
        "TcC": "255 0 255",
        "Chitinase": "0 255 0",
        "region1": "0 255 255",
        "region2": "255 153 255",
        "region3": "204 0 102",
        "region4": "0 0 0",
    }
    for result in hmm_dict:
        print(result)
        i = 0
        for reg in result:
            """ Create a dictionary for key = feature type -> value = location """
            if "forward" in reg:
                location = reg.split("/")[3] + utilities.randstring(5)
                locs[location] = result[reg].split(":")
                strandd = 1
                strand_dict[location] = strandd

            elif "backward" in reg:

                strandd = -1

                location = reg.split("/")[3] + utilities.randstring(5)
                locs[location] = result[reg].split(":")
                strand_dict[location] = strandd

            i += 1

    print("Adding %s genome to diagram" % (seqrecord.name))

    seqrecord.features = []
    for location in locs:
        """ create and add features based on locations """
        color = {"color": colour_dict[location[0:-5]]}
        feature = SeqFeature(
            location=FeatureLocation(
                int(locs[location][0]),
                int(locs[location][1]),
                strand=strand_dict[location],
            ),
            type=location[0:-5],
            qualifiers=color,
        )
        seqrecord.features.append(feature)

        # add_hit_to_database(query_id = query_id, species=species, region=location[0:-5], start=int(locs[location][0]), \
        #                                                                             end=int(locs[
        #                                                                                                       location][
        #                                                                                                       1]),
        #                     expand=expand)

    # sequence = str(seqrecord)[2:-1]
    # seqrecord.seq = Seq(sequence, generic_dna)
    SeqIO.write(seqrecord, output_path, "genbank")


def classify_genomes_original(queries):
    for query in queries:
        region_names = set([hit.region for hit in query.hits])
        print(f"\nClassifying {query.description}")
        print(f"It has the following regions {region_names}")
        if set(["TcB", "TcC"]).issubset(region_names):
            if set(["A1", "A2"]).issubset(region_names):
                if set(["Chitinase"]).issubset(region_names):
                    print(
                        "It had A1 / A2 and chitinases, so we're tagging it as Type 2A"
                    )
                    update_tag(query, "Auto_Type2A")
                else:
                    print(
                        "It had A1 / A2 but no chitinases, so we're tagging it as Type 2B"
                    )
                    update_tag(query, "Auto_Type2B")

            elif set(["TcdA1"]).issubset(region_names):
                if set("Chitinase").issubset(region_names):
                    print(
                        "It had TcdA1 but also chitinase, so we're tagging it for further investigation"
                    )
                    update_tag(query, "Auto_Unsure")

                else:
                    print("It had TcdA1 so we're classifying it as Type 1")
                    update_tag(query, "Auto_Type1")

            else:
                print(
                    "It had a TcB and TcC, but was missing either A1 or TcdA1, so we're classifying it as incomplete"
                )
                update_tag(query, "Auto_Incomplete")

        else:
            print("It was lacking a TcB and TcC, so we're classifying it as incomplete")
            update_tag(query, "Auto_Incomplete")


def classify_genomes(queries):

    for query in queries:
        region_names = set(
            [
                hit.region
                for hit in query.hits
                if "expanded" in hit.region
                if "hidden" not in hit.tags
            ]
        )

        print(f"\nClassifying {query.description}")
        print(f"It has the following regions {region_names}")

        print("Check if it should be Single or Multiple")

        allow_multi = ["Chitinase_expanded", "TcB_expanded", "TcC_expanded"]
        multi_check = [
            hit.region
            for hit in query.hits
            if "expanded" in hit.region and hit.region not in allow_multi
        ]

        print(multi_check)

        if len(multi_check) != len(set(multi_check)):
            count_check = "Multiple"
        else:
            count_check = "Single"

        # Count the number of chitinases
        chitinase_count = 0
        for hit in query.hits:
            if hit.region == "Chitinase_expanded":
                chitinase_count += 1

        # Check if there is at least one A1 and A2 that don't overlap

        a1_starts = []
        a1_ends = []
        a2_starts = []
        a2_ends = []

        if set(["A1_expanded", "A2_expanded"]).issubset(region_names):
            for hit in query.hits:
                if hit.region == "A1_expanded":
                    a1_starts.append(hit.start)
                    a1_ends.append(hit.end)

                if hit.region == "A2_expanded":
                    a2_starts.append(hit.start)
                    a2_ends.append(hit.end)

            print("Here are the A1 start positions ")
            print(a1_starts)

            print("Here are the A1 end positions")
            print(a1_ends)

            print("Here are the A2 start positions ")
            print(a2_starts)

            print("Here are the A2 end positions")
            print(a2_ends)

            start_unique = False
            end_unique = False

            for a1_start in a1_starts:
                if a1_start not in a2_starts:
                    start_unique = True

            for a2_start in a2_starts:
                if a2_start not in a1_starts:
                    start_unique = True

            for a1_end in a1_ends:
                if a1_end not in a2_ends:
                    end_unique = True

            for a2_end in a2_ends:
                if a2_end not in a1_ends:
                    end_unique = True

            if start_unique and end_unique:
                print("There is a non-overlapping A1 and A2 region")

        if set(["TcB_expanded", "TcC_expanded"]).issubset(region_names):
            if (
                set(["A1_expanded", "A2_expanded"]).issubset(region_names)
                and start_unique
                and end_unique
            ):

                # Needs to have more than one chitinase
                if (
                    set(["Chitinase_expanded"]).issubset(region_names)
                    and chitinase_count > 1
                ):
                    print(
                        "It had A1 / A2 and chitinases, so we're tagging it as Type 2A"
                    )
                    update_tag(query, count_check, "Type2a")

                    genome_tag = models.GenomeTags(
                        tag_id=query.name, tags=[count_check, "Type2a"]
                    )
                    genome_tag.save()
                else:
                    print(
                        "It had A1 / A2 but no chitinases, so we're tagging it as Type 2B"
                    )
                    update_tag(query, count_check, "Type2b")
                    genome_tag = models.GenomeTags(
                        tag_id=query.name, tags=[count_check, "Type2b"]
                    )
                    genome_tag.save()

            elif set(["TcdA1_expanded"]).issubset(region_names) or set(
                ["A2_expanded"]
            ).issubset(region_names):
                if set("Chitinase_expanded").issubset(region_names):
                    print(
                        "It had TcdA1 or A2 but also chitinase, so we're tagging it for further investigation"
                    )
                    update_tag(query, count_check, "TcdA1_w_Chitinase")
                    genome_tag = models.GenomeTags(
                        tag_id=query.name, tags=[count_check, "TcdA1_w_Chitinase"]
                    )
                    genome_tag.save()

                else:

                    for hit in query.hits:
                        if hit.region == "TcB_expanded":
                            tcB_start = hit.start
                        elif hit.region == "TcC_expanded":
                            tcC_start = hit.start

                    if tcB_start == tcC_start:
                        print(
                            "It had TcdA1 or A2 and a fused TcB / TcC so we're classifying it as "
                            "Type 3"
                        )
                        update_tag(query, count_check, "Type3")
                        genome_tag = models.GenomeTags(
                            tag_id=query.name, tags=[count_check, "Type3"]
                        )
                        genome_tag.save()
                    else:

                        print(
                            "It had TcdA1 or A2 and a split TcB / TcC so we're classifying it as Type 1"
                        )
                        update_tag(query, count_check, "Type1")
                        genome_tag = models.GenomeTags(
                            tag_id=query.name, tags=[count_check, "Type1"]
                        )
                        genome_tag.save()

            else:
                print(
                    "It had a TcB and TcC, but was missing either A1 or TcdA1, so we're classifying it as incomplete"
                )
                update_tag(query, count_check, "Incomplete")
                genome_tag = models.GenomeTags(
                    tag_id=query.name, tags=[count_check, "Incomplete"]
                )
                genome_tag.save()

        else:
            print("It was lacking a TcB or TcC, so we're classifying it as incomplete")
            update_tag(query, count_check, "Incomplete")
            genome_tag = models.GenomeTags(
                tag_id=query.name, tags=[count_check, "Incomplete"]
            )
            genome_tag.save()

        region_names = [hit.region for hit in query.hits if "expanded" in hit.region]


def delete_genome_tags(queries):
    for query in queries:
        query.update(tags=[])


def update_tag(query, *args):
    for tag in args:
        if tag not in query.tags:
            query.update(push__tags=tag)
