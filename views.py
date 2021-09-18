from phyloisland import app, allfiles
import models
import forms
import genome_overview
import utilities
import getGenomes
import custom_filters
from flask import (
    render_template,
    flash,
    request,
    session,
    send_file,
    has_app_context,
    redirect,
    url_for,
)
from flask_login import login_required, current_user, user_logged_in
from flask_admin import Admin, expose, BaseView
from flask_admin.actions import action
from flask_admin.contrib.mongoengine import ModelView, filters
from Bio.Alphabet import generic_protein
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from markupsafe import Markup
import gettext
import timeit
import time
import os
from urllib.error import HTTPError
import random
import io
import json
import Bio
import warnings
import math
from ete3 import PhyloTree
from wtforms import SelectField
import mongoengine

ref_names = [
    "A1",
    "A2",
    "Chitinase",
    "TcB",
    "TcC",
    "TcdA1",
    "region1",
    "region2",
    "region3",
    "region4",
    "existing_feature",
]

ref_mlgo_dict = {
    "A1": "1",
    "A2": "2",
    "Chitinase": "3",
    "TcB": "4",
    "TcC": "5",
    "TcdA1": "6",
    "Fused_TcB_TcC": "7",
}
mlgo_ref_dict = {
    "1": "A1",
    "2": "A2",
    "3": "Chitinase",
    "4": "TcB",
    "5": "TcC",
    "6": "TcdA1",
    "7": "Fused_TcB_TcC",
}


@user_logged_in.connect_via(app)
def on_user_logged_in(sender, user):
    # Clear any existing session values
    keys = [
        key
        for key in session.keys()
        if key
        not in ["_fresh", "_permanent", "csrf_token", "user_id", "_user_id", "_id"]
    ]

    print(keys)

    for key in keys:
        session.pop(key)

    # Get the current User
    current = models.User.objects().get(username=str(current_user.username))

    session["page_size"] = current.page_size if current.page_size != None else 20
    session["record_size"] = current.record_size if current.record_size != None else 20
    # session['genome'] = None


class UploadView(BaseView):
    """
    View for uploading files
    """

    @login_required
    @expose("/", methods=("GET", "POST"))
    def upload(self):
        form = forms.UploadForm()
        print("is it validated - ", form.validate())

        if request.method == "POST" and not form.validate():
            print("not validated")
            print(form.genome_type)
            print(request.get_json())

        if request.method == "POST" and form.validate():

            if request.method == "POST":
                names = request.get_json()
                print("names coming in")
                print(request)
                print(request.data)
                print(names)
                # for name in names:
                #     print (name)

            start_time = timeit.default_timer()

            # Get the information from the upload form
            filename = allfiles.save(request.files["file"])
            seq_type = form.type.data
            add_seq = form.add_sequence.data
            add_genome = form.add_genome.data
            single = form.single_genome.data
            genome_type = form.genome_type.data
            representative = form.representative.data
            assembly = form.assembly.data
            genbank = form.genbank.data
            failed_genomes = []

            print("seq type is ", seq_type)

            try:

                if seq_type == "protein" or seq_type == "nucleotide":
                    print("here")
                    seq_records = utilities.read_fasta("static/uploads/" + filename)

                    if not seq_records:
                        print("Couldn't find any sequences in the uploaded file.")
                        flash(
                            "Couldn't find any sequences in the uploaded file.",
                            category="error",
                        )

                    else:
                        if add_seq:
                            utilities.addSequence(seq_records)

                        if add_genome:
                            # TODO: Correct this call
                            getGenomes.add_genome(seq_records, seq_type)

                elif seq_type == "species":

                    species_names = utilities.readLinesFromFile(
                        "static/uploads/" + filename
                    )
                    for species_name in species_names:

                        print("Species name is ", species_name)

                        destinations = [genome_type]

                        # TODO: Once the checkbox selection is dynamic we can just add freely here
                        if representative and genome_type not in [
                            "representative genome",
                            "assembly",
                            "genbank",
                        ]:
                            destinations.append("representative genome")

                        if assembly and genome_type not in ["assembly", "genbank"]:
                            destinations.append("assembly")

                        print("Destinations is ", destinations)

                        genome_results = getGenomes.add_genome(
                            species_name, destinations, single=single
                        )

                        if genome_results and genome_results != "Fail":
                            utilities.add_genome(genome_results)
                        else:

                            # Couldn't find it in RefSeq, let's try genbank
                            if genbank:
                                destinations = ["genbank"]
                                genome_results = getGenomes.add_genome(
                                    species_name, destinations, single=single
                                )

                                if genome_results and genome_results != "Fail":

                                    utilities.add_genome(genome_results)
                                else:
                                    failed_genomes.append(species_name)

                            else:
                                failed_genomes.append(species_name)

                elif seq_type == "genome":
                    pass
                    # genome_names = readLinesFromFile("static/uploads/" + filename)
                    #
                    # for name in genome_names:
                    #     genome_ids = {}
                    #     genome_results = {}
                    #     genome_ids = readLinesFromFile("static/uploads/" + filename)
                    #     genome_results = mapToGenome.get_full_genome([name])
                    #     if genome_results and genome_results != 'in_database':
                    #         addGenome(genome_results)
                    #     elif genome_results != 'in_database' and search_shotgun:
                    #         print("\nWe didn't identify any genome records for %s. Attempting to search for shotgun "
                    #               "sequenced genomes \n" % (
                    #                 name))
                    #
                    #         shotgun_id_dict = mapToGenome.get_shotgun_id_dict_from_id(name)
                    #         if shotgun_id_dict:
                    #             # print("Skipping getting shotgun")
                    #             # shotgun_results = None
                    #             genome_results = mapToGenome.get_shotgun_genome(shotgun_id_dict)
                    #             if genome_results:
                    #                 addGenome(genome_results)
                    #
                    #         else:
                    #             print("Couldn't find a shotgun sequence for %s \n" % (name))
                    #
                    #
                    #     elif genome_results != 'in_database' and not search_shotgun:
                    #         print(
                    #             "\nWe didn't identify any genome records for %s. And we are not attempting to search
                    # for shotgun sequenced genomes \n" % (
                    #                 name))

                elif seq_type == "profile":
                    print("it is a profile")
                    print("it is a profile")

                    # hmm_path = "static/uploads/" + filename

                    file = utilities.open_file(filename)
                    utilities.save_profile(file)

                    # while not os.path.exists(hmm_path):
                    #     time.sleep(1)
                    # if os.path.isfile(hmm_path):
                    #     print('path for hmm is ' + hmm_path)
                    #     file = open(hmm_path, 'rb')
                    #
                    #     utilities.save_profile(file)

                elif seq_type == "alignment":

                    with open("static/uploads/" + filename, "rb") as aln_file:
                        aln = models.AlignmentRecords(
                            name=filename.split(".")[0],
                            alignment=aln_file.read(),
                            tool="UploadedAlignment",
                        )
                        aln.save()

                elif seq_type == "tree":
                    print("Uploading a tree")

                    with open("static/uploads/" + filename, "rb") as tree_file:

                        tree = models.TreeRecords(
                            name=filename.split(".")[0],
                            alignment="UploadedTree",
                            tree=tree_file.read(),
                            tool="UploadedTreeToolNotSpecified",
                        )
                        tree.save()

            except HTTPError as error:
                flash("There was a HTTP error. Please try again", category="error")

            seconds = timeit.default_timer() - start_time
            minutes, seconds = divmod(seconds, 60)
            hours, minutes = divmod(minutes, 60)

            periods = [("hours", hours), ("minutes", minutes), ("seconds", seconds)]
            time_string = ", ".join(
                "{} {}".format(value, name) for name, value in periods if value
            )

            print(
                "\nFINISHED GETTING RECORDS: Time taken was {} \n".format(time_string)
            )
            if failed_genomes:
                flash(
                    "The following genomes failed to map - "
                    + " ".join(x for x in failed_genomes),
                    category="error",
                )
                print(
                    "The following genomes failed to map - \n"
                    + " \n".join(x for x in failed_genomes)
                    + "\n"
                )

        elif request.method == "POST":
            print("not validated")
            print(form.genome_type)

        return self.render("upload_admin.html", form=form)

        # def is_accessible(self):
        #     if (not current_user.is_active or not
        #     current_user.is_authenticated):
        #         return False
        #     return True


class SetupView(BaseView):
    @login_required
    @expose("/", methods=("GET", "POST"))
    def setup(self):
        form = forms.SetupForm()
        # del_form = forms.DeleteForm()
        current = models.User.objects().get(username=str(current_user.username))
        prefs = {
            "page_size": current.page_size,
            "record_size": current.record_size,
            "references": current_user.references,
        }

        if request.method == "POST":
            if form.submit.data:
                try:
                    models.User.objects().get(
                        username=str(current_user.username)
                    ).update(page_size=form.page_size.data)
                    print("data")
                    models.User.objects().get(
                        username=str(current_user.username)
                    ).update(record_size=form.record_size.data)

                    session["page_size"] = (
                        form.page_size.data
                        if form.page_size.data != None
                        else session["page_size"]
                    )
                    session["record_size"] = (
                        form.record_size.data
                        if form.record_size.data != None
                        else session["page_size"]
                    )

                    for reference in request.form.getlist("references"):
                        models.User.objects().update(add_to_set__references=reference)

                    flash("User preferences updated", category="success")
                    return self.render("setup.html", form=form, prefs=prefs)

                except Exception as e:
                    print(e)
                    flash("Something went wrong", category="error")
                    return self.render("setup.html", form=form, prefs=prefs)

            else:
                return self.render("setup.html", form=form, prefs=prefs)

        elif request.method == "GET":
            return self.render("setup.html", form=form, prefs=prefs)


class RegionView(BaseView):
    @login_required
    @expose("/", methods=("GET", "POST"))
    def setup(self, outgroup_choices=None):

        if session.get("outgroup_choices") is None:
            outgroup_choices = ["Select tree first"]
            tree_to_reroot = None
        else:
            tree = models.TreeRecords.objects().get(name=session["outgroup_choices"])

            alignment = (
                models.AlignmentRecords.objects()
                .get(name=tree.alignment)
                .alignment.read()
                .decode()
            )

            temp_path = "./tmp/tmp_alignment.aln"

            with open(temp_path, "w+") as temp_align:
                temp_align.write(alignment)

            while not os.path.exists(temp_path):
                time.sleep(1)

            if os.path.isfile(temp_path):
                seqs = utilities.read_fasta(temp_path)
                seq_names = [x for x in seqs.keys()]

            outgroup_choices = seq_names
            tree_to_reroot = tree.name

        print("ogc")
        print(outgroup_choices)

        upload_form = forms.UploadRegion()
        region_form = forms.RegionForm()
        regions_download_form = forms.RegionsDownloadForm()
        alignment_form = forms.AlignmentForm()
        tree_form = forms.MakeTreeForm()
        reroot_tree_form = forms.RerootTreeForm()

        if upload_form.upload_submit.data:
            name = upload_form.name.data

            region = models.RegionRecords(
                name=name, regions=request.files["file"].read()
            )

            region.save()

            flash("Saved " + name + " to Region Records", category="success")

        if region_form.search_regions.data:
            print("Searching for regions with profiles")

            region_to_search = region_form.region.data
            profiles = region_form.profiles.data

            region_dict, domain_dict = utilities.search_regions_with_profiles(
                region_to_search, profiles
            )

            rtp_id = utilities.randstring(5)

            region_to_profile = models.RegionToProfileRecords(
                rtp_id=rtp_id,
                region=region_to_search,
                profiles=profiles,
                region_dict=region_dict,
                domain_dict=domain_dict,
            )

            region_to_profile.save()

            flash(
                "Saved results of searching " + region_to_search + " as " + rtp_id,
                category="success",
            )

        if alignment_form.align.data:
            # We want to make an alignment

            print("Making an alignment \n")

            aln_name = alignment_form.name.data
            region_name = alignment_form.region.data
            tool = alignment_form.tool.data
            regions = models.RegionRecords.objects.get(
                name=region_name
            ).regions.decode()
            aln_path = utilities.make_alignment_from_regions(aln_name, regions, tool)

            with open(aln_path, "rb") as aln_file:
                aln = models.AlignmentRecords(
                    name=aln_name, alignment=aln_file.read(), tool=tool
                )
                aln.save()

            flash(
                "Made alignment "
                + aln_name
                + " based on "
                + region_name
                + " with "
                + tool,
                category="success",
            )

        if tree_form.make_tree.data:
            print("Making a tree \n")

            tree_name = tree_form.name.data
            alignment_name = tree_form.alignment.data
            tool = tree_form.tool.data
            alignment = (
                models.AlignmentRecords.objects.get(name=alignment_name)
                .alignment.read()
                .decode()
            )

            tree_path = utilities.make_tree(tree_name, alignment, tool)

            with open(tree_path, "rb") as tree_file:
                print("here")
                print(tree_name)
                print(alignment_name)
                print(tool)
                tree = models.TreeRecords(
                    name=tree_name,
                    alignment=alignment_name,
                    tree=tree_file.read(),
                    tool=tool,
                )
                tree.save()

            flash(
                "Made tree "
                + tree_name
                + " based on "
                + alignment_name
                + " with "
                + tool,
                category="success",
            )

        if request.method == "POST" and reroot_tree_form.reroot_tree.data:
            print("rerooting tree")
            tree_name = reroot_tree_form.tree.data
            rerooted_tree_name = reroot_tree_form.rerooted_tree_name.data
            reroot_on_seq = reroot_tree_form.seq.data

            print(tree_name)

            tree = models.TreeRecords.objects.get(name=tree_name)

            phylo_tree = PhyloTree(tree.tree.decode())

            phylo_tree.set_outgroup(reroot_on_seq)

            print(phylo_tree)

            print(os.getcwd())

            tree_path = "./tmp/tmptree.nwk"

            phylo_tree.write(outfile=tree_path)

            while not os.path.exists(tree_path):
                time.sleep(1)

            with open(tree_path, "rb") as tree_file:

                rerooted_tree = models.TreeRecords(
                    name=rerooted_tree_name,
                    alignment=tree.alignment,
                    tree=tree_file.read(),
                    tool=tree.tool,
                )
                rerooted_tree.save()

        elif request.method == "POST" and regions_download_form.download_regions.data:

            download_regions = models.RegionRecords.objects().get(
                name=regions_download_form.regions_to_download.data
            )

            regions_path = f"./fasta_folder/{download_regions.name}.fasta"

            with open(regions_path, "w+") as regions_file:
                regions_file.write(download_regions.regions.decode())
            flash(f"Wrote {download_regions.name} to {regions_path} ")

        region_names = [region.name for region in models.RegionRecords.objects()]
        region_to_profile_names = [
            region_to_profile.rtp_id
            + "_"
            + region_to_profile.region
            + " ( "
            + str(len(region_to_profile.profiles))
            + " "
            "profiles "
            ")"
            for region_to_profile in models.RegionToProfileRecords.objects()
        ]

        align_names = [align.name for align in models.AlignmentRecords.objects()]
        tree_names = [tree.name for tree in models.TreeRecords.objects()]

        region_choices = [
            (region.name, region.name) for region in models.RegionRecords.objects()
        ]
        profile_choices = [
            (profile.id, profile.name) for profile in models.Profile.objects()
        ]
        align_choices = [
            (align.name, align.name) for align in models.AlignmentRecords.objects()
        ]
        tree_choices = [(tree.name, tree.name) for tree in models.TreeRecords.objects()]
        outgroup_choices = list(zip(outgroup_choices, outgroup_choices))

        tree_choices.insert(0, (None, "Click to select tree"))

        region_form.region.choices = region_choices
        regions_download_form.regions_to_download.choices = region_choices

        region_form.profiles.choices = profile_choices

        alignment_form.region.choices = region_choices
        tree_form.alignment.choices = align_choices
        reroot_tree_form.tree.choices = tree_choices
        reroot_tree_form.seq.choices = outgroup_choices

        return self.render(
            "regions.html",
            upload_form=upload_form,
            region_form=region_form,
            regions_download_form=regions_download_form,
            alignment_form=alignment_form,
            tree_form=tree_form,
            reroot_tree_form=reroot_tree_form,
            region_names=region_names,
            region_to_profile_names=region_to_profile_names,
            align_names=align_names,
            tree_names=tree_names,
            tree_to_reroot=tree_to_reroot,
        )


class ProfilesView(BaseView):
    @login_required
    @expose("/", methods=("GET", "POST"))
    def region_to_profiles(self):
        view_profiles_on_alignment_form = forms.ViewProfilesOnAlignmentForm()
        form = forms.SelectRegionToProfilesForm()

        align_choices = [
            (align.name, align.name) for align in models.AlignmentRecords.objects()
        ]
        profile_choices = [
            (rtp.id, rtp.rtp_id) for rtp in models.RegionToProfileRecords.objects()
        ]

        view_profiles_on_alignment_form.alignment_name.choices = align_choices
        view_profiles_on_alignment_form.profiles.choices = profile_choices

        region_to_profiles = models.RegionToProfileRecords.objects()[0]

        alignment = None

        if (
            request.method == "POST"
            and view_profiles_on_alignment_form.view_profiles.data
        ):
            print("view on alignment")

            alignment = (
                models.AlignmentRecords.objects()
                .get(name=view_profiles_on_alignment_form.alignment_name.data)
                .alignment
            )

            region_to_profiles = models.RegionToProfileRecords.objects().get(
                id=view_profiles_on_alignment_form.profiles.data
            )

            alignment = utilities.colour_alignment_by_profiles(
                alignment.read().decode(), region_to_profiles
            )

        if request.method == "POST" and form.submit.data:
            region_to_profiles_name = form.name.data
            print(region_to_profiles_name)
            region_to_profiles = models.RegionToProfileRecords.objects().get(
                id=region_to_profiles_name
            )

        form.name.choices = [
            (region_to_profiles.id, region_to_profiles.rtp_id)
            for region_to_profiles in models.RegionToProfileRecords.objects()
        ]

        return self.render(
            "region_to_profiles.html",
            view_profiles_on_alignment_form=view_profiles_on_alignment_form,
            form=form,
            rtp=region_to_profiles,
            alignment=alignment,
        )


class AlignmentsView(BaseView):
    @login_required
    @expose("/", methods=("GET", "POST"))
    def alignments(self):
        form = forms.SelectAlignmentForm()
        alignment_download_form = forms.AlignmentDownloadForm()
        alignment_choices = [
            (aln.name, aln.name) for aln in models.AlignmentRecords.objects()
        ]
        alignment_download_form.alignment.choices = alignment_choices

        if request.method == "POST" and form.submit.data:
            alignment = models.AlignmentRecords.objects().get(id=form.name.data)

        elif (
            request.method == "POST" and alignment_download_form.download_alignment.data
        ):
            alignment = None
            download_alignment = models.AlignmentRecords.objects().get(
                name=alignment_download_form.alignment.data
            )

            aln_path = f"./fasta_folder/{download_alignment.name}.aln"

            with open(aln_path, "w+") as aln_file:
                aln_file.write(download_alignment.alignment.read().decode())
            flash(f"Wrote {download_alignment.name} to {aln_path} ")

        elif request.method == "GET":
            alignment = None

        form.name.choices = [
            (align.id, align.name) for align in models.AlignmentRecords.objects()
        ]

        return self.render(
            "alignments.html",
            form=form,
            alignment_download_form=alignment_download_form,
            align_data=alignment,
        )


class TreeView(BaseView):
    @login_required
    @expose("/", methods=("GET", "POST"))
    def tree(self):
        tree = None
        species_tree = None
        tree_select_form = forms.TreeSelectForm()
        tree_download_form = forms.TreeDownloadForm()
        species_tree_form = forms.SpeciesTreeForm()

        if request.method == "POST" and tree_select_form.submit.data:
            tree = models.TreeRecords.objects().get(
                id=tree_select_form.tree_select_name.data
            )
            tag_dict = getGenomes.get_tags()

            if tree_select_form.profiles.data == "None":
                region_dict = {}
            else:
                region_dict = (
                    models.RegionToProfileRecords.objects()
                    .get(rtp_id=tree_select_form.profiles.data)
                    .region_dict
                )

            if tree_select_form.region_order.data == "None":
                region_order_dict = {}
            else:
                region_order_dict = (
                    models.RegionOrderRecords.objects()
                    .get(name=tree_select_form.region_order.data)
                    .region_order_dict
                )

            if tree_select_form.sequence_content.data == "None":
                sequence_content_dict = {}
            else:
                sequence_content_dict = utilities.get_sequence_content_dict(
                    tree_select_form.sequence_content.data
                )

            full_names = tree_select_form.full_names.data

            collapse_on_genome_tags = tree_select_form.collapse_on_genome_tags.data

            display_circular = tree_select_form.display_circular.data
            display_circular_180 = tree_select_form.display_circular_180.data

            colour_dict = {
                "Type1": "dodgerblue",
                "type1": "dodgerblue",
                "Type1_A2": "white",
                "Type2b": "gold",
                "Type2a": "green",
                "Type2a_force": "green",
                "Type3": "purple",
                "Type3_A2": "purple",
                "Multiple": "red",
                "unknown": "black",
                "Single": "brown",
                "Type?": "pink",
                "Type5": "pink",
            }

        elif request.method == "POST" and tree_download_form.download_tree.data:
            download_tree = models.TreeRecords.objects().get(
                name=tree_download_form.tree.data
            )

            tree_path = f"./fasta_folder/{download_tree.name}.nwk"

            with open(tree_path, "w+") as tree_file:
                tree_file.write(download_tree.tree.decode())
            flash(f"Wrote {download_tree.name} to {tree_path} ")

        # We want to make a species tree and annotate our tree with species info
        elif request.method == "POST" and species_tree_form.make_tree.data:
            species_tree = models.TreeRecords.objects().get(
                name=tree_download_form.tree.data
            )

            # Check if the tree we want to annotate actually has taxids on it
            if "_taxid_" not in species_tree.tree.decode():
                flash(
                    f"{species_tree.name} does not have taxonomic IDs annotated on the newick file. Recreate the "
                    f"Regions file you used to make this tree.",
                    category="error",
                )
                tree_img = ""
                ncbi_tree_img = ""

            else:
                tree_img = utilities.get_species_tree_image(
                    species_tree.tree.decode(), species_tree.name
                )
                tree_img = "/" + tree_img + "#" + utilities.randstring(5)
                ncbi_tree_img = (
                    tree_img.split(".png")[0]
                    + "_ncbi.png"
                    + "#"
                    + utilities.randstring(5)
                )

                print(tree_img)
                print(ncbi_tree_img)

        elif request.method == "GET":
            tree = None
            species_tree = None
            # if models.TreeRecords.objects().count() == 0:
            #     tree = None
            # else:
            #     tree = models.TreeRecords.objects()[0]
            #     tag_dict = {}
            #     region_dict = {}
            #     colour_dict = {}

        if tree:
            tree_img = utilities.get_tree_image(
                tree.tree.decode(),
                tree.name,
                tag_dict,
                region_dict,
                region_order_dict,
                sequence_content_dict,
                colour_dict,
                full_names,
                collapse_on_genome_tags,
                display_circular,
                display_circular_180,
            )

            print(tree_img)
            tree_img = "/" + tree_img + "#" + utilities.randstring(5)

            print(tree_img)
            ncbi_tree_img = ""
            # if tree_img:
            #     tree_img = tree_img.split("static/")[1]
        elif not species_tree:

            tree_img = ""
            ncbi_tree_img = ""

        profile_choices = [
            (rtp.rtp_id, rtp.rtp_id) for rtp in models.RegionToProfileRecords.objects()
        ]

        region_order_choices = [
            (region_order.name, region_order.name)
            for region_order in models.RegionOrderRecords.objects()
        ]
        alignment_choices = [
            (aln.name, aln.name) for aln in models.AlignmentRecords.objects()
        ]

        # Insert a None option in case we don't want to add certain information
        profile_choices.insert(0, (None, None))
        region_order_choices.insert(0, (None, None))
        alignment_choices.insert(0, (None, None))

        tree_choices = [(tree.name, tree.name) for tree in models.TreeRecords.objects()]

        tree_select_form.tree_select_name.choices = [
            (tree.id, tree.name) for tree in models.TreeRecords.objects()
        ]
        tree_select_form.profiles.choices = profile_choices
        tree_select_form.region_order.choices = region_order_choices
        tree_select_form.sequence_content.choices = alignment_choices

        tree_download_form.tree.choices = tree_choices

        species_tree_form.tree.choices = tree_choices

        return self.render(
            "trees.html",
            tree_select_form=tree_select_form,
            tree_download_form=tree_download_form,
            species_tree_form=species_tree_form,
            tree_img=tree_img,
            ncbi_tree_img=ncbi_tree_img,
        )


class DownloadFastaView(BaseView):
    @login_required
    @expose("/", methods=("GET", "POST"))
    def setup(self):
        form = forms.DownloadFastaForm()
        associated_regions = forms.DownloadAssociatedRegions()
        tags = forms.DownloadTags()

        if request.method == "POST":
            if form.submit.data:
                try:
                    print(form.region.data)

                    region = form.region.data
                    translate = form.translate.data
                    filename = form.filename.data
                    include_genome = [
                        x.strip() for x in form.include_genome.data.split(",")
                    ]
                    exclude_genome = [
                        x.strip() for x in form.exclude_genome.data.split(",")
                    ]
                    include_hits = [
                        x.strip() for x in form.include_hits.data.split(",")
                    ]
                    exclude_hits = [
                        x.strip() for x in form.exclude_hits.data.split(",")
                    ]
                    align = form.align.data
                    split_strands = False

                    print("here they come")
                    print("region")
                    print(region)
                    print("filename")
                    print(filename)
                    print("Include genome ")
                    print(include_genome)
                    print("Exclude genome ")
                    print(exclude_genome)
                    print("Include hits ")
                    print(include_hits)
                    print("Exclude hits ")
                    print(exclude_hits)
                    print("translate")
                    print(translate)
                    print("align")
                    print(align)
                    print("split strands")
                    print(split_strands)
                    # if include_genome == [""]:
                    #     include_genome = None
                    #
                    # if exclude_genome == [""]:
                    #     exclude_genome = None
                    #
                    # if include_hits == [""]:
                    #     include_hits = None
                    #
                    # if exclude_hits == [""]:
                    #     exclude_hits = None

                    outpath = getGenomes.download_fasta_regions(
                        region,
                        filename,
                        include_genome,
                        exclude_genome,
                        include_hits,
                        exclude_hits,
                        translate,
                        align,
                        split_strands,
                    )

                    flash(
                        "Downloaded " + region + " file to " + outpath,
                        category="success",
                    )
                    return self.render(
                        "download_fasta.html",
                        form=form,
                        associated_regions=associated_regions,
                        tags=tags,
                    )

                except Exception as e:
                    print(e)
                    flash(e, category="error")
                    return self.render(
                        "download_fasta.html",
                        form=form,
                        associated_regions=associated_regions,
                        tags=tags,
                    )

            elif associated_regions.associated_regions.data:

                outpath = getGenomes.download_associated_regions()

                flash(
                    "Downloaded associated regions file to " + outpath,
                    category="success",
                )

                return self.render(
                    "download_fasta.html",
                    form=form,
                    associated_regions=associated_regions,
                    tags=tags,
                )

            elif tags.tags.data:

                outpath = getGenomes.download_tags()

                flash("Downloaded tags file to " + outpath, category="success")

                return self.render(
                    "download_fasta.html",
                    form=form,
                    associated_regions=associated_regions,
                    tags=tags,
                )

            else:
                return self.render(
                    "download_fasta.html",
                    form=form,
                    associated_regions=associated_regions,
                    tags=tags,
                )

        elif request.method == "GET":
            return self.render(
                "download_fasta.html",
                form=form,
                associated_regions=associated_regions,
                tags=tags,
            )


class DownloadRegionOrderView(BaseView):
    @login_required
    @expose("/", methods=("GET", "POST"))
    def setup(self):
        form = forms.DownloadRegionOrder()

        region_order_names = [
            region_order.name for region_order in models.RegionOrderRecords.objects()
        ]

        if request.method == "POST":
            if form.submit.data:
                try:
                    fasta_dict = {}
                    include_genome = form.include_genome.data.split(",")
                    exclude_genome = form.exclude_genome.data.split(",")
                    include_hits = form.include_hits.data.split(",")
                    exclude_hits = form.exclude_hits.data.split(",")

                    if include_genome == [""]:
                        genomes = models.GenomeRecords.objects()
                    else:
                        genomes = models.GenomeRecords.objects(tags__in=include_genome)

                    getGenomes.write_region_order(
                        genomes,
                        exclude_hits=exclude_hits,
                        save_to_db=form.save_to_db.data,
                    )

                    flash(
                        "Downloaded file to fasta_folder/region_order.txt",
                        category="success",
                    )

                except Exception as e:
                    print(e)
                    flash(e, category="error")

        return self.render(
            "download_region_order.html",
            form=form,
            region_order_names=region_order_names,
        )


class DownloadMLGOView(BaseView):
    @login_required
    @expose("/", methods=("GET", "POST"))
    def setup(self):
        form = forms.DownloadMLGO()
        tree_form = forms.DownloadMLGOTree()

        tree_choices = [(tree.name, tree.name) for tree in models.TreeRecords.objects()]

        tree_form.tree_select_name.choices = tree_choices

        if request.method == "POST" and form.submit.data:
            if form.submit.data:
                try:
                    fasta_dict = {}
                    include_genome = form.include_genome.data.split(",")
                    exclude_genome = form.exclude_genome.data.split(",")
                    include_hits = form.include_hits.data.split(",")
                    exclude_hits = form.exclude_hits.data.split(",")

                    if include_genome == [""]:
                        genomes = models.GenomeRecords.objects(tags__nin=exclude_genome)
                    else:
                        genomes = models.GenomeRecords.objects(
                            tags__in=include_genome, tags__nin=exclude_genome
                        )

                    getGenomes.write_mlgo_order(genomes, ref_mlgo_dict=ref_mlgo_dict)

                    flash(
                        "Downloaded file to fasta_folder/mlgo.txt", category="success"
                    )

                except Exception as e:
                    print(e)
                    flash(e, category="error")
        if request.method == "POST" and tree_form.download_mlgo_tree.data:
            tree = models.TreeRecords.objects.get(name=tree_form.tree_select_name.data)

            tree_path = "fasta_folder/" + tree.name + "_mlgo.nwk"

            print(tree)

            getGenomes.write_mlgo_tree(tree.tree, tree_path)

            flash("Downloaded file to " + tree_path, category="success")

        return self.render("download_MLGO.html", form=form, tree_form=tree_form)


class VisualiseMLGOView(BaseView):
    @login_required
    @expose("/", methods=("GET", "POST"))
    def setup(self):
        upload_form = forms.UploadMLGOTree()
        select_form = forms.SelectMLGOTree()

        tree_img = ""

        tree_choices = [
            (tree.name, tree.name) for tree in models.MLGOTreeRecords.objects()
        ]

        print(tree_choices)

        select_form.select_name.choices = tree_choices

        if request.method == "POST" and upload_form.upload.data:
            # annotated_tree = open(upload_form.annotated_tree.data, "r+")
            # gene_order = open(upload_form.gene_order.data, "r+")
            # mlgo_tree = models.MLGOTreeRecords(upload_form.upload_name.data, annotated_tree,
            #                                    gene_order)

            mlgo_dict = utilities.get_mlgo_dict(
                upload_form.gene_order.data.read().decode()
            )

            mlgo_tree = models.MLGOTreeRecords(
                upload_form.upload_name.data, upload_form.annotated_tree.data, mlgo_dict
            )

            mlgo_tree.save()

            flash(
                "Uploaded ML Gene Order tree "
                + upload_form.upload_name.data
                + " to database"
            )

        if request.method == "POST" and select_form.select.data:
            tree = models.MLGOTreeRecords.objects().get(
                name=select_form.select_name.data
            )

            # print(tree)
            # print (tree.tree.read().decode())

            tree_img = utilities.get_ml_go_tree_image(
                tree.tree.read().decode(),
                tree.name,
                tree.ancestral_orders,
                mlgo_ref_dict,
            )

            tree_img = "/" + tree_img + "#" + utilities.randstring(5)

        return self.render(
            "visualise_MLGO.html",
            upload_form=upload_form,
            select_form=select_form,
            tree_img=tree_img,
        )


class SequenceRecordsView(ModelView):
    edit_modal = True
    can_create = False
    can_view_details = True
    form_edit_rules = ("name", "species", "description", "sequence")

    # list_template = 'list.html'
    # create_template = 'create.html'
    # edit_template = 'edit.html'

    column_formatters = {"sequence": custom_filters.seqdescription_formatter}

    @action("generate_profile", "Generate a profile from these sequences")
    def action_generate_profile(self, ids):
        try:
            query = models.SequenceRecords.objects(id in ids)
            print(query)
            align_list = []
            for record in query.all():
                align_record = SeqRecord(
                    Seq(record.sequence, generic_protein),
                    id=str(record.name) + "_" + "seq",
                )
                align_list.append(align_record)

            utilities.createProfile(align_list)

        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise
            flash(gettext(ex))


class GenomeRecordsView(ModelView):

    form_edit_rules = (
        "name",
        "species",
        "organism",
        "taxid",
        "assembly_type",
        "sequence",
        "hits",
        "release_type",
        "assembly_level",
        "genome_representation",
        "assembly_name",
        "biosample",
        "bioproject",
        "date",
        "wgs_project",
        "genome_coverage",
        "expected_final_version",
        "excluded",
        "genbank_accession_id",
        "refseq_accession_id",
        "r_g_identical",
        "plasmid",
        "description",
    )

    # column_list = (
    #     'name', 'species', 'strain', 'description', 'sequence', 'hits', 'Genome Overview', 'Expanded Genome Overview')
    #
    column_list = (
        "name",
        "species",
        "organism",
        "taxid",
        "assembly_type",
        "sequence",
        "hits",
        "release_type",
        "assembly_level",
        "genome_representation",
        "assembly_name",
        "biosample",
        "bioproject",
        "date",
        "wgs_project",
        "genome_coverage",
        "expected_final_version",
        "excluded",
        "genbank_accession_id",
        "refseq_accession_id",
        "r_g_identical",
        "plasmid",
        "description",
    )

    column_sortable_list = (
        "name",
        "species",
        "organism",
        "taxid",
        "assembly_type",
        "sequence",
        "hits",
        "release_type",
        "assembly_level",
        "genome_representation",
        "assembly_name",
        "biosample",
        "bioproject",
        "date",
        "wgs_project",
        "genome_coverage",
        "expected_final_version",
        "excluded",
        "genbank_accession_id",
        "refseq_accession_id",
        "r_g_identical",
        "plasmid",
        "description",
    )

    # column_searchable_list = ('name')

    @login_required
    @expose("/", methods=("GET", "POST"))
    def index_view(self):
        print(session.get("page_size"))
        # current = models.User.objects().get(username=str(current_user.username))
        # print(current.page_size)
        self.edit_modal = True
        self.can_create = False
        self.can_view_details = True
        self.can_set_page_size = True
        self.search_placeholder()

        self.column_searchable_list = ["species"]

        # self.page_size = current.page_size

        self.page_size = session["page_size"]

        self.form_edit_rules = ("name", "species", "strain", "plasmid", "description")

        self.column_list = (
            "name",
            "species",
            "strain",
            "plasmid",
            "future",
            "description",
            "sequence",
            "hits",
            "Genome Overview",
        )

        return super(ModelView, self).index_view()

    # def is_accessible(self):
    #     if (not current_user.is_active or not
    #     current_user.is_authenticated):
    #         return False
    #     return True

    # def _download_formatter(self, context, model, name):
    #     return Markup(
    #         "<a href='{url}' target='_blank'>View Genome</a>".format(
    #             url=self.get_url('download_genome_overview', id=model.name)))
    #
    # def _expanded_download_formatter(self, context, model, name):
    #     return Markup(
    #         "<a href='{url}' target='_blank'>View Genome</a>".format(
    #             url=self.get_url('download_genome_expanded_overview', id=model.name)))

    column_formatters = {
        "sequence": custom_filters.seqdescription_formatter,
        "hits": custom_filters.hitdescription_formatter,
        # 'Genome Overview': _download_formatter,
        # 'Expanded Genome Overview': _expanded_download_formatter,
    }

    for name in ref_names:
        exec(
            '@action("'
            + name
            + '_check_with_profile", "Check for '
            + name
            + ' region with a profile")\ndef '
            + name
            + '_check_with_profile(self, profile_ids):\n utilities.check_with_profile(profile_ids, "'
            + name
            + '")'
        )


#
# class ProfileView(ModelView):
#     """
#     View of the Profile database for storing HMM Profiles """
#     column_list = (
#         'name', 'a1_profile_ref', 'a2_profile_ref', 'pore_profile_ref', 'chitinase_profile_ref',
# 'region1_profile_ref',
#         'region2_profile_ref',
#         'region3_profile_ref', 'region4_profile_ref', 'download')
#     form_columns = ('name', 'profile')
#
#     form_extra_fields = {'profile': models.BlobUploadField(
#         label='File',
#         allowed_extensions=['pdf', 'doc', 'docx', 'xls', 'xlsx', 'png', 'jpg', 'jpeg', 'gif', 'hmm', 'fa', 'fasta'],
#         size_field='size',
#         filename_field='filename',
#         mimetype_field='mimetype'
#     )}
#
#     def _download_formatter(self, context, model, name):
#         return Markup(
#             "<a href='{url}' target='_blank'>Download</a>".format(url=self.get_url('download_blob', id=model.uid)))
#
#     column_formatters = {
#         'download': _download_formatter,
#     }


class GenomeOverviewView(BaseView):
    @login_required
    @expose("/", methods=("GET", "POST"))
    def genomeoverview(self):
        form = forms.GenomeOverviewSelectForm()
        # form.genome.choices = [(genome.name, genome.species) for genome in models.GenomeRecords.query.all()]

        queries = models.GenomeRecords.objects()
        items = []

        for query in queries:
            query_details = {}
            query_details["src"] = (
                query.name + "_" + query.species.replace(" ", "_") + ".png"
            )
            query_details["srct"] = (
                query.name + "_" + query.species.replace(" ", "_") + ".png"
            )
            query_details["title"] = str(query.name) + "\n" + query.species
            query_details["tags"] = " ".join(tag for tag in query.tags)

            items.append(query_details)

        print(items)

        return self.render("genomeoverview.html", form=form, items=items)


class TempFixView(BaseView):
    @login_required
    @expose("/", methods=("GET", "POST"))
    def temp_fix(self):
        form = forms.TempFixForm()
        return self.render("temp_fix.html", form=form)


class AutomaticTaggingView(BaseView):
    @login_required
    @expose("/", methods=("GET", "POST"))
    def automatic_tagging(self):

        tag_simple_form = forms.TagSimpleForm()
        search_for_promoters_form = forms.SearchForPromoters()
        update_tags_form = forms.UpdateTagsForm()
        auto_hide_form = forms.AutoHideRegionsForm()
        auto_classify_form = forms.AutoClassifyForm()
        auto_classify_test_form = forms.AutoClassifyTestForm()

        unique_tags = models.GenomeRecords.objects().distinct(field="tags")

        regions_list = [
            "TcdA1_expanded",
            "A1_expanded",
            "A2_expanded",
            "TcB_expanded",
            "TcC_expanded",
            "Chitinase_expanded",
        ]

        # The list of tags to choose from
        update_tags_form.old_tag.choices = list(zip(unique_tags, unique_tags))
        auto_hide_form.hide_include_genome.choices = list(zip(unique_tags, unique_tags))
        auto_hide_form.hide_exclude_genome.choices = list(zip(unique_tags, unique_tags))

        auto_hide_form.auto_hide_region.choices = list(zip(regions_list, regions_list))

        unique_tags.insert(0, "Test all")
        auto_classify_test_form.limit_classify_test_tagged.choices = list(
            zip(unique_tags, unique_tags)
        )
        auto_classify_test_form.skip_tags.choices = list(zip(unique_tags, unique_tags))

        test_results = ""

        if request.method == "POST" and tag_simple_form.tag_simple.data:
            include_genome = [
                x.strip() for x in tag_simple_form.include_genome.data.split(",")
            ]
            exclude_hits = [
                x.strip() for x in tag_simple_form.exclude_hits.data.split(",")
            ]

            if include_genome == [""]:
                # specific_queries = ['JPPA01', 'CABPSQ01', 'FNJL01', 'NITZ01', 'RCWL01', 'NC_013216.1', 'VUOC01',
                #                     'FQYP01', 'VSRQ01',
                #                     'NZ_CP010897.2', 'FNKR01', 'CABPSP01', 'BCQP01', 'QGTQ01', 'RCZD01', 'MVGR01',
                #                     'NZ_CP020080.1', 'QHJL01', 'NZ_CP013047.2', 'FMVH01', 'BILZ01', 'MQZX01',
                #                     'CAAJVF01', 'PIQI01', 'NC_009253.1', 'NZ_CP013461.1', 'NZ_CM001559.1', 'SMDG01',
                #                     'JAAQYP01', 'LVTS01', 'SLVP01', 'SSNZ01', 'LYRP01', 'NZ_CP054613.1',
                #                     'NZ_CP038255.1', 'LMXH01', 'SAXA01', 'BJMH01', 'NZ_CP024085.1', 'PQMB01', 'SMSL01',
                #                     'AEDD01', 'RKHU01']
                # genomes = models.GenomeRecords.objects(name__in=specific_queries)

                genomes = models.GenomeRecords.objects()

            else:
                genomes = models.GenomeRecords.objects(tags__in=include_genome)

            for g in genomes:
                print(g.name)

            getGenomes.tag_as_simple(genomes, exclude_hits)

        if (
            request.method == "POST"
            and search_for_promoters_form.search_for_promoters.data
        ):
            mismatch = search_for_promoters_form.mismatch.data
            if not mismatch:
                mismatch = 0

            utilities.search_for_promoters(mismatch=mismatch)

            flash("Finished searching for promoters", category="success")

        if request.method == "POST" and update_tags_form.update_tags.data:
            print("update tags")
            old_tag = update_tags_form.old_tag.data
            new_tag = update_tags_form.new_tag.data
            # db.getCollection("genome_records").updateMany({"tags": "monkey"}, {"$set": {"tags.$": "possum_face"}})
            # BlogPost.objects(tags="mongoEngien").update(set__tags__S="MongoEngine")

            models.GenomeRecords.objects(tags=old_tag).update(set__tags__S=new_tag)

            queries = models.GenomeRecords.objects()

            for query in queries:
                for hit in query.hits:
                    new_tags = [new_tag if x == old_tag else x for x in hit.tags]
                    print("new tags was ")

                    tags_unique = list(set(new_tags))

                    print(new_tags)
                    print(tags_unique)
                    hit.tags = tags_unique

                query.save()

            # models.Hits.objects(tags=old_tag).update(set__tags__S=new_tag)

            # models.GenomeRecords.objects(hits__tags=old_tag).update(set__tags__=new_tag)
            #
            # genomes = models.GenomeRecords.objects().all()
            #
            # mg = models.GenomeRecords.objects(hits__tags=old_tag)
            #
            #
            #
            # print (mg)
            #
            # for g in mg:
            #     g.hits.tags = ['pop']
            #

            # models.Hits.objects(tags=old_tag).update(set__tags__S=new_tag)

            # hits = models.GenomeRecords.objects(field='hits')

            # print (hits)

            # for genome in genomes:
            #     hits = genome.hits
            #
            # #     hits.(tags=old_tag).update(set__tags__S=new_tag)
            # #
            #     for hit in hits:
            #
            #         for idx, tag in enumerate(hit.tags):
            #             if tag == old_tag:
            #                 hit.tags[idx] = new_tag
            #
            # genomes.save()

            #
            # hits = models.GenomeRecords.hits(tags=old_tag).update(set__tags__S=new_tag)

            # hits(tags=old_tag).update(set__tags__S=new_tag)

            # models.GenomeRecords.objects(tags__hits=old_tag).update(set__tags__hits__S=new_tag)

            # print (new_tag)

            models.GenomeTags.objects(tags=old_tag).update(set__tags__S=new_tag)

            flash(
                "Changed all genome tags - " + old_tag + " into " + new_tag,
                category="success",
            )

        if request.method == "POST" and auto_hide_form.hide.data:

            include_genome = auto_hide_form.hide_include_genome.data
            exclude_genome = auto_hide_form.hide_exclude_genome.data
            hide_region = auto_hide_form.auto_hide_region.data

            if include_genome == [""]:
                genomes = models.GenomeRecords.objects(tags__nin=exclude_genome)
            else:
                genomes = models.GenomeRecords.objects(
                    tags__in=include_genome, tags__nin=exclude_genome
                )

            print("check results here")
            print(include_genome)
            print(exclude_genome)

            for g in genomes:
                print(g)
                for hit in g.hits:
                    print(hit)
                    print(hit.region)
                    if hit.region == hide_region:
                        print(hit.region)
                        print(hit.tags)
                        if "hidden" not in hit.tags:
                            print("adding hidden")
                            hit.tags.append("hidden")
                g.save()

        if request.method == "POST" and auto_classify_test_form.auto_classify_test.data:
            limit_classify_test_tagged = (
                auto_classify_test_form.limit_classify_test_tagged.data
            )
            skip_tags = auto_classify_test_form.skip_tags.data

            if limit_classify_test_tagged == "Test all":
                queries = models.GenomeRecords.objects.all().timeout(False)
            else:
                # specific_queries = ['JPPA01', 'CABPSQ01', 'FNJL01', 'NITZ01', 'RCWL01', 'NC_013216.1', 'VUOC01',
                #                     'FQYP01', 'VSRQ01',
                #                     'NZ_CP010897.2', 'FNKR01', 'CABPSP01', 'BCQP01', 'QGTQ01', 'RCZD01', 'MVGR01',
                #                     'NZ_CP020080.1', 'QHJL01', 'NZ_CP013047.2', 'FMVH01', 'BILZ01', 'MQZX01',
                #                     'CAAJVF01', 'PIQI01', 'NC_009253.1', 'NZ_CP013461.1', 'NZ_CM001559.1', 'SMDG01',
                #                     'JAAQYP01', 'LVTS01', 'SLVP01', 'SSNZ01', 'LYRP01', 'NZ_CP054613.1',
                #                     'NZ_CP038255.1', 'LMXH01', 'SAXA01', 'BJMH01', 'NZ_CP024085.1', 'PQMB01', 'SMSL01',
                #                     'AEDD01', 'RKHU01']
                # queries = models.GenomeRecords.objects(name__in=specific_queries, tags=limit_classify_test_tagged).timeout(False)
                queries = models.GenomeRecords.objects(
                    tags=limit_classify_test_tagged
                ).timeout(False)

            test_results = utilities.test_auto_classify(queries, skip_tags=skip_tags)

            print("**")
            print(test_results)

        return self.render(
            "automatic_tagging.html",
            tag_simple_form=tag_simple_form,
            search_for_promoters_form=search_for_promoters_form,
            update_tags_form=update_tags_form,
            auto_hide_form=auto_hide_form,
            auto_classify_form=auto_classify_form,
            auto_classify_test_form=auto_classify_test_form,
            test_results=test_results.replace("\n", "<br />\n"),
        )


@app.route("/temp_assoc_fix", methods=["GET", "POST"])
def temp_assoc_fix():
    print("Fixing the Associated Regions dict")

    queries = models.AssociatedHits.objects()

    for query in queries:
        print(query)
        print(query.region1)
        print(query.region1.split("_information")[0])

        genome_name = query.region1.split("_information")[0]
        genome = models.GenomeRecords.objects().get(name=genome_name)
        genome_id = str(genome.id)
        print("this genome id was ")
        print(genome_id)
        query.genome_id = genome_id
        query.save()

    return redirect("temp_fix")


class TrimRegionsView(BaseView):
    @login_required
    @expose("/", methods=("GET", "POST"))
    def trim_regions(self):
        trim_to_profile_form = forms.TrimToProfileForm()
        trim_around_profile_form = forms.TrimAroundProfileForm()
        region_choices = [
            (region.name, region.name) for region in models.RegionRecords.objects()
        ]
        profile_choices = [
            (profile.name, profile.name) for profile in models.Profile.objects()
        ]
        section_choices = [
            (profile.name, profile.name) for profile in models.Profile.objects()
        ]
        section_choices.insert(0, ("Content", "Content"))

        trim_to_profile_form.trim_to_region.choices = region_choices
        trim_to_profile_form.trim_to_profile.choices = profile_choices

        trim_around_profile_form.trim_around_region.choices = region_choices

        trim_around_profile_form.trim_around_profile.choices = section_choices

        if request.method == "POST" and trim_to_profile_form.trim_to.data:
            regions = (
                models.RegionRecords.objects()
                .get(name=trim_to_profile_form.trim_to_region.data)
                .regions
            )
            profile = models.Profile.objects().get(
                name=trim_to_profile_form.trim_to_profile.data
            )
            trimmed_name = trim_to_profile_form.trim_to_name.data
            failed_seqs = ["fail1", "fail2"]
            failed_seqs = utilities.trim_to_profile(
                regions, profile.profile, trimmed_name
            )

            if failed_seqs:
                print(failed_seqs)
                print(profile.name)
                flash(
                    "The following regions did not have a match for "
                    + profile.name
                    + " and so have not been added to "
                    "the new file - " + " ".join(failed_seqs),
                    category="error",
                )
            flash(
                "Regions have been trimmed to "
                + profile.name
                + " and saved as "
                + trimmed_name,
                category="success",
            )

        if (
            request.method == "POST"
            and trim_around_profile_form.trim_around_submit.data
        ):

            regions = (
                models.RegionRecords.objects()
                .get(name=trim_around_profile_form.trim_around_region.data)
                .regions
            )
            trimmed_name = trim_around_profile_form.trim_around_name.data
            section1 = trim_around_profile_form.section1.data
            section2 = trim_around_profile_form.section2.data

            sections = trim_around_profile_form.trim_around_profile.data

            print("Sections was ")
            print(sections)

            pos1 = None
            pos2 = None

            if len(sections) == 3:

                print("Profile 1 is " + sections[0])
                print("Profile 2 is " + sections[1])

                profile1 = models.Profile.objects().get(name=sections[0]).profile
                profile2 = models.Profile.objects().get(name=sections[2]).profile

                pos1 = "start" if section1 else "end"

                pos2 = "end" if section2 else "start"

            elif len(sections) == 2:
                if sections[0] == "Content":
                    profile1 = None
                    profile2 = models.Profile.objects().get(name=sections[1]).profile

                    print("Profile 1 is None")
                    print("Profile 2 is " + sections[1])
                    pos2 = "end" if section1 else "start"  # Need to check section 1

                elif sections[1] == "Content":
                    profile1 = models.Profile.objects().get(name=sections[0]).profile
                    profile2 = None

                    print("Profile 1 is " + sections[0])
                    print("Profile 2 is None")

                    pos1 = "start" if section1 else "end"  # Need to check section 1

            else:
                flash("Incorrect number of sections chosen")
                return self.render(
                    "trim_regions.html",
                    trim_to_profile_form=trim_to_profile_form,
                    trim_around_profile_form=trim_around_profile_form,
                )

            failed_seqs = utilities.trim_around_profile(
                regions, profile1, profile2, pos1, pos2, trimmed_name
            )
            if failed_seqs:
                print(failed_seqs)
                flash(
                    "The following regions did not have a match for one of the profiles and so have not been added to "
                    "the new file - " + " ".join(failed_seqs),
                    category="error",
                )
            flash(
                "Regions have been trimmed and saved as " + trimmed_name,
                category="success",
            )

        return self.render(
            "trim_regions.html",
            trim_to_profile_form=trim_to_profile_form,
            trim_around_profile_form=trim_around_profile_form,
        )


class ChartView(BaseView):
    @login_required
    @expose("/", methods=("GET", "POST"))
    def chart(self):
        chart_form = forms.ChartsForm()
        unique_tags = models.GenomeRecords.objects().distinct(field="tags")

        labels = []
        values = []

        for tag in unique_tags:
            tag_count = models.GenomeRecords.objects(tags=tag).count()
            labels.append(tag)
            values.append(tag_count)

        # The list of tags to choose from
        chart_form.select_tags.choices = list(zip(labels, labels))
        chart_form.exclude_tags.choices = list(zip(labels, labels))

        colors = [
            "#F7464A",
            "#46BFBD",
            "#FDB45C",
            "#FEDCBA",
            "#ABCDEF",
            "#DDDDDD",
            "#ABCABC",
            "#4169E1",
            "#C71585",
            "#FF4500",
            "#FEDCBA",
            "#46BFBD",
        ]

        if request.method == "POST":
            selected_vals = chart_form.data["select_tags"]
            exclude_vals = chart_form.data["exclude_tags"]

            print("Selected was ")
            print(selected_vals)
            print("Excluded was ")
            print(exclude_vals)

            labels = []
            values = []

            for tag in selected_vals:
                tag_count = models.GenomeRecords.objects(
                    __raw__={"tags": {"$in": [tag], "$nin": exclude_vals}}
                ).count()

                print(tag)
                print(tag_count)
                labels.append(tag)
                values.append(tag_count)

            print("Wow!")

        else:
            selected_vals = labels
            exclude_vals = labels

        if values:
            maxval = max(values) + 10
        else:
            maxval = 10

        return self.render(
            "charts.html",
            chart_form=chart_form,
            title="Unique tags",
            max=maxval,
            labels=labels,
            values=values,
            selected_vals=json.dumps(selected_vals),
            exclude_vals=json.dumps(exclude_vals),
        )


class BatchDeleteView(BaseView):
    @login_required
    @expose("/", methods=("GET", "POST"))
    def batch_delete(self):
        form = forms.BatchDeleteForm()
        return self.render("batch_delete.html", form=form)


class GenomeDetailView(BaseView):
    @login_required
    @expose("/", methods=("GET", "POST"))
    def genomedetail(self):

        # passed_from_page = False

        # print ("Function was passed from " + session['passed_from'])

        select_form = forms.GenomeDiagramSelectForm()
        genome_by_name_form = forms.GenomeByNameForm()
        page_form = forms.GenomeDiagramPageForm()
        region_form = forms.GenomeDiagamShowRegions()
        hit_form = forms.GenomeHitForm()

        # genome = models.GenomeRecords.objects.get(id=session['genome'])
        if session.get("hits") is None:
            session["hits"] = "expanded"

        if session.get("hidden_type") is None:
            session["hidden_type"] = True

        if session.get("show_promoters") is None:
            session["show_promoters"] = False

        if session.get("show_stop_codons") is None:
            session["show_stop_codons"] = False

        if session.get("show_existing_features") is None:
            session["show_existing_features"] = False

        if session.get("hidden_type") is None:
            session["hidden_type"] = True

        if session.get("checked_regions") is None:
            session["checked_regions"] = [
                "A1",
                "A2",
                "TcdA1",
                "TcB",
                "TcC",
                "Chitinase",
                "region1",
            ]

        if session.get("page_choice") is None:
            session["page_choice"] = 0

        if session.get("untagged") is None:
            session["untagged"] = False

        if session.get("limit_genomes") is None:
            session["limit_genomes"] = False

        if session.get("genome_tagged") is None:
            session["genome_tagged"] = [""]

        if session.get("genome") is None:
            print("IT WAS NONE")
            genome = models.GenomeRecords.objects()[0]
            select_form.data["genome"] = [genome.id]
            session["genome"] = str(genome.id)

        # # If it is passed from page we need to reset the genome we're getting
        # if session.get('passed_from') == 'page':
        #     print ('passed from page')
        #     passed_from_page = True

        if session.get("passed_from") == "untagged":
            session["page_choice"] = 0
            session["limit_genomes"] == False

        if session.get("passed_from") == "limit_selection":
            session["page_choice"] = 0

        unique_tags = models.GenomeRecords.objects().distinct(field="tags")
        page_form.genome_tagged.choices = list(zip(unique_tags, unique_tags))

        current = models.User.objects().get(username=str(current_user.username))
        records_per_page = int(session["record_size"])

        # records_per_page = current.record_size if current.record_size != None else 20

        untagged = session["untagged"]
        limit_genomes = session["limit_genomes"]

        genome_tagged = session["genome_tagged"]

        print("ok and now is")
        print(session["show_existing_features"])
        print(session["show_stop_codons"])

        # untagged = False

        if untagged:
            if models.GenomeRecords.objects(tags=[""]).count() != 0:
                genome_count = models.GenomeRecords.objects(tags=[""]).count()
                limit_genomes = False
                session["limit_genomes"] = False

            else:
                genome_count = models.GenomeRecords.objects().count()
                session["untagged"] = False
                untagged = False

        elif limit_genomes:
            genome_count = models.GenomeRecords.objects(tags=genome_tagged[0]).count()

        else:
            genome_count = models.GenomeRecords.objects().count()

        page_count = math.ceil(genome_count / records_per_page)
        page_choices = [(x, "Page " + str(x)) for x in range(page_count)]
        page_choice = int(session["page_choice"])

        # # If we're searching directly by name, just get the genome
        if (
            session.get("passed_from") == "search_by_name"
            and genome_by_name_form.search_by_name.data
        ):
            print("search")

            print(genome_by_name_form.data)
            print(genome_by_name_form.genome_by_name.data)

            session["passed_from"] = None

            try:
                genome = models.GenomeRecords.objects.get(
                    name=genome_by_name_form.genome_by_name.data
                )
                session["genome"] = str(genome.id)
                select_form.genome.choices = [
                    (str(genome.id), genome.name + " " + genome.species)
                ]

            except:
                genome = models.GenomeRecords.objects.get(id=session["genome"])
                session["genome"] = None

                # if genome_by_name_form.genome_by_name.data:

                flash("Genome was not found in the database", category="error")

            tracks, hit_tags, genomesize = utilities.get_genome_items(
                genome,
                hits=session["hits"],
                hidden_type=session["hidden_type"],
                show_promoters=session["show_promoters"],
                show_stop_codons=bool(session["show_stop_codons"][0]),
                show_existing_features=session["show_existing_features"][0] == "True",
                checked_regions=session["checked_regions"],
            )

            associated_dict = utilities.get_associated_dict(genome)

            # session['genome'] = None

            tags = genome["tags"]

            # print('tags was ')
            #
            # print(tags)

            return self.render(
                "genomedetail.html",
                select_form=select_form,
                genome_by_name_form=genome_by_name_form,
                page_form=page_form,
                hit_form=hit_form,
                region_form=region_form,
                tracks=tracks,
                hit_tags=hit_tags,
                associated_dict=associated_dict,
                genome=genome,
                genome_name=genome.name,
                page_selected=page_choice,
                untagged=untagged,
                limit_genomes=limit_genomes,
                genome_tagged=genome_tagged[0],
                genome_tags=genome["tags"],
                hit_type=session["hits"],
                hidden_type=session["hidden_type"],
                show_promoters=session["show_promoters"],
                show_stop_codons=session["show_stop_codons"],
                show_existing_features=session["show_existing_features"],
                checked_regions=session["checked_regions"],
                genomesize=genomesize,
            )

        else:

            print("Total genomes is " + str(genome_count))
            print("Page count is " + str(page_count))
            print("Page choices is " + str(page_choices))
            print("Page choice is " + str(page_choice))
            print("untagged is " + str(untagged))
            print("limit genomes is " + str(limit_genomes))
            # print ('tags to limit to is ' + str(genome_tagged[0]))

            specific_choice = int(page_choice) * records_per_page

            if untagged:
                select_form.genome.choices = [
                    (genome.id, genome.name + " " + genome.species)
                    for genome in models.GenomeRecords.objects(tags=[""])[
                        specific_choice : specific_choice + records_per_page
                    ]
                ]

            elif limit_genomes:
                print("selecing limited genomes")
                selection = models.GenomeRecords.objects(tags=genome_tagged[0])
                print(len(selection))
                select_form.genome.choices = [
                    (genome.id, genome.name + " " + genome.species)
                    for genome in models.GenomeRecords.objects(tags=genome_tagged[0])[
                        specific_choice : specific_choice + records_per_page
                    ]
                ]

            else:

                select_form.genome.choices = [
                    (genome.id, genome.name + " " + genome.species)
                    for genome in models.GenomeRecords.objects()[
                        specific_choice : specific_choice + records_per_page
                    ]
                ]

            # The list of pages to choose from
            page_form.page.choices = page_choices

            # The page that is chosen
            page_form.data["page"] = page_choice

            # If a page change generated this request, we just want to set the genome to whatever is top of that page
            if session.get("passed_from") == "page":

                genome = models.GenomeRecords.objects()[specific_choice]

            elif (
                session.get("passed_from") == "untagged"
                or session.get("passed_from") == "limit_selection"
            ):
                if untagged:
                    genome = models.GenomeRecords.objects(tags=[""])[0]

                elif limit_genomes:
                    print("limiting genomes")
                    genome = models.GenomeRecords.objects(tags=genome_tagged[0])[0]

                else:
                    genome = models.GenomeRecords.objects()[0]

            # Otherwise we want the specific genome we've chosen
            else:
                genome = models.GenomeRecords.objects.get(id=session["genome"])

            # print ('candy')
            #
            # print (session['show_existing_features'])
            # print (bool(session['show_existing_features']))
            # print (session['show_existing_features'][0] == True)
            # print (bool(session['show_existing_features'] == 'True'))

            tracks, hit_tags, genomesize = utilities.get_genome_items(
                genome,
                hits=session["hits"],
                hidden_type=session["hidden_type"],
                show_promoters=session["show_promoters"],
                show_stop_codons=session["show_stop_codons"],
                show_existing_features=session["show_existing_features"],
                checked_regions=session["checked_regions"],
            )

            associated_dict = utilities.get_associated_dict(genome)

            return self.render(
                "genomedetail.html",
                select_form=select_form,
                genome_by_name_form=genome_by_name_form,
                page_form=page_form,
                hit_form=hit_form,
                region_form=region_form,
                tracks=tracks,
                hit_tags=hit_tags,
                associated_dict=associated_dict,
                genome=genome,
                genome_name=genome.name,
                page_selected=page_choice,
                untagged=untagged,
                limit_genomes=limit_genomes,
                genome_tagged=genome_tagged[0],
                genome_tags=genome["tags"],
                hit_type=session["hits"],
                hidden_type=session["hidden_type"],
                show_promoters=session["show_promoters"],
                show_stop_codons=session["show_stop_codons"],
                show_existing_features=session["show_existing_features"],
                checked_regions=session["checked_regions"],
                genomesize=genomesize,
            )


class UserView(ModelView):
    @login_required
    @expose("/", methods=("GET", "POST"))
    def index_view(self):
        self.edit_modal = True
        self.can_view_details = True

        self.can_create = False
        self.column_list = ["username"]

        choice = bool(random.getrandbits(1))
        print("Choice is ", choice)
        if choice:
            self.column_list = ["username"]
        else:
            self.column_list = ["username", "password"]
        return super(ModelView, self).index_view()

    @property
    def _list_columns(self):
        return self.get_list_columns()

    @_list_columns.setter
    def _list_columns(self, value):
        pass


class UserFormView(BaseView):
    @login_required
    @expose("/")
    def userform(self):
        form = forms.UserForm()
        return self.render("user.html", form=form)


class ProfileView(ModelView):
    """
    View of the Profile database for storing HMM Profiles """

    can_create = False

    @login_required
    @expose("/", methods=("GET", "POST"))
    def index_view(self):
        def form_rules():
            action_rules = []

            if not has_app_context():
                ref_names = []
            else:
                current = models.User.objects().get(username=str(current_user.username))
                ref_names = current.references

            print("ref names is ", ref_names)

            for name in ref_names:
                action_rules.append(
                    '@action("'
                    + name
                    + '_set_reference", "Set this profile as the '
                    + name
                    + ' reference")\ndef '
                    + name
                    + '_set_reference(self, ids):\n print("the name here is" + ids)'
                )

            return action_rules

        self.column_list = ("name", "references")
        # 'name', 'references', 'download')

        self.form_columns = ("name", "profile")
        self.form_extra_fields = {
            "profile": models.BlobUploadField(
                label="File",
                allowed_extensions=[
                    "pdf",
                    "doc",
                    "docx",
                    "xls",
                    "xlsx",
                    "png",
                    "jpg",
                    "jpeg",
                    "gif",
                    "hmm",
                    "fa",
                    "fasta",
                ],
                size_field="size",
                filename_field="filename",
                mimetype_field="mimetype",
            )
        }

        @action("doggydog", "Generate a doggydog")
        def kittycat(self, ids):
            print("doggydog")

        self.action_form()

        # self.list_template = 'custom_list.html'
        #
        # ar = form_rules()
        # #
        # ref_names = ['region1']
        # #
        # for name in ref_names:
        #     exec(
        #         '@action("' + name + '_set_reference", "Set this profile as the ' + name + ' reference")\ndef ' + name +
        #         '_set_reference(self, ids):\n print("the name here is" + ids)')
        #
        # print ('here is ar ')
        # print (ar)
        # for rule in ar:
        #     exec(rule)

        def on_model_change(self, form, model, is_created):
            print("change")

        def _download_formatter(self, context, model, name):

            return Markup(
                "<a href='{url}' target='_blank'>Download</a>".format(
                    url=self.get_url("download_blob", id=model.name)
                )
            )

        self.column_formatters = {
            "download": _download_formatter,
            "references": custom_filters.references_formatter,
        }

        return super(ModelView, self).index_view()

    # def is_accessible(self):
    #     if (not current_user.is_active or not
    #     current_user.is_authenticated):
    #         return False
    #     return True

    @property
    def _list_columns(self):
        return self.get_list_columns()

    @_list_columns.setter
    def _list_columns(self, value):
        pass

    @property
    def _action_form_class(self):
        return self.get_action_form()

    @_action_form_class.setter
    def _action_form_class(self, value):
        pass

    for name in ref_names:
        exec(
            '@action("'
            + name
            + '_set_reference", "Set this profile as the '
            + name
            + ' reference")\ndef '
            + name
            + '_set_reference(self, profile_ids):\n utilities.set_profile_as_reference(profile_ids, "'
            + name
            + '")'
        )


class FeatureLogView(BaseView):
    @expose("/")
    def feature_log(self):
        return self.render("features.html")


class DocumentationView(BaseView):
    @expose("/")
    def documentation(self):
        return self.render("documentation.html")


@app.route("/overview_<id>", methods=["GET"])
def download_genome_overview(id):
    """
    Route for downloading profiles from the Profile view
    :param id: Profile to download
    :return:
    """
    print(id)
    record = models.GenomeRecords.objects().get(name=id)
    return send_file(record.genome_overview, attachment_filename=id + ".png")


@app.route("/expanded_overview_<id>", methods=["GET"])
def download_genome_expanded_overview(id):
    """
    Route for downloading profiles from the Profile view
    :param id: Profile to download
    :return:
    """
    print(id)
    record = models.GenomeRecords.objects().get(name=id)
    return send_file(record.genome_expanded_overview, attachment_filename=id + ".png")


# @app.route("/<id>", methods=['GET'])
# def download_blob(id):
#     """
#     Route for downloading profiles from the Profile view
#     :param id: Profile to download
#     :return:
#     """
#     print (id)
#     profile = models.Profile.objects().get(name=id)
#     return send_file(
#         io.BytesIO(profile.profile),
#         attachment_filename=id + '.hmm')


@app.route("/genomeoverview/add_tag", methods=["GET", "POST"])
def add_tag():
    session["genome"] = request.json["genome"]

    for record in request.json["records"]:
        query = models.GenomeRecords.objects().get(
            id=record.split("/")[-1].split("_")[0]
        )
        query.tags.append(request.json["tag"])
        query.save()

    return redirect("genomeoverview")


@app.route("/delete_all_hits", methods=["GET", "POST"])
def delete_all_hits():
    queries = models.GenomeRecords.objects().all()

    queries.update(hits=[])

    models.AssociatedHits.objects().delete()

    flash("Deleted hits", category="success")

    return redirect("batch_delete")


@app.route("/delete_all_tags", methods=["GET", "POST"])
def delete_all_tags():
    print("delete_all_genome_tags")

    queries = models.GenomeRecords.objects().all()

    queries.update(tags=[])

    for query in queries:

        for hit in query.hits:
            hit.tags = []

        query.save()

    models.GenomeTags.objects().delete()

    flash("Deleted tags", category="success")

    return redirect("batch_delete")


@app.route("/genomedetail/tag_hit)", methods=["GET", "POST"])
def tag_hit():
    query = models.GenomeRecords.objects().get(id=request.json["genome"])

    print("query name")

    # print(query.id)

    print(query.name)

    # print(query.hits)

    hits = query.hits

    tag2add = request.json["tag2add"]

    print("poker")

    print(tag2add)

    print("hits")
    print(hits)

    # for hit in hits:
    #     # print(hit.id)

    for hit_id, hit_name in request.json["hits"].items():

        # Add it into the genome's list of hits

        print(hits.get(id=hit_id).tags)

        if hits.get(id=hit_id).tags == [""]:
            print("create new")
            hits.get(id=hit_id).tags = [tag2add]

        else:

            hits.get(id=hit_id).tags.append(tag2add)
        hits.save()
        formatted_hit = hit_name.replace(" ", "_").replace(":", "_")

        # At this stage, I see no reason to write out the hidden tag to the GenomeTags record
        if tag2add != "hidden":

            genome_name = (
                query.name
                + "_information_"
                + query.species.replace(" ", "_")
                + "_region_"
                + formatted_hit
            )

            if not models.GenomeTags.objects(tag_id=genome_name):
                genome_tag = models.GenomeTags(
                    tag_id=genome_name, tags=[request.json["tag2add"]]
                )
                genome_tag.save()

                print("Genome tag does not exist")

            else:

                print("Genome tag exists")
                models.GenomeTags.objects().get(tag_id=genome_name).update(
                    push__tags=request.json["tag2add"]
                )

                # genome_tag = models.GenomeTags(query.name + "_information_" + query.species.replace(" ",
                #                                                                                     "_") + "_region_" +
                #                                formatted_hit,
                #                                tag=request.json[
                #     'tag2add'])
                # genome_tag.save()

                # models.GenomeRecords.objects().get(id=request.json['genome'], hits__id=hit_id).update(push__hits__tags=
                #     request.json['tag2add'])
                #

    return redirect("genomedetail")


@app.route("/genomedetail/update_hit_tags)", methods=["GET", "POST"])
def update_hit_tags():
    genome = request.json["genome"]
    genome_name = request.json["genome_name"]
    genome_species = request.json["genome_species"]

    tags = request.json["hit_tags"].split(",")
    hit_id = request.json["hit_id"]
    hit_name = request.json["hit_name"].replace(" ", "_").replace(":", "_")

    print("hit name is ")
    print(hit_name)

    formatted_hit_name = (
        genome_name
        + "_information_"
        + genome_species.replace(" ", "_")
        + "_region_"
        + hit_name
    )

    print("hits to change is ")
    print(tags)

    print(formatted_hit_name)

    hits = models.GenomeRecords.objects().get(id=genome).hits.get(id=hit_id)

    hits.tags = tags

    hits.save()

    models.GenomeTags.objects().get(tag_id=formatted_hit_name).update(tags=tags)

    return redirect("genomedetail")


@app.route("/genomedetail/tag_genome)", methods=["GET", "POST"])
def tag_genome():
    session["genome"] = request.json["genome"]

    genome_name = request.json["genome_name"]

    tags = request.json["tag2add"].split(",")

    print("tag")

    # Don't add blank tags in
    if tags == [""]:
        tags = []

    # models.GenomeRecords.objects().get(id=request.json['genome']).update(push__tags=
    #                                                                      request.json['tag2add'])

    # models.GenomeRecords.objects().get(id=request.json['genome']).update(tags=[request.json['tag2add']])

    # Add all the tags to GenomeRecords
    models.GenomeRecords.objects().get(id=request.json["genome"]).update(tags=tags)

    # Add all the tags to GenomeTags
    # genome_tag = models.GenomeTags(tag_id=genome_name, tags=tags)
    # genome_tag.save()

    # #TODO: At the moment this just supports one tag per genome
    if not models.GenomeTags.objects(tag_id=genome_name):
        genome_tag = models.GenomeTags(tag_id=genome_name, tags=tags)
        genome_tag.save()

        # print ('Genome tag does not exist')

    else:

        # print ('Genome tag exists')
        models.GenomeTags.objects().get(tag_id=genome_name).update(tags=tags)
        models.GenomeRecords.objects().get(id=request.json["genome"]).update(tags=tags)

    return redirect("genomedetail")


@app.route("/genomedetail/clear_genome_tags)", methods=["GET", "POST"])
def clear_genome_tags():
    session["genome"] = request.json["genome"]
    genome_name = request.json["genome_name"]

    print("from here")
    print("genome name is ")

    print("clear the records")

    models.GenomeRecords.objects().get(id=request.json["genome"]).update(tags=[])

    print("clear the tags")
    models.GenomeTags.objects().get(tag_id=genome_name).delete()

    return redirect("genomedetail")


@app.route("/genomedetail/associate_hits)", methods=["GET", "POST"])
def associate_hits():
    hits = request.json["hits"]

    genome_name = request.json["genome_name"].replace(" ", "_")
    genome_id = request.json["genome"]
    genome_species = request.json["genome_species"].replace(" ", "_")

    # genome_name = models.GenomeRecords.objects().get(id=request.json['genome']).name

    print("here we be ")
    print(hits)

    print("genome name is ")
    print(genome_name)

    print("genome id is ")
    print(genome_id)

    vals = [x for x in hits.values()]

    region1 = sorted(vals)[0].replace(" ", "_").replace(":", "_")
    region2 = sorted(vals)[1].replace(" ", "_").replace(":", "_")

    print(region1)
    print(region2)

    associated_hit = models.AssociatedHits(
        genome_id,
        genome_name + "_information_" + genome_species + "_region_" + region1,
        genome_name + "_information_" + genome_species + "_region_" + region2,
    )
    associated_hit.save()

    # models.AssociatedHits.objects().get(id=request.json['genome']).update(tags=[])

    return redirect("genomedetail")


@app.route("/genomedetail/update_assoc_hits)", methods=["GET", "POST"])
def update_assoc_hits():
    print("assoco")
    print()

    remove_assoc = request.json["remove_assoc"]

    print(remove_assoc)
    print(type(remove_assoc))

    for assoc in remove_assoc:
        models.AssociatedHits.objects(id=assoc).delete()

    # models.AssociatedHits.deleteMany({'_id':{'$in':remove_assoc}})

    # models.GenomeRecords.objects().get(id=request.json['genome']).update(tags=[])
    #
    # print ('clear the tags')
    # models.GenomeTags.objects().get(tag_id=genome_name).delete()

    return redirect("genomedetail")


@app.route("/genomedetail/delete_hit", methods=["GET", "POST"])
def delete_hit():
    query = models.GenomeRecords.objects().get(id=request.json["genome"])

    print("query id")

    print(query.id)

    for hit_id in request.json["hits"].keys():
        print("hop id")
        print(hit_id)
        query.update(pull__hits__id=hit_id)

    print("before")
    print(request.json["genome"])
    print(request.json["hits"])
    print("after")

    return redirect("genomedetail")


@app.route("/genomedetail/show_hits", methods=["GET", "POST"])
def show_hits():
    # Store the information in session

    session["genome"] = request.json["genome"]
    session["hits"] = request.json["hits"]
    session["hidden_type"] = request.json["hidden_type"]
    session["page_choice"] = request.json["page_choice"]
    session["checked_regions"] = request.json["checked_regions"]
    session["untagged"] = request.json["untagged"]
    session["limit_genomes"] = request.json["limit_genomes"]
    session["genome_tagged"] = [request.json["genome_tagged"]]
    session["show_promoters"] = request.json["show_promoters"]
    session["show_stop_codons"] = (request.json["show_stop_codons"],)
    session["show_existing_features"] = (request.json["show_existing_features"],)

    session["passed_from"] = request.json["passed_from"]

    print("got to show hits")

    print(request.json["show_promoters"])

    print(request.json["show_stop_codons"])

    print(session["hidden_type"])

    return redirect("genomedetail")

    # return redirect(url_for('genomedetail.genomedetail'))


@app.route("/regions/update_outgroup_choices", methods=["GET", "POST"])
def update_outgroup_choices():
    # TODO: If we associate trees and regions with sequences correctly, we won't need to write out the alignment to
    # file and parse it back in.

    session["outgroup_choices"] = request.json["tree_choice"]

    print("session outgroup choice here is ")
    print(session["outgroup_choices"])

    return redirect("regions")


@app.route("/regions/update_regions", methods=["GET", "POST"])
def update_regions():
    keep_regions = request.json["regions"]

    models.RegionRecords.objects(name__nin=keep_regions).delete()

    return redirect("regions")


@app.route("/regions/update_region_to_profiles", methods=["GET", "POST"])
def update_region_to_profiles():
    keep_region_to_profiles = request.json["region_to_profiles"]

    print(keep_region_to_profiles)

    keep_region_to_profiles = [x.split("_")[0] for x in keep_region_to_profiles]

    print(keep_region_to_profiles)

    models.RegionToProfileRecords.objects(rtp_id__nin=keep_region_to_profiles).delete()

    return redirect("regions")


@app.route("/regions/update_region_order", methods=["GET", "POST"])
def update_region_order():
    keep_region_order = request.json["region_order"]

    keep_region_order = [x.split("_")[0] for x in keep_region_order]

    models.RegionOrderRecords.objects(name__nin=keep_region_order).delete()

    return redirect("download_region_order")


@app.route("/regions/update_aligns", methods=["GET", "POST"])
def update_aligns():
    keep_aligns = request.json["aligns"]

    models.AlignmentRecords.objects(name__nin=keep_aligns).delete()

    return redirect("regions")


@app.route("/regions/update_trees", methods=["GET", "POST"])
def update_trees():
    keep_trees = request.json["trees"]

    models.TreeRecords.objects(name__nin=keep_trees).delete()

    return redirect("regions")


@app.route("/automatic_tagging/clear_all_promoters", methods=["GET", "POST"])
def clear_all_promoters():
    utilities.clear_all_promoters()

    return redirect("automatic_tagging")


# Auto classify used temporarily to classify the 43 incorrectly tagged Simples.
# @app.route("/automatic_tagging/auto_classify", methods=['GET', 'POST'])
# def auto_classify():
#     print("Delete all existing tags")
#
#     specific_queries = ['JPPA01', 'CABPSQ01', 'FNJL01', 'NITZ01', 'RCWL01', 'NC_013216.1', 'VUOC01', 'FQYP01', 'VSRQ01',
#                         'NZ_CP010897.2', 'FNKR01', 'CABPSP01', 'BCQP01', 'QGTQ01', 'RCZD01', 'MVGR01', 'NZ_CP020080.1',
#                         'QHJL01', 'NZ_CP013047.2', 'FMVH01', 'BILZ01', 'MQZX01', 'CAAJVF01', 'PIQI01', 'NC_009253.1',
#                         'NZ_CP013461.1', 'NZ_CM001559.1', 'SMDG01', 'JAAQYP01', 'LVTS01', 'SLVP01', 'SSNZ01', 'LYRP01',
#                         'NZ_CP054613.1', 'NZ_CP038255.1', 'LMXH01', 'SAXA01', 'BJMH01', 'NZ_CP024085.1', 'PQMB01',
#                         'SMSL01', 'AEDD01', 'RKHU01']
#
#     queries = models.GenomeRecords.objects(name__in=specific_queries).timeout(False)
#
#     queries.update(tags=[])
#
#     for query in queries:
#
#         # By default we keep the tag 'hidden' as this won't get GenomeTags out of sync
#
#         for hit in query.hits:
#             if 'hidden' in hit.tags:
#                 hit.tags = ['hidden']
#             else:
#                 hit.tags = []
#
#         query.save()
#
#     models.GenomeTags.objects().all().delete()
#
#
#     print("Classifying the genomes")
#     genome_overview.classify_genomes(queries)
#
#     flash('Automatically classified')
#
#     return redirect('automatic_tagging')


# Actual auto classify
@app.route("/automatic_tagging/auto_classify", methods=["GET", "POST"])
def auto_classify():
    print("Delete all existing tags")

    utilities.delete_all_tags()
    queries = models.GenomeRecords.objects().timeout(False)

    print("Classifying the genomes")
    genome_overview.classify_genomes(queries)

    flash("Automatically classified")

    return redirect("automatic_tagging")


@app.errorhandler(404)
def page_not_found(e):
    return render_template("404.html")


@app.errorhandler(500)
def error_encountered(e):
    return render_template("500.html", error=e)


class MyModelView(ModelView):
    @property
    def _list_columns(self):
        return self.get_list_columns()

    @_list_columns.setter
    def _list_columns(self, value):
        pass

    @property
    def _action_form_class(self):
        return self.get_action_form()

    @_action_form_class.setter
    def _action_form_class(self, value):
        pass

    @property
    def column_list(self):
        if not has_app_context():

            @action("kittycat", "Generate a kittycat")
            def kittycat(self, ids):
                print("kitty cat")

            return ["team", "project_name", "approve"]
        else:

            @action("doggydog", "Generate a doggydog")
            def kittycat(self, ids):
                print("doggydog")

            return ["team", "project_name"]


admin = Admin(
    app,
    "Phylo Island",
    base_template="layout.html",
    url="/",
    template_mode="bootstrap3",
)

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", "Fields missing from ruleset", UserWarning)

    admin.add_view(UserView(model=models.User, endpoint="user"))
    admin.add_view(SetupView(name="Setup", endpoint="setup"))
    admin.add_view(UploadView(name="Upload", endpoint="upload_admin"))
    # admin.add_view(SequenceRecordsView(model=models.SequenceRecords, endpoint="sequence_records"))
    admin.add_view(
        GenomeRecordsView(model=models.GenomeRecords, endpoint="genome_records")
    )
    admin.add_view(
        ProfileView(model=models.Profile, name="Profile Records", endpoint="profiles")
    )
    # admin.add_view(MyModelView(model=models.Profile, name='Profiles', endpoint='models'))

    admin.add_view(GenomeDetailView(name="Genome Detail", endpoint="genomedetail"))
    # admin.add_view(GenomeOverviewView(name='Genome Overview', endpoint='genomeoverview'))
    admin.add_view(RegionView(name="Region Records", endpoint="regions"))
    admin.add_view(TrimRegionsView(name="Trim Regions", endpoint="trim_regions"))

    admin.add_view(ProfilesView(name="View Profiles", endpoint="region_to_profiles"))
    admin.add_view(ChartView(name="Charts", endpoint="chart"))

    admin.add_view(AlignmentsView(name="Alignments", endpoint="alignments"))

    admin.add_view(TreeView(name="Trees", endpoint="trees"))

    admin.add_view(BatchDeleteView(name="Batch Delete", endpoint="batch_delete"))

    admin.add_view(DownloadFastaView(name="Download FASTA", endpoint="download_fasta"))
    admin.add_view(
        DownloadRegionOrderView(
            name="Download Region Order", endpoint="download_region_order"
        )
    )
    admin.add_view(
        DownloadMLGOView(name="Download ML Gene Order", endpoint="download_mlgo")
    )
    admin.add_view(
        VisualiseMLGOView(name="Visualise ML Gene Order", endpoint="visualise_mlgo")
    )

    # admin.add_view(TempFixView(name='Temp Fix', endpoint='temp_fix'))
    admin.add_view(
        AutomaticTaggingView(name="Automatic Tagging", endpoint="automatic_tagging")
    )

    admin.add_view(FeatureLogView(name="Feature Log", endpoint="features"))
    admin.add_view(
        DocumentationView(name="Documentation & FAQ", endpoint="documentation")
    )
