from ete3 import (
    PhyloTree,
    Tree,
    TreeStyle,
    TextFace,
    add_face_to_node,
    SeqMotifFace,
    NodeStyle,
    faces,
    ImgFace,
    CircleFace,
    AttrFace,
)
import os
import sys
import argparse
import pickle


def load_phylo_tree(tree_path, aln_path=None):
    """
    Load a tree, associate an alignment with it if given
    """
    tree = PhyloTree(tree_path, alignment=aln_path, format=1, alg_format="fasta")
    return tree


def load_tree(tree_path, aln_path=None):
    """
    Load a tree, associate an alignment with it if given
    """
    tree = Tree(tree_path, format=1)
    return tree


###
###
###


def get_region_domains(regions):
    region_domains = []

    region_colour_dict = {
        "TcdA1": "purple",
        "A1": "orange",
        "A2": "red",
        "TcB": "dodgerblue",
        "TcC": "pink",
        "Chi": "green",
        "Chi/Chi": "green",
        "TcB/TcC": "yellow",
        "A2/TcdA1": "red",
        "A1/A2/TcdA1": "red",
        "A1/TcdA1": "orange",
    }

    start = 20
    for region in regions:
        # print ("**")
        # print (region)
        region = region.replace("Chitinase", "Chi")
        region_split = region.split("_joined_")

        region_short = "/".join(
            [region.split("_")[0] for region in region_split]
        ).strip()

        if region_short == "Chitinase":
            region_name = "Chi"
        else:
            region_name = region_short

            # Override the naming of multiple hits for A region here
            if region_short in ["A2/TcdA1", "A1/A2/TcdA1", "A1/A2"]:
                region_name = "A2"
                region_short = "A2"

            if region_name in ["A1/TcdA1"]:
                region_name = "A1"
                region_short = "A1"

            if region_short in region_colour_dict:
                region_domain = [
                    start,
                    start + 120,
                    ">" if "forward" in region else "<",
                    None,
                    40,
                    "black",
                    "rgradient:" + region_colour_dict[region_short],
                    "arial|40|black|" + region_name,
                ]
            else:
                print(
                    "WARNING: Got a value for region colour dict that wasn't in there - "
                    + region_short
                )
                print(region_colour_dict.keys())
                print(region_colour_dict[region_short])
                region_domain = [
                    start,
                    start + 120,
                    ">" if "forward" in region else "<",
                    None,
                    40,
                    "black",
                    "rgradient:black",
                    "arial|40|black|" + region_name,
                ]

        region_domains.append(region_domain)
        start += 130

    return region_domains


def get_ancestor_domains(regions):
    ancestor_domains = []

    region_colour_dict = {
        "TcdA1": "purple",
        "A1": "orange",
        "A2": "red",
        "TcB": "dodgerblue",
        "TcC": "pink",
        "Chi": "green",
        "Chi/Chi": "green",
        "TcB/TcC": "yellow",
        "A2/TcdA1": "red",
        "A1/A2/TcdA1": "red",
        "A1/TcdA1": "orange",
    }

    start = 20
    for region in regions:
        # print ("**")
        print(region)
        region = region.replace("Chitinase", "Chi")

        region_name = region.replace("-", "")
        # region_split = region.split("_joined_")

        # region_short = "/".join([region.split("_")[0] for region in region_split])
        # if region_short == 'Chitinase':
        # #     region_name = 'Chi'
        # else:
        #     region_name = region_short

        if region_name in region_colour_dict:
            region_domain = [
                start,
                start + 30,
                ">" if "forward" in region else "<",
                None,
                20,
                "black",
                "rgradient:" + region_colour_dict[region_name],
                "arial|4|black|" + region_name,
            ]
        else:
            print(
                "WARNING: Got a value for region colour dict that wasn't in there - "
                + region_name
            )
            region_domain = [
                start,
                start + 30,
                ">" if "forward" in region else "<",
                None,
                20,
                "black",
                "rgradient:black",
                "arial|4|black|" + region_name,
            ]

        ancestor_domains.append(region_domain)
        start += 40

    return ancestor_domains


def get_leaf_tag(node, skip_tags=["Single", "Multiple", "Simple"]):
    for descend in node.get_descendants():
        if descend.is_leaf():
            # descend.name

            long_name = descend.name.split("_joined")[0]
            short_name = descend.name.split("_information")[0]

            tag = []

            print(long_name)
            print(short_name)
            if long_name in tag_dict:
                tag = [x for x in tag_dict[long_name] if x not in skip_tags]
            elif short_name in tag_dict:
                tag = [x for x in tag_dict[short_name] if x not in skip_tags]

            if tag:
                return tag[0]
            else:
                return "unknown"


def continue_wo_collapse(node, skip_tags=["Single", "Multiple", "Simple"]):
    """
    If there is still a node in highlight nodes that is a descendant of current node, we want to keep going
    """
    check = False

    leaves = []

    for descend in node.get_descendants():
        if descend.is_leaf():
            leaves.append(descend)

    names = [x.name for x in leaves]

    tag_set = set()

    for name in names:
        long_name = name.split("_joined")[0]
        short_name = name.split("_information")[0]

        if long_name in tag_dict:
            tag = [x for x in tag_dict[long_name] if x not in skip_tags]
        elif short_name in tag_dict:
            tag = [x for x in tag_dict[short_name] if x not in skip_tags]
        else:
            tag = ["missing"]

        if len(tag) > 1:
            print("Incorrectly tagged " + name)

        tag_set.add(tag[0])

    if len(tag_set) > 1:
        return True
    else:
        check = False

    return check


def get_domains(domains):
    pos_dict = {
        "Big_1_Full": 0,
        "VRP1_Full": 1,
        "RBD_A": 2,
        "RBD_C": 3,
        "RBD_B": 4,
        "Neuraminidase": 5,
        "TcA_RBD": 6,
        "TcB_BD_seed": 7,
    }

    domain_list = [
        # seq.start, seq.end, shape, width, height, fgcolor, bgcolor
        [10, 70, "[]", None, 20, "black", "rgradient:red", "arial|40|black|Big_1"],
        [80, 140, "[]", None, 20, "black", "rgradient:pink", "arial|40|black|VRP1"],
        [
            150,
            210,
            "[]",
            None,
            20,
            "black",
            "rgradient:lightgreen",
            "arial|40|black|RBD_A",
        ],
        [220, 280, "[]", None, 20, "black", "rgradient:blue", "arial|40|black|RBD_C"],
        [290, 350, "[]", None, 20, "black", "rgradient:orange", "arial|40|black|RBD_B"],
        [360, 430, "[]", None, 20, "black", "rgradient:purple", "arial|40|black|NMD"],
        [440, 500, "[]", None, 20, "black", "rgradient:gray", "arial|40|black|RBD_D"],
        [
            510,
            560,
            "[]",
            None,
            20,
            "black",
            "rgradient:darkgreen",
            "arial|40|black|TCB_BD",
        ],
    ]

    for k, v in pos_dict.items():
        if k not in domains:
            print(k + "was not there")
            domain_list[v][6] = "white"
            domain_list[v][7] = "arial|3|black|"

    return domain_list


def setup_domain_list(region_dict):
    region_list = []

    print("setup domain list")
    for entry in region_dict.keys():
        region_names = [
            x[0] for x in sorted(region_dict[entry].items(), key=lambda k: k[1])
        ]
        region_list.append(region_names)

    print(region_list)


def get_example_tree(
    tree,
    tag_dict,
    colour_dict,
    region_dict,
    region_order_dict,
    sequence_content_dict,
    skip_list,
    outpath,
    full_names=False,
    highlight_nodes=[],
    collapse_on_genome_tags=True,
):
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.branch_vertical_margin = 15
    ts.layout_fn = lambda x: None

    if region_dict:

        domain_list = setup_domain_list(region_dict)

    #

    # Get the colours for each extant genome
    for node in tree.iter_descendants("postorder"):
        node_style = ""

        # Force internals to be black
        if not node.is_leaf():
            colour = "black"

        if (
            collapse_on_genome_tags and continue_wo_collapse(node) == True
        ) or node.is_leaf():
            if node.is_leaf():
                node.show_leaf_name = True

                long_name = node.name.split("_joined")[0]
                short_name = node.name.split("_information")[0]

                print("raspberries")
                print(tag_dict)
                print(node.name)
                print(long_name)
                print(short_name)

                if long_name in tag_dict:

                    # print ('long name was in tag dict')
                    if tag_dict[long_name][0] in colour_dict:

                        colour = colour_dict[tag_dict[long_name][0]]
                    else:
                        colour = "black"
                        # print(tag_dict[long_name][0] + " wasn't there")

                elif short_name in tag_dict:

                    tags = [x for x in tag_dict[short_name] if x not in skip_list]

                    if len(tags) > 1:
                        print(
                            "\nWARNING: The following genome had multiple tags associated with it"
                        )
                        print(long_name)
                        print("The tags were ")
                        print(tags)
                        colour = "black"
                    else:
                        if tags[0] in colour_dict:
                            colour = colour_dict[tags[0]]

                        else:
                            "WARNING: Colour dict doesn't have an entry for " + tags[0]
                            colour = "black"

                else:
                    print(
                        "\nWARNING: We couldn't find an entry for the following genome in the tag dictionary"
                    )
                    print(node.name)
                    colour = "black"

                if full_names:
                    spaced_name = " ".join(long_name.split("_"))
                else:
                    spaced_name = " ".join(node.name.split("_")[3:5])

                plasmid_colour = "red" if "plasmid_true" in node.name else "black"

                nameFace = TextFace(
                    "  " + spaced_name, fsize=15, fgcolor=plasmid_colour
                )
                node.add_face(nameFace, column=0)

                region_names = []

                cleaned_name = (
                    node.name.replace(".", "***")
                    .replace("__", "_[", 1)
                    .replace("__", "]_", 1)
                )

                # print (cleaned_name)
                # print (region_dict.keys())
                # break

                if cleaned_name.startswith("LVTS"):
                    print("work mork")
                    print(cleaned_name)
                    print(node.name)
                # print (region_dict.keys())

                if cleaned_name in region_dict:
                    # region_names = [x for x in region_dict[cleaned_name].keys()]

                    print("HOGGLES MC BOGGLES")
                    print(cleaned_name)
                    region_names = [
                        x[0]
                        for x in sorted(
                            region_dict[cleaned_name].items(), key=lambda k: k[1]
                        )
                    ]

                    # region_dict = sorted(region_dict.items(), key=lambda k: k[1])

                if region_dict:

                    print("horseradish")
                    print(region_names)

                    box_domains = get_domains(region_names)

                    # box_domains = get_domains([x for x in region_dict[node.name.replace("***", ".")].keys()])

                    seqFace = SeqMotifFace(
                        seq=None, motifs=box_domains, gap_format="line"
                    )
                    node.add_face(seqFace, 0, "aligned")

                elif region_order_dict:

                    # Override QGTQ01 here
                    if node.name.split("_information_")[0] == "QGTQ01":
                        if "QGTQ01" in region_order_dict:
                            print("chocolate bear")
                            print(region_order_dict["QGTQ01"])

                            region_order_dict["QGTQ01"] = (
                                "TcB_expanded_forward_joined_TcC_expanded_forward, "
                                "A2_expanded_forward_joined_TcdA1_expanded_forward"
                            )

                            print(region_order_dict["QGTQ01"])

                    region_order_name = node.name.split("_information_")[0].replace(
                        ".", "***"
                    )

                    if region_order_name in region_order_dict:
                        regions = region_order_dict[region_order_name].split(",")

                        print("regions")
                        print(region_order_name)
                        print(regions)

                        region_domains = get_region_domains(regions)

                        seqFace = SeqMotifFace(
                            seq=None, motifs=region_domains, gap_format="line"
                        )
                        node.add_face(seqFace, 0, "aligned")
                elif sequence_content_dict:

                    if node.name in sequence_content_dict:
                        S = SeqMotifFace(
                            sequence_content_dict[node.name], seq_format="seq", width=6
                        )
                        node.add_face(S, 1, position="aligned")
                    else:
                        print(
                            "\nWARNING: "
                            + node.name
                            + " was not in the sequence content dict"
                        )

            else:
                colour = "black"

            if colour == None:
                colour = "black"

        else:

            if node.name not in highlight_nodes:

                if collapse_on_genome_tags:
                    node_style = "skip"
                    tag = get_leaf_tag(node)

                    format_text = (
                        " extant sequence " if len(node) == 1 else " extant sequences "
                    )
                    N = TextFace(
                        " " + str(len(node)) + format_text, fsize=14, fgcolor="black"
                    )
                    node.add_face(N, 1, position="branch-right")

                    ceiling = 50
                    length = 100

                    if display_circular or display_circular_180:
                        ceiling = 100
                        length = 800

                    wid = min(len(node) * 10, ceiling)

                    #                 node.add_face(ImgFace("./triangle.png",width=wid), 0)
                    box_domains = [
                        10,
                        length,
                        "<",
                        None,
                        wid,
                        "black",
                        "rgradient:" + colour_dict[tag],
                        "arial|30|black|" + tag,
                    ]

                    seqFace = SeqMotifFace(
                        seq=None, motifs=[box_domains], gap_format="line"
                    )
                    seqFace.rotable = False
                    node.add_face(seqFace, 0, "branch-right")
                    node.img_style["draw_descendants"] = False

        if node_style != "skip":
            nstyle = NodeStyle()
            nstyle["fgcolor"] = colour
            nstyle["size"] = 20
            node.set_style(nstyle)

    if display_circular:
        ts.mode = "c"

    if display_circular_180:
        ts.mode = "c"
        ts.arc_start = -180  # 0 degrees = 3 o'clock
        ts.arc_span = 180
    ts.root_opening_factor = 1

    return tree, ts


def get_ancestor_tree(tree, ancestral_order_dict, ref_ml_go_dict, outpath):

    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.branch_vertical_margin = 15
    ts.layout_fn = lambda x: None
    #

    region_domains = False
    ancestral_orders = False

    # Get the colours for each extant genome
    for node in tree.iter_descendants("postorder"):
        if not node.is_leaf():
            # node_style = ""
            print(node.name)
            print(ancestral_order_dict)

            if node.name:

                ancestral_orders = [
                    ref_ml_go_dict[x.replace("-", "")]
                    for x in ancestral_order_dict[node.name]
                    if x != "" and not x.startswith(">")
                ]
            print(ancestral_orders)

            if ancestral_orders:

                region_domains = get_ancestor_domains(ancestral_orders)

            if region_domains:

                seqFace = SeqMotifFace(
                    seq=None, motifs=region_domains, gap_format="line"
                )
                node.add_face(seqFace, 0, "branch-right")

    ts.root_opening_factor = 1

    return tree, ts


def colour_tips(
    tree,
    tag_dict,
    colour_dict,
    region_dict=None,
    region_order_dict=None,
    sequence_content_dict=None,
    skip_list=[],
    display_circular=False,
    display_circular_180=False,
    outpath=None,
    full_names=False,
    highlight_nodes=[],
    custom_layout=False,
    collapse_on_genome_tags=False,
):
    tree, ts = get_example_tree(
        tree,
        tag_dict,
        colour_dict,
        region_dict=region_dict,
        region_order_dict=region_order_dict,
        sequence_content_dict=sequence_content_dict,
        skip_list=skip_list,
        outpath=outpath,
        full_names=full_names,
        collapse_on_genome_tags=collapse_on_genome_tags,
    )
    # if outpath:
    # ts.layout_fn = lambda x: None

    print("outpath is ")
    print(outpath)
    tree.render(outpath, dpi=30, tree_style=ts, w=1600)


def display_ancestors(tree, ancestral_order_dict, ref_ml_go_dict, outpath):
    tree, ts = get_ancestor_tree(
        tree,
        ancestral_order_dict=ancestral_order_dict,
        ref_ml_go_dict=ref_ml_go_dict,
        outpath=outpath,
    )

    tree.render(outpath, dpi=3000, tree_style=ts, w=16000)


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--tree", help="Path to tree", required=True)
    parser.add_argument("-s", "--seqs", help="Path to sequences")
    parser.add_argument("-td", "--tag_dict", help="Path to tag dict")
    parser.add_argument("-rd", "--region_dict", help="Path to region dict")
    parser.add_argument("-rod", "--region_order_dict", help="Path to region order dict")
    parser.add_argument(
        "-scd", "--sequence_content_dict", help="Path to sequence content dict"
    )
    parser.add_argument("-cd", "--colour_dict", help="Path to colour dict")
    parser.add_argument("-o", "--outpath", help="Outpath", default="treegaze.png")
    parser.add_argument(
        "-fn", "--full_names", help="Add full name to leaf node", default=False
    )
    parser.add_argument(
        "-cgt",
        "--collapse_on_genome_tags",
        help="Collapse tree based on genome tags",
        default=False,
    )
    parser.add_argument(
        "-dc",
        "--display_circular",
        help="Display tree as circular",
        action="store_true",
    )
    parser.add_argument(
        "-dco",
        "--display_circular_180",
        help="Display tree as circular (180 degrees)",
        action="store_true",
    )
    parser.add_argument("-ao", "--ancestral_order")
    parser.add_argument("-mlgo", "--ref_ml_go_dict")

    return parser.parse_args(args)


def pickle_open(filename):
    with open(filename, "rb") as handle:
        pickle_dict = pickle.load(handle)
    return pickle_dict


if __name__ == "__main__":

    parser = parse_args(sys.argv[1:])

    loaded_tree = load_tree(parser.tree)
    outpath = parser.outpath

    if parser.ancestral_order != None:
        ancestral_order_dict = pickle_open(parser.ancestral_order)
        ref_ml_go_dict = pickle_open(parser.ref_ml_go_dict)

        display_ancestors(loaded_tree, ancestral_order_dict, ref_ml_go_dict, outpath)

    else:

        tag_dict = pickle_open(parser.tag_dict)
        # tag_dict = {"AP018269.1": ["No Hits"], "AP018270.1": ["No Hits"], "AP018271.1": ["Single", "Type3"], "AP018272.1": ["No Hits"], "AP018273.1": ["No Hits"], "CP034538.1": ["No Hits"], "DAATWL01": ["Single", "Type2b"], "OOHJ01": ["No Hits"], "NZ_CP022961.1": ["Multiple", "Type3"], "JPPA01": ["Single", "Type3"], "NUCS01": ["No Hits"], "LIVZ01": ["Single", "Type3"], "FTPJ01": ["Multiple", "Type1", "Type2b"], "LAVL01": ["No Hits"], "NC_013892.1": ["Multiple", "Type1", "Type2a"], "PCQL01": ["Single", "Type2b"], "JAAMRA01": ["Multiple", "Type2b"], "MDEO01": ["Single", "Type3"], "WIVP01": ["Single", "Type2b"], "VHLG01": ["No Hits"], "MIFU01": ["Single", "Type2b"], "NZ_AP018449.1": ["Single", "Type3"], "CABPSQ01": ["Single", "Type2b"], "CP034537.1": ["No Hits"], "PHFJ01": ["Single", "Type1"], "JAABMF01": ["Single", "Type1"], "RBQC01": ["Single", "Type2b"], "AGJN02": ["Single", "Type3"], "ONZJ01": ["Single", "Type3"], "WWHE01": ["Single", "Type3"], "FMXW01": ["Multiple", "Type2b"], "NGVR01": ["Single", "Type1"], "NZ_CP020038.1": ["No Hits"], "NZ_CP021745.1": ["No Hits"], "NZ_CP021746.1": ["No Hits"], "NZ_CP021747.1": ["No Hits"], "QGAC01": ["No Hits"], "NZ_CP024081.1": ["Single", "Type2b"], "NZ_CP015613.1": ["Single", "Type1"], "VIWL01": ["Single", "Type2b"], "NJAK01": ["Single", "Type2a"], "UUIW01": ["Single", "Type2b"], "NZ_CP012672.1": ["Multiple", "Type1", "Type3"], "NZ_CP047267.1": ["Single", "Type2b"], "PCQC01": ["Single", "Type2b"], "NHML01": ["No Hits"], "FNJL01": ["Single", "Type3"], "NZ_CP012533.1": ["Single", "Type3"], "NZ_CP012534.1": ["No Hits"], "NZ_CP012535.1": ["No Hits"], "NZ_CP012536.1": ["No Hits"], "NZ_CP012537.1": ["No Hits"], "NZ_CP012538.1": ["No Hits"], "NZ_CP012539.1": ["No Hits"], "NZ_CP012540.1": ["No Hits"], "RBSM01": ["Single", "Type2b"], "NZ_CP010029.1": ["Single", "Type2a"], "PHHE01": ["Single", "Type2b"], "NEJT01": ["Single", "Type2b"], "NZ_CP031450.1": ["Single", "Type2b"], "NZ_CP017708.1": ["No Hits"], "NZ_CP017709.1": ["No Hits"], "NZ_CP017710.1": ["No Hits"], "FYEE01": ["Multiple", "Type2b"], "ATXB01": ["Single", "Type3"], "APLI01": ["No Hits"], "AWQP01": ["Single", "Type2b"], "LMTZ01": ["No Hits"], "NCXP01": ["Single", "Type3"], "NZ_CP045799.1": ["Single", "Type2b"], "NZ_CP045800.1": ["No Hits"], "NZ_CP045801.1": ["No Hits"], "AJXJ01": ["Single", "Type2b"], "NZ_CP047073.1": ["Single", "Type2b"], "NVXX01": ["Single", "Type2b"], "MUIN01": ["Multiple", "Type2b"], "VDNE01": ["No Hits"], "LXYR01": ["No Hits"], "NZ_CP022411.1": ["Single", "Type2b"], "WIVR01": ["Multiple", "Type2b"], "WIWH01": ["Single", "Type2b"], "FONE01": ["No Hits"], "MKZS01": ["No Hits"], "QJUG01": ["Single", "Type3"], "LWBP01": ["Single", "Type3"], "QUOK01": ["Multiple", "Type3"], "FQUQ01": ["Single", "Type3"], "VIWO01": ["Multiple", "Type3"], "NITZ01": ["Single", "Type1"], "CBLV01": ["Multiple", "Type2b"], "MASH01": ["No Hits"], "LQOW01": ["No Hits"], "NEHI01": ["Single", "Type2b"], "NZ_CP038254.1": ["Single", "Type3"], "JTBY01": ["No Hits"], "FNTY01": ["Single", "Type2b"], "NZ_CP028826.1": ["Single", "Type2b"], "NIRH01": ["Single", "Type1"], "LVYD01": ["Single", "Type3"], "NZ_CP025800.1": ["Single", "Type1"], "NZ_CP025801.1": ["No Hits"], "NZ_CP025802.1": ["No Hits"], "QLKY01": ["No Hits"], "RCFQ01": ["Single", "Type3"], "SMCH01": ["Single", "Type2b"], "NZ_CP015381.1": ["No Hits"], "NMRE01": ["Single", "Type2b"], "QSNX01": ["Single", "Type2b"], "NZ_CM001558.1": ["Single", "Type2b"], "FNNQ01": ["No Hits"], "FQUS01": ["Multiple", "Type3"], "NZ_LT629762.1": ["Single", "Type2b"], "NZ_CP031065.1": ["No Hits"], "NZ_CP031066.1": ["Single", "Type3"], "NZ_CP031067.1": ["No Hits"], "VSJH01": ["Single", "Type2b"], "FNYO01": ["No Hits"], "NZ_CP036313.1": ["Single", "Type3"], "NZ_CP036314.1": ["No Hits"], "NZ_CP036315.1": ["No Hits"], "NZ_CP041668.1": ["Single", "Type3"], "NZ_CP041669.1": ["No Hits"], "FXYF01": ["Single", "Type3"], "NIRS01": ["Single", "Type2b"], "FNCO01": ["Single", "Type2b"], "RHQN01": ["Multiple", "Type1", "Type2b"], "NZ_CP021983.2": ["No Hits"], "RCWL01": ["Single", "Type3"], "QUMQ01": ["No Hits"], "QAJM01": ["Single", "Type2b"], "NBRZ01": ["Single", "Type2b"], "BJLR01": ["No Hits"], "NZ_CP039291.1": ["No Hits"], "NC_013947.1": ["Multiple", "Type3"], "JMCC02": ["Single", "Type3"], "NBVR01": ["Single", "Type1"], "QAIL01": ["Single", "Type3"], "QWFB01": ["Single", "Type2b"], "NZ_CP031648.1": ["Single", "Type2b"], "MCHY01": ["Single", "Type3"], "NZ_CP041186.1": ["Multiple", "Type3"], "NC_013216.1": ["Single", "Type3"], "JYHW01": ["Single", "Type2b"], "WIVY01": ["Multiple", "Type1", "Type2b"], "FYDX01": ["Single", "Type2b"], "MUGY01": ["No Hits"], "FNYJ01": ["No Hits"], "JJML01": ["Single", "Type3"], "FNTZ01": ["Multiple", "Type2b"], "NZ_CP029064.1": ["Single", "Type3"], "LRUN01": ["Single", "Type2b"], "VIUF01": ["Single", "Type2b"], "VZZK01": ["Single", "Type3"], "AJLJ01": ["Single", "Type3"], "CAADIW01": ["Single", "Type1"], "AXVJ01": ["No Hits"], "VIUC01": ["Multiple", "Type2b"], "AMBZ01": ["Multiple", "Type2b"], "QGGJ01": ["No Hits"], "VUOC01": ["Single", "Type3"], "QAAE01": ["Single", "Type3"], "NCWQ01": ["No Hits"], "PVTU01": ["No Hits"], "BBMZ01": ["Single", "Type2b"], "NZ_CP054043.1": ["Single", "Type1"], "SJSL01": ["Single", "Type3"], "FQYP01": ["Single", "Type3"], "NZ_CP011129.1": ["Single", "Type3"], "NC_012961.1": ["No Hits"], "NC_012962.1": ["Multiple", "Type1", "Type2a"], "BBXD01": ["No Hits"], "NZ_CP029196.1": ["Single", "Type3"], "AKJT01": ["Multiple", "Type2b"], "NVPT01": ["Single", "Type3"], "BBXG01": ["No Hits"], "ALVN01": ["Single", "Type3"], "NJFA02": ["Single", "Type1"], "NC_019738.1": ["No Hits"], "NC_019739.1": ["No Hits"], "NC_019740.1": ["No Hits"], "NC_019741.1": ["No Hits"], "NC_019742.1": ["No Hits"], "NC_019743.1": ["No Hits"], "NC_019760.1": ["No Hits"], "NC_019761.1": ["No Hits"], "NC_019762.1": ["No Hits"], "QPCD01": ["No Hits"], "QTPO01": ["Single", "Type3"], "FOEO01": ["Single", "Type2b"], "QWLL01": ["No Hits"], "QOVA01": ["Single", "Type2b"], "NZ_CP014262.1": ["Multiple", "Type2b"], "FNDJ01": ["No Hits"], "NZ_AP017422.1": ["No Hits"], "SNXZ01": ["Single", "Type3"], "FXWP01": ["Single", "Type1"], "UTBZ01": ["Multiple", "Type2b"], "BCBA01": ["Single", "Type2b"], "VSRQ01": ["No Hits"], "LFWB01": ["No Hits"], "QTUB01": ["Single", "Type2a"], "NZ_CP053584.1": ["Single", "Type2b"], "NZ_CP010897.2": ["Single", "Type2b"], "NZ_CP010898.2": ["No Hits"], "WIVZ01": ["Multiple", "Type1", "Type2b"], "NZ_CP013341.1": ["Single", "Type3"], "JACAQG01": ["Single", "Type2b"], "FNKR01": ["Single", "Type3"], "NZ_CP027723.1": ["Single", "Type2b"], "MDEN01": ["No Hits"], "CVRZ01": ["Single", "Type1"], "NZ_CP038033.1": ["No Hits"], "NZ_CP044217.1": ["No Hits"], "NZ_CP044218.1": ["Single", "Type3"], "PENV01": ["Single", "Type3"], "NRQY01": ["Single", "Type1"], "SISB01": ["Multiple", "Type2b"], "NZ_LT629732.1": ["Single", "Type3"], "AOCZ01": ["Single", "Type1"], "NZ_CP039371.1": ["Single", "Type1"], "NZ_CP039372.1": ["No Hits"], "JAFA01": ["No Hits"], "FNOY01": ["No Hits"], "CABPSP01": ["Single", "Type2b"], "LGSI01": ["Single", "Type2b"], "VZRB01": ["Single", "Type3"], "MKWS01": ["Multiple", "Type2b"], "VIUI01": ["Multiple", "Type2b"], "RXOM01": ["No Hits"], "BCQP01": ["No Hits"], "SMTE01": ["Single", "Type1"], "QMEY01": ["No Hits"], "MBDT01": ["Single", "Type2b"], "LKPJ01": ["Single", "Type3"], "OGTP01": ["Single", "Type3"], "QKTW01": ["Single", "Type3"], "NC_005773.3": ["Multiple", "Type1", "Type2b"], "NC_007274.1": ["No Hits"], "NC_007275.1": ["No Hits"], "NZ_CP048835.1": ["Single", "Type3"], "NC_010162.1": ["Multiple", "Type3"], "NEVM01": ["Single", "Type1"], "FOUX01": ["Single", "Type3"], "NZ_CP023526.1": ["No Hits"], "NZ_CP054422.1": ["Single", "Type2b"], "VOIX01": ["Single", "Type2b"], "VIWA01": ["Single", "Type3"], "VEBC01": ["No Hits"], "WIWK01": ["Multiple", "Type1", "Type2b"], "QREK01": ["Single", "Type3"], "NZ_CM002330.1": ["Single", "Type2b"], "NZ_CM002331.1": ["No Hits"], "BAHC01": ["No Hits"], "NZ_CP042968.1": ["No Hits"], "NZ_CP018049.1": ["Single", "Type2b"], "VZPM01": ["Single", "Type2b"], "QLIN01": ["Single", "Type2b"], "AUYR01": ["No Hits"], "NTYK01": ["No Hits"], "VSFF01": ["Single", "Type3"], "LRTK01": ["No Hits"], "ARBP01": ["Single", "Type3"], "ABCS01": ["Multiple", "Type3"], "BJNF01": ["Single", "Type3"], "VOQD01": ["Single", "Type3"], "VIUL01": ["Multiple", "Type2b"], "WHJD01": ["Multiple", "Type1", "Type2b"], "MLFS01": ["Single", "Type1"], "NZ_CP024900.1": ["Multiple", "Type1", "Type2a", "Type2b"], "NZ_CP009555.1": ["No Hits"], "NZ_CP009556.1": ["No Hits"], "NZ_CP013426.1": ["Single", "Type3"], "NZ_CP013427.1": ["No Hits"], "NZ_CP013428.1": ["No Hits"], "NZ_CP013429.1": ["No Hits"], "POUA01": ["No Hits"], "AJUL01": ["Single", "Type3"], "PCOS01": ["Single", "Type2b"], "QKZA01": ["Single", "Type1"], "FNQW01": ["No Hits"], "JADL01": ["No Hits"], "CABHXE01": ["Single", "Type1"], "VIKS01": ["No Hits"], "MOBX01": ["Single", "Type2b"], "QKLR01": ["Single", "Type3"], "JCLE01": ["Single", "Type1"], "FSRS01": ["Single", "Type3"], "NZ_LR134159.1": ["Multiple", "Type2b"], "VCKW01": ["Multiple", "Type3"], "WTCR01": ["Single", "Type3"], "LLWH01": ["Single", "Type2b"], "NZ_CP027738.1": ["Single", "Type3"], "QKVL01": ["Single", "Type2b"], "NZ_CP033932.1": ["Single", "Type3"], "NZ_CM001441.1": ["No Hits"], "QGTQ01": ["Single", "Type3"], "RCZD01": ["Single", "Type1"], "PYLU01": ["Single", "Type1"], "NZ_CP011288.1": ["Single", "Type2b"], "FPLG01": ["Single", "Type1"], "NZ_CP012371.1": ["No Hits"], "NZ_CP022478.1": ["No Hits"], "NMQR01": ["Multiple", "Type1", "Type2a", "Type2b"], "CTIO01": ["Single", "Type1"], "VCNG01": ["Multiple", "Type2b"], "NZ_CP007410.1": ["Multiple", "Type2b"], "NKHL01": ["No Hits"], "MVGR01": ["Single", "Type3"], "NZ_CP056779.1": ["No Hits"], "NZ_CP056780.1": ["Multiple", "Type1"], "NZ_CP056781.1": ["No Hits"], "NZ_CP056782.1": ["Single", "Type2b"], "NQKQ01": ["Single", "Type2b"], "JOGE01": ["Single", "Type3"], "NZ_CP009533.1": ["Multiple", "Type2b"], "NQKJ01": ["Multiple", "Type1", "Type2b"], "NETK01": ["Single", "Type3"], "NZ_CP031062.1": ["No Hits"], "NZ_CP031063.1": ["Single", "Type3"], "NZ_CP031064.1": ["No Hits"], "NZ_CP004078.1": ["Single", "Type3"], "PJZH01": ["No Hits"], "FNPW01": ["Multiple", "Type1"], "SEUB01": ["Multiple", "Type2b"], "UPHP01": ["Single", "Type3"], "JNGI01": ["Single", "Type3"], "UUFD01": ["No Hits"], "AAWS01": ["No Hits"], "NZ_CP021659.1": ["Multiple", "Type1"], "NZ_CP021660.1": ["No Hits"], "NZ_CP021661.1": ["No Hits"], "NZ_CP021662.1": ["No Hits"], "MOBP01": ["Single", "Type2b"], "OIFR01": ["Single", "Type3"], "JSAL01": ["Multiple", "Type2b"], "NZ_CP011104.1": ["Multiple", "Type1", "Type2b"], "MOBI01": ["Single", "Type2b"], "PUJU01": ["Multiple", "Type1", "Type2a", "Type2b"], "BIFQ01": ["Single", "Type3"], "NZ_CP025035.2": ["No Hits"], "LIUV01": ["Single", "Type2b"], "NC_010830.1": ["Single", "Type3"], "CABPSR01": ["Single", "Type2b"], "CVTM01": ["Single", "Type2b"], "RQJP01": ["Single", "Type3"], "NZ_CP009288.1": ["Single", "Type3"], "NZ_CM001025.1": ["Single", "Type2b"], "MOBT01": ["Single", "Type2b"], "NZ_LR134318.1": ["Multiple", "Type2b"], "ABBN01": ["No Hits"], "NZ_CP039287.1": ["No Hits"], "NZ_CP039288.1": ["Single", "Type3"], "NZ_CP039289.1": ["No Hits"], "LAIJ01": ["Single", "Type3"], "LFCV01": ["Single", "Type3"], "WWJN01": ["Single", "Type3"], "VZPQ01": ["No Hits"], "VOBN01": ["Single", "Type2b"], "QGTG01": ["Single", "Type3"], "AYLO01": ["Single", "Type3"], "NZ_LT707064.1": ["Single", "Type2b"], "NZ_CP020076.1": ["No Hits"], "NZ_CP020077.1": ["No Hits"], "NZ_CP020078.1": ["No Hits"], "NZ_CP020079.1": ["No Hits"], "NZ_CP020080.1": ["Single", "Type1"], "NZ_CP020081.1": ["No Hits"], "AMRI01": ["No Hits"], "NZ_LT629705.1": ["No Hits"], "NRST01": ["Single", "Type2b"], "NZ_CP050291.1": ["Single", "Type1"], "NZ_CP025263.1": ["Single", "Type2b"], "FXWM01": ["Single", "Type3"], "NZ_CP034725.1": ["Multiple", "Type2b"], "MKQR01": ["No Hits"], "FOCT01": ["No Hits"], "NUVY01": ["No Hits"], "MRCJ01": ["Single", "Type3"], "JUQG01": ["Single", "Type1"], "LECZ01": ["Single", "Type1"], "MTHI01": ["Single", "Type2b"], "NZ_CP022121.1": ["Single", "Type3"], "NZ_CM001561.1": ["Single", "Type2b"], "NZ_CP017141.1": ["Multiple", "Type3"], "AZAN01": ["Single", "Type3"], "AGFX01": ["No Hits"], "VDCQ01": ["Single", "Type3"], "QHJL01": ["Single", "Type3"], "QWEX01": ["Multiple", "Type3"], "LMCT01": ["No Hits"], "NTTM01": ["No Hits"], "VZZS01": ["Multiple", "Type1", "Type2b"], "SMFW01": ["Multiple", "Type2b"], "UEXE01": ["Single", "Type1"], "NZ_CP013046.2": ["No Hits"], "NZ_CP013047.2": ["Single", "Type1"], "FRBZ01": ["Single", "Type2b"], "AKJH01": ["Single", "Type2b"], "BBLT01": ["Single", "Type3"], "NBWC01": ["No Hits"], "NZ_CP007039.1": ["Multiple", "Type2b"], "FMCR01": ["Single", "Type3"], "VIUK01": ["Single", "Type2b"], "MVHE01": ["No Hits"], "RCOE01": ["Single", "Type2b"], "QGSY01": ["Single", "Type3"], "AXWS01": ["No Hits"], "AYMJ01": ["Single", "Type2b"], "VOBI01": ["Single", "Type2b"], "AKJK01": ["Single", "Type2b"], "FNUD01": ["Single", "Type2b"], "MOBJ01": ["Multiple", "Type2b"], "CAAKGZ01": ["Single", "Type2b"], "FOUB01": ["Single", "Type3"], "MUXN01": ["Single", "Type3"], "LKBR01": ["Single", "Type2b"], "UTVB01": ["Multiple", "Type1", "Type2b"], "PENZ01": ["Single", "Type1"], "NZ_CP009451.1": ["Single", "Type2b"], "NZ_CP034148.1": ["No Hits"], "NZ_CP034149.1": ["No Hits"], "NZ_CP034150.1": ["No Hits"], "NZ_CP034151.1": ["Single", "Type1"], "NZ_AP017313.1": ["Single", "Type3"], "FAOS01": ["No Hits"], "NZ_CP027727.1": ["Single", "Type2b"], "NZ_CP035319.1": ["No Hits"], "QAIP01": ["Single", "Type2b"], "FNBR01": ["Multiple", "Type3"], "AXDH01": ["Single", "Type1"], "FMVH01": ["Single", "Type1"], "NZ_CP036488.1": ["No Hits"], "NZ_CP036489.1": ["No Hits"], "NZ_CP036490.1": ["Single", "Type2b"], "CAQM01": ["No Hits"], "LOWA01": ["Single", "Type2b"], "NZ_CP049044.1": ["Single", "Type2b"], "NZ_CP010896.1": ["Single", "Type2b"], "NC_017168.1": ["Single", "Type1"], "NC_017169.1": ["No Hits"], "NC_017170.1": ["No Hits"], "NZ_CP031641.1": ["Single", "Type2b"], "VKDC01": ["Multiple", "Type3"], "JOAG01": ["No Hits"], "MWQG01": ["No Hits"], "VDFY01": ["Single", "Type3"], "ALVK01": ["Multiple", "Type3"], "QFRW01": ["Single", "Type3"], "BILZ01": ["Single", "Type3"], "BAXG01": ["Multiple", "Type1"], "MWPQ01": ["Single", "Type3"], "WIWM01": ["Single", "Type2b"], "FOCU01": ["Single", "Type2b"], "MQZX01": ["Single", "Type1"], "RKHS01": ["Single", "Type1"], "QHHZ01": ["No Hits"], "MYFJ01": ["Multiple", "Type1", "Type2a", "Type2b", "Type3"], "NC_016901.1": ["Single", "Type1"], "NC_016905.1": ["No Hits"], "PEIB01": ["No Hits"], "MOBQ01": ["Single", "Type2b"], "NXNJ01": ["Single", "Type2b"], "NZ_CP044407.1": ["No Hits"], "PYBV01": ["No Hits"], "JABTYG01": ["Multiple", "Type2b"], "NZ_CP042468.1": ["Multiple", "Type1", "Type3"], "NZ_CP014135.1": ["Single", "Type2b"], "NC_016818.1": ["Multiple", "Type2b"], "NC_016819.1": ["No Hits"], "NC_016835.1": ["No Hits"], "NC_017092.1": ["No Hits"], "MTAX01": ["Single", "Type3"], "NC_015559.1": ["Single", "Type3"], "LQRT01": ["No Hits"], "NZ_LS999839.1": ["Single", "Type3"], "SOCV01": ["Single", "Type2b"], "ASRX01": ["Single", "Type3"], "NZ_CP044064.1": ["Single", "Type2b"], "AKJM01": ["Single", "Type2b"], "SMKX01": ["Single", "Type3"], "CAAJVF01": ["Single", "Type3"], "VIUJ01": ["Single", "Type2b"], "LGTC01": ["No Hits"], "NZ_CP033893.1": ["No Hits"], "NZ_CP033894.1": ["Single", "Type1"], "NZ_CP033895.1": ["No Hits"], "JXRA01": ["Single", "Type3"], "RQPI01": ["Single", "Type3"], "NZ_CP023695.1": ["Single", "Type3"], "NZ_LR134335.1": ["Single", "Type1"], "SMJU01": ["Single", "Type3"], "LMCV01": ["Single", "Type2b"], "PKNM01": ["Single", "Type1"], "PIQI01": ["Single", "Type1"], "FZPH01": ["Single", "Type3"], "WIWB01": ["Single", "Type2b"], "NC_009253.1": ["Single", "Type3"], "SOZA01": ["Single", "Type2b"], "NZ_LT855380.1": ["No Hits"], "NZ_CP014947.1": ["Multiple", "Type2b"], "ALVJ01": ["No Hits"], "NZ_CP013459.1": ["No Hits"], "NZ_CP013460.1": ["No Hits"], "NZ_CP013461.1": ["Single", "Type3"], "NZ_CP048408.1": ["Multiple", "Type2b"], "NZ_CP003181.2": ["Single", "Type2b"], "VFIO01": ["Single", "Type2b"], "MASS01": ["No Hits"], "NC_020453.1": ["Single", "Type3"], "PYUC01": ["Single", "Type1"], "VEGT01": ["Single", "Type3"], "MKZO01": ["Single", "Type2b"], "WIWE01": ["Single", "Type2b"], "FMWY01": ["Single", "Type3"], "MWQL01": ["No Hits"], "FMVD01": ["No Hits"], "NZ_CP023969.1": ["Single", "Type2b"], "NZ_CP029608.1": ["Single", "Type2b"], "SMKU01": ["No Hits"], "FUKJ01": ["Single", "Type3"], "JONO01": ["Multiple", "Type2a", "Type1"], "RAVW01": ["No Hits"], "PDUD01": ["No Hits"], "MKMC01": ["Single", "Type3"], "NC_017448.1": ["Single", "Type3"], "PVZV01": ["No Hits"], "NZ_CP031069.1": ["No Hits"], "NZ_CP031070.1": ["Single", "Type2b"], "NZ_CP023269.1": ["Multiple", "Type2b"], "VLLP01": ["No Hits"], "NZ_CM001559.1": ["Single", "Type3"], "NZ_CP029983.1": ["Single", "Type2b"], "VHKL01": ["Single", "Type3"], "NZ_CP027218.1": ["Single", "Type2b"], "JPPZ01": ["Single", "Type3"], "AKJD01": ["No Hits"], "VCNJ01": ["Multiple", "Type2b"], "NZ_CP013423.1": ["No Hits"], "NZ_CP013424.1": ["No Hits"], "NZ_CP013425.1": ["No Hits"], "SMKK01": ["Single", "Type3"], "SODH01": ["No Hits"], "AZSS01": ["Multiple", "Type3"], "JFHN01": ["Single", "Type1"], "MUNM01": ["Single", "Type2b"], "NC_021492.1": ["No Hits"], "NC_021500.1": ["Single", "Type1"], "RCSU01": ["Single", "Type3"], "SMOD01": ["Single", "Type3"], "NZ_CP042382.1": ["Single", "Type3"], "NC_008268.1": ["No Hits"], "NC_008269.1": ["No Hits"], "NC_008270.1": ["No Hits"], "NC_008271.1": ["No Hits"], "JYLE01": ["Multiple", "Type1", "Type2b"], "PYMM01": ["Single", "Type1"], "NZ_CP007699.2": ["Single", "Type1"], "QAOU01": ["No Hits"], "WBOI01": ["Multiple", "Type1", "Type2b"], "CAACVJ01": ["No Hits"], "BJMN01": ["Single", "Type3"], "SMDG01": ["Single", "Type1"], "CABIWI01": ["Multiple", "Type1", "Type2b"], "WHZZ01": ["Multiple", "Type1", "Type2a", "Type2b"], "QAOQ01": ["Single", "Type3"], "RCBZ01": ["No Hits"], "NZ_CP022303.1": ["No Hits"], "NZ_CP022304.1": ["No Hits"], "NZ_CP022305.1": ["No Hits"], "NZ_CP022306.1": ["No Hits"], "AQRJ01": ["Single", "Type3"], "FNGP01": ["No Hits"], "RJKM01": ["No Hits"], "PKND01": ["Single", "Type1"], "FOSU01": ["No Hits"], "AWZT01": ["Single", "Type3"], "NZ_CP009458.1": ["Single", "Type1"], "WEGH01": ["No Hits"], "VZZZ01": ["Single", "Type3"], "NZ_CP043060.1": ["Multiple", "Type2b"], "VZZR01": ["No Hits"], "JAAQYP01": ["Single", "Type3"], "MTBD01": ["Single", "Type1"], "NZ_CP019686.1": ["Single", "Type1"], "VUAZ01": ["Single", "Type2b"], "AZXK01": ["No Hits"], "BBIR01": ["Multiple", "Type2b"], "MWLO01": ["No Hits"], "QYZD01": ["Single", "Type2b"], "MOBZ01": ["Single", "Type2b"], "FOBB01": ["Single", "Type3"], "FMDM01": ["Single", "Type3"], "NZ_CP019888.1": ["Single", "Type1"], "AUAX01": ["Single", "Type3"], "AVEF02": ["Single", "Type1"], "FNAD01": ["Single", "Type3"], "BBCC01": ["Single", "Type1"], "QAOV01": ["Single", "Type2b"], "BAZX01": ["No Hits"], "NKQZ01": ["Single", "Type3"], "NZ_CP022960.1": ["Single", "Type2b"], "SHKK01": ["Single", "Type3"], "NZ_CP012673.1": ["Multiple", "Type1", "Type3"], "WBKQ01": ["Single", "Type3"], "LQAL01": ["Single", "Type2b"], "FCNY02": ["Single", "Type3"], "VFEU01": ["Single", "Type2b"], "SHKT01": ["No Hits"], "BAVR01": ["No Hits"], "OAOQ01": ["Single", "Type3"], "MSSW01": ["Single", "Type3"], "JOGR01": ["Single", "Type3"], "NZ_CP034337.1": ["No Hits"], "CAACUY01": ["No Hits"], "QVNU01": ["Single", "Type3"], "ASSC01": ["Single", "Type3"], "NZ_CP028035.1": ["Single", "Type2b"], "NZ_CP028036.1": ["No Hits"], "NZ_CP028037.1": ["No Hits"], "NZ_CP028038.1": ["No Hits"], "NZ_CP028039.1": ["No Hits"], "NZ_CP028040.1": ["No Hits"], "NZ_CP028041.1": ["No Hits"], "NZ_CP028042.1": ["No Hits"], "NZ_CP057330.1": ["Single", "Type2b"], "NZ_CP057331.1": ["No Hits"], "NZ_CP057332.1": ["No Hits"], "NZ_CP057333.1": ["No Hits"], "NZ_CP046052.1": ["No Hits"], "NZ_CP046053.1": ["No Hits"], "NZ_CP046054.1": ["No Hits"], "BFCB01": ["Single", "Type3"], "NZ_CM001489.1": ["No Hits"], "RBOV01": ["Single", "Type2b"], "NZ_CP030750.1": ["Multiple", "Type1", "Type2b"], "SODV01": ["Single", "Type3"], "QEKL01": ["Multiple", "Type2b"], "QJRT01": ["Single", "Type2b"], "NZ_LT629778.1": ["Single", "Type2b"], "NZ_CP024634.1": ["No Hits"], "LVTS01": ["Single", "Type1"], "NZ_AP018150.1": ["Single", "Type2b"], "CPYD01": ["Single", "Type2a"], "RIAR02": ["Multiple", "Type3"], "NZ_CP010407.1": ["Multiple", "Type3"], "NZ_CP010408.1": ["No Hits"], "AUKP01": ["Single", "Type3"], "MTSA01": ["Single", "Type2b"], "MOBL01": ["Single", "Type2b"], "NZ_CP028158.1": ["No Hits"], "MOBY01": ["Multiple", "Type2b"], "NIBV01": ["Multiple", "Type1", "Type2a"], "SJOP01": ["Single", "Type1"], "VEJO01": ["Single", "Type2b"], "FWXB01": ["No Hits"], "SMKT01": ["Single", "Type3"], "FPBO01": ["No Hits"], "NZ_CP012831.1": ["Multiple", "Type2b"], "JAAQYX01": ["Multiple", "Type2b"], "QGSZ01": ["Single", "Type3"], "QARA01": ["Single", "Type2b"], "RXLQ01": ["Single", "Type3"], "ANOR01": ["Single", "Type2b"], "PJBP01": ["Single", "Type3"], "NQKI01": ["Single", "Type2b"], "CBSW01": ["Multiple", "Type1", "Type2a", "Type2b"], "VCKX01": ["No Hits"], "FNGF01": ["Single", "Type3"], "NZ_CP010310.2": ["No Hits"], "WJIE01": ["Multiple", "Type1"], "VAUR01": ["Multiple", "Type2b"], "NC_017956.1": ["No Hits"], "NC_017957.2": ["No Hits"], "NC_017958.1": ["No Hits"], "NC_017959.1": ["No Hits"], "NC_017966.1": ["No Hits"], "NZ_CP013368.1": ["No Hits"], "NZ_CP013369.1": ["Single", "Type2b"], "NZ_CP013370.1": ["No Hits"], "NZ_CP013371.1": ["No Hits"], "VFVY01": ["No Hits"], "NZ_CP023567.1": ["Single", "Type2b"], "NZ_CP023568.1": ["No Hits"], "RJKL01": ["No Hits"], "NZ_CP045767.1": ["Single", "Type2b"], "JOGP01": ["Single", "Type3"], "CABHPT01": ["Single", "Type1"], "UWFA01": ["No Hits"], "QVNQ01": ["Single", "Type3"], "PJBC01": ["Single", "Type2b"], "SHKZ01": ["Single", "Type3"], "JQFM01": ["No Hits"], "SLVP01": ["Single", "Type3"], "NZ_CP005969.1": ["Multiple", "Type1", "Type2b"], "NZ_CP021106.3": ["No Hits"], "JYLH01": ["Single", "Type2b"], "LJCS01": ["Multiple", "Type1"], "NZ_CP017606.1": ["No Hits"], "NZ_CP017607.1": ["Single", "Type1"], "NZ_CP017608.1": ["No Hits"], "NZ_CP017609.1": ["No Hits"], "SOAB01": ["Multiple", "Type3"], "QKLY01": ["Single", "Type2b"], "JAAQWH01": ["Single", "Type2b"], "AKJS01": ["Single", "Type2b"], "QPDS01": ["Multiple", "Type2b"], "NZ_CP023965.1": ["Single", "Type1"], "NZ_CP029231.1": ["No Hits"], "NZ_CP029232.1": ["No Hits"], "NZ_CP029233.1": ["No Hits"], "NZ_CP029234.1": ["No Hits"], "NZ_CP029235.1": ["Single", "Type3"], "LNCD01": ["Single", "Type3"], "NVDM01": ["Single", "Type1"], "NZ_CP011020.1": ["Multiple", "Type2b"], "PZZQ01": ["Multiple", "Type1"], "NZ_CP011807.3": ["Single", "Type2b"], "NZ_CP011808.2": ["No Hits"], "NZ_CP011809.2": ["No Hits"], "ALJC01": ["Single", "Type2b"], "QJRQ01": ["No Hits"], "QEOK01": ["Single", "Type3"], "NZ_CP029618.1": ["Single", "Type3"], "NZ_CP010016.1": ["No Hits"], "NZ_CP010017.1": ["No Hits"], "NZ_CP010018.1": ["No Hits"], "SOCQ01": ["Single", "Type2b"], "RJKE01": ["Multiple", "Type3"], "QAOO01": ["No Hits"], "JMCL01": ["Multiple", "Type2b"], "QBKC01": ["No Hits"], "NZ_CP034335.1": ["No Hits"], "VDLX02": ["No Hits"], "SSMR01": ["Multiple", "Type3"], "NSCM01": ["Multiple", "Type1", "Type2a", "Type2b"], "VMSG01": ["Multiple", "Type2b"], "ABCQ01": ["Single", "Type1"], "OUND01": ["No Hits"], "QAJI01": ["Single", "Type3"], "NZ_CP045761.1": ["Single", "Type2b"], "SUNB01": ["No Hits"], "LJSD01": ["No Hits"], "NZ_CP041692.1": ["Single", "Type3"], "PZZR01": ["No Hits"], "JPMW01": ["No Hits"], "QPIJ01": ["Single", "Type3"], "LYUY01": ["No Hits"], "SMJX01": ["No Hits"], "VATL01": ["Single", "Type1"], "FMYF01": ["Single", "Type3"], "PUWV01": ["Multiple", "Type1", "Type2a", "Type2b", "Type3"], "FQYM01": ["No Hits"], "PJRP01": ["No Hits"], "QRBE01": ["No Hits"], "FOVJ01": ["Multiple", "Type3"], "SOCT01": ["No Hits"], "CABMLW01": ["Multiple", "Type1"], "BDBY01": ["Single", "Type3"], "PYGV01": ["No Hits"], "VRLS01": ["Single", "Type1"], "ASTJ01": ["Single", "Type3"], "LVEJ01": ["Single", "Type2b"], "OUNR01": ["Single", "Type3"], "FPBP01": ["Single", "Type3"], "FSRU01": ["Single", "Type2b"], "SMKN01": ["No Hits"], "ASJB01": ["Single", "Type3"], "VIYH01": ["Single", "Type2b"], "SNZP01": ["Single", "Type3"], "NZ_CP014847.1": ["No Hits"], "NZ_CP014848.1": ["No Hits"], "NZ_CP014849.1": ["No Hits"], "NZ_CP014850.1": ["No Hits"], "NZ_CP014851.1": ["No Hits"], "NZ_CP014852.1": ["No Hits"], "NZ_CP014853.1": ["Single", "Type2b"], "NZ_CP034780.1": ["No Hits"], "NZ_CP034781.1": ["No Hits"], "NZ_CP034782.1": ["No Hits"], "NZ_CP034783.1": ["Single", "Type3"], "PYBJ01": ["Single", "Type3"], "PTJB01": ["No Hits"], "NZ_CP024159.1": ["No Hits"], "JNYY01": ["Single", "Type3"], "NZ_CP027756.1": ["Single", "Type2b"], "SSNZ01": ["Single", "Type3"], "NZ_CP046874.1": ["Multiple", "Type2b"], "WIBD01": ["Single", "Type2b"], "NZ_CP029710.1": ["Single", "Type3"], "RBRE01": ["Multiple", "Type1", "Type2b"], "NZ_CP024866.1": ["Single", "Type2b"], "JAAAGD01": ["Single", "Type2b"], "JAAEFD01": ["Single", "Type1"], "RBUY01": ["No Hits"], "QXQA01": ["Single", "Type3"], "QJRP01": ["No Hits"], "AXBA01": ["Single", "Type3"], "OMPE01": ["No Hits"], "NZ_LT629790.1": ["Multiple", "Type2b"], "LLWI01": ["Multiple", "Type1", "Type2b"], "NZ_LT629746.1": ["Multiple", "Type2b"], "BAOS01": ["Single", "Type1"], "VLPL01": ["Single", "Type1"], "LYRP01": ["Single", "Type1"], "JXDG01": ["Multiple", "Type2b"], "LIPP01": ["Single", "Type3"], "JAAQWI01": ["Multiple", "Type2b"], "NZ_LT629795.1": ["Single", "Type2b"], "LXEN01": ["Single", "Type1"], "NZ_CM001514.1": ["Multiple", "Type2b"], "NC_015731.1": ["Single", "Type3"], "NZ_CP041754.1": ["Multiple", "Type1", "Type2b"], "NZ_CP041755.1": ["No Hits"], "NZ_CP041756.1": ["No Hits"], "NZ_CP045572.1": ["No Hits"], "FOVO01": ["No Hits"], "JMKF01": ["No Hits"], "JACCAY01": ["Single", "Type3"], "MAID01": ["No Hits"], "NZ_CP009727.1": ["No Hits"], "NZ_CP009728.1": ["No Hits"], "AKXV01": ["Single", "Type3"], "CZQA01": ["Single", "Type3"], "NZ_CP054609.1": ["No Hits"], "NZ_CP054610.1": ["No Hits"], "NZ_CP054611.1": ["No Hits"], "NZ_CP054612.1": ["No Hits"], "NZ_CP054613.1": ["Single", "Type3"], "NZ_CP014226.1": ["Single", "Type3"], "CABIVM01": ["Single", "Type2b"], "CDPK01": ["Single", "Type1"], "WIAO01": ["Single", "Type3"], "WSFG01": ["Multiple", "Type1", "Type2a", "Type2b"], "WIVV01": ["Multiple", "Type1", "Type2b"], "NQKN01": ["Single", "Type2b"], "FUWU01": ["No Hits"], "UTLV01": ["Single", "Type2b"], "QOIO01": ["Single", "Type3"], "QGLF01": ["No Hits"], "VSRO01": ["Multiple", "Type2b"], "WIVQ01": ["Single", "Type2b"], "OBDY01": ["No Hits"], "SJZK02": ["Single", "Type1"], "NZ_CP029482.1": ["Single", "Type2b"], "LHVN01": ["Single", "Type2b"], "UEXF01": ["Single", "Type1"], "VIUR01": ["Multiple", "Type2b"], "QBJA02": ["Multiple", "Type1", "Type2b"], "LZEX01": ["Multiple", "Type1", "Type2b"], "JYLN01": ["Single", "Type2b"], "NZ_CP022504.1": ["Single", "Type1"], "LFQK01": ["Single", "Type2b"], "NZ_CP038255.1": ["Single", "Type3"], "NZ_CP024646.1": ["Single", "Type2b"], "CVJX01": ["Single", "Type1"], "NZ_CP031146.1": ["Single", "Type1"], "LACH01": ["Single", "Type2b"], "NZ_CP011253.3": ["Single", "Type2b"], "NZ_CP011518.2": ["No Hits"], "NZ_CP011519.2": ["No Hits"], "VJZE01": ["No Hits"], "QKYK01": ["No Hits"], "NZ_CP023525.1": ["Single", "Type1"], "AZUB01": ["Single", "Type2b"], "JRYA01": ["Multiple", "Type1", "Type2b"], "AUEZ01": ["Single", "Type3"], "VSFG01": ["Single", "Type3"], "QOIL01": ["No Hits"], "SJZD01": ["No Hits"], "VJWF01": ["Single", "Type3"], "FNSJ01": ["Multiple", "Type3"], "AWNH01": ["No Hits"], "NZ_CP016211.1": ["No Hits"], "SMKO01": ["No Hits"], "VXLC01": ["Single", "Type3"], "JFCA01": ["Single", "Type2b"], "PPRY01": ["Single", "Type2b"], "NZ_LT629708.1": ["Single", "Type2b"], "NZ_CP026364.1": ["Single", "Type1"], "NZ_CP015225.1": ["Single", "Type2b"], "MUNY01": ["Single", "Type3"], "MBLO01": ["Multiple", "Type3"], "NZ_CP012159.1": ["Single", "Type3"], "ALBT01": ["Single", "Type3"], "RPND01": ["Single", "Type3"], "LMXH01": ["Single", "Type3"], "CPZI01": ["Single", "Type2b"], "SAXA01": ["Single", "Type3"], "QPAO01": ["Single", "Type1"], "SSMQ01": ["Multiple", "Type3"], "PQKR01": ["Single", "Type2b"], "NC_015510.1": ["Single", "Type3"], "CWJL01": ["Single", "Type1"], "BJMH01": ["Single", "Type3"], "JABWHT01": ["Single", "Type1"], "NZ_CP024085.1": ["Single", "Type3"], "SMFY01": ["No Hits"], "NZ_CP027732.1": ["Single", "Type2b"], "NQKV01": ["Multiple", "Type2b"], "PYGG01": ["No Hits"], "AWXZ01": ["Single", "Type3"], "NZ_LN681227.1": ["Multiple", "Type1", "Type2a", "Type2b"], "PQMB01": ["Single", "Type1"], "FODH01": ["Single", "Type3"], "FOIG01": ["Single", "Type3"], "ARBJ01": ["Single", "Type3"], "NZ_CP029843.1": ["Single", "Type3"], "FOLC01": ["Single", "Type3"], "NZ_LT629691.1": ["Multiple", "Type1", "Type2b"], "FORG01": ["Multiple", "Type1"], "WIWO01": ["Single", "Type2b"], "MUBJ01": ["Multiple", "Type2a", "Type2b"], "QPIP01": ["Single", "Type1"], "QNVV01": ["Single", "Type3"], "LHVM01": ["Single", "Type2b"], "PVZG01": ["Single", "Type3"], "NZ_CP017687.1": ["Single", "Type2b"], "FNQB01": ["Multiple", "Type3"], "NQYG01": ["Single", "Type2b"], "JQNB01": ["No Hits"], "FNVO01": ["No Hits"], "VIVV01": ["Single", "Type3"], "CABLBP01": ["Single", "Type3"], "NZ_CP043422.1": ["Single", "Type1"], "QEOF01": ["Single", "Type2b"], "PVNL01": ["Single", "Type3"], "CQAW01": ["Single", "Type1"], "AUAF01": ["No Hits"], "SZQA01": ["No Hits"], "NC_015379.1": ["Multiple", "Type2b"], "NZ_CP026106.1": ["Single", "Type3"], "VIUG01": ["Multiple", "Type2b"], "SMSC01": ["Single", "Type1"], "MQWC01": ["Single", "Type3"], "NZ_CP042804.1": ["No Hits"], "PENT01": ["Single", "Type1"], "NC_014718.1": ["Multiple", "Type2b"], "NC_014722.1": ["No Hits"], "PENW01": ["Single", "Type1"], "QRAV01": ["Single", "Type2b"], "QEQQ01": ["Single", "Type2b"], "MKCS01": ["Single", "Type3"], "JPIX01": ["Single", "Type1"], "NJAH01": ["No Hits"], "JPLA01": ["No Hits"], "VFPP01": ["No Hits"], "BBXC01": ["No Hits"], "JPOH01": ["Single", "Type3"], "VWSH01": ["Single", "Type3"], "FNON01": ["Single", "Type3"], "FAOZ01": ["Single", "Type3"], "NC_013131.1": ["Single", "Type3"], "FOOS01": ["Single", "Type3"], "NZ_CP053682.1": ["Single", "Type1"], "NZ_CP013949.1": ["Single", "Type3"], "VZPJ01": ["Multiple", "Type2b"], "MIIW01": ["Single", "Type1"], "JXDI01": ["Single", "Type2b"], "NZ_LT629704.1": ["Single", "Type2b"], "NZ_CP016634.1": ["Multiple", "Type1", "Type2b"], "SLYK01": ["No Hits"], "NZ_AP020337.1": ["Multiple", "Type1", "Type2b"], "PYAL01": ["Single", "Type2b"], "FYDV01": ["Single", "Type2b"], "VCBA01": ["No Hits"], "QFZQ01": ["Single", "Type3"], "SMSL01": ["Single", "Type3"], "FODZ01": ["Single", "Type3"], "PYGD01": ["Multiple", "Type3"], "LHVL01": ["Multiple", "Type2b"], "QAIK01": ["Single", "Type2b"], "FODL01": ["Single", "Type2b"], "NZ_LT828648.1": ["Single", "Type3"], "FUXU01": ["Single", "Type1"], "VZII01": ["Single", "Type2b"], "QLLL01": ["Single", "Type3"], "NZ_CP016176.1": ["No Hits"], "FQWR01": ["Multiple", "Type3"], "NZ_CP029197.1": ["Single", "Type3"], "VIUE01": ["Multiple", "Type2b"], "VIVC01": ["Single", "Type2b"], "NEJJ01": ["Single", "Type2b"], "VOQB01": ["Multiple", "Type1"], "NC_007954.1": ["No Hits"], "AEDD01": ["Single", "Type3"], "NC_020209.1": ["Single", "Type2b"], "LKBT01": ["Single", "Type2b"], "NZ_LT629801.1": ["Single", "Type2b"], "JHVK01": ["Single", "Type2b"], "PVNG01": ["No Hits"], "RBIS01": ["No Hits"], "FNTT01": ["Single", "Type2b"], "VCRA01": ["No Hits"], "ABXF01": ["Single", "Type3"], "PJCP01": ["Single", "Type2b"], "NZ_CP045158.1": ["Single", "Type1"], "BBXE01": ["No Hits"], "RKHU01": ["Single", "Type3"], "LIPN01": ["Single", "Type3"], "NZ_LT629788.1": ["Single", "Type2b"], "NIBS01": ["Single", "Type2a"], "JAABNH01": ["Single", "Type1"], "NZ_CP017599.1": ["Multiple", "Type3"], "QKWJ01": ["Single", "Type3"], "FNVU01": ["Single", "Type3"], "RAVX01": ["Single", "Type3"], "NZ_CP024923.1": ["Single", "Type3"], "NKFP01": ["No Hits"], "AEDB02": ["No Hits"], "NZ_CP038613.1": ["No Hits"], "NZ_CP018319.1": ["Multiple", "Type2b"], "SLZA01": ["Single", "Type3"], "WIWL01": ["Single", "Type2b"], "SMFY01_information_Ancylobacter_aquaticus_region_TcB_expanded_[640.9]_4636162_4644109_backward": [""], "SMFY01_information_Ancylobacter_aquaticus_region_A2_expanded_[100.2]_4628824_4636126_backward": [""], "SMFY01_information_Ancylobacter_aquaticus_region_TcC_expanded_[251.6]_4636162_4644109_backward": [""], "QTSW01": ["Single", "Type2b"], "VIVY01": ["No Hits"], "MOBO01": ["Multiple", "Type2b"], "NZ_CP013997.1": ["No Hits"], "RBIO01": ["No Hits"], "SNYE01": ["Single", "Type3"], "NZ_CP038438.1": ["Single", "Type2b"], "SMAS01": ["Single", "Type1"]}
        region_dict = pickle_open(parser.region_dict)
        region_order_dict = pickle_open(parser.region_order_dict)
        sequence_content_dict = pickle_open(parser.sequence_content_dict)

        colour_dict = pickle_open(parser.colour_dict)

        full_names = True if parser.full_names == "True" else False

        collapse_on_genome_tags = (
            True if parser.collapse_on_genome_tags == "True" else False
        )

        display_circular = True if parser.display_circular else False
        display_circular_180 = True if parser.display_circular_180 else False
        skip_list = [
            "Single",
            "Multiple",
            "Simple",
        ]  # Tags to skip when looking for tags to colour on

        print("\nMake tree called with the following dictionaries - ")

        print("tag dict")
        print(tag_dict)
        print("region dict")
        print(region_dict)

        print("region order dict")

        print(region_order_dict)
        print("sequence_content_dict")

        print(sequence_content_dict)
        print("colour dict")
        print(colour_dict)
        print("skip list")
        print(skip_list)

        # Override the tag dict so that Type3_A2 collapses with Type3 and Type1_A2 collapses with Type1
        for name, tags in tag_dict.items():
            if "Type3_A2" in tags:
                tag_dict[name] = [x.replace("Type3_A2", "Type3") for x in tags]

            if "Type1_A2" in tags:
                tag_dict[name] = [x.replace("Type1_A2", "Type1") for x in tags]

            if "Type2a_force" in tags:
                tag_dict[name] = [x.replace("Type2a_force", "Type2a") for x in tags]

        colour_tips(
            loaded_tree,
            tag_dict,
            colour_dict,
            region_dict,
            region_order_dict,
            sequence_content_dict,
            skip_list,
            display_circular,
            display_circular_180,
            outpath,
            full_names=full_names,
            collapse_on_genome_tags=collapse_on_genome_tags,
        )
