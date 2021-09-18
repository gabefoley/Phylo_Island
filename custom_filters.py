
def seqdescription_formatter(view, context, model, name):
    """
    Format a sequence description so that you only visualise the start of a sequence
    :param view:
    :param context:
    :param model:
    :param name:
    :return:
    """
    if model.sequence:
        return model.sequence[:15] + "..."
    else:
        return model.sequence


def hitdescription_formatter(view, context, model, name):

    hit_count = {}
    output_string = ""

    if model.hits:
        for hit in model.hits:
            if hit.region in hit_count:
                hit_count[hit.region] += 1
            else:
                hit_count[hit.region] = 1

    for k, v in hit_count.items():
        output_string += "%s:%s hits " % (k, v)

    return output_string


def references_formatter(view, context, model, name):

    output_string = ""

    if model.references:

        for k in model.references.keys():
            output_string += "%s " % (k)

    return output_string