jQuery(document).ready(function () {
    // delete selection
    jQuery('#btn_del').on('click', function () {
        console.log('was clicked')
        var ngy2data = $("#my_nanogallery2").nanogallery2('data');
        ngy2data.items.forEach(function (item) {
            if (item.selected) {
                console.log(item)
            }
        });
        $("#my_nanogallery2").nanogallery2('resize');
    });

    // switch selection mode on/off
    jQuery('#btn_select_mode').on('click', function () {
        var b = !$('#my_nanogallery2').nanogallery2('option', 'thumbnailSelectable');
        $('#my_nanogallery2').nanogallery2('option', 'thumbnailSelectable', b);
    });

    jQuery('#btnsearch2').on('click', function() {
        jQuery("#my_nanogallery2").nanogallery2('search2', jQuery('#ngy2search').val(), jQuery('#ngy2search2').val() );
        alert( 'found: ' + jQuery("#my_nanogallery2").nanogallery2('search2Execute') );
    });

    jQuery('#tagBtn').on('click', function () {
        var ngy2data = $("#my_nanogallery2").nanogallery2('data');
        tag = jQuery('#tag2add').val()
        records = [];

        ngy2data.items.forEach(function (item) {
            if (item.selected) {
                records.push(item.src)
            }
        });

        console.log (tag)
        console.log(records)
        $.ajax({
            url: Flask.url_for("add_tag"),
            data: JSON.stringify({
                'tag' : tag, 'records' : records
            }),
            contentType: 'application/json;charset=UTF-8',
            type: 'POST',
            error: function(error) {
                console.log(error);
            },
                          success: function (response) {


                        window.location.reload();


                    }
        });

        $("#my_nanogallery2").nanogallery2('resize');
    });
});