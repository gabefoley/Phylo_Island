/* Initialize the layout for the linear plot, the base
   parameters such as the genome size, the height of the
   plot in px, the width of the plot, the div container to
   put the SVG element, and an initial zoom level */

var selected = new Object()
var linearlayout = { genomesize: genomesize,
                     height: 300,
                     width: 500,
                     container: "#linearchart",
                     initStart: 0,
                     initEnd: 200000,
		    };

/* Initialize the layout for the linear brush, including
   what container to put the brush in to */
var contextLayout = { genomesize: 6264404,
                      container: "#brush" };

/* Initialize the plot with the layout and the tracks data,
   then initialize the brush telling it the plot it needs
   to update */
var linearTrack = new genomeTrack(linearlayout, tracks);
var brush = new linearBrush(contextLayout,linearTrack);

/* Tell the linear plot about the brush so when the linear
   track is updated (via a drag event, etc) it will update
   the brush */
linearTrack.addBrushCallback(brush);

/* If we have a circular plot, tell the linear plot
   to let it know about updates, this is needed for
   the combo demo */
window.onload = function() {
  if('undefined' !== typeof cTrack) {
    console.log("Hooking up circular plot callback");
    linearTrack.addBrushCallback(cTrack);
  }
}

/* Callback to demo resizing the linear plot */
function resizeLinearPlot() {
    linearTrack.resize(1000);
}

/* Catch a click callback on the linear plot and show what
   information we're given about the track item */
function linearPopup(trackName, d) {
    console.log(d);
    alert("Received click event from track " + trackName + ", item: " + JSON.stringify(d));
}

/* Callback to demo click functionality on the linear plot,
   the click callback is defined in the data json object for
   each track */
function linearClick(trackName, d) {

    var new_id = d.id;
    console.log('plop')
    console.log(d)
    var new_elem = d.name + " " + d.start + ":" + d.end + " " + d.actual_strand;

    var selected_vals = new Array()

    console.log('selected is ')
    console.log(selected)



    console.log(d)

    console.log('new id is ' + new_id)


    if (selected.hasOwnProperty(new_id)) {
        console.log('it was in it')
        delete selected[new_id]
    }
    else {
        console.log('was not in it')
        selected[new_id] = new_elem
    }


    var list = document.getElementById('myList');
    console.log('selected is now')
    console.log(selected)
    console.log('and the keys')
    console.log(Object.keys(selected))

    for(var index in selected) {
        selected_vals.push(selected[index])
    }

    // console.log('and the selected arr')
    //
    // console.log(selected_arr)
    $('#myList').html(selected_vals.join('<br>'));

}
