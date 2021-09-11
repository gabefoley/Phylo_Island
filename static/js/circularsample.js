// Sample routines to initialize the plot for rendering
// and add some basic functionality.

// Initialize the defaults for the chart such as
// the genome size, the div container to put the
// SVG object in, what function to call during a
// double click and the initial chart size.
var circularlayout = {genomesize: genomesize,
		      container: "#circularchart",
		      dblclick: "doubleClick",
                      w: 550, h: 550
        };


track_names = ['A1', 'A2', 'TcdA1', 'TcB', 'TcC', 'Chitinase']

expanded_track_names = ['A1_expanded', 'A2_expanded', 'TcdA1_expanded', 'TcB_expanded', 'TcC_expanded', 'Chitinase_expanded']

// The actual initialization call which takes two
// parameters, the layout (above) for the plot and
// the dataset to visualize (from data.js, a json
// data structure)
var cTrack = new circularTrack(circularlayout, tracks);

// If we're showing both a circular and linear chart,
// and have a linear brush, attach it (see combo plot demo)
if('undefined' !== typeof linearTrack) {
    console.log("Attaching linear track");
    cTrack.attachBrush(linearTrack);
    cTrack.showBrush();
}

if('undefined' !== typeof brush) {
    console.log("Attaching linear track brush");
    cTrack.attachBrush(brush);
}

// Now some callbacks to make the interactive functionality work.

// Attached to the onchange callback for the GC Plot checkbox,
// call the plot to add/remove the GC Plot as needed
function updateGC(cb) {
    if(cb.checked) {
	cTrack.showTrack("gcplot");
    } else {
	cTrack.hideTrack("gcplot");
    }
}

// Attached to strand track checkbox, call the plot to
// add/remove the inner stranded track
function updateStrand(cb) {
    if(cb.checked) {
	cTrack.showTrack("track1");
    } else {
	cTrack.hideTrack("track1");
    }
}

function updateHits(cb){

    console.log('got to uh')
    if(cb.checked) {
        for (track in track_names) {
            cTrack.showTrack(track);
        }
    }
    else
    {
        for (track in track_names) {

            cTrack.hideTrack(track);
        }

    }

}

function updateExpandedHits(cb){

    console.log("got ot ueh")
    if(cb.checked) {
        for (track in expanded_track_names) {
            cTrack.showTrack(track);
        }
    }
    else
    {
        for (track in expanded_track_names) {

            cTrack.hideTrack(track);
        }

    }

}

// Attached to the contig gap checkbox, call the plot to
// add/remove the contig gap squiggles
function updateGaps(cb) {
    if(cb.checked) {
	cTrack.showTrack("gapTrack");
    } else {
	cTrack.hideTrack("gapTrack");
    }
}

// Attached to the ADB glyph checkbox, call the plot to
// add/remove only the ADB type of glyph
function updateAdb(cb) {
    if(cb.checked) {
	cTrack.showGlyphTrackType("track5", "adb");
    } else {
	cTrack.hideGlyphTrackType("track5", "adb");
    }
}

// Attached to the resize plot button, call the plot to
// resize the plot to 650px diameter
function resizePlot() {
    cTrack.resize(650);
}

function saveImage() {
    cTrack.savePlot(4.0, "islandviewer.png", "tracks.css", 'png');
}

// Demo of the hover over timer, we had to
// do it this way to get around IE <9 not supporting
// parameters to the function called by setTimeout()
//
// If you have over an island, the console log will
// display the callback parameters when the timer expires
//
// The callback for hover (along with click) are defined in
// the data definition for each track in the dataset (data.js)
var timer;
var d_callback;
function islandPopup(d) {
    d_callback = d;
    timer = setTimeout(function() {

        console.log('clicko')

        console.log(d_callback);}, 1000);
}

function islandPopupClear(d) {
    clearTimeout(timer);
}

// Callback defined at the top of this file, for
// double clicks on the plot
function doubleClick(plotid, bp) {
    // If we have an attached linear plot, we're going
    // to refocus the zoomed area, otherwise we'll just
    // alert the user that a double click happened
    if('undefined' !== typeof linearTrack) {
        var halfBP = (cTrack.brushEndBP - cTrack.brushStartBP) /2;

	var newStart = Math.max(0, (bp - halfBP));
	var newEnd = Math.min(genomesize, (bp + halfBP));

        console.log("Moving zoom area to: " + newStart + ", " + newEnd);
        cTrack.moveBrushbyBP(newStart,
                             newEnd);
        linearTrack.update(newStart, newEnd);
    } else {
      alert("Double click! From " + plotid + " at " + bp + " bp" )
      console.log("double click!");
      console.log(plotid);
      console.log(bp);

    }
}

var radius = 40;

window.states = [
    { x : 43, y : 67, label : "first" },
    { x : 340, y : 150, label : "second" },
    { x : 200, y : 250, label : "third" },
    { x : 300, y : 320, label : "fourth" },
    { x : 50, y : 250, label : "fifth" },
    { x : 90, y : 170, label : "last" }
]

window.svg = d3.select("circularchart_svg")
.append("svg")
//.attr("viewBox", "0 0 " + 1000 + " " + 1000 )
//.attr("preserveAspectRatio", "xMinYMin")
.attr("width", "960px")
.attr("height", "500px")
.attr("background", "red");

var gStates = svg.selectAll( "g.state").data( states);

var gState = gStates.enter().append( "g")
    .attr({
        "transform" : function( d) {
            return "translate("+ [d.x,d.y] + ")";
        },
        'class'     : 'state'
    })
;

var drag = d3.behavior.drag()
.on("drag", function( d, i) {
    var selection = d3.selectAll( '.selected');

    if( selection[0].indexOf( this)==-1) {
        selection.classed( "selected", false);
        selection = d3.select( this);
        selection.classed( "selected", true);
    }

    selection.attr("transform", function( d, i) {
        d.x += d3.event.dx;
        d.y += d3.event.dy;
        return "translate(" + [ d.x,d.y ] + ")"
    })
        // reappend dragged element as last
        // so that its stays on top
    this.parentNode.appendChild( this);
    d3.event.sourceEvent.stopPropagation();
});
gState.call( drag);

gState.append( "circle")
    .attr({
        r       : radius + 4,
        class   : 'outer'
    })
;
gState.append( "circle")
    .attr({
        r       : radius,
        class   : 'inner'
    })
    .on( "click", function( d, i) {
        var e = d3.event,
            g = this.parentNode,
            isSelected = d3.select( g).classed( "selected");

        if( !e.ctrlKey) {
            d3.selectAll( 'g.selected').classed( "selected", false);
        }

        d3.select( g).classed( "selected", !isSelected);

            // reappend dragged element as last
            // so that its stays on top
        g.parentNode.appendChild( g);
    })
    .on("mouseover", function(){
        d3.select(this).style( "fill", "aliceblue");
    })
    .on("mouseout", function() {
        d3.select(this).style("fill", "white");
    });
;

gState.append( "text")
    .attr({
        'text-anchor'   : 'middle',
        y               : 4
    })
    .text( function( d) {
        return d.label;
    })
;

gState.append( "title")
    .text( function( d) {
        return d.label;
    })
;

svg
.on( "mousedown", function() {
    if( !d3.event.ctrlKey) {
        d3.selectAll( 'g.selected').classed( "selected", false);
    }

    var p = d3.mouse( this);

    svg.append( "rect")
    .attr({
        rx      : 6,
        ry      : 6,
        class   : "selection",
        x       : p[0],
        y       : p[1],
        width   : 0,
        height  : 0
    })
})
.on( "mousemove", function() {
    var s = svg.select( "rect.selection");

    if( !s.empty()) {
        var p = d3.mouse( this),
            d = {
                x       : parseInt( s.attr( "x"), 10),
                y       : parseInt( s.attr( "y"), 10),
                width   : parseInt( s.attr( "width"), 10),
                height  : parseInt( s.attr( "height"), 10)
            },
            move = {
                x : p[0] - d.x,
                y : p[1] - d.y
            }
        ;

        if( move.x < 1 || (move.x*2<d.width)) {
            d.x = p[0];
            d.width -= move.x;
        } else {
            d.width = move.x;
        }

        if( move.y < 1 || (move.y*2<d.height)) {
            d.y = p[1];
            d.height -= move.y;
        } else {
            d.height = move.y;
        }

        s.attr( d);

            // deselect all temporary selected state objects
        d3.selectAll( 'g.state.selection.selected').classed( "selected", false);

        d3.selectAll( 'g.state >circle.inner').each( function( state_data, i) {
            if(
                !d3.select( this).classed( "selected") &&
                    // inner circle inside selection frame
                state_data.x-radius>=d.x && state_data.x+radius<=d.x+d.width &&
                state_data.y-radius>=d.y && state_data.y+radius<=d.y+d.height
            ) {

                d3.select( this.parentNode)
                .classed( "selection", true)
                .classed( "selected", true);
            }
        });
    }
})
.on( "mouseup", function() {
       // remove selection frame
    svg.selectAll( "rect.selection").remove();

        // remove temporary selection marker class
    d3.selectAll( 'g.state.selection').classed( "selection", false);
})
.on( "mouseout", function() {
    if( d3.event.relatedTarget.tagName=='HTML') {
            // remove selection frame
        svg.selectAll( "rect.selection").remove();

            // remove temporary selection marker class
        d3.selectAll( 'g.state.selection').classed( "selection", false);
    }
});