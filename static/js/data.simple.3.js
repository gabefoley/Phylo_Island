var tracks = [
	      { trackName: "track1",
		trackType: "stranded",
		visible: true,
		inner_radius: 130,
		outer_radius: 185,
		trackFeatures: "complex",
		featureThreshold: 7000000,
		mouseclick: 'islandPopup',
		mouseover_callback: 'islandPopup',
		mouseout_callback: 'islandPopupClear',
		linear_mouseclick: 'linearPopup',
		showLabels: true,
		showTooltip: true,
		linear_mouseclick: 'linearClick',
		items: [
                         {id: 1, start:0, end:30000, name:"island0", strand: -1},
                         {id: 2, start:60000,end:100000, name:"island1", strand: -1},
                         {id: 3, start:800000,end:1000000, name:"island with a very long name", strand: 1},
                         {id: 7, start:1000000,end:1200000, name:"intergenic", strand: 0},
                         {id: 4, start:1200000,end:1500000, name:"island3", strand: 1},
                         {id: 5, start:1500000,end:1700000, name:"island4",strand: -1, extraclass: "cellwall"},
                         {id: 6, start:2000000,end:2100000, name:"the tasty gene", strand: -1, extraclass: "innermembrane"},
                         {id: 5, start:4500000,end:4700000, name:"an interesting region",strand: -1, extraclass: "cellwall"},
                         {id: 6, start:5000000,end:5100000, name:"causes bad things", strand: 1, extraclass: "innermembrane"},
                         {id: 7, start:5100000,end:5300000, name:"cures bad hair days", strand: -1, extraclass: "extracelluar"},
                         {id: 8, start:100000,end:115000, name:"cheeky gene", strand: 1, extraclass: "outermembrane"},
                         {id: 9, start:3000000,end:3100000, name:"meow meow", strand: 1, extraclass: "innermembrane"},
                         {id: 10, start:120000,end:136000, name:"turtle power", strand: -1, extraclass: "extrashelluar"},
                         {id: 11, start:3200000,end:3300000, name:"out of ideas", strand: -1, extraclass: "lostinspace"},

			]
	      }

      ];


