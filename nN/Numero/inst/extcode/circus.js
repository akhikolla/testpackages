
// Global constants.
const SVGNS = "http://www.w3.org/2000/svg";
const DETAILWIDTH = 12;
const HOVERCOLOR = "#000000a0";

// Global variables.
var Origin = {"x":0, "y":0};
var Hovers = [];
var Activity = {};
var Membership = {};
var MenuTimer = null;
var Unit = 10.0;

// Create a new hover element assembly.
function createHoverDistrict(plot, key, dx, dy) {

    // Find target element.
    var elem = document.getElementById(plot + "_paint_" + key);
    if(!elem) {
	window.alert("District not found.")
	return;
    }
    
    // Find calibration element.
    var calibr = document.getElementById(plot + "_calibration");
    if(!calibr) {
	window.alert("Calibration failed.")
	return;
    }

    // Find the drawing area subject to translation.
    var area = document.getElementById(plot + "_contents");
    if(!area) area = document.getElementById("plot_contents");
    if(!area) {
	window.alert("Cannot determine plot contents.");
	return;
    }

    // Find the SVG element within the document.
    var svg = document.getElementById(plot);
    if(!svg) svg = document.getElementById("plot");
    if(!svg) {
	window.alert("Cannot determine svg element.");
	return;
    }

    // Determine scroll offset.
    var bbox = svg.getBoundingClientRect();
    var scroll = {"x":bbox.x, "y":bbox.y};

    // Determine translation of main contents.
    var transl = {"x":area.getAttribute("tfx"),
		  "y":area.getAttribute("tfy")};
    
    // Determine target location.
    var offset = calibr.getBoundingClientRect();
    var x = (dx + offset.x - transl.x - scroll.x);
    var y = (dy + offset.y - transl.y - scroll.y);
    var unit = offset.width;

    // New element for district highlight.
    var clon = elem.cloneNode(false);
    clon.setAttribute("pointer-events", "none");
    clon.style["fill"] = "none";
    clon.style["stroke"] = HOVERCOLOR;
    clon.style["stroke-width"] = "2px";
    clon.style["stroke-linejoin"] = "round";

    // New element for info text.
    var txt = document.createElementNS(SVGNS, "text");
    txt.setAttribute("pointer-events", "none");
    txt.style["font-family"] = "'Arial'";
    txt.style["font-size"] = (Math.round(0.5*unit) + "px");
    txt.style["font-weight"] = "bold";
    txt.style["text-anchor"] = "middle";
    txt.style["fill"] = "white";

    // New element for text background.
    var bg = document.createElementNS(SVGNS, "rect");
    bg.setAttribute("rx", 0.1*unit);
    bg.setAttribute("ry", 0.1*unit);
    bg.setAttribute("pointer-events", "none");
    bg.style["fill"] = HOVERCOLOR;
    bg.style["stroke"] = "none";

    // Set text contents.
    txt.textContent = LABELS[plot][key];

    // Relative position to map center in screen coordinates.
    var dx = (elem.getAttribute("v2"))*unit;
    var dy = (elem.getAttribute("v3"))*unit;
    
    // Set text position.
    txt.setAttribute("x", x);
    txt.setAttribute("y", (y - 0.5*unit));
    
    // Check if district is labeled.
    var label = document.getElementById(plot + "_label_" + key);
    if(label) {
	if(label.getAttribute("visibility") != "hidden") {
	    clon.setAttribute("blocked", label.id);
	    label.setAttribute("visibility", "hidden");
	}
    }

    // Add to parent element.
    area.appendChild(clon);
    area.appendChild(bg);
    area.appendChild(txt);
    
    // Determine text dimensions.
    var tbox = txt.getBoundingClientRect();
    var w = (tbox.width + 0.3*unit);
    var h = (tbox.height + 0.1*unit);
    
    // Finish background box.
    bg.setAttribute("x", (x - 0.5*w));
    bg.setAttribute("y", (y - 0.67*unit - 0.5*h));
    bg.setAttribute("width", w);
    bg.setAttribute("height", h);

    // Add elements to global vector.
    Hovers.push(clon);
    Hovers.push(bg);
    Hovers.push(txt);
}

// Create elements for a symbol assembly.
function createSymbol(x, y, unit, shape, strokeFlag) {

    // Additional modifiers.
    var fontSize = 0.47*unit;
    var strokeWidth = 0.08*unit;
    var radius1 = (0.45*unit - strokeWidth);
    var radius0 = (radius1 - 2.3*strokeWidth);

    // Adjust radius if reduced stroke.
    if(!strokeFlag) radius1 -= 0.3*strokeWidth;
    
    // Halo around activated label and marker for inactivity.
    var halo = undefined;
    var marker = undefined;
    if(shape == "circle") {
	halo = document.createElementNS(SVGNS, "circle");
	halo.setAttribute("cx", x);
	halo.setAttribute("cy", y);
	halo.setAttribute("r", radius1);
	marker = document.createElementNS(SVGNS, "circle");
	marker.setAttribute("cx", x);
	marker.setAttribute("cy", y);
	marker.setAttribute("r", radius0)
    }
    if(shape == "square") {
	halo = document.createElementNS(SVGNS, "rect");
	halo.setAttribute("x", (x - 0.93*radius1));
	halo.setAttribute("y", (y - 0.93*radius1));
	halo.setAttribute("width", 1.86*radius1);
	halo.setAttribute("height", 1.86*radius1);
	marker = document.createElementNS(SVGNS, "rect");
	marker.setAttribute("x", (x - 0.92*radius0));
	marker.setAttribute("y", (y - 0.92*radius0));
	marker.setAttribute("width", 1.84*radius0);
	marker.setAttribute("height", 1.84*radius0);
    }
    if(shape == "diamond") {
	var xy1 = (String(x - 1.2*radius1) + "," + String(y));
	xy1 += (" " + String(x) + "," + String(y + 1.25*radius1));
	xy1 += (" " + String(x + 1.2*radius1) + "," + String(y));
	xy1 += (" " + String(x) + "," + String(y - 1.25*radius1));
	halo = document.createElementNS(SVGNS, "polygon");
	halo.setAttribute("points", xy1);
	var xy0 = (String(x - 1.2*radius0) + "," + String(y));
	xy0 += (" " + String(x) + "," + String(y + 1.3*radius0));
	xy0 += (" " + String(x + 1.2*radius0) + "," + String(y));
	xy0 += (" " + String(x) + "," + String(y - 1.3*radius0));
	marker = document.createElementNS(SVGNS, "polygon");
	marker.setAttribute("points", xy0);
    }
    if(!halo) {
	window.alert("createSymbol: Unknown shape.");
	return null;
    }
    
    // Set halo style and attributes.
    halo.style["stroke"] = "#ffffff";
    if(strokeFlag) halo.style["stroke-width"] = strokeWidth;
    else halo.style["stroke-width"] = 0.6*strokeWidth;
    halo.setAttribute("active", true);
    halo.setAttribute("shape", shape);
	
    // Activated region label.
    var label = document.createElementNS(SVGNS, "text");
    label.style["font-family"] = "'Arial'";
    label.style["fill"] = "#ffffff";
    label.style["font-weight"] = "bold";
    label.style["font-size"] = (Math.round(fontSize) + "px");
    label.style["text-anchor"] = "middle";
    label.setAttribute("pointer-events", "none");
    label.setAttribute("x", x);
    label.setAttribute("y", (y + 0.33*fontSize));
    label.setAttribute("yoffset", 0.33*fontSize);
    label.setAttribute("active", true);

    // Set marker style and attributes.
    marker.style["stroke"] = "#ffffff";
    marker.style["stroke-width"] = 0.6*strokeWidth;
    marker.setAttribute("active", false);
    marker.setAttribute("shape", shape);
    
    // Return results.
    return {"halo":halo, "label":label, "marker":marker};
}

// Respond to clicking on a district.
function downloadRegions() {

    // Column headings.
    var colnames = Object.keys(TOPOLOGY);
    var data = colnames[0];
    for(var j = 1; colnames.length > j; j++)
	data += ("\t" + colnames[j]);
    data += "\n";

    // Copy topology.
    var topo = {};
    for(var j = 0; colnames.length > j; j++) {
	var c = colnames[j];
	topo[c] = TOPOLOGY[c].slice(0);
    }

    // Collect membership data.
    for(var i = 0; NDISTRICTS > i; i++) {
	var rkey = Activity.history[i].current;
	var detl = document.getElementById("detail_" + rkey);
	topo["REGION"][i] = "not_selected";
	topo["REGION.label"][i] = "";
	topo["REGION.color"][i] = "";
	if(1 > rkey.length) continue;
	topo["REGION"][i] = detl.textContent;
	topo["REGION.label"][i] = rkey;
	topo["REGION.color"][i] = REGIONCOLORS[rkey];
    }
   
    // Column values.
    for(var i = 0; NDISTRICTS > i; i++) {
        data += topo[colnames[0]][i];
	for(var j = 1; colnames.length > j; j++)
	    data += ("\t" + topo[colnames[j]][i]);
	data += "\n";
    }
    
    // Prepare output object.
    var output = new Blob([data], {type: 'text/plain'});
    var url = window.URL.createObjectURL(output);

    // Trigger download via a hyperlink.
    var link = document.getElementById("download_data");
    link.setAttribute("href", url);   
    link.click();

    // Release data object.
    window.URL.revokeObjectURL(url);   
}

// Collect district-specific data.
function getDistrictData(elem) {
    if(!elem) return null;
    if((typeof elem) != "object") return null;
    
    // Extract data fields.
    var subplot = elem.getAttribute("v0"); 
    var key = elem.getAttribute("v1");
    var dx = elem.getAttribute("v2");
    var dy = elem.getAttribute("v3");
    if(!subplot) return null;
    if(!key) return null;
    if(!dx) return null;
    if(!dy) return null;

    // Determine activity status and region assignment.
    var active = elem.getAttribute("active"); 
    var region = elem.getAttribute("region"); 
    active = (active == "true");
    
    // Finish results.
    output = {};
    output.subplot = subplot;
    output.key = parseInt(key);
    output.dx = parseFloat(dx);
    output.dy = parseFloat(dy);
    output.active = active;
    output.region = region;
    return output;
}

// Create menu items for subgroup selection.
function initMenu() {

    // Check that regions have been defined.
    if(!REGIONS) {
 	window.alert("initMenu: Regions not defined.");
	return false;
    }
    
    // Additional modifiers.
    const fontSize = 0.47*Unit;
    const textWidth = DETAILWIDTH*fontSize;
    const lineWidth = 0.08*Unit;
  
    // Determine menu block size.
    const width = (textWidth + (1.0 + 0.5)*Unit);
    const height = (Object.keys(REGIONS).length)*Unit;
 
    // Adjust the hovering side block that contains the menu.
    var side = document.getElementById("side");
    side.style.width = (width + "px");
    side.style.height = (height + "px");

    // Top menu layer with editable text and a gap for
    // interacting with symbols on the layer below.
    var top = document.getElementById("top");
    top.style.marginLeft = (1.33*Unit + "px");
  
    // SVG block that contains subgroup details.
    var detl = document.getElementById("detail");
    detl.setAttribute("width", width);
    detl.setAttribute("height", height);

    // SVG block that contains symbols and the background.
    var sym = document.getElementById("symbol");
    sym.setAttribute("width", width);
    sym.setAttribute("height", height);

    // Menu background.
    var menubg = document.getElementById("symbol_background");
    menubg.setAttribute("x", 0);
    menubg.setAttribute("y", 0);
    menubg.setAttribute("width", width);
    menubg.setAttribute("height", height);
    menubg.setAttribute("rx", 0.2*Unit);
    menubg.setAttribute("ry", 0.2*Unit);
    menubg.style["fill"] = "#ffffffa0";
    menubg.style["stroke"] = "none";
 
    // Create menu items.
    var baseline = 0;
    for(var k in REGIONS) {
	var x = (0.5*Unit + lineWidth);
	var y = (baseline + 0.5)*Unit;
	var r = (0.5*Unit - lineWidth);

	// Create a group for elements in activated highlight.
        var hlight = document.createElementNS(SVGNS, "g");
	hlight.id = ("symbol_active_" + k);
	
	// Background for activated detail text area.
	var txtbg = document.createElementNS(SVGNS, "rect");
	txtbg.style["fill"] = REGIONCOLORS[k];
	txtbg.style["fill-opacity"] = 0.4;
	txtbg.setAttribute("pointer-events", "none");
	txtbg.setAttribute("x", x);
	txtbg.setAttribute("y", (y - r + 0.81*lineWidth));
	txtbg.setAttribute("width", (textWidth + 0.2*Unit));
	txtbg.setAttribute("height", 2*(r - 0.81*lineWidth));
	txtbg.setAttribute("rx", 0.2*r);
	txtbg.setAttribute("ry", 0.2*r);
	
        // Add background to highlight group.
	hlight.appendChild(txtbg);

        // Create symbol assembly.
	var shape = "circle";
	if(baseline%3 == 1) shape = "square";
	if(baseline%3 == 2) shape = "diamond";
	var elems = createSymbol(x, y, Unit, shape, true);

	// Set paint color.
	elems.halo.style["fill"] = REGIONCOLORS[k];
 	elems.marker.style["fill"] = REGIONCOLORS[k];

	// Set region identity.
	elems.halo.setAttribute("key", k);
 	elems.marker.setAttribute("key", k);
	elems.label.setAttribute("key", k);
	elems.label.textContent = k;
	
        // Background to prevent text box from
	// showing through transparent halo.
	var halobg = elems.halo.cloneNode(false);
	halobg.setAttribute("pointer-events", "none");
	halobg.style["fill"] = halobg.style["stroke"]

	// Set element identities.
	elems.halo.id = ("symbol_halo_" + k);
 	elems.label.id = ("symbol_label_" + k);
	elems.marker.id = ("symbol_inactive_" + k);
	halobg.id = ("symbol_background_" + k);
	
        // Add highlight elements to group.
	hlight.appendChild(halobg);
	hlight.appendChild(elems.halo);
	hlight.appendChild(elems.label);
	
	// Set initial activity status.
	hlight.style["visibility"] = "hidden";
	elems.marker.style["visibility"] = "visible";

        // Add activity items to the symbol block.
        sym.appendChild(hlight);
        sym.appendChild(elems.marker);
        baseline++;
    }
        
    // Add menu-specific event listeners.
    var middle = document.getElementById("middle");
    if(window.PointerEvent) {
	middle.addEventListener("pointerover", pointerOverSymbol, false); 
	middle.addEventListener("pointerdown", pointerDownSymbol, false);
    }
    else {
	middle.addEventListener("mouseover", pointerOverSymbol, false); 
	middle.addEventListener("mousedown", pointerDownSymbol, false);
    }
    
    // Launch background refresh for updates.
    if(!MenuTimer) {
	refreshMenu(); // initial set up
	MenuTimer = window.setInterval(refreshMenu, 100);
    }
    return true;
}

// Launch widget.
function initPage(plot, guiFlag) {
    
    // Check that there is at least one plot.
    if(!SUBPLOTS) {
 	window.alert("initPage: Subplots not defined.");
	return;
    }
    if(1 > SUBPLOTS.length) {
 	window.alert("initPage: No subplots.");
	return;
    }
    
    // Determine window content offset.
    var conts = document.getElementById(plot + "_contents")
    var item0 = conts.transform.baseVal.getItem(0);
    if(item0.type == SVGTransform.SVG_TRANSFORM_TRANSLATE) {
        Origin.x = item0.matrix.e;
        Origin.y = item0.matrix.f;
    }
    
    // No interactive elements.
    if(!guiFlag) {
	document.addEventListener("pointerover", pointerOverDistrict, false);
	return;
    }

    // Clear any shadow elements.
    for(var j = 0; SUBPLOTS.length > j; j++) {
	var sGroup = document.getElementById(SUBPLOTS[j] + "_shadow");
	if(sGroup) sGroup.parentNode.removeChild(sGroup);
    }
    
    // Find calibration element.
    var id = (SUBPLOTS[0] + "_calibration");
    var calibr = document.getElementById(id);
    if(!calibr) {
 	window.alert("initPage: Calibration failed.");
	return;
    }

    // Tables to keep track of activity.
    Activity.hold = false;
    Activity.recent = "";
    Activity.regions = {};
    Activity.history = [];
    for(var i = 0; NDISTRICTS > i; i++)
	Activity.history.push({"current":"", "previous":""});
    
    // Set up membership matrix for regions.
    Membership = {};
    for(var k in REGIONS)
	Membership[k] = {"length":0};

    // Set distance unit.
    var bbox = calibr.getBoundingClientRect();
    Unit = bbox.width;

    // Set up subgroup selector and define symbol shapes
    // for regions, must be before initRegions().
    if(!initMenu()) {
 	window.alert("initPage: Menu failed.");
	return;
    }
    
    // Set up regions.
    if(!initRegions(TOPOLOGY)) {
 	window.alert("initRegions: Highlights failed.");
	return;
    }
    
    // Add content-specific event listeners.
    var bottom = document.getElementById("bottom");
    if(window.PointerEvent) {
	bottom.addEventListener("pointerover", pointerOverDistrict, false);
	bottom.addEventListener("pointerdown", pointerDownDistrict, false);
    }
    else {
	bottom.addEventListener("mouseover", pointerOverDistrict, false);
	bottom.addEventListener("mousedown", pointerDownDistrict, false);
    }
}

// Set up region elements.
function initRegions(topology) {

    // Clear any existing highlights.
    for(var j = 0; SUBPLOTS.length > j; j++) {
	var sGroup = document.getElementById(SUBPLOTS[j] + "_region");
        if(sGroup) sGroup.parentNode.removeChild(sGroup);
    }
    
    // Create highlights according to topology.
    for(var i = 0; NDISTRICTS > i; i++) {
        var rkey = topology["REGION.label"][i];
	if(1 > rkey.length) continue;
	 
        // District details.
        var key = (SUBPLOTS[0] + "_paint_" + i);
	var target = document.getElementById(key)
        var dst = getDistrictData(target);
	var history = Activity.history[dst.key];

        // Update membership.
        history.current = rkey;
	Membership[history.current][i] = target;
	Membership[history.current].length += 1;
        updateHighlights(dst, history.current);
    }
    return true;
}

// Location and scale of a subplot.
function locateSubplot(sub) {

    // Find the drawing area subject to translation.
    var area = document.getElementById(sub + "_contents");
    if(!area) area = document.getElementById("plot_contents");
    if(!area) {
	window.alert("locateSubplot: Cannot find contents.");
	return null;
    }
	
    // Find the SVG element within the document.
    var svg = document.getElementById(sub);
    if(!svg) svg = document.getElementById("plot");
    if(!svg) {
	window.alert("locateSubplot: Cannot find svg element.");
	return null;
    }

    // Find calibration element.
    var calibr = document.getElementById(sub + "_calibration");
    if(!calibr) calibr = document.getElementById("plot_calibration");
    if(!calibr) {
	window.alert("locateSubplot: Calibration failed.");
	return null;
    }
	
    // Determine scroll offset.
    var bbox = svg.getBoundingClientRect();
    var scroll = {"x":bbox.x, "y":bbox.y};
	
    // Determine translation of main contents.
    var transl = {"x":parseFloat(area.getAttribute("tfx")),
		  "y":parseFloat(area.getAttribute("tfy"))};

    // Finish calibration box.
    var offset = calibr.getBoundingClientRect();
    offset.x -= (transl.x + scroll.x);
    offset.y -= (transl.y + scroll.y);

    // Return results.
    var output = {};
    output.x = offset.x;
    output.y = offset.y;
    output.unit = offset.width;
    output.area = area;
    return output;
}

// Respond to clicking on a district.
function pointerDownDistrict(evt) {
    
    // Clear any hover elements.
    removeHovers();

    // Check if over a district.
    var dst = getDistrictData(evt.target);
    if(!dst) {
	if(Activity.hold) return;
	if(window.confirm("Download results?"))
	    downloadRegions();
	return;
    }

    // Check that there is a single active region.
    var nregs = 0;
    var rkey = "";
    for(var k in Activity.regions) {
	if(Activity.regions[k]) {
	    rkey = k;
 	    nregs++;
	}
	if(nregs > 1) {
	    window.alert("Please select a single region.");
	    return;
	}
    }
    if(1 > nregs) {
	window.alert("Please select a region first.");
	return;
    }
    
    // Current and previous assignments.
    var history = Activity.history[dst.key];

    // Determine replacement region.
    if(Activity.hold) {

        // Check if clearance is needed.
	if(1 > Activity.recent.length) {
	    if(history.current != rkey) return;
	    rkey = history.previous;
	}
	else {
	    rkey = Activity.recent;
	}
	
	// Check if already up-to-date.
	if(history.current == rkey) return;
    }

    // Remove current assignment.
    if(history.current.length > 0) {
	Membership[history.current][dst.key] = null;
	Membership[history.current].length -= 1;
    }
	
    // If repeated click, restore previous assignment.
    // If new click, update assignment.
    if(history.current == rkey) {
	history.current = history.previous;
	history.previous = "";
	if(!Activity.hold) Activity.recent = "";
    }
    else {
	history.previous = history.current;
	history.current = rkey;
        if(!Activity.hold) Activity.recent = rkey;
    }
    
    // Set new assignment.
    if(history.current.length > 0) {
	Membership[history.current][dst.key] = evt.target;
	Membership[history.current].length += 1;
    }

    // Update visuals.
    updateHighlights(dst, history.current);
 }

// Respond to clicking on a symbol area.
function pointerDownSymbol(evt) {
    
    // Check if valid target.
    var rkey = evt.target.getAttribute("key");
    if(!rkey) return;

    // Check activity status.
    var active = evt.target.getAttribute("active");
    if(!active) {
	window.alert("pointerDownSymbol: Unusable activity status.");
	return;
    }

    // Check that region is defined in the symbol menu.
    var target0 = document.getElementById("symbol_inactive_" + rkey);
    var target1 = document.getElementById("symbol_active_" + rkey);
    if(!target0 || !target1) {
	window.alert("pointerDownSymbol: Region not defined.");
	return;
    }
    
    // Activate a region that is currently inactive.
    if(active == "false") {
	target0.style["visibility"] = "hidden";
	target1.style["visibility"] = "visible";
	Activity.regions[rkey] = true;
    }

    // De-activate a region that is currently active.
    if(active == "true") {
	target0.style["visibility"] = "visible";
	target1.style["visibility"] = "hidden";
	Activity.regions[rkey] = false;
    }
    
    // Update districts.
    var memb = Membership[rkey];
    for(k in memb) {
	var dst = getDistrictData(memb[k]);
        updateHighlights(dst, rkey);
    }
}

// Respond to entering an element's area.
function pointerOverDistrict(evt) {
    
    // Check if pointer is down and over a district.
    Activity.hold = (evt.buttons != 0);
    if(Activity.hold) {
	pointerDownDistrict(evt);
	return;
    }

    // Remove existing hover elements.
    removeHovers();

    // Check if over a district.
    var dst = getDistrictData(evt.target);
    if(!dst) return;
    
    // Create new hover elements.
    for(var j = 0; SUBPLOTS.length > j; j++)
	createHoverDistrict(SUBPLOTS[j], dst.key, dst.dx, dst.dy);
}

// Respond to pointer within a menu symbol.
function pointerOverSymbol(evt) {
    
    // Remove any previous hover elements.
    removeHovers();

    // Check if pointer is down.
    if(evt.buttons != 0) return;
    
    // Find the active version of the symbol.
    var key = evt.target.getAttribute("key");
    var elem = document.getElementById("symbol_halo_" + key);
    if(!elem) return;
   
    // Create a new hover element.
    var halo = elem.cloneNode(false);
    var lw = parseFloat(halo.style["stroke-width"])
    halo.style["stroke-width"] = lw;
    halo.style["stroke"] = HOVERCOLOR;
    halo.style["fill"] = "none";
    halo.setAttribute("pointer-events", "none");
    
    // Add element to the figure.
    evt.target.parentElement.appendChild(halo);
    Hovers.push(halo);
}

// Document-wide subgroup activation.
function refreshMenu() {
    const fontSize = 0.55*Unit;

    // Check that menu is available.
    var detl = document.getElementById("detail");
    if(!detl) {
	window.alert("refreshMenu: Details not defined.");
	return false;
    }
    
    // Make sure text is up-to-date.
    var reserved = {};
    for(k in REGIONS) {
	var txt = document.getElementById("detail_" + k);

	// Create a text element for region detail.
	if(!txt) {
	    var label = document.getElementById("symbol_label_" + k);
	    var halo = document.getElementById("symbol_halo_" + k);
            var halobox = halo.getBoundingClientRect();
	    var y = (0.5*(halobox.top + halobox.bottom) - 0.25*fontSize);
	    txt = document.createElementNS(SVGNS, "text");
	    txt.style["fill"] = "#a0a0a0";
            txt.style["font-family"] = label.style["font-family"];
            txt.style["font-size"] = (Math.round(fontSize) + "px");
            txt.setAttribute("x", 0);
            txt.setAttribute("y", y);
            txt.setAttribute("key", k);
            txt.id = ("detail_" + k);
	    detl.appendChild(txt);
	}

        // Check if first refresh.
	if(!MenuTimer) txt.textContent = REGIONS[k];
	
	// Make sure name is usable.
	while(2 > txt.textContent.length) {
	    var msg = ("Please enter a name for Region " + k + ":");
	    removeHovers(); // hovers may get stuck after prompt
	    txt.textContent = window.prompt(msg, "");
	}
	
	// Include region name in the lookup table.
	reserved[txt.textContent] = k;
    }

    // Update text style according to memberships and activity.
    for(k in REGIONS) {
	var txt = document.getElementById("detail_" + k);
	var memb = Membership[k];
	if(Activity.regions[k] == true) {
	    if(memb.length > 0) {
		txt.style["font-weight"] = "bold";
		txt.style["fill"] = "#000000a0";
	    }
	    else {
		txt.style["font-weight"] = "normal";
		txt.style["fill"] = "#00000080";
	    }
	}
	else {
	    if(memb.length > 0) {
		txt.style["font-weight"] = "bold";
		txt.style["fill"] = "#505050";
	    }
	    else {
		txt.style["font-weight"] = "normal";
		txt.style["fill"] = "#a0a0a0";
	    }
	}

	// Check if the name is duplicated.
	if(reserved[txt.textContent] != k)
	    txt.style["fill"] = "#fa1010";
	reserved[txt.textContent] = k;
    }
}

// Clear all hover elements.
function removeHovers() {
    for(var i = 0; Hovers.length > i; i++) {
	var parent = Hovers[i].parentElement;
	var blocked = Hovers[i].getAttribute("blocked");
        parent.removeChild(Hovers[i]);
	if(!blocked) continue;
	var label = document.getElementById(blocked);
	label.setAttribute("visibility", "visible");
    }
    Hovers = [];
}

// Update region indicators on a district.
function updateHighlights(dst, rkey) {
    if(!dst) return;
    for(var j = 0; SUBPLOTS.length > j; j++) {
	var sub = SUBPLOTS[j];

	// Check if valid district.
	var elem = document.getElementById(sub + "_paint_" + dst.key);
	if(!elem) {
	    window.alert("updateDistrict: Unusable element.");
	    return;
	}

	// Find highlight elements.
        const hlkey1 = (sub + "_active_" + dst.key);
	const hlkey0 = (sub + "_inactive_" + dst.key);
	var hlight1 = document.getElementById(hlkey1);
	var hlight0 = document.getElementById(hlkey0);

        // Determine district location.
        var offset = locateSubplot(sub);
	var x = (offset.x + dst.dx);
	var y = (offset.y + dst.dy);
	
        // Remove previous highlights if region changed.
	var prev = getDistrictData(elem);
	if(prev.region != rkey) {
	    if(hlight1) offset.area.removeChild(hlight1);
	    if(hlight0) offset.area.removeChild(hlight0);
	    hlight1 = null;
	    hlight2 = null;
	}

	// Update district data.
	elem.setAttribute("region", rkey);
	
	// Check if district is labeled.
	var label = document.getElementById(sub + "_label_" + dst.key);
	if(1 > rkey.length) {
	    if(label) label.setAttribute("visibility", "visible");
	    continue;
	}

	// Hide label to prevent overlap with highlight.
	if(label) label.setAttribute("visibility", "hidden");	

	// Create new highlight assemblies.
	if(!hlight1 || !hlight0) {

            // Determine highlight type and color.
	    var outline = document.getElementById("symbol_halo_" + rkey);
            if(!outline) window.alert("updateHighlights: Missing symbol.");
	    var shape = outline.getAttribute("shape");
	
            // Create a new symbol assembly.
	    var r = 0.87*(offset.unit);
	    var symbol = createSymbol(x, y, r, shape, false);
	
	    // Set active halo style and attributes.
 	    symbol.halo.style["fill"] = outline.style["fill"];
	    symbol.halo.setAttribute("pointer-events", "none");

	    // Set label style and attributes.
	    symbol.label.textContent = rkey;
 	    symbol.label.setAttribute("pointer-events", "none");
	
	    // Create a group for active highlight elements.
            hlight1 = document.createElementNS(SVGNS, "g");
	    hlight1.appendChild(symbol.halo);
	    hlight1.appendChild(symbol.label);
	    hlight1.id = hlkey1;

	    // Set inactive marker style and attributes.
	    hlight0 = symbol.marker;
	    hlight0.style["fill"] = outline.style["fill"];
	    hlight0.setAttribute("pointer-events", "none");
	    hlight0.id = hlkey0;
	    
	    // Add highlight on the district.
	    offset.area.appendChild(hlight0);
	    offset.area.appendChild(hlight1);
	}
	    
	// Set visibility based on activity status.
	if(Activity.regions[rkey]) {
	    hlight1.style["visibility"] = "visible";
	    hlight0.style["visibility"] = "hidden";
	} else {
	    hlight1.style["visibility"] = "hidden";
	    hlight0.style["visibility"] = "visible";
	}
    }
}
