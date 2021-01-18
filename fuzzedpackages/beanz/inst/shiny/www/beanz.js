// $(document).delegate("#uiMdlInst",'DOMSubtreeModified','DOMNodeInserted',
// function(event) {
//     MathJax.Hub.Queue(["Typeset",MathJax.Hub,"uiMdlInst"]);
// });


function toggleDisplay(event, divid) {
    document.getElementById(divid).style.visibility = "visible";

    if(document.getElementById(divid).style.display == "none" ) {
	document.getElementById(divid).style.display = "";
    }
    else {
	document.getElementById(divid).style.display = "none";
    }

    event.preventDefault();
}

function toggleChkbox(event, divid, divid2) {

    if(document.getElementById(divid).checked) {
	document.getElementById(divid).checked = false;
	document.getElementById(divid2).style.display = "none";
    }
    else {
	document.getElementById(divid).checked = true;
	document.getElementById(divid2).style.visibility = "visible";
	document.getElementById(divid2).style.display = "";
    }

    event.preventDefault();
}

