function openAdvanced() {
    var x = document.getElementById("advanced");
    if (x.style.display === "none") {
        x.style.display = "block";
    } else {
        x.style.display = "none";
    }
}


function validate_read_files() {
    var gff = document.getElementById("gff").value;
    var genome = document.getElementById("genome").value;
    
    if (gff != '' && (!gff.endsWith(".zip"))) {
        alert("GFF3 file must be in a zip archive!");
        return false;
    }
    if (genome != '' && (!genome.endsWith(".zip"))) {
        alert("full genome file must be in a zip archive!");
        return false;
    }
    if (gff == '' && genome != '') {
        alert("GFF3 file is missing!");
        return false;
    } else if (gff != '' && genome == '') {
        alert("full genome file is missing!");
        return false;
    }
        
    return true;
}

function openCIS() {
    var gff = document.getElementById("gff").value;
    var genome = document.getElementById("genome").value;
    
    if (gff != '' && genome != '') {
        var x = document.getElementById("CIS");
        if (x.style.display === "none") {
            x.style.display = "block";
        } else {
            x.style.display = "none";
        }
    } else {
        var x = document.getElementById("CIS");
        document.getElementById("PIP-box").checked = false;
        document.getElementById("hrp-box").checked = false;
        document.getElementById("mxiE-box").checked = false;
        document.getElementById("exs-box").checked = false;
        x.style.display = "none";
    }
}
