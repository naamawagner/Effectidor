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
    
    if (gff == '' && genome != '') {
        alert("GFF3 file is missing!");
        return false;
    } else if (gff != '' && genome == '') {
        alert("full genome FASTA file is missing!");
        return false;
    }
        
    return true;
}
