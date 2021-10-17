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
    var host = document.getElementById("host").value;
    var no_T3SS = document.getElementById("no_T3SS").value;
    
    if (host != '' && (!host.endsWith('.zip'))) {
        alert("Host proteome must be given as a zip archive!");
        return false;
    }
    if (no_T3SS != '' && (!no_T3SS.endsWith('.zip'))) {
        alert("Proteomes of close bacteria without T3SS must be given as a zip archive!");
        return false;
    }
    
    if (gff != '' && (!gff.endsWith(".zip")) && (!gff.endsWith(".txt")) && (!gff.endsWith(".gff3"))) {
        alert("GFF3 file must be gff3 file or a zip archive! Only gff3/txt/zip file types are acceptable.");
        return false;
    }
    if (genome != '' && (!genome.endsWith(".zip")) && (!genome.endsWith(".txt")) && (!genome.endsWith(".fasta")) && (!genome.endsWith(".fna"))) {
        alert("Full genome file must be a fasta file or a zip archive! Only fasta/fna/txt/zip file types are acceptable.");
        return false;
    }
    if (gff == '' && genome != '') {
        alert("GFF3 file is missing!");
        return false;
    } else if (gff != '' && genome == '') {
        alert("Full genome file is missing!");
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

function ShowHideDiv() {
    var full_genome_yes = document.getElementById("full_genome_yes");
    var promoters = document.getElementById("promoters");
    promoters.style.display = full_genome_yes.checked ? "block" : "none";
}

function openOption() {
    var x = document.getElementById("homology_effectors");
    if (x.style.display === "none") {
        x.style.display = "block";
    } else {
        x.style.display = "none";
    }
    
}