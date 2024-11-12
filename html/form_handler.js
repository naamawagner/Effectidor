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
    var job_title = document.getElementById("job_title").value;
    var ORFs = document.getElementById("ORFs").value;
    var effectors = document.getElementById("effectors").value;
    var T3Es = document.getElementById("T3Es").value;
    
    if (gff != '' && (!/^[\x00-\x7F]*$/.test(gff))) {
        alert("GFF3 file name contains non ASCII characters!")
        return false;
    }
    
    if (genome != '' && (!/^[\x00-\x7F]*$/.test(genome))) {
        alert("Full genome file name contains non ASCII characters!")
        return false;
    }
    
    if (host != '' && (!/^[\x00-\x7F]*$/.test(host))) {
        alert("Host proteome file name contains non ASCII characters!")
        return false;
    }
    
    if (no_T3SS != '' && (!/^[\x00-\x7F]*$/.test(no_T3SS))) {
        alert("Proteomes of closely related bacteria without T3SS file name contains non ASCII characters!")
        return false;
    }
    
    if (ORFs != '' && (!/^[\x00-\x7F]*$/.test(ORFs))) {
        alert("ORFs file name contains non ASCII characters!")
        return false;
    }
    
    if (effectors != '' && (!/^[\x00-\x7F]*$/.test(effectors))) {
        alert("effectors file name contains non ASCII characters!")
        return false;
    }
    
    if (T3Es != '' && (!/^[\x00-\x7F]*$/.test(T3Es))) {
        alert("Effectors for homology search file name contains non ASCII characters!")
        return false;
    }
    
    if (job_title != '' && (!/^[\x00-\x7F]*$/.test(job_title))) {
        alert("Job title contains non ASCII characters!")
        return false;
    }
    
    if (host != '' && (!host.endsWith('.zip'))) {
        alert("Host proteome must be given as a zip archive!");
        return false;
    }
    if (no_T3SS != '' && (!no_T3SS.endsWith('.zip'))) {
        alert("Proteomes of close bacteria without T3SS must be given as a zip archive!");
        return false;
    }
    
    if (gff != '' && (!gff.endsWith(".zip")) && (!gff.endsWith(".txt")) && (!gff.endsWith(".gff3")) && (!gff.endsWith(".gff"))) {
        alert("GFF3 file must be gff file or a zip archive! Only gff/txt/zip file types are acceptable.");
        return false;
    }
    if (genome != '' && (!genome.endsWith(".zip")) && (!genome.endsWith(".txt")) && (!genome.endsWith(".fasta")) && (!genome.endsWith(".fna")) && (!genome.endsWith(".fsa")) && (!genome.endsWith(".fa"))) {
        alert("Full genome file must be a fasta file or a zip archive! Only fasta/fna/txt/zip file types are acceptable.");
        return false;
    }
    if (gff == '' && genome != '') {
        alert("GFF3 file is missing!");
        return false;
    }
    if (ORFs.endsWith(".zip")) {
        if (gff != '' && (!gff.endsWith(".zip"))) {
            alert("For pan-genome analysis gff input must be given in a zip archive!");
            return false;
        } else if (genome != '' && (!genome.endsWith(".zip"))) {
            alert("For pan-genome analysis full genome input must be given in a zip archive!");
            return false;
        }
    } else {
        if (gff != '' && (gff.endsWith(".zip"))) {
            alert("For a single genome analysis gff input cannot be given in a zip archive!");
            return false;
        } else if (genome != '' && (genome.endsWith(".zip"))) {
            alert("For a single genome analysis full genome input cannot be given in a zip archive!");
            return false;
        }
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
    var gff = document.getElementById("gff").value;
    var promoters = document.getElementById("promoters");
    promoters.style.display = gff != '' ? "block" : "none";
}

function openOption() {
    var x = document.getElementById("homology_effectors");
    if (x.style.display === "none") {
        x.style.display = "block";
    } else {
        x.style.display = "none";
    }
    
}