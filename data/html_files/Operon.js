function Toggle_all() {
    var checkbox = document.getElementById("show_all");
    var domaintypes = ['kripp','peptidase','transporter','regulator','biosyn','other'];
    var descriptions = document.getElementById("descriptions");
    var filler = document.getElementById("filler")
    var domain_div = document.getElementById("domains")
    var j; var i;
    Toggle_class("gene_sub_domain");
    if (checkbox.checked === true) {
//        Show all the domains 
        domain_div.style.visiblity = 'visible'
        var genes = document.getElementsByTagName("polygon");
        domains_found = ShowDomains(domaintypes,genes);

        filler.innerHTML = "Showing all operon domains."
    } else {
//        Clear all domains and show the one of the gene selected if any
        for (j = 0; j < domaintypes.length; j ++) {
            var outp_box = document.getElementById(domaintypes[j]);
            outp_box.innerHTML = '';
            outp_box.parentElement.style.display = "none";
        }
        var broken = false;
        var all_outer = document.getElementsByClassName("outer");
        for (i = 0; i < all_outer.length; i++) {
            if (all_outer[i].style.visibility === 'visible') {
                all_outer[i].style.visibility = 'hidden';
                var id = all_outer[i].attributes.getNamedItem("id").value;
                var id_arrow = id.split('.')[0];
                Toggle(id_arrow);
                broken = true;
                break;
            }
        } if (broken === false) {
            filler.innerHTML = "Select a gene to show domains.";
        }
    }
}

function Toggle_class(arg) {
    console.log("Toggle sub");
    var subs = document.getElementsByClassName(arg);
    console.log("subs: " + subs);
    console.log("length: " + subs.length);
    for (i = 0; i < subs.length; i++) {
        console.log(subs[i]);
        console.log(subs[i].style.visibility)
        if (subs[i].style.display === "none") {
            subs[i].style.display = "block";
        }
        else {
           subs[i].style.display = "none";
        }
    }

}


function Toggle(arg) {
    console.log("Toggle function");
    var outer = document.getElementById(arg + ".2");
    var inp = document.getElementById(arg);
    var filler = document.getElementById("filler");
    var domaintypes = ['precursor','kripp','peptidase','transporter','regulator','biosyn','other'];
    var geneinfo = ['gene','locus_tag','protein_id','name','cog','coordinates','pseudo','svm_hit','transl','dna','rre','rre_loc','rre_prob','rre_hit'];
    var descriptions = document.getElementById("descriptions");
    var checkbox = document.getElementById("show_all");
    
//    Clear the previous boxes
    var all_outer = document.getElementsByClassName("outer");
    var i;
    for (i = 0; i < all_outer.length; i++) {
        if (all_outer[i] != outer) {
            all_outer[i].style.visibility = "hidden";
            }
        }
//    Clear all domain lists
    var j;
    if (checkbox.checked === false) {
        for (j = 0; j < domaintypes.length; j++) {
            var domaintype = domaintypes[j];
            var outp_box = document.getElementById(domaintype);
            outp_box.innerHTML = '';
            outp_box.parentElement.style.display = "none";
        }
    }
    
//    Toggle the box around the clicked gene and show/hide domains (if any)
    if (outer.style.visibility === 'hidden') {
        outer.style.visibility = 'visible';
        ShowGeneInfo(inp,geneinfo,true);
        if (checkbox.checked === false) {
            domains_found = ShowDomains(domaintypes,[inp]);
            if (domains_found === false) {
                filler.innerHTML = 'No domains found.'
            } else {
                filler.innerHTML = ''
            }
        }
    
    } else {
        outer.style.visibility = 'hidden';
        ShowGeneInfo(inp,geneinfo,false);
        if (checkbox.checked === false) {
            filler.innerHTML = 'Select a gene to show domains.';
        }
    }
}

function ShowGeneInfo(inp,geneinfo,show) {
    var k;
    console.log("Running ShowGeneInfo");
    for (k = 0; k < geneinfo.length; k++) {
        var geneheader = geneinfo[k];
        console.log(geneheader);
        var outp_box = document.getElementById("gene." + geneheader);
        if (outp_box === null) {
            console.log("No place found to put data for " + geneheader)
            continue
        }
        outp_box.innerHTML = '';
        if (show === true) {
            if (inp.hasAttribute(geneheader)) {
                var value = inp.attributes.getNamedItem(geneheader).value;
                if (geneheader === 'transl' || geneheader ==='dna') {
                    var newv = ''
                    var l;
                    for (l = 1; (l-1)*50 < value.length; l++) {
                        newv += value.slice((l-1)*50,l*50) + '<br>';
                    }
                    value=newv;
                }
            } else {
                var value = 'N/A';
            }
            outp_box.innerHTML = value;
        } else {
            outp_box.innerHTML = '';
        }
    }
}

function ShowDomains(domaintypes,genes) {
//Display all domains in their appropriate slots
    var domains_found = false;
    var pfamlink = "https://pfam.xfam.org/family/"
    var tigrfamlink = "http://www.jcvi.org/cgi-bin/tigrfams/HmmReportPage.cgi?acc="
    var j;
    var k;
    for (j = 0; j < domaintypes.length; j++) {
        var domaintype = domaintypes[j];
        var outp_box = document.getElementById(domaintype);
        var text = '';
        var domains_seen = [];
        outp_box.innerHTML = '';
        outp_box.parentElement.style.display = "none";
        for (k = 0; k < genes.length; k++) {
            var inp = genes[k]
            if (inp.hasAttribute(domaintype)) {
                var domains = inp.attributes.getNamedItem(domaintype).value;
                var dom_split = domains.split(' ');
                for (i = 0; i < dom_split.length; i++) { 
                    if (domains_seen.includes(dom_split[i])) {
                        continue
                    } else {
                        domains_seen.push(dom_split[i]);
                    }
                    if (dom_split[i].slice(0,2) === 'PF') {
                        var baselink = pfamlink; 
                    } else {
                        var baselink = tigrfamlink;
                    }
                    console.log(dom_split[i]);
                    var descr = descriptions.attributes.getNamedItem(dom_split[i]).value;
                    text += "<a href='" + baselink + dom_split[i] + "'target='_blank'>" + dom_split[i] + "</a>" + '  ' + descr + "<br>";
                }
                if (dom_split[0] !== "") {
                     outp_box.innerHTML = text
                     outp_box.parentElement.style.display="block"
                     domains_found = true
                }
            }
        }
    }
    return domains_found;
}
//Set the collapsible elements
var coll = document.getElementsByClassName("collapsible");
var i;

for (i = 0; i < coll.length; i++) {
  coll[i].addEventListener("click", function() {
    this.classList.toggle("active");
    var content = this.nextElementSibling;
    if (content.style.display === "block") {
      content.style.display = "none";
    } else {
      content.style.display = "block";
    }
  });
};
//Set the precursors sequences
var precursors_box = document.getElementById("precursors.sequence");
var genes = document.getElementsByTagName("polygon");
var text = '';
var i;

for (i = 0; i < genes.length; i++) {
    if (genes[i].hasAttribute('svm_hit')) {
        if (genes[i].attributes.getNamedItem('svm_hit').value === "True" ) {
            var name = genes[i].attributes.getNamedItem('name').value;
            var transl = genes[i].attributes.getNamedItem('transl').value;
            transl = transl.replace(/T/g, '<tspan class="hydroxyl">T</tspan>');
            transl = transl.replace(/S/g, '<tspan class="hydroxyl">S</tspan>');
            transl = transl.replace(/C/g, '<tspan class="cysteine">C</tspan>');
            text += "<b>" + name + "</b><br>";
            text += transl + "<br>";
        }
    }
}
precursors_box.innerHTML=text;

//Uncheck the show all domains box 
document.getElementById("show_all").checked=false;




