// Script to parse requirement and filter the given operons
// 1) Function to convert passed requirements to js dictionary
// 2) Set the dictionary to localStorage (for other webpages)
// 3) Figure out which operons are filtered, store them per group
// 4) Figure out which groups are filtered and hide them

function standardDeviation(values){
  var avg = average(values);
  
  var squareDiffs = values.map(function(value){
    var diff = value - avg;
    var sqrDiff = diff * diff;
    return sqrDiff;
  });
  
  var avgSquareDiff = average(squareDiffs);

  var stdDev = Math.sqrt(avgSquareDiff);
  return [avg, stdDev];
}

function average(data){
  var sum = data.reduce(function(sum, value){
    return sum + value;
  }, 0);

  var avg = sum / data.length;
  return avg;
}

function ParseRequirements_test() {
    var requirements = { 
        "length_core": {
            "type": "equal",
            "value": 10
            }
        }
    return requirements
}

function ParseRequirements(requirements) {
    length_req = ParseRequirement("length_core");
    cog_req = ParseRequirement("COG_avg");
    domain_req = ParseRequirementDomainType();
    domain_spec_req = ParseRequirementSpecificDomain();
    reqs = [length_req,cog_req,domain_req,domain_spec_req];
    for (i = 0; i < reqs.length; i ++ ) {
        req = reqs[i];
        console.log(req);
        if (! Array.isArray(req[2]) && req[2] !== "" && ! isNaN(req[2]))  {
            //Overwrite the current values (not for arrays)
            requirements[req[0]] = {
                "type": req[1],
                "value": req[2]
            };
        } else if (Array.isArray(req[2]) && req[2].length > 0) {
            //Update the Array currently used if there is already one
            key = req[0]
            if (key in requirements) {
                new_value = requirements[key]["value"].concat(req[2]);
                // Get only the unique elements
                temp = {};
                for (j = 0; j < new_value.length; j++) {
                    temp[new_value[j]] = true
                };
                new_value = Object.keys(temp);
            } else {
                new_value = req[2];
            }
            requirements[key] = {
                "type": req[1],
                "value": new_value
            };
        } else {
            delete requirements[req[0]];
        };
    };
    console.log(requirements);
};

function ParseRequirement(name) {
    var name_type = "type_" + name;
    var name_value = "value_" + name;
    var type = document.getElementById(name_type).value;
    var value = parseFloat(document.getElementById(name_value).value);
    return([name,type,value]);
};

function ParseRequirementDomainType() {
    var domaintype = document.getElementById("domain_select_domaintype").value;
    var domainlocation = document.getElementById("domain_select_location").value;
    var filtertype = document.getElementById("domain_select_filtertype").value;
    var domainnr = parseInt(document.getElementById("domain_select_number").value);
    if (domainlocation === "all") {
        domaintype += '_all';
    };
    return([domaintype,filtertype,domainnr])
}

function ParseRequirementSpecificDomain() {
    var text = document.getElementById("specific_domains").value;
    if (text.includes(" ") || text.includes(",") ) {
        var domains = text.split(/[\s,]+/);
    } else {
        var domains = [text]
    }
    domains_clean = [];
    for (i=0; i<domains.length; i++) {
        domain = domains[i].replace(" ","");
        if (domain !== "") {
            domains_clean.push(domain);
        }
    }
    if (domains_clean.length === 1 && domains_clean[0] === "") {
        domains_clean = "";
    }
    return(["all_domains","contains",domains_clean]);
}

function FilterOperons(operons,requirements) {
    console.log("FilterOperons")
    var filtered_operons = [];
    for (var operon in operons) {
        console.log("FilterOperons " + operon);
        var attributes = operons[operon];
        var all_passed = true;
        for (var attr in requirements) {
            var type_req = requirements[attr]["type"];
            var value_req = requirements[attr]["value"];
            if (attr in attributes) {
                value_operon = attributes[attr];
            } else {
                if (type_req === "contains") {
                    value_operon = [];
                } else {
                    value_operon = 0;
                };
            };
            // Comparing ints to arrays --> get lengths
            if (Array.isArray(value_operon) && type_req !== "contains") {
                console.log("using length");
                value_operon_real = value_operon.length;
            } else {
                value_operon_real = value_operon;
            }
            if (type_req === "equal") {
                if (value_operon_real !== value_req) {
                    all_passed = false;
                }
            } else if (type_req === "min") {
                if (value_operon_real < value_req) {
                    all_passed = false;
                }
            } else if (type_req === "max") {
                if (value_operon_real > value_req) {
                    all_passed = false;
                }
            } else if (type_req === "contains") {
                var found_all = true ;
                for (i = 0; i < value_req.length; i ++) {
                    item_req = value_req[i]
                    if (value_operon_real.includes(item_req) === false ) {
                        found_all = false;
                        break;
                    }
                }
                if (found_all === false) {
                    all_passed = false;
                }
            }
        }
        if (all_passed === true) {
            console.log("operon passed");
            filtered_operons.push(operon);
        }
    }
    return filtered_operons
}

function ChangeTableDataEntry(nr) {
    nr_operons_display = document.getElementById("nr_operons")
    nr_operons_display.innerText = nr
}


function FilterOperonsEntry(operons_dict,operons,requirements) {
    var filtered_operons = [];
    for (i=0 ; i<operons.length ;i++) {
        operon = operons[i];
        all_passed = true;
        var attributes = operons_dict[operon];
        for (var attr in requirements) {
            var type_req = requirements[attr]["type"];
            var value_req = requirements[attr]["value"];
            if (attr in attributes) {
                value_operon = attributes[attr];
                // Comparing ints to arrays --> get lengths
                if (Array.isArray(value_operon) && type_req !== "contains") {
                    console.log("using length");
                    value_operon_real = value_operon.length;
                } else {
                    value_operon_real = value_operon;
                }
                if (type_req === "equal") {
                    if (value_operon_real !== value_req) {
                        all_passed = false;
                    }
                } else if (type_req === "min") {
                    if (value_operon_real < value_req) {
                        all_passed = false;
                    }
                } else if (type_req === "max") {
                    if (value_operon_real > value_req) {
                        all_passed = false;
                    }
                } else if (type_req === "contains") {
                    var found_all = true ;
                    for (j = 0; j < value_req.length; j ++) {
                        item_req = value_req[j]
                        if (value_operon_real.includes(item_req) === false ) {
                            found_all = false;
                            break;
                        }
                    }
                    if (found_all === false) {
                        all_passed = false;
                    }
                }
            }
        }
        if (all_passed === true) {
            console.log("operon passed");
            filtered_operons.push(operon);
        }
    }
    return filtered_operons
}

function FilterGroups(groups,operons) {
    var filtered_groups = {};
    for (groupname in groups) {
        operons_group = groups[groupname];
        operons_filtered = [];
        for (i = 0; i < operons_group.length; i ++) {
            operon = operons_group[i];
            if (operons.includes(operon)) {
                operons_filtered.push(operon);
            }
        }
        if (operons_filtered.length > 0) {
            filtered_groups[groupname] = operons_filtered;
        }
    }
    return filtered_groups
}

function ChangeTableData(groups,operons) {
    console.log("ChangeTableData with ");
    console.log(operons);
    for (groupname in groups) {
        operons_group = groups[groupname];
        console.log(groupname);
        all_COGs = [];
        all_scores = [];
        total_mibig = 0;
        total_antismash = 0;
        total_filtered = 0;
        total_kripp = 0;
        for (i=0; i<operons_group.length; i++) {
            operon = operons_group[i];
            console.log(operon)
            data = operons[operon];
            for (j = 0; j < data['COG'].length; j++) {
                all_COGs.push(data['COG'][j]);
            }
            if (data['score'] !== '') { 
                all_scores.push(data['score']);
            }
            if (data['mibig'] !== "") {
                total_mibig += 1;
            }
            if (data['filtered']) {
                total_filtered += 1
            }
            if ("kripp" in data && data["kripp"].length > 0) {
                total_kripp += 1
            }
        }
        // Calculate the new values
        res = standardDeviation(all_scores);
        average_score = Math.round(res[0]*100)/100;
        std_score = Math.round(res[1]*100)/100;
        res = standardDeviation(all_COGs);
        average_COG = Math.round(res[0]*100)/100;
        std_COG = Math.round(res[1]*100)/100;
        if (isNaN(average_score)) {
            score_text = "N\\A";
        } else {
            score_text = average_score + " +- " + std_score;
        }
        if (isNaN(average_COG)) {
            cog_text = "N\\A";
        } else {
            cog_text = average_COG + " +- " + std_COG;
        }
        // Get the relevant table row
        tr = document.getElementById(groupname);
        tds = tr.getElementsByTagName("td");
        tds[1].innerText = score_text;
        tds[2].innerText = operons_group.length;
        tds[3].innerText = total_mibig;
        tds[4].innerText = total_antismash;
        tds[5].innerText = total_kripp;
        tds[6].innerText = cog_text;
        tds[7].innerText = total_filtered;
        
    }
}

function GetCommonDomains() {
    var common_domains = {};
}

function HideRows(groups) {
    rows = document.getElementsByTagName("tr");
    if (Array.isArray(groups)) {
        groupnames = groups;
    } else {
        groupnames = Object.keys(groups);
    }
    for (i = 0; i < rows.length; i ++) {
        row = rows[i];
        if (groupnames.includes(row.id) === false && row.id !== "") {
            row.style.display = "none";
        }
        else {
            row.style.display = "";
        }
    }
}

function DisplayRequirements(requirements) {
    textbox = document.getElementById("requirements_textbox");
    text = '';
    for (req in requirements) {
        value_req = requirements[req]["value"];
        type_req = requirements[req]["type"];
        text += '<li>' + req + ' ' + type_req + ' ' + value_req.toString() + '</li>'
    }
    textbox.innerHTML = text;
}

function ToggleFilter(groups_dict,Operons,requirements) {
    HideRows(groups_dict);
    ChangeTableData(groups_dict,Operons);
    DisplayRequirements(requirements);
    localStorage.setItem("requirements", JSON.stringify(requirements));
}

function ApplyFilter(requirements,Operons,Groups) {
    ParseRequirements(requirements);
    if (jQuery.isEmptyObject(requirements)) {
        filtered_operons = base_operons;
    } else {
        filtered_operons = FilterOperons(Operons,requirements);
    };
    filtered_groups_dict = FilterGroups(Groups,filtered_operons);
    ToggleFilter(filtered_groups_dict,Operons,requirements);
}

function ResetFilter(Groups,Operons) {
    requirements = {};
    ToggleFilter(Groups,Operons,requirements);
}

