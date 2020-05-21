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

//Get the requirements from the local storage

var requirements = localStorage.getItem("requirements")
if (requirements === null) {
    requirements = {}
} else {
    requirements = JSON.parse(requirements);
}

//Filter operons
groupname = document.getElementsByTagName("html_id")[0].innerText;
all_operons = Groups[groupname];
if (jQuery.isEmptyObject(requirements)) {
    filtered_operons = all_operons;
} else {
    filtered_operons = FilterOperonsEntry(Operons,all_operons,requirements);
    ChangeTableDataEntry(filtered_operons.length)
};
HideRows(filtered_operons);

