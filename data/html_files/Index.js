function opentab(evt, id) {
// Declare all variables
var i, tabcontent, tablinks;

// Get all elements with class="tabcontent" and hide them
tabcontent = document.getElementsByClassName("tabcontent");
for (i = 0; i < tabcontent.length; i++) {
    tabcontent[i].style.display = "none";
}

// Get all elements with class="tablinks" and remove the class "active"
tablinks = document.getElementsByClassName("tablinks");
for (i = 0; i < tablinks.length; i++) {
    tablinks[i].className = tablinks[i].className.replace(" active", "");
}

// Show the current tab, and add an "active" class to the button that opened the tab
document.getElementById(id).style.display = "block";
evt.currentTarget.className += " active";
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

//Set the initial requirements for the filters (load requirements if on the system)
var requirements = localStorage.getItem("requirements")
if (requirements === null) {
    var requirements = {}
} else {
    var requirements = JSON.parse(requirements);
    localStorage.setItem("requirements", JSON.stringify(requirements));
}
var base_groups = Object.keys(Groups);
var base_operons = Object.keys(Operons);

if (! jQuery.isEmptyObject(requirements) ) {
    ApplyFilter(requirements,Operons,Groups);
}
