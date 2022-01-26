kseHover = function() {
	var sfEls = document.getElementById("horznav").getElementsByTagName("LI");
	for (var i=0; i<sfEls.length; i++) {
		sfEls[i].onmouseover=function() {
			this.className+=" ksehover";
		}
		sfEls[i].onmouseout=function() {
			this.className=this.className.replace(new RegExp(" ksehover\\b"), "");
		}
	}
}
if (window.attachEvent) window.attachEvent("onload", kseHover);