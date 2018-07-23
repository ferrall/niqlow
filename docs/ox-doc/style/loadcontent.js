function loadContent(defaultSrc) {
    var src;
    var s = location.search;
    var pos = s.indexOf("?content=");
    if (pos >= 0) {
        src = s.substring(pos + 9) + location.hash;
    }
	else {
		src = defaultSrc || "main.html";
	}
    document.getElementById("content-frame").src = src;
}
