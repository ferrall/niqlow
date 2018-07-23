var sContent_url;
var sContent_href = "ox.html";
var sContent_hash = "";
var bIsWebKit = navigator.userAgent.indexOf('WebKit') > -1;

function setContent() {
    var s = location.search;
    var pos = s.indexOf("?content=");
    if (pos >= 0) {
        sContent_href = s.substring(pos + 9);
    }
    sContent_hash = location.hash;
    sContent_url = sContent_href + sContent_hash;
}
function loadDoc(url) {
	document.open("text/html","replace");
	document.write(
	'<!DOCTYPE html><html><frameset rows="50,*" border="0" framespacing="0" id="frameset1" framespacing="0">' +
	 '<frame scrolling="no" frameborder="0" name="toolbar" src="toolbar.html" marginwidth="20" marginheight="2" noresize>' +
	  '<frameset cols="200,*" border="0" frameborder="0" framespacing="0">' +
	   '<frame src="oxmenu.html" frameborder="0" name="index" scrolling="auto" marginheight="0" marginwidth="4" noresize>' +
	   '<frame src="' + url + '" frameborder="0" name="content" id="content" marginheight="12" marginwidth="12">' +
	  '<\/frameset>' +
	'<\/frameset></html>'
	);
	document.close();
}

// parse the url
setContent();

// hack: delay for webkit based browsers (Chrome, Safari)
// to ensure that the anchor is located properly
if (bIsWebKit && sContent_hash != "")
{
	loadDoc(sContent_href);
	setTimeout("loadDoc(sContent_url);", 30);
}
else
{
	loadDoc(sContent_url);
}
