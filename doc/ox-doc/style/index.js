var html_main = "main.html";
var html_toc = "";
var html_idx = "";

var pos_width = ($("#menubar").width() + 2);
var neg_width = -pos_width;

function loadContent(defaultSrc) {
    var src;
    var s = location.search;
    var pos = s.indexOf("?content=");
    if (pos >= 0) {
        src = s.substring(pos + 9) + location.hash;
    }
	else {
		src = defaultSrc || html_main;
	}
    $("#content-frame").attr("src", src);
}

function toggleMenubarStatus() {
	$("#menubar").toggleClass("active");
}
function toggleMenubar() {
	if ($("#menubar").hasClass("active")) {
		$("#menubar").show();
		$("#header").css("left", pos_width);
		$("#content").css("left", pos_width);
	}
	else {
		$("#menubar").hide();
		$("#header").css("left", 0);
		$("#content").css("left", 0);
	}
}
function pushMenubar() {
	$("#content").width($(window).width());
	if ($("#menubar").hasClass("active")) {
		$("#menubar").show();
		$("#header").css("marginLeft", 0);
		$("#content").css("marginLeft", 0);
		$("#navicon-text").text("x");
		$("#navicon-text").css("font-weight","Bold");
	}
	else {
		$("#menubar").hide();
		$("#header").css("marginLeft", neg_width);
		$("#content").css("marginLeft", neg_width);
		$("#navicon-text").text("Menu");
		$("#navicon-text").css("font-weight","Normal");
	}
}
function showMenubar() {
	if ($("#menubar").hasClass("push")) {
		pushMenubar();
	}
	else {
		toggleMenubar();
	}
}

function initContent() {
	if ($(window).width() < 500) {
		$("#menubar").addClass("push");
		$("#menubar").removeClass("active");
	}
	else if ($(window).width() < 700) {
		$("#menubar").removeClass("push");
		$("#menubar").removeClass("active");
	}
	else {
		$("#menubar").addClass("active");
	}
	
	$(window).resize(function() {
		if ($("#menubar").hasClass("push")) {
			$("#content").width($(window).width());
		}
	});
	
	showMenubar();
	loadContent();

	if ($("#load-toc").length) {
		$("#load-toc").on("click", function(event, ui) { $("#menu-frame").attr("src", html_toc); } );
	}
	if ($("#load-idx").length) {
		$("#load-idx").on("click", function(event, ui) { $("#menu-frame").attr("src", html_idx);  } );
	}
	$("#navicon").on("click", function(event, ui) { toggleMenubarStatus(); showMenubar(); } );
}
function setContent(main, toc, idx) {
	html_main = main;
	html_toc = toc;
	html_idx = idx;
	initContent();
}
