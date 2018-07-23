
if (self.location.href == top.location.href) {
    top.location.replace("index.html?content=" + location.href);
    document.close();
}
else {
    top.document.title = document.title;
}
