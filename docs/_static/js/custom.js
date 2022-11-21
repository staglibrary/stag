window.onload = () => {
    // Add custom html to the documentation pages.
    var target = document.body;
    var newElement = document.createElement("div");
    newElement.innerHTML = "<!-- ======= Header ======= -->\n" +
        "<header id=\"header\" class=\"fixed-top d-flex align-items-center\">\n" +
        "    <div class=\"container d-flex align-items-center justify-content-between\">\n" +
        "\n" +
        "        <div class=\"logo\">\n" +
        "            <h1 class=\"text-light\"><a href=\"https://staglibrary.io/index.html\"><span>STAG Library</span></a></h1>\n" +
        "            <!-- Uncomment below if you prefer to use an image logo -->\n" +
        "            <!-- <a href=\"index.html\"><img src=\"assets/img/logo.png\" alt=\"\" class=\"img-fluid\"></a>-->\n" +
        "        </div>\n" +
        "\n" +
        "        <nav id=\"navbar\" class=\"navbar\">\n" +
        "            <ul>\n" +
        "                <li><a class=\"nav-link scrollto\" href=\"https://staglibrary.io/index.html\">Home</a></li>\n" +
        "                <li class=\"dropdown\"><a href=\"#\"><span>Components</span> <i class=\"bi bi-chevron-down\"></i></a>\n" +
        "                    <ul>\n" +
        "                        <li class=\"dropdown\"><a href=\"\"><span> Clustering</span> <i\n" +
        "                                class=\"bi bi-chevron-right\"></i></a>\n" +
        "                            <ul>\n" +
        "                                <li><a href=\"https://staglibrary.io/index.html#local\">Local Graph Clustering</a></li>\n" +
        "                                <li><a href=\"https://staglibrary.io/clustering-demo.html\">Demo</a></li>\n" +
        "                                <li>\n" +
        "                            </ul>\n" +
        "                        </li>\n" +
        "                    </ul>\n" +
        "                </li>\n" +
        "                <li class=\"dropdown active\"><a href=\"#\"><span>Documentation</span> <i class=\"bi bi-chevron-down\"></i></a>\n" +
        "                            <ul>\n" +
        "                                <li><a href=\"index.html\">C++</a></li>\n" +
        // "                                <li><a href=\"https://staglibrary.io/docs/python\">Python</a></li>\n" +
        "                                <li>\n" +
        "                            </ul>\n" +
        "                </li>\n" +
        "                <li><a class=\"nav-link scrollto\" href=\"https://staglibrary.io/index.html#about\">About Us</a></li>\n" +
        "                <li><a class=\"nav-link scrollto\" href=\"https://staglibrary.io/index.html#faq\">FAQ</a></li>\n" +
        "                <li><a class=\"nav-link scrollto\" href=\"https://staglibrary.io/index.html#contact\">Contact</a></li>\n" +
        "            </ul>\n" +
        "            <i class=\"bi bi-list mobile-nav-toggle\"></i>\n" +
        "        </nav><!-- .navbar -->\n" +
        "\n" +
        "    </div>\n" +
        "</header><!-- End Header -->\n"

    target.insertBefore(newElement, target.firstChild);

}
