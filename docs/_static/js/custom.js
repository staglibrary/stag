var dark_background = "rgb(26, 28, 30)";
var light_background = "rgb(248, 249, 251)";
var dark_color= "#9ca0a5";
var light_color= "#000";

var dom_event_fired = false;

document.addEventListener("DOMContentLoaded", function (){
    console.log("DOM");

    // If this is the first dom event fired, add html to the window
    if (!dom_event_fired) {
        dom_event_fired = true;

        update_header();

        // Fire another DOM event
        window.document.dispatchEvent(new Event("DOMContentLoaded", {
            bubbles: true,
            cancelable: true
        }));
    }
});

function update_header() {
    // Add custom html to the documentation pages.
    var target = document.body;
    var newElement = document.createElement("div");
    newElement.innerHTML = "<!-- ======= Header ======= -->\n" +
        "<header id=\"header\" class=\"fixed-top d-flex align-items-center\">\n" +
        "    <div class=\"container d-flex align-items-center justify-content-between\">\n" +
        "\n" +
        "        <div class=\"logo\">\n" +
        "            <h1 class=\"text-light\"><a href=\"https://staglibrary.io/docs/python/index.html\"><span>STAG Python</span></a></h1>\n" +
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
        "                                <li><a href=\"https://staglibrary.io/docs/cpp/index.html\">C++</a></li>\n" +
        "                                <li><a href=\"index.html\">Python</a></li>\n" +
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

    // Add favicon code to the header
    if (!location.hostname.includes("localhost")) {
        target = document.head;
        target.innerHTML += "<link rel=\"apple-touch-icon\" sizes=\"180x180\" href=\"/apple-touch-icon.png\">\n" +
            "<link rel=\"icon\" type=\"image/png\" sizes=\"32x32\" href=\"/favicon-32x32.png\">\n" +
            "<link rel=\"icon\" type=\"image/png\" sizes=\"16x16\" href=\"/favicon-16x16.png\">\n" +
            "<link rel=\"manifest\" href=\"/site.webmanifest\">\n" +
            "<link rel=\"mask-icon\" href=\"/safari-pinned-tab.svg\" color=\"#5bbad5\">\n" +
            "<meta name=\"msapplication-TileColor\" content=\"#2b5797\">\n" +
            "<meta name=\"theme-color\" content=\"#ffffff\">";
    }
}

window.onload = () => {
    // Set the header background to the side-bar background color
    $("#header").css('background', $(".sidebar-drawer").css('background'));

    var current_background = $(".sidebar-drawer").css('background');

    // Set the new background and colors accordingly
    if (current_background.includes(dark_background)) {
        $(".navbar a").css('color', dark_color);
    } else {
        $(".navbar a").css('color', light_color);
    }
}

function cycleHeaderThemeOnce() {
    const currentTheme = localStorage.getItem("theme") || "auto";
    const prefersDark = window.matchMedia("(prefers-color-scheme: dark)").matches;

    if (prefersDark) {
        // Auto (dark) -> Light -> Dark
        if (currentTheme === "auto") {
            $("#header").css('background', light_background);
            $(".navbar a").css('color', light_color);
        } else if (currentTheme == "light") {
            $("#header").css('background', dark_background);
            $(".navbar a").css('color', dark_color);
        } else {
            $("#header").css('background', dark_background);
            $(".navbar a").css('color', dark_color);
        }
    } else {
        // Auto (light) -> Dark -> Light
        if (currentTheme === "auto") {
            $("#header").css('background', dark_background);
            $(".navbar a").css('color', dark_color);
        } else if (currentTheme == "dark") {
            $("#header").css('background', light_background);
            $(".navbar a").css('color', light_color);
        } else {
            $("#header").css('background', light_background);
            $(".navbar a").css('color', light_color);
        }
    }
}

// Update the header colors on toggle
$(".theme-toggle").click(function() {
    cycleHeaderThemeOnce();
});
