# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'STAG'
copyright = '2022, The STAG Team'
author = 'The STAG Team'
release = '0.3.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['breathe']

# Breathe Configuration
breathe_default_project = "STAG"

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_logo = '_static/img/logo/ginger.png'
html_theme_options = {
    'logo_only': True,
    'sticky_navigation': False,
}

# Add all vendor CSS and javascript
import os
html_css_files = []
html_js_files = []
for root, dirs, files in os.walk("_static/vendor"):
    for file in files:
        if file.endswith(".css"):
            html_css_files.append(os.path.join(root[8:], file))
        if file.endswith(".js"):
            html_js_files.append(os.path.join(root[8:], file))

# Finally, add custom scripts for the documentation
html_css_files.append('css/custom.css')
html_css_files.append('css/style.css')
html_js_files.append('js/custom.js')

# Needed to fix the scrolling bug
html_js_files.append('https://cdn.jsdelivr.net/npm/jquery.scrollto@2.1.3/jquery.scrollTo.min.js')

