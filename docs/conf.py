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

html_theme = 'furo'
html_static_path = ['_static']
html_logo = '_static/img/logo/ginger.png'
html_theme_options = {
    "sidebar_hide_name": True,
}

# Add all vendor CSS and javascript
html_css_files = []
html_js_files = []
html_css_files.append('vendor/bootstrap-icons/bootstrap-icons.css')

# Finally, add custom scripts for the documentation
html_css_files.append('css/custom.css')
html_css_files.append('css/style.css')
html_js_files.append('js/custom.js')
