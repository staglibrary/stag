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
    # Since Furo doesn't allow us to disable dark mode, we make dark mode
    # equivalent to light mode by overriding all colors back to their light value.
    # See: https://github.com/pradyunsg/furo/issues/28
    'dark_css_variables': {
        # Taken from: https://github.com/pradyunsg/furo/blob/c682d5d3502f3fa713c909eebbf9f3afa0f469d9/src/furo/assets/styles/variables/_colors.scss
        'color-problematic': '#b30000',

        # Base Colors
        'color-foreground-primary': 'black', # for main text and headings
        'color-foreground-secondary': '#5a5c63', # for secondary text
        'color-foreground-muted': '#646776', # for muted text
        'color-foreground-border': '#878787', # for content borders

        'color-background-primary': 'white', # for content
        'color-background-secondary': '#f8f9fb', # for navigation + ToC
        'color-background-hover': '#efeff4ff', # for navigation-item hover
        'color-background-hover--transparent': '#efeff400',
        'color-background-border': '#eeebee', # for UI borders
        'color-background-item': '#ccc', # for "background" items (eg: copybutton)

        # Announcements
        'color-announcement-background': '#000000dd',
        'color-announcement-text': '#eeebee',

        # Brand colors
        'color-brand-primary': '#2962ff',
        'color-brand-content': '#2a5adf',

        # Highlighted text (search)
        'color-highlighted-background': '#ddeeff',

        # GUI Labels
        'color-guilabel-background': '#ddeeff80',
        'color-guilabel-border': '#bedaf580',

        # API documentation
        'color-api-keyword': 'var(--color-foreground-secondary)',
        'color-highlight-on-target': '#ffffcc',

        # Admonitions
        'color-admonition-background': 'transparent',

        # Cards
        'color-card-border': 'var(--color-background-secondary)',
        'color-card-background': 'transparent',
        'color-card-marginals-background': 'var(--color-background-hover)',

        # Code blocks
        'color-code-foreground': 'black',
        'color-code-background': '#f8f9fb',
    }
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
