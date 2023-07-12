
function check_version(docs_version) {
    fetch("https://staglibrary.io/docs/cpp/version.txt")
        .then(response => response.text())
        .then(data => {
            if (data.trim() !== docs_version) {
                // Create the warning message
                const warn_elem = document.createElement('div');
                warn_elem.innerHTML = "<dl class=\"section note\"><dt>Note</dt><dd>" +
                    "This documentation is not for the latest library version. " +
                    "Documentation for the latest version is available <a href='https://staglibrary.io/docs/cpp/'>here</a>." +
                    "</dd></dl>";

                // Find the contents element.
                const contents_elem = document.querySelector('.contents');

                // If the contents element has display: flex, then find
                // the first textblock in the content.
                if (contents_elem){
                    // Get the first child of the contents block
                    const first_child = contents_elem.firstElementChild;

                    if (first_child) {
                        if (first_child.className == 'textblock') {
                            // If the first child is a textblock element, insert the
                            // warning inside the textblock element
                            const next_child = first_child.firstChild;
                            first_child.insertBefore(warn_elem, next_child);
                        }  else if (first_child.className == 'toc') {
                            // If the first child is a toc, find the first textblock
                            // element and add there.
                            const textblock_elem = document.querySelector('.textblock');
                            if (textblock_elem) {
                                const next_child = textblock_elem.firstChild;
                                textblock_elem.insertBefore(warn_elem, next_child);
                            }
                        } else {
                            // The first child is not a textblock or toc - insert
                            // the warning directly into the content element.
                            contents_elem.insertBefore(warn_elem, first_child);
                        }
                    }
                }
            }
        })
}
