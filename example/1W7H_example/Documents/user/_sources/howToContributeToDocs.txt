######################################
How to contribute to the documentation
######################################

If you want to contribute to the documentation of PELE, you are welcome.

PELE documentation uses the `Sphinx <http://sphinx-doc.org/>`__ tool, which is based in the `reStructuredText <http://docutils.sf.net/rst.html>`__ markup language (very similar to the `Markdown <http://daringfireball.net/projects/markdown/>`__ markup language previously used in PELE). A quick introduction to reStructuredText with Sphinx is at http://sphinx-doc.org/rest.html

For any page where you will contribute, you can obtain the original document for that page on the bar on the left, following the `This Page > Show Source`. Download that file, make the modifications you want, and send it to the PELE support team, which will review it and merge it with the official documentation as suitable.

If you want to add a completely new page, proceed as above, but build the complete page structure, using one of the existing pages as a template. Then send it to the PELE support team and tell them where in the documentation structure the page fits.

Below are some general instructions on how to write the documentation. For all the details about the Sphinx tool, check `Shpinx <http://sphinx-doc.org/>`__.

****************************************
Notes on writing a page of documentation
****************************************

Structure
=========

Follow the sectioning structure that Sphinx uses as a convention (http://sphinx-doc.org/rest.html#sections). This is:

* ``#`` with overline, for parts
* ``*`` with overline, for chapters
* ``=``, for sections
* ``-``, for subsections
* ``^``, for subsubsections
* ``"``, for paragraphs

Depending on where your documentation page fits into the PELE documentation structure, your page top level header can be any of the previous ones.

Cross-references
================

To help the user connecting ideas, and also to avoid redundance by keeping the information in a single place, use cross-references (see http://sphinx-doc.org/markup/inline.html#cross-referencing-arbitrary-locations).

Labels for section will start with ``sec-``, while labels for figures start with ``fig-`` and labels for tables with ``table-``.

For sections, the main header of a page will take the name of the file, but with the first character in lower case (unless the file starts with an acronym as PELE). For example, the label for the main header of file :file:`SelectionExamples.rst` is ``sec-selectionExamples``; however, for :file:`PELESimulation.rst`, it is ``sec-PELESimulation``. Subsections in a file will start with the same label as the main header of the file, plus additional text, using ``-`` as a separator. For example, the label for the ``Parameters`` section in the :file:`Anm.rst` file is ``sec-anm-parameters``.

Footnotes
=========

The text corresponding to footnotes will appear in the same page as the footnote. See Sphinx for reference: http://sphinx-doc.org/rest.html#footnotes

Footnotes will be labeled as ``#f<footnote_number>`` consecutively in the file. The footnote section will be defined with the `rubric` directive. For example:

.. code-block:: rest

    A peleSimulation has a PeleTask block. It is an array comprised of all
    the tasks that the program has to carry out before a command (usually a
    simulation) is considered to be finished [#f1]_. Each task is a block, so it
    should be enclosed in curly braces.

    .. rubric:: Footnotes

    .. [#f1] It can also finish because numberOfPeleSteps is reached.

Citations
=========

Each part of the documentation (for example, the ``Reference Manual`` or the ``Tutorial``) may have its own references, which should be gathered in a page of their own. Each reference will be labeled by ``[firstauthor:year]``, and follows the format: Authors "Article title" *Journal* volume, pp. 1-n, Year. URL.

For example:

.. code-block:: rest

    References
    ==========

    .. [atilgan:2001] Atilgan, A.R. and Durell, S.R. and Jernigan, R.L. and Demirel, M.C. and Keskin, O. and Bahar, I. "Anisotropy of Fluctuation Dynamics of Proteins with an Elastic Network Model" *Biophysical Journal* 80, pp. 505--515, 2001. http://dx.doi.org/10.1016/S0006-3495(01)76033-X

For citing those references, follow the instructions given by Sphinx (http://sphinx-doc.org/rest.html#citations). For example, the previous paper will be cited as:

.. code-block:: rest

  The original form of the potential ([atilgan:2001]_) is ...

Mathematics
===========

The documentation uses MathJax to display mathematical formulas, symbols and equations (see http://sphinx-doc.org/ext/math.html). This is basically LaTeX input, though only basic mathematical commands are available (for a complete documentation, see `MathJax <https://www.mathjax.org/>`__).

For example, the Pythagoras equation, inline, would be coded as:

.. code-block:: rest

  Since Pythagoras, we know that :math:`a^2 + b^2 = c^2`.

And it will appear as :math:`a^2 + b^2 = c^2`.

Though the ``\AA`` command is not defined in MathJax, the PELE documentation configures MathJax (in :file:`_templates/layout.html`) so that you can use it to represent :math:`\AA{}`.

Notice that the HTML representation of mathematics through MathJax requires an Internet connection, since the MathJax Javascript code is always retrieved from the Internet.

File names
==========

When providing a file name or path, use the ``:file:`` role of Sphinx (http://sphinx-doc.org/markup/inline.html#role-file):

.. code-block:: rest

  Open file :file:`README.txt`.


Examples and source code
========================

For displaying the contents (or extacts) of files, possibly source code files, use the ``code-block`` directive (http://sphinx-doc.org/markup/code.html#showing-code-examples). For example, to show python code, use:

.. code-block:: rest

  .. code-block:: python

    for i in [1,2,3]:
       print i

Which results in:

.. code-block:: python

  for i in [1,2,3]:
    print i

Hyperlinks
==========

For referencing external links, use the inline format given in http://sphinx-doc.org/rest.html#hyperlinks

Notice that external links will appear with discontinued underlining in the current PELE theme, but internal cross-references will simply appear as colored cursive text, with no underlining.

Copyright notice
================

Since the copyright notice is included in the general template for all pages, you do not have to add any copyright text to the documentation pages. Sphinx will take care of it.

