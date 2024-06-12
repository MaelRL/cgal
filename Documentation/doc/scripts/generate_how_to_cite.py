#!/usr/bin/env python3
# coding: utf8

import re
import codecs
from sys import argv
from sys import stderr
import os

### Constants ###

RESULT_BIB = r"""@book{ cgal:eb-${CGAL_RELEASE_YEAR_ID}
, title        = "{CGAL} User and Reference Manual"
, author      = "{The CGAL Project}"
, publisher     = "{CGAL Editorial Board}"
, edition     = "{${CGAL_CREATED_VERSION_NUM}}"
, year         = ${CGAL_BUILD_YEAR4}
, url =    "https://doc.cgal.org/${CGAL_CREATED_VERSION_NUM}/Manual/packages.html"
}"""

RESULT_TXT = r"""// This file was generated by generate_how_to_cite.py. You shouldn't modify it directly.

/*!
\page how_to_cite_cgal Acknowledging %CGAL

\details \cgal is implemented for the most part by researchers. The
academic world evaluates and rewards researchers for a good part by
the analysis of the number of published papers and the number of
citations of their papers, which measures their impact. In order to
make the \cgal project attractive for researchers to contribute their
work (which allows users to benefit from new contributions), we are
pushing a model where the \cgal manual chapters are considered like
publications, and can be cited in articles as such.

We therefore kindly ask users to cite \cgal as appropriately as
possible in their papers, and to mention the use of \cgal on the web
pages of their projects using \cgal and provide us with links to these
web pages. Feel free to contact us in case you have any question or
remark on this topic.

We provide bibtex entries for the chapters of the User and Reference
Manual, as well as for publications directly related to the \cgal
software.

## Citing the %CGAL Library or the %CGAL project ##

If you want to cite the \cgal Library or project as a whole, please

- cite: \cgal, Computational Geometry Algorithms Library, https://www.cgal.org
- use the first bibtex entry from the file <a href="how_to_cite_cgal.bib">how_to_cite_cgal.bib</a>.

## Citing the User and Reference Manual ##

If you want to refer to \cgal manual, please cite the appropriate
  entry from the bibliographic entries for individual chapters listed
  in the table below.

<table>

<tr valign="top">
<td align="right" class="bibtexnumber">
[<a name="cgal:eb-${CGAL_RELEASE_YEAR_ID}">1</a>]
</td>
<td class="bibtexitem">
The \cgal Project.
 <em>\cgal User and Reference Manual</em>.
 \cgal Editorial Board, ${CGAL_CREATED_VERSION_NUM} edition, ${CGAL_BUILD_YEAR4}.
[&nbsp;<a href="how_to_cite.html#cgal:eb-${CGAL_RELEASE_YEAR_ID}">bib</a>&nbsp;|
<a href="packages.html">http</a>&nbsp;]

</td>
</tr>


"""
RESULT_TXT_FOOTER = r"""</td>
</tr>
</table><hr>
*/
"""

PRE_HTML = r"""<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "https://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<link rel="icon" type="image/png" href="../Manual/g-196x196-doc.png"/>
<meta http-equiv="Content-Type" content="text/xhtml;charset=UTF-8"/>
<meta http-equiv="X-UA-Compatible" content="IE=9"/>
<link href="cgal_stylesheet.css" rel="stylesheet" type="text/css" />
<title>CGAL ${CGAL_CREATED_VERSION_NUM} - Manual: Acknowledging CGAL</title>
</head>
<body>
"""
POST_HTML = r"""</body>
</html>
"""

RESULT_HTML = r"""<h1>how_to_cite_cgal.bib</h1><a name="cgal:eb-${CGAL_RELEASE_YEAR_ID}"></a><pre>
@book{<a href="how_to_cite_cgal.html#cgal:eb-${CGAL_RELEASE_YEAR_ID}">cgal:eb-${CGAL_RELEASE_YEAR_ID}</a>,
  title = {{CGAL} User and Reference Manual},
  author = {{The CGAL Project}},
  publisher = {{CGAL Editorial Board}},
  edition = {{${CGAL_CREATED_VERSION_NUM}}},
  year = ${CGAL_BUILD_YEAR4},
  url = {<a href="https://doc.cgal.org/${CGAL_CREATED_VERSION_NUM}/Manual/packages.html">https://doc.cgal.org/${CGAL_CREATED_VERSION_NUM}/Manual/packages.html</a>}
}
</pre>

"""

### Functions ###


def gen_bib_entry(title, authors, bib, anchor):
    res = (
        "\n\
@incollection{"
        + bib
        + '-${CGAL_RELEASE_YEAR_ID}\n\
, author =  "'
        + authors
        + '"\n\
, title =   "'
        + title
        + '"\n\
, publisher =  "{CGAL Editorial Board}"\n\
, edition =     "{${CGAL_CREATED_VERSION_NUM}}"\n\
, booktitle =   "{CGAL} User and Reference Manual"\n\
, url = "https://doc.cgal.org/${CGAL_CREATED_VERSION_NUM}/Manual/packages.html#'
        + anchor
        + '"\n\
, year =        ${CGAL_BUILD_YEAR4}\n\
}\n'
    )
    return res


def gen_txt_entry(title, authors, bib, anchor, k):
    title_r = (
        title.replace("Kernel", "%Kernel")
        .replace("Interval", "%Interval")
        .replace("Matrix", "%Matrix")
        .replace("Kinetic", "%Kinetic")
        .replace("CGAL", r"\cgal")
        .replace("Range", "%Range")
    )
    authors = authors.replace("CGAL", r"\cgal")
    res = (
        '<tr valign="top">\n\
<td align="right" class="bibtexnumber">\n\
[<a name="'
        + bib
        + '-${CGAL_RELEASE_YEAR_ID}">'
        + str(k)
        + '</a>]\n\
</td>\n\
<td class="bibtexitem">\n '
        + authors
        + ".\n "
        + title_r
        + '.\n\
 In <em>\\cgal User and Reference Manual</em>. \\cgal Editorial Board,\n\
  ${CGAL_CREATED_VERSION_NUM} edition, ${CGAL_BUILD_YEAR4}.\n\
[&nbsp;<a href="how_to_cite.html#'
        + bib
        + '-${CGAL_RELEASE_YEAR_ID}">bib</a>&nbsp;| \n\
<a href="packages.html#'
        + anchor
        + '">http</a>&nbsp;]\n\
\n\
</td>\n\
</tr>\n\n\n'
    )
    return res


def gen_html_entry(title, authors, bib, anchor):
    res = (
        '<a name="'
        + bib
        + '-${CGAL_RELEASE_YEAR_ID}"></a><pre>\n\
@incollection{<a href="how_to_cite_cgal.html#'
        + bib
        + '-${CGAL_RELEASE_YEAR_ID}">'
        + bib
        + "-${CGAL_RELEASE_YEAR_ID}</a>,\n\
  author = {"
        + authors
        + "},\n\
  title = {"
        + title
        + '},\n\
  publisher = {{CGAL Editorial Board}},\n\
  edition = {{${CGAL_CREATED_VERSION_NUM}}},\n\
  booktitle = {{CGAL} User and Reference Manual},\n\
  url = {<a href="https://doc.cgal.org/${CGAL_CREATED_VERSION_NUM}/Manual/packages.html#'
        + anchor
        + '">https://doc.cgal.org/${CGAL_CREATED_VERSION_NUM}/Manual/packages.html#'
        + anchor
        + "</a>},\n\
  year = ${CGAL_BUILD_YEAR4}\n\
}\n\
</pre>\n\n"
    )
    return res


def protect_upper_case(title):
    return (
        title.replace("dD", "{dD}")
        .replace("2D", "{2D}")
        .replace("3D", "{3D}")
        .replace("CGAL", "{CGAL}")
        .replace("Qt", "{Qt}")
        .replace("Boost", "{Boost}")
    )


def protect_accentuated_letters(authors):
    res = (
        authors.replace("é", r"{\'e}")
        .replace("è", r"{\`e}")
        .replace("É", r"{\'E}")
        .replace("ä", r"{\"a}")
        .replace("ö", r"{\"o}")
        .replace("ñ", r"{\~n}")
        .replace("ã", r"{\~a}")
        .replace("ë", r"{\"e}")
        .replace("ı", r"{\i}")
        .replace("Ş", r"{\c{S}}")
        .replace("ş", r"{\c{s}}")
        .replace("%", "")
        .replace("đ", r"{\-d}")
    )
    try:
        res.encode("ascii")
    except UnicodeEncodeError:
        stderr.write(
            "WARNING: a non-ASCII character has been found in author string for bibtex (probably a non-handled accentuated letter)."
            "Check the new package added and update the function protect_accentuated_letters in Documentation/scripts/generate_how_to_cite.py\n\n"
        )
    return res


def make_doc_path(pkg, arg, source_dir, branch_build):
    if branch_build:
        return os.path.join(source_dir, pkg, "doc", pkg, arg)
    else:
        return os.path.join(source_dir, "doc", pkg, arg)

    ### Start of the main function ###


def main():
    assert (
        len(argv) == 4
    ), "require exactly three arguments: source_dir, build_dir, branch_build"

    source_dir = argv[1]
    build_dir = argv[2]
    branch_build = argv[3]

    result_txt = RESULT_TXT
    result_bib = RESULT_BIB
    result_html = RESULT_HTML

    pattern = re.compile(r"\\package_listing{([^}]*)}")

    pattern_title_and_anchor = re.compile(
        r"\\cgalPkgDescriptionBegin{([^}]*),\s?([^}]*)}"
    )
    pattern_author = re.compile(r"\\cgalPkgAuthors?{([^}]*)}")
    pattern_bib = re.compile(r"\\cgalPkgBib{([^}]*)}")

    f = codecs.open(
        make_doc_path("Documentation", "packages.txt", source_dir, branch_build),
        "r",
        encoding="utf-8",
    )
    k = 2
    for line in f:
        match = pattern.match(line)
        if match:
            pkg = match.group(1)
            filename = make_doc_path(
                pkg, "PackageDescription.txt", source_dir, branch_build
            )
            pkgdesc = codecs.open(filename, "r", encoding="utf-8")
            authors = ""
            bib = ""
            anchor = ""
            for pkg_line in pkgdesc:
                match = pattern_title_and_anchor.match(pkg_line)
                if match:
                    title = match.group(1).replace(r"\,", ",")
                    anchor = match.group(2)
                    continue
                match = pattern_author.match(pkg_line)
                if match:
                    authors = (
                        match.group(1).replace(", and", " and").replace(",", " and")
                    )
                    continue
                match = pattern_bib.match(pkg_line)
                if match:
                    bib = match.group(1)
                    continue
            assert len(bib) > 0, "Did you forget a \\cgalPkgBib{} in %r?" % filename
            assert len(authors) > 0, (
                "Did you forget a \\cgalPkgAuthors{} in %r?" % filename
            )
            assert len(anchor) > 0, (
                "Did you forget the anchor in \\cgalPkgDescriptionBegin{} in %r?"
                % filename
            )
            result_txt += gen_txt_entry(title, authors, bib, anchor, k)
            # convert title and author to bibtex format
            title = protect_upper_case(title)
            authors = protect_accentuated_letters(authors)
            result_bib += gen_bib_entry(title, authors, bib, anchor)
            result_html += gen_html_entry(title, authors, bib, anchor)
            k += 1

    f = codecs.open(build_dir + "/how_to_cite_cgal.bib.in", "w", encoding="utf-8")
    f.write(result_bib)
    result_txt += RESULT_TXT_FOOTER
    f = codecs.open(build_dir + "/how_to_cite_cgal.txt.in", "w", encoding="utf-8")
    f.write(result_txt)
    f = codecs.open(build_dir + "/how_to_cite.html.in", "w", encoding="utf-8")
    f.write(PRE_HTML + result_html + POST_HTML)


if __name__ == "__main__":
    main()
