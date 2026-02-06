
See the link:

https://github.com/lsst
https://github.com/LSSTDESC

for examples of how to document.


See for Sphinx documentation making:

https://www.sphinx-doc.org/en/master/tutorial/index.html
https://www.youtube.com/watch?v=nZttMg_n_s0
https://docs.readthedocs.io/en/stable/intro/getting-started-with-sphinx.html
https://sphinx-themes.org/

To create a html version of man page do from cTreeBalls directory:

man doc/cballs.m | man2html -topm 0 -botm 0 -cgiurl \$title.html doc/cballs.m > doc/cballs.html

If we follow tutorial en youtube:

sphinx-quickstart

$ pip install sphinx

Install additional themes or MkDocs
pip install sphinx-rtd-theme
pip install mkdocs


Other useful links

https://sphinx-tutorial.readthedocs.io/step-1/
https://github.com/kiith-sa/RestructuredText-tutorial/tree/master

$ make

$ make html
$ make latexpdf
$ make man
$ make 
