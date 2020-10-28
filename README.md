# nbis-epigenomics-workshop

This is a repository for the NBIS Epigenomics Workshop website.

Current `rst` version in `master`; older `mkdocs` version in a separate branch




## Modifications and contributions

To modify the webpage:

1. Clone this repo;

3. Edit the documents, add files etc. If the changes are to be viewed locally, see below.

4. Commit and push to this repo (to master). Readthedocs will automatically build the site upon new commits. In case of issues, contact Agata or Olga for building the site on Readthedocs.


### Guidelines for contributions

Please keep the tutorials content in their own folder under `docs/content/tutorials`. Please add the appropriate link to appropriate contents file in `docs/content/tutorials/*.rst`.


### Viewing changes locally


#### Install Sphix and readthedocs theme

```
conda install sphinx
pip install sphinx-rtd-theme

# to support Markdown
pip install recommonmark
```

Versions currently used:

```
sphinx                    3.2.1
sphinx-rtd-theme          0.5.0
```

#### Edit the docs




#### Build the local html

In project root (i.e. `/docs`):

```
make html
```

The html to view is `_build/html/index.html`


## Formats

### Markdown `.md`

It is supported. However, the tables are not rendered properly in html.


### Restructured text `.rst`


A useful primer on syntax:

https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html


And... Better tables (finally)!



