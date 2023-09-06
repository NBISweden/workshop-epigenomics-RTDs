# Contributing teaching material

The website is rendered to <https://nbis-workshop-epigenomics.readthedocs.io/en/latest/> and build with Sphinx. To add or modify the content:

## Editing workflow

### Clone this repository

``` bash
# clone the repo
git clone https://github.com/NBISweden/workshop-epigenomics-RTDs
cd workshop-epigenomics-RTDs

# checkout feature branch to work on
git checkout -b session-example
```

### Edit files

Add and edit files at will and note that:

* **Supported formats**: `.rst`, `.md`
* **Tutorial files**: go under `docs/content/tutorials`
* **Links to tutorials**: add to  `docs/content/tutorials/*.rst`

### Push feature branch to repo

```bash
# code & commit changes while working on the materials
git add session-feature.md
git commit -m "commit message"

# push to feature when ready
git push
```

#### Git commits good practices

* Commit messages should contain relevant information regarding the feature(s) you add, what type of analyses they can be used for, *etc.*.
* The subject line should be written in an imperative, e.g. *Fix typos* and be 50 characters or less
* More about [good commit messages][git-commits]

### Make a pull request to master branch when ready

* Go to course repository [https://github.com/NBISweden/workshop-epigenomics-RTDs](https://github.com/NBISweden/workshop-epigenomics-RTDs) and create a pull request.
* **OK** to merge with the main branch, if you know what you're doing and have technical permission to do so.

> Readthedocs will automatically build the site upon new commits

## Viewing changes locally

You can view changes locally while working on. To do so:

### Install Sphinx and readthedocs theme

```bash
conda install sphinx
pip install sphinx-rtd-theme
pip install recommonmark      # to support Markdown
pip install sphinx_togglebutton
pip install sphinx_copybutton
```
Or via environment:

```bash
conda create -n sphinx7 -c conda-forge python=3.11 pip sphinx=7.2.4 sphinx_rtd_theme=1.3.0 docutils sphinx-togglebutton sphinx-copybutton pygments recommonmark 
```

Versions currently used:

```tsv
sphinx                    7.2.4
sphinx-rtd-theme          1.3.0
sphinx-copybutton         0.5.2
sphinx-togglebutton       0.3.2
pygments                  2.16.1
docutils                  0.18.1
```



### Build the local html

In project root (i.e. `/docs`):

```bash
make html
```

The html to view is `_build/html/index.html`

## More of file formats

### Markdown `.md`

It is supported. However, the tables are not rendered properly in html.

### Restructured text `.rst`

The preferred format. A useful primer on syntax:

<https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>

### Issues

Contact Agata

[git-commits]: https://chris.beams.io/posts/git-commit/
