# latex Makefile
ifndef texpath
texpath=/usr/texbin/
endif
PDFLATEX=${texpath}pdflatex -halt-on-error -synctex=1 --interaction=nonstopmode
SKIPERR=${texpath}pdflatex --interaction=nonstopmode
LATEX=${PDFLATEX}
BIBTEX=bibtex
DVIPS=dvips
PS2PDF=ps2pdf
SHELL=/bin/bash

all: gitstuff.tex cfe

.PHONY: cfe
cfe: 
	echo "texpath: ${texpath}"
	#python make_apjform.py
	python make.py --texpath=${texpath}
	gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=cfe.pdf cfe_compressed.pdf
	python paper2arxiv.py --apj
	python paper2arxiv.py

cfe.tex: cfe

.PHONY: diff
diff:
	python parse_macros.py cfe.tex cfe.tex
	python parse_macros.py submitted.tex original_for_diff.tex
	#python parse_macros.py cfe.tex cfe.tex
	#python2 parse_macros.py cfe.tex cfe.tex
	latexdiff original_for_diff.tex cfe.tex > diff.tex
	${SKIPERR} diff.tex
	${BIBTEX} diff
	${SKIPERR} diff.tex
	${BIBTEX} diff
	${SKIPERR} diff.tex
	gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=diff_compressed.pdf diff.pdf

.PHONY: referee
referee:
	echo "texpath: ${texpath}"
	python make.py --referee --texpath=${texpath}


gitstuff.tex: ../.git/logs/HEAD
	echo "%%% This file is generated by the Makefile." > gitstuff.tex
	git log -1 --date=short --format="format:\\newcommand{\\githash}{%h}\\newcommand{\\gitdate}{%ad\\xspace}\\newcommand{\\gitauthor}{%an\\xspace}" >> gitstuff.tex

