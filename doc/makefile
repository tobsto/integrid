#integrid_documentation.pdf: integrid_documentation.dvi integrid_documentation.ps integrid_documentation_*.tex integrid_documentation.bib
integrid_documentation.dvi: integrid_documentation_*.tex integrid_documentation.bib
	latex integrid_documentation
	bibtex integrid_documentation
	latex integrid_documentation.tex
	makeindex integrid_documentation.idx
	latex integrid_documentation.tex

pdf: integrid_documentation_*.tex integrid_documentation.bib
	latex integrid_documentation
	bibtex integrid_documentation
	latex integrid_documentation.tex
	makeindex integrid_documentation.idx
	latex integrid_documentation.tex
	dvips -o integrid_documentation.ps integrid_documentation.dvi
	ps2pdf integrid_documentation.ps
clean:
	rm -f *.out *.aux *.log *.bak *.backup *.toc *~ *.ps *.dvi *.bbl *.blg *.idx *.ilg *.ind
