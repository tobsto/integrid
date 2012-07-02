#multigrid_documentation.pdf: multigrid_documentation.dvi multigrid_documentation.ps multigrid_documentation_*.tex multigrid_documentation.bib
multigrid_documentation.dvi: multigrid_documentation_*.tex multigrid_documentation.bib
	latex multigrid_documentation
	bibtex multigrid_documentation
	latex multigrid_documentation.tex
	makeindex multigrid_documentation.idx
	latex multigrid_documentation.tex

pdf: multigrid_documentation.pdf multigrid_documentation_*.tex multigrid_documentation.bib
	latex multigrid_documentation
	bibtex multigrid_documentation
	latex multigrid_documentation.tex
	makeindex multigrid_documentation.idx
	latex multigrid_documentation.tex
	dvips -o multigrid_documentation.ps multigrid_documentation.dvi
	ps2pdf multigrid_documentation.ps
clean:
	rm -f *.out *.aux *.log *.bak *.backup *.toc *~ *.ps *.dvi *.bbl *.blg *.idx *.ilg *.ind
