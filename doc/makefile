multigrid_documentation.pdf: multigrid_documentation.tex
	latex multigrid_documentation.tex
	latex multigrid_documentation.tex
	dvips -o multigrid_documentation.ps multigrid_documentation.dvi
	ps2pdf multigrid_documentation.ps
clean:
	rm -f *.aux *.log *.backup *.toc *.tex~ *.ps *.dvi

	
