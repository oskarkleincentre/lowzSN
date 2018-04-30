name=references
refs :
	pdflatex $(name)
	bibtex $(name)
	pdflatex $(name)
	pdflatex $(name)

clean:
	rm -f *.aux *.log *.pdf *out *synctex.gz *.blg *.bbl

default:  refs
