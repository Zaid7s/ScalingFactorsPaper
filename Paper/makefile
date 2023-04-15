DOCNAME=ScalingFactors

all: report

.PHONY: clean

report:
	pdflatex $(DOCNAME).tex
	bibtex8 $(DOCNAME).aux
	pdflatex $(DOCNAME).tex
	pdflatex $(DOCNAME).tex

view: report
	open $(DOCNAME).pdf

clean:
	rm *.blg *.bbl *.aux *.log *.out
