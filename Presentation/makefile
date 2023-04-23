PROJECT = 'ScalingFactors_Presentation' 

LOG_File 		= $(shell mkdir -p LogFiles)

MAIN = ScalingFactors_Presentation

SECTIONS := $(shell find Inputs/* -iname '*.tex' && find Inputs/* -iname '*.cls')

TEX_SOURCES = 		Makefile \
              		$(MAIN).tex \
                	$(SECTIONS) 

FIGURES :=  $(shell find Figures/* -iname '*.{png,jpg,pdf}')

SHELL=/bin/bash

OPT = -interaction=nonstopmode 

$(MAIN).pdf: $(LOG_File) $(TEX_SOURCES) $(FIGURES)
	pdflatex $(MAIN).tex 
	pdflatex $(MAIN).tex $(MAIN).aux $(MAIN).toc
	@mv *.{aux,log,out,toc,nav,snm} LogFiles/
	clear
	@echo "				 											"
	@echo "			!!!!!!!!!!!!!!!!!!!!!!!!!!!			"
	@echo "			!!!!!Tex file compiled!!!!!			"
	@echo "			!!!!!!!!!!!!!!!!!!!!!!!!!!!			"
	@echo "															"
	@echo "															"
	@echo "			!!!!!!!!!!!!!!!!!!!!!!!!!!!			"
	@echo "			!!!!!!!!Check .pdf!!!!!!!!!			"
	@echo "			!!!!!!!!!!!!!!!!!!!!!!!!!!!			"
	@echo "															"

clean: 
	-rm -f *.{nav,snm,txt,tgz,zip,dvi,aux,fdb_latexmk,fls,idx,ilg,ind,lof,log,lot,out,pdf,toc} 
	-rm -rf LogFiles/
	-rm -f Inputs/*.{xdv,tgz,zip,dvi,aux,fdb_latexmk,fls,idx,ilg,ind,lof,log,lot,out,pdf,toc}

once: 
	pdflatex $(MAIN).tex

notes:
	pdflatex Notes.tex
	@mv *.{aux,log,out} LogFiles/
	clear
	@echo "				 											"
	@echo "			!!!!!!!!!!!!!!!!!!!!!!!!!!!			"
	@echo "			!!!!!Tex file compiled!!!!!			"
	@echo "			!!!!!!!!!!!!!!!!!!!!!!!!!!!			"
	@echo "															"
	@echo "															"
	@echo "			!!!!!!!!!!!!!!!!!!!!!!!!!!!			"
	@echo "			!!!!!!!!Check .pdf!!!!!!!!!			"
	@echo "			!!!!!!!!!!!!!!!!!!!!!!!!!!!			"
	@echo "															"

zip:
	$(MAKE) clean
	zip -q $(MAIN).zip $(TEX_SOURCES) $(FIGURES)

live: 
	@latexmk -pdf -pvc -f -pdflatex="pdflatex $(OPT)" $(MAIN).tex  
	@mv *.{aux,fdb_latexmk,fls,idx,ilg,ind,lof,log,lot,out,toc} LogFiles/
