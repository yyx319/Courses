#SLIDES = "Numbers and their Representation"
#SLIDES = "Python numbers"
#SLIDES = "Relative Error"
#SLIDES = "Decimal Rounding"
SLIDES = "Rounding"


#SLIDES = "Stability"

.IGNORE: clean

slides:
	
	jupyter nbconvert --to markdown ${SLIDES}.ipynb
	
	pandoc -t beamer -N --slide-level=2 -o ${SLIDES}.pdf -H header.tex ${SLIDES}.md

clean:
	@-rm --force ${SLIDES}.md ${SLIDES}.pdf
