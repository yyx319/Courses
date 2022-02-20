#SLIDES = DirectSolvers
#SLIDES = Pivoting
SLIDES = JacobiGS
slides:
	jupyter nbconvert --to markdown ${SLIDES}.ipynb
	pandoc -t beamer -N --slide-level=2 -o ${SLIDES}.pdf -H header.tex ${SLIDES}.md
clean:
	rm ${SLIDES}.md ${SLIDES}.pdf

${SLIDES}.md: 
	jupyter nbconvert --to markdown ${SLIDES}.ipynb

${SLIDES}.html: ${SLIDES}.md
	pandoc ${SLIDES}.md -t html -s -o ${SLIDES}.html --mathjax=https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML
