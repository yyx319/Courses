#SLIDES = DirectSolvers
SLIDES = Pivoting
#SLIDES = JacobiGS
#SLIDES = Gradients
#SLIDES = Bisection
#SLIDES = FixedPoint
#SLIDES = Newton
#SLIDES = SOR_IR
#SLIDES = FFT

slides: ${SLIDES}.md
	pandoc +RTS -V0 -RTS -t beamer -N -fmarkdown-implicit_figures --slide-level=2 -o ${SLIDES}.pdf -H header.tex ${SLIDES}.md
clean:
	-rm -f ${SLIDES}.md ${SLIDES}.pdf

${SLIDES}.md: ${SLIDES}.ipynb
	jupyter nbconvert --to markdown ${SLIDES}.ipynb

html: ${SLIDES}.md
	pandoc +RTS -V0 -RTS ${SLIDES}.md -t html -s -o ${SLIDES}.html --mathjax=https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML

tex: ${SLIDES}.md
	pandoc +RTS -V0 -RTS -t beamer -N --slide-level=2 -s -o ${SLIDES}.tex -H header.tex ${SLIDES}.md

notebook:
	jupyter notebook --no-browser
