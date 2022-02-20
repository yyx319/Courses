#SLIDES = ODEs
#SLIDES = OneStep
#SLIDES = Runge-Kutta
SLIDES = A-Stability

pdfslides: ${SLIDES}.md
	pandoc +RTS -V0 -RTS -t beamer -N -fmarkdown-implicit_figures --slide-level=2 -o ${SLIDES}.pdf -H header.tex ${SLIDES}.md

${SLIDES}.md: ${SLIDES}.ipynb
	jupyter nbconvert --to markdown $<

clean:
	-rm -f ${SLIDES}.md ${SLIDES}.html ${SLIDES}.slides.html

reveal: ${SLIDES}.md
	pandoc -t revealjs $< -o ${SLIDES}.html -S -s --mathjax=https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML --slide-level=2 -c white.css

nbslides: ${SLIDES}.ipynb
	jupyter nbconvert --to slides ${SLIDES}.ipynb

html: ${SLIDES}.md
	pandoc +RTS -V0 -RTS ${SLIDES}.md -t html -s -o ${SLIDES}.html --mathjax=https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML

tex: ${SLIDES}.md
	pandoc +RTS -V0 -RTS -t beamer -N --slide-level=2 -s -o ${SLIDES}.tex -H header.tex ${SLIDES}.md

notebook:
	jupyter notebook --no-browser
