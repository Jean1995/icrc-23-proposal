texoptions = -lualatex \
			 -output-directory=build \
			 -interaction=nonstopmode \
			 -halt-on-error


all: build/poster.pdf

build/poster.pdf: poster.tex beamerthemetudoposter.sty code/example.txt | build
	latexmk $(texoptions) poster.tex

build:
	mkdir -p build

clean:
	rm -rf build
