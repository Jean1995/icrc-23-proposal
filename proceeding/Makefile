all: proposal_proceeding

proposal_proceeding: %:
	latexmk -pdf -outdir=build $@

clean:
	latexmk -C -outdir=build

.PHONY: proposal_proceeding clean
