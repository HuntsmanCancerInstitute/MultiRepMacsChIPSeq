# MultiRepMacsChIPSeq - render

|[Home](Readme.md)|[Overview](Overview.md)|[Usage](Usage.md)|[Variations](Variations.md)|[Applications](applications.md)|[Install](Install.md)|

## render.pl

A script to render a markdown report into a self-contained HTML file.
The HTML file will be written in the same directory as the input file.
This uses [Pandoc](https://pandoc.org), which must be in the environment
`PATH` or explicitly provided.

VERSION: 0.1

USAGE:

	render.pl report.md 

OPTIONS:

	-i --in <file>           The markdown report
	--pandoc <path>          Path to pandoc ($app)
	-h --help                Print documentation


