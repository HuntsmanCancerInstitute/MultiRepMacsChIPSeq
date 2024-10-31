# MultiRepMacsChIPSeq - generate\_mean\_bedGraph

|[Home](Readme.md)|[Overview](Overview.md)|[Usage](Usage.md)|[Variations](Variations.md)|[Applications](applications.md)|[Install](Install.md)|

## generate\_mean\_bedGraph.pl
  
A script to generate a chromosomal mean coverage bedGraph track 
to be used in Macs2 as the global control track when there is 
no input for generating a control_lambda chromatin bias track.
This uses a bigWig or bedGraph coverage file to calculate a global 
mean. Intervals without coverage are not included in the calculation.
It will write out a simple bedGraph representing the genome 
with the respective mean for each chromosome. 

Provide an input file in either bigWig or bedGraph format. 

USAGE:

	generate_mean_bedGraph.pl -i file.bw -o mean.bdg

OPTIONS:

    -i --in <file>      The input bigWig or bedGraph file
    -o --out <file>     The output bedGraph file (input.global_mean.bdg)
    -h --help           Display help
        


