# Pysano Alignment Templates

These are `cmd.txt` alignment templates for users of the HCI 
[pysano](https://healthcare.utah.edu/huntsmancancerinstitute/research/shared-resources/center-managed/bioinformatics/pysano/) 
system for executing jobs at [CHPC](https://www.chpc.utah.edu). If you are not a HCI 
user, feel free to use this as a shell script template for aligning and processing 
your reads.

## Using the command templates

- Modify the template as necessary

    Add your email address, cluster name if so desired, adjust the path to your 
    organism index file, and set the effective size of your organism genome. For 
    optical duplicate size checking, set the distance: 100 pixels for older HiSeq 
    runs, and 2500 for newer NovaSeq runs. This is customization common to all 
    of your sample alignment jobs.
    
    If necessary, tweak the commands as appropriate. In general, they probably 
    don't need much changing.

- Set the sample name

    For each of your sample alignment jobs, set the name of your sample. This can 
    be edited in batch using `sed`. For example, to edit the name of the sample and 
    write the resulting file into a job directory named after GNomEx identifier:
    
        $ sed 's/MYNAME/sample1/' cmd_template.txt > 1234X1/cmd.txt
    
    For lots of jobs, use Gnu `parallel` with a `list.txt` file, comprised of two 
    columns, GNomEx ID and sample name.
    
        $ parallel -a list.txt -C '\s' sed 's/MYNAME/{2}/' cmd_template.txt '>' {1}/cmd.txt

- Run your jobs

    Start your jobs using `pstart`.

- Collect metrics

    After all of your jobs are finished, you can use the included script to collate the 
    results
    
        $ combine_std_chipstats.pl 1234R_stats.txt 1234X1/ 1234X2/ 1234X3/ 1234X4/
    
    This script will descend into job directories, collect various statistic numbers 
    from the `stdout.txt` and `stderr.txt` files, collate them into a single file, and 
    write a tab-delimited text output file, suitable for loading into a spreadsheet 
    program or R and examining or plotting. It will collect statistics from Novoalign 
    alignment, `bam_partial_dedup`, `bam2wig`, and `Macs predictd`. Missing information 
    is recorded as empty columns.

- Plot Macs2 prediction PDFs

    Macs2 writes out its prediction as a R script, which can then be executed to 
    generate a 2-page PDF report of the profile and shift correlations. Look at the 
    second page to determine the optimum shift size, which is sometimes accurately 
    identified, but sometimes an alternate shift size is more suitable or desirable.
    
        $ Rscript 1234X1/sample1.predictd
    
- Plot bam2wig shift models

    The [bam2wig.pl](https://metacpan.org/pod/bam2wig.pl) script will write out two 
    text files: `$NAME_model.txt` shows the forward, reverse, and shifted mean coverage 
    at peaks, and `$NAME_correlations.txt` shows the mean correlations at different 
    shift values. These can be plotted with the included script.
    
        $ plot_shift_models.R --input 1234X1/sample1
    
    Provide the path and base output name (`$NAME`) that you indicated in your `cmd.txt`
    file. It will write out two corresponding PDF files. *NOTE*: the script requires 
    `ggplot2` module.


 
