# MultiRepMacsChIPSeq - Installation

|[Home](Readme.md)|[Overview](Overview.md)|[Usage](Usage.md)|[Variations](Variations.md)|[Examples](Examples.md)|[Applications](applications.md)|[Install](Install.md)|


## Installation

The MultiRepMacsChIPSeq package itself is primarily Perl scripts, but relies heavily
on Python, R, and compiled C code as well. It assumes installation in a Unix-like
environment (Linux and macOS) and runs via the command line. A standard build
environment (`make`, `gcc`, etc) is required and should be available from your OS
package manager. It may be installed in a custom package directory, suitable for
multiple-user environments with tools such as
[Modules](https://modules.readthedocs.io/en/latest), or in a container.

Since this requires multiple programming language support, a single installation
utility would be overly complex and likely redundant in most typical scientific
computing environments. This document will try to guide most of package-specific
installation stuff. Programming skills are not required, but familiarity with
package installation and execution is assumed. 

**NOTE on existing applications:** This pipeline uses many common bioinformatic tools
that may already be available on your compute system. If that is the case, you do not
necessarily need to re-install. Simply point, or symlink, the tool into a convenient
location, or otherwise ensure that it is in your `PATH`. 

**NOTE on installation PREFIX:** The following examples show an installation into a
custom directory, noted below by the environment variable `$DEST`. Adjust accordingly,
or leave out if installing to standard locations, i.e. `/usr/local`. 

## List of external library packages
This is a list of the primary Perl and R modules that need to be installed. These, 
of course, require various prerequisites, so the list is not to be considered
absolute.

- Perl [Bio::ToolBox](http://tjparnell.github.io/biotoolbox)
  
- Perl [Bio::DB::HTS](https://metacpan.org/pod/Bio::DB::HTS)

- Perl [Bio::DB::Big](https://metacpan.org/pod/Bio::DB::Big)

- Python [Macs2](http://github.com/taoliu/MACS/)

- R [ggplot2](https://ggplot2.tidyverse.org)

- R [pheatmap](https://cran.r-project.org/web/packages/pheatmap)

- R [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)

- R [normR](https://bioconductor.org/packages/release/bioc/html/normr.html)

- External [Samtools](https://github.com/samtools/samtools)

- External [BedTools](https://github.com/arq5x/bedtools2). 

- External [UCSC bigWig utilities](http://hgdownload.soe.ucsc.edu/downloads.html#utilities_downloads)

- External [Pandoc](https://pandoc.org)


## Language-specific installations

These are typically already installed on most research compute environments, either
with your OS or separately for other purposes. However, some modern OS do not
include modern releases, so check your versions. OS package managers may have later
versions available. MacOS users should consider [Homebrew](https://brew.sh). 
Similarly, Linux (and evidently Microsoft Windows Subsystem for Linux) users could
use [Linuxbrew](http://linuxbrew.sh).

In all cases, the executable is assumed to be in the user's `PATH`.

### Perl

Perl should be available on most Linux and macOS OS versions and distributions.
You can check your Perl version by running

	perl --version

A modern Perl should be at least `5.3x`, with the latest being `5.40` as of this
writing. If your version is < `5.30`, strongly consider installing a newer version. You
can install Perl from [Perl.org](https://www.perl.org). If you're managing multiple
Perl versions, look at a Perl manager such as [PerlBrew](https://perlbrew.pl).


### Python

Python should be available on most Linux and macOS versions and distributions.
You can check your Python version by running

	python3 --version

or simply

	python --version

A modern Python is at least `3.8` and should be available in your OS package manager.
You can install current versions of Python from [Python.org](https://www.python.org).

### R

You can check your version of R by running

	R --version

A modern R version is at least `4.x`. You can obtain the latest version at
[R-project](https://www.r-project.org). 


## C libraries for necessary Perl modules

Install the [HTSlib](https://github.com/samtools/htslib) library for working with Bam
files. This may or may not already be present on your system. It may also be
available through your package manger. Version `1.9` is officially recommended by the
Bio::DB::HTS authors; however, later versions appear to work just fine and should
probably be preferred. [Version
1.19](https://github.com/samtools/htslib/releases/download/1.19/htslib-1.19.tar.bz2)
has been used successfully by this author. 

Note that the `htslib` package is included with the `samtools` source code packages
and is compiled therein, but is not installed by default. The code below compiles and
installs both `samtools` and the `htslib` library (since we need `samtools` as well).

	curl -o samtools-1.19.2.tar.bz2 \
	  -L https://github.com/samtools/htslib/releases/download/1.19/samtools-1.19.2.tar.bz2
	tar -xf samtools-1.19.2.tar.bz2
	cd samtools-1.19.2
	sh ./configure --prefix=$DEST
	make all all-htslib
	make install install-htslib

Install the [libBigWig](https://github.com/dpryan79/libBigWig) library for working with
bigWig files. 

	curl -o libBigWig-0.4.7.tar.gz -L https://github.com/dpryan79/libBigWig/archive/0.4.7.tar.gz
	tar -xf libBigWig-0.4.7.tar.gz
	cd libBigWig-0.4.7

The `Makefile` installs into `/usr/local` by default. Leave as is or edit this file to
change it to your target directory. For example,

	perl -i.bak -p -e "s|/usr/local|$DEST|" Makefile

Then compile and install

	make && make install

## Perl Modules

A modern Perl package manager, such as [CPAN
Minus](https://metacpan.org/pod/App::cpanminus), or any other package manager, is
highly recommended and used in the examples below. Follow the link to install it. The
`cpanm` utility will conveniently install required prerequisites with little
additional effort on the user. Most Perl packages require numerous dependencies.

Where the Perl modules will be installed depends on how Perl is installed and whether
you have write access.

- Private Perl installation

	You should be able to write directly under the Perl installation; no adjustments
	should need to be made.

- System Perl with specific destination directory

	You can use `cpanm -l $DEST` to target the files into your destination prefix, by
	default `$DEST/lib/perl5` for modules. You will need to manually set your
	environment variable `PERL5LIB` to this directory and append `$DEST/bin` to your
	`$PATH` environment variables as appropriate.

- System Perl with `local::lib`

	The [local::lib](https://metacpan.org/pod/local::lib) is a convenience Perl
	module that automatically sets your private Perl environment, usually in your
	`$HOME` folder. An example installation and modifying your shell `.profile` to
	modify your environment would be

		curl -L https://cpanmin.us | perl - -l $DEST local::lib App::cpanminus
		echo 'eval "$(perl -I$HOME/lib/perl5 -Mlocal::lib)"' >> ~/.profile
		source ~/.profile

- System Perl with root access

	If you have root privileges and can write to the system Perl libraries, you may
	use `cpanm --sudo` during installation to elevate privileges. This is generally
	not recommended except possibly personal workstations. If building a container,
	you are likely already `root` and can install with impunity.
	
Once location is determined, you can begin installing Perl modules. If you're
installing on a macOS system and need additional help, take a look at the
[MacOS Notes](http://tjparnell.github.io/biotoolbox/MacOSNotes.html) for BioToolBox.

[BioPerl](https://metacpan.org/pod/BioPerl) is a large bundle of bioinformatics
related Perl modules and is required by Bio::DB::HTS. Unless you have other software
needs for BioPerl, install an unofficial, custom, stripped-down [minimal
version](https://github.com/tjparnell/bioperl-live/tree/minimal-tjparnell), which is
considerably faster and easier to install and has a footprint one-third the size of
the full distribution.

	curl -O -L https://github.com/tjparnell/bioperl-live/releases/download/minimal-v1.7.8/Minimal-BioPerl-1.7.8.tar.gz
	cpanm Minimal-BioPerl-1.7.8.tar.gz

The [Bio::DB::HTS](https://metacpan.org/pod/Bio::DB::HTS) package requires the C
`htslib` library that was compiled and installed above. If this was installed in a
standard prefix location, such as `/usr/local`, then it should be found
automatically. Otherwise you need to specify the location with an option to the
`Build.PL` installer, which you can specify in `cpanm` as below

	cpanm --configure-args="--htslib $DEST" Bio::DB::HTS

Likewise, the [Bio::DB::Big](https://metacpan.org/pod/Bio::DB::Big) needs the
location of the `libBigWig` C library installed above, particularly if you've
installed into a custom path. The version on CPAN is currently out-of-date, so
install from the GitHub repository. It also requires `Template::Tiny`, which was
(mistakenly?) not specified in the original manifest.

	cpanm Template::Tiny
	curl -o bio-db-big-master.tar.gz -L https://github.com/Ensembl/Bio-DB-Big/tarball/master
	cpanm --configure-args="--libbigwig $DEST" bio-db-big-master.tar.gz

The last package is the
[MultiRepMacsChIPSeq](https://github.com/HuntsmanCancerInstitute/MultiRepMacsChIPSeq)
package itself. This should install a number of other prerequisites, including
[Bio::ToolBox](https://metacpan.org/pod/Bio::ToolBox). The code below extracts the
tar ball and builds the package in situ. This source directory is required to run
example pipelines and test the installation; see below and the [Examples
page](Examples.md).

	curl -o multirepchipseq-master.tar.gz -L https://github.com/HuntsmanCancerInstitute/MultiRepMacsChIPSeq/tarball/master
	tar -xf multirepchipseq-master.tar.gz
	cd HuntsmanCancerInstitute-MultiRepMacsChIPSeq-*
	cpanm .



## Python Applications

The [MACS2](https://pypi.org/project/MACS2/) application can be installed from PyPi.
It's best to install this into a separate environment. You can use a virtual
environment in python3, placing it in its own directory, and then link the main
executable to a convenient location if desired. The MultiRepMacsChIPSeq pipeline will
need to call this executable directly without having to explicitly activate a virtual
environment. 

	python3 -m venv $DEST/macs2
	source $DEST/macs2/activate
	pip3 install MACS2
	deactivate
	ln -s $DEST/macs2/bin/macs2 $DEST/bin

## R applications

You should be able to install the modules into your personal R library in your home
directory, where they can be found when you execute an R script. Open an R terminal
by executing `R` and run the following commands to install the required packages.
Note that some of them may already be installed if you've worked with R before.

	install.packages(c("optparse","RColorBrewer","reshape2","ggplot2","pheatmap"), repos="http://cran.r-project.org")
	quit()

These R packages are the essentials for generating plots. Additional packages may be
installed as desired, for example
[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) from
BioConductor.


## Additional applications

Finally, there are additional helper applications that are needed. These may or may not
be available already on your system; adjust accordingly. They should be placed in the
working environment `PATH` as appropriate.

The first is [Samtools](https://github.com/samtools/samtools). This may already be
installed from steps above.

The second are helper [bigWig
utilities](http://hgdownload.soe.ucsc.edu/downloads.html#utilities_downloads) 
from UCSC. The following command collects Linux executables and copies them
into your destination `bin` directory. MacOS executables are also available.

	for name in wigToBigWig bigWigToBedGraph bigWigToWig bedToBigBed; \
	do curl -o $DEST/bin/$name -L http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/$name \
	&& chmod +x $DEST/bin/$name; done;

The third is [BedTools](https://github.com/arq5x/bedtools2), which you may already
have on your system. It may be available through your OS package manager; if not, you
can build from source as below.

	curl -o bedtools-2.29.2.tar.gz -L \
	  https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools-2.29.2.tar.gz
	tar -zxvf bedtools-2.29.2.tar.gz
	cd bedtools2
	make
	mv bedtools $DEST/bin

Last is [Pandoc](https://pandoc.org), which is used to generate a self-contained HTML
document from the Markdown report and PNG images generated by the pipeline. This is
usually available from your OS package manager. While not critical for the execution
of the pipeline, it is nevertheless convenient for interpreting results.


## Finished

You should now have a functional environment for MultiRepMacsChIPSeq. 

To test whether all helper applications are available and found by the pipeline,
run the following command

	multirep_macs2_pipeline.pl --help

In the output under the section `Application Paths`, the full path is printed
for each helper application found in your environment `PATH`. If the path is
empty, the application was not found.



## Tests

There are no official tests for the MultiRepMacsChIPSeq package, primarily because of
so many external dependencies. However, there are example pipeline scripts and dummy
data files located in the
[examples directory](https://github.com/HuntsmanCancerInstitute/MultiRepMacsChIPSeq/blob/master/examples).

To run these tests, change into the `MultiRepMacsChIPSeq` source folder (see above)
and verify that the `blib` directory exists. If not, then the package will have to be
(re)built. Run the example scripts, for example

	cd examples
	./run_se.sh

These example pipelines should only take 1-3 minutes to complete. If successful, the
standard output should say `Finished ChIPSeq multi-replicate pipeline`. The output
files will be in a subdirectory, and include multiple sub directories, a
`<NAME>_job_output_logs.txt` log file, a `<NAME>_summary.tsv` table file, a
`<NAME>.md` report markdown file, and if `pandoc` is installed, a `<NAME>.html`
report file.

Look at the [Examples page](Examples.md) for additional details of the different
example scripts.





