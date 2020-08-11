# MultiRepMacsChIPSeq Manual Installation

The MultiRepMacsChIPSeq package is primarily Perl scripts, but relies heavily on 
Python, R, and compiled C code as well. It assumes installation in a Unix-like environment 
(Linux and macOS) and runs via the command line. It may be installed in a custom 
working directory, suitable for multiple-user environments with tools such as 
[Modules](https://modules.readthedocs.io/en/latest) or in a Docker container.

Since this requires multiple programming language support, a single installation 
utility would be overly complex and likely redundant in most typical scientific 
computing environments. This document will try to guide most of package-specific 
installation stuff. Programming skills are not required, but familiarity with 
package installation and execution is assumed. 

**NOTE on existing applications:** This pipeline uses many common bioinformatic tools 
that may already be available on your compute system. If that is the case, you do not 
necessarily need to re-install. Simply point, or symlink, the tool into a convenient 
location, or otherwise ensure that it is in your `PATH`. 

**NOTE on installations:** The following examples show an installation into your 
home directory, noted by the environment variable `$HOME`. This may or may not be 
desirable. Adjust accordingly. 

## List of external library packages

This is a list of the primary Perl and R modules that need to be installed. These, 
of course, require various prerequisites, so the list is not to be considered 
absolute.

- Perl [Bio::ToolBox](https://metacpan.org/pod/Bio::ToolBox)
  
- Perl [Bio::DB::HTS](https://metacpan.org/pod/Bio::DB::HTS)

- Perl [Bio::DB::Big](https://metacpan.org/pod/Bio::DB::Big)

- Perl [Set::IntervalTree](https://metacpan.org/pod/Set::IntervalTree)

- Perl [Parallel::ForkManager](https://metacpan.org/pod/Parallel::ForkManager)

- Python [Macs2](http://github.com/taoliu/MACS/)

- R [ggplot2](https://ggplot2.tidyverse.org)

- R [pheatmap](https://cran.r-project.org/web/packages/pheatmap)

- R [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)

- R [normR](https://bioconductor.org/packages/release/bioc/html/normr.html)

- External [BedTools](https://github.com/arq5x/bedtools2). 

- External [UCSC bigWig utilities](http://hgdownload.soe.ucsc.edu/downloads.html#utilities_downloads)


## Language-specific installations

These are typically already installed on most research compute environments, either
with your OS or separately for other purposes. However, some modern OS do not 
include modern releases, so check your versions. OS package managers may have later 
versions available. MacOS users should consider [Homebrew](https://brew.sh). 
Similarly, Linux (and evidently Microsoft Windows Subsystem for Linux) users could 
use [Linuxbrew](http://linuxbrew.sh).

In all cases, the executable is assumed to be in the users' `PATH`.

### Perl

You can check your Perl version by running

	$ perl --version

A modern Perl version is at least `5.22`, with the latest being `5.30`. Some OS 
distributions include `5.16` or `5.18` as part of the base installation. If your 
version is that old, strongly consider installing a newer version. You can install 
Perl from [Perl.org](https://www.perl.org). An easy installation might be 

	$ sh ./Configure -d -e -s -Dprefix=$HOME
	$ make
	$ make install

This will accept standard configuration (64 bit, thread-enabled) and set the 
installation directory to your home directory, or wherever your specify.

### Python

You can check your Python version by running

	$ python3 --version

or simply

	$ python --version

A modern Python is at least `3.5`, with the latest being `3.8`. Python `2.7` or earlier 
is deprecated and you should seriously upgrade. You can install current versions of 
Python from [Python.org](https://www.python.org). 

### R

You can check your version of R by running

	$ R --version

A modern R version is at least `3.5`, with the latest being `4.0`. You can obtain the 
latest version at [R-project](https://www.r-project.org). 


## C libraries for necessary Perl modules

Install the [HTSlib](https://github.com/samtools/htslib) library for working with Bam 
files. Version `1.9` is recommended at this time. This may or may not already be 
present on your system. We need both the `samtools` executable and explicitly the 
`htslib` library.

	$ curl -o samtools-1.9.tar.bz2 \
      -L https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
    $ tar -xf samtools-1.9.tar.bz2 && cd samtools-1.9
    $ sh ./configure --prefix=$HOME && make all all-htslib && make install install-htslib

Install the [libBigWig](https://github.com/dpryan79/libBigWig) library for working with 
bigWig files. 

	$ curl -o libBigWig-0.4.4.tar.gz \
	  -L https://github.com/dpryan79/libBigWig/archive/0.4.4.tar.gz
	$ tar -xf libBigWig-0.4.4.tar.gz && cd libBigWig-0.4.4

The `Makefile` installs into `/usr/local` by default. Edit this file to change it to 
your target directory. For example, this will change to our `$HOME` directory.

	$ perl -i.bak -p -e "s|/usr/local|$HOME|" Makefile 

Then compile and install

	$ make && make install


## Perl Modules

I strongly recommend the use of a modern Perl package manager like [CPAN
Minus](https://metacpan.org/pod/App::cpanminus), used 
below. If you have installed your own Perl and have write permissions to the directory, 
then you do not need to set a `-l` prefix. Otherwise, if you have a custom working 
directory, use that. 

You may also consider installing [local::lib](https://metacpan.org/pod/local::lib) for 
automatically enabling your Perl environment; this is not necessary if you use some 
other means to export custom `PERL5LIB` and `PATH` environment variables.

The `cpanm` utility will conveniently install required prerequisites. There will be 
dozens of additional modules installed along the way. In some cases, we are explicitly 
installing specific versions to avoid complications.

    $ curl -L https://cpanmin.us | perl - -l $HOME local::lib App::cpanminus \
    && echo 'eval "$(perl -I$HOME/lib/perl5 -Mlocal::lib)"' >> ~/.profile \
    && . ~/.profile


Then install additional modules, including 

	$ cpanm Module::Build
	$ cpanm --notest Data::Swap Data::Types
	$ cpanm Parallel::ForkManager Set::IntervalTree Set::IntSpan::Fast::XS Template::Tiny Archive::Zip

Install specific version of BioPerl (we're trying to avoid installing unnecessary stuff 
here). We run the installer interactively; answer `n` to all of the interactive questions. 

	$ curl -o BioPerl-1.007002.tar.gz -L \
	  https://cpan.metacpan.org/authors/id/C/CJ/CJFIELDS/BioPerl-1.007002.tar.gz
	$ cpanm --interactive BioPerl-1.007002.tar.gz

We can then install additional modules

	$ curl -o Bio-DB-Big.zip -L \
	  https://github.com/Ensembl/Bio-DB-Big/archive/master.zip
	$ cpanm --notest --configure-args="--libbigwig $HOME" Bio::DB::Big.zip
	$ cpanm --configure-args="--htslib $HOME" Bio::DB::HTS
	$ curl -o Bio-ToolBox.zip -L \
	  https://github.com/tjparnell/biotoolbox/archive/master.zip
	$ cpanm Bio-ToolBox.zip
	$ curl -o MultiRepMacsChIPSeq.zip -L \
	  https://github.com/HuntsmanCancerInstitute/MultiRepMacsChIPSeq/archive/master.zip
	$ cpanm MultiRepMacsChIPSeq.zip


## Python Applications

The [Macs2](https://pypi.org/project/MACS2/) application can be installed from PyPi. It's 
best to install this into a separate environment. You can use a virtual environment in 
python3, placing it in its own directory, and then link the main executable to a 
convenient location if desired. The MultiRepMacsChIPSeq pipeline will need to call this 
executable directly without having to explicitly activate a virtual environment.

	$ python3 -m venv $HOME/macs2
	$ source $HOME/macs2/activate
	$ pip3 install macs2
	$ deactivate
	$ ln -s $HOME/macs2/bin/macs2 $HOME/bin


## R applications

You should be able to install the modules into your personal R library in your home 
directory, where they can be found when you execute an R script. Open an R terminal 
as below and execute the command to install the following packages. Note that some of 
them may already be installed if you've worked with R before.

	$ R
	> install.packages(c("optparse","RColorBrewer","reshape2","ggplot2","pheatmap"), repos="http://cran.r-project.org")
	> quit()
	$


## Additional C applications

Finally, there are additional helper applications that are needed. These may or may not 
be available already on your system; adjust accordingly.

First are helper bigWig [utilities](http://hgdownload.soe.ucsc.edu/downloads.html#utilities_downloads) 
from UCSC. The following command is for Linux.

    $ for name in wigToBigWig bedGraphToBigWig bigWigToWig bedToBigBed; \
    do curl -o $HOME/bin/$name -L http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/$name \
    && chmod +x $HOME/bin/$name; done;

Second is [BedTools](https://github.com/arq5x/bedtools2). 

	$ curl -o bedtools-2.29.2.tar.gz -L \
	  https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools-2.29.2.tar.gz
	$ tar -zxvf bedtools-2.29.2.tar.gz
	$ cd bedtools2
	$ make
	$ mv bedtools $HOME/bin


## Finished

You should now have a functional environment for MultiRepMacsChIPSeq.


