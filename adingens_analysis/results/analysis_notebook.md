# Mutational antigenic profiling of broad human IgG
### Adam S. Dingens, Jesse Bloom
### In collabortation with Florian Klein's group, particularly Philipp Schommers and Henning Gruell.

We are performing mutational antigenic profiling of IgG purified from individuals in Florian Klein's cohort, using the BG505.T332N mutant Env libraries, first described and characterized in [Haddox, Dingens et al 2018](https://elifesciences.org/articles/34420). The triplicate mutant libraries examined here correspond to the three BG505 replicates in this paper. Replicates annotated with an additional letter (i.e. "rep 1b") were done on an indepdent day and have their own mock selected control. 

The _fraction surviving_ statistic used in this analysis is explained in detail [here](https://jbloomlab.github.io/dms_tools2/fracsurvive.html). Our first mutational antigenic profiling analysis of escape from PGT151 using the BF520 env libraries was published [here](http://dx.doi.org/10.1016/j.chom.2017.05.003) in June 2017, and this original analysis is located [in this ipython notebook](https://github.com/adingens/BF520_MutationalAntigenicProfiling_PGT151).

We use [dms_tools2](https://jbloomlab.github.io/dms_tools2/) to analyze these data. This notebook processes the Illumina deep sequencing data software package, and then analyzes the selection in the context of the antibody. Experiments and analysis performed by Adam Dingens in the [Bloom lab](http://research.fhcrc.org/bloom/en.html) in early 2019. 

## Import `Python` modules, define widely used functions
Set `use_existing` to `yes` below if you want to use any existing output, and `no` if you want to re-generate all output.

### Import Python modules / packages
Import modules / packages.
In particular, we use:

 - [plotnine](https://plotnine.readthedocs.io) for ggplot2-like plotting syntax
 - [dmslogo](https://jbloomlab.github.io/dmslogo/) to draw sequence logo plots
 - [dms_tools2](https://jbloomlab.github.io/dms_tools2/) for much of the analysis


```python
import os
import math
import collections
import itertools
import warnings
import subprocess

import yaml
import pandas as pd

from IPython.display import display, HTML

import matplotlib
matplotlib.use('webagg')

import matplotlib.pyplot as plt

from plotnine import *

import dms_tools2
from dms_tools2.ipython_utils import showPDF
from dms_tools2.plot import COLOR_BLIND_PALETTE_GRAY as PALETTE

import dmslogo
%matplotlib inline


print("Using dms_tools2 version {0}".format(dms_tools2.__version__))

# results will go in this directory
resultsdir = './results/' 
if not os.path.isdir(resultsdir):
    os.mkdir(resultsdir)
    
# CPUs to use, -1 means all available
ncpus = 14

# do we use existing results or generate everything new?
use_existing = "yes"
```

    Using dms_tools2 version 2.6.6


    /fh/fast/bloom_j/software/conda_v2/envs/BloomLab/lib/python3.6/site-packages/dmslogo/logo.py:40: MatplotlibDeprecationWarning: 
    The createFontList function was deprecated in Matplotlib 3.2 and will be removed two minor releases later. Use FontManager.addfont instead.
      matplotlib.font_manager.findSystemFonts(_FONT_PATH)))


## UPDATE UPON PUB - Download the sequencing data from the Sequence Read Archive
Here we download sequencing data for each sample from the [Sequence Read Archive](http://www.ncbi.nlm.nih.gov/sra). Sample names specificy if they are mutant (mut) or wildtype (wt), as well as virus (virus) or DNA (DNA). The _rep#_ is which replicate.  If a sample is antibody selected, it also has the antibody and `ug/mL` concentration. There is metadata available from the SRA, and the dictionary below also specifies the accession number of each sample.

The 1-18 reads were submitted as SRA # [PRJNA561627](https://www.ncbi.nlm.nih.gov/sra/PRJNA561627) with SRA accession numbers `SRX6752366-SRX6752371`  on Aug 22, 2019. The BioProject ID is `PRJNA561627`. 

We will also re-analyze data on VRC01 and 3BNC117 from [Dingens et al Immunity 2019](https://www.cell.com/immunity/fulltext/S1074-7613(18)30565-X). These reads were submitted as SRA submission [SRP157948](https://www.ncbi.nlm.nih.gov/sra/SRP157948) on Aug 14 and Aug 27, 2018. The BioProject ID is `PRJNA486029`. 

We download these files using the [dms_tools2.sra.fastqFromSRA](https://jbloomlab.github.io/dms_tools2/dms_tools2.sra.html) function from the dms_tools2 Python API. Note that the call to this function below uses two external programs that are not part of dms_tools2, and which you therefore must install externally on the computer that you are using:
1. The fastq-dump program from the SRA Toolkit. If you do not already have this toolkit installed, you will need to install a relatively recent version.
2. The Aspera Connect program for rapid downloads. You need both the executable ascp and an Aspera Key file. Installing Aspera Connect and a key can be somewhat complex, so if you do not want to do this then just set aspera=None in the command below and fastq-dump will do the downloads (albeit more slowly).



samples = pd.DataFrame.from_records(
#new data published here for 1-18
        [('BG505_mut_virus_rep3d_1-18-4ug','SRR10014244'),
         ('BG505_mut_virus_rep3d_1-18-8ug','SRR10014243'),
         ('BG505_mut_virus_rep2d_1-18-4ug','SRR10014242'),
         ('BG505_mut_virus_rep2d_1-18-8ug','SRR10014241'),
         ('BG505_mut_virus_rep2d','SRR10014240'),
         ('BG505_mut_virus_rep3d','SRR10014239'),
#Data on VRC01 and 3BNC117 from Dingens et al Immunity 2019 
#Here, I do not download the data on 10-1074 and pooled 10-1074/3BNC117. While I look at this data briefly for one analysis, I simply download the analyzed files from github rather than redoing all analyses. However these can be downloaded and anaylzed in parallel by uncommentin the relevant lines below. 
         ('BG505_mut-virus-rep1b','SRR7693968'),
         ('BG505_mut-virus-rep2b-3BNC117-4ug','SRR7693969'),
#         ('BG505_wt-virus-rep2b','SRR7693970'),
         ('BG505_mut-virus-rep1-VRC01-11ug','SRR7693971'),
#         ('BG505_wt-virus-rep3','SRR7693972'),
#         ('BG505_mut-virus-rep2b-101074-3ug','SRR7693974'),
         ('BG505_mut-virus-rep3','SRR7693976'),
         ('BG505_mut-virus-rep3-3BNC117-1-1ug','SRR7693977'),
#         ('BG505_mut-virus-rep2b-3BN-1074-pool-4ug','SRR7693978'),
#         ('BG505_mut-virus-rep1b-3BN-1074-pool-3ug','SRR7693979'),
#         ('BG505_wt-virus-rep1b','SRR7693980'),
#         ('BG505_mut-virus-rep1b-3BN-1074-pool-4ug','SRR7693981'),
#         ('BG505_mut-virus-rep3b-3BN-1074-pool-4ug','SRR7693982'),
#         ('BG505_wt-virus-rep1','SRR7693983'),
#         ('BG505_wt-virus-rep3b','SRR7693984'),
#         ('BG505_mut-virus-rep2b-101074-4ug','SRR7693985'),
         ('BG505_mut-DNA-rep1','SRR7693986'),
         ('BG505_mut-DNA-rep3','SRR7694021'),
         ('BG505_mut-virus-rep2b','SRR7694018'),
         ('BG505_wt-DNA-rep3','SRR7694017'),
         ('BG505_mut-virus-rep2-VRC01-8ug','SRR7694015'),
#         ('BG505_mut-virus-rep2b-3BN-1074-pool-3ug','SRR7694014'),
#         ('BG505_wt-virus-rep2','SRR7694013'),
#         ('BG505_mut-virus-rep3b-101074-4ug','SRR7694011'),
#         ('BG505_mut-virus-rep3b-3BN-1074-pool-3ug','SRR7694009'),
         ('BG505_mut-virus-rep1b-3BNC117-4ug','SRR7694006'),
         ('BG505_mut-virus-rep3b-3BNC117-4ug','SRR7694005'),
#         ('BG505_mut-virus-rep1b-101074-4ug','SRR7694004'),
         ('BG505_mut-virus-rep3b','SRR7694003'),
         ('BG505_mut-DNA-rep2','SRR7694002'),
         ('BG505_wt-DNA-rep2','SRR7693998'),
         ('BG505_mut-virus-rep2','SRR7693997'),
         ('BG505_mut-virus-rep3-VRC01-tr2-11ug','SRR7693996'),
         ('BG505_mut-virus-rep3-VRC01-11ug','SRR7693992'),
         ('BG505_mut-virus-rep2-VRC01-11ug','SRR7693991'),
         ('BG505_mut-virus-rep2b-3BNC117-3ug','SRR7693990'),
         ('BG505_wt-DNA-rep1','SRR7693989'),
         ('BG505_mut-virus-rep1','SRR7693987')],
        columns=['name', 'run']
        )


fastqdir = './results/FASTQ_files/'
if not os.path.isdir(fastqdir):
    os.mkdir(fastqdir)
print("Downloading FASTQ files from the SRA...")
dms_tools2.sra.fastqFromSRA(
        samples=samples,
        fastq_dump='fastq-dump', # valid path to this program on the Hutch server
        fastqdir=fastqdir,
        aspera=(
            '/app/aspera-connect/3.5.1/bin/ascp', # valid path to ascp on Hutch server
            '/app/aspera-connect/3.5.1/etc/asperaweb_id_dsa.openssh' # Aspera key on Hutch server
            ),
        )
print("Here are the names of the downloaded files now found in {0}".format(fastqdir))
display(HTML(samples.to_html(index=False)))

## Define samples from FASTQ_files


```python
R1fastqfilelist_df = pd.read_csv("./data/samples.csv", header =0)
display(HTML(R1fastqfilelist_df.to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>name</th>
      <th>R1</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>mtvir-rep3d</td>
      <td>../../2019/MAP_118/results/FASTQ_files/BG505_mut_virus_rep3d_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep2d</td>
      <td>../../2019/MAP_118/results/FASTQ_files/BG505_mut_virus_rep2d_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep2d-118-4ug</td>
      <td>../../2019/MAP_118/results/FASTQ_files/BG505_mut_virus_rep2d_1-18-4ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep2d-118-8ug</td>
      <td>../../2019/MAP_118/results/FASTQ_files/BG505_mut_virus_rep2d_1-18-8ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3d-118-4ug</td>
      <td>../../2019/MAP_118/results/FASTQ_files/BG505_mut_virus_rep3d_1-18-4ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3d-118-8ug</td>
      <td>../../2019/MAP_118/results/FASTQ_files/BG505_mut_virus_rep3d_1-18-8ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>wt-DNA-rep2</td>
      <td>../../2019/MAP_118/results/FASTQ_files/BG505_wt-DNA-rep2_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>wt-DNA-rep3</td>
      <td>../../2019/MAP_118/results/FASTQ_files/BG505_wt-DNA-rep3_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>wt-DNA-rep1</td>
      <td>../../2019/MAP_118/results/FASTQ_files/BG505_wt-DNA-rep1_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-IDC0136-5500ug</td>
      <td>./FASTQ_files/BG-M12b-IDC0136-5500ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-IDC0136-7000ug</td>
      <td>./FASTQ_files/BG-M12b-IDC0136-7000ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-IDC0136-8500ug</td>
      <td>./FASTQ_files/BG-M12b-IDC0136-8500ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-IDC0337-4000ug</td>
      <td>./FASTQ_files/BG-M12b-IDC0337-4000ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-IDC0337-4500ug</td>
      <td>./FASTQ_files/BG-M12b-IDC0337-4500ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-IDC403-7500ug</td>
      <td>./FASTQ_files/BG-M12b-IDC403-7500ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-IDC403-6000ug</td>
      <td>./FASTQ_files/BG-M12b-IDC403-6000ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-IDC208-6000ug</td>
      <td>./FASTQ_files/BG-M12b-IDC208-6000ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-IDC208-5000ug</td>
      <td>./FASTQ_files/BG-M12b-IDC208-5000ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-IDC208-4000ug</td>
      <td>./FASTQ_files/BG-M12b-IDC208-4000ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-IDC0337-6000ug</td>
      <td>./FASTQ_files/BG-M12b-IDC0337-6000ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b</td>
      <td>./FASTQ_files/BG-M12b-mock_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-561-IgG-2800ug</td>
      <td>./FASTQ_files/BG-M12b-561-IgG-2800ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-561-IgG-2200ug</td>
      <td>./FASTQ_files/BG-M12b-561-IgG-2200ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-212-6ug</td>
      <td>./FASTQ_files/BG-M12b-2-12-6ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-QA013-2282dpi-d40</td>
      <td>./FASTQ_files/BG-M12b-QA013-2282dpi-d40_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-QA013-2282dpi-d25</td>
      <td>./FASTQ_files/BG-M12b-QA013-2282dpi-d25_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-IDC403-8500ug</td>
      <td>./FASTQ_files/BG-M12b-IDC403-8500ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3e-212-8ug</td>
      <td>../../2019/CD4bs_bnAb_comparison/FASTQ_files/BG505_BG5-M3E-2-12-8ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3e</td>
      <td>../../2019/CD4bs_bnAb_comparison/FASTQ_files/BG505_BG5-M3E-mock_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3e-QA013-2282dpi-d7-5</td>
      <td>../../2019/CD4bs_bnAb_comparison/FASTQ_files/BG505_BG5-M3E-QA013-2282dpi-1-7-5_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3e-561-IgG-2500ug</td>
      <td>../../2019/CD4bs_bnAb_comparison/FASTQ_files/BG505_BG5-M3E-561-IgG-2500ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2a-QA013-2282dpi-d15</td>
      <td>../../2019/CD4bs_bnAb_comparison/FASTQ_files/BG505_BG5-M1v2A-QA013-2282dpi-1-15_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2a-QA013-2282dpi-d10</td>
      <td>../../2019/CD4bs_bnAb_comparison/FASTQ_files/BG505_BG5-M1v2A-QA013-2282dpi-1-10_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2a</td>
      <td>../../2019/CD4bs_bnAb_comparison/FASTQ_files/BG505_BG5-M1v2A-mock_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2a-561-IgG-3500ug</td>
      <td>../../2019/CD4bs_bnAb_comparison/FASTQ_files/BG505_BG5-M1v2A-561-IgG-3500ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2a-561-IgG-3000ug</td>
      <td>../../2019/CD4bs_bnAb_comparison/FASTQ_files/BG505_BG5-M1v2A-561-IgG-3000ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2a-212-10ug</td>
      <td>../../2019/CD4bs_bnAb_comparison/FASTQ_files/BG505_BG5-M1v2A-2-12-10ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2a-212-7ug</td>
      <td>../../2019/CD4bs_bnAb_comparison/FASTQ_files/BG505_BG5-M1v2A-2-12-7ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep2-QA013-2-2ug</td>
      <td>../../2017/BG505_MAP/FASTQ_files/BG505_mut_virus_rep2_QA013-2_2ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep2-QA013-2-4ug</td>
      <td>../../2017/BG505_MAP/FASTQ_files/BG505_mut_virus_rep2_QA013-2_4ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3-QA013-2-4ug</td>
      <td>../../2017/BG505_MAP/FASTQ_files/BG505_mut_virus_rep3_QA013-2_4ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep2</td>
      <td>../../2017/BG505_MAP/FASTQ_files/BG505_mut_virus_rep2_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3</td>
      <td>../../2017/BG505_MAP/FASTQ_files/BG505_mut_virus_rep3_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3v2a</td>
      <td>./FASTQ_files/BG-M32a-mock_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3v2a-IDC0191-6000ug</td>
      <td>./FASTQ_files/BG-M32a-IDC0191-6000ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3v2a-IDC0245-8000ug</td>
      <td>./FASTQ_files/BG-M32a-IDC0245-8000ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3v2a-IDC0396-9000ug</td>
      <td>./FASTQ_files/BG-M32a-IDC0396-9000ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3v2a-IDC0513-5000ug</td>
      <td>./FASTQ_files/BG-M32a-IDC0513-5000ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3v2a-IDC0518-8000ug</td>
      <td>./FASTQ_files/BG-M32a-IDC0518-8000ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3v2a-IDC0546-8500ug</td>
      <td>./FASTQ_files/BG-M32a-IDC0546-8500ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3v2a-IDF033-1400ug</td>
      <td>./FASTQ_files/BG-M32a-IDF033-1400ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3v2a-IDF033-1800ug</td>
      <td>./FASTQ_files/BG-M32a-IDF033-1800ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3v2a-IDF065-1400ug</td>
      <td>./FASTQ_files/BG-M32a-IDF065-1400ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3v2a-IDF141-6000ug</td>
      <td>./FASTQ_files/BG-M32a-IDF141-6000ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3v2a-IDF147-5300ug</td>
      <td>./FASTQ_files/BG-M32a-IDF147-5300ug_R1.fastq.gz</td>
    </tr>
  </tbody>
</table>


## Process the FASTQ files to count the mutations for each sample
We used [barcoded-subamplicon](https://jbloomlab.github.io/dms_tools2/bcsubamp.html) sequencing to obtain high accuracy during the Illumina deep sequencing. We therefore use the [dms2_batch_bcsubamp](https://jbloomlab.github.io/dms_tools2/dms2_batch_bcsubamp.html#dms2-batch-bcsubamp) program to analyze these data.

Running that program requires specifying a `--batchfile` that lists the samples, a wildtype `--refseq` to which we make alignments, and --alignspecs that tell us where the subamplicons should align. 
The batch file that we specify is printed by the cell below. The alignment specs need to be exactly correct for the subamplicons to align. We also do some trimming of the reads using the `--R1trim` and `--R2trim` parameters. 

The wildtype sequence of the BG505.W6.C2.T332N Env used in this experiment is in the file [./data/BG505.W6.C2.T332N_env.fasta](./data/BG505.W6.C2.T332N_env.fasta).
This sequence is based on GenBank accession number [DQ208458.1](https://www.ncbi.nlm.nih.gov/nucleotide/77025198?report=genbank&log$=nuclalign&blast_rank=1&RID=WMZ5XNUG014), with the introduced T332N mutation to knock in the glycan commonly targeted by this class of bnAbs. 


```python
# counts and alignments placed in this directory
countsdir = os.path.join(resultsdir, 'codoncounts')
if not os.path.isdir(countsdir):
    os.mkdir(countsdir)
```


```python
refseq = './data/BG505.W6.C2.T332N_env.fasta'

#fastqdir ='./FASTQ_files/'
# define subamplicon alignment specifications
alignspecs = ' '.join(['87,375,39,36', 
                       '376,666,36,39',
                       '663,954,33,41',
                       '955,1228,33,37',
                       '1228,1527,34,35',
                       '1527,1815,32,39',
                       '1816,2098,36,41'])


    
# write sample information to a batch file for dms2_batch_bcsubamplicons
countsbatchfile = os.path.join(countsdir, 'batch.csv')
print("Here is the batch file that we write to CSV format to use as input:")
display(HTML(R1fastqfilelist_df[['name', 'R1']].to_html(index=False)))
R1fastqfilelist_df[['name', 'R1']].to_csv(countsbatchfile, index=False)

print('\nNow running dms2_batch_bcsubamp...')
log = !dms2_batch_bcsubamp \
        --batchfile {countsbatchfile} \
        --refseq {refseq} \
        --alignspecs {alignspecs} \
        --outdir {countsdir} \
        --summaryprefix summary \
        --R1trim 200 \
        --R2trim 170 \
        --minq 15 \
        --ncpus {ncpus} \
        --use_existing {use_existing} 
print("Completed dms2_batch_bcsubamp.")
```

    Here is the batch file that we write to CSV format to use as input:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>name</th>
      <th>R1</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>mtvir-rep3d</td>
      <td>../../2019/MAP_118/results/FASTQ_files/BG505_mut_virus_rep3d_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep2d</td>
      <td>../../2019/MAP_118/results/FASTQ_files/BG505_mut_virus_rep2d_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep2d-118-4ug</td>
      <td>../../2019/MAP_118/results/FASTQ_files/BG505_mut_virus_rep2d_1-18-4ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep2d-118-8ug</td>
      <td>../../2019/MAP_118/results/FASTQ_files/BG505_mut_virus_rep2d_1-18-8ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3d-118-4ug</td>
      <td>../../2019/MAP_118/results/FASTQ_files/BG505_mut_virus_rep3d_1-18-4ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3d-118-8ug</td>
      <td>../../2019/MAP_118/results/FASTQ_files/BG505_mut_virus_rep3d_1-18-8ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>wt-DNA-rep2</td>
      <td>../../2019/MAP_118/results/FASTQ_files/BG505_wt-DNA-rep2_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>wt-DNA-rep3</td>
      <td>../../2019/MAP_118/results/FASTQ_files/BG505_wt-DNA-rep3_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>wt-DNA-rep1</td>
      <td>../../2019/MAP_118/results/FASTQ_files/BG505_wt-DNA-rep1_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-IDC0136-5500ug</td>
      <td>./FASTQ_files/BG-M12b-IDC0136-5500ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-IDC0136-7000ug</td>
      <td>./FASTQ_files/BG-M12b-IDC0136-7000ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-IDC0136-8500ug</td>
      <td>./FASTQ_files/BG-M12b-IDC0136-8500ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-IDC0337-4000ug</td>
      <td>./FASTQ_files/BG-M12b-IDC0337-4000ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-IDC0337-4500ug</td>
      <td>./FASTQ_files/BG-M12b-IDC0337-4500ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-IDC403-7500ug</td>
      <td>./FASTQ_files/BG-M12b-IDC403-7500ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-IDC403-6000ug</td>
      <td>./FASTQ_files/BG-M12b-IDC403-6000ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-IDC208-6000ug</td>
      <td>./FASTQ_files/BG-M12b-IDC208-6000ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-IDC208-5000ug</td>
      <td>./FASTQ_files/BG-M12b-IDC208-5000ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-IDC208-4000ug</td>
      <td>./FASTQ_files/BG-M12b-IDC208-4000ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-IDC0337-6000ug</td>
      <td>./FASTQ_files/BG-M12b-IDC0337-6000ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b</td>
      <td>./FASTQ_files/BG-M12b-mock_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-561-IgG-2800ug</td>
      <td>./FASTQ_files/BG-M12b-561-IgG-2800ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-561-IgG-2200ug</td>
      <td>./FASTQ_files/BG-M12b-561-IgG-2200ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-212-6ug</td>
      <td>./FASTQ_files/BG-M12b-2-12-6ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-QA013-2282dpi-d40</td>
      <td>./FASTQ_files/BG-M12b-QA013-2282dpi-d40_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-QA013-2282dpi-d25</td>
      <td>./FASTQ_files/BG-M12b-QA013-2282dpi-d25_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2b-IDC403-8500ug</td>
      <td>./FASTQ_files/BG-M12b-IDC403-8500ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3e-212-8ug</td>
      <td>../../2019/CD4bs_bnAb_comparison/FASTQ_files/BG505_BG5-M3E-2-12-8ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3e</td>
      <td>../../2019/CD4bs_bnAb_comparison/FASTQ_files/BG505_BG5-M3E-mock_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3e-QA013-2282dpi-d7-5</td>
      <td>../../2019/CD4bs_bnAb_comparison/FASTQ_files/BG505_BG5-M3E-QA013-2282dpi-1-7-5_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3e-561-IgG-2500ug</td>
      <td>../../2019/CD4bs_bnAb_comparison/FASTQ_files/BG505_BG5-M3E-561-IgG-2500ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2a-QA013-2282dpi-d15</td>
      <td>../../2019/CD4bs_bnAb_comparison/FASTQ_files/BG505_BG5-M1v2A-QA013-2282dpi-1-15_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2a-QA013-2282dpi-d10</td>
      <td>../../2019/CD4bs_bnAb_comparison/FASTQ_files/BG505_BG5-M1v2A-QA013-2282dpi-1-10_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2a</td>
      <td>../../2019/CD4bs_bnAb_comparison/FASTQ_files/BG505_BG5-M1v2A-mock_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2a-561-IgG-3500ug</td>
      <td>../../2019/CD4bs_bnAb_comparison/FASTQ_files/BG505_BG5-M1v2A-561-IgG-3500ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2a-561-IgG-3000ug</td>
      <td>../../2019/CD4bs_bnAb_comparison/FASTQ_files/BG505_BG5-M1v2A-561-IgG-3000ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2a-212-10ug</td>
      <td>../../2019/CD4bs_bnAb_comparison/FASTQ_files/BG505_BG5-M1v2A-2-12-10ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep1v2a-212-7ug</td>
      <td>../../2019/CD4bs_bnAb_comparison/FASTQ_files/BG505_BG5-M1v2A-2-12-7ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep2-QA013-2-2ug</td>
      <td>../../2017/BG505_MAP/FASTQ_files/BG505_mut_virus_rep2_QA013-2_2ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep2-QA013-2-4ug</td>
      <td>../../2017/BG505_MAP/FASTQ_files/BG505_mut_virus_rep2_QA013-2_4ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3-QA013-2-4ug</td>
      <td>../../2017/BG505_MAP/FASTQ_files/BG505_mut_virus_rep3_QA013-2_4ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep2</td>
      <td>../../2017/BG505_MAP/FASTQ_files/BG505_mut_virus_rep2_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3</td>
      <td>../../2017/BG505_MAP/FASTQ_files/BG505_mut_virus_rep3_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3v2a</td>
      <td>./FASTQ_files/BG-M32a-mock_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3v2a-IDC0191-6000ug</td>
      <td>./FASTQ_files/BG-M32a-IDC0191-6000ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3v2a-IDC0245-8000ug</td>
      <td>./FASTQ_files/BG-M32a-IDC0245-8000ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3v2a-IDC0396-9000ug</td>
      <td>./FASTQ_files/BG-M32a-IDC0396-9000ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3v2a-IDC0513-5000ug</td>
      <td>./FASTQ_files/BG-M32a-IDC0513-5000ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3v2a-IDC0518-8000ug</td>
      <td>./FASTQ_files/BG-M32a-IDC0518-8000ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3v2a-IDC0546-8500ug</td>
      <td>./FASTQ_files/BG-M32a-IDC0546-8500ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3v2a-IDF033-1400ug</td>
      <td>./FASTQ_files/BG-M32a-IDF033-1400ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3v2a-IDF033-1800ug</td>
      <td>./FASTQ_files/BG-M32a-IDF033-1800ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3v2a-IDF065-1400ug</td>
      <td>./FASTQ_files/BG-M32a-IDF065-1400ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3v2a-IDF141-6000ug</td>
      <td>./FASTQ_files/BG-M32a-IDF141-6000ug_R1.fastq.gz</td>
    </tr>
    <tr>
      <td>mtvir-rep3v2a-IDF147-5300ug</td>
      <td>./FASTQ_files/BG-M32a-IDF147-5300ug_R1.fastq.gz</td>
    </tr>
  </tbody>
</table>


    
    Now running dms2_batch_bcsubamp...
    Completed dms2_batch_bcsubamp.


Now we look at the summary plots created by `dms2_batch_bcsubamp`. 
All of these files are found in the directory specified by --outdir, and with the prefix specified by summaryprefix. 
So we define them using this plot prefix plus the suffix for each plot.
Note that these files all refer to sites in sequential 1, 2, ... numbering of the BF520 sequence.


```python
countsplotprefix = os.path.join(countsdir, 'summary')
```

The `*_readstats.pdf` plot below shows the statistics on the reads. This plot shows that most of the reads were retained, and a small fraction discarded because of low-quality barcodes. None failed the Illumina filter as those were already filtered out those when we downloaded from the SRA.


```python
showPDF(countsplotprefix + '_readstats.pdf', width=600)
```


![png](analysis_notebook_files/analysis_notebook_14_0.png)


The `*_readsperbc.pdf` plot below shows how many times different barcodes were observed for each sample. Barcodes need to be observed multiple times to be useful for barcoded-subamplicon sequencing error correction.


```python
showPDF(countsplotprefix + '_readsperbc.pdf')
```


![png](analysis_notebook_files/analysis_notebook_16_0.png)


The `*_bcstats.pdf` plot below shows statistics on the barcodes. Some of the barcodes had to be discarded because they had too few reads (these are the single-read barcodes in the plot above), a small fraction with adequate reads were not alignable, and the rest aligned to the Env gene properly.
This plot and the one above suggest that probably could have gotten additional depth by sequencing more, since then we would have had more barcodes with multiple reads.


```python
showPDF(countsplotprefix + '_bcstats.pdf', width=600)
```


![png](analysis_notebook_files/analysis_notebook_18_0.png)


The `*_depth.pdf` plot below shows the depth (number of called codons) at each site in the gene. 
For most of the samples, the depth across the gene is fairly uniform, indicating that the subamplicons were pooled fairly evenly.  
Note that some samples (in particular the mock samples) were intentionally sequenced to higher depth than the antibody-selected samples, as we expect to see more diversity in the mock samples. 
Note that the gene was not sequenced past codon site 691, and so there is no coverage there.


```python
showPDF(countsplotprefix + '_depth.pdf')
```


![png](analysis_notebook_files/analysis_notebook_20_0.png)


The `*_mutfreq.pdf` plot below shows the per-codon frequency of mutations at each site. 
For each antibody-selected sample, we see a few sites of clear peaks in mutation frequency. 
These peaks tend to occur at the same sites in different replicates, and so presumably represent the sites where antibody-escape mutations are selected. 
There are no such peaks for the mock sample since there is no antibody selection to favor specific mutations.
Note also that the gene was not mutagenized or sequenced past codon site 691, so there are no mutations there.


```python
showPDF(countsplotprefix + '_mutfreq.pdf')
```


![png](analysis_notebook_files/analysis_notebook_22_0.png)


The `*_codonmuttypes.pdf` plot below shows the per-codon frequency of nonsynonymous, synonymous, and stop codon mutations across the entire gene. 


```python
showPDF(countsplotprefix + '_codonmuttypes.pdf', width=600)
```


![png](analysis_notebook_files/analysis_notebook_24_0.png)


The `*_codonntchanges.pdf` plot below shows same data as above but categorizes codon mutations by the number of nucleotides that are changed (e.g., ATG to AAG changes 1 nucleotide, ATG to AAC changes 2 nucleotides, and ATG to CAC changes 3 nucleotides).


```python
showPDF(countsplotprefix + '_codonntchanges.pdf', width=600)
```


![png](analysis_notebook_files/analysis_notebook_26_0.png)


The `*_singlentchanges.pdf` plot below shows the frequency of each type of nucleotide change among only codon mutations with one nucleotide change. This plot is mostly useful to check if there is a large bias in which mutations appear. In particular, if you are getting oxidative damage (which causes G to T mutations) during the library preparation process, you will see a large excess of C to A or G to T mutations (or both). There is not much oxidative damage in the samples plotted below, which mostly have a fairly even distribution of nucleotide changes.

We do see that transitions (G <-> A and C <-> T) are a bit more common than most of the other types of mutations (all of which are transversions). This is expected, since [PCR based sources of mutations are expected to preferentially introduce transitions](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3720931/). 

This plot would also be important to examine any sign of APOBEC hypermutation (also G <-> A, with preferences for specific motifs) occuring in the HIV genome. One of the reasons we selected the SupT1.R5 cell line was because [there is very little of the APOBEC proteins expressed in  a subclone of SupT1 cells](https://www.ncbi.nlm.nih.gov/pubmed/20308164). In the viral samples, there does appear to be an increase in G to A and C to T mutations. However, we passaged WT virus in these cells and confirmed there were not signs of extensive APOBEC hypermutation in hotspot motifs via deep sequencing of env after passaging. Thus, we do not think there is extensive APOBEC hypermutation in our data, and this increase certain types of mutations in the viral samples is likely due to bias in the viral reverse transcription.


```python
showPDF(countsplotprefix + '_singlentchanges.pdf', width=600)
```


![png](analysis_notebook_files/analysis_notebook_28_0.png)


## Renumber codon counts to HXB2 numbering 
The standard numbering scheme for HIV Env is the [HXB2 numbering scheme](https://www.hiv.lanl.gov/content/sequence/HIV/REVIEWS/HXB2.html).
The file [./data/BF520c2_to_HXB2.csv](./data/BF520c2_to_HXB2.csv) gives the mapping from sequential 1, 2, ... numbering of the BF520 protein sequence to the HXB2 numbering scheme. 
This file was generated by aligning the HXB2 sequence [taken from Genbank](http://www.ncbi.nlm.nih.gov/protein/1906385) with the BF520c2 sequence using the [LANL alignment interface](http://www.hiv.lanl.gov/cgi-bin/VIRALIGN/viralign.cgi) at the protein sequence level. 
Insertions relative to HXB2 are given letter suffixes as [described here](http://www.hiv.lanl.gov/content/sequence/HIV/REVIEWS/HXB2.html).

Additionally, not all residues in BF520 Env were mutagenized. 
The N-terminal signal peptide and the C-terminal cytoplasmic tail were excluded because they seem likely to affect expression level. 
These sites are not listed in the renumbering file, and so are dropped when we do the re-numbering.

To do the re-numbering, we use the [dms_tools2.utils.renumberSites](https://jbloomlab.github.io/dms_tools2/dms_tools2.utils.html#dms_tools2.utils.renumberSites) function from the [dms_tools Python API](https://jbloomlab.github.io/dms_tools2/api.html) to create a new directory that contains all of the re-numbered files with the same name as in the original codon counts directory produced above.


```python
renumberedcountsdir = os.path.join(resultsdir, 'renumberedcounts')
os.makedirs(renumberedcountsdir, exist_ok=True)
```


```python
renumberfile = './data/BG505_to_HXB2_new.csv'

# renumbered counts will go here
renumberedcountsdir = os.path.join(resultsdir, 'renumberedcounts')
import glob
# counts files to renumber
countsfiles = glob.glob('{0}/*codoncounts.csv'.format(countsdir))

dms_tools2.utils.renumberSites(renumberfile, countsfiles, missing='drop', 
        outdir=renumberedcountsdir)
```


### Compute immune selection
Now we run [dms2_batch_diffsel](https://jbloomlab.github.io/dms_tools2/dms2_batch_diffsel.html) to compute the immune selection.
We then add to our `selections` data frame the name of the files holding the computed site (*site*) and mutation (*mut*) level selection for each sample.


```python
diffseldir = os.path.join(resultsdir, 'diffsel')
os.makedirs(diffseldir, exist_ok=True)

# write batch file used by program
batchdf = pd.read_csv("./data/diffselbatch.csv")
display(HTML(batchdf.to_html(index=False)))

batchfile = os.path.join(diffseldir, 'batch.csv')
(batchdf.to_csv(batchfile, index=False)
 )

```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>group</th>
      <th>name</th>
      <th>sel</th>
      <th>mock</th>
      <th>err</th>
      <th>mds_names</th>
      <th>libfracsurvive</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>118</td>
      <td>4ug-rep2d</td>
      <td>mtvir-rep2d-118-4ug</td>
      <td>mtvir-rep2d</td>
      <td>wt-DNA-rep2</td>
      <td>118-4ug-rep2d</td>
      <td>0.005177</td>
    </tr>
    <tr>
      <td>118</td>
      <td>8ug-rep2d</td>
      <td>mtvir-rep2d-118-8ug</td>
      <td>mtvir-rep2d</td>
      <td>wt-DNA-rep2</td>
      <td>118-8ug-rep2d</td>
      <td>0.000168</td>
    </tr>
    <tr>
      <td>IDC561</td>
      <td>2500ug-rep3e</td>
      <td>mtvir-rep3e-561-IgG-2500ug</td>
      <td>mtvir-rep3e</td>
      <td>wt-DNA-rep3</td>
      <td>IDC561-2500ug-rep3e</td>
      <td>0.023230</td>
    </tr>
    <tr>
      <td>118</td>
      <td>4ug-rep3d</td>
      <td>mtvir-rep3d-118-4ug</td>
      <td>mtvir-rep3d</td>
      <td>wt-DNA-rep3</td>
      <td>118-4ug-rep3d</td>
      <td>0.005182</td>
    </tr>
    <tr>
      <td>118</td>
      <td>8ug-rep3d</td>
      <td>mtvir-rep3d-118-8ug</td>
      <td>mtvir-rep3d</td>
      <td>wt-DNA-rep3</td>
      <td>118-8ug-rep3d</td>
      <td>0.000673</td>
    </tr>
    <tr>
      <td>212</td>
      <td>8ug-rep3e</td>
      <td>mtvir-rep3e-212-8ug</td>
      <td>mtvir-rep3e</td>
      <td>wt-DNA-rep3</td>
      <td>212-8ug-rep3e</td>
      <td>0.006215</td>
    </tr>
    <tr>
      <td>QA013-2282dpi</td>
      <td>d7-5-rep3e</td>
      <td>mtvir-rep3e-QA013-2282dpi-d7-5</td>
      <td>mtvir-rep3e</td>
      <td>wt-DNA-rep3</td>
      <td>QA013-2282dpi-d7-5-rep3e</td>
      <td>0.000388</td>
    </tr>
    <tr>
      <td>QA013-2282dpi</td>
      <td>d10-rep1v2a</td>
      <td>mtvir-rep1v2a-QA013-2282dpi-d10</td>
      <td>mtvir-rep1v2a</td>
      <td>wt-DNA-rep1</td>
      <td>QA013-2282dpi-d10-rep1v2a</td>
      <td>0.000317</td>
    </tr>
    <tr>
      <td>QA013-2282dpi</td>
      <td>d15-rep1v2a</td>
      <td>mtvir-rep1v2a-QA013-2282dpi-d15</td>
      <td>mtvir-rep1v2a</td>
      <td>wt-DNA-rep1</td>
      <td>QA013-2282dpi-d15-rep1v2a</td>
      <td>0.001444</td>
    </tr>
    <tr>
      <td>IDC561</td>
      <td>3000ug-rep1v2a</td>
      <td>mtvir-rep1v2a-561-IgG-3000ug</td>
      <td>mtvir-rep1v2a</td>
      <td>wt-DNA-rep1</td>
      <td>IDC561-3000ug-rep1v2a</td>
      <td>0.002478</td>
    </tr>
    <tr>
      <td>IDC561</td>
      <td>3500ug-rep1v2a</td>
      <td>mtvir-rep1v2a-561-IgG-3500ug</td>
      <td>mtvir-rep1v2a</td>
      <td>wt-DNA-rep1</td>
      <td>IDC561-3500ug-rep1v2a</td>
      <td>0.001820</td>
    </tr>
    <tr>
      <td>212</td>
      <td>10ug-rep1v2a</td>
      <td>mtvir-rep1v2a-212-10ug</td>
      <td>mtvir-rep1v2a</td>
      <td>wt-DNA-rep1</td>
      <td>212-10ug-rep1v2a</td>
      <td>0.002501</td>
    </tr>
    <tr>
      <td>212</td>
      <td>7ug-rep1v2a</td>
      <td>mtvir-rep1v2a-212-7ug</td>
      <td>mtvir-rep1v2a</td>
      <td>wt-DNA-rep1</td>
      <td>212-7ug-rep1v2a</td>
      <td>0.007341</td>
    </tr>
    <tr>
      <td>IDC0136</td>
      <td>5500ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC0136-5500ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC0136-5500ug-rep1v2b</td>
      <td>0.136001</td>
    </tr>
    <tr>
      <td>IDC0136</td>
      <td>7000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC0136-7000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC0136-7000ug-rep1v2b</td>
      <td>0.062973</td>
    </tr>
    <tr>
      <td>IDC0136</td>
      <td>8500ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC0136-8500ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC0136-8500ug-rep1v2b</td>
      <td>0.041266</td>
    </tr>
    <tr>
      <td>IDC0337</td>
      <td>4000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC0337-4000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC0337-4000ug-rep1v2b</td>
      <td>0.102009</td>
    </tr>
    <tr>
      <td>IDC0337</td>
      <td>4500ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC0337-4500ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC0337-4500ug-rep1v2b</td>
      <td>0.091161</td>
    </tr>
    <tr>
      <td>IDC403</td>
      <td>7500ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC403-7500ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC403-7500ug-rep1v2b</td>
      <td>0.000893</td>
    </tr>
    <tr>
      <td>IDC403</td>
      <td>6000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC403-6000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC403-6000ug-rep1v2b</td>
      <td>0.001969</td>
    </tr>
    <tr>
      <td>IDC208</td>
      <td>6000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC208-6000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC208-6000ug-rep1v2b</td>
      <td>0.068419</td>
    </tr>
    <tr>
      <td>IDC208</td>
      <td>5000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC208-5000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC208-5000ug-rep1v2b</td>
      <td>0.103581</td>
    </tr>
    <tr>
      <td>IDC208</td>
      <td>4000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC208-4000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC208-4000ug-rep1v2b</td>
      <td>0.150962</td>
    </tr>
    <tr>
      <td>IDC0337</td>
      <td>6000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC0337-6000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC0337-6000ug-rep1v2b</td>
      <td>0.059010</td>
    </tr>
    <tr>
      <td>IDC561</td>
      <td>2800ug-rep1v2b</td>
      <td>mtvir-rep1v2b-561-IgG-2800ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC561-2800ug-rep1v2b</td>
      <td>0.007367</td>
    </tr>
    <tr>
      <td>IDC561</td>
      <td>2200ug-rep1v2b</td>
      <td>mtvir-rep1v2b-561-IgG-2200ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC561-2200ug-rep1v2b</td>
      <td>0.012112</td>
    </tr>
    <tr>
      <td>212</td>
      <td>6ug-rep1v2b</td>
      <td>mtvir-rep1v2b-212-6ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>212-6ug-rep1v2b</td>
      <td>0.034943</td>
    </tr>
    <tr>
      <td>QA013-2282dpi</td>
      <td>d40-rep1v2b</td>
      <td>mtvir-rep1v2b-QA013-2282dpi-d40</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>QA013-2282dpi-d40-rep1v2b</td>
      <td>0.470669</td>
    </tr>
    <tr>
      <td>QA013-2282dpi</td>
      <td>d25-rep1v2b</td>
      <td>mtvir-rep1v2b-QA013-2282dpi-d25</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>QA013-2282dpi-d25-rep1v2b</td>
      <td>0.278701</td>
    </tr>
    <tr>
      <td>IDC403</td>
      <td>8500ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC403-8500ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC403-8500ug-rep1v2b</td>
      <td>0.000359</td>
    </tr>
    <tr>
      <td>QA013-2</td>
      <td>2ug-rep2</td>
      <td>mtvir-rep2-QA013-2-2ug</td>
      <td>mtvir-rep2</td>
      <td>wt-DNA-rep2</td>
      <td>QA013-2-2ug-rep2</td>
      <td>0.032488</td>
    </tr>
    <tr>
      <td>QA013-2</td>
      <td>4ug-rep2</td>
      <td>mtvir-rep2-QA013-2-4ug</td>
      <td>mtvir-rep2</td>
      <td>wt-DNA-rep2</td>
      <td>QA013-2-4ug-rep2</td>
      <td>0.002637</td>
    </tr>
    <tr>
      <td>QA013-2</td>
      <td>4ug-rep3</td>
      <td>mtvir-rep3-QA013-2-4ug</td>
      <td>mtvir-rep3</td>
      <td>wt-DNA-rep3</td>
      <td>QA013-2-4ug-rep3</td>
      <td>0.001661</td>
    </tr>
    <tr>
      <td>IDC0191</td>
      <td>6000ug-rep3v2a</td>
      <td>mtvir-rep3v2a-IDC0191-6000ug</td>
      <td>mtvir-rep3v2a</td>
      <td>wt-DNA-rep3</td>
      <td>IDC0191-6000ug-rep3v2a</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>IDC0245</td>
      <td>8000ug-rep3v2a</td>
      <td>mtvir-rep3v2a-IDC0245-8000ug</td>
      <td>mtvir-rep3v2a</td>
      <td>wt-DNA-rep3</td>
      <td>IDC0245-8000ug-rep3v2a</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>IDC0396</td>
      <td>9000ug-rep3v2a</td>
      <td>mtvir-rep3v2a-IDC0396-9000ug</td>
      <td>mtvir-rep3v2a</td>
      <td>wt-DNA-rep3</td>
      <td>IDC0396-9000ug-rep3v2a</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>IDC0513</td>
      <td>5000ug-rep3v2a</td>
      <td>mtvir-rep3v2a-IDC0513-5000ug</td>
      <td>mtvir-rep3v2a</td>
      <td>wt-DNA-rep3</td>
      <td>IDC0513-5000ug-rep3v2a</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>IDC0518</td>
      <td>8000ug-rep3v2a</td>
      <td>mtvir-rep3v2a-IDC0518-8000ug</td>
      <td>mtvir-rep3v2a</td>
      <td>wt-DNA-rep3</td>
      <td>IDC0518-8000ug-rep3v2a</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>IDC0546</td>
      <td>8500ug-rep3v2a</td>
      <td>mtvir-rep3v2a-IDC0546-8500ug</td>
      <td>mtvir-rep3v2a</td>
      <td>wt-DNA-rep3</td>
      <td>IDC0546-8500ug-rep3v2a</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>IDF033</td>
      <td>1400ug-rep3v2a</td>
      <td>mtvir-rep3v2a-IDF033-1400ug</td>
      <td>mtvir-rep3v2a</td>
      <td>wt-DNA-rep3</td>
      <td>IDF033-1400ug-rep3v2a</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>IDF033</td>
      <td>1800ug-rep3v2a</td>
      <td>mtvir-rep3v2a-IDF033-1800ug</td>
      <td>mtvir-rep3v2a</td>
      <td>wt-DNA-rep3</td>
      <td>IDF033-1800ug-rep3v2a</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>IDF065</td>
      <td>1400ug-rep3v2a</td>
      <td>mtvir-rep3v2a-IDF065-1400ug</td>
      <td>mtvir-rep3v2a</td>
      <td>wt-DNA-rep3</td>
      <td>IDF065-1400ug-rep3v2a</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>IDF141</td>
      <td>6000ug-rep3v2a</td>
      <td>mtvir-rep3v2a-IDF141-6000ug</td>
      <td>mtvir-rep3v2a</td>
      <td>wt-DNA-rep3</td>
      <td>IDF141-6000ug-rep3v2a</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>IDF147</td>
      <td>5300ug-rep3v2a</td>
      <td>mtvir-rep3v2a-IDF147-5300ug</td>
      <td>mtvir-rep3v2a</td>
      <td>wt-DNA-rep3</td>
      <td>IDF147-5300ug-rep3v2a</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>



```python
cmds = ['dms2_batch_diffsel',
        '--summaryprefix', 'summary',
        '--batchfile', batchfile,
        '--outdir', diffseldir,
        '--indir', renumberedcountsdir,
        '--use_existing', use_existing, #config['use_existing'],
        '--ncpus', "{0}".format(ncpus),
        ]

print(f"Computing diffsel using dms2_batch_diffsel with command:\n{' '.join(cmds)}")
subprocess.check_output(cmds)
```

    Computing diffsel using dms2_batch_diffsel with command:
    dms2_batch_diffsel --summaryprefix summary --batchfile ./results/diffsel/batch.csv --outdir ./results/diffsel --indir ./results/renumberedcounts --use_existing yes --ncpus 14





    b''




```python
ind_names = list(set(batchdf["mds_names"].tolist()))
```


```python
for antibody in ind_names:
    mutdiffsel = os.path.join(diffseldir, '{0}_mutdiffsel.csv'.format(antibody))
    
    #scale bar unit is maximum effect
    mutdiffseldf = pd.read_csv(mutdiffsel)
    scaleunit = '{0:.1g}'.format(mutdiffseldf['mutdiffsel'].max())
    #print(scaleunit)
    scalelabel = '"differential selection = {0}"'.format(scaleunit)
    logoplot = os.path.join(diffseldir, '{0}_diffsel.pdf'.format(antibody))
    #print(logoplot)
    logoname = '{0}'.format(antibody)
    print("\nCreating logo plot for {0} from {1}".format(antibody, mutdiffsel))
    log = !dms2_logoplot \
            --diffsel {mutdiffsel} \
            --name {logoname} \
            --outdir {diffseldir} \
            --restrictdiffsel positive \
            --sepline no \
            --nperline 84 \
            --overlay1 {mutdiffsel} wildtype wildtype \
            --scalebar {scaleunit} {scalelabel} \
            --underlay yes \
            --use_existing {use_existing}
    #showPDF(logoplot)
```

    
    Creating logo plot for QA013-2282dpi-d10-rep1v2a from ./results/diffsel/QA013-2282dpi-d10-rep1v2a_mutdiffsel.csv
    
    Creating logo plot for QA013-2-2ug-rep2 from ./results/diffsel/QA013-2-2ug-rep2_mutdiffsel.csv
    
    Creating logo plot for IDC403-6000ug-rep1v2b from ./results/diffsel/IDC403-6000ug-rep1v2b_mutdiffsel.csv
    
    Creating logo plot for 212-7ug-rep1v2a from ./results/diffsel/212-7ug-rep1v2a_mutdiffsel.csv
    
    Creating logo plot for IDC0136-7000ug-rep1v2b from ./results/diffsel/IDC0136-7000ug-rep1v2b_mutdiffsel.csv
    
    Creating logo plot for QA013-2282dpi-d7-5-rep3e from ./results/diffsel/QA013-2282dpi-d7-5-rep3e_mutdiffsel.csv
    
    Creating logo plot for IDC561-3500ug-rep1v2a from ./results/diffsel/IDC561-3500ug-rep1v2a_mutdiffsel.csv
    
    Creating logo plot for 118-8ug-rep3d from ./results/diffsel/118-8ug-rep3d_mutdiffsel.csv
    
    Creating logo plot for IDC0337-6000ug-rep1v2b from ./results/diffsel/IDC0337-6000ug-rep1v2b_mutdiffsel.csv
    
    Creating logo plot for IDC561-3000ug-rep1v2a from ./results/diffsel/IDC561-3000ug-rep1v2a_mutdiffsel.csv
    
    Creating logo plot for IDC561-2800ug-rep1v2b from ./results/diffsel/IDC561-2800ug-rep1v2b_mutdiffsel.csv
    
    Creating logo plot for IDC208-4000ug-rep1v2b from ./results/diffsel/IDC208-4000ug-rep1v2b_mutdiffsel.csv
    
    Creating logo plot for IDF065-1400ug-rep3v2a from ./results/diffsel/IDF065-1400ug-rep3v2a_mutdiffsel.csv
    
    Creating logo plot for IDF033-1400ug-rep3v2a from ./results/diffsel/IDF033-1400ug-rep3v2a_mutdiffsel.csv
    
    Creating logo plot for 212-6ug-rep1v2b from ./results/diffsel/212-6ug-rep1v2b_mutdiffsel.csv



    ---------------------------------------------------------------------------

    FileNotFoundError                         Traceback (most recent call last)

    <ipython-input-21-9a31de43bfd9> in <module>
          3 
          4     #scale bar unit is maximum effect
    ----> 5     mutdiffseldf = pd.read_csv(mutdiffsel)
          6     scaleunit = '{0:.1g}'.format(mutdiffseldf['mutdiffsel'].max())
          7     #print(scaleunit)


    /fh/fast/bloom_j/software/conda_v2/envs/BloomLab/lib/python3.6/site-packages/pandas/io/parsers.py in parser_f(filepath_or_buffer, sep, delimiter, header, names, index_col, usecols, squeeze, prefix, mangle_dupe_cols, dtype, engine, converters, true_values, false_values, skipinitialspace, skiprows, skipfooter, nrows, na_values, keep_default_na, na_filter, verbose, skip_blank_lines, parse_dates, infer_datetime_format, keep_date_col, date_parser, dayfirst, cache_dates, iterator, chunksize, compression, thousands, decimal, lineterminator, quotechar, quoting, doublequote, escapechar, comment, encoding, dialect, error_bad_lines, warn_bad_lines, delim_whitespace, low_memory, memory_map, float_precision)
        683         )
        684 
    --> 685         return _read(filepath_or_buffer, kwds)
        686 
        687     parser_f.__name__ = name


    /fh/fast/bloom_j/software/conda_v2/envs/BloomLab/lib/python3.6/site-packages/pandas/io/parsers.py in _read(filepath_or_buffer, kwds)
        455 
        456     # Create the parser.
    --> 457     parser = TextFileReader(fp_or_buf, **kwds)
        458 
        459     if chunksize or iterator:


    /fh/fast/bloom_j/software/conda_v2/envs/BloomLab/lib/python3.6/site-packages/pandas/io/parsers.py in __init__(self, f, engine, **kwds)
        893             self.options["has_index_names"] = kwds["has_index_names"]
        894 
    --> 895         self._make_engine(self.engine)
        896 
        897     def close(self):


    /fh/fast/bloom_j/software/conda_v2/envs/BloomLab/lib/python3.6/site-packages/pandas/io/parsers.py in _make_engine(self, engine)
       1133     def _make_engine(self, engine="c"):
       1134         if engine == "c":
    -> 1135             self._engine = CParserWrapper(self.f, **self.options)
       1136         else:
       1137             if engine == "python":


    /fh/fast/bloom_j/software/conda_v2/envs/BloomLab/lib/python3.6/site-packages/pandas/io/parsers.py in __init__(self, src, **kwds)
       1915         kwds["usecols"] = self.usecols
       1916 
    -> 1917         self._reader = parsers.TextReader(src, **kwds)
       1918         self.unnamed_cols = self._reader.unnamed_cols
       1919 


    pandas/_libs/parsers.pyx in pandas._libs.parsers.TextReader.__cinit__()


    pandas/_libs/parsers.pyx in pandas._libs.parsers.TextReader._setup_parser_source()


    FileNotFoundError: [Errno 2] File b'./results/diffsel/IDC0513-5000ug-rep3v2a_mutdiffsel.csv' does not exist: b'./results/diffsel/IDC0513-5000ug-rep3v2a_mutdiffsel.csv'



```python
names = list(set(batchdf["group"].tolist()))
names
```




    ['IDC0513',
     'IDF033',
     'IDC0546',
     'IDF141',
     'QA013-2',
     'IDF147',
     'IDC0396',
     '118',
     'IDF065',
     '212',
     'IDC0518',
     'QA013-2282dpi',
     'IDC0136',
     'IDC0337',
     'IDC208',
     'IDC0191',
     'IDC561',
     'IDC0245',
     'IDC403']




```python
use_existing
```


```python
for antibody in names:
    mutdiffsel = os.path.join(diffseldir, 'summary_{0}-medianmutdiffsel.csv'.format(antibody))
    
    #scale bar unit is maximum effect
    mutdiffseldf = pd.read_csv(mutdiffsel)
    scaleunit = '{0:.1g}'.format(mutdiffseldf['mutdiffsel'].max())
    #print(scaleunit)
    scalelabel = '"differential selection = {0}"'.format(scaleunit)
    logoplot = os.path.join(diffseldir, '{0}-median_diffsel.pdf'.format(antibody))
    #print(logoplot)
    logoname = '{0}-median'.format(antibody)
    print("\nCreating logo plot for {0} from {1}".format(antibody, mutdiffsel))
    log = !dms2_logoplot \
            --diffsel {mutdiffsel} \
            --name {logoname} \
            --outdir {diffseldir} \
            --restrictdiffsel positive \
            --sepline no \
            --nperline 84 \
            --overlay1 {mutdiffsel} wildtype wildtype \
            --scalebar {scaleunit} {scalelabel} \
            --underlay yes \
            --use_existing {use_existing}
    #showPDF(logoplot)
```


```python
for antibody in names:
    mutdiffsel = os.path.join(diffseldir, 'summary_{0}-meanmutdiffsel.csv'.format(antibody))
    
    #scale bar unit is maximum effect
    mutdiffseldf = pd.read_csv(mutdiffsel)
    scaleunit = '{0:.1g}'.format(mutdiffseldf['mutdiffsel'].max())
    #print(scaleunit)
    scalelabel = '"differential selection = {0}"'.format(scaleunit)
    logoplot = os.path.join(diffseldir, '{0}-mean_diffsel.pdf'.format(antibody))
    #print(logoplot)
    logoname = '{0}-mean'.format(antibody)
    print("\nCreating logo plot for {0} from {1}".format(antibody, mutdiffsel))
    log = !dms2_logoplot \
            --diffsel {mutdiffsel} \
            --name {logoname} \
            --outdir {diffseldir} \
            --restrictdiffsel positive \
            --sepline no \
            --nperline 84 \
            --overlay1 {mutdiffsel} wildtype wildtype \
            --scalebar {scaleunit} {scalelabel} \
            --underlay yes \
            --use_existing {use_existing}
    #showPDF(logoplot)
```

Put it all together and plot...
Get data in df
plot epitope zooms by region. 



```python
batchdf
```


```python
selections = batchdf.copy()
#add % from fracsurvive
selections['percent_infectivity'] =  selections['libfracsurvive'] * 100
selections.sort_values(by=['group', 'percent_infectivity'], inplace=True)#)
selections
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>group</th>
      <th>name</th>
      <th>sel</th>
      <th>mock</th>
      <th>err</th>
      <th>mds_names</th>
      <th>libfracsurvive</th>
      <th>percent_infectivity</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1</th>
      <td>118</td>
      <td>8ug-rep2d</td>
      <td>mtvir-rep2d-118-8ug</td>
      <td>mtvir-rep2d</td>
      <td>wt-DNA-rep2</td>
      <td>118-8ug-rep2d</td>
      <td>0.000168</td>
      <td>0.016777</td>
    </tr>
    <tr>
      <th>4</th>
      <td>118</td>
      <td>8ug-rep3d</td>
      <td>mtvir-rep3d-118-8ug</td>
      <td>mtvir-rep3d</td>
      <td>wt-DNA-rep3</td>
      <td>118-8ug-rep3d</td>
      <td>0.000673</td>
      <td>0.067347</td>
    </tr>
    <tr>
      <th>0</th>
      <td>118</td>
      <td>4ug-rep2d</td>
      <td>mtvir-rep2d-118-4ug</td>
      <td>mtvir-rep2d</td>
      <td>wt-DNA-rep2</td>
      <td>118-4ug-rep2d</td>
      <td>0.005177</td>
      <td>0.517724</td>
    </tr>
    <tr>
      <th>3</th>
      <td>118</td>
      <td>4ug-rep3d</td>
      <td>mtvir-rep3d-118-4ug</td>
      <td>mtvir-rep3d</td>
      <td>wt-DNA-rep3</td>
      <td>118-4ug-rep3d</td>
      <td>0.005182</td>
      <td>0.518168</td>
    </tr>
    <tr>
      <th>11</th>
      <td>212</td>
      <td>10ug-rep1v2a</td>
      <td>mtvir-rep1v2a-212-10ug</td>
      <td>mtvir-rep1v2a</td>
      <td>wt-DNA-rep1</td>
      <td>212-10ug-rep1v2a</td>
      <td>0.002501</td>
      <td>0.250140</td>
    </tr>
    <tr>
      <th>5</th>
      <td>212</td>
      <td>8ug-rep3e</td>
      <td>mtvir-rep3e-212-8ug</td>
      <td>mtvir-rep3e</td>
      <td>wt-DNA-rep3</td>
      <td>212-8ug-rep3e</td>
      <td>0.006215</td>
      <td>0.621478</td>
    </tr>
    <tr>
      <th>12</th>
      <td>212</td>
      <td>7ug-rep1v2a</td>
      <td>mtvir-rep1v2a-212-7ug</td>
      <td>mtvir-rep1v2a</td>
      <td>wt-DNA-rep1</td>
      <td>212-7ug-rep1v2a</td>
      <td>0.007341</td>
      <td>0.734148</td>
    </tr>
    <tr>
      <th>26</th>
      <td>212</td>
      <td>6ug-rep1v2b</td>
      <td>mtvir-rep1v2b-212-6ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>212-6ug-rep1v2b</td>
      <td>0.034943</td>
      <td>3.494317</td>
    </tr>
    <tr>
      <th>15</th>
      <td>IDC0136</td>
      <td>8500ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC0136-8500ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC0136-8500ug-rep1v2b</td>
      <td>0.041266</td>
      <td>4.126587</td>
    </tr>
    <tr>
      <th>14</th>
      <td>IDC0136</td>
      <td>7000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC0136-7000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC0136-7000ug-rep1v2b</td>
      <td>0.062973</td>
      <td>6.297330</td>
    </tr>
    <tr>
      <th>13</th>
      <td>IDC0136</td>
      <td>5500ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC0136-5500ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC0136-5500ug-rep1v2b</td>
      <td>0.136001</td>
      <td>13.600096</td>
    </tr>
    <tr>
      <th>23</th>
      <td>IDC0337</td>
      <td>6000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC0337-6000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC0337-6000ug-rep1v2b</td>
      <td>0.059010</td>
      <td>5.900967</td>
    </tr>
    <tr>
      <th>17</th>
      <td>IDC0337</td>
      <td>4500ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC0337-4500ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC0337-4500ug-rep1v2b</td>
      <td>0.091161</td>
      <td>9.116124</td>
    </tr>
    <tr>
      <th>16</th>
      <td>IDC0337</td>
      <td>4000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC0337-4000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC0337-4000ug-rep1v2b</td>
      <td>0.102009</td>
      <td>10.200934</td>
    </tr>
    <tr>
      <th>20</th>
      <td>IDC208</td>
      <td>6000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC208-6000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC208-6000ug-rep1v2b</td>
      <td>0.068419</td>
      <td>6.841890</td>
    </tr>
    <tr>
      <th>21</th>
      <td>IDC208</td>
      <td>5000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC208-5000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC208-5000ug-rep1v2b</td>
      <td>0.103581</td>
      <td>10.358083</td>
    </tr>
    <tr>
      <th>22</th>
      <td>IDC208</td>
      <td>4000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC208-4000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC208-4000ug-rep1v2b</td>
      <td>0.150962</td>
      <td>15.096161</td>
    </tr>
    <tr>
      <th>29</th>
      <td>IDC403</td>
      <td>8500ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC403-8500ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC403-8500ug-rep1v2b</td>
      <td>0.000359</td>
      <td>0.035906</td>
    </tr>
    <tr>
      <th>18</th>
      <td>IDC403</td>
      <td>7500ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC403-7500ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC403-7500ug-rep1v2b</td>
      <td>0.000893</td>
      <td>0.089306</td>
    </tr>
    <tr>
      <th>19</th>
      <td>IDC403</td>
      <td>6000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC403-6000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC403-6000ug-rep1v2b</td>
      <td>0.001969</td>
      <td>0.196913</td>
    </tr>
    <tr>
      <th>10</th>
      <td>IDC561</td>
      <td>3500ug-rep1v2a</td>
      <td>mtvir-rep1v2a-561-IgG-3500ug</td>
      <td>mtvir-rep1v2a</td>
      <td>wt-DNA-rep1</td>
      <td>IDC561-3500ug-rep1v2a</td>
      <td>0.001820</td>
      <td>0.181995</td>
    </tr>
    <tr>
      <th>9</th>
      <td>IDC561</td>
      <td>3000ug-rep1v2a</td>
      <td>mtvir-rep1v2a-561-IgG-3000ug</td>
      <td>mtvir-rep1v2a</td>
      <td>wt-DNA-rep1</td>
      <td>IDC561-3000ug-rep1v2a</td>
      <td>0.002478</td>
      <td>0.247840</td>
    </tr>
    <tr>
      <th>24</th>
      <td>IDC561</td>
      <td>2800ug-rep1v2b</td>
      <td>mtvir-rep1v2b-561-IgG-2800ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC561-2800ug-rep1v2b</td>
      <td>0.007367</td>
      <td>0.736652</td>
    </tr>
    <tr>
      <th>25</th>
      <td>IDC561</td>
      <td>2200ug-rep1v2b</td>
      <td>mtvir-rep1v2b-561-IgG-2200ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC561-2200ug-rep1v2b</td>
      <td>0.012112</td>
      <td>1.211230</td>
    </tr>
    <tr>
      <th>2</th>
      <td>IDC561</td>
      <td>2500ug-rep3e</td>
      <td>mtvir-rep3e-561-IgG-2500ug</td>
      <td>mtvir-rep3e</td>
      <td>wt-DNA-rep3</td>
      <td>IDC561-2500ug-rep3e</td>
      <td>0.023230</td>
      <td>2.322992</td>
    </tr>
    <tr>
      <th>32</th>
      <td>QA013-2</td>
      <td>4ug-rep3</td>
      <td>mtvir-rep3-QA013-2-4ug</td>
      <td>mtvir-rep3</td>
      <td>wt-DNA-rep3</td>
      <td>QA013-2-4ug-rep3</td>
      <td>0.001661</td>
      <td>0.166100</td>
    </tr>
    <tr>
      <th>31</th>
      <td>QA013-2</td>
      <td>4ug-rep2</td>
      <td>mtvir-rep2-QA013-2-4ug</td>
      <td>mtvir-rep2</td>
      <td>wt-DNA-rep2</td>
      <td>QA013-2-4ug-rep2</td>
      <td>0.002637</td>
      <td>0.263700</td>
    </tr>
    <tr>
      <th>30</th>
      <td>QA013-2</td>
      <td>2ug-rep2</td>
      <td>mtvir-rep2-QA013-2-2ug</td>
      <td>mtvir-rep2</td>
      <td>wt-DNA-rep2</td>
      <td>QA013-2-2ug-rep2</td>
      <td>0.032488</td>
      <td>3.248800</td>
    </tr>
    <tr>
      <th>7</th>
      <td>QA013-2282dpi</td>
      <td>d10-rep1v2a</td>
      <td>mtvir-rep1v2a-QA013-2282dpi-d10</td>
      <td>mtvir-rep1v2a</td>
      <td>wt-DNA-rep1</td>
      <td>QA013-2282dpi-d10-rep1v2a</td>
      <td>0.000317</td>
      <td>0.031665</td>
    </tr>
    <tr>
      <th>6</th>
      <td>QA013-2282dpi</td>
      <td>d7-5-rep3e</td>
      <td>mtvir-rep3e-QA013-2282dpi-d7-5</td>
      <td>mtvir-rep3e</td>
      <td>wt-DNA-rep3</td>
      <td>QA013-2282dpi-d7-5-rep3e</td>
      <td>0.000388</td>
      <td>0.038845</td>
    </tr>
    <tr>
      <th>8</th>
      <td>QA013-2282dpi</td>
      <td>d15-rep1v2a</td>
      <td>mtvir-rep1v2a-QA013-2282dpi-d15</td>
      <td>mtvir-rep1v2a</td>
      <td>wt-DNA-rep1</td>
      <td>QA013-2282dpi-d15-rep1v2a</td>
      <td>0.001444</td>
      <td>0.144448</td>
    </tr>
    <tr>
      <th>28</th>
      <td>QA013-2282dpi</td>
      <td>d25-rep1v2b</td>
      <td>mtvir-rep1v2b-QA013-2282dpi-d25</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>QA013-2282dpi-d25-rep1v2b</td>
      <td>0.278701</td>
      <td>27.870066</td>
    </tr>
    <tr>
      <th>27</th>
      <td>QA013-2282dpi</td>
      <td>d40-rep1v2b</td>
      <td>mtvir-rep1v2b-QA013-2282dpi-d40</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>QA013-2282dpi-d40-rep1v2b</td>
      <td>0.470669</td>
      <td>47.066878</td>
    </tr>
  </tbody>
</table>
</div>




```python
selfilecols = []
for selfile in ['mutdiffsel', 'sitediffsel']:
    selfilecol = selfile + '_file'
    selfilecols.append(selfilecol)
    selections[selfilecol] = (diffseldir + '/' + selections['mds_names'] + '_' +
                              selfile + '.csv')
    #print(selections['serum'])
    assert all(selections[selfilecol].map(os.path.isfile)), 'missing files'
    print(f"Created {len(selections[selfilecol])} {selfile} files, adding to "
          f"`selections` data frame in column {selfilecol}")
```

    Created 33 mutdiffsel files, adding to `selections` data frame in column mutdiffsel_file
    Created 33 sitediffsel files, adding to `selections` data frame in column sitediffsel_file



```python
selections
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>group</th>
      <th>name</th>
      <th>sel</th>
      <th>mock</th>
      <th>err</th>
      <th>mds_names</th>
      <th>libfracsurvive</th>
      <th>percent_infectivity</th>
      <th>mutdiffsel_file</th>
      <th>sitediffsel_file</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1</th>
      <td>118</td>
      <td>8ug-rep2d</td>
      <td>mtvir-rep2d-118-8ug</td>
      <td>mtvir-rep2d</td>
      <td>wt-DNA-rep2</td>
      <td>118-8ug-rep2d</td>
      <td>0.000168</td>
      <td>0.016777</td>
      <td>./results/diffsel/118-8ug-rep2d_mutdiffsel.csv</td>
      <td>./results/diffsel/118-8ug-rep2d_sitediffsel.csv</td>
    </tr>
    <tr>
      <th>4</th>
      <td>118</td>
      <td>8ug-rep3d</td>
      <td>mtvir-rep3d-118-8ug</td>
      <td>mtvir-rep3d</td>
      <td>wt-DNA-rep3</td>
      <td>118-8ug-rep3d</td>
      <td>0.000673</td>
      <td>0.067347</td>
      <td>./results/diffsel/118-8ug-rep3d_mutdiffsel.csv</td>
      <td>./results/diffsel/118-8ug-rep3d_sitediffsel.csv</td>
    </tr>
    <tr>
      <th>0</th>
      <td>118</td>
      <td>4ug-rep2d</td>
      <td>mtvir-rep2d-118-4ug</td>
      <td>mtvir-rep2d</td>
      <td>wt-DNA-rep2</td>
      <td>118-4ug-rep2d</td>
      <td>0.005177</td>
      <td>0.517724</td>
      <td>./results/diffsel/118-4ug-rep2d_mutdiffsel.csv</td>
      <td>./results/diffsel/118-4ug-rep2d_sitediffsel.csv</td>
    </tr>
    <tr>
      <th>3</th>
      <td>118</td>
      <td>4ug-rep3d</td>
      <td>mtvir-rep3d-118-4ug</td>
      <td>mtvir-rep3d</td>
      <td>wt-DNA-rep3</td>
      <td>118-4ug-rep3d</td>
      <td>0.005182</td>
      <td>0.518168</td>
      <td>./results/diffsel/118-4ug-rep3d_mutdiffsel.csv</td>
      <td>./results/diffsel/118-4ug-rep3d_sitediffsel.csv</td>
    </tr>
    <tr>
      <th>11</th>
      <td>212</td>
      <td>10ug-rep1v2a</td>
      <td>mtvir-rep1v2a-212-10ug</td>
      <td>mtvir-rep1v2a</td>
      <td>wt-DNA-rep1</td>
      <td>212-10ug-rep1v2a</td>
      <td>0.002501</td>
      <td>0.250140</td>
      <td>./results/diffsel/212-10ug-rep1v2a_mutdiffsel.csv</td>
      <td>./results/diffsel/212-10ug-rep1v2a_sitediffsel...</td>
    </tr>
    <tr>
      <th>5</th>
      <td>212</td>
      <td>8ug-rep3e</td>
      <td>mtvir-rep3e-212-8ug</td>
      <td>mtvir-rep3e</td>
      <td>wt-DNA-rep3</td>
      <td>212-8ug-rep3e</td>
      <td>0.006215</td>
      <td>0.621478</td>
      <td>./results/diffsel/212-8ug-rep3e_mutdiffsel.csv</td>
      <td>./results/diffsel/212-8ug-rep3e_sitediffsel.csv</td>
    </tr>
    <tr>
      <th>12</th>
      <td>212</td>
      <td>7ug-rep1v2a</td>
      <td>mtvir-rep1v2a-212-7ug</td>
      <td>mtvir-rep1v2a</td>
      <td>wt-DNA-rep1</td>
      <td>212-7ug-rep1v2a</td>
      <td>0.007341</td>
      <td>0.734148</td>
      <td>./results/diffsel/212-7ug-rep1v2a_mutdiffsel.csv</td>
      <td>./results/diffsel/212-7ug-rep1v2a_sitediffsel.csv</td>
    </tr>
    <tr>
      <th>26</th>
      <td>212</td>
      <td>6ug-rep1v2b</td>
      <td>mtvir-rep1v2b-212-6ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>212-6ug-rep1v2b</td>
      <td>0.034943</td>
      <td>3.494317</td>
      <td>./results/diffsel/212-6ug-rep1v2b_mutdiffsel.csv</td>
      <td>./results/diffsel/212-6ug-rep1v2b_sitediffsel.csv</td>
    </tr>
    <tr>
      <th>15</th>
      <td>IDC0136</td>
      <td>8500ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC0136-8500ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC0136-8500ug-rep1v2b</td>
      <td>0.041266</td>
      <td>4.126587</td>
      <td>./results/diffsel/IDC0136-8500ug-rep1v2b_mutdi...</td>
      <td>./results/diffsel/IDC0136-8500ug-rep1v2b_sited...</td>
    </tr>
    <tr>
      <th>14</th>
      <td>IDC0136</td>
      <td>7000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC0136-7000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC0136-7000ug-rep1v2b</td>
      <td>0.062973</td>
      <td>6.297330</td>
      <td>./results/diffsel/IDC0136-7000ug-rep1v2b_mutdi...</td>
      <td>./results/diffsel/IDC0136-7000ug-rep1v2b_sited...</td>
    </tr>
    <tr>
      <th>13</th>
      <td>IDC0136</td>
      <td>5500ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC0136-5500ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC0136-5500ug-rep1v2b</td>
      <td>0.136001</td>
      <td>13.600096</td>
      <td>./results/diffsel/IDC0136-5500ug-rep1v2b_mutdi...</td>
      <td>./results/diffsel/IDC0136-5500ug-rep1v2b_sited...</td>
    </tr>
    <tr>
      <th>23</th>
      <td>IDC0337</td>
      <td>6000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC0337-6000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC0337-6000ug-rep1v2b</td>
      <td>0.059010</td>
      <td>5.900967</td>
      <td>./results/diffsel/IDC0337-6000ug-rep1v2b_mutdi...</td>
      <td>./results/diffsel/IDC0337-6000ug-rep1v2b_sited...</td>
    </tr>
    <tr>
      <th>17</th>
      <td>IDC0337</td>
      <td>4500ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC0337-4500ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC0337-4500ug-rep1v2b</td>
      <td>0.091161</td>
      <td>9.116124</td>
      <td>./results/diffsel/IDC0337-4500ug-rep1v2b_mutdi...</td>
      <td>./results/diffsel/IDC0337-4500ug-rep1v2b_sited...</td>
    </tr>
    <tr>
      <th>16</th>
      <td>IDC0337</td>
      <td>4000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC0337-4000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC0337-4000ug-rep1v2b</td>
      <td>0.102009</td>
      <td>10.200934</td>
      <td>./results/diffsel/IDC0337-4000ug-rep1v2b_mutdi...</td>
      <td>./results/diffsel/IDC0337-4000ug-rep1v2b_sited...</td>
    </tr>
    <tr>
      <th>20</th>
      <td>IDC208</td>
      <td>6000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC208-6000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC208-6000ug-rep1v2b</td>
      <td>0.068419</td>
      <td>6.841890</td>
      <td>./results/diffsel/IDC208-6000ug-rep1v2b_mutdif...</td>
      <td>./results/diffsel/IDC208-6000ug-rep1v2b_sitedi...</td>
    </tr>
    <tr>
      <th>21</th>
      <td>IDC208</td>
      <td>5000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC208-5000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC208-5000ug-rep1v2b</td>
      <td>0.103581</td>
      <td>10.358083</td>
      <td>./results/diffsel/IDC208-5000ug-rep1v2b_mutdif...</td>
      <td>./results/diffsel/IDC208-5000ug-rep1v2b_sitedi...</td>
    </tr>
    <tr>
      <th>22</th>
      <td>IDC208</td>
      <td>4000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC208-4000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC208-4000ug-rep1v2b</td>
      <td>0.150962</td>
      <td>15.096161</td>
      <td>./results/diffsel/IDC208-4000ug-rep1v2b_mutdif...</td>
      <td>./results/diffsel/IDC208-4000ug-rep1v2b_sitedi...</td>
    </tr>
    <tr>
      <th>29</th>
      <td>IDC403</td>
      <td>8500ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC403-8500ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC403-8500ug-rep1v2b</td>
      <td>0.000359</td>
      <td>0.035906</td>
      <td>./results/diffsel/IDC403-8500ug-rep1v2b_mutdif...</td>
      <td>./results/diffsel/IDC403-8500ug-rep1v2b_sitedi...</td>
    </tr>
    <tr>
      <th>18</th>
      <td>IDC403</td>
      <td>7500ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC403-7500ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC403-7500ug-rep1v2b</td>
      <td>0.000893</td>
      <td>0.089306</td>
      <td>./results/diffsel/IDC403-7500ug-rep1v2b_mutdif...</td>
      <td>./results/diffsel/IDC403-7500ug-rep1v2b_sitedi...</td>
    </tr>
    <tr>
      <th>19</th>
      <td>IDC403</td>
      <td>6000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC403-6000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC403-6000ug-rep1v2b</td>
      <td>0.001969</td>
      <td>0.196913</td>
      <td>./results/diffsel/IDC403-6000ug-rep1v2b_mutdif...</td>
      <td>./results/diffsel/IDC403-6000ug-rep1v2b_sitedi...</td>
    </tr>
    <tr>
      <th>10</th>
      <td>IDC561</td>
      <td>3500ug-rep1v2a</td>
      <td>mtvir-rep1v2a-561-IgG-3500ug</td>
      <td>mtvir-rep1v2a</td>
      <td>wt-DNA-rep1</td>
      <td>IDC561-3500ug-rep1v2a</td>
      <td>0.001820</td>
      <td>0.181995</td>
      <td>./results/diffsel/IDC561-3500ug-rep1v2a_mutdif...</td>
      <td>./results/diffsel/IDC561-3500ug-rep1v2a_sitedi...</td>
    </tr>
    <tr>
      <th>9</th>
      <td>IDC561</td>
      <td>3000ug-rep1v2a</td>
      <td>mtvir-rep1v2a-561-IgG-3000ug</td>
      <td>mtvir-rep1v2a</td>
      <td>wt-DNA-rep1</td>
      <td>IDC561-3000ug-rep1v2a</td>
      <td>0.002478</td>
      <td>0.247840</td>
      <td>./results/diffsel/IDC561-3000ug-rep1v2a_mutdif...</td>
      <td>./results/diffsel/IDC561-3000ug-rep1v2a_sitedi...</td>
    </tr>
    <tr>
      <th>24</th>
      <td>IDC561</td>
      <td>2800ug-rep1v2b</td>
      <td>mtvir-rep1v2b-561-IgG-2800ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC561-2800ug-rep1v2b</td>
      <td>0.007367</td>
      <td>0.736652</td>
      <td>./results/diffsel/IDC561-2800ug-rep1v2b_mutdif...</td>
      <td>./results/diffsel/IDC561-2800ug-rep1v2b_sitedi...</td>
    </tr>
    <tr>
      <th>25</th>
      <td>IDC561</td>
      <td>2200ug-rep1v2b</td>
      <td>mtvir-rep1v2b-561-IgG-2200ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC561-2200ug-rep1v2b</td>
      <td>0.012112</td>
      <td>1.211230</td>
      <td>./results/diffsel/IDC561-2200ug-rep1v2b_mutdif...</td>
      <td>./results/diffsel/IDC561-2200ug-rep1v2b_sitedi...</td>
    </tr>
    <tr>
      <th>2</th>
      <td>IDC561</td>
      <td>2500ug-rep3e</td>
      <td>mtvir-rep3e-561-IgG-2500ug</td>
      <td>mtvir-rep3e</td>
      <td>wt-DNA-rep3</td>
      <td>IDC561-2500ug-rep3e</td>
      <td>0.023230</td>
      <td>2.322992</td>
      <td>./results/diffsel/IDC561-2500ug-rep3e_mutdiffs...</td>
      <td>./results/diffsel/IDC561-2500ug-rep3e_sitediff...</td>
    </tr>
    <tr>
      <th>32</th>
      <td>QA013-2</td>
      <td>4ug-rep3</td>
      <td>mtvir-rep3-QA013-2-4ug</td>
      <td>mtvir-rep3</td>
      <td>wt-DNA-rep3</td>
      <td>QA013-2-4ug-rep3</td>
      <td>0.001661</td>
      <td>0.166100</td>
      <td>./results/diffsel/QA013-2-4ug-rep3_mutdiffsel.csv</td>
      <td>./results/diffsel/QA013-2-4ug-rep3_sitediffsel...</td>
    </tr>
    <tr>
      <th>31</th>
      <td>QA013-2</td>
      <td>4ug-rep2</td>
      <td>mtvir-rep2-QA013-2-4ug</td>
      <td>mtvir-rep2</td>
      <td>wt-DNA-rep2</td>
      <td>QA013-2-4ug-rep2</td>
      <td>0.002637</td>
      <td>0.263700</td>
      <td>./results/diffsel/QA013-2-4ug-rep2_mutdiffsel.csv</td>
      <td>./results/diffsel/QA013-2-4ug-rep2_sitediffsel...</td>
    </tr>
    <tr>
      <th>30</th>
      <td>QA013-2</td>
      <td>2ug-rep2</td>
      <td>mtvir-rep2-QA013-2-2ug</td>
      <td>mtvir-rep2</td>
      <td>wt-DNA-rep2</td>
      <td>QA013-2-2ug-rep2</td>
      <td>0.032488</td>
      <td>3.248800</td>
      <td>./results/diffsel/QA013-2-2ug-rep2_mutdiffsel.csv</td>
      <td>./results/diffsel/QA013-2-2ug-rep2_sitediffsel...</td>
    </tr>
    <tr>
      <th>7</th>
      <td>QA013-2282dpi</td>
      <td>d10-rep1v2a</td>
      <td>mtvir-rep1v2a-QA013-2282dpi-d10</td>
      <td>mtvir-rep1v2a</td>
      <td>wt-DNA-rep1</td>
      <td>QA013-2282dpi-d10-rep1v2a</td>
      <td>0.000317</td>
      <td>0.031665</td>
      <td>./results/diffsel/QA013-2282dpi-d10-rep1v2a_mu...</td>
      <td>./results/diffsel/QA013-2282dpi-d10-rep1v2a_si...</td>
    </tr>
    <tr>
      <th>6</th>
      <td>QA013-2282dpi</td>
      <td>d7-5-rep3e</td>
      <td>mtvir-rep3e-QA013-2282dpi-d7-5</td>
      <td>mtvir-rep3e</td>
      <td>wt-DNA-rep3</td>
      <td>QA013-2282dpi-d7-5-rep3e</td>
      <td>0.000388</td>
      <td>0.038845</td>
      <td>./results/diffsel/QA013-2282dpi-d7-5-rep3e_mut...</td>
      <td>./results/diffsel/QA013-2282dpi-d7-5-rep3e_sit...</td>
    </tr>
    <tr>
      <th>8</th>
      <td>QA013-2282dpi</td>
      <td>d15-rep1v2a</td>
      <td>mtvir-rep1v2a-QA013-2282dpi-d15</td>
      <td>mtvir-rep1v2a</td>
      <td>wt-DNA-rep1</td>
      <td>QA013-2282dpi-d15-rep1v2a</td>
      <td>0.001444</td>
      <td>0.144448</td>
      <td>./results/diffsel/QA013-2282dpi-d15-rep1v2a_mu...</td>
      <td>./results/diffsel/QA013-2282dpi-d15-rep1v2a_si...</td>
    </tr>
    <tr>
      <th>28</th>
      <td>QA013-2282dpi</td>
      <td>d25-rep1v2b</td>
      <td>mtvir-rep1v2b-QA013-2282dpi-d25</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>QA013-2282dpi-d25-rep1v2b</td>
      <td>0.278701</td>
      <td>27.870066</td>
      <td>./results/diffsel/QA013-2282dpi-d25-rep1v2b_mu...</td>
      <td>./results/diffsel/QA013-2282dpi-d25-rep1v2b_si...</td>
    </tr>
    <tr>
      <th>27</th>
      <td>QA013-2282dpi</td>
      <td>d40-rep1v2b</td>
      <td>mtvir-rep1v2b-QA013-2282dpi-d40</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>QA013-2282dpi-d40-rep1v2b</td>
      <td>0.470669</td>
      <td>47.066878</td>
      <td>./results/diffsel/QA013-2282dpi-d40-rep1v2b_mu...</td>
      <td>./results/diffsel/QA013-2282dpi-d40-rep1v2b_si...</td>
    </tr>
  </tbody>
</table>
</div>



### Get all selection information in one data frame

For further processing, we want to create a dataframe that holds all of the selection information at the site and mutation levels for all samples. We create such a dataframe, sel_df, by reading the files in selections into the data frame using dms_tools2.diffsel.df_read_filecols:



```python
sel_df = (dms_tools2.diffsel.df_read_filecols(selections, selfilecols)
          .drop(columns=selfilecols)
          )
```

Now sel_df is a very large data frame, but it has all the information we want to plot. Here are the first few rows:


```python
print(f"sel_df has {len(sel_df)} rows. Here are the last few:")
display(HTML(sel_df.tail(n=5).to_html(index=False)))
```

    sel_df has 442200 rows. Here are the last few:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>group</th>
      <th>name</th>
      <th>sel</th>
      <th>mock</th>
      <th>err</th>
      <th>mds_names</th>
      <th>libfracsurvive</th>
      <th>percent_infectivity</th>
      <th>site</th>
      <th>wildtype</th>
      <th>mutation</th>
      <th>mutdiffsel</th>
      <th>abs_diffsel</th>
      <th>positive_diffsel</th>
      <th>negative_diffsel</th>
      <th>max_diffsel</th>
      <th>min_diffsel</th>
      <th>isite</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>QA013-2282dpi</td>
      <td>d40-rep1v2b</td>
      <td>mtvir-rep1v2b-QA013-2282dpi-d40</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>QA013-2282dpi-d40-rep1v2b</td>
      <td>0.470669</td>
      <td>47.0669</td>
      <td>533</td>
      <td>A</td>
      <td>I</td>
      <td>-0.000345</td>
      <td>1.281547</td>
      <td>0.0</td>
      <td>-1.281547</td>
      <td>0.0</td>
      <td>-0.624449</td>
      <td>500</td>
    </tr>
    <tr>
      <td>QA013-2282dpi</td>
      <td>d40-rep1v2b</td>
      <td>mtvir-rep1v2b-QA013-2282dpi-d40</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>QA013-2282dpi-d40-rep1v2b</td>
      <td>0.470669</td>
      <td>47.0669</td>
      <td>533</td>
      <td>A</td>
      <td>T</td>
      <td>-0.213414</td>
      <td>1.281547</td>
      <td>0.0</td>
      <td>-1.281547</td>
      <td>0.0</td>
      <td>-0.624449</td>
      <td>500</td>
    </tr>
    <tr>
      <td>QA013-2282dpi</td>
      <td>d40-rep1v2b</td>
      <td>mtvir-rep1v2b-QA013-2282dpi-d40</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>QA013-2282dpi-d40-rep1v2b</td>
      <td>0.470669</td>
      <td>47.0669</td>
      <td>533</td>
      <td>A</td>
      <td>G</td>
      <td>-0.438162</td>
      <td>1.281547</td>
      <td>0.0</td>
      <td>-1.281547</td>
      <td>0.0</td>
      <td>-0.624449</td>
      <td>500</td>
    </tr>
    <tr>
      <td>QA013-2282dpi</td>
      <td>d40-rep1v2b</td>
      <td>mtvir-rep1v2b-QA013-2282dpi-d40</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>QA013-2282dpi-d40-rep1v2b</td>
      <td>0.470669</td>
      <td>47.0669</td>
      <td>533</td>
      <td>A</td>
      <td>V</td>
      <td>-0.624449</td>
      <td>1.281547</td>
      <td>0.0</td>
      <td>-1.281547</td>
      <td>0.0</td>
      <td>-0.624449</td>
      <td>500</td>
    </tr>
    <tr>
      <td>QA013-2282dpi</td>
      <td>d40-rep1v2b</td>
      <td>mtvir-rep1v2b-QA013-2282dpi-d40</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>QA013-2282dpi-d40-rep1v2b</td>
      <td>0.470669</td>
      <td>47.0669</td>
      <td>533</td>
      <td>A</td>
      <td>A</td>
      <td>NaN</td>
      <td>1.281547</td>
      <td>0.0</td>
      <td>-1.281547</td>
      <td>0.0</td>
      <td>-0.624449</td>
      <td>500</td>
    </tr>
  </tbody>
</table>


#plot for all samples 
STOP
sel_df = sel_df.assign(    # add informative names for serum and samples
        serum_name_formatted=lambda x: x['serum_name'],
        name_formatted=lambda x:
            x['library'] + ', ' + x['percent_infectivity'].apply(
                dms_tools2.utils.sigFigStr, nsig=2) + '% infectivity, ' + "1:" + x["serum_dilution"].apply(
                dms_tools2.utils.sigFigStr, nsig=1) + " serum dilution",  
        name=lambda x:
            x['library'] + '-' + x['percent_infectivity'].apply(
                dms_tools2.utils.sigFigStr, nsig=2)
        )


```python
#sel_df["percent_infectivity_rnd"] = sel_df["percent_infectivity"].round(2)
sel_df["percent_infectivity_rnd"] = sel_df["percent_infectivity"].apply(lambda x: round(x, 2))
#sel_df
```


```python

```


```python
#STOP
sel_df["name_formatted"] = (sel_df["mds_names"].astype(str) + " % infectivity:"+ sel_df["percent_infectivity_rnd"].astype(str))
#display(HTML(sel_df.head(n=5).to_html(index=False)))
```


```python

import seaborn as sns
cm = sns.light_palette("green", as_cmap=True)
for serum_name, sub_sel_df in sel_df.groupby('group'):

    print(f"\n\n******************* {serum_name} *******************")

    fig, ax = dmslogo.facet_plot(
                sub_sel_df,
                x_col='isite',
                gridcol_col='name_formatted',
                show_col=None,
                wspace=0.6,
                draw_line_kwargs=dict(
                        xtick_col='site',
                        height_col='positive_diffsel',
                        ylabel='positive_diffsel',
                        )
                )
    display(fig)
    plt.close(fig)

    corr_df = (sub_sel_df
               #.rename(columns={'name_formatted': 'sample'})
               .pivot_table(values='positive_diffsel',
                            columns='name_formatted',
                            index=['site'])
               .corr()
               .round(3)
               )
    display(corr_df.style.background_gradient(cmap=cm, low=0, high=1))

```

    
    
    ******************* 118 *******************



![png](analysis_notebook_files/analysis_notebook_54_1.png)



<style  type="text/css" >
    #T_50b4d416_5e61_11ea_aee3_0cc47abc358drow0_col0 {
            background-color:  #72bf72;
            color:  #000000;
        }    #T_50b4d416_5e61_11ea_aee3_0cc47abc358drow0_col1 {
            background-color:  #d0f3d0;
            color:  #000000;
        }    #T_50b4d416_5e61_11ea_aee3_0cc47abc358drow0_col2 {
            background-color:  #d5f6d5;
            color:  #000000;
        }    #T_50b4d416_5e61_11ea_aee3_0cc47abc358drow0_col3 {
            background-color:  #ddfbdd;
            color:  #000000;
        }    #T_50b4d416_5e61_11ea_aee3_0cc47abc358drow1_col0 {
            background-color:  #d3f5d3;
            color:  #000000;
        }    #T_50b4d416_5e61_11ea_aee3_0cc47abc358drow1_col1 {
            background-color:  #72bf72;
            color:  #000000;
        }    #T_50b4d416_5e61_11ea_aee3_0cc47abc358drow1_col2 {
            background-color:  #d8f8d8;
            color:  #000000;
        }    #T_50b4d416_5e61_11ea_aee3_0cc47abc358drow1_col3 {
            background-color:  #dffcdf;
            color:  #000000;
        }    #T_50b4d416_5e61_11ea_aee3_0cc47abc358drow2_col0 {
            background-color:  #ddfbdd;
            color:  #000000;
        }    #T_50b4d416_5e61_11ea_aee3_0cc47abc358drow2_col1 {
            background-color:  #defbde;
            color:  #000000;
        }    #T_50b4d416_5e61_11ea_aee3_0cc47abc358drow2_col2 {
            background-color:  #72bf72;
            color:  #000000;
        }    #T_50b4d416_5e61_11ea_aee3_0cc47abc358drow2_col3 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_50b4d416_5e61_11ea_aee3_0cc47abc358drow3_col0 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_50b4d416_5e61_11ea_aee3_0cc47abc358drow3_col1 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_50b4d416_5e61_11ea_aee3_0cc47abc358drow3_col2 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_50b4d416_5e61_11ea_aee3_0cc47abc358drow3_col3 {
            background-color:  #72bf72;
            color:  #000000;
        }</style><table id="T_50b4d416_5e61_11ea_aee3_0cc47abc358d" ><thead>    <tr>        <th class="index_name level0" >name_formatted</th>        <th class="col_heading level0 col0" >118-4ug-rep2d % infectivity:0.52</th>        <th class="col_heading level0 col1" >118-4ug-rep3d % infectivity:0.52</th>        <th class="col_heading level0 col2" >118-8ug-rep2d % infectivity:0.02</th>        <th class="col_heading level0 col3" >118-8ug-rep3d % infectivity:0.07</th>    </tr>    <tr>        <th class="index_name level0" >name_formatted</th>        <th class="blank" ></th>        <th class="blank" ></th>        <th class="blank" ></th>        <th class="blank" ></th>    </tr></thead><tbody>
                <tr>
                        <th id="T_50b4d416_5e61_11ea_aee3_0cc47abc358dlevel0_row0" class="row_heading level0 row0" >118-4ug-rep2d % infectivity:0.52</th>
                        <td id="T_50b4d416_5e61_11ea_aee3_0cc47abc358drow0_col0" class="data row0 col0" >1</td>
                        <td id="T_50b4d416_5e61_11ea_aee3_0cc47abc358drow0_col1" class="data row0 col1" >0.617</td>
                        <td id="T_50b4d416_5e61_11ea_aee3_0cc47abc358drow0_col2" class="data row0 col2" >0.572</td>
                        <td id="T_50b4d416_5e61_11ea_aee3_0cc47abc358drow0_col3" class="data row0 col3" >0.538</td>
            </tr>
            <tr>
                        <th id="T_50b4d416_5e61_11ea_aee3_0cc47abc358dlevel0_row1" class="row_heading level0 row1" >118-4ug-rep3d % infectivity:0.52</th>
                        <td id="T_50b4d416_5e61_11ea_aee3_0cc47abc358drow1_col0" class="data row1 col0" >0.617</td>
                        <td id="T_50b4d416_5e61_11ea_aee3_0cc47abc358drow1_col1" class="data row1 col1" >1</td>
                        <td id="T_50b4d416_5e61_11ea_aee3_0cc47abc358drow1_col2" class="data row1 col2" >0.558</td>
                        <td id="T_50b4d416_5e61_11ea_aee3_0cc47abc358drow1_col3" class="data row1 col3" >0.528</td>
            </tr>
            <tr>
                        <th id="T_50b4d416_5e61_11ea_aee3_0cc47abc358dlevel0_row2" class="row_heading level0 row2" >118-8ug-rep2d % infectivity:0.02</th>
                        <td id="T_50b4d416_5e61_11ea_aee3_0cc47abc358drow2_col0" class="data row2 col0" >0.572</td>
                        <td id="T_50b4d416_5e61_11ea_aee3_0cc47abc358drow2_col1" class="data row2 col1" >0.558</td>
                        <td id="T_50b4d416_5e61_11ea_aee3_0cc47abc358drow2_col2" class="data row2 col2" >1</td>
                        <td id="T_50b4d416_5e61_11ea_aee3_0cc47abc358drow2_col3" class="data row2 col3" >0.499</td>
            </tr>
            <tr>
                        <th id="T_50b4d416_5e61_11ea_aee3_0cc47abc358dlevel0_row3" class="row_heading level0 row3" >118-8ug-rep3d % infectivity:0.07</th>
                        <td id="T_50b4d416_5e61_11ea_aee3_0cc47abc358drow3_col0" class="data row3 col0" >0.538</td>
                        <td id="T_50b4d416_5e61_11ea_aee3_0cc47abc358drow3_col1" class="data row3 col1" >0.528</td>
                        <td id="T_50b4d416_5e61_11ea_aee3_0cc47abc358drow3_col2" class="data row3 col2" >0.499</td>
                        <td id="T_50b4d416_5e61_11ea_aee3_0cc47abc358drow3_col3" class="data row3 col3" >1</td>
            </tr>
    </tbody></table>


    
    
    ******************* 212 *******************



![png](analysis_notebook_files/analysis_notebook_54_4.png)



<style  type="text/css" >
    #T_5362899c_5e61_11ea_aee3_0cc47abc358drow0_col0 {
            background-color:  #73c073;
            color:  #000000;
        }    #T_5362899c_5e61_11ea_aee3_0cc47abc358drow0_col1 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_5362899c_5e61_11ea_aee3_0cc47abc358drow0_col2 {
            background-color:  #d4f6d4;
            color:  #000000;
        }    #T_5362899c_5e61_11ea_aee3_0cc47abc358drow0_col3 {
            background-color:  #e3fee3;
            color:  #000000;
        }    #T_5362899c_5e61_11ea_aee3_0cc47abc358drow1_col0 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_5362899c_5e61_11ea_aee3_0cc47abc358drow1_col1 {
            background-color:  #73c073;
            color:  #000000;
        }    #T_5362899c_5e61_11ea_aee3_0cc47abc358drow1_col2 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_5362899c_5e61_11ea_aee3_0cc47abc358drow1_col3 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_5362899c_5e61_11ea_aee3_0cc47abc358drow2_col0 {
            background-color:  #c6eec6;
            color:  #000000;
        }    #T_5362899c_5e61_11ea_aee3_0cc47abc358drow2_col1 {
            background-color:  #d5f6d5;
            color:  #000000;
        }    #T_5362899c_5e61_11ea_aee3_0cc47abc358drow2_col2 {
            background-color:  #72bf72;
            color:  #000000;
        }    #T_5362899c_5e61_11ea_aee3_0cc47abc358drow2_col3 {
            background-color:  #d4f6d4;
            color:  #000000;
        }    #T_5362899c_5e61_11ea_aee3_0cc47abc358drow3_col0 {
            background-color:  #d8f8d8;
            color:  #000000;
        }    #T_5362899c_5e61_11ea_aee3_0cc47abc358drow3_col1 {
            background-color:  #dbf9db;
            color:  #000000;
        }    #T_5362899c_5e61_11ea_aee3_0cc47abc358drow3_col2 {
            background-color:  #dbf9db;
            color:  #000000;
        }    #T_5362899c_5e61_11ea_aee3_0cc47abc358drow3_col3 {
            background-color:  #72bf72;
            color:  #000000;
        }</style><table id="T_5362899c_5e61_11ea_aee3_0cc47abc358d" ><thead>    <tr>        <th class="index_name level0" >name_formatted</th>        <th class="col_heading level0 col0" >212-10ug-rep1v2a % infectivity:0.25</th>        <th class="col_heading level0 col1" >212-6ug-rep1v2b % infectivity:3.49</th>        <th class="col_heading level0 col2" >212-7ug-rep1v2a % infectivity:0.73</th>        <th class="col_heading level0 col3" >212-8ug-rep3e % infectivity:0.62</th>    </tr>    <tr>        <th class="index_name level0" >name_formatted</th>        <th class="blank" ></th>        <th class="blank" ></th>        <th class="blank" ></th>        <th class="blank" ></th>    </tr></thead><tbody>
                <tr>
                        <th id="T_5362899c_5e61_11ea_aee3_0cc47abc358dlevel0_row0" class="row_heading level0 row0" >212-10ug-rep1v2a % infectivity:0.25</th>
                        <td id="T_5362899c_5e61_11ea_aee3_0cc47abc358drow0_col0" class="data row0 col0" >1</td>
                        <td id="T_5362899c_5e61_11ea_aee3_0cc47abc358drow0_col1" class="data row0 col1" >0.372</td>
                        <td id="T_5362899c_5e61_11ea_aee3_0cc47abc358drow0_col2" class="data row0 col2" >0.548</td>
                        <td id="T_5362899c_5e61_11ea_aee3_0cc47abc358drow0_col3" class="data row0 col3" >0.449</td>
            </tr>
            <tr>
                        <th id="T_5362899c_5e61_11ea_aee3_0cc47abc358dlevel0_row1" class="row_heading level0 row1" >212-6ug-rep1v2b % infectivity:3.49</th>
                        <td id="T_5362899c_5e61_11ea_aee3_0cc47abc358drow1_col0" class="data row1 col0" >0.372</td>
                        <td id="T_5362899c_5e61_11ea_aee3_0cc47abc358drow1_col1" class="data row1 col1" >1</td>
                        <td id="T_5362899c_5e61_11ea_aee3_0cc47abc358drow1_col2" class="data row1 col2" >0.465</td>
                        <td id="T_5362899c_5e61_11ea_aee3_0cc47abc358drow1_col3" class="data row1 col3" >0.432</td>
            </tr>
            <tr>
                        <th id="T_5362899c_5e61_11ea_aee3_0cc47abc358dlevel0_row2" class="row_heading level0 row2" >212-7ug-rep1v2a % infectivity:0.73</th>
                        <td id="T_5362899c_5e61_11ea_aee3_0cc47abc358drow2_col0" class="data row2 col0" >0.548</td>
                        <td id="T_5362899c_5e61_11ea_aee3_0cc47abc358drow2_col1" class="data row2 col1" >0.465</td>
                        <td id="T_5362899c_5e61_11ea_aee3_0cc47abc358drow2_col2" class="data row2 col2" >1</td>
                        <td id="T_5362899c_5e61_11ea_aee3_0cc47abc358drow2_col3" class="data row2 col3" >0.517</td>
            </tr>
            <tr>
                        <th id="T_5362899c_5e61_11ea_aee3_0cc47abc358dlevel0_row3" class="row_heading level0 row3" >212-8ug-rep3e % infectivity:0.62</th>
                        <td id="T_5362899c_5e61_11ea_aee3_0cc47abc358drow3_col0" class="data row3 col0" >0.449</td>
                        <td id="T_5362899c_5e61_11ea_aee3_0cc47abc358drow3_col1" class="data row3 col1" >0.432</td>
                        <td id="T_5362899c_5e61_11ea_aee3_0cc47abc358drow3_col2" class="data row3 col2" >0.517</td>
                        <td id="T_5362899c_5e61_11ea_aee3_0cc47abc358drow3_col3" class="data row3 col3" >1</td>
            </tr>
    </tbody></table>


    
    
    ******************* IDC0136 *******************



![png](analysis_notebook_files/analysis_notebook_54_7.png)



<style  type="text/css" >
    #T_55df1fc8_5e61_11ea_aee3_0cc47abc358drow0_col0 {
            background-color:  #72bf72;
            color:  #000000;
        }    #T_55df1fc8_5e61_11ea_aee3_0cc47abc358drow0_col1 {
            background-color:  #dbf9db;
            color:  #000000;
        }    #T_55df1fc8_5e61_11ea_aee3_0cc47abc358drow0_col2 {
            background-color:  #dcfadc;
            color:  #000000;
        }    #T_55df1fc8_5e61_11ea_aee3_0cc47abc358drow1_col0 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_55df1fc8_5e61_11ea_aee3_0cc47abc358drow1_col1 {
            background-color:  #72bf72;
            color:  #000000;
        }    #T_55df1fc8_5e61_11ea_aee3_0cc47abc358drow1_col2 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_55df1fc8_5e61_11ea_aee3_0cc47abc358drow2_col0 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_55df1fc8_5e61_11ea_aee3_0cc47abc358drow2_col1 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_55df1fc8_5e61_11ea_aee3_0cc47abc358drow2_col2 {
            background-color:  #72bf72;
            color:  #000000;
        }</style><table id="T_55df1fc8_5e61_11ea_aee3_0cc47abc358d" ><thead>    <tr>        <th class="index_name level0" >name_formatted</th>        <th class="col_heading level0 col0" >IDC0136-5500ug-rep1v2b % infectivity:13.6</th>        <th class="col_heading level0 col1" >IDC0136-7000ug-rep1v2b % infectivity:6.3</th>        <th class="col_heading level0 col2" >IDC0136-8500ug-rep1v2b % infectivity:4.13</th>    </tr>    <tr>        <th class="index_name level0" >name_formatted</th>        <th class="blank" ></th>        <th class="blank" ></th>        <th class="blank" ></th>    </tr></thead><tbody>
                <tr>
                        <th id="T_55df1fc8_5e61_11ea_aee3_0cc47abc358dlevel0_row0" class="row_heading level0 row0" >IDC0136-5500ug-rep1v2b % infectivity:13.6</th>
                        <td id="T_55df1fc8_5e61_11ea_aee3_0cc47abc358drow0_col0" class="data row0 col0" >1</td>
                        <td id="T_55df1fc8_5e61_11ea_aee3_0cc47abc358drow0_col1" class="data row0 col1" >0.584</td>
                        <td id="T_55df1fc8_5e61_11ea_aee3_0cc47abc358drow0_col2" class="data row0 col2" >0.581</td>
            </tr>
            <tr>
                        <th id="T_55df1fc8_5e61_11ea_aee3_0cc47abc358dlevel0_row1" class="row_heading level0 row1" >IDC0136-7000ug-rep1v2b % infectivity:6.3</th>
                        <td id="T_55df1fc8_5e61_11ea_aee3_0cc47abc358drow1_col0" class="data row1 col0" >0.584</td>
                        <td id="T_55df1fc8_5e61_11ea_aee3_0cc47abc358drow1_col1" class="data row1 col1" >1</td>
                        <td id="T_55df1fc8_5e61_11ea_aee3_0cc47abc358drow1_col2" class="data row1 col2" >0.54</td>
            </tr>
            <tr>
                        <th id="T_55df1fc8_5e61_11ea_aee3_0cc47abc358dlevel0_row2" class="row_heading level0 row2" >IDC0136-8500ug-rep1v2b % infectivity:4.13</th>
                        <td id="T_55df1fc8_5e61_11ea_aee3_0cc47abc358drow2_col0" class="data row2 col0" >0.581</td>
                        <td id="T_55df1fc8_5e61_11ea_aee3_0cc47abc358drow2_col1" class="data row2 col1" >0.54</td>
                        <td id="T_55df1fc8_5e61_11ea_aee3_0cc47abc358drow2_col2" class="data row2 col2" >1</td>
            </tr>
    </tbody></table>


    
    
    ******************* IDC0337 *******************



![png](analysis_notebook_files/analysis_notebook_54_10.png)



<style  type="text/css" >
    #T_584f8838_5e61_11ea_aee3_0cc47abc358drow0_col0 {
            background-color:  #72bf72;
            color:  #000000;
        }    #T_584f8838_5e61_11ea_aee3_0cc47abc358drow0_col1 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_584f8838_5e61_11ea_aee3_0cc47abc358drow0_col2 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_584f8838_5e61_11ea_aee3_0cc47abc358drow1_col0 {
            background-color:  #dffcdf;
            color:  #000000;
        }    #T_584f8838_5e61_11ea_aee3_0cc47abc358drow1_col1 {
            background-color:  #73c073;
            color:  #000000;
        }    #T_584f8838_5e61_11ea_aee3_0cc47abc358drow1_col2 {
            background-color:  #ddfbdd;
            color:  #000000;
        }    #T_584f8838_5e61_11ea_aee3_0cc47abc358drow2_col0 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_584f8838_5e61_11ea_aee3_0cc47abc358drow2_col1 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_584f8838_5e61_11ea_aee3_0cc47abc358drow2_col2 {
            background-color:  #72bf72;
            color:  #000000;
        }</style><table id="T_584f8838_5e61_11ea_aee3_0cc47abc358d" ><thead>    <tr>        <th class="index_name level0" >name_formatted</th>        <th class="col_heading level0 col0" >IDC0337-4000ug-rep1v2b % infectivity:10.2</th>        <th class="col_heading level0 col1" >IDC0337-4500ug-rep1v2b % infectivity:9.12</th>        <th class="col_heading level0 col2" >IDC0337-6000ug-rep1v2b % infectivity:5.9</th>    </tr>    <tr>        <th class="index_name level0" >name_formatted</th>        <th class="blank" ></th>        <th class="blank" ></th>        <th class="blank" ></th>    </tr></thead><tbody>
                <tr>
                        <th id="T_584f8838_5e61_11ea_aee3_0cc47abc358dlevel0_row0" class="row_heading level0 row0" >IDC0337-4000ug-rep1v2b % infectivity:10.2</th>
                        <td id="T_584f8838_5e61_11ea_aee3_0cc47abc358drow0_col0" class="data row0 col0" >1</td>
                        <td id="T_584f8838_5e61_11ea_aee3_0cc47abc358drow0_col1" class="data row0 col1" >0.666</td>
                        <td id="T_584f8838_5e61_11ea_aee3_0cc47abc358drow0_col2" class="data row0 col2" >0.645</td>
            </tr>
            <tr>
                        <th id="T_584f8838_5e61_11ea_aee3_0cc47abc358dlevel0_row1" class="row_heading level0 row1" >IDC0337-4500ug-rep1v2b % infectivity:9.12</th>
                        <td id="T_584f8838_5e61_11ea_aee3_0cc47abc358drow1_col0" class="data row1 col0" >0.666</td>
                        <td id="T_584f8838_5e61_11ea_aee3_0cc47abc358drow1_col1" class="data row1 col1" >1</td>
                        <td id="T_584f8838_5e61_11ea_aee3_0cc47abc358drow1_col2" class="data row1 col2" >0.67</td>
            </tr>
            <tr>
                        <th id="T_584f8838_5e61_11ea_aee3_0cc47abc358dlevel0_row2" class="row_heading level0 row2" >IDC0337-6000ug-rep1v2b % infectivity:5.9</th>
                        <td id="T_584f8838_5e61_11ea_aee3_0cc47abc358drow2_col0" class="data row2 col0" >0.645</td>
                        <td id="T_584f8838_5e61_11ea_aee3_0cc47abc358drow2_col1" class="data row2 col1" >0.67</td>
                        <td id="T_584f8838_5e61_11ea_aee3_0cc47abc358drow2_col2" class="data row2 col2" >1</td>
            </tr>
    </tbody></table>


    
    
    ******************* IDC208 *******************



![png](analysis_notebook_files/analysis_notebook_54_13.png)



<style  type="text/css" >
    #T_5ab6b240_5e61_11ea_aee3_0cc47abc358drow0_col0 {
            background-color:  #72bf72;
            color:  #000000;
        }    #T_5ab6b240_5e61_11ea_aee3_0cc47abc358drow0_col1 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_5ab6b240_5e61_11ea_aee3_0cc47abc358drow0_col2 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_5ab6b240_5e61_11ea_aee3_0cc47abc358drow1_col0 {
            background-color:  #e4fee4;
            color:  #000000;
        }    #T_5ab6b240_5e61_11ea_aee3_0cc47abc358drow1_col1 {
            background-color:  #72bf72;
            color:  #000000;
        }    #T_5ab6b240_5e61_11ea_aee3_0cc47abc358drow1_col2 {
            background-color:  #cbf1cb;
            color:  #000000;
        }    #T_5ab6b240_5e61_11ea_aee3_0cc47abc358drow2_col0 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_5ab6b240_5e61_11ea_aee3_0cc47abc358drow2_col1 {
            background-color:  #cdf2cd;
            color:  #000000;
        }    #T_5ab6b240_5e61_11ea_aee3_0cc47abc358drow2_col2 {
            background-color:  #72bf72;
            color:  #000000;
        }</style><table id="T_5ab6b240_5e61_11ea_aee3_0cc47abc358d" ><thead>    <tr>        <th class="index_name level0" >name_formatted</th>        <th class="col_heading level0 col0" >IDC208-4000ug-rep1v2b % infectivity:15.1</th>        <th class="col_heading level0 col1" >IDC208-5000ug-rep1v2b % infectivity:10.36</th>        <th class="col_heading level0 col2" >IDC208-6000ug-rep1v2b % infectivity:6.84</th>    </tr>    <tr>        <th class="index_name level0" >name_formatted</th>        <th class="blank" ></th>        <th class="blank" ></th>        <th class="blank" ></th>    </tr></thead><tbody>
                <tr>
                        <th id="T_5ab6b240_5e61_11ea_aee3_0cc47abc358dlevel0_row0" class="row_heading level0 row0" >IDC208-4000ug-rep1v2b % infectivity:15.1</th>
                        <td id="T_5ab6b240_5e61_11ea_aee3_0cc47abc358drow0_col0" class="data row0 col0" >1</td>
                        <td id="T_5ab6b240_5e61_11ea_aee3_0cc47abc358drow0_col1" class="data row0 col1" >0.671</td>
                        <td id="T_5ab6b240_5e61_11ea_aee3_0cc47abc358drow0_col2" class="data row0 col2" >0.665</td>
            </tr>
            <tr>
                        <th id="T_5ab6b240_5e61_11ea_aee3_0cc47abc358dlevel0_row1" class="row_heading level0 row1" >IDC208-5000ug-rep1v2b % infectivity:10.36</th>
                        <td id="T_5ab6b240_5e61_11ea_aee3_0cc47abc358drow1_col0" class="data row1 col0" >0.671</td>
                        <td id="T_5ab6b240_5e61_11ea_aee3_0cc47abc358drow1_col1" class="data row1 col1" >1</td>
                        <td id="T_5ab6b240_5e61_11ea_aee3_0cc47abc358drow1_col2" class="data row1 col2" >0.742</td>
            </tr>
            <tr>
                        <th id="T_5ab6b240_5e61_11ea_aee3_0cc47abc358dlevel0_row2" class="row_heading level0 row2" >IDC208-6000ug-rep1v2b % infectivity:6.84</th>
                        <td id="T_5ab6b240_5e61_11ea_aee3_0cc47abc358drow2_col0" class="data row2 col0" >0.665</td>
                        <td id="T_5ab6b240_5e61_11ea_aee3_0cc47abc358drow2_col1" class="data row2 col1" >0.742</td>
                        <td id="T_5ab6b240_5e61_11ea_aee3_0cc47abc358drow2_col2" class="data row2 col2" >1</td>
            </tr>
    </tbody></table>


    
    
    ******************* IDC403 *******************



![png](analysis_notebook_files/analysis_notebook_54_16.png)



<style  type="text/css" >
    #T_5d2f9690_5e61_11ea_aee3_0cc47abc358drow0_col0 {
            background-color:  #72bf72;
            color:  #000000;
        }    #T_5d2f9690_5e61_11ea_aee3_0cc47abc358drow0_col1 {
            background-color:  #dcfadc;
            color:  #000000;
        }    #T_5d2f9690_5e61_11ea_aee3_0cc47abc358drow0_col2 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_5d2f9690_5e61_11ea_aee3_0cc47abc358drow1_col0 {
            background-color:  #d5f6d5;
            color:  #000000;
        }    #T_5d2f9690_5e61_11ea_aee3_0cc47abc358drow1_col1 {
            background-color:  #73c073;
            color:  #000000;
        }    #T_5d2f9690_5e61_11ea_aee3_0cc47abc358drow1_col2 {
            background-color:  #defbde;
            color:  #000000;
        }    #T_5d2f9690_5e61_11ea_aee3_0cc47abc358drow2_col0 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_5d2f9690_5e61_11ea_aee3_0cc47abc358drow2_col1 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_5d2f9690_5e61_11ea_aee3_0cc47abc358drow2_col2 {
            background-color:  #72bf72;
            color:  #000000;
        }</style><table id="T_5d2f9690_5e61_11ea_aee3_0cc47abc358d" ><thead>    <tr>        <th class="index_name level0" >name_formatted</th>        <th class="col_heading level0 col0" >IDC403-6000ug-rep1v2b % infectivity:0.2</th>        <th class="col_heading level0 col1" >IDC403-7500ug-rep1v2b % infectivity:0.09</th>        <th class="col_heading level0 col2" >IDC403-8500ug-rep1v2b % infectivity:0.04</th>    </tr>    <tr>        <th class="index_name level0" >name_formatted</th>        <th class="blank" ></th>        <th class="blank" ></th>        <th class="blank" ></th>    </tr></thead><tbody>
                <tr>
                        <th id="T_5d2f9690_5e61_11ea_aee3_0cc47abc358dlevel0_row0" class="row_heading level0 row0" >IDC403-6000ug-rep1v2b % infectivity:0.2</th>
                        <td id="T_5d2f9690_5e61_11ea_aee3_0cc47abc358drow0_col0" class="data row0 col0" >1</td>
                        <td id="T_5d2f9690_5e61_11ea_aee3_0cc47abc358drow0_col1" class="data row0 col1" >0.591</td>
                        <td id="T_5d2f9690_5e61_11ea_aee3_0cc47abc358drow0_col2" class="data row0 col2" >0.523</td>
            </tr>
            <tr>
                        <th id="T_5d2f9690_5e61_11ea_aee3_0cc47abc358dlevel0_row1" class="row_heading level0 row1" >IDC403-7500ug-rep1v2b % infectivity:0.09</th>
                        <td id="T_5d2f9690_5e61_11ea_aee3_0cc47abc358drow1_col0" class="data row1 col0" >0.591</td>
                        <td id="T_5d2f9690_5e61_11ea_aee3_0cc47abc358drow1_col1" class="data row1 col1" >1</td>
                        <td id="T_5d2f9690_5e61_11ea_aee3_0cc47abc358drow1_col2" class="data row1 col2" >0.553</td>
            </tr>
            <tr>
                        <th id="T_5d2f9690_5e61_11ea_aee3_0cc47abc358dlevel0_row2" class="row_heading level0 row2" >IDC403-8500ug-rep1v2b % infectivity:0.04</th>
                        <td id="T_5d2f9690_5e61_11ea_aee3_0cc47abc358drow2_col0" class="data row2 col0" >0.523</td>
                        <td id="T_5d2f9690_5e61_11ea_aee3_0cc47abc358drow2_col1" class="data row2 col1" >0.553</td>
                        <td id="T_5d2f9690_5e61_11ea_aee3_0cc47abc358drow2_col2" class="data row2 col2" >1</td>
            </tr>
    </tbody></table>


    
    
    ******************* IDC561 *******************



![png](analysis_notebook_files/analysis_notebook_54_19.png)



<style  type="text/css" >
    #T_60e9611c_5e61_11ea_aee3_0cc47abc358drow0_col0 {
            background-color:  #72bf72;
            color:  #000000;
        }    #T_60e9611c_5e61_11ea_aee3_0cc47abc358drow0_col1 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_60e9611c_5e61_11ea_aee3_0cc47abc358drow0_col2 {
            background-color:  #ddfbdd;
            color:  #000000;
        }    #T_60e9611c_5e61_11ea_aee3_0cc47abc358drow0_col3 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_60e9611c_5e61_11ea_aee3_0cc47abc358drow0_col4 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_60e9611c_5e61_11ea_aee3_0cc47abc358drow1_col0 {
            background-color:  #e1fde1;
            color:  #000000;
        }    #T_60e9611c_5e61_11ea_aee3_0cc47abc358drow1_col1 {
            background-color:  #72bf72;
            color:  #000000;
        }    #T_60e9611c_5e61_11ea_aee3_0cc47abc358drow1_col2 {
            background-color:  #ddfbdd;
            color:  #000000;
        }    #T_60e9611c_5e61_11ea_aee3_0cc47abc358drow1_col3 {
            background-color:  #e2fde2;
            color:  #000000;
        }    #T_60e9611c_5e61_11ea_aee3_0cc47abc358drow1_col4 {
            background-color:  #d8f8d8;
            color:  #000000;
        }    #T_60e9611c_5e61_11ea_aee3_0cc47abc358drow2_col0 {
            background-color:  #d2f4d2;
            color:  #000000;
        }    #T_60e9611c_5e61_11ea_aee3_0cc47abc358drow2_col1 {
            background-color:  #d5f6d5;
            color:  #000000;
        }    #T_60e9611c_5e61_11ea_aee3_0cc47abc358drow2_col2 {
            background-color:  #72bf72;
            color:  #000000;
        }    #T_60e9611c_5e61_11ea_aee3_0cc47abc358drow2_col3 {
            background-color:  #e0fce0;
            color:  #000000;
        }    #T_60e9611c_5e61_11ea_aee3_0cc47abc358drow2_col4 {
            background-color:  #d7f7d7;
            color:  #000000;
        }    #T_60e9611c_5e61_11ea_aee3_0cc47abc358drow3_col0 {
            background-color:  #dffcdf;
            color:  #000000;
        }    #T_60e9611c_5e61_11ea_aee3_0cc47abc358drow3_col1 {
            background-color:  #e1fde1;
            color:  #000000;
        }    #T_60e9611c_5e61_11ea_aee3_0cc47abc358drow3_col2 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_60e9611c_5e61_11ea_aee3_0cc47abc358drow3_col3 {
            background-color:  #72bf72;
            color:  #000000;
        }    #T_60e9611c_5e61_11ea_aee3_0cc47abc358drow3_col4 {
            background-color:  #d6f7d6;
            color:  #000000;
        }    #T_60e9611c_5e61_11ea_aee3_0cc47abc358drow4_col0 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_60e9611c_5e61_11ea_aee3_0cc47abc358drow4_col1 {
            background-color:  #dcfadc;
            color:  #000000;
        }    #T_60e9611c_5e61_11ea_aee3_0cc47abc358drow4_col2 {
            background-color:  #e4fee4;
            color:  #000000;
        }    #T_60e9611c_5e61_11ea_aee3_0cc47abc358drow4_col3 {
            background-color:  #dcfadc;
            color:  #000000;
        }    #T_60e9611c_5e61_11ea_aee3_0cc47abc358drow4_col4 {
            background-color:  #72bf72;
            color:  #000000;
        }</style><table id="T_60e9611c_5e61_11ea_aee3_0cc47abc358d" ><thead>    <tr>        <th class="index_name level0" >name_formatted</th>        <th class="col_heading level0 col0" >IDC561-2200ug-rep1v2b % infectivity:1.21</th>        <th class="col_heading level0 col1" >IDC561-2500ug-rep3e % infectivity:2.32</th>        <th class="col_heading level0 col2" >IDC561-2800ug-rep1v2b % infectivity:0.74</th>        <th class="col_heading level0 col3" >IDC561-3000ug-rep1v2a % infectivity:0.25</th>        <th class="col_heading level0 col4" >IDC561-3500ug-rep1v2a % infectivity:0.18</th>    </tr>    <tr>        <th class="index_name level0" >name_formatted</th>        <th class="blank" ></th>        <th class="blank" ></th>        <th class="blank" ></th>        <th class="blank" ></th>        <th class="blank" ></th>    </tr></thead><tbody>
                <tr>
                        <th id="T_60e9611c_5e61_11ea_aee3_0cc47abc358dlevel0_row0" class="row_heading level0 row0" >IDC561-2200ug-rep1v2b % infectivity:1.21</th>
                        <td id="T_60e9611c_5e61_11ea_aee3_0cc47abc358drow0_col0" class="data row0 col0" >1</td>
                        <td id="T_60e9611c_5e61_11ea_aee3_0cc47abc358drow0_col1" class="data row0 col1" >0.353</td>
                        <td id="T_60e9611c_5e61_11ea_aee3_0cc47abc358drow0_col2" class="data row0 col2" >0.443</td>
                        <td id="T_60e9611c_5e61_11ea_aee3_0cc47abc358drow0_col3" class="data row0 col3" >0.363</td>
                        <td id="T_60e9611c_5e61_11ea_aee3_0cc47abc358drow0_col4" class="data row0 col4" >0.324</td>
            </tr>
            <tr>
                        <th id="T_60e9611c_5e61_11ea_aee3_0cc47abc358dlevel0_row1" class="row_heading level0 row1" >IDC561-2500ug-rep3e % infectivity:2.32</th>
                        <td id="T_60e9611c_5e61_11ea_aee3_0cc47abc358drow1_col0" class="data row1 col0" >0.353</td>
                        <td id="T_60e9611c_5e61_11ea_aee3_0cc47abc358drow1_col1" class="data row1 col1" >1</td>
                        <td id="T_60e9611c_5e61_11ea_aee3_0cc47abc358drow1_col2" class="data row1 col2" >0.444</td>
                        <td id="T_60e9611c_5e61_11ea_aee3_0cc47abc358drow1_col3" class="data row1 col3" >0.383</td>
                        <td id="T_60e9611c_5e61_11ea_aee3_0cc47abc358drow1_col4" class="data row1 col4" >0.408</td>
            </tr>
            <tr>
                        <th id="T_60e9611c_5e61_11ea_aee3_0cc47abc358dlevel0_row2" class="row_heading level0 row2" >IDC561-2800ug-rep1v2b % infectivity:0.74</th>
                        <td id="T_60e9611c_5e61_11ea_aee3_0cc47abc358drow2_col0" class="data row2 col0" >0.443</td>
                        <td id="T_60e9611c_5e61_11ea_aee3_0cc47abc358drow2_col1" class="data row2 col1" >0.444</td>
                        <td id="T_60e9611c_5e61_11ea_aee3_0cc47abc358drow2_col2" class="data row2 col2" >1</td>
                        <td id="T_60e9611c_5e61_11ea_aee3_0cc47abc358drow2_col3" class="data row2 col3" >0.397</td>
                        <td id="T_60e9611c_5e61_11ea_aee3_0cc47abc358drow2_col4" class="data row2 col4" >0.409</td>
            </tr>
            <tr>
                        <th id="T_60e9611c_5e61_11ea_aee3_0cc47abc358dlevel0_row3" class="row_heading level0 row3" >IDC561-3000ug-rep1v2a % infectivity:0.25</th>
                        <td id="T_60e9611c_5e61_11ea_aee3_0cc47abc358drow3_col0" class="data row3 col0" >0.363</td>
                        <td id="T_60e9611c_5e61_11ea_aee3_0cc47abc358drow3_col1" class="data row3 col1" >0.383</td>
                        <td id="T_60e9611c_5e61_11ea_aee3_0cc47abc358drow3_col2" class="data row3 col2" >0.397</td>
                        <td id="T_60e9611c_5e61_11ea_aee3_0cc47abc358drow3_col3" class="data row3 col3" >1</td>
                        <td id="T_60e9611c_5e61_11ea_aee3_0cc47abc358drow3_col4" class="data row3 col4" >0.418</td>
            </tr>
            <tr>
                        <th id="T_60e9611c_5e61_11ea_aee3_0cc47abc358dlevel0_row4" class="row_heading level0 row4" >IDC561-3500ug-rep1v2a % infectivity:0.18</th>
                        <td id="T_60e9611c_5e61_11ea_aee3_0cc47abc358drow4_col0" class="data row4 col0" >0.324</td>
                        <td id="T_60e9611c_5e61_11ea_aee3_0cc47abc358drow4_col1" class="data row4 col1" >0.408</td>
                        <td id="T_60e9611c_5e61_11ea_aee3_0cc47abc358drow4_col2" class="data row4 col2" >0.409</td>
                        <td id="T_60e9611c_5e61_11ea_aee3_0cc47abc358drow4_col3" class="data row4 col3" >0.418</td>
                        <td id="T_60e9611c_5e61_11ea_aee3_0cc47abc358drow4_col4" class="data row4 col4" >1</td>
            </tr>
    </tbody></table>


    
    
    ******************* QA013-2 *******************



![png](analysis_notebook_files/analysis_notebook_54_22.png)



<style  type="text/css" >
    #T_63403530_5e61_11ea_aee3_0cc47abc358drow0_col0 {
            background-color:  #72bf72;
            color:  #000000;
        }    #T_63403530_5e61_11ea_aee3_0cc47abc358drow0_col1 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_63403530_5e61_11ea_aee3_0cc47abc358drow0_col2 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_63403530_5e61_11ea_aee3_0cc47abc358drow1_col0 {
            background-color:  #cef2ce;
            color:  #000000;
        }    #T_63403530_5e61_11ea_aee3_0cc47abc358drow1_col1 {
            background-color:  #72bf72;
            color:  #000000;
        }    #T_63403530_5e61_11ea_aee3_0cc47abc358drow1_col2 {
            background-color:  #c1ebc1;
            color:  #000000;
        }    #T_63403530_5e61_11ea_aee3_0cc47abc358drow2_col0 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_63403530_5e61_11ea_aee3_0cc47abc358drow2_col1 {
            background-color:  #d5f6d5;
            color:  #000000;
        }    #T_63403530_5e61_11ea_aee3_0cc47abc358drow2_col2 {
            background-color:  #72bf72;
            color:  #000000;
        }</style><table id="T_63403530_5e61_11ea_aee3_0cc47abc358d" ><thead>    <tr>        <th class="index_name level0" >name_formatted</th>        <th class="col_heading level0 col0" >QA013-2-2ug-rep2 % infectivity:3.25</th>        <th class="col_heading level0 col1" >QA013-2-4ug-rep2 % infectivity:0.26</th>        <th class="col_heading level0 col2" >QA013-2-4ug-rep3 % infectivity:0.17</th>    </tr>    <tr>        <th class="index_name level0" >name_formatted</th>        <th class="blank" ></th>        <th class="blank" ></th>        <th class="blank" ></th>    </tr></thead><tbody>
                <tr>
                        <th id="T_63403530_5e61_11ea_aee3_0cc47abc358dlevel0_row0" class="row_heading level0 row0" >QA013-2-2ug-rep2 % infectivity:3.25</th>
                        <td id="T_63403530_5e61_11ea_aee3_0cc47abc358drow0_col0" class="data row0 col0" >1</td>
                        <td id="T_63403530_5e61_11ea_aee3_0cc47abc358drow0_col1" class="data row0 col1" >0.877</td>
                        <td id="T_63403530_5e61_11ea_aee3_0cc47abc358drow0_col2" class="data row0 col2" >0.845</td>
            </tr>
            <tr>
                        <th id="T_63403530_5e61_11ea_aee3_0cc47abc358dlevel0_row1" class="row_heading level0 row1" >QA013-2-4ug-rep2 % infectivity:0.26</th>
                        <td id="T_63403530_5e61_11ea_aee3_0cc47abc358drow1_col0" class="data row1 col0" >0.877</td>
                        <td id="T_63403530_5e61_11ea_aee3_0cc47abc358drow1_col1" class="data row1 col1" >1</td>
                        <td id="T_63403530_5e61_11ea_aee3_0cc47abc358drow1_col2" class="data row1 col2" >0.895</td>
            </tr>
            <tr>
                        <th id="T_63403530_5e61_11ea_aee3_0cc47abc358dlevel0_row2" class="row_heading level0 row2" >QA013-2-4ug-rep3 % infectivity:0.17</th>
                        <td id="T_63403530_5e61_11ea_aee3_0cc47abc358drow2_col0" class="data row2 col0" >0.845</td>
                        <td id="T_63403530_5e61_11ea_aee3_0cc47abc358drow2_col1" class="data row2 col1" >0.895</td>
                        <td id="T_63403530_5e61_11ea_aee3_0cc47abc358drow2_col2" class="data row2 col2" >1</td>
            </tr>
    </tbody></table>


    
    
    ******************* QA013-2282dpi *******************



![png](analysis_notebook_files/analysis_notebook_54_25.png)



<style  type="text/css" >
    #T_66f52302_5e61_11ea_aee3_0cc47abc358drow0_col0 {
            background-color:  #72bf72;
            color:  #000000;
        }    #T_66f52302_5e61_11ea_aee3_0cc47abc358drow0_col1 {
            background-color:  #d4f6d4;
            color:  #000000;
        }    #T_66f52302_5e61_11ea_aee3_0cc47abc358drow0_col2 {
            background-color:  #defbde;
            color:  #000000;
        }    #T_66f52302_5e61_11ea_aee3_0cc47abc358drow0_col3 {
            background-color:  #d9f8d9;
            color:  #000000;
        }    #T_66f52302_5e61_11ea_aee3_0cc47abc358drow0_col4 {
            background-color:  #d9f8d9;
            color:  #000000;
        }    #T_66f52302_5e61_11ea_aee3_0cc47abc358drow1_col0 {
            background-color:  #d2f4d2;
            color:  #000000;
        }    #T_66f52302_5e61_11ea_aee3_0cc47abc358drow1_col1 {
            background-color:  #72bf72;
            color:  #000000;
        }    #T_66f52302_5e61_11ea_aee3_0cc47abc358drow1_col2 {
            background-color:  #dcfadc;
            color:  #000000;
        }    #T_66f52302_5e61_11ea_aee3_0cc47abc358drow1_col3 {
            background-color:  #d5f6d5;
            color:  #000000;
        }    #T_66f52302_5e61_11ea_aee3_0cc47abc358drow1_col4 {
            background-color:  #d1f4d1;
            color:  #000000;
        }    #T_66f52302_5e61_11ea_aee3_0cc47abc358drow2_col0 {
            background-color:  #e0fce0;
            color:  #000000;
        }    #T_66f52302_5e61_11ea_aee3_0cc47abc358drow2_col1 {
            background-color:  #e2fde2;
            color:  #000000;
        }    #T_66f52302_5e61_11ea_aee3_0cc47abc358drow2_col2 {
            background-color:  #72bf72;
            color:  #000000;
        }    #T_66f52302_5e61_11ea_aee3_0cc47abc358drow2_col3 {
            background-color:  #acdfac;
            color:  #000000;
        }    #T_66f52302_5e61_11ea_aee3_0cc47abc358drow2_col4 {
            background-color:  #dbf9db;
            color:  #000000;
        }    #T_66f52302_5e61_11ea_aee3_0cc47abc358drow3_col0 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_66f52302_5e61_11ea_aee3_0cc47abc358drow3_col1 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_66f52302_5e61_11ea_aee3_0cc47abc358drow3_col2 {
            background-color:  #b2e3b2;
            color:  #000000;
        }    #T_66f52302_5e61_11ea_aee3_0cc47abc358drow3_col3 {
            background-color:  #72bf72;
            color:  #000000;
        }    #T_66f52302_5e61_11ea_aee3_0cc47abc358drow3_col4 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_66f52302_5e61_11ea_aee3_0cc47abc358drow4_col0 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_66f52302_5e61_11ea_aee3_0cc47abc358drow4_col1 {
            background-color:  #e1fde1;
            color:  #000000;
        }    #T_66f52302_5e61_11ea_aee3_0cc47abc358drow4_col2 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_66f52302_5e61_11ea_aee3_0cc47abc358drow4_col3 {
            background-color:  #e5ffe5;
            color:  #000000;
        }    #T_66f52302_5e61_11ea_aee3_0cc47abc358drow4_col4 {
            background-color:  #72bf72;
            color:  #000000;
        }</style><table id="T_66f52302_5e61_11ea_aee3_0cc47abc358d" ><thead>    <tr>        <th class="index_name level0" >name_formatted</th>        <th class="col_heading level0 col0" >QA013-2282dpi-d10-rep1v2a % infectivity:0.03</th>        <th class="col_heading level0 col1" >QA013-2282dpi-d15-rep1v2a % infectivity:0.14</th>        <th class="col_heading level0 col2" >QA013-2282dpi-d25-rep1v2b % infectivity:27.87</th>        <th class="col_heading level0 col3" >QA013-2282dpi-d40-rep1v2b % infectivity:47.07</th>        <th class="col_heading level0 col4" >QA013-2282dpi-d7-5-rep3e % infectivity:0.04</th>    </tr>    <tr>        <th class="index_name level0" >name_formatted</th>        <th class="blank" ></th>        <th class="blank" ></th>        <th class="blank" ></th>        <th class="blank" ></th>        <th class="blank" ></th>    </tr></thead><tbody>
                <tr>
                        <th id="T_66f52302_5e61_11ea_aee3_0cc47abc358dlevel0_row0" class="row_heading level0 row0" >QA013-2282dpi-d10-rep1v2a % infectivity:0.03</th>
                        <td id="T_66f52302_5e61_11ea_aee3_0cc47abc358drow0_col0" class="data row0 col0" >1</td>
                        <td id="T_66f52302_5e61_11ea_aee3_0cc47abc358drow0_col1" class="data row0 col1" >0.347</td>
                        <td id="T_66f52302_5e61_11ea_aee3_0cc47abc358drow0_col2" class="data row0 col2" >0.245</td>
                        <td id="T_66f52302_5e61_11ea_aee3_0cc47abc358drow0_col3" class="data row0 col3" >0.204</td>
                        <td id="T_66f52302_5e61_11ea_aee3_0cc47abc358drow0_col4" class="data row0 col4" >0.204</td>
            </tr>
            <tr>
                        <th id="T_66f52302_5e61_11ea_aee3_0cc47abc358dlevel0_row1" class="row_heading level0 row1" >QA013-2282dpi-d15-rep1v2a % infectivity:0.14</th>
                        <td id="T_66f52302_5e61_11ea_aee3_0cc47abc358drow1_col0" class="data row1 col0" >0.347</td>
                        <td id="T_66f52302_5e61_11ea_aee3_0cc47abc358drow1_col1" class="data row1 col1" >1</td>
                        <td id="T_66f52302_5e61_11ea_aee3_0cc47abc358drow1_col2" class="data row1 col2" >0.259</td>
                        <td id="T_66f52302_5e61_11ea_aee3_0cc47abc358drow1_col3" class="data row1 col3" >0.233</td>
                        <td id="T_66f52302_5e61_11ea_aee3_0cc47abc358drow1_col4" class="data row1 col4" >0.266</td>
            </tr>
            <tr>
                        <th id="T_66f52302_5e61_11ea_aee3_0cc47abc358dlevel0_row2" class="row_heading level0 row2" >QA013-2282dpi-d25-rep1v2b % infectivity:27.87</th>
                        <td id="T_66f52302_5e61_11ea_aee3_0cc47abc358drow2_col0" class="data row2 col0" >0.245</td>
                        <td id="T_66f52302_5e61_11ea_aee3_0cc47abc358drow2_col1" class="data row2 col1" >0.259</td>
                        <td id="T_66f52302_5e61_11ea_aee3_0cc47abc358drow2_col2" class="data row2 col2" >1</td>
                        <td id="T_66f52302_5e61_11ea_aee3_0cc47abc358drow2_col3" class="data row2 col3" >0.557</td>
                        <td id="T_66f52302_5e61_11ea_aee3_0cc47abc358drow2_col4" class="data row2 col4" >0.191</td>
            </tr>
            <tr>
                        <th id="T_66f52302_5e61_11ea_aee3_0cc47abc358dlevel0_row3" class="row_heading level0 row3" >QA013-2282dpi-d40-rep1v2b % infectivity:47.07</th>
                        <td id="T_66f52302_5e61_11ea_aee3_0cc47abc358drow3_col0" class="data row3 col0" >0.204</td>
                        <td id="T_66f52302_5e61_11ea_aee3_0cc47abc358drow3_col1" class="data row3 col1" >0.233</td>
                        <td id="T_66f52302_5e61_11ea_aee3_0cc47abc358drow3_col2" class="data row3 col2" >0.557</td>
                        <td id="T_66f52302_5e61_11ea_aee3_0cc47abc358drow3_col3" class="data row3 col3" >1</td>
                        <td id="T_66f52302_5e61_11ea_aee3_0cc47abc358drow3_col4" class="data row3 col4" >0.105</td>
            </tr>
            <tr>
                        <th id="T_66f52302_5e61_11ea_aee3_0cc47abc358dlevel0_row4" class="row_heading level0 row4" >QA013-2282dpi-d7-5-rep3e % infectivity:0.04</th>
                        <td id="T_66f52302_5e61_11ea_aee3_0cc47abc358drow4_col0" class="data row4 col0" >0.204</td>
                        <td id="T_66f52302_5e61_11ea_aee3_0cc47abc358drow4_col1" class="data row4 col1" >0.266</td>
                        <td id="T_66f52302_5e61_11ea_aee3_0cc47abc358drow4_col2" class="data row4 col2" >0.191</td>
                        <td id="T_66f52302_5e61_11ea_aee3_0cc47abc358drow4_col3" class="data row4 col3" >0.105</td>
                        <td id="T_66f52302_5e61_11ea_aee3_0cc47abc358drow4_col4" class="data row4 col4" >1</td>
            </tr>
    </tbody></table>



```python
print(f"sel_df has {len(sel_df)} rows. Here are the last few:")

```

    sel_df has 442200 rows. Here are the last few:


### More closely plot the correlation between replicates
Now, I will show the actual correlation plots at the site and mutation level. 



```python
names = batchdf["group"].astype(str)
names = names.tolist()
diffselprefix = "./results/diffsel/summary_"
```


```python
for seltype in ['mutdiffsel', 'positivesitediffsel']:
    print("\n{0} correlations:".format(seltype))
    plots = []
    for g in names:
        plot = diffselprefix + g + '-' + seltype + 'corr.pdf'
        if os.path.isfile(plot):
            plots.append(plot)
        else:
            print("{0} does not exist.".format(plot))
    showPDF(plots, width=1800)
```


```python
#We will use axes with shared ylimits across rows for all plots except for the antibody serum group:



#I will by hand decide what sites to zoom in on. These will include both the 241/289 and C3/465 epitopes of interest. 
ZoomSitesList = ["301", "302", "303", '320', '321','322','322a','323', '324', '325', '326', '327', '328', '329', '330', '331', '332', '333', '334', '415', '416', '417', '441']


#add a column to the sel_df to indicated if the site should be plotted.

sel_df["zoom_site"] = False
sel_df.loc[sel_df['site'].isin(ZoomSitesList), 'zoom_site'] = True
#sel_df

sel_df = sel_df.assign(site_label=lambda x: x['wildtype'] +
                               x['site'].astype('str'))


```


```python
share_ylim_across_rows = {group: ('antibody' not in group)
                          for group in sel_df['group'].unique()}
```


```python
#sel_df["verydiscrip_name"] = (sel_df["serum_alldils_name"].astype(str) + "-vacc, " + sel_df["name_formatted"].astype(str))
#display(HTML(sel_df.head(n=5).to_html(index=False)))
```


```python
os.makedirs("./results/diffsel/zoom/", exist_ok=True) #"./results/diffsel/zoom/"

for group, df in sel_df.groupby('group'):

    plotfile = os.path.join("./results/diffsel/zoom/",
                            f"{group}.pdf")
    print(f"\n\n{'*' * 72}\n {group}, saving to {plotfile}\n")

    fig, axes = dmslogo.facet_plot(
            data=df,#.query('library == @avg_type'),
            x_col='isite',
            show_col='zoom_site',
            gridrow_col='name_formatted',
            share_xlabel=True,
            share_ylabel=True,
            share_ylim_across_rows=share_ylim_across_rows[group],
            wspace=0.6,
            draw_line_kwargs=dict(
                    height_col='positive_diffsel',
                    xtick_col='site',
                    ylabel='immune selection',
                    ),
            draw_logo_kwargs=dict(
                    letter_col='mutation',
                    letter_height_col='mutdiffsel',
                    xtick_col='site_label',
                    xlabel='site',
                    ylabel='immune selection',
                    clip_negative_heights=True,
                    ),
            )
    display(fig)
    fig.savefig(plotfile)
    plt.close(fig)
```

    
    
    ************************************************************************
     118, saving to ./results/diffsel/zoom/118.pdf
    



![png](analysis_notebook_files/analysis_notebook_62_1.png)


    
    
    ************************************************************************
     212, saving to ./results/diffsel/zoom/212.pdf
    



![png](analysis_notebook_files/analysis_notebook_62_3.png)


    
    
    ************************************************************************
     IDC0136, saving to ./results/diffsel/zoom/IDC0136.pdf
    



![png](analysis_notebook_files/analysis_notebook_62_5.png)


    
    
    ************************************************************************
     IDC0337, saving to ./results/diffsel/zoom/IDC0337.pdf
    



![png](analysis_notebook_files/analysis_notebook_62_7.png)


    
    
    ************************************************************************
     IDC208, saving to ./results/diffsel/zoom/IDC208.pdf
    



![png](analysis_notebook_files/analysis_notebook_62_9.png)


    
    
    ************************************************************************
     IDC403, saving to ./results/diffsel/zoom/IDC403.pdf
    



![png](analysis_notebook_files/analysis_notebook_62_11.png)


    
    
    ************************************************************************
     IDC561, saving to ./results/diffsel/zoom/IDC561.pdf
    



![png](analysis_notebook_files/analysis_notebook_62_13.png)


    
    
    ************************************************************************
     QA013-2, saving to ./results/diffsel/zoom/QA013-2.pdf
    



![png](analysis_notebook_files/analysis_notebook_62_15.png)


    
    
    ************************************************************************
     QA013-2282dpi, saving to ./results/diffsel/zoom/QA013-2282dpi.pdf
    



![png](analysis_notebook_files/analysis_notebook_62_17.png)



```python

```


```python

```

### make avgselections df


```python
mediansiteselfiledict={}
medianmutselfiledict = {}
for name in names:
    mediansiteselfiledict[name] = "./results/diffsel/summary_{0}-mediansitediffsel.csv".format(name)
    medianmutselfiledict[name] = "./results/diffsel/summary_{0}-medianmutdiffsel.csv".format(name)

mediansiteselfile_df = pd.DataFrame(list(mediansiteselfiledict.items()), columns=['group','sitediffsel_file'])
medianmutselfile_df = pd.DataFrame(list(medianmutselfiledict.items()), columns=['group','mutdiffsel_file'])
avg_selections = avg_selections.merge(mediansiteselfile_df, left_on='group', right_on='group')
avg_selections = avg_selections.merge(medianmutselfile_df, left_on='group', right_on='group')
avg_selections
```


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-75-34d104e6a7b4> in <module>
          7 mediansiteselfile_df = pd.DataFrame(list(mediansiteselfiledict.items()), columns=['group','sitediffsel_file'])
          8 medianmutselfile_df = pd.DataFrame(list(medianmutselfiledict.items()), columns=['group','mutdiffsel_file'])
    ----> 9 avg_selections = avg_selections.merge(mediansiteselfile_df, left_on='group', right_on='group')
         10 avg_selections = avg_selections.merge(medianmutselfile_df, left_on='group', right_on='group')
         11 avg_selections


    NameError: name 'avg_selections' is not defined



```python
#avg_selections['rabbit_id'] = avg_selections['rabbit_id'].astype(str)
#share_ylim_across_rows = {serum_group: ('antibody' not in serum_group)
#                          for serum_group in avg_selections['rabbit_id'].unique()}
```


```python
selfilecols = ['mutdiffsel_file', 'sitediffsel_file']
avg_sel_df = (dms_tools2.diffsel.df_read_filecols(avg_selections, selfilecols)
          .drop(columns=selfilecols)
          )
```


```python
ZoomSitesList = ["84", "85", "86", "87", "88", "89", "90", "91",  "229", "230", "231", "232", "240","241", "242", "243", "268",  "289", "290", "291","347", "350","351","352", "353", "354", "355", "356", "357", "358", "359", "360", "396", "459", "460", "461", "462", "463", "464", "465", "466", "467", "629"]
hole241 = ["84", "85", "86", "87", "88", "89", "90", "91", "229", "230", "231", "232", "240","241", "242", "243", "268", "289", "290", "291", "347", "629"]
c3465 = ["350","351","352", "353", "354", "355", "356", "357", "358", "359", "360", "396", "459", "460", "461", "462", "463", "464", "465", "466", "467"]

#ZoomSitesList = ["350","351","352", "353", "354", "355", "356", "357", "358", "359", "360",  "396", "459", "460", "461", "462", "463", "464", "465", "466", "467","518", "625", "629"]
avg_sel_df["zoom_site"] = False
avg_sel_df.loc[avg_sel_df['site'].isin(ZoomSitesList), 'zoom_site'] = True

avg_sel_df = avg_sel_df.assign(site_label=lambda x: x['wildtype'] +
                               x['site'].astype('str'))

avg_sel_df["color"] = "gray"
avg_sel_df.loc[avg_sel_df['site'].isin(hole241), 'color'] = "green"
avg_sel_df.loc[avg_sel_df['site'].isin(c3465), 'color'] = "blue"



```


```python

```


```python
os.makedirs("./results/diffsel/AvgAcrossDils_individual/", exist_ok=True)

for rabbit_id, df in avg_sel_df.groupby('serum_alldils_name'):

    plotfile = os.path.join("./results/diffsel/AvgAcrossDils_individual/",
                            f"{rabbit_id}_NewPNGRed.pdf")
    print(f"\n\n{'*' * 72}\nSerum group {rabbit_id}, saving to {plotfile}\n")
    fig, axes = dmslogo.facet_plot(
            data=df,#.query('library == @avg_type'),
            x_col='isite',
            show_col='zoom_site',
            gridrow_col='serum_alldils_name',
            share_xlabel=True,
            share_ylabel=True,
            #share_ylim_across_rows=share_ylim_across_rows[rabbit_id],
            wspace=0.6,
            draw_line_kwargs=dict(
                    height_col='positive_diffsel',
                    xtick_col='site',
                    ylabel='immune selection',
                    ),
            draw_logo_kwargs=dict(
                    letter_col='mutation',
                    letter_height_col='mutdiffsel',
                    xtick_col='site_label',
                    xlabel='site',
                    ylabel='immune selection',
                    clip_negative_heights=True,
                    color_col='color'
                    ),
            )
    display(fig)
    fig.savefig(plotfile)
    plt.close(fig)
```


```python

```


```python

```

## Compute the fraction surviving 
Now we compute the [fraction surviving](https://jbloomlab.github.io/dms_tools2/fracsurvive.html).This caluclation takes into account the level of antibody selection, found in the input file below. We will use mutliple different controls to estimate the error rates to  correct fo, and put the output in its own subdirectory, named according to its control sample. 

This [csv file](/data/BG505_qPCR_master.csv) contains the fraction remaining infectivity for each antibody selected sample, as quatified using pol qPCR and computed based on a standard curve of infecting cells with dilutions of mutant virus (library and experiment specific), with dilutiuons ranging from 0.1 to .0001. 

We first create a batch file to use with [dms2_batch_fracsurvive](https://jbloomlab.github.io/dms_tools2/dms2_batch_fracsurvive.html). 
Note we make the **group** argument the antibody, the **name** argument the replicate, and assign the **sel**, **mock**, and **err** arguments based on the names used for the batch file when generating the counts files above with `dms2_batch_bcsubamp`.
By grouping replicates for the same antibody in the batch file, we tell [dms2_batch_fracsurvive](https://jbloomlab.github.io/dms_tools2/dms2_batch_fracsurvive.html) to analyze these together and take their mean and median.

### Average fraction surviving across multiple antibody dilutions
For select antibodies, we have escape profiles at numerous dilutions. Oftentimes, the antibody concentrations are quite similar (e.g. 3 vs 4 ug/mL), and the single replicate escape profiles look very similar as well. I am averaging across these additional dilutions. 



```python
fracsurvivedir = os.path.join(resultsdir, 'fracsurvive')
if not os.path.isdir(fracsurvivedir):
    os.mkdir(fracsurvivedir)
    
fracsurviveaboveavgdir = os.path.join(resultsdir, 'fracsurviveaboveavg')
if not os.path.isdir(fracsurviveaboveavgdir):
    os.mkdir(fracsurviveaboveavgdir) 
```


```python
fracsurvivebatchavg = pd.read_csv("./data/diffselbatch.csv", header =0)
fracsurvivebatchavg = fracsurvivebatchavg.sort_values(by='group')
display(HTML(fracsurvivebatchavg.to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>group</th>
      <th>name</th>
      <th>sel</th>
      <th>mock</th>
      <th>err</th>
      <th>mds_names</th>
      <th>libfracsurvive</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>118</td>
      <td>4ug-rep2d</td>
      <td>mtvir-rep2d-118-4ug</td>
      <td>mtvir-rep2d</td>
      <td>wt-DNA-rep2</td>
      <td>118-4ug-rep2d</td>
      <td>0.005177</td>
    </tr>
    <tr>
      <td>118</td>
      <td>8ug-rep2d</td>
      <td>mtvir-rep2d-118-8ug</td>
      <td>mtvir-rep2d</td>
      <td>wt-DNA-rep2</td>
      <td>118-8ug-rep2d</td>
      <td>0.000168</td>
    </tr>
    <tr>
      <td>118</td>
      <td>4ug-rep3d</td>
      <td>mtvir-rep3d-118-4ug</td>
      <td>mtvir-rep3d</td>
      <td>wt-DNA-rep3</td>
      <td>118-4ug-rep3d</td>
      <td>0.005182</td>
    </tr>
    <tr>
      <td>118</td>
      <td>8ug-rep3d</td>
      <td>mtvir-rep3d-118-8ug</td>
      <td>mtvir-rep3d</td>
      <td>wt-DNA-rep3</td>
      <td>118-8ug-rep3d</td>
      <td>0.000673</td>
    </tr>
    <tr>
      <td>212</td>
      <td>6ug-rep1v2b</td>
      <td>mtvir-rep1v2b-212-6ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>212-6ug-rep1v2b</td>
      <td>0.034943</td>
    </tr>
    <tr>
      <td>212</td>
      <td>8ug-rep3e</td>
      <td>mtvir-rep3e-212-8ug</td>
      <td>mtvir-rep3e</td>
      <td>wt-DNA-rep3</td>
      <td>212-8ug-rep3e</td>
      <td>0.006215</td>
    </tr>
    <tr>
      <td>212</td>
      <td>10ug-rep1v2a</td>
      <td>mtvir-rep1v2a-212-10ug</td>
      <td>mtvir-rep1v2a</td>
      <td>wt-DNA-rep1</td>
      <td>212-10ug-rep1v2a</td>
      <td>0.002501</td>
    </tr>
    <tr>
      <td>212</td>
      <td>7ug-rep1v2a</td>
      <td>mtvir-rep1v2a-212-7ug</td>
      <td>mtvir-rep1v2a</td>
      <td>wt-DNA-rep1</td>
      <td>212-7ug-rep1v2a</td>
      <td>0.007341</td>
    </tr>
    <tr>
      <td>IDC0136</td>
      <td>8500ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC0136-8500ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC0136-8500ug-rep1v2b</td>
      <td>0.041266</td>
    </tr>
    <tr>
      <td>IDC0136</td>
      <td>7000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC0136-7000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC0136-7000ug-rep1v2b</td>
      <td>0.062973</td>
    </tr>
    <tr>
      <td>IDC0136</td>
      <td>5500ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC0136-5500ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC0136-5500ug-rep1v2b</td>
      <td>0.136001</td>
    </tr>
    <tr>
      <td>IDC0337</td>
      <td>6000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC0337-6000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC0337-6000ug-rep1v2b</td>
      <td>0.059010</td>
    </tr>
    <tr>
      <td>IDC0337</td>
      <td>4500ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC0337-4500ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC0337-4500ug-rep1v2b</td>
      <td>0.091161</td>
    </tr>
    <tr>
      <td>IDC0337</td>
      <td>4000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC0337-4000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC0337-4000ug-rep1v2b</td>
      <td>0.102009</td>
    </tr>
    <tr>
      <td>IDC208</td>
      <td>6000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC208-6000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC208-6000ug-rep1v2b</td>
      <td>0.068419</td>
    </tr>
    <tr>
      <td>IDC208</td>
      <td>5000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC208-5000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC208-5000ug-rep1v2b</td>
      <td>0.103581</td>
    </tr>
    <tr>
      <td>IDC208</td>
      <td>4000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC208-4000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC208-4000ug-rep1v2b</td>
      <td>0.150962</td>
    </tr>
    <tr>
      <td>IDC403</td>
      <td>8500ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC403-8500ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC403-8500ug-rep1v2b</td>
      <td>0.000359</td>
    </tr>
    <tr>
      <td>IDC403</td>
      <td>7500ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC403-7500ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC403-7500ug-rep1v2b</td>
      <td>0.000893</td>
    </tr>
    <tr>
      <td>IDC403</td>
      <td>6000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC403-6000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC403-6000ug-rep1v2b</td>
      <td>0.001969</td>
    </tr>
    <tr>
      <td>IDC561</td>
      <td>3000ug-rep1v2a</td>
      <td>mtvir-rep1v2a-561-IgG-3000ug</td>
      <td>mtvir-rep1v2a</td>
      <td>wt-DNA-rep1</td>
      <td>IDC561-3000ug-rep1v2a</td>
      <td>0.002478</td>
    </tr>
    <tr>
      <td>IDC561</td>
      <td>2500ug-rep3e</td>
      <td>mtvir-rep3e-561-IgG-2500ug</td>
      <td>mtvir-rep3e</td>
      <td>wt-DNA-rep3</td>
      <td>IDC561-2500ug-rep3e</td>
      <td>0.023230</td>
    </tr>
    <tr>
      <td>IDC561</td>
      <td>3500ug-rep1v2a</td>
      <td>mtvir-rep1v2a-561-IgG-3500ug</td>
      <td>mtvir-rep1v2a</td>
      <td>wt-DNA-rep1</td>
      <td>IDC561-3500ug-rep1v2a</td>
      <td>0.001820</td>
    </tr>
    <tr>
      <td>IDC561</td>
      <td>2800ug-rep1v2b</td>
      <td>mtvir-rep1v2b-561-IgG-2800ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC561-2800ug-rep1v2b</td>
      <td>0.007367</td>
    </tr>
    <tr>
      <td>IDC561</td>
      <td>2200ug-rep1v2b</td>
      <td>mtvir-rep1v2b-561-IgG-2200ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC561-2200ug-rep1v2b</td>
      <td>0.012112</td>
    </tr>
    <tr>
      <td>QA013-2</td>
      <td>2ug-rep2</td>
      <td>mtvir-rep2-QA013-2-2ug</td>
      <td>mtvir-rep2</td>
      <td>wt-DNA-rep2</td>
      <td>QA013-2-2ug-rep2</td>
      <td>0.032488</td>
    </tr>
    <tr>
      <td>QA013-2</td>
      <td>4ug-rep3</td>
      <td>mtvir-rep3-QA013-2-4ug</td>
      <td>mtvir-rep3</td>
      <td>wt-DNA-rep3</td>
      <td>QA013-2-4ug-rep3</td>
      <td>0.001661</td>
    </tr>
    <tr>
      <td>QA013-2</td>
      <td>4ug-rep2</td>
      <td>mtvir-rep2-QA013-2-4ug</td>
      <td>mtvir-rep2</td>
      <td>wt-DNA-rep2</td>
      <td>QA013-2-4ug-rep2</td>
      <td>0.002637</td>
    </tr>
    <tr>
      <td>QA013-2282dpi</td>
      <td>d40-rep1v2b</td>
      <td>mtvir-rep1v2b-QA013-2282dpi-d40</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>QA013-2282dpi-d40-rep1v2b</td>
      <td>0.470669</td>
    </tr>
    <tr>
      <td>QA013-2282dpi</td>
      <td>d25-rep1v2b</td>
      <td>mtvir-rep1v2b-QA013-2282dpi-d25</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>QA013-2282dpi-d25-rep1v2b</td>
      <td>0.278701</td>
    </tr>
    <tr>
      <td>QA013-2282dpi</td>
      <td>d15-rep1v2a</td>
      <td>mtvir-rep1v2a-QA013-2282dpi-d15</td>
      <td>mtvir-rep1v2a</td>
      <td>wt-DNA-rep1</td>
      <td>QA013-2282dpi-d15-rep1v2a</td>
      <td>0.001444</td>
    </tr>
    <tr>
      <td>QA013-2282dpi</td>
      <td>d10-rep1v2a</td>
      <td>mtvir-rep1v2a-QA013-2282dpi-d10</td>
      <td>mtvir-rep1v2a</td>
      <td>wt-DNA-rep1</td>
      <td>QA013-2282dpi-d10-rep1v2a</td>
      <td>0.000317</td>
    </tr>
    <tr>
      <td>QA013-2282dpi</td>
      <td>d7-5-rep3e</td>
      <td>mtvir-rep3e-QA013-2282dpi-d7-5</td>
      <td>mtvir-rep3e</td>
      <td>wt-DNA-rep3</td>
      <td>QA013-2282dpi-d7-5-rep3e</td>
      <td>0.000388</td>
    </tr>
  </tbody>
</table>



```python
fracsurvivebatchavgcopy = fracsurvivebatchavg.copy()


names = fracsurvivebatchavg["group"].astype(str) + "-" + fracsurvivebatchavg["name"].astype(str)
names = names.tolist()
```

We now `run dms2_batch_survive` twice with the following difference:
    1. First we run it simply computing the fraction surviving for each mutation.
    2. Then we run it with the --aboveavg yes option to compute the fraction surviving for each mutation above the overall library average.

Note how the results for these two different runs are output to two different subdirectories.



```python
fracsurvivebatch = fracsurvivebatchavg.copy()
fracsurvivebatchfile = os.path.join(fracsurvivedir, 'batch.csv')
print("Here is the batch input that we write to the CSV file {0}:".format(fracsurvivebatchfile))
display(HTML(fracsurvivebatch.to_html(index=False)))
fracsurvivebatch.to_csv(fracsurvivebatchfile, index=False, encoding='utf-8')


for (arg_aboveavg, outdir) in [('', fracsurvivedir), ('--aboveavg yes', fracsurviveaboveavgdir)]:
    print("\nRunning dms2_batch_fracsurvive {0}and writing output to {1}".format(
            {'':'', '--aboveavg yes':'with `--aboveavg yes` '}[arg_aboveavg], outdir))
    log = !dms2_batch_fracsurvive \
            --summaryprefix summary \
            --batchfile {fracsurvivebatchfile} \
            --outdir {outdir} \
            --indir {renumberedcountsdir} \
            --use_existing {use_existing} \
            {arg_aboveavg} 
    print("Completed run.")
```

    Here is the batch input that we write to the CSV file ./results/fracsurvive/batch.csv:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>group</th>
      <th>name</th>
      <th>sel</th>
      <th>mock</th>
      <th>err</th>
      <th>mds_names</th>
      <th>libfracsurvive</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>118</td>
      <td>4ug-rep2d</td>
      <td>mtvir-rep2d-118-4ug</td>
      <td>mtvir-rep2d</td>
      <td>wt-DNA-rep2</td>
      <td>118-4ug-rep2d</td>
      <td>0.005177</td>
    </tr>
    <tr>
      <td>118</td>
      <td>8ug-rep2d</td>
      <td>mtvir-rep2d-118-8ug</td>
      <td>mtvir-rep2d</td>
      <td>wt-DNA-rep2</td>
      <td>118-8ug-rep2d</td>
      <td>0.000168</td>
    </tr>
    <tr>
      <td>118</td>
      <td>4ug-rep3d</td>
      <td>mtvir-rep3d-118-4ug</td>
      <td>mtvir-rep3d</td>
      <td>wt-DNA-rep3</td>
      <td>118-4ug-rep3d</td>
      <td>0.005182</td>
    </tr>
    <tr>
      <td>118</td>
      <td>8ug-rep3d</td>
      <td>mtvir-rep3d-118-8ug</td>
      <td>mtvir-rep3d</td>
      <td>wt-DNA-rep3</td>
      <td>118-8ug-rep3d</td>
      <td>0.000673</td>
    </tr>
    <tr>
      <td>212</td>
      <td>6ug-rep1v2b</td>
      <td>mtvir-rep1v2b-212-6ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>212-6ug-rep1v2b</td>
      <td>0.034943</td>
    </tr>
    <tr>
      <td>212</td>
      <td>8ug-rep3e</td>
      <td>mtvir-rep3e-212-8ug</td>
      <td>mtvir-rep3e</td>
      <td>wt-DNA-rep3</td>
      <td>212-8ug-rep3e</td>
      <td>0.006215</td>
    </tr>
    <tr>
      <td>212</td>
      <td>10ug-rep1v2a</td>
      <td>mtvir-rep1v2a-212-10ug</td>
      <td>mtvir-rep1v2a</td>
      <td>wt-DNA-rep1</td>
      <td>212-10ug-rep1v2a</td>
      <td>0.002501</td>
    </tr>
    <tr>
      <td>212</td>
      <td>7ug-rep1v2a</td>
      <td>mtvir-rep1v2a-212-7ug</td>
      <td>mtvir-rep1v2a</td>
      <td>wt-DNA-rep1</td>
      <td>212-7ug-rep1v2a</td>
      <td>0.007341</td>
    </tr>
    <tr>
      <td>IDC0136</td>
      <td>8500ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC0136-8500ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC0136-8500ug-rep1v2b</td>
      <td>0.041266</td>
    </tr>
    <tr>
      <td>IDC0136</td>
      <td>7000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC0136-7000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC0136-7000ug-rep1v2b</td>
      <td>0.062973</td>
    </tr>
    <tr>
      <td>IDC0136</td>
      <td>5500ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC0136-5500ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC0136-5500ug-rep1v2b</td>
      <td>0.136001</td>
    </tr>
    <tr>
      <td>IDC0337</td>
      <td>6000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC0337-6000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC0337-6000ug-rep1v2b</td>
      <td>0.059010</td>
    </tr>
    <tr>
      <td>IDC0337</td>
      <td>4500ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC0337-4500ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC0337-4500ug-rep1v2b</td>
      <td>0.091161</td>
    </tr>
    <tr>
      <td>IDC0337</td>
      <td>4000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC0337-4000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC0337-4000ug-rep1v2b</td>
      <td>0.102009</td>
    </tr>
    <tr>
      <td>IDC208</td>
      <td>6000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC208-6000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC208-6000ug-rep1v2b</td>
      <td>0.068419</td>
    </tr>
    <tr>
      <td>IDC208</td>
      <td>5000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC208-5000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC208-5000ug-rep1v2b</td>
      <td>0.103581</td>
    </tr>
    <tr>
      <td>IDC208</td>
      <td>4000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC208-4000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC208-4000ug-rep1v2b</td>
      <td>0.150962</td>
    </tr>
    <tr>
      <td>IDC403</td>
      <td>8500ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC403-8500ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC403-8500ug-rep1v2b</td>
      <td>0.000359</td>
    </tr>
    <tr>
      <td>IDC403</td>
      <td>7500ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC403-7500ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC403-7500ug-rep1v2b</td>
      <td>0.000893</td>
    </tr>
    <tr>
      <td>IDC403</td>
      <td>6000ug-rep1v2b</td>
      <td>mtvir-rep1v2b-IDC403-6000ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC403-6000ug-rep1v2b</td>
      <td>0.001969</td>
    </tr>
    <tr>
      <td>IDC561</td>
      <td>3000ug-rep1v2a</td>
      <td>mtvir-rep1v2a-561-IgG-3000ug</td>
      <td>mtvir-rep1v2a</td>
      <td>wt-DNA-rep1</td>
      <td>IDC561-3000ug-rep1v2a</td>
      <td>0.002478</td>
    </tr>
    <tr>
      <td>IDC561</td>
      <td>2500ug-rep3e</td>
      <td>mtvir-rep3e-561-IgG-2500ug</td>
      <td>mtvir-rep3e</td>
      <td>wt-DNA-rep3</td>
      <td>IDC561-2500ug-rep3e</td>
      <td>0.023230</td>
    </tr>
    <tr>
      <td>IDC561</td>
      <td>3500ug-rep1v2a</td>
      <td>mtvir-rep1v2a-561-IgG-3500ug</td>
      <td>mtvir-rep1v2a</td>
      <td>wt-DNA-rep1</td>
      <td>IDC561-3500ug-rep1v2a</td>
      <td>0.001820</td>
    </tr>
    <tr>
      <td>IDC561</td>
      <td>2800ug-rep1v2b</td>
      <td>mtvir-rep1v2b-561-IgG-2800ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC561-2800ug-rep1v2b</td>
      <td>0.007367</td>
    </tr>
    <tr>
      <td>IDC561</td>
      <td>2200ug-rep1v2b</td>
      <td>mtvir-rep1v2b-561-IgG-2200ug</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>IDC561-2200ug-rep1v2b</td>
      <td>0.012112</td>
    </tr>
    <tr>
      <td>QA013-2</td>
      <td>2ug-rep2</td>
      <td>mtvir-rep2-QA013-2-2ug</td>
      <td>mtvir-rep2</td>
      <td>wt-DNA-rep2</td>
      <td>QA013-2-2ug-rep2</td>
      <td>0.032488</td>
    </tr>
    <tr>
      <td>QA013-2</td>
      <td>4ug-rep3</td>
      <td>mtvir-rep3-QA013-2-4ug</td>
      <td>mtvir-rep3</td>
      <td>wt-DNA-rep3</td>
      <td>QA013-2-4ug-rep3</td>
      <td>0.001661</td>
    </tr>
    <tr>
      <td>QA013-2</td>
      <td>4ug-rep2</td>
      <td>mtvir-rep2-QA013-2-4ug</td>
      <td>mtvir-rep2</td>
      <td>wt-DNA-rep2</td>
      <td>QA013-2-4ug-rep2</td>
      <td>0.002637</td>
    </tr>
    <tr>
      <td>QA013-2282dpi</td>
      <td>d40-rep1v2b</td>
      <td>mtvir-rep1v2b-QA013-2282dpi-d40</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>QA013-2282dpi-d40-rep1v2b</td>
      <td>0.470669</td>
    </tr>
    <tr>
      <td>QA013-2282dpi</td>
      <td>d25-rep1v2b</td>
      <td>mtvir-rep1v2b-QA013-2282dpi-d25</td>
      <td>mtvir-rep1v2b</td>
      <td>wt-DNA-rep1</td>
      <td>QA013-2282dpi-d25-rep1v2b</td>
      <td>0.278701</td>
    </tr>
    <tr>
      <td>QA013-2282dpi</td>
      <td>d15-rep1v2a</td>
      <td>mtvir-rep1v2a-QA013-2282dpi-d15</td>
      <td>mtvir-rep1v2a</td>
      <td>wt-DNA-rep1</td>
      <td>QA013-2282dpi-d15-rep1v2a</td>
      <td>0.001444</td>
    </tr>
    <tr>
      <td>QA013-2282dpi</td>
      <td>d10-rep1v2a</td>
      <td>mtvir-rep1v2a-QA013-2282dpi-d10</td>
      <td>mtvir-rep1v2a</td>
      <td>wt-DNA-rep1</td>
      <td>QA013-2282dpi-d10-rep1v2a</td>
      <td>0.000317</td>
    </tr>
    <tr>
      <td>QA013-2282dpi</td>
      <td>d7-5-rep3e</td>
      <td>mtvir-rep3e-QA013-2282dpi-d7-5</td>
      <td>mtvir-rep3e</td>
      <td>wt-DNA-rep3</td>
      <td>QA013-2282dpi-d7-5-rep3e</td>
      <td>0.000388</td>
    </tr>
  </tbody>
</table>


    
    Running dms2_batch_fracsurvive and writing output to ./results/fracsurvive
    Completed run.
    
    Running dms2_batch_fracsurvive with `--aboveavg yes` and writing output to ./results/fracsurviveaboveavg
    Completed run.


Running [dms2_batch_fracsurvive](https://jbloomlab.github.io/dms_tools2/dms2_batch_fracsurvive.html) creates plots showing how the fraction surviving estimates correlate among replicates. 
These plots have names like `summary_*-avgfracsurvivecorr.pdf` and `summary_*-maxfracsurvivecorr.pdf`. 
Note that the plots show the correlations between all pairs, and also on the diagonal show the density of the different selection values for each replicates (most of them are close to zero).


```python
fracsurviveprefix = os.path.join(fracsurvivedir, 'summary_')
groups = fracsurvivebatchavg['group'].unique()

```


```python
for seltype in ['avgfracsurvive', 'maxfracsurvive']:
    print("\n{0} correlations:".format(seltype))
    plots = []
    for g in groups:
        plot = fracsurviveprefix + g + '-' + seltype + 'corr.pdf'
        if os.path.isfile(plot):
            plots.append(plot)
        else:
            print("{0} does not exist.".format(plot))
    showPDF(plots, width=1800)
```

    
    avgfracsurvive correlations:



![png](analysis_notebook_files/analysis_notebook_82_1.png)


    
    maxfracsurvive correlations:



![png](analysis_notebook_files/analysis_notebook_82_3.png)


Next, we can look at the correlation for the `fraction surviving above average` values.


```python
fracsurviveaboveavgprefix = os.path.join(fracsurviveaboveavgdir, 'summary_')

for seltype in ['avgfracsurvive', 'maxfracsurvive']:
    print("\n{0} correlations:".format(seltype))
    plots = []
    for g in groups:
        plot = fracsurviveaboveavgprefix + g + '-' + seltype + 'corr.pdf'
        if os.path.isfile(plot):
            plots.append(plot)
        else:
            print("{0} does not exist.".format(plot))
    showPDF(plots, width=1800)
```

    
    avgfracsurvive correlations:



![png](analysis_notebook_files/analysis_notebook_84_1.png)


    
    maxfracsurvive correlations:



![png](analysis_notebook_files/analysis_notebook_84_3.png)


Now, lets look at the median `average fraction suurviving above average` for each antibody.


```python
phiaboveavg_sub_tuple = [
 ('1-18', './results/fracsurviveaboveavg/summary_118-mediansitefracsurvive.csv'),#]
 ('VRC01', './results/fracsurviveaboveavg/summary_VRC01-mediansitefracsurvive.csv'),
 ('3BNC117', './results/fracsurviveaboveavg/summary_3BNC117-mediansitefracsurvive.csv')]

names = []
diffselfiles = []
phiaboveavg_dict = {}
for tup in phiaboveavg_sub_tuple:
    names.append(tup[0])
    diffselfiles.append(tup[1])
    phiaboveavg_dict[tup[0]] = tup[1]
```


```python
import warnings
warnings.filterwarnings('ignore')
```


```python
diffseltype = "avgfracsurvive"
plotfile = "./results/fracsurviveaboveavg/CD4bsAb_subset_median_avgsitefracsurvive.pdf"
dms_tools2.plot.plotSiteDiffSel(names, diffselfiles, plotfile, diffseltype, maxcol=1, white_bg=True)
showPDF(plotfile)
```


![png](analysis_notebook_files/analysis_notebook_88_0.png)



```python
groups = fracsurvivebatch['group'].unique()
print(groups)
```

    ['118' '3BNC117' 'VRC01']



```python
for antibody in groups:
    mutdiffsel = os.path.join(fracsurviveaboveavgdir, 'summary_{0}-medianmutfracsurvive.csv'.format(antibody))
    logoplot = os.path.join(fracsurviveaboveavgdir, '{0}-median_fracsurvive.pdf'.format(antibody))
        #scale bar unit is maximum effect
    mutdiffseldf = pd.read_csv(mutdiffsel)
    scaleunit = '{0:.1g}'.format(mutdiffseldf['mutfracsurvive'].max())
    scalelabel = '"fraction surviving above avg = {0}"'.format(scaleunit)
    logoname = '{0}-median'.format(antibody)
    print("Creating logo plot for {0} from {1}".format(antibody, mutdiffsel))
    log = !dms2_logoplot \
            --fracsurvive {mutdiffsel} \
            --name {logoname} \
            --outdir {fracsurviveaboveavgdir} \
            --sepline no \
            --nperline 84 \
            --overlay1 {mutdiffsel} wildtype wildtype \
            --scalebar {scaleunit} {scalelabel} \
            --underlay yes \
            --use_existing {use_existing}
    #showPDF(logoplot)
    
    mutdiffsel = os.path.join(fracsurvivedir, 'summary_{0}-medianmutfracsurvive.csv'.format(antibody))
    logoplot = os.path.join(fracsurvivedir, '{0}-median_fracsurvive.pdf'.format(antibody))
        #scale bar unit is maximum effect
    mutdiffseldf = pd.read_csv(mutdiffsel)
    scaleunit = '{0:.1g}'.format(mutdiffseldf['mutfracsurvive'].max())
    scalelabel = '"fraction surviving = {0}"'.format(scaleunit)
    logoname = '{0}-median'.format(antibody)
    print("Creating logo plot for {0} from {1}".format(antibody, mutdiffsel))
    log = !dms2_logoplot \
            --fracsurvive {mutdiffsel} \
            --name {logoname} \
            --outdir {fracsurvivedir} \
            --sepline no \
            --nperline 84 \
            --overlay1 {mutdiffsel} wildtype wildtype \
            --scalebar {scaleunit} {scalelabel} \
            --underlay yes \
            --use_existing {use_existing}
    #showPDF(logoplot)
```

    Creating logo plot for 118 from ./results/fracsurviveaboveavg/summary_118-medianmutfracsurvive.csv
    Creating logo plot for 118 from ./results/fracsurvive/summary_118-medianmutfracsurvive.csv
    Creating logo plot for 3BNC117 from ./results/fracsurviveaboveavg/summary_3BNC117-medianmutfracsurvive.csv
    Creating logo plot for 3BNC117 from ./results/fracsurvive/summary_3BNC117-medianmutfracsurvive.csv
    Creating logo plot for VRC01 from ./results/fracsurviveaboveavg/summary_VRC01-medianmutfracsurvive.csv
    Creating logo plot for VRC01 from ./results/fracsurvive/summary_VRC01-medianmutfracsurvive.csv


The overall effect size of mutations at the site under strongest selection is larger for 1-18 than 3BNC117 or VRC01 when examining differential selection. However, there is still lesse selection at the cannonically defined CD4bs. 

I will focus onescess fraction surviving in the paper, as there is not an a priori reason to use one or the other Also, for 3BNC117 and VRC01, there were replicates that were not neutralized as stringently as the 1-18 replicates.  

## Identifying significant escape sites.
As explained by Jesse [here](https://github.com/jbloomlab/computational_notebooks/blob/master/jbloom/2018/MAP_Identifying_Significant_Escape/analysis_notebook.ipynb), we will use [dms_tools2.plot.findSigSel](https://jbloomlab.github.io/dms_tools2/dms_tools2.plot.html#dms_tools2.plot.findSigSel) to fit a [gamma distribution](https://en.wikipedia.org/wiki/Gamma_distribution) to binned `avgfracsurvive` values for each antibody using robust regression. We then identify sites that are larger than expected from this distribution at an FDR=0.01. These are significantly higher than expected from a gamma-fit distribution. 


```python
phiaboveavg_sub_tuple = [
 ('118', './results/fracsurviveaboveavg/summary_118-mediansitefracsurvive.csv'),
 ('VRC01', './results/fracsurviveaboveavg/summary_VRC01-mediansitefracsurvive.csv'),
 ('3BNC117', './results/fracsurviveaboveavg/summary_3BNC117-mediansitefracsurvive.csv')]
```


```python
antibodies = []
datafiles = []
for tup in phiaboveavg_sub_tuple:
    antibodies.append(tup[0])
    datafiles.append(tup[1])

fitresultsdir = './results/fracsurvive/SignificantEscape/'
if not os.path.isdir(fitresultsdir):
    os.mkdir(fitresultsdir)
    

def sigSitesDisplay(valcol, fdrs, antibodies, datafiles):
    """Display the significant sites at some FDRs.
    
    Args:
        `valcol` (str): column to use, `maxfracsurvive` or `avgfracsurvive`
        `fdr` (float): FDR
        `antibodies` (list): antibody names
        `datafiles` (list): sitefracsurvive files
    """
    
    for (a, datafile) in zip(antibodies, datafiles):
        print("\n*************************\nAnalyzing antibody {0}".format(a))
        plots = []
        df_printed = True
        for fdr in sorted(fdrs, reverse=True):
            plot = os.path.join(fitresultsdir, '{0}_{1}_{2}_fit.pdf'.format(a, valcol, fdr))
            df = dms_tools2.plot.findSigSel(
                    pd.read_csv(datafile),
                    valcol,
                    plot,
                    fdr=fdr,
                    title='{0} at FDR {1}'.format(a, fdr),
                    )[0]
            plots.append(plot)
            if not df_printed:
                print(df.query('sig').sort_values('Q')[['site', valcol, 'Q']]
                        .reset_index(drop=True))
                df_printed = True # only print for lowest FDR

        showPDF(plots, width=400)

sigsites = {}
def sigSitesDisplaySave(valcol, fdrs, antibodies, datafiles, sigsitedict):
    """Display the significant sites at a single FDR, and saves that list of sites to s dictionary keyed by antibody
    
    Args:
        `valcol` (str): column to use, `maxfracsurvive` or `avgfracsurvive`
        `fdr` (float): FDR
        `antibodies` (list): antibody names
        `datafiles` (list): sitefracsurvive files
    """
    
    for (a, datafile) in zip(antibodies, datafiles):
        print("\n*************************\nAnalyzing antibody {0}".format(a))
        plots = []
        df_printed = False
        for fdr in sorted(fdrs, reverse=True):
            if  len(fdrs)>1:
                raise NameError('only 1 FDR allowed for this function')
            else:
                plot = os.path.join(fitresultsdir, '{0}_{1}_{2}_fit.pdf'.format(a, valcol, fdr))
                df = dms_tools2.plot.findSigSel(
                        pd.read_csv(datafile),
                        valcol,
                        plot,
                        fdr=fdr,
                        title='{0} at FDR {1}'.format(a, fdr),
                        )[0]
                plots.append(plot)
                if not df_printed:
                    print(df.query('sig').sort_values('Q')[['site', valcol, 'Q']]
                            .reset_index(drop=True))
                    df_printed = True # only print for lowest FDR
                    siglist = df.query('sig')['site'].tolist()
                    
                    print(siglist)
                    sigsitedict[a] = siglist
        showPDF(plots, width=400)
```


```python
notaboveavgdatafiles =[]
for file in datafiles:
    mutfile = file.replace("aboveavg","")
    notaboveavgdatafiles.append(mutfile)
```

First we look fo sites of significant escape using the `avgfracsurvive` for a site.


We also print the Q-values for all sites found at even the most lenient FDR: the Q values are the minimum FDR at which each site would be called:



```python
sigsites = {}
fdrs = [0.01]
sigSitesDisplaySave('avgfracsurvive', fdrs, antibodies, notaboveavgdatafiles, sigsites)
```

    
    *************************
    Analyzing antibody 118
      site  avgfracsurvive             Q
    0  304        0.010021  4.097188e-16
    1  119        0.006338  9.844564e-06
    2  207        0.006232  1.234732e-05
    3  318        0.005310  1.667233e-03
    ['304', '119', '207', '318']



![png](analysis_notebook_files/analysis_notebook_97_1.png)


    
    *************************
    Analyzing antibody VRC01
      site  avgfracsurvive         Q
    0  279        0.013127  0.000001
    1  326        0.012705  0.000004
    2  369        0.012131  0.000036
    3  197        0.011616  0.000227
    4  209        0.011280  0.000693
    ['279', '326', '369', '197', '209']



![png](analysis_notebook_files/analysis_notebook_97_3.png)


    
    *************************
    Analyzing antibody 3BNC117
       site  avgfracsurvive             Q
    0   207        0.020230  1.441091e-14
    1   304        0.019376  7.440499e-13
    2   197        0.017588  4.181925e-09
    3   318        0.017310  1.170851e-08
    4   471        0.016126  1.884601e-06
    5   279        0.015686  9.857346e-06
    6   119        0.015012  1.216440e-04
    7   308        0.014894  1.597877e-04
    8   209        0.014873  1.597877e-04
    9   206        0.014474  6.223021e-04
    10  182        0.014292  1.075935e-03
    11  120        0.014017  2.489072e-03
    12  204        0.013997  2.489072e-03
    13  369        0.013947  2.734806e-03
    14  274        0.013707  5.587751e-03
    ['207', '304', '197', '318', '471', '279', '119', '308', '209', '206', '182', '120', '204', '369', '274']



![png](analysis_notebook_files/analysis_notebook_97_5.png)


Now, lets look at the significant sites of escape in a single table.


```python
print(sigsites)
```

    {'118': ['304', '119', '207', '318'], 'VRC01': ['279', '326', '369', '197', '209'], '3BNC117': ['207', '304', '197', '318', '471', '279', '119', '308', '209', '206', '182', '120', '204', '369', '274']}



```python
from natsort import natsorted
sigsitestabledict = {}
for ab in sigsites.keys():
    sitelist = sigsites[ab]
    sortedsitelist = natsorted(sitelist)
    sitestring = ', '.join(sortedsitelist)
    sigsitestabledict[ab] = sitestring


sigsitestable = pd.DataFrame.from_dict(sigsitestabledict, orient='index')
with pd.option_context('display.max_colwidth', -1): 
    display(HTML(sigsitestable.to_html(index=True)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>0</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>118</th>
      <td>119, 207, 304, 318</td>
    </tr>
    <tr>
      <th>VRC01</th>
      <td>197, 209, 279, 326, 369</td>
    </tr>
    <tr>
      <th>3BNC117</th>
      <td>119, 120, 182, 197, 204, 206, 207, 209, 274, 279, 304, 308, 318, 369, 471</td>
    </tr>
  </tbody>
</table>



```python
fracsurviveaboveavg_dict = {}
fracsurviveaboveavg_dict["118"] = "./results/fracsurviveaboveavg/summary_118-medianmutfracsurvive.csv"

fracsurviveaboveavg_dict["VRC01"] = "./results/fracsurviveaboveavg/summary_VRC01-medianmutfracsurvive.csv"
fracsurviveaboveavg_dict["3BNC117"] = "./results/fracsurviveaboveavg/summary_3BNC117-medianmutfracsurvive.csv" 
fracsurviveaboveavg_dict["101074"] = "./results/fracsurviveaboveavg/summary_101074-medianmutfracsurvive.csv"

sitefracsurviveaboveavg_dict = {}
sitefracsurviveaboveavg_dict["118"] = "./results/fracsurviveaboveavg/summary_118-mediansitefracsurvive.csv"
sitefracsurviveaboveavg_dict["VRC01"] = "./results/fracsurviveaboveavg/summary_VRC01-mediansitefracsurvive.csv" 
sitefracsurviveaboveavg_dict["3BNC117"] = "./results/fracsurviveaboveavg/summary_3BNC117-mediansitefracsurvive.csv"
sitefracsurviveaboveavg_dict["101074"] = "./results/fracsurviveaboveavg/summary_101074-mediansitefracsurvive.csv"

```


```python
import seaborn as sns
from scipy import stats
%matplotlib inline
```



## Generate figure ready plots 
Now, I will replot the data in a number of ways for paper figures. 


```python
summplotdir = "./results/fracsurviveaboveavg/SummaryPlots/"
if not os.path.isdir(summplotdir):
    os.mkdir(summplotdir)
```


```python
import re
import os
import math
import natsort
#import pandas
import numpy
import scipy.stats
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib import rc
from plotnine import *
# set ggplot theme
theme_set(theme_bw(base_size=12)) 

import seaborn


from dms_tools2 import CODONS, AAS, AAS_WITHSTOP
import dms_tools2.utils
COLOR_BLIND_PALETTE = ["#000000", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]

def plotSiteDiffSel_YLimDiff(names, diffselfiles, plotfile, 
        diffseltype, y_lim=None, highlighted_sites =[], underlay_1 = [], point = [], maxcol=2, white_bg=False):
    """Plot site diffsel or fracsurvive along sequence.
    Despite the function name, this function can be used to
    plot either differential selection or fraction surviving.
    Args:
        `names` (list or series)
            Names of samples for which we plot statistics.
        `diffselfiles` (list or series)
            ``*sitediffsel.csv`` files from ``dms2_diffsel`` or
            ``*sitefracsurvive.csv`` files from ``dms2_fracsurvive``.
        `plotfile` (str)
            Name of created PDF plot file.
        `diffseltype` (str)
            Type of diffsel or fracsurvive to plot:
                - `positive`: positive sitediffsel
                - `total`: positive and negative sitediffsel
                - `max`: maximum mutdiffsel
                - `minmax`: minimum and maximum mutdiffsel
                - `avgfracsurvive`: total site fracsurvive
                - `maxfracsurvive`: max mutfracsurvive at site
        `maxcol` (int)
            Number of columns in faceted plot.
        `white_bg` (bool)
            Plots will have a white background with limited other formatting.
            
        DOCUMENT CHANGES
    """
    assert len(names) == len(diffselfiles) == len(set(names)) > 0
    assert os.path.splitext(plotfile)[1].lower() == '.pdf'

    diffsels = [pd.read_csv(f).assign(name=name) for (name, f) 
            in zip(names, diffselfiles)]
    assert all([set(diffsels[0]['site']) == set(df['site']) for df in 
            diffsels]), "diffselfiles not all for same sites"
    diffsel = pd.concat(diffsels, ignore_index=True)

    ylabel = 'differential selection'
    if diffseltype == 'positive':
        rename = {'positive_diffsel':'above'}
    elif diffseltype == 'total':
        rename = {'positive_diffsel':'above',
                  'negative_diffsel':'below'}
    elif diffseltype == 'max':
        rename = {'max_diffsel':'above'}
    elif diffseltype == 'minmax':
        rename = {'max_diffsel':'above',
                  'min_diffsel':'below'}
    elif diffseltype in ['avgfracsurvive', 'maxfracsurvive']:
        ylabel = 'fraction surviving'
        rename = {diffseltype:'above'}
    else:
        raise ValueError("invalid diffseltype {0}".format(diffseltype))
    diffsel = (diffsel.rename(columns=rename)
                      .melt(id_vars=['site', 'name'], 
                            value_vars=list(rename.values()),
                            value_name='diffsel',
                            var_name='direction')
                      )


    # natural sort by site: https://stackoverflow.com/a/29582718
    diffsel = diffsel.reindex(index=natsort.order_by_index(
            diffsel.index, natsort.index_natsorted(diffsel.site)))#,
            #signed=True)))
    # now some manipulations to make site str while siteindex is int
    diffsel['site'] = diffsel['site'].apply(str)
    diffsel['siteindex'] = pd.Categorical(diffsel['site'],
            diffsel['site'].unique()).codes
    
    ncol = min(maxcol, len(names))
    nrow = math.ceil(len(names) / float(ncol))
    y_lim_min = y_lim - (float(y_lim/30))

    # make name a category to preserve order
    diffsel['name'] = diffsel['name'].astype('category')#, 
            #categories=names)
    
    (xbreaks, xlabels) = dms_tools2.plot.breaksAndLabels(diffsel['siteindex'].unique(), 
            diffsel['site'].unique(), n=6)
    
    #diffsel['highlight'] = diffsel['site'].isin(highlighted_sites) 
    #diffsel['highlight'] =np.where(diffsel['highlight']==True , -0.0092, -1)#-0.0092, -1)#-1)
    diffsel['highlight'] = diffsel['site'].isin(highlighted_sites) 
    diffsel['highlight'] =np.where(diffsel['highlight']==True , y_lim, -1)
    #diffsel['underlay_1'] = 0.16
    diffsel['underlay_1'] = diffsel['site'].isin(underlay_1)
    #worksforunderlaydiffsel['underlay_1'] =np.where(diffsel['underlay_1']==True , -0.003, -1) #0, -1) #last two things are height and bottom
    ymin = float(y_lim) *(float(.005/.15))*float(-1)
    diffsel['underlay_1'] =np.where(diffsel['underlay_1']==True , ymin, -1)
    diffsel['point'] = diffsel['site'].isin(point)
    diffsel['point'] =np.where(diffsel['point']==True , -0.0069, -1)
    #diffsel['point'] =np.where(diffsel['point']==True , diffsel['diffsel'], -1)
    #this will do an overlay diffsel['underlay_2'] =np.where(diffsel['underlay_2']==True , y_lim, -1) #last two things are height and bottom

    
    

    if white_bg:
        p = (ggplot(diffsel, aes(x='siteindex', y='diffsel',
                    color='direction', fill='direction'))             
             + geom_bar(aes(y='highlight'), alpha=0.0, stat="identity", color="#d8d8d8", size=0.5, show_legend=False) #pink was #ECBABA
             + geom_step(size=0.5) #was .3
             + geom_point(aes(y='underlay_1'), fill="#0000e1", shape="|", size=3, show_legend=False) #light blue - #0072b2 old color: "#440154ff"             
             + scale_y_continuous(limits=(ymin, y_lim))#-0.0092, y_lim))
             + xlab('site')
             + ylab(ylabel)
             + scale_x_continuous(breaks=xbreaks, labels=xlabels)
             + scale_color_manual(COLOR_BLIND_PALETTE)
             + scale_fill_manual(COLOR_BLIND_PALETTE)
             + guides(color=False)
             + theme(#panel_background=element_rect(fill='white'),
                     axis_line_x=element_line(color='black'),
                     axis_line_y=element_line(color='black'),
                     axis_title_x=element_blank(),
                     axis_title_y=element_blank(),
                     panel_grid=element_blank(),
                     panel_border=element_blank(),
                     strip_background=element_blank()
                     )
            )
    else:
        p = (ggplot(diffsel, aes(x='siteindex', y='diffsel', color='direction'))
             + geom_step(size=0.4)
             + xlab('site')
             + scale_y_continuous(limits=(0, y_lim))
             + ylab(ylabel)
             + scale_x_continuous(breaks=xbreaks, labels=xlabels)
             + scale_color_manual(COLOR_BLIND_PALETTE)
             + guides(color=False)
             )
    if not ((len(names) == 1) and ((not names[0]) or names[0].isspace())):
        p += facet_wrap('~name', ncol=ncol)
    p += theme(figure_size=(8 * (0.3 + ncol), 1.5 * (0.2 + nrow))) #WAS 4.6  and 1.9
    p.save(plotfile, verbose=False, transparent=True)
plt.close()
```


```python
CD4bs_sites_sub = ['119', '120', '121', "195", "196",'197', '198', '199', '206','207', '208', '209', '275', '276', '277', '278', '279', '280', '281', '282', '283', '304', '305', '306', '307', '308', '309', '310', '311','312','313','314','315', '316', '317', '318', '319', '320', '366', '367', '368', '369', '370', '371', '372', '373', '456', '457', '458', '459', '460', '461']
antibodies = ['1-18', 'VRC01', '3BNC117']

diffseltype = "avgfracsurvive"
plotfile = "./results/fracsurviveaboveavg/CD4bsAb_overlay_median_avgsitefracsurvive.pdf"
ylimit=0.01
plotSiteDiffSel_YLimDiff(antibodies, diffselfiles, plotfile, diffseltype, highlighted_sites = CD4bs_sites_sub, y_lim=ylimit, maxcol=1, white_bg=True)
showPDF(plotfile)
```


![png](analysis_notebook_files/analysis_notebook_107_0.png)


## Make logoplots of key sites
Here, I will plot two sets of sites for each antibody. 

First, I plot the `fraction surviving above average`, `BG505 rescaled preferencs`, and the `frequency in nature` for **structurally defined contact sites**. These figures are not in the paper, but I speculate that some people (like me) will love to examine how these three metrics compare in the epitope of different antibodies. 

Next, I will plot the `fraction surviving above average` for **"sites of interest"** for each epitope. These groups of sites are defined above loosely based on contact sites, sites of significant escape, and prior literature for each antibody that targets that epitope. They are the sites shown in the paper figures. 


```python
epitopelogodir = os.path.join(resultsdir, 'fracsurviveaboveavg/EpitopeLogoplots/')
if not os.path.isdir(epitopelogodir):
    os.mkdir(epitopelogodir)
```


```python
import phydmslib.weblogo
AA_COLORS_FG = phydmslib.weblogo.FunctionalGroupColorMapping()[1]

def siteSubsetGGSeqLogo(logodata, chars, plotfile, width, height,
        yname='', char_colors=AA_COLORS_FG, ylimits=None):
    """Creates one-row logo plot with subset of sites.
    Designed to show logo plot for a subset of sites. This
    is useful when you have data for many sites, but only
    want to look at a few of them. 
    Args:
        `logodata` (pandas DataFrame)
            Contains data to plot. Should have the columns
            `site`, `show`, and a column giving the height
            height of each char in `chars`. Only sites
            where `show` is `True` are shown. Sites are 
            shown in the order they occur in this dataframe,
            with spaces every time there is an interspersed
            site with `show` being `False`. 
        `chars` (list)
            Letters for which we plot heights.
        `plotfile` (str)
            Name of created plot.
        `width` (float)
            Width of plot in inches.
        `height` (float)
            Height of plot in inches.
        `yname` (str)
            If set to a non-empty string, is the y-axis label
            and yticks are drawn.
        `char_colors` (dict)
            Values give color for every character in `chars`.
        `ylimits` (`None` or 2-tuple)
            If not `None`, should give the ylimits for the plot
            as `(ymin, ymax)`
    Here is an example that creates a plot for a subset of
    sites for two characters:
    >>> logodata = pandas.read_csv(io.StringIO(
    ...     '''site show    A    C
    ...        A101 True  0.8  0.2
    ...        N102 True  0.7  0.3
    ...        K103 False 0.1  0.9
    ...        L104 True  0.8  0.2
    ...        S105 True  0.5  0.5
    ...        T106 False 0.2  0.8
    ...        G107 False 0.4  0.6
    ...        L108 True  0.7  0.3'''),
    ...     delim_whitespace=True, index_col=False)
    >>> plotfile = '_siteSubsetGGSeqLogo_test_plot.png'
    >>> siteSubsetGGSeqLogo(logodata,
    ...         chars=['A', 'C'],
    ...         plotfile=plotfile,
    ...         width=3.5, height=2
    ...         )
    >>> os.path.isfile(plotfile)
    True
    Here is the plot created by the code block above:
    .. image:: _static/_siteSubsetGGSeqLogo_test_plot.png
       :width: 55%
       :align: center
    """
    if os.path.isfile(plotfile):
        os.remove(plotfile)

    assert set(chars) <= set(char_colors.keys()), \
            "`char_colors` not defined for all chars"

    expectcol = ['site', 'show'] + chars
    assert set(logodata.columns) >= set(expectcol), \
            "`logodata` needs these column: {0}".format(expectcol)

    assert logodata['show'].any(), "no sites to show"

    # for each consecutive set of rows not to show, keep just one
    logodata = logodata[expectcol]
    logodata['keeprow'] = (
            ((logodata['show']) | 
                (logodata['show'] != logodata['show'].shift(1)))
            )
    logodata = logodata.query('keeprow').reset_index()

    # trim first and last row if they are not to be shown
    if not logodata.iloc[0]['show']:
        logodata = logodata.iloc[1 : ].reset_index()
    if not logodata.iloc[-1]['show']:
        logodata = logodata.iloc[ : -1]

    # set site label to empty and data to zero for rows not to show
    logodata.loc[~logodata['show'], 'site'] = ''
    logodata.loc[~logodata['show'], chars] = 0
    vertlines = logodata.query('~show').index.values + 1

    # generate matrix to plot
    sites = logodata['site']
    matrix = r.matrix(logodata.set_index('site')[chars].values.ravel(),
            ncol=len(sites),
            dimnames=[chars, sites]
            )

    if ylimits is None:
        ylimits = rinterface.NULL
    else:
        ylimits = FloatVector(ylimits)

    # make the plot
    with warnings.catch_warnings():
        warnings.simplefilter(SHOW_WARNINGS)
        _RFUNCS.siteSubsetGGSeqLogo(
                mat=matrix,
                plotfile=plotfile,
                width=width,
                height=height,
                xlabels=list(map(str, sites)),
                vertlines=vertlines,
                yname=yname,
                chars=StrVector(chars),
                char_colors=StrVector([char_colors[x] for x in chars]),
                ylimits=ylimits
                )

    if not os.path.isfile(plotfile):
        raise RuntimeError("failed to create {0}".format(plotfile))
```

### munge the LANL sequence alignment
`HIV1_FLT_2016_env_PRO.fasta` is the filtered web alignment of all HIV-1 sequences, downoaded Jan 3, 2018 from [`LANL`](https://www.hiv.lanl.gov/cgi-bin/NEWALIGN/align.cgi). We convert the alignment with the function below to amino acid frequencies. 
[./data/HXB2_HIV1_FLT_2016_env_PRO_numbering.csv](./data/HXB2_HIV1_FLT_2016_env_PRO_numbering.csv) is a csv file that contains the HXB2 sequence from the alignment (containing the insertions), with the site labeled according to it's site in the alingment. 
["./data/HXB2_alignment_to_HXB2.csv](./data/HXB2_alignment_to_HXB2.csv) contains a conversion between the alignment site and the HXB2 site. note that insertions sites relative to HXB2 are labled with an x rather than the actual notation (i.e. a, b etc). 


```python
#here, I am redefining dms_tools2.prefs.aafreqsFromAlignment such that it gets rid of sequences that:
#are not the same length 
#mask sites that not IUPAC AAs. 
import Bio.SeqIO
def aafreqsFromAlignment_AD(alignmentfile, codon_to_aa,
        ignore_gaps=True, ignore_stop=True):
    """Get amino-acid frequencies at each site in alignment.

    Args:
        `alignmentfile` (str)
            FASTA file with alignment of proteins or coding sequences.
        `codon_to_aa` (bool)
            If `True`, translate codon alignment to amino acids.
        `ignore_gaps` (bool)
            Ignore gaps when calculating frequencies.
        `ignore_stop` (bool)
            Ignore stop codons when calculating frequencies.

    Returns:
        A `pandas.DataFrame` with columns being `site` (1, 2, ...
        numbering) and other columns being amino acids and values
        giving frequencies in alignment.

    >>> with tempfile.NamedTemporaryFile(mode='w') as f:
    ...     x = f.write('>seq1\\n'
    ...                 'ATGGGGCAG\\n'
    ...                 '>seq2\\n'
    ...                 '---AGGCAG\\n'
    ...                 '>seq3\\n'
    ...                 'ATGTGACAG')
    ...     f.flush()
    ...     aafreqs = aafreqsFromAlignment(f.name, codon_to_aa=True)
    >>> aas_counts = ['M', 'G', 'R', 'Q']
    >>> aas_nocounts = [a for a in dms_tools2.AAS if a not in aas_counts]
    >>> (0 == aafreqs[aas_nocounts].values).all()
    True
    >>> expected_counts = pandas.DataFrame.from_items([
    ...         ('site', [1, 2, 3]), ('M', [1.0, 0.0, 0.0]),
    ...         ('G', [0.0, 0.5, 0]), ('R', [0.0, 0.5, 0.0]),
    ...         ('Q', [0.0, 0.0, 1.0])])
    >>> expected_counts.equals(aafreqs[['site'] + aas_counts])
    True
    """
    # read sequences
    seqs = [s.seq for s in Bio.SeqIO.parse(alignmentfile, 'fasta')]
    if codon_to_aa:
        seqs = [s.translate(gap='-', stop_symbol='*') for s in seqs]
    assert seqs, "No sequences"
    seqlen = len(seqs[0])
    #assert seqlen, "sequences have no length" #commented this and line below out. 
    #assert all([seqlen == len(s) for s in seqs]), "seqs not same length"
    # get character sets
    aas = dms_tools2.AAS.copy()
    skipchars = ["#","?","X","$"] #here I added these other characters
    if ignore_gaps:
        skipchars.append('-')
    else:
        aas.append('-')
    if ignore_stop:
        skipchars.append('*')
    else:
        aas.append('*')
    # tally amino-acid frequencies
    aafreqs = dict([(col, [0] * seqlen) for col in aas])
    aafreqs['site'] = list(range(1, seqlen + 1))
    for s in seqs:
        for (r, aa) in enumerate(s):
            if aa in skipchars:
                continue
            else:
                aafreqs[aa][r] += 1
    # convert to dataframe and change counts to freqs
    aafreqs = pd.DataFrame(aafreqs)
    ncounts = aafreqs[aas].sum(axis=1).astype('float')
    for aa in aas:
        aafreqs[aa] = aafreqs[aa] / ncounts
    return aafreqs[['site'] + aas].fillna(0)
```


```python
groupM_alignmentfile = "./data/HIV1_FLT_2016_env_PRO.fasta" 
groupM_df = aafreqsFromAlignment_AD(groupM_alignmentfile , codon_to_aa=False, ignore_gaps=True, ignore_stop=True)
alginment_to_HXB2 = pd.read_csv("./data/HXB2_alignment_to_HXB2.csv")
AnnotatedAlignment_df = alginment_to_HXB2.merge(groupM_df, left_on='alignment_site', right_on='site')
AnnotatedAlignment_df = AnnotatedAlignment_df.drop('site', 1)
AnnotatedAlignment_df = AnnotatedAlignment_df.drop('HXB2_AA', 1)
AnnotatedAlignment_df = AnnotatedAlignment_df.drop('alignment_site', 1)
AnnotatedAlignment_df = AnnotatedAlignment_df.rename(columns={'HXB2_site': 'site'})
AnnotatedAlignment_DropInserts_df = AnnotatedAlignment_df[AnnotatedAlignment_df.site.str.contains("x") == False].copy()
```


```python
import natsort
def EpitopeFracsurviveLogoplot(mutfracsurvivefile, keysites, outfile, max_sitefracsurvive):
    pandadf = pd.read_csv(mutfracsurvivefile)
    df = dms_tools2.diffsel.tidyToWide(pandadf, valuecol='mutfracsurvive')

    #sort sites 
    df = df.reindex(index=natsort.order_by_index(df.index,natsort.index_natsorted(df.site)))
    df['show'] = df['site'].isin(keysites)
    df['site'] = df['wildtype'] + df['site']
    width = len(keysites) / 4
    dms_tools2.rplot.siteSubsetGGSeqLogo(
        logodata=df,
        chars=dms_tools2.AAS,
        plotfile=outfile,
        width= width,
        height=2.5,
        yname='frac surviving',
        ylimits=(0,max_sitefracsurvive),
    )
    
def EpitopePrefLogoplot(preffile, keysites, outfile):
    df = pd.read_csv(preffile)
    df = df.reindex(index=natsort.order_by_index(df.index,natsort.index_natsorted(df.site)))
    df['show'] = df['site'].isin(keysites)
    #df['site'] = df['wildtype'] + df['site']
    width = len(keysites) / 4
    dms_tools2.rplot.siteSubsetGGSeqLogo(
            logodata=df,
            chars=dms_tools2.AAS,
            plotfile=outfile,
            width=width,
            height=2,
            yname='preference',
            ylimits=(0,1.1),
            )
    
def EpitopeNatSeqLogoplot(HXB2annotated_df, keysites, outfile):
    df = HXB2annotated_df
    df = df.reindex(index=natsort.order_by_index(df.index,natsort.index_natsorted(df.site)))
    df['show'] = df['site'].isin(keysites)
    width = len(keysites) / 4
    dms_tools2.rplot.siteSubsetGGSeqLogo(
            logodata=df,
            chars=dms_tools2.AAS,
            plotfile=outfile,
            width=width,
            height=2,
            yname='freq in nature',
            ylimits=(0,1.1),
            )
```


```python
preffile = "./data/BG505-avg-rescaled-prefs_ADrealigned.csv" 
```


```python
epitopelogodir = os.path.join(resultsdir, 'fracsurviveaboveavg/EpitopeLogoplots/SitesOfInterest/')
if not os.path.isdir(epitopelogodir):
    os.mkdir(epitopelogodir)
    
fracsurviveaboveavg_dict_ug = {}
wtDNActrldict = "./results/fracsurviveaboveavg"
fracsurviveaboveavg_dict_ug["VRC01"] = "./results/fracsurviveaboveavg/summary_VRC01-medianmutfracsurvive.csv" 
fracsurviveaboveavg_dict_ug["3BNC117"] = "./results/fracsurviveaboveavg/summary_3BNC117-medianmutfracsurvive.csv" 
fracsurviveaboveavg_dict_ug["1-18"] = "./results/fracsurviveaboveavg/summary_118-medianmutfracsurvive.csv" 


max_sitefracsurvive_dict = {}
max_sitefracsurvive_dict["VRC01"] = "0.2"
max_sitefracsurvive_dict["3BNC117"] = "0.2"
max_sitefracsurvive_dict["1-18"] = "0.2"



natseq_df = AnnotatedAlignment_DropInserts_df.copy()
for ab in antibodies:
    epitopesites = CD4bs_sites_sub #epitope_sites_of_interest[ab]
    print(ab)
    max_sitefracsurvive = max_sitefracsurvive_dict[ab]
    fracsurvivefile= fracsurviveaboveavg_dict_ug[ab]
    fs_outfile = "{0}/{1}_epitope_fracsurviveaboveavg.pdf".format(epitopelogodir, ab)
    pref_outfile = "{0}/{1}_epitope_prefs.pdf".format(epitopelogodir, ab)
    natseq_outfile = "{0}/{1}_NatSeq.pdf".format(epitopelogodir, ab)
    EpitopeFracsurviveLogoplot("{0}".format(fracsurvivefile), epitopesites, fs_outfile, max_sitefracsurvive)
    EpitopePrefLogoplot(preffile, epitopesites, pref_outfile)
    EpitopeNatSeqLogoplot(natseq_df, epitopesites, natseq_outfile)
    print(ab)
    print("Here is the fraction surviving {0} at sites of interest".format(ab, ab))
    pdflist = []
    showPDF(fs_outfile)

```

    1-18
    1-18
    Here is the fraction surviving 1-18 at sites of interest



![png](analysis_notebook_files/analysis_notebook_116_1.png)


    VRC01
    VRC01
    Here is the fraction surviving VRC01 at sites of interest



![png](analysis_notebook_files/analysis_notebook_116_3.png)


    3BNC117
    3BNC117
    Here is the fraction surviving 3BNC117 at sites of interest



![png](analysis_notebook_files/analysis_notebook_116_5.png)



```python
escapedir = './results/fracsurviveaboveavg/escapability_plots/'
if not os.path.isdir(escapedir):
    os.mkdir(escapedir)
```

Now, lets compare the largest effect size mutations across antibodies. For each antibody, I will just plot the effect size of the top 100 mutations...

First, I need to get the median excess mut frac survive files from [Dingens et al Immunity 2019](https://www.cell.com/immunity/fulltext/S1074-7613(18)30565-X) for 10-1074 and pooled 10-1074/3BNC117. I did not reanalyze these samples start to finish because they are not included in other samples. 


```python
import requests
 
url = 'https://raw.githubusercontent.com/jbloomlab/EnvsAntigenicAtlas/master/results/fracsurviveaboveavg/concavg_wtDNA_ctrl/summary_101074-medianmutfracsurvive.csv'
ab101074file = requests.get(url)
open('./results/fracsurviveaboveavg/summary_101074-medianmutfracsurvive.csv', 'wb').write(ab101074file.content)

url = 'https://raw.githubusercontent.com/jbloomlab/EnvsAntigenicAtlas/master/results/fracsurviveaboveavg/concavg_wtDNA_ctrl/summary_3BN-1074-pool-medianmutfracsurvive.csv'
ab3BN101074file = requests.get(url)
open('./results/fracsurviveaboveavg/summary_3BN-1074-pool-medianmutfracsurvive.csv', 'wb').write(ab3BN101074file.content)

```




    183206




```python
colorcycle = ["#279846", "#293381", "#00a2d9",
             #"#7f2d19", "#c04526",
             #"#eb5e18", "#f5ad8a",
             #"#1da19d", "#53dfdb",
             "#e05400", "#808080"]
phiaboveavg_sub_tuple_mut = [('1-18', './results/fracsurviveaboveavg/summary_118-medianmutfracsurvive.csv'),
    ('VRC01', './results/fracsurviveaboveavg/summary_VRC01-medianmutfracsurvive.csv'),
    ('3BNC117', './results/fracsurviveaboveavg/summary_3BNC117-medianmutfracsurvive.csv'),
    ('101074', './results/fracsurviveaboveavg/summary_101074-medianmutfracsurvive.csv'),
    ('3BNC117/101074', './results/fracsurviveaboveavg/summary_3BN-1074-pool-medianmutfracsurvive.csv')]

    

def PlotRankedMutationsFrac(labels_files, outfileprefix, include_stdev=False, rank_lims=False, y_lims=False, title=False,  
                       colorcycle=colorcycle, alpha=0.6, make_legend=True, figsize=(5,4), ylabel='fraction surviving', convert_to_enrichment=False):
    '''labels_files is a list of tuples of (label, diffsel_file).
    Those diffsel_files must have a stdev column if include_stdev=True.
    
    To keep files friendly to existing dms_tools programs like dms_merge and dms_logoplot, phi files use the header `diffsel` for phi (but are named appropriately)
    
    the prefix in outfile will be saved in the plots directory with .pdf added.
    rank_lims sets the x-axis limit (mutation ranks)'''
    
    fig = plt.figure(figsize=figsize)
    
    for i, (difflabel, difffile) in enumerate(labels_files):
        df = pd.read_csv(difffile).dropna()
        if not convert_to_enrichment:
            plt.plot(df['mutfracsurvive'], marker='.', label = difflabel, linewidth=0.6, mew=0, color=colorcycle[i]) # i think this plots rank technically as 0, 1, 2 instead of 1, 2, 3...
        else:
            plt.plot(2**df['mutfracsurvive'], marker='.', label = difflabel, linewidth=0.6, mew=0, color=colorcycle[i]) # i think this plots rank technically as 0, 1, 2 instead of 1, 2, 3...

    plt.xlabel('mutation rank')
    plt.ylabel(ylabel) # default is differential selection, but can change to plot rank-ordered phis.
    
    spineOffset = {'left': 4, 'bottom': 4}    
    [spine.set_position(('outward',spineOffset[loc])) if loc in ['left','bottom'] else spine.set_color('none') for loc, spine in plt.gca().spines.items() ] 
    plt.gca().tick_params(axis='x', direction='out')
    plt.gca().tick_params(axis='y', direction='out')
    plt.gca().get_xaxis().tick_bottom()
    plt.gca().get_yaxis().tick_left()
    
    if rank_lims:
        plt.xlim(rank_lims)
    if y_lims:
        plt.ylim(y_lims)
    if title:
        plt.title(title)
    if make_legend:
        legend = plt.legend(fontsize=10.5, fancybox = True)
        
    if include_stdev:
        for i, (difflabel, difffile) in enumerate(labels_files):
            df = pd.read_csv(difffile).dropna()
            plt.gca().errorbar(range(0, len(df.index)),
                               df['diffsel'],
                               yerr=df['stdev'], 
                               marker=None, color=colorcycle[i], 
                               alpha=alpha, capsize=0, elinewidth=0.9)
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()
    



    
tempoutdir = './results/fracsurviveaboveavg/escapability_plots/'
if not os.path.isdir(tempoutdir):
    os.mkdir(tempoutdir)

outfile = "./results/fracsurviveaboveavg/escapability_plots/BG505_All_fracsurvive.pdf"
print("all mutations, ranked phi:")                           
PlotRankedMutationsFrac(phiaboveavg_sub_tuple_mut, outfile, title='all mutations', make_legend=True, y_lims=(0, .4), figsize=(9,4.5), ylabel='fracsurvive aboveavg', colorcycle=colorcycle)
showPDF([outfile], width=700)

outfile = "./results/fracsurviveaboveavg/escapability_plots/BG505_top1000_fracsurvive.pdf"
print("top 1000, ranked phi:")                           
PlotRankedMutationsFrac(phiaboveavg_sub_tuple_mut, outfile, rank_lims=(-10, 1000), y_lims=(0, .4), title='top 1000', make_legend=True, figsize=(9,4.5), ylabel='fracsurvive aboveavg', colorcycle=colorcycle)
#showPDF([outfile], width=700)

outfile = "./results/fracsurviveaboveavg/escapability_plots/BG505_top100_fracsurvive.pdf"
print("top 100, ranked phi:")                           
PlotRankedMutationsFrac(phiaboveavg_sub_tuple_mut, outfile, rank_lims=(-2, 40), y_lims=(0, .4), title=' ', make_legend=True, figsize=(6,4), ylabel='fracsurvive aboveavg', colorcycle=colorcycle)
showPDF([outfile], width=700)
```

    all mutations, ranked phi:



![png](analysis_notebook_files/analysis_notebook_120_1.png)


    top 1000, ranked phi:
    top 100, ranked phi:



![png](analysis_notebook_files/analysis_notebook_120_3.png)


Now, lets color that plot by sequence accessible vs non accessible mutations


```python
from dms_tools2.utils import codonEvolAccessibility
#below functions edited from Shirleen Soh, https://github.com/jbloomlab/computational_notebooks/blob/master/yqsoh/2018/PB2-DMS-full-A549-CCL141/AnalyzeAdaptiveSites.ipynb

def calc_accessibility(row):
    if row['min Subst']==0:
        return '0'
    elif row['min Subst']<=1.1: # Picked this cut off to accommodate values slightly >1 when considering many sequences
        return '1'
    else:
        return '>1'
def seqsToCodonAcc(filename, seqlen):
    allowed_chars = set('ATCG')
    seqs = []
    for seq_record in Bio.SeqIO.parse(filename, 'fasta'):
        if len(seq_record.seq)==seqlen:
            if set(seq_record.seq).issubset(allowed_chars):
                seqs.append(str(seq_record.seq))
    #         else:
    #             print('Invalid bases:', seq_record.id, len(seq_record.seq))
        else:
            print('Not full length:', seq_record.id, len(seq_record.seq))
    accessibilitydf = (codonEvolAccessibility(seqs)
                     .melt(id_vars='site', value_vars=dms_tools2.AAS_WITHSTOP, 
                           var_name='toAA', value_name='min Subst')
                    )
    accessibilitydf['accessibility'] = accessibilitydf.apply(lambda row: calc_accessibility(row), axis=1)
    return accessibilitydf

accessibilityBG505 = seqsToCodonAcc('./data/BG505.W6.C2.T332N_env.fasta', 2583)

#now, I need to convert to HXB2 numbering, and be able to merge into mutfracsurvive df. 
convert_df = pd.read_csv("./results/HXB2_numbering/BG505_to_HXB2.csv")
convert_df.drop('N-glycan', axis=1, inplace=True)
accessibilityBG505_convert = accessibilityBG505.merge(convert_df, left_on = "site", right_on = "original")
accessibilityBG505_convert = accessibilityBG505_convert.rename(columns={'new': 'HXB2_site'})
accessibilityBG505_convert["mut"] = accessibilityBG505_convert["wildtype"].astype(str) + accessibilityBG505_convert["HXB2_site"].astype(str) + accessibilityBG505_convert["toAA"].astype(str) 
```


```python
def PlotRankedMutationsFracAccessible(labels_files, outfileprefix, include_stdev=False, rank_lims=False, y_lims=False, title=False,  
                       colorcycle=colorcycle, alpha=0.6, make_legend=True, figsize=(5,4), ylabel='fraction surviving', convert_to_enrichment=False):
    '''labels_files is a list of tuples of (label, diffsel_file).
    Those diffsel_files must have a stdev column if include_stdev=True.
    
    To keep files friendly to existing dms_tools programs like dms_merge and dms_logoplot, phi files use the header `diffsel` for phi (but are named appropriately)
    
    the prefix in outfile will be saved in the plots directory with .pdf added.
    rank_lims sets the x-axis limit (mutation ranks)'''
    
    fig = plt.figure(figsize=figsize)
    
    for i, (difflabel, difffile) in enumerate(labels_files):
        df = pd.read_csv(difffile).dropna()
        if not convert_to_enrichment:
            #merge in an drop!
            df["mut"] = df["wildtype"].astype(str) + df["site"].astype(str) + df["mutation"].astype(str) 
            accessibilityBG505_convert_copy = accessibilityBG505_convert.copy()
            df = df.merge(accessibilityBG505_convert_copy, left_on = "mut", right_on = "mut")
            
            #only plot if accesible! I will simply drop if not accessible
            dfaccess = df[df.accessibility != ">1"]
            dfnonaccess = df[df.accessibility == ">1"]
            
            plt.plot(df['mutfracsurvive'], marker=None, label = difflabel, linewidth=0.6, color=colorcycle[i]) # i think this plots rank technically as 0, 1, 2 instead of 1, 2, 3...
    
            plt.plot(dfnonaccess['mutfracsurvive'], marker='o', markersize=4, fillstyle="none", markerfacecolor="white", label = difflabel, linewidth=0, color=colorcycle[i]) # i think this plots rank technically as 0, 1, 2 instead of 1, 2, 3...
            plt.plot(dfaccess['mutfracsurvive'], marker='o', markersize=4,  label = difflabel, linewidth=0, color=colorcycle[i]) # i think this plots rank technically as 0, 1, 2 instead of 1, 2, 3...

            
        else:
            #merge in an drop!
            df["mut"] = df["wildtype"].astype(str) + df["site"].astype(str) + df["mutation"].astype(str) 
            accessibilityBG505_convert_copy = accessibilityBG505_convert.copy()
            df = df.merge(accessibilityBG505_convert_copy, left_on = "mut", right_on = "mut")
            #only plot if accesible! I will simply drop if not accessible
            df = df[df.accessibility != ">1"]
            print(df)
            plt.plot(2**df['mutfracsurvive'], marker='.', label = difflabel, linewidth=0.6, mew=0, color=colorcycle[i]) # i think this plots rank technically as 0, 1, 2 instead of 1, 2, 3...

    plt.xlabel('mutation rank')
    plt.ylabel(ylabel) # default is differential selection, but can change to plot rank-ordered phis.
    
    spineOffset = {'left': 4, 'bottom': 4}    
    [spine.set_position(('outward',spineOffset[loc])) if loc in ['left','bottom'] else spine.set_color('none') for loc, spine in plt.gca().spines.items() ] 
    plt.gca().tick_params(axis='x', direction='out')
    plt.gca().tick_params(axis='y', direction='out')
    plt.gca().get_xaxis().tick_bottom()
    plt.gca().get_yaxis().tick_left()
    
    if rank_lims:
        plt.xlim(rank_lims)
    if y_lims:
        plt.ylim(y_lims)
    if title:
        plt.title(title)
    if make_legend:
        legend = plt.legend(fontsize=10.5, fancybox = False)
        
    if include_stdev:
        for i, (difflabel, difffile) in enumerate(labels_files):
            df = pd.read_csv(difffile).dropna()
            plt.gca().errorbar(range(0, len(df.index)),
                               df['diffsel'],
                               yerr=df['stdev'], 
                               marker=None, color=colorcycle[i], 
                               alpha=alpha, capsize=0, elinewidth=0.9)

    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()
```


```python
tempoutdir = './results/fracsurviveaboveavg/escapability_plots/Accessible/'
if not os.path.isdir(tempoutdir):
    os.mkdir(tempoutdir)


outfile = "./results/fracsurviveaboveavg/escapability_plots/Accessible/AccOnly_BG505_top100_fracsurvive.pdf"
print("top 40, ranked phi:")                           
PlotRankedMutationsFracAccessible(phiaboveavg_sub_tuple_mut, outfile, rank_lims=(-2, 40), y_lims=(0, .15), title='', make_legend=False, figsize=(6,4), ylabel='excess fracsurvive', colorcycle=colorcycle)
showPDF([outfile], width=700)
```

    top 40, ranked phi:



![png](analysis_notebook_files/analysis_notebook_124_1.png)



```python
import pandas as pd
from colour import Color
import os

pymoldir = './results/fracsurviveaboveavg/pymol/'
if not os.path.isdir(pymoldir):
    os.mkdir(pymoldir)

def MapFracSurvtoPDB(infile, 
                     scriptfile, 
                     colors = ['#fafafa', '#ff0000'], 
                     map_type = 'site_fracsurv', 
                     restrict_to_chain = False, 
                     script_preamble = None,
                     script_postamble = None,
                     abname = None):
    '''Writes a colormapping script to be run in pymol; the colormapping is based on fracsurvive 
    to color a structure'''
    df = pd.read_csv(infile)
    df = df.dropna()
    column_names = list(df)
    
    # establish the color spectrum in hex and rgb.
    n_subdivisions = 500 # the color spectrum will be divided into this many discrete colors
    color1 = Color(colors[0])
    color2 = Color(colors[1])
    rgb_spectrum = [c.rgb for c in color1.range_to(color2, n_subdivisions)]
    rgb_spectrum_dict = dict([(i, rgb_spectrum[i]) for i in range(len(rgb_spectrum))])
    
    if map_type == 'site_fracsurv':
        assert 'avgfracsurvive' in column_names
        min_avg = df.min()['avgfracsurvive']  
        max_avg = df.max()['avgfracsurvive']  # the min and max will be mapped to color1 and color2, respectively
        range_avg = max_avg - min_avg
        df['colorindex'] =  (df.avgfracsurvive - min_avg)/range_avg*(n_subdivisions-1)
        
    elif map_type == 'max_fracsurv':
        assert 'maxfracsurvive' in column_names
        min_val = df.min()['maxfracsurvive']  
        max_val = df.max()['maxfracsurvive']  # the min and max will be mapped to color1 and color2, respectively
        range_val = max_val - min_val
        df['colorindex'] =  (df.maxfracsurvive - min_val)/range_val*(n_subdivisions-1)
    
    df['colorindex'] = df['colorindex'].astype(int) # round to nearest index
    df['rgb'] = df['colorindex'].map(rgb_spectrum_dict)        
    site_color_mapping = pd.concat([df['site'], df['rgb']], axis=1)
    
    # write out the script to *scriptfile*:
    f = open(scriptfile, 'w')
    
    if script_preamble:
        preamblef = open(script_preamble, 'r')
        for line in preamblef.readlines():
            f.write(line)
        f.write('\n\n')
        preamblef.close()
    
    for i in range(len(df.index)):
        rgblist = [min(1, c) for c in site_color_mapping.iloc[i]['rgb']]
        r = site_color_mapping.iloc[i]['site']
        
        f.write("cmd.set_color(\'color{0}\', \'{1}\')\n".format(r, rgblist))
        f.write("cmd.color(\'color{0}\', \'resi {0}\')\n".format(r))
    if script_postamble:
        postamblef = open(script_postamble, 'r')
        f.write('abname = "{0}"'.format(abname))
        for line in postamblef.readlines():
            f.write(line)
        f.write('\n\n')
        postamblef.close()
    f.close()
```


```python
for ab in groups:
    print (ab)
    MapFracSurvtoPDB("./results/fracsurviveaboveavg/summary_{0}-mediansitefracsurvive.csv".format(ab), 
                 '{0}/{1}-sitefracsurvive_avgsitefracsurvive.py'.format(pymoldir, ab), 
                 map_type = 'site_fracsurv',
                 abname = ab)
    MapFracSurvtoPDB("./results/fracsurviveaboveavg/summary_{0}-mediansitefracsurvive.csv".format(ab), 
                 '{0}/{1}-sitefracsurvive_maxsitefracsurvive.py'.format(pymoldir, ab), 
                 map_type = 'max_fracsurv',
                 abname = ab)
```

    118
    3BNC117
    VRC01

