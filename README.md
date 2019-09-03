# CWL Definition files for Bioinformatics NGS data analysis

## Requirements
 
 * Python 3.x
 * CWLtools
 * Docker  

## Tools 
Defined in folder *tools*

## Workflows
Defined in *workflows*

## Bioconda environment

We created some requirement files for installing in **conda** environment. 

Please, see https://bioconda.github.io/ for Bioconda instalation.

### Creating a *~/.condarc* file

Create a *.condarc* in your home directory and add:

    channels:
        - conda-forge
        - bioconda
        - defaults
    ssl_verify: true
    
### Bioconda environment for RNA-Seq

#### Creating a *conda* environment

    conda create -n rnaseq
    
#### Installing the *Bioconda* packages

    conda install -n rnaseq --file https://raw.githubusercontent.com/ncbi/cwl-ngs-workflows-cbb/master/requirements/conda-rnaseq.txt

Check the https://github.com/ncbi/cwl-ngs-workflows-cbb/tree/master/requirements folder 
for more environments files. 

#### Activating the *rnaseq* env

    source activate rnaseq
    
Then, all the required programs and tools will be available in the user **$PATH**

### Bioconda environment for ChIP-Seq

#### Creating a *conda* environment

    conda create -n chipseq
    
#### Installing the *Bioconda* packages

    conda install -n chipseq --file https://raw.githubusercontent.com/ncbi/cwl-ngs-workflows-cbb/master/requirements/conda-chipseq.txt

Check the https://github.com/ncbi/cwl-ngs-workflows-cbb/tree/master/requirements folder 
for more environments files. 

#### Activating the *chipseq* env

    source activate chipseq
    
Then, all the required programs and tools will be available in the user **$PATH**

# Public Domain notice

National Center for Biotechnology Information.

This software is a "United States Government Work" under the terms of the United States
Copyright Act. It was written as part of the authors' official duties as United States
Government employees and thus cannot be copyrighted. This software is freely available
to the public for use. The National Library of Medicine and the U.S. Government have not
 placed any restriction on its use or reproduction.

Although all reasonable efforts have been taken to ensure the accuracy and reliability
of the software and data, the NLM and the U.S. Government do not and cannot warrant the
performance or results that may be obtained by using this software or data. The NLM and
the U.S. Government disclaim all warranties, express or implied, including warranties
of performance, merchantability or fitness for any particular purpose.

Please cite NCBI in any work or product based on this material.
