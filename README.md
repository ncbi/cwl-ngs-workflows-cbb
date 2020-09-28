# CWL Definition files for Bioinformatics NGS data analysis

## Requirements
 
 * Python 3.x
 * [CWLtool](https://github.com/common-workflow-language/cwltool)
 * Node.JS
 * Conda ([Bioconda](https://bioconda.github.io/)) or Docker  

## CWL Tools 
Defined in folder *tools*

## CWL Workflows
Defined in *workflows*

## Conda execution

Install and configure Conda and [Bioconda](https://bioconda.github.io/)

### Creating a *~/.condarc* file

Create a *.condarc* in your home directory and add:

    channels:
        - conda-forge
        - bioconda
        - defaults
    ssl_verify: true
    
### Bioconda environment cwltool

#### Creating a *conda* environment

    conda create -n cwltool
    
#### Installing the packages

    conda install -n rnaseq cwltool nodejs

#### Activating the *cwltool* env

    conda activate cwltool
    
Then, all the required programs and tools will be available in the user **$PATH**

#### Running cwltool

Using Conda as runtime environment the **cwltool** should be executed with the following examples:

    cwltool --no-container --beta-conda-dependencies --beta-dependencies-directory /path-to-directory/storing-conda-envs

## Docker execution

Execute **cwltool** with these options if Docker is installed and user has permissions to pull images:

    cwltool  --no-read-only --beta-use-biocontainers

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
