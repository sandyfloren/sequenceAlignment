# sequenceAlignment
Tools for personalizing Multiple Sequence Alignments


#### Update: June 2020
As of now the NCBI SNP GeneView page no longer supports hg37. 
As a large part of this project required data to be scraped from this specific page, and since 23andMe data still uses the hg37 build, 
it will be necessary to completely reimplement this part of the code in order to ensuree compatibility with NCBI's Variation Viewer page,
which is currently the only place where human genetic variation on this older assembly is still reported by NCBI.

#### Usage

    python3 sequenceAligment.py -gI 79068 -acc NM_001080432.2 -f phylip
    exampleGenomeFile.txt outputFile.txt
    
    python3 sequenceAlignment.py
    
### Important:
This program requires an email address and NCBI API key. To obtain a key and make an account,
visit: https://www.ncbi.nlm.nih.gov/account/, then change the placeholders for Entrez.api_key and Entrez.email, near the top of 
sequenceAlignment.py to your ownn personal API key and email address.
    
## sequenceAlignment
Currently, only genomic 23andMe tab-delimited files are suitable input data.
23andMe is still using human genome build 37, which makes it necessary to
convert NCBI genomic SNP locations data from hg38 to hg37. This is done with
the MyVariant tool: https://pypi.org/project/myvariant/

Other important information is retrieved from NCBI databases using
Biopython's Entrez, and sequence manipulation and alignment uses other Bio
modules, including SeqIO, AlignIO, etc: https://biopython.org,
Copyright 1999-2003 by Jeffrey Chang.

All the SNPs in a given mRNA sequence can be found on dbSNP's GeneView page,
which provides a batch query service that seems to be broken. Because of
this, it was necessary to scrape the data directly using Selenium. Thus in
order to obtain the data it is important to install the latest version of
the Google Chrome Webdriver and browser:
https://chromedriver.chromium.org/getting-started
https://www.google.com/chrome/
