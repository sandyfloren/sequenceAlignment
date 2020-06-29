#!/usr/bin/env python3
# Name: Sandy Floren (afloren)
# Group Members: None
"""This module provides tools for personalizing Multiple Sequence Alignments.

Currently, only genomic 23andMe tab-delimited files are suitable input data.
23andMe is still using human genome build 37, which makes it necessary to
convert NCBI genomic SNP locations data from hg38 to hg37. This is done with
the MyVariant tool: https://pypi.org/project/myvariant/

Other important information is retrieved from NCBI databases using
Biopython's Entrez, and sequence manipulation and alignment uses other Bio
modules, including SeqIO, AlignIO, etc: https://biopython.org,
Copyright 1999-2003 by Jeffrey Chang.



This program requires an email address and
NCBI API key. Please use your own. To obtain a key and make an account,
visit: https://www.ncbi.nlm.nih.gov/account/

All the SNPs in a given mRNA sequence can be found on dbSNP's GeneView page,
which provides a batch query service that seems to be broken. Because of
this, it was necessary to scrape the data directly using Selenium. Thus in
order to obtain the data it is important to install the latest version of
the Google Chrome Webdriver and browser:
https://chromedriver.chromium.org/getting-started
https://www.google.com/chrome/



Classes:
    GenomeReader: Defines objects to read 23andMe files.
        -Adapted from FastAreader by David Bernick for BME160, UCSC, Feb 2020
    SequenceAlignment: Define objects to create a personalized alignments.
    CommandLine: Handle the command line, usage and help requests.
        -Adapted from CommandLine by David Bernick for BME160, UCSC, Mar 2020

Methods:
    getChromosome: Find the chromosome number for a given gene ID.
    getAccessionsFromGene: Find all accessions in GeneView for a gene.
    getGeneFromAccession: Find the gene ID from an accesssion number.
    rsidsToHG37Positions: Find the positions on hg37 from a list of rsids.
    geneIDtoHG37Positions: Find the positions on hg37 of all SNPS in a gene.

Example command line usage:
    python3 sequenceAligment.py -gI 79068 -acc NM_001080432.2 -f phylip
        exampleGenomeFile.txt outputFile.txt

    python3 sequenceAlignment.py
"""
import sys
import pandas as pd
from Bio import Entrez, SeqIO, AlignIO, Align
from Bio.Alphabet import IUPAC
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from urllib.error import HTTPError

# this try/except block was adapted from Bio.NaiveBayes, Copyright 2000 by
# Jeffrey Chang
try:
    from selenium import webdriver
except ImportError:
    raise ImportError("sequenceAlignment.py requires installation of the "
                      "Selenium module.")
try:
    import myvariant
except ImportError:
    raise ImportError("sequenceAlignment.py requires installation of the "
                      "myvariant module.")

# change these for personal use

Entrez.api_key = 'YOUR KEY HERE'
Entrez.email = 'YOUR EMAIL HERE'
Entrez.tool = "sequenceAlignment.py"


def getChromosome(geneID):
    """Find what chromosome a gene is on.

        args:
            geneID (str or int): the gene ID number

        returns:
            chrom (int): the chromosome number
    """
    handle = Entrez.efetch(db='gene', id=geneID, rettype='xml')
    result = Entrez.read(handle)
    chrom = int(
        result[0]['Entrezgene_locus'][0]['Gene-commentary_accession'][-2:])
    return chrom


def getAccessionsFromGene(geneID):
    """Return a list of valid mRNA accession numbers from GeneView.

        args:
            geneID (str or int): the gene ID number

        returns:
            accessions (list): the valid accession numbers

    Does the same thing as SequenceAlignment.validAccessions(), but
    without creating a SequenceAlignment object.
    """
    url = f'https://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?locusId={geneID}'
    driver = webdriver.Chrome("/usr/local/bin/chromedriver")
    driver.get(url)

    mrnaElements = driver.find_elements_by_class_name('gm_mrna')
    accessions = list(
        elem.text for elem in mrnaElements)

    driver.close()

    return accessions


def getGeneFromAccession(accession):
    """Find the gene ID for a given nucleotide accession number.

        args:
            accession (str): the accession number

        returns:
            geneID (str): the gene ID
    """
    try:
        handle = Entrez.efetch(db='nuccore', id=accession, rettype='gb',
                               retmode='xml')
        result = Entrez.read(handle)
    except HTTPError:
        raise ValueError("Invalid accession number.")

    # find the geneID from the EFetch XML result
    for dictElement in result[0]['GBSeq_feature-table']:
        if dictElement['GBFeature_key'] == 'CDS':
            for elem in dictElement['GBFeature_quals']:
                if elem['GBQualifier_name'] == 'db_xref':
                    if 'GeneID' in elem['GBQualifier_value']:
                        stringElement = elem['GBQualifier_value']
                        geneID = stringElement[stringElement.find(':') + 1:]
                        break
    return geneID


def geneLocationHG37(geneID):
    """Find the genomic position of a gene on hg37."""
    url = f'https://www.ncbi.nlm.nih.gov/gene/?term={geneID}'
    driver = webdriver.Chrome("/usr/local/bin/chromedriver")
    driver.get(url)
    locXpath = '//*[@id="ui-ncbigrid-11"]/tbody/tr[2]/td[5]'
    chromXpath = '//*[@id="ui-ncbigrid-11"]/tbody/tr[2]/td[4]'
    text = driver.find_element_by_xpath(locXpath).text
    text = text[text.find('(') + 1: text.find(')') + 1]
    text = text[0:text.find(',')]
    location = text.split('..')

    chrom = driver.find_element_by_xpath(chromXpath).text
    driver.close()
    startPos = location[0]
    endPos = location[1]
    return chrom, startPos, endPos


def getHG37PositionsInRange(chromosome, startPos, endPos):
    """Return a DataFrame containing hg37 positions for all rsids in a range.

        args:
            chromosome (int or str): the chromosome number
            startPos (int or str): the start position on the chromosome
            endPos (int or str): the end position on the chromosome

        returns:
            df (DataFrame): all the rsids found in the genomic range
                between startPos and endPos, indexed by rsid
            chromosome (int or str): the chromosome number
    """
    queryString = f'chr{chromosome}:{startPos}-{endPos}'

    mv = myvariant.MyVariantInfo()
    gen = mv.query(queryString, scopes='dbsnp.rsid',
                   fields='dbsnp.rsid, dbsnp.hg19.start', fetch_all=True,
                   assembly='hg37')

    rsids = {}
    for row in gen:
        try:
            rsid = (row['dbsnp']['rsid'])
            start = (row['dbsnp']['hg19']['start'])
            rsids[rsid] = start
        except KeyError:
            continue
    df = pd.DataFrame.from_dict(rsids, orient='index')
    return df, chromosome


def rsidsToHG37Positions(rsidList):
    """Return a DataFrame containing hg37 positions for a list of rsids.

        args:
            rsidList (list of str): the rsids

        returns:
            df (DataFrame): all the rsids found in the genomic range
                between startPos and endPos, indexed by rsid
    """

    mv = myvariant.MyVariantInfo()
    gen = mv.querymany(rsidList, scopes='dbsnp.rsid',
                       fields='dbsnp.rsid, dbsnp.hg19.start', fetch_all=True,
                       assembly='hg37')

    rsids = {}
    for row in gen:
        try:
            rsid = (row['dbsnp']['rsid'])
            start = (row['dbsnp']['hg19']['start'])
            rsids[rsid] = start
        except KeyError:
            continue
    df = pd.DataFrame.from_dict(rsids, orient='index')
    return df


def geneIDtoHG37Positions(geneID):
    """Return a DataFrame containing hg37 positions for all rsids in a gene.

           args:
               geneID (int or str): the geneID

           returns:
               df (DataFrame): all the rsids found in the genomic range
                   between startPos and endPos, indexed by rsid
               chromosome (int or str): the chromosome number
       """
    return getHG37PositionsInRange(geneLocationHG37(geneID))


class GenomeReader:
    """Define objects to read genome files.

        instantiation:
            thisReader = GenomeReader ('genome_John_Cleese_v5_full.txt')

        usage:
            for rsid, genotype, chrom, pos in thisReader.readGenome():
                print (rsid, genotype, chrom, pos)
    """

    def __init__(self, filename=''):
        """Contructor: save attribute filename."""
        self.filename = filename

    def openFile(self):
        """Handle file opens, allowing stdin."""
        if self.filename is '':
            return sys.stdin
        else:
            return open(self.filename)

    def readGenome(self):
        """Parse the data from a 23andMe genome file."""
        count = 0
        with self.openFile() as handle:
            line = handle.readline()
            while line:
                # get rid of comments at top of file and internal 23andMe
                # markers.
                if not line.startswith('#') and not line.startswith('i'):
                    entry = line.split()
                    rsid = entry[0]
                    location = entry[2]
                    chromosome = entry[1]
                    genotype = entry[3]  # save genotype
                    count += 1
                    yield rsid, genotype, chromosome, location
                line = handle.readline()


class SequenceAligment:
    """Define objects to create a personal alignment to the reference genome.

        instantiation:
            sa = SequenceAlignment('79068,' 'NM_001080432.2',
            genome_file='genome_John_Cleese_v5_full.txt', indels=True)

        usage:
            from Bio import AlignIO
            AlignIO.write([sa.get_alignment()], handle=outfile.aln,
                  format='clustal')

    Given a gene ID, an accession number, and a file of SNP data,
    a SequenceAlignment object will scrape the dbSNP GeneView webpage using a
    Chrome webdriver from the selenium module. In the future, I would like to
    find a faster way to collect SNP data for a specific transcript sequence,
    but as of yet I have not found one.

    The SNP data scraped from GeneView are represented as a pandas DataFrame
    object, which is subsequently joined with another DataFrame containing
    all the corresponding rsids (and their genotypes) from the input genome
    file. After dropping NaN values, what is left is a much smaller DataFrame
    containing only SNPs that are present in both the sequence of interest
    and the genome file.

    This small subset of SNPs is then reduced to an even smaller DataFrame
    containing only non-wildtype (mutant) alleles. A Seq object is then
    created for the reference sequence, from which a MutableSeq object is
    then created and modified according to the mutations found in the genome
    file.

    Finally, the two Seqs are used to create a MultipleSeqAlignment object,
    comparing the personalized sequence to the reference sequence.
    """
    # allowable accession number prefixes
    NUCLEOTIDE_PREFIXES = {'NM', 'NR', 'XM', 'XR'}
    PROTEIN_PREFIXES = {'AP', 'NP', 'YP', 'XP', 'WP'}

    # ambiguous bases. Used to convert a heterozygous genotype to a Seq object,
    # e.g. an 'AG' genotype would be represented as 'R'
    AMBIGUOUS_DNA_LETTERS = {
        'AA': 'A',
        'CC': 'C',
        'TT': 'T',
        'GG': 'G',
        'AC': 'M',
        'CA': 'M',
        'AG': 'R',
        'GA': 'R',
        'AT': 'W',
        'TA': 'W',
        'CG': 'S',
        'GC': 'S',
        'CT': 'Y',
        'TC': 'Y',
        'GT': 'K',
        'TG': 'K'
    }

    def __init__(self, geneID, accession, genomeFile='', indels=True):
        """Initialize a new SequenceAnalysis object.
            Args:
                geneID (str or int): the gene ID number
                accession (str): the accession number
                genomeFile (str): the output filename/path. Defaults to ""
                indels (bool): whether or not to include insertions/deletions
                    in the alignment. Defaults to True

        """
        self.geneID = geneID
        self.accessionNumber = accession
        self.genomeFile = genomeFile

        # if indels = True, treat insertions/deletions as
        # homozygous. Heterozygous indels are difficult to represent with
        # ambiguous bases in the way that, for example, a CT genotype can be
        # represented with a 'Y'.
        self.indels = indels
        self.chromosome = getChromosome(geneID)
        self.driver = self.webdriver()
        self.data = self.getDataFrame()
        self.referenceSeqRecord = self.getReferenceSeqRecord()
        self.driver.close()

    def webdriver(self):
        """Return a Google Chrome webdriver object."""
        url = f'https://www.ncbi.nlm.' \
              f'nih.gov/SNP/snp_ref.cgi?locusId={self.geneID} '
        driver = webdriver.Chrome("/usr/local/bin/chromedriver")
        driver.get(url)
        return driver

    def getValidAccessions(self):
        """Return a list of valid mRNA accession numbers from GeneView."""
        mrnaElements = self.driver.find_elements_by_class_name('gm_mrna')
        accessions = list(
            elem.text for elem in mrnaElements)

        return accessions

    def geneViewBodyText(self):
        """Get the SNP data out of the SNP GeneView page."""
        names = self.getValidAccessions()
        index = names.index(self.accessionNumber)
        sequenceXpath = f'/html/body/table[2]/tbody/tr/td[3]/div/table[' \
                        f'2]/tbody/tr/td/table[3]/tbody/tr[{index + 3}] '
        viewButtonXpath = sequenceXpath + '/td[7]/a'

        element = self.driver.find_element_by_xpath(sequenceXpath)

        # click on the sequence of interest if not currently shown
        if 'currently shown' not in element.text:
            viewButton = self.driver.find_element_by_xpath(viewButtonXpath)
            viewButton.click()

        # now access the SNP data:
        body = self.driver.find_element_by_xpath(
            '/html/body/table[2]/tbody/tr/td[3]/div/table['
            '2]/tbody/tr/td/table[6]/tbody/tr/td/table/tbody')

        return body.text

    def getSNPs(self):
        """Return a DataFrame of all the SNPs found on SNP GeneView."""
        text = self.geneViewBodyText()
        SNPsInSequence = {}

        snp = {
            'chromosome_position': None,
            'mrna_position': None,
            'rsid': None,
            'amino_acid_position': None,
            'missense': [],
            'nonsense': [],
            'synonymous': [],
            'shift': [],
            'reference': None,
            'A': None,
            'C': None,
            'T': None,
            'G': None,
            '-': None
        }

        typeset = {'missense',
                   'nonsense',
                   'synonymous',
                   'shift',
                   'reference'}
        for line in text.splitlines():

            lineList = line.split()
            # check if in the first line of a new SNP entry
            if lineList[0].isnumeric():
                # add the previous SNP to SNPsinSequence
                if snp['rsid']:
                    SNPsInSequence[rsid] = snp

                # create a new SNP
                snp = {
                    'chromosome_position': None,
                    'mrna_position': None,
                    'rsid': None,
                    'amino_acid_position': None,
                    'missense': [],
                    'nonsense': [],
                    'synonymous': [],
                    'shift': [],
                    'reference': None,
                    'A': None,
                    'C': None,
                    'T': None,
                    'G': None,
                    '-': None
                }
                chromosomePosition = int(lineList[0])
                mrnaPosition = int(lineList[1])
                rsid = lineList[2]
                snp['chromosome_position'] = chromosomePosition
                snp[
                    'mrna_position'] = mrnaPosition
                snp['rsid'] = rsid

            for elem in lineList:
                index = lineList.index(elem)

                if elem in 'ACTG-':
                    SNPType = lineList[index - 1]

                    if SNPType not in typeset:  # handle missing data
                        snp['rsid'] = None
                        break

                    if SNPType == 'reference':
                        snp['amino_acid_position'] = lineList[-1]
                        snp[elem] = SNPType
                        snp[SNPType] = elem
                        break

                    snp[elem] = SNPType
                    snp[SNPType].append(elem)

        mvQuery = []
        for rsid in SNPsInSequence.keys():
            mvQuery.append(rsid)

        # convert chromosome positions to build hg37
        df1 = rsidsToHG37Positions(mvQuery)

        df2 = pd.DataFrame.from_dict(SNPsInSequence, orient='index')
        df3 = df1.join(df2).drop('chromosome_position', axis=1).rename(
            columns={0: 'location'})
        return df3

    def getDataFrame(self):
        """Return a DataFrame with only genotyped SNPs."""
        df = self.getSNPs()
        alleles = {'rsid': [], 'genotype': []}
        g = GenomeReader(self.genomeFile)

        for rsid, genotype, chromosome, location in g.readGenome():
            if rsid in list(df['rsid']):

                # don't include ungenotyped SNPS in 23andMe data
                if '-' not in genotype:
                    alleles['rsid'].append(rsid)
                    alleles['genotype'].append(genotype)

            # handle old merges, don't include mtDNA
            elif chromosome.isnumeric() and int(chromosome) == self.chromosome:
                if int(location) in list(df['location']):

                    # don't include ungenotyped SNPS in 23andMe data
                    if '-' not in genotype:
                        alleles['rsid'].append(rsid)
                        alleles['genotype'].append(genotype)

        df2 = pd.DataFrame(data=alleles)
        df3 = df.join(df2.set_index('rsid')).dropna(subset=['genotype'])
        return df3

    def getMutants(self):
        """Return a DataFrame of all non-wildtype SNPs found in genomeFile."""
        df = self.data
        mutants = []
        for row in df.iterrows():

            # handle indels
            if row[1]['shift']:
                if row[1]['shift'][0] == '-':  # if deletion
                    if 'I' in row[1]['genotype']:
                        mutants.append(row)

                else:  # if insertion
                    if 'D' in row[1]['genotype']:
                        mutants.append(row)


            elif row[1]['reference'] * 2 != row[1]['genotype']:
                mutants.append(row)

        if mutants:
            return pd.DataFrame(data=mutants)
        else:
            return None

    def getReferenceSeqRecord(self):
        """Return a reference Seq object for the provided accession number."""

        handle = Entrez.efetch(db='nuccore', id=self.accessionNumber,
                               rettype='gb',
                               retmode='text')

        seqrec = SeqIO.read(handle, "gb", alphabet=IUPAC.ambiguous_dna)
        handle.close()

        return seqrec

    def getSeq(self):
        """Return a Seq object representing a personalized sequence."""
        refSeqRecord = self.referenceSeqRecord
        refseq = refSeqRecord.seq
        mutants = self.getMutants()

        # if genotypes are all wild-type, return the reference Seq object
        if mutants is None or mutants.empty:
            return refseq

        # create a MutableSeq object
        myseq = refseq.tomutable()
        offset = -1
        for row in mutants.iterrows():
            base = ''
            index = row[1][1]['mrna_position'] + offset
            genotype = row[1][1]['genotype']

            # check if homozygous:
            if genotype[0] == genotype[1]:
                base = genotype[0]

            # handle indels
            shift = row[1][1]['shift']
            if shift:  # if frameshift:

                if shift[0] == '-':  # if deletion
                    if base:  # if homozygous
                        base = row[1][1]['shift']

                        # make sure that the sequence to be deleted is present
                        assert myseq[index:index + len(base)] == base

                        # delete the sequence
                        for i in range(len(base)):
                            myseq.remove(index)
                            offset -= 1

                    else:
                        if self.indels:
                            base = row[1][1]['shift']

                            # make sure the sequence to be deleted is present
                            assert myseq[index:index + len(base)] == base

                            # delete the sequence
                            for i in range(len(base)):
                                myseq.remove(index)
                                offset -= 1

                else:  # if insertion
                    if base:  # if homozygous:
                        base = shift[0]
                        myseq.insert(index, base)
                        offset += len(base)

                    else:  # if heterozygous:
                        if self.indels:
                            base = shift[0]
                            myseq.insert(index,
                                         base)
                            offset += len(base)

            # handle snvs
            else:
                base = SequenceAligment.AMBIGUOUS_DNA_LETTERS[genotype]
                myseq[index] = base

        myseq.toseq()
        return myseq

    def getMultipleSeqAlignment(self):
        """Return a MultipleSeqAlignment object."""
        refseq = self.referenceSeqRecord
        myseq = self.getSeq()
        align = MultipleSeqAlignment(
            [refseq, SeqRecord(myseq, id='My Sequence')])
        return align

    def getPairwiseAlignment(self):
        """Return the formatted pairwise alignment."""
        from Bio import pairwise2
        refseq = self.referenceSeqRecord.seq
        myseq = self.getSeq()

        # globalms alignment using match/mismatch scores of 5/-4 and gap
        # penalties (open/extend) of 2/0.5. Code adapted from BioPython Tutorial
        # and Cookbook, Ch. 6.5.1
        alignments = pairwise2.align.globalms(refseq, myseq, 5, -4, -2, -0.5)
        return pairwise2.format_alignment(*alignments[0], full_sequences=True)


class CommandLine:
    """Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. it implements
    a standard command line argument parser with various argument options,
    a standard usage and help.

    Attributes: all arguments received from the commandline using
    .add_argument will be available within the .args attribute of object
    instantiated from CommandLine. For example, if myCommandLine is an object
    of the class, and requiredbool was set as an option using add_argument,
    then myCommandLine.args.requiredbool will name that option.

    Adapted from CommandLine by David Bernick, BME160, UCSC.
    March 8 2020
    """

    def __init__(self, inOpts=None):
        """Implement a parser to interpret the command line argv string."""

        import argparse
        self.parser = argparse.ArgumentParser(
            description='This program creates a personalized alignment for a '
                        'sequence of interest using 23andMe genome data, '
                        'comparing to the human reference genome.',
            epilog='If not specified on the command line, the gene ID, '
                   'accession number, and genome filename can all be supplied'
                   'by user input. Do not redirect your genome file to stdin.'
                   'You can specify an output filename if desired, otherwise'
                   'the alignment will be printed to stdout.'
                   '\n'
                   'Data is collected from NCBI\'s dbSNP GeneView, '
                   'dbNucleotide, and dbGene using myVariant, Bio.Entrez, '
                   'Selenium and Google\'s Chromedriver.',
            add_help=True,  # default is True
            prefix_chars='-',

        )
        self.parser.add_argument('infile', nargs='?', action='store',
                                 default='')
        self.parser.add_argument('outfile', nargs='?',
                                 type=argparse.FileType('w'),
                                 default=sys.stdout)
        self.parser.add_argument('-gI', '--geneID', action='store', nargs='?',
                                 help='the gene ID (e.g. \'79068\')',
                                 default=None)
        self.parser.add_argument('-acc', '--accession', action='store',
                                 nargs='?', default=None,
                                 help='the accession number (e.g. '
                                      '\'NM_001080432.2\')')
        self.parser.add_argument('-ind', '--indels', action='store',
                                 nargs='?', const=True, default=False,
                                 help='whether to include indels. '
                                      'default=False')
        self.parser.add_argument('-f', '--format', action='store',
                                 nargs='?', default='clustal',
                                 help='the multiple sequence alignment file '
                                      'format (e.g. \'phylip\'). '
                                      'default=\'clustal\'')
        self.parser.add_argument('-p', '--pairwise', action='store',
                                 nargs='?', const=True, default=False,
                                 help='whether to use a pairwise alignment. If '
                                      'false, use a MultipleSeqAlignment. '
                                      'Default=False')
        self.parser.add_argument('-v', '--version', action='version',
                                 version='%(prog)s 0.1')

        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


def main(inCL=None):
    """Get a personalized alignment."""
    import os

    if inCL is None:
        myCommandLine = CommandLine()
    else:
        myCommandLine = CommandLine(inCL)

    genomefile = myCommandLine.args.infile
    outfile = myCommandLine.args.outfile
    geneID = myCommandLine.args.geneID
    accession = myCommandLine.args.accession
    indels = myCommandLine.args.indels
    msaFormat = myCommandLine.args.format
    pairwise = myCommandLine.args.pairwise

    mrnaNames = []
    if not geneID:
        if not accession:
            while True:
                geneID = input('Enter a gene ID.').strip()
                mrnaNames = getAccessionsFromGene(geneID)
                if mrnaNames:
                    break
                else:
                    print(f'No variants found in GeneView for gene {geneID}.')

        else:
            while True:

                # if accession is invalid
                if accession[0:2] not in {'NM', 'NR', 'XM', 'XR'}:
                    accession = input(
                        'Invalid accession number. Enter an accession number '
                    ).strip()

                try:
                    geneID = getGeneFromAccession(accession)

                # if the accession is invalid
                except:
                    accession = input(
                        'Invalid accession number. Enter an accession number.'
                    ).strip()
                    continue

                mrnaNames = getAccessionsFromGene(geneID)
                if mrnaNames:
                    # a valid geneID and accession have been found
                    break

                # if no variants are displayed in GeneView
                else:
                    print(f'No variants found in GeneView for gene {geneID}.')
                    accession = input('Enter an accession number.').strip()

        # at this point a valid geneID has been assigned.
    if not mrnaNames:
        mrnaNames = getAccessionsFromGene(geneID)

    if accession not in mrnaNames:
        accession = None

    if not accession:

        # display available accession numbers
        print(f'Available mRNAs: {repr(mrnaNames)}')
        while True:
            accession = input('Enter an accession number.').strip()
            if accession not in mrnaNames:
                print(
                    f'Accession number {accession} not '
                    f'found in SNP GeneView Report.\n'
                )

            else:
                try:  # test if accession is found in Entrez Nucleotide
                    handle = Entrez.efetch(db='nuccore',
                                           id=accession,
                                           rettype='gb',
                                           retmode='text')
                    break

                except:
                    print(
                        f'Accession number {accession} not '
                        f'found in Entrez Nucleotide database.\n'
                    )
                    continue

    # if no genome_file was provided, ask user for input file name:
    if not genomefile:
        while True:
            genomefile = input('Genome file name:')
            if os.path.isfile(genomefile):  # check that input is a valid file
                break
            else:
                print('File not found.')

    sa = SequenceAligment(geneID, accession, genomeFile=genomefile,
                          indels=indels)

    if not pairwise:
        AlignIO.write([sa.getMultipleSeqAlignment()], handle=outfile,
                      format=msaFormat)
        outfile.close()
    else:
        if outfile:
            outfile.write(sa.getPairwiseAlignment())
        else:
            print(sa.getPairwiseAlignment())


if __name__ == "__main__":
    main()
