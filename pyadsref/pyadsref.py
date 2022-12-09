import pybtex
from pybtex.database import parse_file
from pybtex.database import parse_string

import requests

from bs4 import BeautifulSoup

import string

try:
    import pyperclip
except ImportError:
    msg = "\n Module 'pyperclip' not installed. Function flag \n "
    msg += "to copy output bibtex to clipboard will not function."
    print(msg)

__author__ = "Steven Gough-Kelly"
__copyright__ = "Copyright 2022, UCLan Galaxy Dynamics"
__credits__ = ["Steven Gough-Kelly"]
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Steven Gough-Kelly"
__email__ = "sgoughkelly@gmail.com"
__status__ = "Production"

def ads_lookup(doi):
    """
    Takes a paper doi and searches NASA ADS abstract service for the bibtex
    citation and uses BeautifulSoup to scrape the html.

    Parameters
    ----------
    doi : string
        The doi of the paper of interest.
    Returns
    -------
    bibtex : string of the bibtex citation
    """

    adswebadd = 'https://ui.adsabs.harvard.edu/abs/'
    expcite = '/exportcitation'

    doi = adswebadd + doi + expcite

    page = requests.get(url)
    soup = BeautifulSoup(page.content, "html.parser")
    elements = soup.find_all("textarea",
                             class_="export-textarea form-control")

    return elements[0].text

def rm_punc(text):
    """
    Removes the punctuation from a string except dashes for double-barrelled
    surnames

    Parameters
    ----------
    text : string
        The input text to remove punctuation from.
    Returns
    -------
    output string : text with punctuation removed.
    """

    punc = string.punctuation
    punc = punc.replace('-', '')

    return text.translate(str.maketrans('', '', punc))

def gen_ref(refbib):
    """
    Generates a pybtex BibliographyData object with one new entry created
    from a bibtex string.

    Parameters
    ----------
    refbib : string
        The input bibtex entry.
    Returns
    -------
    new entry : BibliographyData object with one entry.
    """

    tmpent = parse_string(refbib, 'bibtex')

    key = list(tmpent.entries.keys())[0]
    tmpauth = tmpent.entries[key].persons['author']
    author_list = tmpent.entries[key].persons['author']
    ln = author_list[0].last_names[0].split('{')[1].split('}')[0].lower()
    ln = rm_punc(ln)

    if len(author_list)>0:
        ln += '+'
    else:
        ln += '_'

    ln += tmpent.entries[key].fields['year']

    lc = refbib.find('{')+1
    rc = refbib.find(',')

    newrefbib = refbib[:lc]+ln+refbib[rc:]

    return parse_string(newrefbib, 'bibtex')

def sort_db(db):
    """
    Sorts a pybtex BibliographyData object by creating a new object and
    building it from the sorted keys of input object

    Parameters
    ----------
    db : BibliographyData object
        The input Bibliography database to sort.
    Returns
    -------
    sorted db : BibliographyData object of same size as db.
    """

    keys = list(db.entries.keys())
    keys.sort()

    sortout = parse_string(allrefs.entries[keys[0]].to_string('bibtex'),
                          'bibtex')
    for key in keys[1:]:
        sortout.add_entry(key,db.entries[key])

    return sortout

def add_new_entry(db,doi):
    """
    Appends a new entry to a  pybtex BibliographyData object from a paper doi.

    Parameters
    ----------
    db : BibliographyData object
        The input Bibliography database to sort.
    doi : string
        The doi of the paper of interest.
    Returns
    -------
    output db : sorted BibliographyData object of db with new entry
    """

    new_entry = gen_ref(ads_lookup(doi))
    # check db for key and figure out repeats.
    nkey = list(newent.entries.keys())[0]
    db.add_entry(nkey,newent.entries[nkey])

    return sort_db(db)

def save_bib(db,filename='./db.bib',clipboard=False):
    """
    Saves a BibliographyData object db to a file and optional clipboard.

    Parameters
    ----------
    db : BibliographyData object
        The input Bibliography database to save.
    filename : string
        The location and filename to save the bibtex file to.
    clipboard: bool
        Flag to copy output bibtex to clipboard for easy pasting into overleaf
        without having upload new file.
    Returns
    -------
    None
    """

    with open(filename, 'w') as f:
        f.write(db.to_string('bibtex'))

    if clipboard:
        pyperclip.copy(db.to_string('bibtex'))

    return None
