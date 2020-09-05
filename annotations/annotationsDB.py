import requests
import re
import json
from bs4 import BeautifulSoup

def isDB(cp):

    """
    given a compound name
    see if it doesnt exist in DB, return None
    otherwise, return the response
    """

    searchurl = f"https://www.drugbank.ca/unearth/q?utf8=%E2%9C%93&query={cp}&searcher=drugs"
    response = requests.get(searchurl)
    response.raise_for_status()    
    #if cp exists in db, redirects to another url
    f =  response.url.startswith("https://www.drugbank.ca/drugs/")
    if f:

        return response

    return


def searchDBClasses(cp):

    """
    given a compound name
    if it doesnt exist in DB, return None
    otherwise, return its identifier and classes
    """
    response = isDB(cp)
    if response is not None:

        DB  = response.url.split("/")[-1]
        regexp = r"<a [^>]*categories[^>]*>[^<]*"
        classifications = re.findall(regexp,response.text.replace("\n", ""))

        ret = [0 for _ in range(len(classifications))]
        ret[0] = DB

        for i,s in enumerate(classifications[1:]): #1st element is always blank
            ret[i+1] = re.sub("<[^>]*>","",s)


        return ret
    
    return

def _add_edges(edgeList, data):

    compound, description = BeautifulSoup(data[0]), data[1]
    compound_id = compound.a["href"].split("/")[-1]
    compound_name = compound.a.text
    edgeList.append((compound_id,compound_name,description))

def _getDBInteractionResponse(drugBankID, start, length):

    url = (f"https://www.drugbank.ca/drugs/{drugBankID}/drug_interactions.json?draw=6"
   "&columns%5B0%5D%5Bdata%5D=0&columns%5B0%5D%5Bname%5D=&columns%5B0%5D%5Bsearchable%5D=true"
   "&columns%5B0%5D%5Borderable%5D=true&columns%5B0%5D%5Bsearch%5D%5Bvalue%5D=&columns"
   "%5B0%5D%5Bsearch%5D%5Bregex%5D=false&columns%5B1%5D%5Bdata%5D=1&columns%5B1%5D%5Bname%5D=&columns"
   "%5B1%5D%5Bsearchable%5D=true&columns%5B1%5D%5Borderable%5D=true&columns%5B1%5D%5Bsearch%5D%5Bvalue%5D"
   f"=&columns%5B1%5D%5Bsearch%5D%5Bregex%5D=false&start={start}&length={length}&search%5Bvalue%5D=&search"
   "%5Bregex%5D=false&_=1592990200279")
    r =requests.get(url)
    r.raise_for_status()
    return json.loads(r.text)

def searchDBInteractions(cp):

    response = isDB(cp)

    if response is not None:

        edgeList = []
        drugBankID =  response.url.split("/")[-1]
        start = 0
        length = 100

        j = _getDBInteractionResponse(drugBankID, start, length)
        tot_records = j["recordsTotal"]
        
        for data in j["data"]:

            _add_edges(edgeList,data)

        for i in range(1, tot_records//1+1):

            start = 100*i
            j = _getDBInteractionResponse(drugBankID, start, length)
            
            for data in j["data"]:

                _add_edges(edgeList, data)

        return edgeList

    return
