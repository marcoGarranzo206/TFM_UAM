from xml.etree import ElementTree as ET
import requests
from Bio import Entrez
from xml.parsers import expat

def get_IDS(url):
    
    Esearch = requests.get(url) 
    root = ET.fromstring(Esearch.text)
    Idset = set()

    for child in root:
        
        if child.tag == "Count" and child.text == "0":
            
            print("No results found")
            return Idset

        if child.tag == "IdList":

            for ID in child:

                if ID.text[0] == "2":
                    
                    Idset.add(ID.text)

            break

    return Idset

def write_GSE_info(dictionary, attributes, Filter, filename, logfile = "log.txt" ,sep = "\t"):

    with open(filename, "a+") as f, open(logfile, "a+") as l:
    
        for element in Filter.keys():

            if Filter[element] != dictionary[element]:

                l.write( f"{dictionary['Accession']}:{str(element)} is {str(dictionary[element])}, not {Filter[element]}\n" )
                break

            else:

                l.write( f"{dictionary['Accession']}:{str(element)} not found\n" )
                 
        f.write( sep.join([str(dictionary[element]) for element in attributes ]) + "\n" ) 

def write_GSEs(gseList, filename):

    length = len(gseList)
    maxSearch = 300
    Entrez.email ="eigtw59tyjrt403@gmail.com"
    #with possibly hundreds or thousands of GSEs
    #we dont want to stop the whole proccess if one batch fails
    #if it does, log it for later

    keyErrorList = []
    notReadErrorList = []
    otherErrorsList = []
    for i in range(length//maxSearch + 1):

        idx = i*maxSearch
        to_search = gseList[idx:idx+maxSearch]
        postRequest = Entrez.epost("gds", id = ",".join(to_search))

        try:

            postResult = Entrez.read(postRequest)
            w = postResult["WebEnv"]
            q = postResult["QueryKey"]
            summaryRequest = Entrez.esummary(db = "gds", query_key= q, WebEnv=w)
            summaryResult = Entrez.read(summaryRequest)
            for result in summaryResult:

                write_GSE_info(result, ["Id", "Accession", "title", "summary", "GPL", "GSE", "taxon", "gdsType"],\
                          {"taxon": "Homo sapiens"}, filename = filename + ".txt", logfile=filename + "Log.txt")
        except KeyError:

            keyErrorList.extend(gseList[idx:idx+maxSearch])
        except (expat.ExpatError, ValueError):

            notReadErrorList.extend(gseList[idx:idx+maxSearch])
        except Exception as e:

            otherErrorsList.extend([ str(e) + ":" + gse for gse in gseList[idx:idx+maxSearch] ])
    if keyErrorList:

        with open(filename + "_key_error.txt", "w+") as f:

            for gse in keyErrorList:

                f.write(gse+"\n")

    if notReadErrorList:

        with open(filename + "_not_read_error.txt", "w+") as f:

            for gse in notReadErrorList:

                f.write(gse+"\n")

    if otherErrorsList:

        with open(filename + "_other_errors.txt", "w+") as f:

            for gse in otherErrorsList:

                f.write(gse+"\n")

def get_column_info(html_text, colname):

    '''
    Extract from GEO sample html page information relating to a column
    Example url: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2674900
    Example Columns: Source name, Characteristics
    
    Source name: simple column, one value "Pca_22Rv1"
    Characteristics: more complicated. Column itself is a table. Each line, separated by an <br>, is a key:value
    
    For simple columns, extract the value
    For complicated columns, extract key and value pair
    
    '''

    opening_td = r"<td[^>]*>"
    #opening column tag, we dont care for its characteristics [^>]* allows anything but a closing bracket
    closing_td = r'</td>' # simple closing tag   
    regexp_row = re.compile(opening_td+colname+closing_td + opening_td + r"([^<]*(<br>|<a[^>]*>[^<]*</a>)*)*" + closing_td) 
    #extract entire row
    stripped = html_text.replace("\n", "")
    stripped = re.sub("http[s]{0,1}:", "" ,stripped) # remove https: to avoid errors when splitting by :

    e,s = list(regexp_row.finditer(stripped))[0].span()
    regexp_col = re.compile(opening_td+colname+closing_td)
    s_noncol = stripped[e:s]
    s_noncol = s_noncol[regexp_col.match(stripped[e:s]).span()[1]:] #extract entire column content
    splitted = re.findall( r"<td[^>]*>(.*?)</td>", s_noncol)[0].split("<br>") #extract column values, separated by <br>

    return [i.split(":")  for i in splitted][:-1]            
        

