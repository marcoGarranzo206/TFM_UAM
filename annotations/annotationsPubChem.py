import requests
class pubChemRet:

    """
    to search if compound_name is in pubChem:

        domain = "compound", namespace = "name" in __init__
        search: compound_name, operation = "cids"

    to search for a classification of a CID:

        domain = compound, namespace = CID, output = JSON in __init__
        search: CID, operation = "classification"
    
    self._content contains the results from these operations
    more info at https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest¶
    """

    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

    def __init__(self, domain, namespace,output = "JSON"):

        self.domain = domain
        self.namespace = namespace
        self.output = output

    def modDomain(self, newDomain):

        self.domain = newDomain

    def modNamespace(self, newNamespace):

        self.namespace = newNamespace

    def modOutput(self, newOutput):

        self.output = newOutput

    def search(self, ID, params = {}, operation = None):

        self.ID = ID
        if operation is None:

            url = pubChemRet.base_url + "/" + self.domain + "/"+ self.namespace + "/" + self.ID + "/"  + self.output

        else:

            url =  pubChemRet.base_url + "/" + self.domain + "/"+ self.namespace + "/" + self.ID + "/" + operation  + "/" + self.output

        response = requests.get(url, params=params)
        response.raise_for_status()
        self.__dict__.update(response.__dict__)

