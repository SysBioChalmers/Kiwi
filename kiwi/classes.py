class Metabolome(object):
    """An object of class Metabolome contains a list of Metabolite objects."""
    def __init__(self):
        self.metaboliteList = []
        
    def __str__(self):
        s = ""
        for met in self.metaboliteList:
            s = s + met.name + "\n"
        return s
        
    def size(self):
        """It returns the number of Metabolite objects."""
        return len(self.metaboliteList)
    
    def addMetabolite(self,metabolite):
        """It adds a Metabolite object to the Metabolome.
        
        :param metabolite: the metabolite to be added.
        :type metabolite: a Metabolite object"""
        self.metaboliteList.append(metabolite)
        
    def removeMetabolite(self,metabolite):
        """It removes a Metabolite object to the Metabolome.
        
        :param metabolite: the metabolite to be removed.
        :type metabolite: a Metabolite object"""
        self.metaboliteList.remove(metabolite)
    
    def getMetabolite(self,name):
        """Given a metabolite name, it returns the corresponding Metabolite object.
        If no object is found with that name, it returns an empty list.
        
        :param name: a metabolite name.
        :type name: string
        :returns: a Metabolite object that matches the name or an empty list if no match was found."""
        try:
            met = [m for m in self.metaboliteList if m.name==name][0]
        except:
            met = []
        return met
                
    def removeNotSignificantMetabolites(self,pcutoff,nonPiano):
        """It removes all Metabolite objects whose *pNonDirectional* attribute is
        above pcutoff.
        
        :param pcutoff: a minimum threshold for the *pNonDirectional* or *pValue*.
        :param nonPiano: a logical to determine whether to use *pNonDirectional* or *pValue*.
        :type pcutoff: float
        :type nonPiano: logical"""
        if nonPiano:
            argcheck = [hasattr(met,'pValue') for met in self.metaboliteList]
            if not all(argcheck):
                raise NameError("pValue attribute is not assigned to Metabolite in the Metabolome")   
            mlOld = self.metaboliteList
            self.metaboliteList = [met for met in mlOld if met.pValue < pcutoff]    
        else:
            argcheck = [hasattr(met,'pNonDirectional') for met in self.metaboliteList]
            if not all(argcheck):
                raise NameError("pNonDirectional attribute is not assigned to Metabolite in the Metabolome")   
            mlOld = self.metaboliteList
            self.metaboliteList = [met for met in mlOld if met.pNonDirectional < pcutoff]
        
    def removeMetabolitesNotInMetNet(self,MN):
        """It removes all Metabolite objects that are not nodes in the input
        metabolic network Graph.
        
        :param MN: a graph that contains Metabolite objects as nodes.
        :type MN: a networkx Graph object."""
        mlOld = self.metaboliteList
        self.metaboliteList = [met for met in mlOld if met.name in MN.node]
         
    def removeHighDegreeMetabolites(self,MN,degreecutoff):
        """It removes all Metabolite objects if the degree as nodes in the input
        metabolic network Graph is higher than degreecutoff.
        
        :param MN: a graph that contains Metabolite objects as nodes.
        :param degreecutoff: a minimum threshold for the degree of a metabolite in the metabolic graph.
        :type degreecutoff: float
        :type MN: a networkx Graph object."""
        import networkx as nx
        mlOld = self.metaboliteList
        self.metaboliteList = [met for met in mlOld if nx.degree(MN,met.name) < degreecutoff]  
            
class Metabolite(object):
    """An object of class Metabolite represents a metabolite to which gene-set
    statistics can be assigned and genes be associated."""
    
    def __init__(self,name):
        self.name = name
        self.label = name
            
    def __str__(self):
        return self.name
        
    def addGeneSetStats(self,stats,header,adj):
        """It adds gene-set statistics as attributes. The statistics values are 
        contained in stats and the mapping to the correct attribute is enforced 
        by the header that matches stats. It supports the use of adjusted statistics.
        
        :param stats: a list of gene-set statistics.
        :param header: a list of names that correspond to the different statistics in stats.
        :param adj: either empty or "adj" if adjusted statistics should be used.
        :type stats: a numerical array
        :type header: a string array
        :type adj: string  
        """
        import numpy as np
        try:
            self.pNonDirectional = float(stats[np.where(header=='p '+adj+'(non-dir.)')[0]])
        except:
            self.pNonDirectional = np.nan
        try:
            self.pMixDirUp = float(stats[np.where(header=='p '+adj+'(mix.dir.up)')[0]])
        except:
            self.pMixDirUp = np.nan
        try:
            self.pMixDirDn = float(stats[np.where(header=='p '+adj+'(mix.dir.dn)')[0]])
        except:
            self.pMixDirDn = np.nan
        try:
            self.pDistDirUp = float(stats[np.where(header=='p '+adj+'(dist.dir.up)')[0]])
        except:
            self.pDistDirUp = np.nan
        try:
            self.pDistDirDn = float(stats[np.where(header=='p '+adj+'(dist.dir.dn)')[0]])
        except:
            self.pDistDirDn = np.nan
        #If non piaNo file, the sixth column contains the gene-set p-values
        if len(adj)>0:
            try:
                self.pValue = float(stats[np.where(header=='p-adj')[0]])
            except:
                self.pValue = np.nan
        else:
            try:
                self.pValue = float(stats[np.where(header=='p-value')[0]])
            except:
                self.pValue = np.nan
            
    def addGene(self,gene):
        """It adds a Gene object to the gene list associated with the metabolite.
        
        :param gene: a gene associated with the metabolite.
        :type gene: a Gene object."""
        if not isinstance(gene,Gene):
            raise NameError("Object must be of class Gene")
        else:
            if hasattr(self, 'geneList'):
                self.geneList.append(gene)
            else:
                self.geneList = [gene]

class Genome(object):
    """An object of class Genome is a list of Gene objects."""
    def __init__(self):
        self.geneList = []
        
    def __str__(self):
        s = ""
        for gene in self.geneList:
            s = s + gene.name + "\n"
        return s
        
    def size(self):
        """It returns the number of Gene objects."""
        return len(self.geneList)
    
    def addGene(self,gene):
        """It adds a Gene object to the list of gene in the Genome.
        
        :param gene: a gene.
        :type gene: a Gene object."""
        self.geneList.append(gene)   

    def getGene(self,name):
        """Given a gene name, it returns the corresponding Gene object.
        If no object is found with that name, it returns an empty list.
        
        :param name: a gene name.
        :type name: string
        :returns: a Gene object that matches the name or an empty list if no match was found."""
        try:
            gene = [g for g in self.geneList if g.name==name][0]
        except:
            gene = []
        return gene
                
class Gene(object):
    """An object of class Gene represents a gene to which statistics can be assigned."""
    def __init__(self,name):
        self.name = name
        
    def __str__(self):
        return self.name
        
    def addGeneStats(self,p,FC):
        """It adds *p*-value and fold-change statistics as gene attributes.
        
        :param p: the gene *p*-value.
        :param FC: the gene fold-change.
        :type p: float
        :type FC: float"""
        self.p = p
        self.FC = FC