def getEdgeProperty(G,eprop):
    """Given an object of class Graph, it returns an ordered array of edge property
    values.
    
    :param G: a graph that contains at least one edge.
    :param eprop: a property of the edges in the graph.
    :type G: an Graph() object from the networkx module
    :type eprop: string
    :returns: a list of ordered edge property values."""
    argcheck = [eprop in G[e[0]][e[1]].keys() for e in G.edges()]
    if not all(argcheck):
        raise NameError("Selected edge property is not assigned to every edge in the graph")    
    p = [G[e[0]][e[1]][eprop] for e in G.edges()]
    return p
    
def updateLabel(metabolite):
    """Given an object of class Metabolite, it updates the attribute *label*. If
    updated for the first time, it adds the suffix "\*2"; otherwise, it adds +1 to
    the numeric part of the suffix.
    
    :param metabolite: a metabolite.
    :type metabolite: a Metabolite() object""" 
    label = metabolite.label
    if "*" in label:
        ind = int(label[label.index("*")+1:]) + 1
        label = label[0:label.index("*")+1] + str(ind)
        metabolite.label = label
    else:
        metabolite.label = label + "*2"
    
def drawHeatmap(PG,hmType,filename,pzero,convertGeneNamesFromType="",
                    convertGeneNamesToType="",convertGeneNamesSpecies=""):
    """Given an object of class Graph where nodes are of class Metabolite, it constructs
    and display a matrix where rows are the metabolites and the columns are the overall set of
    genes associated with the metabolites. Each entry indicates if the gene is associated
    with the metabolites. The entry can be binary if hmType = "Binary", otherwise
    it gets the -log10 value for the attribute *p* of the gene. In the latter case,
    the minimum value the logarithm can get is set by pzero. The matrix is displayed
    as a heatmap, where rows and columns are hierarchically clustered. Gene labels can be 
    interconverted using different annotation sources.     
    
    :param PG: a graph that contains Metabolite objects as nodes.
    :param hmType: either "Binary" or "Numeric", according to whether the heatmap entries are binary or -log10(p) of each gene.
    :param filename: the name of the output figure.
    :param pzero: the minimum value the logarithm can accept.
    :param convertGeneNamesFromType: the annotation source for the genes in each Metabolite object in the graph. It should match the name in mygene.info. *Default: ""*.
    :param convertGeneNamesToType: the annotation source preferred to be shown in the heatmap for the genes in each Metabolite object in the graph. It should match the name in mygene.info. *Default: ""*.
    :param convertGeneNamesSpecies: the species the gene names refer to. It should match the name in mygene.info. *Default: ""*.
    :type PG: an Graph() object from the networkx module
    :type hmType: string
    :type filename: string
    :type pzero: float
    :type convertGeneNamesFromType: string
    :type convertGeneNamesToType: string
    :type convertGeneNamesSpecies: string"""
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('PDF')
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import scipy.spatial.distance as distance
    import scipy.cluster.hierarchy as sch
    
    # Helper for cleaning up axes by removing ticks, tick labels, frame, etc:
    def clean_axis(ax):
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
            
    # Get gene and metabolites from PG (skip genes without data):
    listOfGeneLists = [m.geneList for m in PG.nodes()]
    PGgeneList = np.unique([g for gList in listOfGeneLists for g in gList if not np.isnan(g.p)]) 
    PGmetList = [m for m in PG.nodes()]
    
    # Construct heatmap data matrix:
    hmMatBinary = np.zeros([len(PGmetList),len(PGgeneList)])
    hmMatValues = np.zeros([len(PGmetList),len(PGgeneList)])*np.nan
    mInd = 0
    for met in PGmetList:
        # Get gene index for genes in met (skip genes without data):
        gInd = [int(np.where(PGgeneList==g)[0]) for g in met.geneList if not np.isnan(g.p)] 
        hmMatBinary[mInd,gInd] = 1
        hmMatValues[mInd,gInd] = \
           [-np.sign(g.FC)*np.log10(g.p+pzero) for g in PGgeneList[gInd]] 
        mInd = mInd + 1
    
    # Convert to dataframe:
    hmMatBinary = pd.DataFrame(hmMatBinary)   
    hmMatBinary.index = [m.label for m in PGmetList]  
    hmMatValues = pd.DataFrame(hmMatValues)   
    hmMatValues.index = [m.label for m in PGmetList] 
    columnLabs =  [g.name for g in PGgeneList]  
    if not (len(convertGeneNamesFromType) == 0 or len(convertGeneNamesToType) == 0 or len(convertGeneNamesSpecies) == 0):
        columnLabs = convertGeneNames(columnLabs,convertGeneNamesFromType,convertGeneNamesToType,convertGeneNamesSpecies)
    hmMatBinary.columns = columnLabs
    hmMatValues.columns = columnLabs
    
    # Initialize plot:
    plt.rcParams.update({'font.size': 8})
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
    fig = plt.figure(figsize=(8,8))
    hm = gridspec.GridSpec(2,2,wspace=0.0,hspace=0.0,width_ratios=[0.25,1],height_ratios=[0.1,1])
    
    # Alt A. Order the columns by gene stats:
    #geneLogP = np.zeros([1,sum([1 for g in PGgeneList])])
    #geneLogP[0,:] = np.asarray([-np.sign(g.FC)*np.log10(g.p+0.0000001) for g in PGgeneList])
    #colorder = np.argsort(geneLogP)[0]
    # Alt B. Or cluster the columns:
    col_pairwise_dists = distance.squareform(distance.pdist(hmMatBinary.T))
    col_clusters = sch.linkage(col_pairwise_dists,method='complete')
    col_denD = sch.dendrogram(col_clusters,no_plot=True)
    colorder = col_denD['leaves']
    
    # Cluster the metabolites and make dendrogram (rows):
    pairwise_dists = distance.squareform(distance.pdist(hmMatBinary))
    row_clusters = sch.linkage(pairwise_dists,method='complete')
    row_denAX = fig.add_subplot(hm[1,0])
    sch.set_link_color_palette(['black'])
    row_denD = sch.dendrogram(row_clusters,color_threshold=np.inf,orientation='right',
                              color_list=['black'])
    clean_axis(row_denAX)
    
    # The heatmap:
    heatmapAX = fig.add_subplot(hm[1,1])
    if hmType=='binary':
        plotmat = hmMatBinary
        c_map = plt.cm.Greys
        vmin = 0
        vmax = 1
    else:
        plotmat = hmMatValues
        c_map = plt.cm.RdBu_r
        vmax=plotmat.abs().max().max()
        vmin=-vmax
    axi = heatmapAX.imshow(plotmat.ix[row_denD['leaves'],colorder],interpolation='nearest',
                           aspect='auto', origin='lower',cmap=c_map, 
                           vmin=vmin,vmax=vmax)
    clean_axis(heatmapAX)
    
    # Add labels:
    heatmapAX.set_yticks(np.arange(hmMatValues.shape[0]))
    heatmapAX.yaxis.set_ticks_position('right')
    heatmapAX.set_yticklabels(hmMatValues.index[row_denD['leaves']])
    
    heatmapAX.set_xticks(np.arange(hmMatValues.shape[1]))
    xlabelsL = heatmapAX.set_xticklabels(hmMatValues.columns[col_denD['leaves']])
    for label in xlabelsL:
        label.set_rotation(90)
    for l in heatmapAX.get_xticklines() + heatmapAX.get_yticklines(): 
        l.set_markersize(0)
        
    # Colorkey:
    if hmType=='values':
        scale_cbGSSS = gridspec.GridSpecFromSubplotSpec(2,1,subplot_spec=hm[0,1],wspace=0.0,hspace=0.0,height_ratios=[0.3,0.7])
        scale_cbAX = fig.add_subplot(scale_cbGSSS[0,0]) # colorbar for scale in upper left corner
        cb = fig.colorbar(axi,scale_cbAX,orientation='horizontal') # note that we tell colorbar to use the scale_cbAX axis
        cb.set_label('log10 p-values')
        cb.ax.yaxis.set_ticks_position('left') # move ticks to left side of colorbar to avoid problems with tight_layout
        cb.ax.yaxis.set_label_position('left') # move label to left side of colorbar to avoid problems with tight_layout
        cb.outline.set_linewidth(0)
    
    # Draw the figure:
    fig.tight_layout()
    if len(filename)>0: 
        plt.savefig(filename, bbox_inches='tight')
        plt.close(fig)
    else:
        plt.show()

def convertGeneNames(geneList, fromtype, totype,species,keepFirst=False):
    """ It uses information for gene and species annotations from mygene.info to 
    convert a list of gene names in a different annotation. I{e.g. fromtype="'ensembl.gene', 
    totype = "'symbol'", species="human" converts the gene IDs from ENSEMBL ID 
    to Official Gene Symbol for human genes.} The returned converted list of gene
    names is ordered according to the input list. NOTE: if no hit is found, the 
    original gene name is returned; if multiple hits are found, only the first hit
    is returned. 
    
    :param geneList: a list of gene names.
    :param fromtype: the annotation source for the genes in each Metabolite object in the graph. It should match the name in mygene.info.
    :param totype: the annotation source preferred to be shown in the heatmap for the genes in each Metabolite object in the graph. It should match the name in mygene.info.
    :param species: the species the gene names refer to. It should match the name in mygene.info.
    :type geneList: list
    :type fromtype: string
    :type totype: string
    :type species: string
    :returns: an ordered list of gene names that match the input gene list."""
    
    import mygene
    mg  = mygene.MyGeneInfo()
    raw = mg.querymany(geneList, scopes=fromtype, fields=totype, species=species,returnall=True)
    out = raw['out'] #The actual dictionary with the mappings
    #If a query returns no hit, the mapping is the query itself
    for diction in out:
        if "notfound" in diction.keys() and diction["notfound"]:
            diction[totype]=diction["query"]   
            print("Gene id unchanged for: "+diction["query"]) 
    #If a query returns multiple hits, either...
    if keepFirst:
    # ...keep the first one
        alldupl = [d[0] for d in raw["dup"]]
        for duplo in alldupl:
            count = 0
            for dic in out:
                if dic["query"] == duplo:
                    count = count + 1
                    dic["isdupl"] = count
        outf = [dic for dic in out if "isdupl" not in dic.keys() or dic["isdupl"] <= 1]
    else:
    # ...concatenate all possible mappings in a string separated by "*"
        alldupl = [d[0] for d in raw["dup"]]
        for duplo in alldupl:
            #First build the name by concatening all instances of the duplicate mapping
            name = ""
            for dic in out:
                if dic["query"] == duplo:
                    name = name + dic[totype] +"*"
            name = name[:-1]
            #Then assign to all duplicates this name
            count = 0
            for dic in out:
                if dic["query"] == duplo:
                    count = count + 1
                    dic[totype] = name
                    dic["isdupl"] = count   
            #Finally remove all duplicates but one
        outf = [dic for dic in out if "isdupl" not in dic.keys() or dic["isdupl"] <= 1]
    #Return a list that index the mapping sorted as the geneList
    id_out = [gene[totype] for gene in outf]   
    return id_out