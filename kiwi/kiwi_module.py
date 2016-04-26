# -*- coding: utf-8 -*-
def plot(gsn, gss, gsc='', gls='', pc=0.01, adj=True, gsz=10**6,
         deg=10**6, spl=2, sae=False, lmp=False, nsz=300, eds=5, lbs=9, hmt='binary', 
         cft='', ctt='', csp='', nwf='', hmf='', gml=''):
    """Kiwi: a tool to combine gene-set analyses with biological networks.
    
    Produces a network plot and a heatmap given a gene-set interaction network, 
    gene-set statistics, gene-set collection and gene-level statistics.
    
    Parameters
    ----------
    
    INPUT
    
    gsn : The geneset-geneset interaction file path. It is a 2 column table where 
          each row contains two interacting genesets, each one per column, with 
          no header. [required]
          
    gss : The geneset statistics result file path. It is, in its minimum 
          implementation, a 2 or more column table where each row must contain 
          as first column the genesets and in any column a statistics value for 
          its significance. The header must identify the first column as Name 
          and the geneset statistics column as p-value. If -ad is set to True, 
          adjusted p-value are expected to be found in a column identified by 
          the header as p-adj. Alternatively, a full geneset result file 
          generated using piano can be supplemented. If generated using the 
          function writeFilesForKiwi.R in piano, the file is parsed 
          automatically. [required]

    gsc : The geneset collection (geneset-gene) file path. It is a two column 
          table with genesets - genes associations in each row and no header.
          [default = '']

    gls : The gene-level statistics file path. It is a three column table with 
          gene - p-value - fold-change in each row. The genes must be as first 
          column, while p values and fold-changes must be identified in the 
          headers as p-value and FC. [default = '']
          
    NETWORK ANALYSIS PARAMETERS

    pc  : This flag controls the maximum threshold for the significance of a 
          geneset, after which it is discarded from the results. [default = 0.01]

    adj : This flag controls if p-values adjusted for multiple testing should be 
          used. [default = True]

    gsz : This flag control the maximum number of genes for a geneset after which 
          it is discarded from the results. It can be used to exclude high gene 
          count genesets, whose interpretability would be difficult and possibly 
          messy to display in the heatmap. [default = 10**6]

    deg : This flag controls the maximum degree of a geneset in the interaction 
          network after which it is discarded from the results. It can be used 
          to exclude very connected genesets (like ATP in a metabolic network) 
          from the plots. [default = 10**6]

    spl : This flag controls the maximum shortest path length between two genesets 
          in the interaction network after which no edge would connect them. This 
          can be interpreted as the threshold after which two genesets are 
          considered unrelated. [default = 2]
          
    sae : Show all edges. This flag controls if all but the best edges (defined
          by having the shortest path length) for each node should be removed 
          before plotting. If set to True, all edges (passing the spl cutoff)
          will be drawn. If set to False, only the best edges for each node will
          be included. [default = False]
          
    lmp : Lump genesets that share an identical gene list. This flag controls if 
          only one geneset (randomly chosen) among all genesets that share an identical
          list of genes should be kept in the geneset-geneset interaction network. 
          If set to True, all but one geneset among those that share an identical
          gene lsit will be removed from the interaction network. If set to False, 
          all genesets are treated as independent entities in the network even if
          they share an identical gene list. [default = False]

    PLOTTING PARAMETERS

    nsz : This flag defines the node size of genesets with the highest p value 
          to be plotted (in the extreme scenario, this equals pCutoff).
          [default = 300]

    eds : This flag control the scaling of edge widths with one increment in the 
          shortest path length between two connected genesets. [default = 5]

    lbs : This flag sets the label sizes in the plots. [default = 9]

    hmt : This flag sets the color of the entries in the heatmap to black or 
          white when set to “binary” or to a blue-to-red colormap according to 
          the fold-change when set to “values”. [default = 'binary']

    cft : This flag reports the annotation source for the genes in the GSA run 
          (e.g. Ensembl). It should match the annotation sources listed in 
          mygene.info (e.g. ensembl.gene). [default = '']

    ctt : This flag reports the annotation desired to plot the gene names in the 
          heatmap (e.g. HUGO). It should match the annotation sources listed in 
          mygene.info (e.g. symbol).  [default = ''] :

    csp : This flag reports the species for the gene annotation (e.g. Homo 
          sapiens). It should match the species listed in mygene.info (e.g. 
          human). [default = '']

    EXPORT OPTIONS

    nwf : This flag defines the name of the file where the network plot is saved 
          (as PDF). If empty, it displays in the current device. [default = '']

    hmf : This flag defines the name of the file where the heatmap is saved (as 
          PDF). If empty, it displays in the current device. [default = '']

    gml : This flag defines the name of the file where the graph shown in the 
          network plot is saved (as GraphML).
    """
    # Import stuff: 
    import matplotlib
    matplotlib.use('PDF')
    import itertools
    import matplotlib.pyplot as plt
    import numpy as np
    import networkx as nx
    import os
    import classes as kiwiC
    import functions as kiwiF
    import warnings
    
    MNfile      = gsn
    GSfile      = gss
    GMfile      = gsc
    Gfile       = gls
    pcutoff     = float(pc)
    splcutoff   = int(spl)
    dcutoff     = int(deg)
    minNodeSize = float(nsz)
    eScaleFac   = float(eds)
    labSize     = float(lbs)
    nwPlotFile  = nwf
    hmType      = hmt
    hmPlotFile  = hmf
    graphMLFile = gml
    if adj:
        adj = 'adj '
    else:
        adj = ''
    maxGSsize   = float(gsz)
    cgnFromType = cft
    cgnToType   = ctt
    cgnSpecies  = csp
    lump        = lmp
    arbitrary_small = 1e-30
    
    # Check for bad input stuff:
    if not os.access(MNfile,os.R_OK): raise IOError("Specified gene-set network file cannot be accessed")
    if not os.access(GSfile,os.R_OK): raise IOError("Specified gene-set statistics file cannot be accessed")
    if len(Gfile)==0: warnings.warn('No gene-level statistics files is provided: no heatmap will be generated',RuntimeWarning)
    if len(GMfile)==0: warnings.warn('No gene-geneset association file is provided: no heatmap will be generated',RuntimeWarning)
    if len(GMfile)==0 and lump: warnings.warn('No gene-geneset association file is provided: no lumping of overlapping gene-sets will be performed\n nor heatmap will be generated',RuntimeWarning)
    if not os.access(GMfile,os.R_OK) and len(GMfile) > 0: raise IOError("Specified gene-geneset association file cannot be accessed")
    if not os.access(Gfile,os.R_OK) and len(Gfile) > 0: raise IOError("Specified gene-level statistics file cannot be accessed")
    if pcutoff >= 1: raise NameError("pc must be lower than 1")
    if splcutoff < 0: warnings.warn('spl is negative and it has been set to 0',RuntimeWarning)
    if minNodeSize <= 0: raise NameError("nsz must be a positive value")
    if eScaleFac <= 0: raise NameError("eds must be a positive value")
    if labSize <= 0: raise NameError("nsz must be a positive value")
    if minNodeSize <= 0: raise NameError("lbs must be a positive value")
    if not hmType in ["binary","values"]: raise NameError("hmt must be either 'binary' or 'values'")
    if maxGSsize <= 0: raise NameError("gsz must be a positive value")
    if len(GMfile)==0 and maxGSsize < 10**6: raise NameError("gsz has been supplied but no gene-geneset association file is provided")
    if lump: warnings.warn('Lump genesets that share an identical gene list is set to True. This will remove genesets from the interaction network',RuntimeWarning)

    # A general comment: note that in the code below it is assumed that the
    # gene-set interaction network is a metabolic network and that gene-sets
    # are metabolites, hence the nomenclature. Of course the code still works
    # for the general case!
    
    # Make metabolic network MN:
    MN = nx.read_edgelist(MNfile,nodetype=str,delimiter='\t')
    
    # Initialize metabolome M:
    M = kiwiC.Metabolome()
    
    # Read the metabolites from stats file and add to M:
    content = np.genfromtxt(GSfile,dtype=None,delimiter='\t')
    header = content[0]
    if header[0]!='Name': raise NameError("Gene-set statistic file has invalid header: first column should be named 'Name'.")
    if 'p (non-dir.)' in header and 'p (mix.dir.up)' in header and 'p (mix.dir.dn)' in header and 'p (dist.dir.up)' in header and 'p (dist.dir.dn)' in header:
        if len(adj)>0 and 'p adj (non-dir.)' not in header: raise NameError("Gene-set statistics file has invalid header: one column should be named 'p adj (non-dir).' or pAdjusted should be set as 'False'")
        if len(adj)>0 and 'p adj (mix.dir.up)' not in header: raise NameError("Gene-set statistics file has invalid header: one column should be named 'p adj (mix.dir.up).' or pAdjusted should be set as 'False'")
        if len(adj)>0 and 'p adj (mix.dir.dn)' not in header: raise NameError("Gene-set statistics file has invalid header: one column should be named 'p adj (mix.dir.dn).' or pAdjusted should be set as 'False'")
        if len(adj)>0 and 'p adj (dist.dir.up)' not in header: raise NameError("Gene-set statistics file has invalid header: one column should be named 'p adj (dist.dir.up).' or pAdjusted should be set as 'False'")
        if len(adj)>0 and 'p adj (dist.dir.dn)' not in header: raise NameError("Gene-set statistics file has invalid header: one column should be named 'p adj (dist.dir.dn).' or pAdjusted should be set as 'False'")
        nonPiano = False
    elif 'p (non-dir.)' in header:
        header[np.where(header=='p (non-dir.)')] = 'p-value'
        if len(adj)>0 and 'p adj (non-dir.)' not in header: raise NameError("Gene-set statistics file has invalid header: one column should be named 'p adj (non-dir).' or pAdjusted should be set as 'False'")
        if 'p adj (non-dir.)' in header:
            header[np.where(header=='p adj (non-dir.)')] = 'p-adj'
        nonPiano = True
    elif 'p-value' in header:
        if len(adj)>0 and 'p-adj' not in header: raise NameError("Gene-set statistics file has invalid header: one column should be named 'p-adj' or pAdjusted should be set as 'False'")
        nonPiano = True
    else:
        raise NameError("Gene-set statistics file has invalid header.")
    
    for stats in content[1:]:
        stats = tuple(stats)
        m = kiwiC.Metabolite(stats[0])
        m.addGeneSetStats(stats,header,adj)
        M.addMetabolite(m)
    if not any([m.pNonDirectional!=np.nan for m in M.metaboliteList]): raise NameError("Invalid data type for gene-set p-value statistic: all values are NaN")
    
    ## Retrieve minimum p-value in the GSS statistics dataset
    p_all = np.array(list(itertools.chain.from_iterable((met.pNonDirectional,met.pMixDirUp,met.pMixDirDn,met.pDistDirUp,met.pDistDirDn,met.pValue) for met in M.metaboliteList)))
    p_all[np.isnan(p_all)] = 1
    ### The minimum p_value of the dataset represents the extreme scenario. This is used for two purposes in computing the directionality score.
    ### First, it is used to add a small quantity before taking the log values (and only if there is a p-value
    ### equal to 0 in the dataset). Second, it is used to normalize the directionality score, so that +1/-1
    ### represent the extreme scenarios for the geneset (minimum pvalues in all classes in one direction, max 
    ### in the other direction). If the minimum is actually 0, a slightly higher number is used not to divide by zero.
    ### This is either an arbitrary small number or the smallest non-zero p_value in the dataset
    if p_all.min() == 0:        
        p_min_tostabilize = p_all[p_all!=0].min()
        if p_min_tostabilize > arbitrary_small:
            p_min_tostabilize = arbitrary_small
        p_min_tonormalize = p_min_tostabilize        
    else:
        p_min_tostabilize = 0
        p_min_tonormalize = p_all.min()
    
    # Remove non-significant metabolites and metabolites not in MN and metabolites w high degree:
    M.removeNotSignificantMetabolites(pcutoff,nonPiano)    
    if len(M.metaboliteList)==0:
        raise NameError('No metabolites passed the pCutoff')
    M.removeMetabolitesNotInMetNet(MN)
    if len(M.metaboliteList)==0:
        raise NameError('No more metabolites from the gene-set statistics file are present in the metabolite-metabolite network')
    if dcutoff < max(nx.degree(MN).values()):
        M.removeHighDegreeMetabolites(MN,dcutoff)
    if len(M.metaboliteList)==0:
        raise NameError('No metabolites passed the dCutoff')
    
    # Import the gene-level statistics:
    G = kiwiC.Genome()
    if len(Gfile)>0:
        content = np.genfromtxt(Gfile,dtype=None,delimiter='\t')
        header = content[0]
        if ('p' not in header) or ('FC' not in header): raise NameError("Gene-level statistics file should contain 'p' and 'FC' as column headers")
        for glstat in content[1:]:
            glstat = tuple(glstat)
            g = kiwiC.Gene(glstat[0])
            p = float(glstat[np.where(header == 'p')[0]])
            FC = float(glstat[np.where(header == 'FC')[0]])
            g.addGeneStats(p,FC)
            G.addGene(g)
    
    if len(GMfile)>0: 
        # Import the gene-metabolite information:
        content      = np.genfromtxt(GMfile,dtype=str,delimiter='\t')
        metList      = np.copy(M.metaboliteList)
        genenameList = []
        for met in metList:
            #Take the gene list from the GSC file for a metabolite
            genenames = np.unique(np.squeeze(np.asarray((content[np.where(content[:,0]==met.name)[0],1]))))
            isremoved = False   
            if lump:
            #If lumping, mets that share the same gene list will be regarded as the same met (specifically the first one
            #to be encountered). Then it checks only for the first met if it complies other stuff    
                genenamesAsList = [g for g in genenames]
                if genenamesAsList in genenameList:
                    #Check if there exist a previous metabolite that has the same gene list as the current one. If so,
                    #remove it and update label
                    M.removeMetabolite(met)
                    genenameList.append([np.nan])
                    metInd = np.where([gl == genenamesAsList for gl in genenameList])[0][0]
                    kiwiF.updateLabel(metList[metInd])
                    isremoved = True
                else:
                    genenameList.append(genenamesAsList)
            # Check if some genes were excluded by the upstream GSA from a metabolite due to missing stats
            if not isremoved:
                if len(genenames) > maxGSsize:
                    M.removeMetabolite(met)
                else:
                    for genename in genenames:
                        gene = G.getGene(genename)
                        if isinstance(gene,kiwiC.Gene):
                            met.addGene(gene)
                        else:
                            newGene = kiwiC.Gene(genename)
                            newGene.addGeneStats(p=np.nan,FC=np.nan)
                            G.addGene(newGene)
                            met.addGene(newGene)     
    
    # Check if there are at least two metabolites to draw a network
    if M.size()<=1: 
        raise NameError('Less than two metabolites passed the cutoffs. No network can be generated')

    # Construct a dense plot graph:
    PG = nx.Graph()
    PG.add_nodes_from(M.metaboliteList)
    PG.add_edges_from(itertools.combinations(PG.nodes(),2))
    
    # Calculate distance and add to edge property. Add edge weight
    for e in PG.edges():
        try:
            PG[e[0]][e[1]]['shortest_path_length'] = nx.shortest_path_length(MN,e[0].name,e[1].name)
        except nx.NetworkXNoPath:
            PG[e[0]][e[1]]['shortest_path_length'] = float('Inf')
        PG[e[0]][e[1]]['weight'] = eScaleFac/PG[e[0]][e[1]]['shortest_path_length']
                
    # Remove edges according to shortest path length:
    edges_to_remove = [e for e in PG.edges() if PG[e[0]][e[1]]['shortest_path_length']>splcutoff]
    PG.remove_edges_from(edges_to_remove)
    
    # Keep only the best edge/s for each node (if not argument sae=True):
    if not sae:
        all_edges_to_save = []
        for met in PG.nodes():
            minspl = float('Inf')
            for e in nx.edges(PG,met):
                minspl = min(minspl,PG[e[0]][e[1]]['shortest_path_length'])
            edges_to_save = [e for e in nx.edges(PG,met) if PG[e[0]][e[1]]['shortest_path_length'] == minspl]
            for e in edges_to_save:
                all_edges_to_save.append(e)
        PG.remove_edges_from([e for e in PG.edges() if e not in all_edges_to_save 
            and (e[1],e[0]) not in all_edges_to_save])
    
    # Get edge width for plotting:
    edge_width = kiwiF.getEdgeProperty(PG,'weight')
    
    # Get node attribute for plotting:
    p          = np.array([[node.pNonDirectional,node.pMixDirUp,node.pDistDirUp,
                            node.pMixDirDn,node.pDistDirDn,node.pValue] for node in PG.nodes()])
    p_stable   = p
    p_stable[np.isnan(p_stable)] = 1 # nan are transformed to 1 and in the next lines into something a bit higher than
                                     # 1, so that the -log10 is almost zero
    p_stable   = p_stable + p_min_tostabilize # Add a small number that is at most 
                                              # as high as the smallest non-zero number in the dataset
    p_stable   = -np.log10(p_stable)
    if nonPiano:
        node_size  = minNodeSize*(p_stable[:,5]+np.log10(pcutoff))+minNodeSize # p_stable[:,5] is the non-piano p-value called pValue
        color_code = p_stable[:,5]
    else:
        node_size  = minNodeSize*(p_stable[:,0]+np.log10(pcutoff))+minNodeSize # p_stable[:,0] is the piano p-value called pNonDirectional
        color_code = ((p_stable[:,1]*(p_stable[:,0]+p_stable[:,2]) - 
                p_stable[:,3]*(p_stable[:,0]+p_stable[:,4])) /
                (2*((np.log10(p_min_tonormalize))**2)))
    
    # Assign plot attributes for a node as node attributes:
    k = 0
    for node in PG.nodes(): 
        PG.node[node]['directionalityScore'] = float(color_code[k])
        if nonPiano:
            PG.node[node]['-log10p'] = float(p_stable[k,5])
        else:
            PG.node[node]['-log10p'] = float(p_stable[k,0])
        k = k+1
    
    # Plot network:
    fig_nw = plt.figure(figsize=(8,8))
    pos=nx.spring_layout(PG,iterations=50,scale=5)
    if nonPiano:
        colormap = plt.cm.Greens
        v_min = 0
    else:
        colormap = plt.cm.RdBu_r
        v_min = -abs(color_code).max()
    nx.draw(PG, pos, width=edge_width, node_size=node_size, node_color=color_code, cmap=colormap,
            vmin=v_min, vmax=abs(color_code).max(), with_labels=False)
    nx.draw_networkx_labels(PG,pos,dict([[n,n.label] for n in PG.nodes()]), font_size=labSize)
    sm = plt.cm.ScalarMappable(cmap=colormap, 
         norm=plt.Normalize(vmin=v_min, vmax=abs(color_code).max()))
    sm._A = []
    plt.colorbar(sm)
    if len(nwPlotFile)>0: 
        plt.savefig(nwPlotFile, bbox_inches='tight')
        plt.close(fig_nw)
    else:
        plt.show()
        
    # Heatmap:
    if len(GMfile)>0 and len(Gfile)>0:
        #In the same fashion for genesets, calculate a small number to stabilize the p-values upon log
        #transformation.
        p_all_gl = np.array([gene.p for gene in G.geneList])
        if p_all_gl.min() == 0:        
            pzero = p_all_gl[p_all_gl!=0].min()
            if pzero > arbitrary_small:
                pzero = arbitrary_small
        else:
            pzero = 0      
        kiwiF.drawHeatmap(PG,hmType,hmPlotFile,pzero,cgnFromType,cgnToType,cgnSpecies)
            
    # Export to graphML:
    if len(graphMLFile) > 0:
        nx.write_graphml(PG,graphMLFile)