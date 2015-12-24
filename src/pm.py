#!/usr/bin/env python
#  Name:  
#    pm.py
#
#  Purpose:  
#    Perform acquisition path analysis on the physical model of a nuclear fuel cycle  
#
#  Usage:        
#        Calculate optimal inspection strategies 
#        from acquisition path analysis
#        
#        python %s [-h] [-d] [-b value] Inspection_effort
#        -h this usage message
#        -d draw diversion path graph
#        -b set perceived sanction parameter
#           (incentive d is 9 for most attractive path)
#
#  Copyright (c) 2015, Mort Canty, Clemens Listner
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
import networkx as nx
import matplotlib.pyplot as plt
import itertools, time
import numpy as np
from itertools import combinations
from nash.nasheq import nashEquilibria

plt.rcParams['figure.figsize'] = 12, 10

def undominated(K,Ks):
    for KK in Ks:
        if set(K) <= set(KK):
            return False 
    return True

def APA(ps,ls,betas,alphas,ws,W,b=1,f=1,a=10,c=100,e=1,select='one',s1=0,verbose=True):  
    '''   Acquisition path analysis game theoretical effectiveness evaluation
    Usage:
     python APA(ps,ls,betas,alphas,ws,W,b=10,f=1,a=10,c=100,e=1) 
     ps =    list of acquisition paths. Each path consists 
             of a list of edges (illegal activities connecting two nodes or material types)
     ls =    list of corresponding path lengths
     betas = list of non-detection probabilities in event of  
            violation and inspection at a location/activity
     alphas = false alarm probability in event of inspection of 
              location/activity where State has behaved legally
     ws =    list of inspection efforts (person-days) for each 
             activity/location
     W  =    available inspection effort (person-days)
     b =     perceived sanction in event of detection relative 
             to the payoff 10 for undetected violation 
             along the shortest path.   
     f =     false alarm resolution cost for any location/activity 
             on same scale as above
     k =     number of locations/activities inspected
     a,c,e = IAEA utilities (-a = loss for detection, 
             -c = loss for undetected violation, 
             -e loss incurred by false alarm)
     s1 = use s1+1 almost complementary path for Lehmke-Howson
     All utilities are normalized to the utility 0 for 
     legal behavior (deterrence) for both players     
     Returns equilibria as a list of tuples [(P*,H1*,Q*,H2*), ...]
     (c) Mort Canty 2015'''

    ls = np.asarray(ls)
    ps = np.asarray(ps)
    ws = np.asarray(ws) 
    betas = np.asarray(betas)
    alphas = np.asarray(alphas)             
    n = len(ps)
    m = len(betas)
#  pure strategies of Inspectorate    
    Ks = []
    if verbose:
        print 'Enumerating feasible undominated inspection strategies ...'
    idx = np.argsort(ws)
    cs = np.cumsum(ws[idx])
#  largest possible number of activities   
    kl = 0
    while kl < len(cs) and cs[kl] <= W:
        kl += 1  
    if verbose:      
        print 'Largest inspection strategy covers %i activities'%kl       
    for i in range(kl,0,-1):
        cm = combinations(range(m),i)
        while True:
            try:
                K = list(cm.next())
                if (sum(ws[K]) <= W) and undominated(K,Ks):
                    Ks.append(K)
            except:
                break  
        if verbose:  
            print 'length: %i cumulative inspectorate strategies: %i' \
                                                        %(i,len(Ks))               
    mk = len(Ks)
#  path preferences inversely proportional to path lengths    
    ds = 1./ls
#  path preferences in decreasing order beginning with 9    
    ds = 9*ds/ds[0]
#  construct bimatrix
    if verbose:
        print 'Building bimatrix...'
    A = np.zeros((mk,n+1))
    B = np.zeros((mk,n+1))
    i = 0
    for K in Ks:
        K = set(K)
        for j in range(n+1):
            if j<n: 
                P = set(ps[j])
            else:
                P = set([])    
            beta = 1
            for ell in K & P:
                beta *= betas[ell]
            tmp = 1
            for ell in K - (K & P):
                tmp *= 1 - alphas[ell]
            alpha = 1 - tmp 
            if j < n:         
                A[i,j] = -c*beta - a*(1-beta) - e*alpha
                B[i,j] = (ds[j])*beta - b*(1-beta) -f*alpha
            else:
                A[i,j] = -e*alpha
                B[i,j] = -f*alpha         
        i += 1
#  solve the game     
    if verbose:
        print 'Calling nashEquilibria on %i by %i bimatrix ...'%(mk,n+1)          
    return (Ks,nashEquilibria(A.tolist(),B.tolist(),select=select,s1=s1))


class Pm(object):
    '''
    A Physical Model class for acquisition path analysis. Initialize
    with 'fuelcycle', the list of edges appropriate for the safeguarded
    nuclear fuel cycle of a specific NPT member state.
    '''
    def __init__(self,fuelcycle,name='Physical Model'):
        self.template = {
                                'Origin':(1,10),
                                           'SourceMaterialResources':(2,9),
                                           'SourceMaterial':(2,8),'EnrichmentFeed':(3,7.5),
                                                                  'IUEnrichmentProduct':(3,6.5),       'DUEnrichmentProduct':(5,6),
                                           'NUFuelFeed':(2,6),     'IUFuelFeed':(3,5.5),            'DUFuelFeed':(5,5),  
                                           'NUFuel':(2,5),         'IUFuel':(3,4.5),                'DUFuel':(5,4), 
                                                                   'IrradiatedFuel':(3,3.5),
                                         'IUReprecessedMaterial':(2,3),                 'DUReprocessedMaterial':(4,2.5)                                               
                         }
        self.name = name
        self.G = nx.DiGraph()
        self.G.add_nodes_from(self.template.keys())
        numedges = len(fuelcycle)
        self.G.add_edges_from(fuelcycle)
        paths1 = nx.all_simple_paths(self.G,'Origin','DUEnrichmentProduct')
        paths2 = nx.all_simple_paths(self.G,'Origin','DUFuelFeed')
        paths3 = nx.all_simple_paths(self.G,'Origin','DUFuel')
        paths4 = nx.all_simple_paths(self.G,'Origin','DUReprocessedMaterial')
        pathiterator = itertools.chain(paths1,paths2,paths3,paths4)
        paths = []
        pathlengths = []
        for path in pathiterator:
            paths.append(path)
            edges = zip(path,path[1:])
            length = 0        
            for e in edges:
                length += self.G.edge[e[0]][e[1]]['cost']
            pathlengths.append(length) 
        tmp = zip(paths,pathlengths) 
        tmp.sort(key=lambda p:p[1]) 
        self.paths = tmp 
        self.ps = []
        self.ls = []
        for p in self.paths:
            self.ls.append(p[1])
            edges = zip(p[0],p[0][1:])
            idxpath = []
            for e in edges:
                idxpath.append(self.G.edge[e[0]][e[1]]['idx']) 
            self.ps.append(idxpath) 
        usededges = set([item for sublist in self.ps for item in sublist])                 
        self.betas = np.zeros(numedges); self.alphas = np.zeros(numedges); self.ws = np.zeros(numedges); self.act = []
        for e in self.G.edges():
            idx = self.G.edge[e[0]][e[1]]['idx']
            self.betas[idx] = self.G.edge[e[0]][e[1]]['beta']
            self.alphas[idx] = self.G.edge[e[0]][e[1]]['alpha']
            self.ws[idx] = self.G.edge[e[0]][e[1]]['w']
            self.act.append(self.G.edge[e[0]][e[1]]['act']+' - ')
        self.ws[list(set(range(numedges))-usededges)] = 1000000
        
    def listedges(self):
        idxs = []
        acts = []
        for e in self.G.edges():
            idxs.append(self.G.edge[e[0]][e[1]]['idx'])
            acts.append(self.G.edge[e[0]][e[1]]['act'])
        tmp = zip(idxs,acts)
        tmp.sort(key=lambda p:p[0])
        print 'Illegal activities for %s'%self.name
        for t in tmp:
            print '%i:%s'%t
            
            
    def _display(self):
        nx.draw_networkx(self.G,self.template,edge_color='lightgray',node_color='lightgray')
        plt.title(self.name)
        
    def display(self):
        self._display()
        plt.show()
        
    def displaypath(self,path,color='red'):
        self._display()
        pathedges = zip(path[0],path[0][1:])
        nx.draw_networkx_edges(self.G,self.template,edgelist=pathedges,edge_color=color)
        plt.xlabel('path length: %s'%str(path[1]))
        plt.show()
        
    def solve(self,W,b=1,c=100,verbose=True):
        start = time.time()
        Ks, eqs = APA(self.ps,self.ls,self.betas,self.alphas,self.ws,W,b=b,c=c,verbose=verbose)
        if verbose:
            print 'elapsed time: %s'%str(time.time()-start)
            eq = eqs[0]
            P = np.array(eq[0])*1000
            H1 = eq[1]
            Q = np.array(eq[2])*1000
            H2 = eq[3]
            print 'P:'
            P = np.array(map(round,P))/1000.
            for i in range(len(Ks)):
                if P[i] != 0:
                    print P[i], Ks[i]
            print 'Q:'
            Q = np.array(map(round,Q))/1000.
            for i in range(len(self.ps)):
                if Q[i] != 0:               
                    print  Q[i], self.ps[i]
            print 'H1 = %f, H2 = %f'%(H1,H2)   
        return (Ks,eqs)
        
    def simulate(self):
    # simulation   
        np.random.seed(222)
        start = time.time()   
        edges = 10
        betas = np.random.rand(edges)*0.5 
        alphas = np.zeros(edges)+0.05
        ws = np.random.random(edges)
        W = 1.4
        ps = [(1,2,5,8),(1,3),(2,6,7),(2,3,7,9),(2,4,6,8),(1,4,6,7),(4,5,6)]
        ls = np.array([1,3,4,3,4,1,2])
        print 'Vertex Enumeration:--------------'
        Ks, eqs = APA(ps,ls,betas,alphas,ws,W,select='all')
        print 'There are %i extreme equilibria in all'%len(eqs)
        k = 1
        for eq in eqs:
            print 'equlibrium %i --------------'%k
            P = eq[0]
            H1 = eq[1]
            Q = np.array(eq[2])*100
            H2 = eq[3]
            print 'P:'
            for i in range(len(Ks)):
                if P[i] != 0:
                    print Ks[i], P[i]
            print 'Q:'
            print np.array(map(round,Q))/100.
            print 'H1 = %f, H2 = %f'%(H1,H2)
            k += 1
        print 'elapsed time: %s'%str(time.time()-start) 
        print ''
        print 'Lemke-Howson:------------------' 
        etimes = []
        for s1 in range(len(Ks)+len(ps)+1):
            start = time.time()
            Ks, eqs = APA(ps,ls,betas,alphas,ws,W,select='one',s1=s1,verbose=False)
            etime = time.time()-start
            etimes.append(etime)
            print '%i-almost complementary start'%(s1+1)
            eq = eqs[0]
            P = eq[0]
            H1 = eq[1]
            Q = np.array(eq[2])*100
            H2 = eq[3]
            print 'P:'
            for i in range(len(Ks)):
                if P[i] != 0:
                    print Ks[i], P[i]
            print 'Q:'
            print np.array(map(round,Q))/100.
            print 'H1 = %f, H2 = %f'%(H1,H2)   
            print 'elapsed time: %s'%str(time.time()-start)  
        ya,xa = np.histogram(np.array(etimes),bins=50)
        print 'min: %f'%np.min(etimes)
        print 'max: %f'%np.max(etimes)
        print 'mean: %f'%np.mean(etimes)
#        ya[0] = 0    
#        plt.plot(xa[0:-1],ya)
#        plt.show() 

if __name__ == '__main__':   
# test
    edges = [
              ('Origin','SourceMaterial',                   {'cost':1,'beta':0.05,'alpha':0,'w':10,'idx':0,'act':'SourceMaterialDiversion'}),
              ('Origin','IrradiatedFuel',                   {'cost':1,'beta':0.05,'alpha':0,'w':10,'idx':1,'act':'IrradiatedFuelDiversion'}),
              ('Origin','IUEnrichmentProduct',              {'cost':1,'beta':0.05,'alpha':0,'w':10,'idx':2,'act':'IUEnrichmentProductlDiversion'}),
              ('Origin','EnrichmentFeed',                   {'cost':1,'beta':0.05,'alpha':0,'w':10,'idx':3,'act':'EnrichmentFeedlDiversion'}),
              ('Origin','SourceMaterialResources',          {'cost':1,'beta':0.05,'alpha':0,'w':10,'idx':4,'act':'SourceMaterialResourcesDiversion'}),
              ('Origin','IUFuelFeed',                       {'cost':1,'beta':0.05,'alpha':0,'w':10,'idx':5,'act':'IUFuelFeedDiversion'}),
              ('Origin','IUFuel',                           {'cost':1,'beta':0.05,'alpha':0,'w':10,'idx':6,'act':'IUFuelDiversion'}),
              ('Origin','DUFuel',                           {'cost':1,'beta':0.05,'alpha':0,'w':10,'idx':7,'act':'DUFuelDiversion'}),
              ('SourceMaterialResources','SourceMaterial',  {'cost':1,'beta':0.1,'alpha':0,'w':50,'idx':8,'act':'MiningMillingMisuse'}),
              ('SourceMaterial','NUFuelFeed',               {'cost':1,'beta':0.1,'alpha':0,'w':50,'idx':9,'act':'ConversionMisuse'}),
              ('SourceMaterial','EnrichmentFeed',           {'cost':1,'beta':0.1,'alpha':0,'w':50,'idx':10,'act':'Conversion1Misuse'}),
              ('EnrichmentFeed','IUEnrichmentProduct',      {'cost':1,'beta':0.1,'alpha':0,'w':50,'idx':11,'act':'EnrichmentPlantMisuse'}),
              ('EnrichmentFeed','DUEnrichmentProduct',      {'cost':1,'beta':0.1,'alpha':0,'w':200,'idx':12,'act':'CascadeMisuse'}),
              ('EnrichmentFeed','NUFuelFeed',               {'cost':1,'beta':0.1,'alpha':0,'w':200,'idx':13,'act':'ReconversionMisuse'}),
              ('IUEnrichmentProduct','IUFuelFeed',          {'cost':1,'beta':0.1,'alpha':0,'w':100,'idx':14,'act':'FuelFabMisuse'}),
              ('IUEnrichmentProduct','DUEnrichmentProduct', {'cost':1,'beta':0.1,'alpha':0,'w':1000,'idx':15,'act':'EnrichmentPlantMisuse1'}),
              ('NUFuelFeed','NUFuel',                       {'cost':1,'beta':0.1,'alpha':0,'w':100,'idx':16,'act':'FuelFabMisuse1'}), 
              ('NUFuelFeed','EnrichmentFeed',               {'cost':1,'beta':0.1,'alpha':0,'w':100,'idx':17,'act':'NUFuelFeedEnrichment'}), 
              ('IUFuelFeed','IUEnrichmentProduct',          {'cost':1,'beta':0.1,'alpha':0,'w':100,'idx':18,'act':'IUFuelFeedEnrichment'}),
              ('IUFuelFeed','IUFuel',                       {'cost':1,'beta':0.1,'alpha':0,'w':100,'idx':19,'act':'FuelFabMisuse'}),
              ('NUFuel','NUFuelFeed',                       {'cost':1,'beta':0.1,'alpha':0,'w':100,'idx':20,'act':'NUFuelFabMisuse'}),
              ('IUFuel','IUFuelFeed',                       {'cost':1,'beta':0.1,'alpha':0,'w':100,'idx':21,'act':'IUFuelFabMisuse'}),
              ('IUFuel','IrradiatedFuel',                   {'cost':1,'beta':0.1,'alpha':0,'w':100,'idx':22,'act':'ReactorMisuse'}),
          ]        
    pm = Pm(edges,'German Fuel Cycle')
#    pm.listedges()
#    pm.display()
#    Ks, eqs = pm.solve(50)
    print '---------------Simulation --------------'
    pm.simulate()
    