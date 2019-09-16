#!/bin/python3
# -*- coding: utf-8 -*-
"""
Chronicle 
"""

import warnings
import numpy as np
import sys, getopt
import warnings
import scipy.sparse.csgraph
from lark import Lark
from math import exp
#from db_generator import pattern

def resize(l, n, d = None):
    while len(l) < n:
        l.append(d)

class EventMapper:
    def __init__(self):
        self.__event_map={}
        self.__rev_event_map={}
        
    def id(self, event):
        """
        return a unique identifier corresponding to the event label
        """
        idv = self.__event_map.setdefault(event, len(self.__event_map))
        self.__rev_event_map[idv]= event
        #print("create "+str(event)+" with id "+str(idv))
        return idv 
        
    def event(self, idv):
        """
        return the name of an event by its identifier
        """
        if not idv in self.__rev_event_map:
            raise KeyError('EventMapper error: unknown event with id '+str(idv)+". Known events: "+str(list(self.__rev_event_map))+".")
        else:
            return self.__rev_event_map[idv]

class Chronicle:
    """Class for a chronicle pattern modeling
    
    """
    
    lamda=0.04
    w=150
    npat = 0
    
    """
    CRS_grammar is a grammar for parsing CRS files
    """
    CRS_grammar=r"""start: chronicle+

chronicle: "chronicle" NAME "()" "{" event+ constraint* "}"

event: "event" "(" NAME "," ID ")"
constraint: ID "-" ID "in" INTERVAL

INTERVAL: "[" NUMBER "," NUMBER "]"
ID: "t" NUMBER
NAME: CNAME ["[]"]
WHITESPACE: (" " | "\t" | "\n")+
%ignore WHITESPACE
%import common.SIGNED_NUMBER    -> NUMBER
%import common.CNAME
"""

    def __init__(self, emapper=None):
        """
        - emapper is an event mapper, if not provided, a new one is created
        """
        
        self.tconst={}  #temporal constraints,
                        # keys: couple (ei,ej) where ei is a index in the item
                        #   in the multiset
                        # values; couple (lb,ub)
        self.inconsistent = False
        self.name = ""
        self.sequence = {}      # description of the pattern events
        self.tidl=[]    # list of the transactions in which the pattern occurs
        self.pid=Chronicle.npat   # pattern id
        Chronicle.npat += 1
	    
        
        if not emapper:
            self.emapper = EventMapper()
        else:
            self.emapper = emapper

    def add_event(self, pos, event):
        """Add an event to the chronicle multiset
        Contrary to add_item, an integer is not required to denote an event!
        """
        self.__add_item(pos, self.emapper.id(event) )
        
    def __add_item(self, pos, item):
        """Add an item to the chronicle
        The function creates all infinite constraints, without variability
        - the id of the event correspond to the order of added items
        """
        self.sequence[pos] = item
        for i in range(pos):
            if not (i,pos) in self.tconst:
                if i in self.sequence and self.sequence[i]==item:
                    self.tconst[(i,pos)]= (1,float("inf")) #here: 1 means that the same items must occur after!
                else:
                    self.tconst[(i,pos)]= (-float("inf"),float("inf"))
        for i in range(pos+1,max(self.sequence.keys())+1):
            if not (pos,i) in self.tconst:
                if i in self.sequence and self.sequence[i]==item:
                    self.tconst[(pos,i)]= (1,float("inf"))
                else:
                    self.tconst[(pos,i)]= (-float("inf"),float("inf"))
        
    def add_constraint(self, ei, ej, constr):
        """Add a constraint-template to the chronicle pattern
        - ei, ej: index of the events in the multiset
        - constr: a 2-tuple (min,max)
        """
        if not type(constr) is tuple:
            print ("error: constraint must be a tuple (=> constraint not added)")
            return
            
        if len(constr)!=2:
            print ("error: constraint must have 2 values (=> constraint not added)")
            return
            
        try:
            self.tconst[(ei,ej)] = constr
        except IndexError:
            print ("add_constraint: index_error (=> constraint not added)")
        
    def __getitem__(self, i):
        """return the item at position i in the multiset if i is an integer
        and return the constraint between i[0] and i[1] if i is a couple
        """
        if not type(i) is tuple:
            return self.sequence[i]
        else:
            try:
                return self.tconst[(min(i[0],i[1]),max(i[0],i[1]))]
            except KeyError:
                return (-float("inf"),float("inf"))
        
            
    def __len__(self):
        """ Length of the patterns (number of items)
        """
        if not self.sequence:
            return 0
        return max(self.sequence.keys())+1
        
    def __str__(self):
        """
        s = "C"+str(self.pid)+": {"+str(self.sequence) + "}\n"
        s += "\t "+str(self.tconst) + "\n"
        """
        s = "C"+str(self.pid)+"\t {{"+ ','.join([str(self.emapper.event(v)) for k,v in self.sequence.items()]) + "}}\n"
        for i in self.sequence.keys():
            for j in [v for v in self.sequence.keys() if v>i]:
                s += str(i) + "," + str(j) + ": " + str(self.tconst[(i,j)])+"\n"
        s +="tid: "+ str(self.tidl)
        return s
        
    def delete(self, itempos):
        del self.sequence[ itempos ]
        todelete=[]
        for k in self.tconst:
            if k[0]==itempos or k[1]==itempos :
                todelete.append(k)
        for k in todelete:
            del self.tconst[k]

    def minimize(self):
        #construction of distance graph
        mat=np.matrix( np.zeros( (max(self.sequence.keys())+1,max(self.sequence.keys())+1) ))
        for i in range(max(self.sequence.keys())+1):
            for j in range(i+1,max(self.sequence.keys())+1):
                if (i,j) in self.tconst:
                    mat[i,j] = self.tconst[ (i,j) ][1]
                    mat[j,i] = -self.tconst[ (i,j) ][0]
                else:
                    mat[i,j] = float("inf")
                    mat[j,i] = -float("inf")
        try:
            matfw = scipy.sparse.csgraph.floyd_warshall( mat )
            #construction of simplified chronicle
            for i in range(max(self.sequence.keys())+1):
                for j in range(i+1,max(self.sequence.keys())+1):
                    self.tconst[ (i,j) ] = (- int(matfw[j,i]), int(matfw[i,j]))
        except:
            warnings.warn("*** Minimisation: Inconsistent chronicle ***")
            self.inconsistent = True
    ################
    
    def __CRS_read_tree(tree, chronicle=None, id_map={}):
        if tree.data =="start":
            return Chronicle.__CRS_read_tree(tree.children[0], chronicle, id_map)
        elif tree.data == "chronicle":
            if not chronicle:
                C = Chronicle()
            else:
                C = chronicle
            print(id_map)
            C.name = str(tree.children[0][:-2]) #remove the last two characters '[]'
            for i in range(1,len(tree.children)):
                Chronicle.__CRS_read_tree(tree.children[i],C, id_map)
            return C
        elif tree.data=="event":
            event = str(tree.children[0])
            event = event.strip("[]") #remove the '[]' if necessary
            eid = id_map.setdefault(str(tree.children[1]), len(id_map))
            chronicle.add_event(eid, event)
        elif tree.data=="constraint":
            eid1=id_map[str(tree.children[0])]
            eid2=id_map[str(tree.children[1])]
            interval=str(tree.children[2]).strip('[]').split(',')
            if eid1<eid2 :
                chronicle.add_constraint(eid1,eid2, (-float(interval[1]), -float(interval[0])))
            else:
                chronicle.add_constraint(eid2,eid1, (float(interval[0]), float(interval[1])))
        
    def load(crs, emapper=None):
        """Load a chronicle from a string in the CRS format.
        Note that the all brackets ("[]" in chronicle or events names; and "()") are assumed to be empty in this function !!!
        
        This is a class-function.
        
        parameters:
        - crs: string describing a string in a CRS format
        - emapper (optional): an external event mapper
        
        return the newly instantiated chronicle
        """
        chro_parser = Lark(Chronicle.CRS_grammar)
        tree= chro_parser.parse(crs)
        if not emapper:
            return Chronicle.__CRS_read_tree(tree, id_map={})
        else:
            C = Chronicle(emapper)
            return Chronicle.__CRS_read_tree(tree, C, {})
            
            
    def to_crs(self):
        """Generate a string representing the chronicle in the CRS format.
        
        Unnamed events (must be figures) are called "E"+str(X) in the event description to avoid events name starting with figures (CNAME conventions)
        Infinite intervals are not printed out, but semi-infinite intervals will generate an description like '[-inf,23]', or '[34,inf]' : do not know whether it is sound or not!
        
        return a string
        """
        s="chronicle "
        if self.name!="":
            s+=str(self.name)
        else:
            s+="C"+str(self.pid)
        s+="[]()\n{\n"

        for pos,e in self.sequence.items():
            if self.emapper:
                evt = self.emapper.event(e)
                if isinstance(evt, str):
                    s+="\tevent("+evt+"[], t{:03d})\n".format(pos)
                else:
                    s+="\tevent("+str(evt)+"[], t{:03d})\n".format(pos)
            else:
                s+="\tevent("+str(e)+"[], t{:03d})\n".format(pos)
        s+="\n"
        
        for events,interval in self.tconst.items():
            if interval[0]!=float("-inf") or interval[1]!=float("inf"): #infinite intervals are not printed out
                s+="\tt{:03d}-t{:03d} in [{},{}]\n".format(events[0],events[1],interval[0],interval[1])
        s+="}"
        return s 
                
    
    ################
    def find(self,a,p,i):
        l=[item for item in a if item == (p,i) ]
        return len(l)!=0

    def __complete_recognition__(self, occurrence, item_index, sequence, proo):
        """
        return a list of occurrences that add the description of the matching of the item_index-th item of the chronicle
        to the occurrence
        """
        
        if not item_index in self.sequence:
            return [occurrence],proo
        
        #if(occurrence[item_index][0]> occurrence[item_index][1]):
        #    return [],0
        item=self.sequence[item_index]
        tmax=sequence[len(sequence)-1][0]     
        if occurrence[item_index][0]==occurrence[item_index][1]:
            if occurrence[item_index][0]<tmax and self.find(sequence,occurrence[item_index][0],item)==True:
                if not(type(proo) is list):
                    proo=[proo]
                return [occurrence], proo
            else:
                return [],[]
        
        occurrences = []
        prob=[]
        #assert(occurrence[item_index][1]<len(sequence))
        w=Chronicle.w
        interval=[occurrence[item_index][0],occurrence[item_index][1]]
        pos=occurrence[:][item_index-1][0]
        borne=self.__borne(pos,interval,w,tmax)
        i=self.__startReaserch(sequence,borne[0])
        if(i is None):
            return [],[]
        while i< len(sequence ) and sequence[i][0] <= borne[1] :
            if sequence[i][1]==item:
                p=sequence[i][0]
                #create a new occurrence to be modified
                pp=self.__prob_likelihood([occurrence[item_index][0],occurrence[item_index][1]],p)
                new_occ = occurrence[:]
                new_occ[item_index] = (p,p)
                pp=proo*pp
                #if(p< occurrence[item_index][0] or p>occurrence[item_index][1] ):
                #    prob=prob*self.pp
                
                satisfiable=True
                #propagate chronicle constraints
                for k in self.tconst:
                    v = self.tconst[k]
                    if (k[0]==item_index) and (k[1] in self.sequence):
                        #new_occ[ k[1] ] = (max(new_occ[ k[1] ][0], p+v[0]), min(new_occ[ k[1] ][1], p+v[1]))
                        new_occ[ k[1] ] = ( p+v[0], p+v[1])
                        if new_occ[ k[1] ][0]>new_occ[ k[1] ][1]: #if empty interval, it is not satisfiable
                            #new_occ[ k[1] ] = ( p+v[0], p+v[1]))
                            satisfiable=False
                            break
                        
                if satisfiable:
                    #add the occurrence to the list
                    occurrences.append( new_occ )
                    prob.append(pp)
            i +=1
        return occurrences, prob
    
    
    def __recrecognize__(self, occurrence, last_item_index, sequence,prob_pre,seuil):
        """
        recursive call for occurrence recognition
        return a list of occurrences recognized from the last_item_index of the chronicle until its last item
        """
        chro_size=max( self.sequence.keys() )
        if last_item_index==chro_size:
            return [occurrence], prob_pre
        
        item_index=last_item_index+1
        proba =[]
        occurrences = []
        loc_occs , loc_prob= self.__complete_recognition__(occurrence, item_index, sequence, prob_pre)
        loc_occs =[loc_occs[i] for i in range(len(loc_occs)) if loc_prob[i]>seuil]
        for occ in loc_occs:
           reoccs, prob = self. __recrecognize__(occ, item_index, sequence, loc_prob[loc_occs.index(occ)],seuil) 
           occurrences.extend(reoccs)
           if(type(prob) is list):
                proba.extend(prob)
           else:
                proba.append(prob) 
        #self.prob /= self.pp
        #if(len(occurrences)==0):
        #    self.prob=0
        
        return occurrences ,proba
                
    def recognize(self,sequence,seuil=0):
        """
        Method that checks whether the chronicle occurs in the sequence 
        sequence: list of events
        Return a list of occurrences
        """
        return self.__recognize([(timestamp,self.emapper.id(event))  for timestamp, event in sequence],seuil)
    
    def __recognize(self,sequence,seuil):
        """
        Method that checks whether the chronicle occurs in the sequence 
        sequence: list of event identifiers
        Return a list of occurrences
        """
        if(len(sequence)==0):
            return [],0
        occurrences = [] #list of occurrences
        proba=[]
        
        chro_size=max( self.sequence.keys() )+1
        if chro_size==0 :
            return occurrences
        
        item_index = 0
        item=self.sequence[item_index]
        #seq_len = len(sequence)
        #self.minimize()
        tmax=sequence[len(sequence)-1][0]
        for p,e in sequence:
            if e==item:
                #create a new occurrence
                new_occ = []
                resize(new_occ, chro_size, (0,tmax))
                new_occ[item_index] = (p,p)

                #propagate chronicle constraints
                for k in self.tconst:
                    v = self.tconst[k]
                    if (k[0]==item_index) and (k[1] in self.sequence):
                        new_occ[ k[1] ] = (max(0,p+v[0]), min(p+v[1],tmax))
                
                #ajouter l'occurrence à la liste des occurrences
                #self.prob=1
                loc_occ , loc_prob  = self.__recrecognize__(new_occ, item_index, sequence,1,seuil)
                occurrences.extend( loc_occ )
                if(type(loc_prob) is list):
                    proba.extend(loc_prob)
                else:
                    proba.append(loc_prob) 
        if(len(proba)!=0):
            pmax=max(proba) 
            occurrences=[occurrences[i] for i in range(len(occurrences)) if proba[i]==pmax]
        else:
            pmax=0
        return occurrences,pmax
        #occ=proba.index(pmax)
        #return occurrences[occ], pmax  
    

    def start_event(self):
        evnt=[events[0]  for events,interval in self.tconst.items() if  interval[0]!=float("-inf") or interval[1]!=float("inf")]    
        return max(evnt,key=evnt.count)
        
    def __prob_likelihood(self,interval,t):
        lamda=Chronicle.lamda
        if(interval[0] <=t and interval[1] >=t ):
            p=1
        elif (interval[0]>interval[1]):
            p=exp(-lamda*(abs(t-interval[0])))
        else:
            p=exp(-lamda*(min(abs(t-interval[0]),abs(t-interval[1]))))
        return p
    
    def __borne(self,pos,interval,w,len_seq):
        borne=interval[:]
        if(interval[0]>pos):
            #borne[0]=max(interval[0]-w,pos +1)
            borne[0]=max(interval[0]-w,pos)
        elif(interval[0]<pos):
            borne[0]=max(interval[0]- w,0)
        else:
            borne[0]=pos
            
        if(interval[1]<pos):
            borne[1]=min(pos -1,interval[1]+w)
        elif(interval[1]>pos):
            borne[1]=min(interval[1]+1+w,len_seq)
        else:
            borne[1]=pos
            
        return borne
            
    def __ToList(self,liste):
        nelement = []
        for element in liste :
            nelement.extend(element)
        return nelement
        
    def BestRecognition(self,occ):
        mx=max(occ[1]) #prob_occ
        return occ[0][occ[1].index(mx)],mx
    
    def __startReaserch(self,sequence,tmin):
        for t,e in sequence:
            if(t>=tmin):
                return sequence.index((t,e))
        
                
if __name__ == "__main__":
    """
    seq=[(160, 16), (308, 2), (556, 19)]
    print("sequence: "+str(seq))
    c=Chronicle()
    #print(c)

    c.add_event(0,16)
    c.add_event(1,16)
    #c.add_event(2,5)
    c.add_constraint(0,1, (-float('inf'), float('inf')))
    #c.add_constraint(0,2, (-40.07730049601006, 109.52508766955438))
    #c.add_constraint(1,2, (67.4333242689155, 258.7202212149524))
    print(c)
    occs=c.recognize(seq)
    print("reco: "+str(occs[0]))
    print(occs[1])

    seq=[(1,1),(0,2),(2,3),(0,4),(0,5),(0,6)]
    c=Chronicle()
    c.add_event(0,1)
    c.add_event(1,2)
    c.add_event(2,2)
    c.add_constraint(0,1, (1,2))
    c.add_constraint(1,2, (-float('inf'),float('inf')))
    #c.add_constraint(1,2, (67.4333242689155, 258.7202212149524))
    print(c)
    occs=c.recognize(seq)


    seq=[(173, 2), (294, 0), (376, 5), (429, 13)] 
    c=Chronicle()
    c.add_event(0,5)
    c.add_event(1,13)
    c.add_constraint(0,1,(58.61928024080464, 114.51122038387246))
    occ=c.recognize(seq)
    print(occ)

    seq= [(245, 2), (413, 1), (493, 5)] 
    c=Chronicle()
    c.add_event(0,1)
    c.add_event(1,2)
    c.add_event(2,5)
    c.add_constraint(0,2,(-60.181163021013596, -5.822358418911001))
    print(c)
    occ=c.recognize(seq)
    print(occ)
    """
  
