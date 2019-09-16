from  recognition.chronicle_recognition import Chronicle
from seq_generation.chronicle_generator import *
from math import pi,sqrt,exp
from random import uniform,sample,shuffle,gauss
import queue
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import as_completed
import concurrent.futures
from sklearn.model_selection import train_test_split
import pandas as pd
import numpy as np
from sklearn import metrics

def disturbance(c,sd):
    for k in c.tconst:
        v=c.tconst[k]
        if(k[0] != k[1] and v[0] != -float('inf') and v[1]!= float('inf')):
            #mean=(k[1]-k[0])/2
            #sd=5 #k[1]-k[0]
            #b1=gauss(v[0],sd)  #gaussian_probability(v[0]-sd+2*random.random()*sd,v[0],sd)
            #b2=gauss(v[1],sd)  #gaussian_probability(v[1]-sd+2*random.random()*sd,v[1],sd)
            #mean=int((v[0]+v[1])/2)
            b1=uniform(v[0]-sd,v[0])
            b2=uniform(v[1],v[1]+sd)
            c.add_constraint(k[0],k[1],(b1,b2))
    return c

def  make_noise(DB,per=0.3,sd=300):
    #c,p =train_test_split(DB,test_size=per)
    c,p=split_db_ch(DB,per)
    for i in range(len(p)):
        p[i]=disturbance(p[i],sd)
    return p,c


def affectation(ch_gen):
    c=Chronicle()
    i=0
    for e in ch_gen.sequence:
        c.add_event(i,e)
        i+=1
    for k in ch_gen.tconst:
        if(k[0] != k[1]):
            v=ch_gen.tconst[k]
            c.add_constraint(k[0],k[1],v)
    c.tidl=ch_gen.tidl
    return c


def split_db_ch(data,per=0.3):
    l=int(len(data)*per)
    if(l==0):
        l=1
    db=data[:]
    disturb_data=[]
    shuffle(db)
    i=0
    while len(disturb_data)<l and i<len(data):
        if (not check_bounds(db[i])):
            disturb_data.append(db[i])
            db.remove(db[i])
        i +=1
    nn_disturb_data = db[:]
    return nn_disturb_data,disturb_data

def check_bounds(ch):
     return {e for e in ch.tconst.values()}=={(-float('inf') , float('inf'))}

def nsplit_db(data,n):
    """
     Yield successive n-sized chunks from data.
     """
    k, m = divmod(len(data), n)
    return list(data[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))


def sort_db_ch(DB_ch):
    DB_ch.sort(key = lambda v: len(v.tidl),reverse=True)
    return DB_ch
       

def KNN(DB_ch,seq):
    occs=0
    chro=[]
    for ch in DB_ch:
        v=ch.recognize(seq,occs)
        if(v[1]>occs):
            occs=v[1]
            chro=ch
        if(occs==1):
            return 1,chro,seq 
    return round(occs,4),chro,seq


def KNN_op(DB_ch,seq,d):
    occs=0
    chro=[]
    for ch in DB_ch:
        if len(d.queue)!=0 and max(list(d.queue))==1 :
            return 0,
        v=ch.recognize(seq,occs)
        if(v[1]>occs):
            occs=v[1]
            chro=ch
        if(occs==1):
            return 1,chro,seq 
    return occs,chro,seq

def isTreadAlive(threads):
  for t in threads:
    if t.isAlive():
      return 1
  return 0

def worker(DB_ch,seq,d):
    d.put(KNN_op(DB_ch,seq,d)[0])

def KNN_Multi(DB_ch,seq,nth=5):
    """
    Parallel research for a sequence in the chronicle database 
    """
    d=queue.Queue()
    threads = []
    num_worker_threads=nth
    chunks=nsplit_db(DB_ch,num_worker_threads)
    for i in range(num_worker_threads):
        t = threading.Thread(target=worker,args=[chunks[i],seq,d])
        t.start()
        threads.append(t)
    #alive=0
    #while(max(list(d.queue))!=1 and alive!= 0) :
    #    alive=isTreadAlive(threads)
    for t in threads:
        t.join()
    return max(list(d.queue))


def query_consumer(queries,d,stop_flag,DB_ch):
    while not queries.empty():
        
        seq = queries.get()
        d.put(KNN(DB_ch,seq)[0])
        #d.put(KNN_Multi(DB_ch,seq))
        if stop_flag or seq==None :
            break
        queries.task_done() 

        
   
def parallel_query_consumer(DB_ch,queries,i=10):
    d=queue.Queue()
    threads = []
    stop_flag=False
    num_worker_threads=i
    for i in range(num_worker_threads):
        t = threading.Thread(target=query_consumer,args=[queries,d,stop_flag,DB_ch])
        t.start()
        threads.append(t)


    # block until all tasks are done
    queries.join()
    # stop workers
    stop_flag=True
    for t in threads:
        t.join()
    
    return list(d.queue)

def syn(DB_ch,DB_seq,i):
    with ThreadPoolExecutor(max_workers = i) as executor:
        results = [executor.submit(KNN,DB_ch,seq) for seq in DB_seq]
        concurrent.futures.wait(results)
    return [ r.result()[0] for r in  results]


def producer(DB_seq):
    q = queue.Queue()
    source=DB_seq
    for item in source:
        q.put(item)
    return q


def predict(DB_ch,l):
    pred=[KNN(DB_ch,seq)[0] for seq in l]
    return pred
def decision(l,th=0.9):
    return [int(i>=th) for i in l]


if __name__ == "__main__":
    """

    print("====================\nChronicle DB generation\n====================\n")
    generator=chrodb_generator(nbitems=20,l=3, lp=3)
    sequences = generator.generate(nb=6, l=6, npat=3, th=.1)

    print("======== PATTERNS =======")
    DB_c=[ch for ch in generator.all_patterns()]
    DB_ch=[affectation(ch) for ch in generator.all_patterns()]
    DB_seq=[s.seq for s in sequences]
    for i in DB_ch:
        print(i)
    #DB_ch=sort_db_ch(DB_ch)    
    print("======== make_noise in CH DB =======")
    #p,c=make_noise(DB_c,0.3,200)
    c,p =split_db_ch(DB_c,per=.3)

    print("======== Chronicle p generation =======")
    p_seq=[s.seq for s in generator.generate(nb=6, l=6, npat=3, th=.1,patterns=p,pert=True)]

    print("======== SEQUENCES =======")
    print("======== label seq ==========")
    seq_df=pd.DataFrame({'sequence':[], 'label':[]})
    seq_df.sequence=pd.Series(DB_seq+p_seq)
    seq_df.label=pd.Series(len(DB_seq)*[1]+len(p_seq)*[0])
    print(len(p))
    print(len(p_seq))
    print(len(DB_seq))
    print(seq_df)
    print("======== Searching for a sequence in the chronicle database I =======")
    Bseq=(seq_df.sequence).tolist()
    predict=predict(DB_ch,Bseq)
    print(predict)
    y_pred=decision(predict,1)
    print(y_pred)
    """
