
from __future__ import division

import copy
import csv
import time
from collections import Counter
from datetime import datetime, timedelta
from math import sqrt

import distance
import matplotlib.pyplot as plt
import numpy as np
import unicodecsv
from jellyfish._jellyfish import damerau_levenshtein_distance
from keras.models import load_model
from sklearn import metrics
import functools 
from math import exp
import sys
import pickle
import category_encoders as ce
import pandas as pd

def loadobj2(filename):
    with open("./LSTM/output_files/data/%s"% filename, 'rb') as config_file:
        return pickle.load(config_file)

def anomalydetect(DB_seq,param,mod):


    def extract_feature(DB_seq):
        lines = [] #these are all the activity seq
        timeseqs = [] #time sequences (differences between two events)
        timeseqs2 = [] #time sequences (differences between the current and first)
        caseids = []
        codelines=[]

        times = []
        times2 = []
        evnts=[]

        numlines= len(DB_seq)
        for i,seq in enumerate(DB_seq): #the rows are "ChID,sequence,TC"
            if len(seq)==0:
                continue
            caseids.append(i)
            lastevnettime=seq[0][0]
            firsteventtime=seq[0][0]
			
            
 
            for t,e in seq:
                evnts.append(e)
                times.append(t-lastevnettime)
                times2.append(t-firsteventtime)
                lastevnettime=t
            b=ce_bin.transform(pd.DataFrame(evnts,columns=['evt']))
            codelines.append(b.values.tolist())
            lines.append(evnts)
            timeseqs.append(times)
            timeseqs2.append(times2)
            return  lines,timeseqs,timeseqs2,caseids,codelines




    #######################################################################
    divisor,divisor2,maxlen,target_chars,target_char_indices,target_indices_char,char_indices,chars=param.divisor,param.divisor2,param.maxlen,param.target_chars,param.target_char_indices,param.target_indices_char,param.char_indices,param.chars
    label=["evt"]
    df = pd.DataFrame(chars,columns=label)
    ce_bin = ce.BinaryEncoder(cols=['evt'],drop_invariant=True)
    r=ce_bin.fit_transform(df.evt.to_frame())
    lines,timeseqs,timeseqs2,caseids,codelines=extract_feature(DB_seq)




    """
    X=loadobj2("X.dictionary")
    y_a=loadobj2("y_a.dictionary")
    y_t=loadobj2("y_t.dictionary")
    """


    # load model, set this to the model generated by train.py
    #model = load_model('./output_files/model_94-0.83.h5')
    model =mod


    ######################################################################
    def feature_modilisation(sentences,sentences_t,codelines,maxlen=maxlen):
        num_features =  r.shape[1]+2
        #print('num features: {}'.format(num_features))
        X = np.zeros((len(sentences), maxlen, num_features), dtype=np.float32)
        for i, sentence in enumerate(sentences):
            leftpad = maxlen-len(sentence)
            sentence_t = sentences_t[i]
            bb=codelines[i]

            for t, char in enumerate(sentence):
                #multiset_abstraction = Counter(sentence[:t+1])
                X[i, t+leftpad ]=bb[t]+[t+1 , sentence_t[t]/divisor] 
        return X
    # define helper functions
    def encode(sentence, times,code, maxlen=maxlen):
        num_features = r.shape[1]+2
        X = np.zeros((1, maxlen, num_features), dtype=np.float32)
        leftpad = maxlen-len(sentence)
        #times2 = np.cumsum(times)
        for t, char in enumerate(sentence):
            #multiset_abstraction = Counter(sentence[:t+1])
            X[i, t+leftpad ]=code[t]+[t+1 , times[t]/divisor] 
        return X

    def getSymbol(predictions):
        maxPrediction = 0
        symbol = ''
        i=np.where(predictions == np.amax(predictions))[0]
        symbol = target_indices_char[i[0]]
        return symbol
    def prob_likelihood(t1,t2):
        lamda=0.0000015
        #lamda=0.0045
        p=exp(-lamda*(abs(t1-t2)))
        return p
    
    #start_time0 = time.perf_counter()
    Sim=[] #event matching
    Sim_t=[] #time matching
    predicted = []

    #set parameters
    predict_size = 1

    # make predictions
    line=lines[0]
    caseid=caseids[0]
    times=timeseqs[0]
    codeline=codelines[0]
    sentences=[]
    sentence_t=[]
    codes=[]
    ground_truth=[]
    ground_truth_t=[]
    for prefix_size in range(2,len(line)):
        #print(prefix_size)
        cropped_line = line[:prefix_size]
        cropped_times = times[:prefix_size]
        code=codeline[:prefix_size]
        #### truth
        ground_truth.append(line[prefix_size:prefix_size+predict_size])
        ground_truth_t.append(times[prefix_size:prefix_size+predict_size])
        #######
        sentences.append(cropped_line)
        sentence_t.append(cropped_times)
        codes.append(code)
        #start_time = time.perf_counter()
        enc=feature_modilisation(sentences,sentence_t,codes)
        end1 = time.perf_counter()
       ############## 
    predicted = []
    predicted_t = []
    for i in range(predict_size):
        if len(ground_truth)<=i:
            continue
            ############
        y = model.predict(enc, verbose=0)
        #end2= time.perf_counter()
        #print("encode",end1-start_time,"model",end2-end1)
        y_char = y[0]
        y_t = y[1]
        for char in y_char:
            prediction = getSymbol(char)              
            predicted.append(prediction)
  
        predicted_t =[max(t[0],0)* divisor for t in y_t ]
        #print(predicted_t,times[2:])
    if len(line)>0:
        Sim = 1 - distance.nlevenshtein(str(predicted), str(line[2:]))
        #print(str(predicted), str(line[2:]),Sim)
        #for t,d in zip(times[2:],predicted_t):
        #    print(t,d)
        Sim_t= [ prob_likelihood(t,d)   for t,d in zip(times[2:],predicted_t) ]
        #print(Sim_t)
    #start_time2 = time.perf_counter()
    #print("end seq",start_time2-start_time0)
    if  Sim>.8:
        return functools.reduce (lambda a, b: a * b ,Sim_t,1)
    else:
        return 0

