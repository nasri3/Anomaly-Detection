from keras.models import Sequential, Model
from keras.layers.core import Dense
from keras.layers.recurrent import LSTM, GRU, SimpleRNN
from keras.layers import Input
from keras.utils.data_utils import get_file
from keras.optimizers import Nadam
from keras.callbacks import EarlyStopping, ModelCheckpoint, ReduceLROnPlateau
from keras.layers.normalization import BatchNormalization
from collections import Counter
import unicodecsv
import numpy as np
import random
import sys
import os
import copy
import csv
import time
import pickle
import category_encoders as ce
import pandas as pd

from datetime import datetime
from math import log

from keras import backend as K
print((K.tensorflow_backend._get_available_gpus()))

import keras
import tensorflow as tf
config = tf.ConfigProto(device_count={'GPU':1})
sess = tf.Session(config = config)
K.set_session(sess)

if tf.test.gpu_device_name():
    print('GPU found')
else:
    print("No GPU found")



if __name__ == "__main__":
    from LSTMpredict import *
    from LSTMmonotoring import *
else:
    from .LSTMpredict import *
    from .LSTMmonotoring import *
import re 


class modelparameter:
    def __init__(self,divisor,divisor2,maxlen,target_chars,target_char_indices,target_indices_char,char_indices,chars):
        self.divisor,self.divisor2,self.maxlen,self.target_chars,self.target_char_indices,self.target_indices_char,self.char_indices,self.chars=divisor,divisor2,maxlen,target_chars,target_char_indices,target_indices_char,char_indices,chars
    def __str__(self):
        return '{0}{1}{2}{3}{4}{5}{6}{7}'.format(self.divisor,self.divisor2,self.maxlen,self.target_chars,self.target_char_indices,self.target_indices_char,self.char_indices,self.chars)
        
  
def loadobj1(filename):
    with open("./LSTM/output_files/data/%s"% filename, 'rb') as config_file:
        return pickle.load(config_file)

def extract(input):  
     # \d+ is a regular expression which means 
     # one or more digit 
     # output will be like ['100','564','365'] 
    numbers = re.findall('\d+',input) 

     # now we need to convert each number into integer 
     # int(string) converts string into integer 
    numbers = list(map(int,numbers))
    l=[]
    for i in range(0,len(numbers)-1,2):
        l.append((numbers[i],numbers[i+1]))
    return [ s for s in l if len(s)!=0]

def read_text_file(filename):
    print('Reading file ' + filename + "...")
    with open(filename, "r", encoding='utf8') as textfile:
        L = []
        for line in textfile:
            L.append(line.strip())
        print('File contains ', len(L), "lines.\n")
        return L

#lignes =seq ch+!
#timeseqs= TC0
#timeseqs2=TC1
def train(DB_seq):

    lines = [] #these are all the activity seq
    timeseqs = [] #time sequences (differences between two events)
    timeseqs2 = [] #time sequences (differences between the current and first)
        
    times = []
    times2 = []
    evnts=[]

    numlines= len(DB_seq)
    for seq in DB_seq: #the rows are "ChID,sequence,TC"
        lastevnettime=seq[0][0]
        firsteventtime=seq[0][0]

        times = []
        times2 = []
        evnts=[] 
        for t,e in seq:
            evnts.append(e)
            times.append(t-lastevnettime)
            times2.append(t-firsteventtime)
            lastevnettime=t
        lines.append(evnts)
        timeseqs.append(times)
        timeseqs2.append(times2)
        

    ########################################

    divisor = np.mean([item for sublist in timeseqs for item in sublist]) #average time between events
    print('divisor: {}'.format(divisor))
    divisor2 = np.mean([item for sublist in timeseqs2 for item in sublist]) #average time between current and first events
    print('divisor2: {}'.format(divisor2))



    #########################################################################################################


    step = 1
    sentences = []
    softness = 0
    next_chars = []
    #lines = [x+'!' for x in lines]
    maxlen = max([len(x) for x in lines]) #find maximum line size

    # next lines here to get all possible characters for events and annotate them with numbers
    chars = [set(x) for x in lines]
    chars = list(set().union(*chars))
    chars.sort()
    target_chars = copy.copy(chars)
    #chars.remove('!')
    print('total chars: {}, target chars: {}'.format(len(chars), len(target_chars)))
    char_indices = dict((c, i) for i, c in enumerate(chars))
    indices_char = dict((i, c) for i, c in enumerate(chars))
    target_char_indices = dict((c, i) for i, c in enumerate(target_chars))
    target_indices_char = dict((i, c) for i, c in enumerate(target_chars))
    print(indices_char)

    ######### save parameters ###########
    print(divisor,divisor2,maxlen,target_chars,target_char_indices,target_indices_char,char_indices,chars)
    m=modelparameter(divisor,divisor2,maxlen,target_chars,target_char_indices,target_indices_char,char_indices,chars)
    with open('./LSTM/output_files/config.dictionary', 'wb') as config_file:
        pickle.dump(m, config_file)
    
    #####################################
    
    label=["evt"]
    df = pd.DataFrame(chars,columns=label)
    ce_bin = ce.BinaryEncoder(cols=['evt'],drop_invariant=True)
    r=ce_bin.fit_transform(df.evt.to_frame())
    ##############################
    codelines=[]
    sentences_t = []
    next_chars_t = []
    lines_t=timeseqs
    for line, line_t in  zip(lines, lines_t):
        b=ce_bin.transform(pd.DataFrame(line,columns=['evt']))
        bb=b.values.tolist()
        for i in range(0, len(line), step):
            if i==0:
                continue

            #we add iteratively, first symbol of the line, then two first, three...
            codelines.append(bb[0:i])
            sentences.append(line[0: i])
            sentences_t.append(line_t[0:i])
            next_chars.append(line[i])
            if i==len(line)-1: # special case to deal time of end character
                next_chars_t.append(0)
            else:
                next_chars_t.append(line_t[i])

    print('nb sequences:', len(sentences))

    print('Vectorization...')
    num_features = r.shape[1]+2
    print('num features: {}'.format(num_features))
    
    X = np.zeros((len(sentences), maxlen, num_features), dtype=np.float32)
    y_a = np.zeros((len(sentences), len(target_chars)), dtype=np.float32)
    y_t = np.zeros((len(sentences)), dtype=np.float32)
    for i, sentence in enumerate(sentences):
        leftpad = maxlen-len(sentence)
        next_t = next_chars_t[i]

        sentence_t = sentences_t[i]
        bb=codelines[i]
        #b=ce_bin.transform(pd.DataFrame(sentence,columns=['evt']))
        #bb=b.values.tolist()
        for t, char in enumerate(sentence):
            #multiset_abstraction = Counter(sentence[:t+1])
            X[i, t+leftpad ]=bb[t]+[t+1 , sentence_t[t]/divisor] 

        y_a[i, target_char_indices[next_chars[i]]] = 1-softness
        y_t[i] = next_t/divisor
        np.set_printoptions(threshold=sys.maxsize)
    """
    with open('./LSTM/output_files/data/X.dictionary', 'wb') as config_file:
        pickle.dump(X, config_file)
    
    with open('./LSTM/output_files/data/y_a.dictionary', 'wb') as config_file:
        pickle.dump(y_a, config_file)
    
    with open('./LSTM/output_files/data/y_t.dictionary', 'wb') as config_file:
        pickle.dump(y_t, config_file)
    
    X=loadobj1("X.dictionary")
    y_a=loadobj1("y_a.dictionary")
    y_t=loadobj1("y_t.dictionary")
    """
    # build the model: 
    print('Build model...')
    print(maxlen)
    #keras.backend.get_session().run(tf.global_variables_initializer())
    main_input = Input(shape=(maxlen, num_features), name='main_input')
    # train a 2-layer LSTM with one shared layer
    l1 = LSTM(100, implementation=2, kernel_initializer='glorot_uniform', return_sequences=True, dropout=0.2)(main_input) # the shared layer
    b1 = BatchNormalization()(l1)
    l2 = LSTM(100, implementation=2, kernel_initializer='glorot_uniform', return_sequences=True, dropout=0.2)(b1) # the shared layer
    b2 = BatchNormalization()(l2)
    l3 = LSTM(100, implementation=2, kernel_initializer='glorot_uniform', return_sequences=True, dropout=0.2)(b2) # the shared layer
    b3 = BatchNormalization()(l3)
    l2_1 = LSTM(100, implementation=2, kernel_initializer='glorot_uniform', return_sequences=False, dropout=0.2)(b3) # the layer specialized in activity prediction
    b2_1 = BatchNormalization()(l2_1)
    l2_2 = LSTM(100, implementation=2, kernel_initializer='glorot_uniform', return_sequences=False, dropout=0.2)(b3) # the layer specialized in time prediction
    b2_2 = BatchNormalization()(l2_2)

    act_output = Dense(len(target_chars), activation='softmax', kernel_initializer='glorot_uniform', name='act_output')(b2_1)
    time_output = Dense(1, kernel_initializer='glorot_uniform', name='time_output')(b2_2)

    model = Model(inputs=[main_input], outputs=[act_output, time_output])

    opt = Nadam(lr=0.002, beta_1=0.9, beta_2=0.999, epsilon=1e-08, schedule_decay=0.004, clipvalue=3)

    model.compile(loss={'act_output':'categorical_crossentropy', 'time_output':'mae'}, optimizer=opt)
    early_stopping = EarlyStopping(monitor='val_loss', patience=42)
    model_checkpoint = ModelCheckpoint('./LSTM/output_files/model.h5', monitor='val_loss', verbose=0, save_best_only=True, save_weights_only=False, mode='auto')
    lr_reducer = ReduceLROnPlateau(monitor='val_loss', factor=0.5, patience=10, verbose=0, mode='auto', min_delta=0.0001, cooldown=0, min_lr=0)

    h=model.fit(X, {'act_output':y_a, 'time_output':y_t}, validation_split=0.2, verbose=2, callbacks=[early_stopping, model_checkpoint, lr_reducer], batch_size=maxlen, epochs=500)

    return h

if __name__ == "__main__":
    """
    p_seq=[[(115, 3), (116, 19), (446, 16), (538, 5), (669, 1), (735, 5), (779, 13), (801, 14), (819, 0)],
    [(148, 3), (148, 8), (162, 12), (169, 19), (316, 6), (655, 12), (895, 6), (984, 2)],
    [(63, 15), (135, 3), (313, 3), (371, 19), (434, 14), (840, 7), (882, 11), (902, 3)],
    [(84, 18), (109, 12), (243, 7), (381, 11), (396, 3), (465, 19), (484, 2), (530, 2), (673, 6), (784, 15), (992, 5)],
    [(59, 9), (74, 7), (410, 3), (459, 19), (510, 11), (526, 9), (576, 11), (755, 19), (906, 12)]]

    #divisor,divisor2,maxlen,target_chars,target_char_indices,target_indices_char,char_indices,chars=train(p_seq)
    #print(divisor,divisor2,maxlen,target_chars,target_char_indices,target_indices_char,char_indices,chars)

    divisor=91.86666666666666 
    divisor2=392.73333333333335
    maxlen=45
    target_chars=[0, 1, 2, 3, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 18, 19]
    target_char_indices = {0: 0, 1: 1, 2: 2, 3: 3, 5: 4, 6: 5, 7: 6, 8: 7, 9: 8, 11: 9, 12:10, 13: 11, 14: 12, 15: 13, 16: 14, 18: 15, 19: 16} 
    target_indices_char={0: 0, 1: 1, 2: 2, 3: 3, 4: 5, 5: 6, 6: 7, 7: 8, 8: 9, 9: 11, 10: 12, 11: 13, 12: 14, 13: 15, 14: 16, 15: 18, 16: 19} 
    char_indices={0: 0, 1: 1, 2: 2, 3: 3, 5: 4, 6: 5, 7: 6, 8: 7, 9: 8, 11: 9, 12: 10, 13: 11, 14: 12, 15: 13, 16: 14, 18: 15, 19: 16}
    chars= [0, 1, 2, 3, 5, 6, 7, 8, 9, 11, 12, 13, 14,15, 16, 18, 19]
    

    #predict(p_seq,divisor,divisor2,maxlen,target_chars,target_char_indices,target_indices_char,char_indices,chars)
    SS=[anomalydetect([seq],divisor,divisor2,maxlen,target_chars,target_char_indices,target_indices_char,char_indices,chars) for seq in p_seq]
    print(SS)
    
    DB=read_text_file("./LSTM/input_data/data.txt")
    trainseq=[ extract(line) for line in DB ]
    divisor,divisor2,maxlen,target_chars,target_char_indices,target_indices_char,char_indices,chars=train(trainseq)
    #divisor,divisor2,maxlen,target_chars,target_char_indices,target_indices_char,char_indices,chars=104.10958904109589, 357.027397260274, 73, [0, 2, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], {0: 0, 2: 1, 5: 2, 6: 3, 7: 4, 8: 5, 10: 6, 11: 7, 12: 8, 13: 9, 14: 10, 15: 11, 16: 12, 17: 13, 18: 14, 19: 15}, {0: 0, 1: 2, 2: 5, 3: 6, 4: 7, 5: 8, 6: 10, 7: 11, 8: 12, 9: 13, 10: 14, 11: 15, 12: 16, 13: 17, 14: 18, 15: 19}, {0: 0, 2: 1, 5: 2, 6: 3, 7: 4, 8: 5, 10: 6, 11: 7, 12: 8, 13: 9, 14: 10, 15: 11, 16: 12, 17: 13, 18: 14, 19: 15}, [0, 2, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
    print(divisor,divisor2,maxlen,target_chars,target_char_indices,target_indices_char,char_indices,chars)
    SS=[anomalydetect([seq],divisor,divisor2,maxlen,target_chars,target_char_indices,target_indices_char,char_indices,chars) for seq in trainseq]
    print(SS)





    l1 = LSTM(200, implementation=2, kernel_initializer='glorot_uniform', return_sequences=True, dropout=0.2)(main_input) # the shared layer
    b1 = BatchNormalization()(l1)
    l2 = LSTM(200, implementation=2, kernel_initializer='glorot_uniform', return_sequences=True, dropout=0.2)(b1) # the shared layer
    b2 = BatchNormalization()(l2)
    l3 = LSTM(200, implementation=2, kernel_initializer='glorot_uniform', return_sequences=True, dropout=0.2)(b2) # the shared layer
    b3 = BatchNormalization()(l3)
    l2_10 = LSTM(400, implementation=2, kernel_initializer='glorot_uniform', return_sequences=True, dropout=0.2)(b3)
    b2_10 = BatchNormalization()(l2_10)
    l2_20 = LSTM(400, implementation=2, kernel_initializer='glorot_uniform', return_sequences=True, dropout=0.2)(b2_10)
    b2_20 = BatchNormalization()(l2_20)
    l2_1 = LSTM(400, implementation=2, kernel_initializer='glorot_uniform', return_sequences=False, dropout=0.2)(b2_20) # the layer specialized in activity prediction
    b2_1 = BatchNormalization()(l2_1)
    l2_2 = LSTM(200, implementation=2, kernel_initializer='glorot_uniform', return_sequences=False, dropout=0.2)(b3) # the layer specialized in time prediction
    b2_2 = BatchNormalization()(l2_2)
    """
