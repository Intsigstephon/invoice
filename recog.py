# -*- coding: utf-8 -*-
from keras.models import Model
import numpy as np
from PIL import Image
import keras.backend  as K
import os
from keras.models import load_model
from charset import charset
from keras.backend.tensorflow_backend import set_session
import tensorflow as tf
import cv2
import time
from glob import glob

def decode(pred):
    batch_size = pred.shape[0]   
    length = pred.shape[1]    
    t = pred.argmax(axis=2)   
    char_lists = []
    n = len(charset)

    for i in range(batch_size):
        char_list = ''
        for ii in range(length):
            c = t[i]   
            if c[ii] != n and (not (ii > 0 and c[ii - 1] == c[ii])):
               char_list = char_list + charset[c[ii]]  
        char_lists.append(char_list)

    return char_lists

def init(gpu):
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu)
    config = tf.ConfigProto()
    config.gpu_options.per_process_gpu_memory_fraction = 0.06
    set_session(tf.Session(config=config))
    modelPath = "./model/basemodel509.h5"
    basemodel = load_model(modelPath)
    return basemodel

def predict(im, basemodel):
    out = basemodel.predict(im)
    y_pred = out[:,2:,:]   

    numlabel = [26,93,25,94,632,631,933,29,27,1109,5530,5531]
    a = np.zeros((5532))  
    a[numlabel] = 1
    y_pred[0] = np.multiply(a, y_pred[0])   
    y_pred[1] = np.multiply(a, y_pred[1])

    #---------------------------------------
    out = decode(y_pred)          
    return out

def date_postprocess(tmp):
    localtime = time.localtime(time.time())
    def check_valid(tmp):
        if tmp[0:4].isdigit() and tmp[5:7].isdigit() and tmp[8:10].isdigit() and tmp[4] == u'年' and tmp[7] == u'月' and tmp[10] == u'日':
            return 1
        return 0
    if check_valid(tmp) == 1:
          tmp = tmp.replace(u"年","-")
          tmp = tmp.replace(u"月","-")
          tmp = tmp.replace(u"日","")
          return tmp
    tt = tmp
    tt = tt.replace(u" ","/")
    tt = tt.replace(u"年","/")
    tt = tt.replace(u"月","/")
    tt = tt.replace(u"日","/")
    tt = tt.split("/")
    try:
        if int(tt[0])> 2030 or int(tt[0]) < 1990:
            tt[0] = localtime.tm_year
        if int(tt[1])> 12 or int(tt[0]) == 0:
            if localtime.tm_mon == 1:
                tt[1] = 12
            else:
                tt[1] = localtime.tm_mon - 1
        if len(tt[2]) > 2:
            tt[2] = tt[2][0:2]
        tt = [str(i) for i in tt]
        if len(tt[1]) == 1:
            tt[1] = "0" + tt[1]
        tmp = tt[0]+"-"+tt[1]+"-"+tt[2]
    except Exception:
        tmp = str(localtime.tm_year) + "-" + str(localtime.tm_mon-1) + "-" + str(localtime.tm_mday)

    return tmp 

def money_postprocess(tmp):
    if len(tmp) == 0:
       return 0.00
    if tmp[-1] == ".":
        tmp = tmp[0:-1]
    index = tmp.find('.')
    if(index == -1):  #not found; add .
        if(len(tmp) > 2):
            tmpl = list(tmp)
            tmpl.insert(len(tmpl) - 2, '.')
            tmp = "".join(tmpl)
    elif len(tmp) - 1 - index < 2: #add 0
        tmp = tmp + '0'
    elif len(tmp) - 1 - index > 2: #erase 
        tmp = tmp[0:index+3]
    else:
        pass
    index = tmp.find(u"¥")
    if not(index == -1):
       tmp = tmp[index+1:]
    return tmp

def forward(imgs, model):    
    Max = 0
    for im in imgs:
        if im is None:
            return

        h = im.shape[0]
        w = im.shape[1]
        w = int(float(w)/(float(h)/32.0))
        if w > Max:
           Max = w
    
    im_group = None
    for im in imgs:
        h = im.shape[0]
        w = im.shape[1]
        w = int(float(w)/(float(h) / 32.0))
        im = cv2.resize(im,(w,32))
        im = cv2.copyMakeBorder(im,0,0,0,Max-w,cv2.BORDER_CONSTANT,value=250)
        im = im.astype(np.float32)
        im = ((im/255.0)-0.5)*2
        X  = im.reshape((32, Max, 1))
        X = np.array([X])
        if im_group is None:
           im_group = X
           continue

        im_group = np.concatenate((im_group, X), axis=0) 

    #ocr
    result = predict(im_group, model)       #handle parallel 

    return result

if __name__ == '__main__':

    #load basemodel
    basemodel = init()   #last 0.9s   #no training configure saved in *.h5

    #prepare data
    imagedir = "/home/stephon/下载/Invoice_Image/img/4"
    imgs = glob(imagedir + "/*.jpeg")
    imgs.extend(glob(imagedir + "/*.jpg"))

    batch_size = 60
    imgs = imgs[0: batch_size]
    print(len(imgs))
    for img in imgs:
        print(img)

    #resize and concatenate
    aa = None
    Max = 0

    #get the max len after resize
    for i in range(batch_size):
        path_1 = imgs[i]
        im = cv2.imread(path_1)
        if im is None:
            exit()

        h = im.shape[0]
        w = im.shape[1]
        w = int(float(w)/(float(h)/32.0))
        if w > Max:
           Max = w
    
    print("Max length is %d" % Max)
    
    #concatenate data
    for i in range(batch_size):
        path_1 = imgs[i]
        im = cv2.imread(path_1)
        h = im.shape[0]
        w = im.shape[1]
        w = int(float(w)/(float(h)/32.0))
        im = cv2.resize(im, (w,32))
        im = cv2.cvtColor(im,cv2.COLOR_BGR2GRAY)
        im = cv2.copyMakeBorder(im,0,0,0,Max-w,cv2.BORDER_CONSTANT,value=250)
        im = im.astype(np.float32)
        im = ((im/255.0)-0.5)*2   #-1~1
        X  = im.reshape((32, Max, 1))   #
        X = np.array([X])
        if aa is None:
           aa = X
           continue

        aa = np.concatenate((aa, X), axis=0)   #last 0.007s

    #do ocr
    start = time.time()
    out = basemodel.predict(aa)  #1: 0.8s   #10: 0.84s  #60: 1.29s; need to parallel
    end = time.time()

    numlabel = [26,93,25,94,632,631,933,29,27,1109,5201,5530,5531]
    numlabel.extend([1, 466])       
    a = np.zeros((5532))
    a[numlabel] = 1
    for i in range(len(out)):
        out[i] = np.multiply(a, out[i])
    # print(out.shape)
    y_pred = out[:,2:,:]    
    result = decode(y_pred)       #handle parallel   #last 0.004s
    for ii in result:
        print ii
    #time cost
    print("time cost:" + str(end - start))
    print("-------------------")



