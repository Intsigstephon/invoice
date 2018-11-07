#coding=utf-8

import cv2
import ctypes
import numpy as np
import time
from recog import forward, init
from preprocess import preprocess
from glob import glob
import os
from detect import loadlib, detect

debug = 1
def crop(img, boxes):
    imgs = []
    for i in range(len(boxes)):
        box = boxes[i]
        imgs.append(img[box[1]:box[3], box[0]:box[2]])
    return imgs

def init_model(gpu):
    lib = loadlib()
    basemodel = init(gpu)
    return lib, basemodel

def invoice_ocr(img, lib, basemodel):
    #1. detect vertex
    rslt = []

    dst, boxes, flag = detect(img, lib)

    if flag < 0:
        return rslt

    if debug:
        img_show = dst
        for tmp in boxes:
            cv2.rectangle(img_show, (tmp[0], tmp[1]), (tmp[2], tmp[3]), (0,0,255), 2)
        cv2.imwrite("./tmp/a.jpg", img_show)

    #2.crop image:
    imgs = crop(dst, boxes)
  
    #3.do preprocess
    imgs = preprocess(imgs)

    if debug:
        for i in range(len(boxes)):
            cv2.imwrite("./tmp/" + str(i) + ".jpg", imgs[i])
    
    #4.recog
    rslt = forward(imgs, basemodel)

    return rslt

if __name__ == '__main__':

    #imgdir = "/media/stephon/_E/OCR_Invoice/Invoice_Image/img"
    imgdir = "./test"
    paths = glob(imgdir + "/*.jpeg")
    paths.extend(glob(imgdir + "/*.jpg"))
    lib, basemodel = init_model(0)

    for path in paths:
        print(path)
        start = time.time() 
        img = cv2.imread(path, 0)
        rslt = invoice_ocr(img, lib, basemodel) 
        end = time.time()
        if len(rslt) < 1:
            print("bad case")

        for x in rslt:
            print(x)
            
        print("time cost:" + str(end - start))

        #show img
        img = cv2.imread(path, 1)
        # cv2.imshow("src", img)
        # cv2.waitKey(0)


    print("-------------------\n")
