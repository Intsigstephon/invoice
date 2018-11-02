#coding=utf-8

import cv2
import ctypes
import numpy as np
import time
import os
from glob import glob

width  = 2592
height = 1656

lu = [107,  356]
ru = [2496, 353]
ld = [108,  1486] 
rd = [2496, 1481] 

def loadlib():
    so = ctypes.cdll.LoadLibrary
    lib = so("./libLine.so")
    return lib

def vertex(img, lib):
    """
    img(one channel, gray) and lib
    output vertex: left-up, right-up, right-down, left-down
    """
    pyvertex = [0] * 8
    cvertex = (ctypes.c_int * len(pyvertex))(*pyvertex)
    flag = lib.getVertex(img.ctypes.data_as(ctypes.POINTER(ctypes.c_ubyte)), img.shape[1], img.shape[0], img.shape[1], cvertex)
    rslt = np.array(list(cvertex))
    return rslt, flag

def detect(img, lib, Debug=False):
    """
    img and lib
    output: perspective image and rect list
    """
    boxes = []
    rslt, flag = vertex(img, lib)
    if flag < 0:
        return img, boxes, flag

    #do perpective image
    pts1 = np.float32([[rslt[0], rslt[1]], [rslt[2], rslt[3]], [rslt[6], rslt[7]], [rslt[4], rslt[5]]])
    pts2 = np.float32([lu, ru, ld, rd])

    M =   cv2.getPerspectiveTransform(pts1, pts2)
    dst = cv2.warpPerspective(img, M, (width, height))

    #area need recognize
    #solid area
    box = [374, 123, 805, 199]
    boxes.append(box)
    box = [1888, 112, 2270, 199]
    boxes.append(box)

    if Debug:
        cv2.imshow("dst", dst)
        cv2.waitKey(0)

    return dst, boxes, flag

def test_vertex(imgdir):
    """
    input a imgdir, test all the images 
    """
    testdir = imgdir + "/test" 
    if os.path.exists(testdir):
        for file in os.listdir(testdir):
            os.remove(os.path.join(testdir, file))
        os.rmdir(testdir)

    imgs = glob(imgdir + "/*.jpeg")
    imgs.extend(glob(imgdir + "/*.jpg"))    

    os.mkdir(testdir)
    lib = loadlib()

    for j in range(len(imgs)):
        imgpath = os.path.join(imgdir, imgs[j])
        print("{}: {}".format(j, imgpath))
        img = cv2.imread(imgpath, 0)
        rslt, flag = vertex(img, lib)                    #get rslt and flag
        print("result is {}: {}".format(flag, rslt))

        if(flag < 0):
            continue

        #show border line
        img_rgb = cv2.imread(imgpath, 1)
        for i in range(3):
            cv2.circle(img_rgb, (rslt[2 * i + 0], rslt[2 * i + 1]), 15, (255,0,0), -1)
            cv2.line(img_rgb, (rslt[2 * i + 0], rslt[2 * i + 1]), (rslt[2 * i + 2], rslt[2 * i + 3]), (0,0,255), 3)

        cv2.circle(img_rgb, (rslt[6], rslt[7]), 10, (255,0,0), -1)
        cv2.line(img_rgb, (rslt[6], rslt[7]), (rslt[0], rslt[1]), (0,0,255), 3)   

        # cv2.imshow("rslt", img_rgb)
        # cv2.waitKey(0)

        cv2.imwrite(os.path.join(testdir, imgs[j]), img_rgb)

def test_detect(imgdir):
    """
    input a imgdir, test all the images 
    """
    testdir = imgdir + "/test" 
    if os.path.exists(testdir):
        for file in os.listdir(testdir):
            os.remove(os.path.join(testdir, file))
        os.rmdir(testdir)

    imgs = glob(imgdir + "/*.jpeg")
    imgs.extend(glob(imgdir + "/*.jpg"))

    os.mkdir(testdir)
    lib = loadlib()

    for i in range(len(imgs)):
        imgpath = os.path.join(imgdir, imgs[i])
        print("{}: {}".format(i, imgpath))
        img = cv2.imread(imgpath, 0)
        rslt = detect(img, lib)
        if rslt[2] < 0:
            print("detect error\n")
            continue

        img_dst = rslt[0]
        boxes   = rslt[1]
        print("succeed! result is {}".format(boxes))

        #show rect
        for tmp in boxes:
            cv2.rectangle(img_dst, (tmp[0], tmp[1]), (tmp[2], tmp[3]), (0,0,255), 2)

        # cv2.namedWindow("rslt img", 0)
        # cv2.imshow("rslt img", dst)
        # cv2.waitKey(0)

        cv2.imwrite(os.path.join(testdir, imgs[i]), img_dst)

if __name__ == "__main__":
    
    #debug for one
    # imgpath = "./test/bad0918/S28BW-418090316460_0006.jpg"
    # lib = loadlib()
    # img = cv2.imread(imgpath, 0)
    # detect(img, lib, 1)

    imgdir = "/media/stephon/_E/OCR_Invoice/Invoice_Image/img"
    #test_vertex(imgdir)
    test_detect(imgdir)
    #test_crop(imgdir, 3)

    # start = time.time()
    # end = time.time()
    # print("%f.s" % (end - start))