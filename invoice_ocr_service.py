#!/usr/bin/env python  
#encoding: utf-8  
#Author:lubin
#Date:2018.04.11
#tips:
#if the socket was locked, then do the below command to release the socked file:
#rm -rf /tmp/service_name.sock

import atexit, os, sys, time, signal
from deamon import cdaemon
import socket
import json
import base64
import cv2
import redis
import cPickle as pickle
from StringIO import StringIO
from PIL import Image
import struct
import numpy as np

service_name ='invoice_ocr_service'
class ClientDaemon(cdaemon):
    def __init__(self, name, gpu, save_path, stdin=os.devnull, stdout=os.devnull, stderr=os.devnull, home_dir='.', umask=022, verbose=1,sock=None,port=None):  
        cdaemon.__init__(self, save_path, stdin, stdout, stderr, home_dir, umask, verbose)  
        self.name = name 
        self.sock = sock
        self.port = port
        self.gpu = gpu

    def readb64(self, base64_string):
        buf = base64.b64decode(base64_string)
        buf = np.asarray(bytearray(buf), dtype="uint8")
        image = cv2.imdecode(buf, cv2.IMREAD_GRAYSCALE)
        return image

    def parse(self, data):
        app_data=json.loads(data["data"])
        request_id = app_data["request_id"]
        image_io = app_data["img"]

        image = self.readb64(image_io)
        return request_id, image
 
    def recv(self, connection):
        """
        two stage receive
        """
        data = connection.recv(4)
        fileLen = struct.unpack('i', data)[0]      
        data=''
        while fileLen > 0:
            readLen = 1024
            tmp = connection.recv(readLen)
            data += tmp
            #print("********************************{}".format(len(tmp)))
            fileLen -= len(tmp)
        return data

    def run(self, output_fn, **kwargs):  
        import invoice_ocr
        self.lib, self.basemodel= invoice_ocr.init_model(self.gpu)
        print("init done")
        server = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        sock="/tmp/%s"%(self.sock)
        if os.path.exists(sock):
            os.unlink(sock)
        server.bind(sock)
        print(sock)
        server.listen(99)

        while True:
            connection, address = server.accept() 
            try:
                connection.settimeout(20000)
                a={"code":"0", "rslt":''}
                data = self.recv(connection)
                print("receive done")

                data=json.loads(data)
                request_id, img = self.parse(data)
                ret_v = invoice_ocr.invoice_ocr(img, self.lib, self.basemodel)
                a={"code":"0", "rslt":ret_v}
                connection.send(json.dumps(a))
            except socket.timeout:
                a={"code":"1","rslt":''}
                print('socket time out')
                connection.send(json.dumps(a))
            except Exception,e:  
                print("=-------------------->{}".format(str(e)))
                a={"code":"1","rslt":''}
                connection.send(json.dumps(a))
            finally:  
                connection.close()
        self.sess.close()


if __name__ == '__main__':  
    help_msg = 'Usage: python %s <start|stop|restart|status>' % sys.argv[0]  
    if len(sys.argv) != 4:
        print(help_msg)
        sys.exit(1)

    port=sys.argv[2]
    gpu= sys.argv[3]
    tag='{}_{}'.format(service_name, port)
    sock = tag +".sock"

    #r = redis.StrictRedis(host='localhost', port=6379, db=0)
    #IDENTITY_CARD_MOSAIC_ALGO = pickle.loads(r.hget('identity_card_mosaic_algo','socket_port_map'))
    #sock=IDENTITY_CARD_MOSAIC_ALGO[port][service_name]

    pid_fn = '/tmp/%s.pid'%(tag)     #process pid
    log_fn = '/tmp/%s.log'%(tag)     #process log absolute path 
    err_fn = '/tmp/%s.err.log'%(tag) #process error log
    cD = ClientDaemon('%s_%s'%(tag,"deamon"), gpu, pid_fn, stdout=log_fn, stderr=err_fn, verbose=1, sock=sock, port=port)
  
    if sys.argv[1] == 'start':  
        cD.start(log_fn)  
    elif sys.argv[1] == 'stop':  
        cD.stop()  
    elif sys.argv[1] == 'restart':  
        cD.restart(log_fn)  
    elif sys.argv[1] == 'status':  
        alive = cD.is_running()
        if alive:  
            print('process [%s] is running ......' % cD.get_pid())
        else:  
            print('daemon process [%s] stopped' %cD.name)
    else:  
        print('invalid argument!')
        print(help_msg)
