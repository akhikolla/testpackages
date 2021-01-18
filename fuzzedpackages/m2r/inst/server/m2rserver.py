import socket
import sys
import time
import random
import string
import os
import inspect

from subprocess import Popen
from subprocess import call
from contextlib import closing

MAX_PORTS = 20
LISTEN_PORT = 27435
FIRST_PORT = (LISTEN_PORT+1)


def m2PortForRPort(rport):
	return 2*LISTEN_PORT - rport

def generateLogFileName(prefix = "session", N = 10):
	randpart = ''.join(random.SystemRandom().choice(string.ascii_uppercase + string.digits) for _ in range(N))
	
	ret = "m2logs/" + prefix + "-"
	ret = ret + time.strftime("%Y-%m-%d-%H-%M-%S-", time.gmtime())
	ret = ret + randpart + ".txt"
	
	return ret

def pathOfFile():
	return os.path.dirname(inspect.stack()[0][1])


with closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as serversocket:
	serversocket.bind((socket.gethostname(), LISTEN_PORT))
	serversocket.listen(5)
	
	while True:
		# Wait for a connection
		print("Waiting for a connection")
		connection, client_address = serversocket.accept()
		
		try:
			print("Connection from " + str(client_address))
			
			rport = 0
			containername = ""
			
			portrange = range(FIRST_PORT, FIRST_PORT+MAX_PORTS)
			
			# find the first open port
			for port in portrange:
				m2port = m2PortForRPort(port)
				m2port = port
				with closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as sock:
					if sock.connect_ex(("localhost",port)) != 0:
						# create M2 docker instance
						containername = "m2rserver" + str(port)
						retval = call(["docker", "create", "-it", "--expose=" + str(m2port), "--publish=" + str(m2port) + ":" + str(m2port), "--name=" + containername, "sommars/m2r"])
						if retval == 0:
							rport = port
							break
			
			if rport != 0:
				logfilename = generateLogFileName()
				
				# launch M2 docker instance
				# /bin/sh -c docker start m2rserver27436 &>/dev/null ; 
				#   docker exec m2rserver27436 M2 --script /m2rserverscript.m2 27434 &>m2logs/session-2017-07-27-14-20-31-LC8QEU72VQ.txt ; 
				#   docker stop m2rserver27436 &>/dev/null ; 
				#   docker rm m2rserver27436 &>/dev/null ; 
				#   aws s3 cp m2logs/session-2017-07-27-14-20-31-LC8QEU72VQ.txt s3://m2r-dev-logs/logfiles/
				cmd = ["docker", "start", containername, "&>/dev/null"]
				cmd = cmd + [";", "docker", "exec", containername, "M2", "--script", "/m2rserverscript.m2", str(m2port), "&>" + logfilename]
				cmd = cmd + [";", "docker", "stop", containername, "&>/dev/null"]
				cmd = cmd + [";", "docker", "rm", containername, "&>/dev/null"]
				cmd = cmd + [";", "aws", "s3", "cp", logfilename, "s3://m2r-dev-logs/logfiles/"]
				Popen(" ".join(cmd), shell=True)
				# p.terminate()
				
				# # launch python middleman
				# cmd = ["python", "/home/ec2-user/m2rserver/m2rdockerserver.py", str(rport), str(m2port), "&>dockerserverlog" + str(rport) + ".txt"]# "&>/dev/null"]
				# Popen(" ".join(cmd), shell=True)
				
				print("Session logging to " + logfilename)
			else:
				print("No session started, all ports taken")
			
			# send open port back to the client
			connection.sendall(str(rport) + "\n")
		
		finally:
			# Clean up the connection
			connection.close()
		
		print("")





