#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Copyright (C) 2022 Priyanka Golia, Brendan Juba, and Kuldeep Meel

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
'''

from __future__ import print_function
import sys
import os
import math
import random
import argparse
import copy
import tempfile
import numpy as np
import statistics
from waps import sampler as samp
import time
from fractions import Fraction
from gmpy2 import mpq,mpfr,log2


import signal
from contextlib import contextmanager

TimeoutError = "Timeout ..need to stop"

@contextmanager

def timeout(time):
    # Register a function to raise a TimeoutError on the signal.
	signal.signal(signal.SIGALRM, raise_timeout)
    # Schedule the signal to be sent after ``time``.
	signal.alarm(time)
	try:
		yield
	except Exception as e:
		print(e)
		if e == TimeoutError:
			pass
	finally:
        # Unregister the signal so it won't be triggered
        # if the timeout is not reached.
		signal.signal(signal.SIGALRM, signal.SIG_IGN)

def raise_timeout(signum, frame):
    raise Exception(TimeoutError)

SAMPLER_WAPS = 1
SAMPLER_SPUR = 2

COUNTER_GANAK = 1

class Experiment:
	def __init__(self, Xlist, Ylist, Zlist ):
		self.Xlist = Xlist
		self.Ylist = Ylist
		self.Zlist = Zlist
		
	def doSmallCheck(self, modelcountX, content_count, content_sample, seed):
		
		m = len(self.Ylist)
		n = len(self.Xlist)+len(self.Zlist)
		
		yassign_list, yclause_list = GenerateSamples(1, self.Xlist, self.Zlist, self.Ylist, content_sample, seed)
		lowEnt = False
		entropy = 0
		T = 0
		
		for i in range(1):
			
			yassign = yassign_list[i]
			yclause = yclause_list[i]
			count = ComputeCount(content_count.strip("\n")+"\n"+yassign)
			count = mpq(count)
			modelcountX = mpq(modelcountX)
			
			r = mpfr(count/modelcountX,100)
			
			if args.verbose:
				print("modelcount for fix Y",count)
				print("value of r ",r)
				
			if r > 0.5:
				print("found high r, calculating new value of t....")
				T = (6 / (math.pow(args.epsilon,2)))  * (min  ((n / math.pow(math.log(1/(1-r),2),2)), (m + math.log(m + math.log(m,2)+2.5,2))))
				if args.verbose:
					print("new value of T",T)
				
				lowEnt = True
		return lowEnt, T, yassign, r
	
	def oneExperiment(self, modelcountX, content_count, content_sample, lowEnt, r, T, seed):		
		entropy = 0
		yassign_list, yclause_list = GenerateSamples(int(math.ceil(T)), self.Xlist, self.Zlist, self.Ylist, content_sample, seed)
		for i in range(math.ceil(T)):
			yassign = yassign_list[i]
			yclause = yclause_list[i]
			count = ComputeCount(content_count.strip("\n")+"\n"+yassign)
			count = mpq(count)
			modelcountX = mpq(modelcountX)
			
			r = mpfr(count/modelcountX,100)
			r_1 = mpfr(modelcountX/count,100)
			entropy += log2(r_1)
			if args.verbose == 2:
				print("modelcount for fix Y",count)
				print("value of r ",r)
		if lowEnt:
			entropy = ((1-r)*mpfr(entropy/T,100))+(r*log2(r_1))
			print("entropy",entropy)
		else:
			print("total entropy",mpfr(entropy/T,100))
			entropy = mpfr(entropy/T,100)
		return entropy




def ComputeCount(content):
	inputfile = args.input
	inputfile = inputfile.split("/")[-1]
	tempfilename = tempfile.gettempdir() + '/' + inputfile + "_countfile"
	#tempfilename = inputfile + "_countfile"
	outfile = tempfile.gettempdir() + '/' + inputfile + "_count"

	
	f = open(tempfilename,"w")
	f.write(content)
	f.close()

	
	#print("using ganak to do projected counting..\n")
	cmd = "./dependencies/run_ganak.sh  %s > %s " %(tempfilename,outfile)
	find_string = "s mc"


	os.system(cmd)

	with open(outfile) as f:
		lines = f.readlines()
	f.close()

	



	for line in lines:
		if line.startswith(find_string):
			modelcount = int(line.strip(find_string).strip(" ").strip("\n"))
			break
	os.unlink(outfile)
	os.unlink(tempfilename)

	return modelcount



def ParserInput(inputfile):
	Xlist = []
	Ylist = []
	with open(inputfile) as f:
		lines = f.readlines()
	f.close()

	'''
	c max represents X variables and c ind presents Y variables
	'''

	for line in lines:
		if line.startswith("c max"):
			Xlist.append(line.strip("c max").strip(" ").strip("\n").split(" ")[:-1])
			continue
		if line.startswith("c ind"):
			Ylist.append(line.strip("c ind").strip(" ").strip("\n").split(" ")[:-1])
			continue
		if line.startswith("p cnf"):
			totalvars = int(line.split(" ")[2])
			totalclause = int(line.split(" ")[3])
			continue

	
	Xlist =  [int(y) for x in Xlist for y in x]
	Ylist =  [int(y) for x in Ylist for y in x]
	Zlist =  []
	Zind = 'c ind '
	for var in range(totalvars):
		var = var + 1
		if (var not in Xlist) and (var not in Ylist):
			Zlist.append(var)
			Zind += "%d " %(var)

	Zind += "0\n"

	if len(Zlist) == 0:
		Zind =''


	if args.verbose:
		print("len of X",len(Xlist))
		print("len of Y",len(Ylist))
		print("len of Z",len(Zlist))
		
		if args.sampler == 1:
			print("for sampling going to use WAPS")
		if args.sampler == 2:
			print("for sampling going to use SPUR")

		
		if args.counter ==1:
			print("for counting going to use Ganak")

	if args.verbose==2:
		print("Xlist",Xlist)
		print("Ylist",Ylist)
		print("Zlist",Zlist)
	

	return Xlist, Ylist, Zlist, Zind, totalclause

def GenerateSamples(num_sample, Xlist, Zlist, Ylist, content, seed):
	inputfile = args.input
	inputfile = inputfile.split("/")[-1]
	tempfilename = tempfile.gettempdir() + '/' + inputfile+"_samples"
	outfile = tempfile.gettempdir() + '/' + inputfile+"_samplesout"

	f = open(tempfilename,"w")
	f.write(content)
	f.close()

	yassign_list = []
	yclause_list = []
	

	
	if args.sampler == 1:
		samples = []
		sampler = samp(cnfFile=tempfilename)
		sampler.compile()
		sampler.parse()
		sampler.annotate()
		solList = sampler.sample(totalSamples=num_sample) 
		solList = [i.strip().split() for i in list(solList)]
		for i in range(len(solList)):
			solList_one = sorted(list(map(int,solList[i])),key=abs)
			solList_x = []
			y_assignment = ''
			yclause = ''
			for var in solList_one :
				if (abs(var) in Xlist) or (abs(var) in Zlist):
					solList_x.append(var)
				if abs(var) in Ylist:
					y_assignment += "%d 0\n" %(var)
					yclause += "%d " %(-1*var)
			samples.append(solList_x)
			yassign_list.append(y_assignment)
			yclause_list.append(yclause+"0\n")

	if args.sampler == 2:
		cmd = './dependencies/spur -seed %d -q -s %d -out %s -cnf %s' % (seed, num_sample, outfile, tempfilename)
		os.system(cmd)
		with open(outfile, 'r') as f:
			lines = f.readlines()

		solList = []
		samples = []
		startParse = False


		for line in lines:
			if (line.startswith('#START_SAMPLES')):
				startParse = True
				continue
			if (not(startParse)):
				continue
			if (line.startswith('#END_SAMPLES')):
				startParse = False
				continue

			fields = line.strip().split(',')
			solCount = int(fields[0])
			sol = ''
			i = 1
			for x in list(fields[1]):
				if (x == '0'):
					sol += '-'+str(i)+" "
				if (x == '1'):
					sol += str(i)+" "
				i += 1

			sol = sol.strip(" ")

			for i in range(solCount):
				solList.append(sol)
		


		for i in range(len(solList)):
			solList_one = sorted(list(map(int,solList[i].split(" "))),key=abs)
			solList_x = []
			y_assignment = ''
			yclause = ''
			for var in solList_one :
				if (abs(var) in Xlist) or (abs(var) in Zlist):
					solList_x.append(var)
				if abs(var) in Ylist:
					y_assignment += "%d 0\n" %(var)
					yclause += "%d " %(-1*var)
			samples.append(solList_x)
			yassign_list.append(y_assignment)
			yclause_list.append(yclause+"0\n")

		#os.unlink(outfile)

	
	#os.unlink(tempfilename)
	return yassign_list,yclause_list


def nCr(n,r):
	f = math.factorial
	return f(n) // f(r) // f(n-r)

def compute_itr():
	itr = 9999
	delta = 0.9 * args.delta
	for i in range(4,101):
		tot = 0
		for j in range(i//2+1,i+1,1):
			tot+= nCr(i,j)*(args.prob)**(j)*(1-args.prob)**(i-j)
	    #print(tot)
		if 1-delta < tot:
			itr = i
			break
	print("number of iteration %d with delta %f" %(itr, 0.9*args.delta))


	return itr




def estimate():

	# As per benchmarks
	# Y is unique with respect to X and Z

	start = time.time()
	Xlist, Ylist, Zlist, Zind, num_clause = ParserInput(args.input)

	m = len(Ylist)
	n = len(Xlist) + len(Zlist)

	with open(args.input) as f:
		content = f.read()
	f.close()

	### projected model count on X and Z ###

	content = content.replace('c ind','c y')
	content = content.replace('c max','c ind')
	content = Zind + content
	
	modelcountX = ComputeCount(content)
	if args.verbose:
		print("total projected count on X and Z ",modelcountX)

	

	content_sample = content
	### Ganak considers the header: update the header with clause = clause + numberof Y
	n_var = len(Xlist)+len(Ylist)+len(Zlist)
	content_count = content.replace('p cnf %d %d' %(n_var,num_clause),'p cnf %d %d' %(n_var,num_clause+len(Ylist)))

	
	entropylist = []

	t = math.ceil( math.log(10/args.delta,2))

	print("checking lowEnt for %d itr " %(t))

	exp = Experiment(Xlist = Xlist, Ylist = Ylist, Zlist = Zlist)


	for i in range(t):
		seed = (args.seed) * i
		lowEnt, T, yclause, r = exp.doSmallCheck(modelcountX, content_count, content_sample, seed)
		if lowEnt:
			content_sample += yclause
			content_count  += yclause
			content_count = content_count.replace('p cnf %d %d' %(n_var,num_clause+len(Ylist)),'p cnf %d %d' %(n_var,num_clause+len(Ylist)+1))
			break

	if not lowEnt:
		
		T = (6 / (math.pow(args.epsilon,2)))  * ( min  (n, (m + math.log(m + math.log(m,2)+1.1,2))) - 1)
		
	
	itr = compute_itr()

	if args.verbose:
		print("number of itr",itr)
		print("Is low ent",lowEnt)
		print("number of samples",T) 

	


	for oneitr in range(itr):
		seed = (args.seed) * oneitr		
		one_entropy = exp.oneExperiment(modelcountX, content_count, content_sample, lowEnt, r, T, seed)
		entropylist.append(one_entropy)

	entropylist = np.array(entropylist)
	print(entropylist)
	entropy = np.median(entropylist)
	print("final entropy", entropy)
	print("total time taken",time.time()-start)



	

if __name__ == "__main__":

	
	samplers = str(SAMPLER_WAPS) +" for WAPS\n"
	samplers += str(SAMPLER_SPUR) +" for SPUR\n"
	counter = str(COUNTER_GANAK) + " for GANAK\n"
	parser = argparse.ArgumentParser()
	parser.add_argument('--seed', type=int, default=10, dest='seed')
	parser.add_argument('--verb', type=int, help="0 ,1 ,2", dest='verbose')
	parser.add_argument('--epsilon', type=float, default=0.8, dest='epsilon')
	parser.add_argument('--delta', type=float, default=0.09, dest='delta')
	parser.add_argument('--sampler', type=int, help=samplers, default=SAMPLER_SPUR, dest='sampler')
	parser.add_argument('--counter', type=int, help=counter, default=COUNTER_GANAK, dest='counter')
	parser.add_argument('--prob',type=float,default=0.84,dest='prob')
	parser.add_argument("input", help="input file")
	parser.add_argument('--timeout', type=int, default=3000, help = "timeout; default 3000s", dest='timeout')
	args = parser.parse_args()
	with timeout(args.timeout):
		estimate()