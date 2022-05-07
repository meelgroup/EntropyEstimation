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
import numpy as median
import statistics
from waps import sampler as samp
import time
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
	
	def oneExperiment(self, modelcountX, content_count, content_sample,  T , seed):		
		entropy = 0
		samples, yassign_list, yclause_list = GenerateSamples(T, self.Xlist, self.Zlist, self.Ylist, content_sample, seed)
		yassign = yassign_list[0]
		yclause = yclause_list[0]
		count = ComputeCount(content_count.strip("\n")+"\n"+yassign)
		count = mpq(count)
		modelcountX = mpq(modelcountX)
		r = mpfr(count/modelcountX)
		r_1 = mpfr(modelcountX/count)
		entropy += log2(r_1)
		if args.verbose == 2:
			print("modelcount for fix Y",count)
			print("value of r ",r)
		return entropy, content_sample+yclause




def ComputeCount(content):
	inputfile = args.input
	inputfile = inputfile.split("/")[-1]
	tempfilename = tempfile.gettempdir() + '/' + inputfile + "_countfile"
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

	''''
	Y is unique in terms of X and Z variables
	'''
	
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


		if args.counter == 1:
			print("for counting going to use Ganak")

	if args.verbose==2:
		print("Xlist",Xlist)
		print("Ylist",Ylist)
		print("Zlist",Zlist)
	

	return Xlist, Ylist, Zlist, Zind, totalclause

def GenerateSamples(num_sample, Xlist, Zlist, Ylist, content, seed):
	inputfile = args.input
	inputfile = inputfile.split("/")[-1]
	tempfilename = tempfile.gettempdir() + '/' + inputfile + "_samples"
	outfile = tempfile.gettempdir() + '/' + inputfile + "_samplesout"
	
	#print("num_sample",num_sample)

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
				else:
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

		os.unlink(outfile)

	
	os.unlink(tempfilename)
	return np.array(samples), yassign_list,yclause_list



def estimate():

	# As per benchmarks
	# Y is unique with respect to X and Z

	start = time.time()
	
	Xlist, Ylist, Zlist, Zind, num_clause = ParserInput(args.input)


	with open(args.input) as f:
		content = f.read()
	f.close()

	### projected model count on Y ###
	content_Y = content.replace('c max','c ')	
	T = ComputeCount(content_Y)
	print("Required number of samples",math.ceil(T))


	### projected model count on X and Z ###

	content = content.replace('c ind','c y')
	content = content.replace('c max','c ind')
	content = Zind + content
	
	modelcountX = ComputeCount(content)

	content_sample = content.replace('c y','c ind')
	

	### Ganak considers the header: update the header with clause = clause + numberof Y
	n_var = len(Xlist)+len(Ylist)+len(Zlist)
	content_count = content.replace('p cnf %d %d' %(n_var,num_clause),'p cnf %d %d' %(n_var,num_clause+len(Ylist)))

	if args.verbose:
		print("total projected count on X and Z ",modelcountX)


	exp = Experiment(Xlist = Xlist, Ylist = Ylist, Zlist = Zlist)

	entropy = 0

	for i in range(T):
		one_entropy, content_sample = exp.oneExperiment(modelcountX, content_count, content_sample, 1, args.seed)
		entropy += one_entropy
	
	
	entropy = mpq(entropy)
	final_entropy = mpfr(entropy/T,100)
	print("final entropy", final_entropy)
	print("total time taken",time.time()-start)



	

if __name__ == "__main__":
	samplers = str(SAMPLER_WAPS) +" for WAPS\n"
	samplers += str(SAMPLER_SPUR) +" for SPUR\n"
	
	counter = str(COUNTER_GANAK) + " for GANAK\n"
	parser = argparse.ArgumentParser()
	parser.add_argument('--seed', type=int, default=10, dest='seed')
	parser.add_argument('--verb', type=int, help="0 ,1 ,2", dest='verbose')
	parser.add_argument('--epsilon', type=float, default=0.3, dest='epsilon')
	parser.add_argument('--delta', type=float, default=0.1, dest='delta')
	parser.add_argument('--sampler', type=int, help=samplers, default=SAMPLER_SPUR, dest='sampler')
	parser.add_argument('--counter', type=int, help=counter, default=COUNTER_GANAK, dest='counter')
	parser.add_argument('--timeout', type=int, default=3000, help = "timeout; default 3000s", dest='timeout')
	parser.add_argument("input", help="input file")
	args = parser.parse_args()
	with timeout(args.timeout):
		estimate()
		