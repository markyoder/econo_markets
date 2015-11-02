from math import *	
from scipy import *
import scipy
from pylab import *
from matplotlib import *
import numpy.fft as nft
import scipy.optimize as spo
from matplotlib import pyplot as plt

import string
import sys
#from matplotlib import *
#from pylab import *
import os
import random
import time
#
# gamma function lives here:
#import scipy.special
from scipy.special import gamma
#from scipy.optimize import leastsq
from matplotlib import axis as aa
from threading import Thread
#
#
import datetime
import calendar

#####################
## utilities:		
def datetimeFromStrings(strDt, strTm, dtdelim='/'):
	# looks like date[time].strptime(date_string, format) does this...
	if strTm=='': strTm='00:00:00.0'
	#
	ldt=strDt.split(dtdelim)
	ltm=strTm.split(':')
	
	lsecs=ltm[2].split('.')
	secs = long(lsecs[0])
	if secs==60:
		secs=59
		# simple work-around. 
	msecs=0
	if len(lsecs)>1:
		msecs = long(float("."+lsecs[1])*10**6)
	#
	#return datetime.datetime(long(ldt[0]), long(ldt[1]), long(ldt[2]), long(ltm[0]), long(ltm[1]), long(ltm[2]))
	return datetime.datetime(long(ldt[0]), long(ldt[1]), long(ldt[2]), long(ltm[0]), long(ltm[1]), secs, msecs)

def floatToDateTime(fdate):
	# this might be off by a few hundred microseconds. in one example, i get 890000-> 889993
	datepart=datetime.date.fromordinal(long(fdate))
	hrs=24.0*(fdate-long(fdate))
	hr=long(hrs)
	mins=(hrs-hr)*60
	mn=long(mins)
	secs=(mins-mn)*60
	sec=long(secs)
	microsecs=long((secs-sec)*10**6)
	
	timepart=datetime.time(hr,mn,sec, microsecs)
	
	return datetime.datetime.combine(datepart, datetime.time(hr, mn, sec, microsecs))

def datetimeToFloat(dtm):
	# return number of days, including fractional bit.
	return float(dtm.toordinal()) + float(dtm.hour)/24.0 + float(dtm.minute)/(24.0*60.0) + float(dtm.second)/(24.0*60.0*60.0) + float(dtm.microsecond)/(24.0*60*60*10**6)

################
################


def plotRecords(recfile="econoRB/djia10/rb-djiaClose-ts128.dat"):
	# plot records with date axis. maybe just use the log-values for the ratios?
	# data come as : fdt, nEvents, Nbig, Nsmall, Npois-big, Npois-small.
	# note the fdt, nEvents may be reversed from the file header
	# also note that since these files are likely in 1-day resolution it's not too important (there will be a minor weekend error).
	#
	f=open(recfile, 'r')
	recs=[[],[],[],[],[],[]]
	# and i guess we know what we want to plot already:
	Xdt=[]
	Y=[]
	Xdjia=[]
	Ydjiaclose=[]	# this is really a RB ratio, not the close values.
	#
	# and these will be for our trading game...
	djiaOp=[]	# open, closing values for djia. we will win or lose P*(close-open)/open
	djiaCl=[]
	minRatioDate=None
	for rw in f:
		if rw[0]=="#" or rw[0]==" " or rw[0]=="\n" or rw[0]=="\t": continue
		rws=rw.split("\t")
		
		if minRatioDate==None or float(rws[0])<minRatioDate: minRatioDate=float(rws[0])
		
		for i in xrange(len(rws)):
			if i >= len(recs): continue
			recs[i]+=[float(rws[i])]
			#
		#
		Xdt+=[floatToDateTime(float(rws[0]))]
		Y+=[log10(float(rws[2])/float(rws[3]))]
	f.close()
	fdj=open("stockdata/djia-history.csv")
	for rw in fdj:
		if rw[0]=="#" or rw[0]==" " or rw[0]=="\n" or rw[0]=="\t" or rw[0]=="D": continue
		rws=rw.split(",")
		dY=10*((float(rws[4])-float(rws[1]))/float(rws[1]))
		djiaOp.insert(0,float(rws[1]))
		djiaCl.insert(0,float(rws[4]))
		if dY<.5: continue
		#print "lens: %d, %d" % (len(Xdjia), len(Ydjiaclose))
		Xdjia+=[datetimeFromStrings(rws[0], "00:00:00.0", dtdelim='-')]
		Ydjiaclose+=[dY]
	fdj.close()
	#
	#
	plt.figure(0)
	plt.clf()
	plt.plot_date(Xdt, Y, fmt='b-', xdate=True, ydate=False)
	plt.plot_date(Xdjia, Ydjiaclose, fmt='go', xdate=True, ydate=False)
	#plt.show()
	#
	#
	# now, do a trading experiment:
	# start with $10,000 in stock, $10,000 in cash.
	# r() is our rb indicator, NRBbig/NRBsmall
	# if r()>.5, buy 100 shares. if r()<.5 sell 100 shares (if we have the cash/shared. we may need to adjust the trading volume. maybe 10% of what we've got, etc.)...
	# yeah, let's to that. buy/sell 10% of cash/stocks.
	#
	cash=301.49
	stocks=cash
	cashes=[cash/2.0]
	stockses=[stocks/2.0]
	investors=[cash]
	totals=[cash]
	# now, spin through our stocks and RB data sets. trade on the above rules.
	
	tradeArrays=[[],[],[],[]]		#dt, r, open, close
	djiaStart=0	# starting positon for the djia (aka, where our djia date=min(r-date)
	#print minRatioDate
	
	for i in xrange(len(Xdt)):
		thisfdt=datetimeToFloat(Xdt[i])
		if thisfdt<minRatioDate: continue
		#
		tradeArrays[0]+=[thisfdt]
		tradeArrays[1]+=[Y[djiaStart]]	# aka, the log(rb ratio)
		tradeArrays[2]+=[djiaOp[i]]
		tradeArrays[3]+=[djiaCl[i]]
		
		# i guess we don't really need the arrays, but maybe later...
		# do the trading in the same swoop:
		buyLim=.5
		sellLim=.2
		tradeFactor=.1	# how much of the stock/cash portfolio to trade
		commissionFactor=0
		dX=0
		dY=(tradeArrays[3][-1]-tradeArrays[2][-1])/tradeArrays[2][-1]		# percent change in index value from open to close.
		if tradeArrays[1][-1]>=buyLim:
			# buy:
			dX=cashes[-1]*tradeFactor

		if tradeArrays[1][-1]<sellLim:
			# sell:
			dX=-stockses[-1]*tradeFactor
		#dX=0
		#print "transaction: %f, %f, %f" % (stockses[-1], cashes[-1], dX)
		cashes+=[cashes[-1]-dX]
		stockses+=[stockses[-1]+dX]
		if abs((cashes[-2]+stockses[-2])-(cashes[-1]+stockses[-1]))>.000001: print "cassh conservation error... %f" % ((cashes[-2]+stockses[-2])-(cashes[-1]+stockses[-1]))
		investors+=[investors[-1]*(1.0+dY)]
		#investors+=[djiaCl[i]]
		#investors+=[tradeArrays[3][-1]]
		#
		# and stocks do the market:
		stockses[-1]*=(1.0+dY)
		# and pay a commission:
		cashes[-1]-=abs(dX)*commissionFactor
		totals+=[stockses[-1]+cashes[-1]]
		djiaStart+=1
		#print "final: %f, %f/%f, %f, %f, %f gain: %f" % (tradeArrays[1][-1], stockses[-1], stockses[-1]/(1.0+dY), cashes[-1], 1.0+dY, investors[-1], (stockses[-1]-stockses[-2] + cashes[-1]-cashes[-2]))
		#if i>=500: break
	#
	# now, let's look at stockses and cashes:
	plt.figure(1)
	plt.clf()
	plt.plot_date(tradeArrays[0], stockses[1:], "b-", xdate=True, ydate=False)
	plt.plot_date(tradeArrays[0], cashes[1:], "g-", xdate=True, ydate=False)
	plt.plot_date(tradeArrays[0], investors[1:], "r-", xdate=True, ydate=False)
	plt.plot_date(tradeArrays[0], totals[1:], "y-", xdate=True, ydate=False)
	plt.show()
	
	print "final: stock: %f, cash: %f, index: %f/%f, gain: %f" % (stockses[-1], cashes[-1], tradeArrays[3][-1],investors[-1], (stockses[-1]+cashes[-1])/tradeArrays[3][-1])
	#print "start: %f, %f, %f, %f" % (tradeArrays[0][0], tradeArrays[1][0], tradeArrays[2][0], tradeArrays[3][0])
		
	#print "i, i`: %d, %d" % (len(Xdt), djiaStart)
	
	
			


