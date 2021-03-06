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

import rbIntervals as rbi
import recordBreaker as rbb

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


def plotRecords(recfile="econoRB/djia10/rb-djiaClose-ts64.dat"):
	# plot records with date axis. maybe just use the log-values for the ratios?
	# data come as : fdt, nEvents, Nbig, Nsmall, Npois-big, Npois-small.
	# note the fdt, nEvents may be reversed from the file header
	# also note that since these files are likely in 1-day resolution it's not too important (there will be a minor weekend error).
	#
	outputFolder = "econoRB/djia10"
	f=open(recfile, 'r')	# given a file of record-breaking data. we're going to plot this and use it as a guide for trading.
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
	djiaHi=[]
	djiaLo=[]
	minRatioDate=None
	for rw in f:
		if rw[0]=="#" or rw[0]==" " or rw[0]=="\n" or rw[0]=="\t": continue
		rws=rw.split("\t")
		#
		# get the minimum date in the RB data.
		if minRatioDate==None or float(rws[0])<minRatioDate: minRatioDate=float(rws[0])
		#
		for i in xrange(len(rws)):
			if i >= len(recs): continue
			recs[i]+=[float(rws[i])]
			#
		#
		Xdt+=[floatToDateTime(float(rws[0]))]		# x axis dates
		Y+=[log10(float(rws[2])/float(rws[3]))]	# ratio of nrb_big/nrb_small (trend)
	f.close()
	fdj=open("stockdata/djia-history.csv")
	Xdjia2=[]
	nrows=0
	for rw in fdj:
		nrows+=1
		if rw[0]=="#" or rw[0]==" " or rw[0]=="\n" or rw[0]=="\t" or rw[0]=="D": continue
		rws=rw.split(",")
		dY=10*((float(rws[4])-float(rws[1]))/float(rws[1]))
		#djiaOp.insert(0,float(rws[1]))
		#djiaCl.insert(0,float(rws[4]))
		djiaOp+=[float(rws[1])]
		djiaCl+=[float(rws[4])]
		djiaHi+=[float(rws[2])]
		djiaLo+=[float(rws[3])]
		Xdjia2+=[datetimeFromStrings(rws[0], "00:00:00.0", dtdelim='-')]
		if dY<.5: continue
		#print "lens: %d, %d" % (len(Xdjia), len(Ydjiaclose))
		Xdjia+=[datetimeFromStrings(rws[0], "00:00:00.0", dtdelim='-')]
		Ydjiaclose+=[dY]
	fdj.close()
	djiaOp.reverse()
	djiaCl.reverse()
	Xdjia2.reverse()
	#
	#
	plt.figure(0)
	plt.clf()
	plt.title("Ratios of nrb_big/nrb_small")
	plt.ylabel("r(t)=nrb_big/nrb_small")
	plt.xlabel("date")
	plt.plot_date(Xdt, Y, fmt='b-', xdate=True, ydate=False)
	plt.plot_date(Xdjia, Ydjiaclose, fmt='go', xdate=True, ydate=False)
	#plt.show()
	#
	#
	# now, do our first trading experiment:
	# start with some cash and some stocks. we typically start evenly distributed.
	# arguably, we need to experiment with different initial conditions. the 1929 crash occurs shortly after our experiment, which might skew our results.
	# r() is our rb indicator, NRBbig/NRBsmall
	# if r()>.5, buy 100 shares. if r()<.5 sell 100 shares (if we have the cash/shared. we may need to adjust the trading volume. maybe 10% of what we've got, etc.)...
	# yeah, let's to that. buy/sell 10% of cash/stocks.
	#
	cash=301.0
	stocks=cash
	cashes=[cash/1000.0]
	stockses=[stocks/1.0]
	investors=[cash]
	totals=[cash]
	tradeGains=[]	# how does our method compare to the market? (not implemented yet).
	#marketGains=[]
	# now, spin through our stocks and RB data sets. trade on the above rules.
	
	tradeArrays=[[],[],[],[]]		#dt, r, open, close
	djiaStart=0	# starting positon for the djia (aka, where our djia date=min(r-date)
	#print minRatioDate
	
	#for i in xrange(len(Xdt)):
	#i=0
	for i in xrange(len(Xdjia2)):
	#while i<xrange(nrows-1):
		#thisfdt=datetimeToFloat(Xdt[i])
		thisfdt=datetimeToFloat(Xdjia2[i])
		if thisfdt<=minRatioDate: continue
		#
		tradeArrays[0]+=[thisfdt]
		tradeArrays[1]+=[Y[djiaStart]]	# aka, the log(rb ratio)
		tradeArrays[2]+=[djiaOp[i]]
		tradeArrays[3]+=[djiaCl[i]]
		
		# i guess we don't really need the arrays, but maybe later...
		# do the trading in the same swoop:
		buyLim=.7
		sellLim=.6
		buyFactor=.1	# how much of the stock/cash portfolio to trade
		sellFactor=.1
		commissionFactor=0
		dX=0
		dY=(tradeArrays[3][-1]-tradeArrays[2][-1])/tradeArrays[2][-1]		# percent change in index value from open to close.
		if tradeArrays[1][-1]>=buyLim:
			# buy:
			dX=cashes[-1]*buyFactor

		if tradeArrays[1][-1]<sellLim:
			# sell:
			dX=-stockses[-1]*sellFactor
		#dX=0
		#print "transaction: %f, %f, %f" % (stockses[-1], cashes[-1], dX)
		cashes+=[cashes[-1]-dX]
		stockses+=[stockses[-1]+dX]
		if abs((cashes[-2]+stockses[-2])-(cashes[-1]+stockses[-1]))>.000001: print "cassh conservation error... %f" % ((cashes[-2]+stockses[-2])-(cashes[-1]+stockses[-1]))
		investors+=[investors[-1]*(1.0+dY)]
		marketGain=(investors[-1]-investors[-2])/investors[-2]
		#investors+=[djiaCl[i]]
		#investors+=[tradeArrays[3][-1]]
		#
		# and stocks do the market:
		stockses[-1]*=(1.0+dY)
		# and pay a commission:
		cashes[-1]-=abs(dX)*commissionFactor
		totals+=[stockses[-1]+cashes[-1]]
		tradeGains+=[((totals[-1]-totals[-2])/totals[-2])-marketGain]
		
		djiaStart+=1
		#print "final: %f, %f/%f, %f, %f, %f gain: %f" % (tradeArrays[1][-1], stockses[-1], stockses[-1]/(1.0+dY), cashes[-1], 1.0+dY, investors[-1], (stockses[-1]-stockses[-2] + cashes[-1]-cashes[-2]))
		#if i>=500: break
		
		#i+=1
	#
	# now, let's look at stockses and cashes:
	plt.figure(1)
	plt.clf()
	plt.plot_date(tradeArrays[0], stockses[1:], "b-", xdate=True, ydate=False)
	plt.plot_date(tradeArrays[0], cashes[1:], "g-", xdate=True, ydate=False)
	plt.plot_date(tradeArrays[0], investors[1:], "r-", xdate=True, ydate=False)
	plt.plot_date(tradeArrays[0], totals[1:], "y-", xdate=True, ydate=False)
	
	plt.figure(2)
	plt.clf()
	plt.plot_date(tradeArrays[0], tradeGains, "b-", xdate=True, ydate=False)
	#plt.show()
	
	print "final: stock: %f, cash: %f, index: %f/%f, gain: %f" % (stockses[-1], cashes[-1], tradeArrays[3][-1],investors[-1], (stockses[-1]+cashes[-1])/tradeArrays[3][-1])
	#print "start: %f, %f, %f, %f" % (tradeArrays[0][0], tradeArrays[1][0], tradeArrays[2][0], tradeArrays[3][0])
		
	#######################################################
	#######################################################
	#
	# a second experiment (because the first is a big loser):
	# theory:
	# - nrb_large is actually slightly less than poisson up to somewhere between 64<nN128. after n=128, nrb_large appers to converge to poisson.
	# - nrb_small deviates significantly from poisson through n=128.
	# gains: large one time gains are rare but long runs of record highs are disproportionately common. gains occur in long strings of small gains.
	# losses: losses are catastrophic. strings of low records are more common than random; large losses are more common than random. losses occur suddenly.
	# - so, develop separate buy/sell rules to take advantage of this asymmetry. let's try to buy on a trend, but use stop-loss
	# rules for selling. in otherwords, we buy when we detect a positive trend over Nevents. each time we break a record value, we set a
	# sell price. if the buy/sell instructions conflict... well, then i don't know what. we'll play with that.
	#
	# note: this experiement depends on the data file provided in the parameter. the affirmative results obtained earlier used Closing values. what about ClOpen values?
	#
	
	startPos=10000	# skip this many rows before trading begins.
	#startPos=5000
	
	cash=301.0
	stocks=cash
	cashes=[cash/2.0]
	stockses=[stocks/2.0]
	investors=[cash]
	totals=[cash]
	#rbWinLen=64		# (parameter) for now, set this value. we will need to keep track of our own record values, as our data file only contains nrb.
	rbWinLen=len(Xdjia2)-len(tradeArrays[0])		#... or actually, i think if we use the already built in minRatioDate thing, which accounts for the lag in the raw data
	print "auto-rbLen: %d" % rbWinLen
							# and the RB data.
	#rbVals=[]	# large record-breaking prices.
	rbVal=djiaCl[0]	# rbVal represents our minimum-hold price. aka, below this price, we start selling (agressively). rbVal goes up with a record price, down with a sell.
	# now, spin through our stocks and RB data sets. trade on the above rules.
	
	tradeArrays=[[],[],[],[],[],[]]		#dt, r, open, close, hi, lo
	djiaStart=0	# starting positon for the djia (aka, where our djia date=min(r-date)
	#print minRatioDate
	
	#for i in xrange(len(Xdt)):
	#i=0
	nTradeInstructionConflicts=0
	nBuys=0
	nSells=0
	#
	avLen=int(2*math.log(rbWinLen))
	#avLen=10
	meanvals=[0,0,0,0]	# so these will be 0,1,2,3 according ot tradeArrays[0],...,[3]
	buyLim=.1
	sellLim=-.5			#RB-rule limit
	sellRBfactor=.95	# aka, sell at price<rbVal*sellRBfactor
	buyFactor=.2		# how much of the stock/cash portfolio to trade
	sellFactor=.90 	# how much of the stock/cash portfolio to trade
	commissionFactor=0.000000001	# just an estimate. i know, most commissions these days are fixed.

	for i in xrange(len(Xdjia2)):
	#while i<xrange(nrows-1):
		#thisfdt=datetimeToFloat(Xdt[i])
		thisfdt=datetimeToFloat(Xdjia2[i])
		#if djiaCl[i]>rbVal: rbVal=djiaCl[i]
		#oldrbVal=rbVal
		if thisfdt<minRatioDate: continue
		if i<startPos: continue
		#
		# now, smoothe r(t) over avlen preceeding elements. let's just leave the first aveLen elements rough (average over whatever we've got).
		aveStart=djiaStart-avLen
		if aveStart<0: aveStart=0
		aveRatio=sum(Y[aveStart:djiaStart+1])/float(len(Y[aveStart:djiaStart+1]))	# do the float conversion just in case...
		
		#print "ratio vals: %s, %f\n" % (Y[aveStart:djiaStart+1], aveRatio)
		
		tradeArrays[0]+=[thisfdt]
		#tradeArrays[1]+=[Y[djiaStart]]	# aka, the log(rb ratio). r typicall ranges from about .05 to 50, and we look at the log of that value.
		tradeArrays[1]+=[aveRatio]
		tradeArrays[2]+=[djiaOp[i]]
		tradeArrays[3]+=[djiaCl[i]]
		tradeArrays[4]+=[djiaHi[i]]
		tradeArrays[5]+=[djiaLo[i]]
		
		if len(tradeArrays[0])<2: continue	# we're looking at close to prev-close, so we need a buffer row.
		
		# i guess we don't really need the arrays, but maybe later...
		# do the trading in the same swoop:
		
		if tradeArrays[3][-2]>rbVal: rbVal=tradeArrays[3][-2]		# if yesterday's close or today's open is higher than the current highest value...
		if tradeArrays[2][-1]>rbVal: rbVal=tradeArrays[2][-1]		# note, we can't use yesterday's high unless we know if a sell occured after that high occured (which resets
																					# the rbVal). maybe we can catch the high value (in this limited res simulation) in the BUY rule?
		oldrbVal=rbVal
		didSell=0 # note: our sell rule will trump our buy rule, so be careful with this flag.
		#
		dX=0
		#dY=(tradeArrays[3][-1]-tradeArrays[2][-1])/tradeArrays[2][-1]		# percent change in index value from open to close. noe: we may need a separate dY for the trader...
		#dY=tradeArrays[3][-1]/tradeArrays[2][-1]
		dY=tradeArrays[3][-1]/tradeArrays[3][-2]	# percent change in the market...
		dYtrader=dY
		# stocks do the market:
		# this may be a problem. are we trading on the morning's price when we sell?
		investors+=[investors[-1]*dY]		# isn't open/close a more efficient way of doing this?
		#stockses[-1]*=(1.0+dY)
		
		# buy rule:
		#if tradeArrays[1][-1]>=buyLim or tradeArrays[3][-1]/tradeArrays[3][-avLen]>=(1.0/sellRBfactor):
		if tradeArrays[1][-1]>=buyLim or tradeArrays[3][-1]/tradeArrays[3][-avLen]>=1.05:
		#if tradeArrays[1][-1]>=buyLim:
			# (we'll buy at the very end of the day based on today's r(t) value). note, we buy at today's close.
			# buy:
			dX=cashes[-1]*buyFactor		# note: we'll buy at the end of the day; this is the number of $1 shares to buy.
			nBuys+=1	#... unless we get trumped by the sell.
		#
		# note that as per sequence, "sell" will trump buy, but let's throw a note if they contradict:
		#if dX>0:
		#	 #print "buy and sell order; default to sell(%d): r(t)=%f, open=%f, rbVal=%f" % (nTradeInstructionConflicts, tradeArrays[1][-1], tradeArrays[3][-1], rbVal)
		#	 nTradeInstructionConflicts+=1
		#
		# sell rule:
		sellPrice=rbVal*sellRBfactor
		if tradeArrays[3][-1]<sellPrice:
		#if tradeArrays[5][-1]<sellPrice:
			# we're ready to sell all day. if the price dips below our sell factor, we sell. this could be good or bad for us really.
			didSell=1
			# if close value<sellLim 
			# what is the open price? if the open price is less than our sell price, that is our trading price.
			#if tradeArays[2][-1]<sellPrice: sellPrice=tradeArays[2][-1]	# but for now, assume continuity...
			# sell:
			# stocks adjust to new sell value:
			dYtrader=sellPrice/tradeArrays[3][-2]	# percent change for trader; we sell at our sell-limit price, not the close price; use yesterday's close.
			#
			# sell x fraction of holdings at the time of the sell:
			dX=-stockses[-1]*sellFactor*dYtrader	# value of stock sold.
			#dX=-stockses[-1]*sellFactor				# number of shares of stock (valued at $1) sold (sort of).
			#stockToCash=dX*dYtrader						# cash value of stock when we sell.
			#dX=-stockses[-1]*sellFactor
			# and reset rbVal:
			#rbVal=tradeArrays[3][-1]/(1+sellFactor)/2	# trading induced rbVal is more sensitive; crashes beget crashes.
			oldrbVal=rbVal
			#rbVal=sellPrice*(3-sellFactor)/2			# aka, sellPrice*(1 + (1-sellFactor)/2), so that we set the effective next sell price closer than sellFactor to current price.
			rbVal=sellPrice
			nSells+=1
		# and if we didn't sell today, we can set the rbVal to today's high:
		#print "tA[4]: %d, %d" % (len(tradeArrays[4]), i)
		
		#
		# we can't go negative from commission:
	#	fixedCommission=investors[-1]*.1*commissionFactor
	#	#fixedCommission=1.0
	#	if fixedCommission > (cashes[-1]+stockses[-1]):
	#		# no money left; we can't trade.
	#		dX=0
	#		fixedCommission=0	# if we're trading with the fixed type commission
		
		# eliminate small transactions (commissions can kill us):
		if abs(dX)<.01*(cashes[-1]+stockses[-1]):
			dX=0
			fixedCommission=0	# if we're trading with the fixed type commission
			#rbVal=oldrbVal
		
		# if there was no sell yesterday, set our max-val to yesterday's high price (if it is a record value):
	#	if didSell==0 and tradeArrays[4][-1]>rbVal:
	#		# use rbVal is set to today's high (not close) value. 
	#		#print "HiVal rbVal change: %f, %f" % (rbVal, tradeArrays[4][-1])
	#		rbVal=tradeArrays[4][-1]	
		#
		cashes+=[cashes[-1]-dX]
		if dX>=0:
			#a buy at the end of the day:
			stockses+=[(stockses[-1]*dY) + dX]	# stock changes with the market and then we buy. our contribution does not see the market.
		if dX<0:
			# this one is a little bit more complicated.
			# a sell sometime during the day:
			#stockses+=[(stockses[-1]*dYtrader)+dX]
			stockses+=[dY*stockses[-1]*(1-sellFactor)]		# what's left in the portfolio declines by dY. the rest was converted to cash.
			#cashes+=[cashes[-1]+stockses[-2]*dYtrader*sellFactor	# in contrast to stockses[]; note this is the same as dX.
			#print "sell transaction (%s): r(t): %f open: %f, close: %f, sellPrice: %f, stockStart: %f, cashStart: %f, stockFinish: %f, cashFinish: %f" % (datetime.datetime.fromordinal(tradeArrays[0][-1]), tradeArrays[1][-1], tradeArrays[5][-2], tradeArrays[5][-1], sellPrice, stockses[-2], cashes[-2], stockses[-1], cashes[-1])
		
		# since we're no longer trading at the same times, this will almost always go nuts unless we get a little smarter...
		#if abs((cashes[-2]+stockses[-2])-(cashes[-1]+stockses[-1]))>.000001: print "cassh conservation error... %f" % ((cashes[-2]+stockses[-2])-(cashes[-1]+stockses[-1]))
		#

		
		# and pay a commission:
		#print "commission: %f, %f" % (investors[-1], fixedCommission)
		cashes[-1]-=abs(dX)*commissionFactor
		#cashes[-1]-=fixedCommission
		totals+=[stockses[-1]+cashes[-1]]
		djiaStart+=1
		#print "final: %f, %f/%f, %f, %f, %f gain: %f" % (tradeArrays[1][-1], stockses[-1], stockses[-1]/(1.0+dY), cashes[-1], 1.0+dY, investors[-1], (stockses[-1]-stockses[-2] + cashes[-1]-cashes[-2]))
		#if i>=500: break
		
		#i+=1
	#
	print "nBuys, nSells: %d, %d" % (nBuys, nSells)
	# now, let's look at stockses and cashes:
	fig=plt.figure(2)
	plt.clf()
	plt.plot_date(tradeArrays[0], tradeArrays[1], "b-", xdate=True, ydate=False, label="r(t)_ave")
	fig=plt.figure(3)
	plt.clf()
	plt.plot_date(tradeArrays[0], stockses[0:], "b-", xdate=True, ydate=False, label="stocks")
	plt.plot_date(tradeArrays[0], cashes[0:], "g-", xdate=True, ydate=False, label="cash")
	plt.plot_date(tradeArrays[0], investors[0:], "r-", xdate=True, ydate=False, label="DJIA")
	plt.plot_date(tradeArrays[0], totals[0:], "y-", xdate=True, ydate=False, label="trader")
	fig.autofmt_xdate()
	fig.gca().set_yscale('log')
	plt.legend(loc="upper left")
	plt.savefig("%s/rbTradeModel.pdf" % (outputFolder))
	plt.show()
	
	print "final(2): stock: %f, cash: %f, index: %f/%f, gain: %f" % (stockses[-1], cashes[-1], tradeArrays[3][-1],investors[-1], (stockses[-1]+cashes[-1])/tradeArrays[3][-1])
	#print "start: %f, %f, %f, %f" % (tradeArrays[0][0], tradeArrays[1][0], tradeArrays[2][0], tradeArrays[3][0])

def getRecordStocks(rbLen=128, mode="close"):
	# (taken from recordBreaker.py)
	# basically, copy calcFullSet()
	# start with some known data, then generalize...
	# nominally, we should pull directly from sql...
	import _mysql
	import MySQLdb
		# make a sql connection:
	sqlHost='localhost'
	sqlPort=3306
	thisdb='markets'
	#sqlUser='defaultInter'
	#sqlPW = 'Int3r@ct1Ve'
	sqlUser='myoder'
	sqlPW='yoda'
	#mycon = _mysql.connect(host=sqlHost, user=sqlUser, passwd=sqlPW, port=sqlPort, db='QuakeData')
	con1 = MySQLdb.connect(host=sqlHost, user=sqlUser, passwd=sqlPW, port=sqlPort, db=thisdb)
	#con2  = MySQLdb.connect(host=sqlHost, user=sqlUser, passwd=sqlPW, port=sqlPort, db='markets')
	c1=con1.cursor()
	#c2=con2.cursor()
	#
	isNT=True
	#rbLen=128
	rbStep=1
	#outputFolder = "econoPhys/djia10"
	outputFolder = "econoRB/djia10"
	strNT='time'
	if isNT: strNT='NT'
	djia1="select cast(dt as datetime), high/low as hilo, `close`/`open` as closopen from stocks where ticker=\'djia\'"	# nominally, close-open should be close_i - close_{i-1}
	djia2="select  cast(dt as datetime), `open`, `close` from stocks where ticker=\'djia\'"
	
	if mode=="close":
		djiaSelect=djia1
		catstring="djiaClose"
	if mode=="clopen":
		djiaSelect=djia2
		catstring="djiaClopen"
	c1.execute(djiaSelect)
	#
	rb=rbb.recordbreaker()
	#
	# this goes directly to getRecordArrays() which uses getIntervals() which returns, [[evNum], [fdate], [interval]]
	djiaCat=[[],[],[]]
	nrows=1
	print 'spin djia cursor...'
	while (1==1):	# this strange notation in case we get 0 rows...
		row = c1.fetchone()		# or use .fetchall()
		if row==None:
			break		# end of cursor
		#
		djiaCat[0]+=[nrows]
		djiaCat[1]+=[long(row[0].toordinal())]
		#djiaCat[2]+=[row[1]]
		djiaCat[2]+=[(row[2])]
		nrows+=1
	#
	'''
	plt.figure(0)
	plt.clf()
	plt.plot_date(djiaCat[1], djiaCat[2], '-')
	plt.show()
	'''
	
	print "do record-breaking..."
	djiaRecords=rb.getRecordArrays(djiaCat, isNT, rbLen, rbStep)
	print "records acquired."
	# this returns: [[biggers], [smallers], [bigMags], [smallMags], [dateVectors]]
	# the first 4 lists are records like: [[x_<=1, x_2, x_4, x_8, x_16]], or in other words the log-"bins" (not really bins because it's cumulative)
	# values might be like [[1,1,2,3,3]], which is the number of records broken at t<=1,2,4,8,16. note t can be natural time.
	#
	# now, make average and ts plots...
	#
	c1.close()
	con1.close()
	#return djiaRecords
	
	print "get poisson records"
	recsPoisson=rb.getRecordArrays(rb.getPoissonIntervals(len(djiaCat[0])), isNT, rbLen, rbStep)
	print "poisson cat len: %d" % len(recsPoisson[0])
	meanvals=[[],[],[],[], [],[],[],[]]	# aveNbigs, aveNsmalls, aveMagBig, aveMagSmall, {same for poisson cat}
	#
	meanvals[0]=rb.averageRecordCols(djiaRecords[0])
	meanvals[1]=rb.averageRecordCols(djiaRecords[1])
	meanvals[2]=rb.averageRecordCols(djiaRecords[2])
	meanvals[3]=rb.averageRecordCols(djiaRecords[3])
	#
	meanvals[4]=rb.averageRecordCols(recsPoisson[0])
	meanvals[5]=rb.averageRecordCols(recsPoisson[1])
	meanvals[6]=rb.averageRecordCols(recsPoisson[2])
	meanvals[7]=rb.averageRecordCols(recsPoisson[3])
	#
	# now, make a gnuplottable file:
	#outfname="rbdata/rb-%s-gnu.dat" % catname.split('/')[-1]
	outfname="%s/rb-%s-ints-gnu.dat" % (outputFolder, catstring)
	
	fout=open(outfname, 'w')
	fout.write("#djia %s" % catstring)
	fout.write("#record breaking summary data output\n")
	fout.write("# (catname, minmag, minLat, maxLat, minLon, maxLon, isNT, winLen, winStep)\n")
	fout.write("#calcFullSet(%s, %d, %d, %d)\n" % (catstring, isNT, rbLen, rbStep))
	fout.write("#>=log2(nDays)\t<nBig>\tsdevnBig\t<nSmall>\tsdevnSmall\t<magBig>\tsdevmagBig\t<magSmall>\tsdevmagSmall\n")
	for i in xrange(len(meanvals[0][0])):
		fout.write("\n%d\t" % i)
			# each bin-val:
		for ary in meanvals:
			if i>=len(ary[0]): continue
			fout.write("%f\t%f\t" % (ary[0][i], ary[1][i]))
			#if ary[0][i]<ary[0][i-1]: print "anomoly: %f, %f" % (ary[0][i], ary[0][i-1])
	fout.write("\n")
	#
	fout.close()
	
	# let's get some time series data as well.
	#outfname="rbdata/rb-%s-tsBig.dat" % catname.split('/')[-1]
	outfname="%s/rb-%s-ts%d.dat" % (outputFolder, catstring, rbLen)
	fout=open(outfname,'w')
	nbins=len(djiaRecords[0][0])
	fout.write("#djia %s" % catstring)
	fout.write("#record breaking time series: number of large records broken before Nmax=2**%d days\n" % nbins)
	fout.write("#NT\tfTime\tNbig\tNsmall\tpoisNbig\tpoisNsmall\n")
	#
	rbtsBig=[]
	rbtsSmall=[]
	rbTime=djiaRecords[4][1]		# these time/NT columns might be reversed.
	rbNT=djiaRecords[4][0]
	lenBig=len(djiaRecords[0][0])
	lenSmall=len(djiaRecords[1][0])	# lenBig and lenSmall should almost always have the same value.
	#print "lenPlots: %d, %d, %d" % (len(rbNT), len(recset[0]), len(recsPoisson[0]))
	for rownum in xrange(len(djiaRecords[0])):
		if len(djiaRecords[0][rownum])!=lenBig: continue	# we could break...
		# if len(recset[0][rownum])==lenBig: 
		rbtsBig+=[djiaRecords[0][rownum][-1]]
		#if len(recset[1][rownum])==lenSmall: 
		rbtsSmall+=[djiaRecords[1][rownum][-1]]
		# datafile:
		fout.write("%f\t%f" % (rbTime[rownum], rbNT[rownum]))
		fout.write("\t%d\t%d\t%d\t%d\n" % (djiaRecords[0][rownum][-1], djiaRecords[1][rownum][-1], recsPoisson[0][rownum][-1], recsPoisson[1][rownum][-1]))
		
		#for val in recset[0][rownum]:
		#	fout.write("\t%f" % val)
	fout.close()
	#
	# now, pyplot the mean values. we assume/know that the values are log_2 binned; each meanvals[i] is a 3-tuple: [[means], [sdevs], [denoms]]. we only want the first two cols.
	X=[[],[]]
	for i in xrange(len(meanvals[0][0])):
		X[0]+=[i]
		X[1]+=[2**i]
	#
	plt.figure(0)
	plt.title("Record Breaking Intervals\n(DJIA %s)" % catstring)
	plt.xlabel("Number of Events (Natural Time)")
	plt.ylabel("Number of Record Breaking Events")
	plt.errorbar(X[0],meanvals[0][0], yerr=meanvals[0][1], fmt='.-', label="DJIA Large")
	plt.errorbar(X[0],meanvals[1][0], yerr=meanvals[1][1], fmt='.-', label="DJIA Small")
	plt.errorbar(X[0],meanvals[4][0], yerr=meanvals[4][1], fmt='.-', label="Poisson Large")
	plt.errorbar(X[0],meanvals[5][0], yerr=meanvals[5][1], fmt='.-', label="Poisson Small")
	plt.xticks(X[0], X[1])
	plt.legend(loc='upper left')
	#plt.show()
	plt.savefig("%s/rb-%s-ints%s-nrbBigSmall%d.pdf" % (outputFolder, catstring, strNT, rbLen))
	plt.savefig("%s/rb-%s-ints%s-nrbBigSmall%d.png" % (outputFolder, catstring, strNT, rbLen))
	plt.clf()
	
	plt.figure(1)
	plt.title("Record Breaking Intervals\n(DJIA %s)" % catstring)
	plt.xlabel("Number of Events (Natural Time)")
	plt.ylabel("Duration of Record Breaking Interval (days)")
	plt.errorbar(X[0],meanvals[2][0], yerr=meanvals[2][1], fmt='.-', label="DJIA Large")
	plt.errorbar(X[0],meanvals[3][0], yerr=meanvals[3][1], fmt='.-', label="DJIA Small")
	plt.errorbar(X[0],meanvals[6][0], yerr=meanvals[6][1], fmt='.-', label="Poisson Large")
	plt.errorbar(X[0],meanvals[7][0], yerr=meanvals[7][1], fmt='.-', label="Poisson Small")
	plt.xticks(X[0], X[1])
	plt.legend(loc='upper left')
	#plt.show()
	plt.savefig("%s/rb-%s-ints%s-magsBigSmall%d.pdf" % (outputFolder, catstring, strNT, rbLen))
	plt.savefig("%s/rb-%s-ints%s-magsBigSmall%d.png" % (outputFolder, catstring, strNT, rbLen))
	plt.clf()


def getRecordStocksCat():
	#
	# (taken from recordBreaker.py)
	# but, this doesn't seem to work very well because it wants to calculate intervals and all that crap. back to plan A...
	
	# def calcFullSet(self, catname='cats/cmt7705.cat', minmag=5.5, minLat=-90, maxLat=90, minLon=-180, maxLon=360, isNT=True, winLen=128, winStep=1, outputFolder='paperOutput/cmt/intsNT'):
	import _mysql
	import MySQLdb
		# make a sql connection:
	sqlHost='localhost'
	sqlPort=3306
	thisdb='markets'
	#sqlUser='defaultInter'
	#sqlPW = 'Int3r@ct1Ve'
	sqlUser='myoder'
	sqlPW='yoda'
	#mycon = _mysql.connect(host=sqlHost, user=sqlUser, passwd=sqlPW, port=sqlPort, db='QuakeData')
	con1 = MySQLdb.connect(host=sqlHost, user=sqlUser, passwd=sqlPW, port=sqlPort, db=thisdb)
	#con2  = MySQLdb.connect(host=sqlHost, user=sqlUser, passwd=sqlPW, port=sqlPort, db='markets')
	c1=con1.cursor()
	#c2=con2.cursor()
	#
	#djiaSelect="select dt, adjClose, volume, high-low as hilo, `adjClose`-`open` as closopen from stocks where ticker=\'djia\'"	# nominally, close-open should be close_i - close_{i-1}
	# nominally, getting some volume, close, etc. data would be a good idea, but for the time being,
	# we have to trick the system into thinking this is an earthquake catalog.
	djiaSelect="select dt, high-low as hilo, `adjClose`-`open` as closopen from stocks where ticker=\'djia\'"
	c1.execute(djiaSelect)
	#
	rb=recordbreaker()
	#
	djiaCat=[[],[],[],[],[],[]]
	nrows=0
	djiafname="econoPhys/djia10/djia.cat"
	fout=open(djiafname, 'w')
	fout.write("#note: we always use adjusted Close values\n")
	#fout.write("#djia record-breaking catalog\n#rownum\tdate\tadjClose\tvolume\thilo\tclosopen\n")
	fout.write("#djia record-breaking catalog\n#rownum\tdate\tlat\tlon\thilo\tclosopen\n")
	print 'spin djia cursor and make a catalog......'
	while (1==1):	# this strange notation in case we get 0 rows...
		row = c1.fetchone()		# or use .fetchall()
		if row==None:
			break		# end of cursor
		nrows+=1
		#fout.write("%d\t%s\t%f\t%d\t%f\t%f\n" % (nrows, row[0], row[1], row[2], row[3], row[4]))
		fout.write("%s\t%f\t%f\t%f\t%f\n" % (("%s\t00:00:00.0" % row[0]).replace("-", "/"), 42, 42, 5.5, row[2]))	# eventually, maybe we use hi/lo as magnitude, but there's this negative sign problem...
	fout.close()
	c1.close()
	#
	print "catalog written: %s" % djiafname
	con1.close()
	#
	# now, we have a catalog. do the full set:
	rb.calcFullSet(djiafname, 0, -90, 90, -360, 360, True, 10, 1, 'econoPhys/djia10')
	# is that is?
	return 0
				

	
def getRecordStockses():
	getRecordStocks(64, "close")
	getRecordStocks(64, "clopen")

