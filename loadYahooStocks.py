#import math
import string
#import pylab
#import os
import time
import datetime as dtm
import urllib


########################################################################
# production:
#####################
def yahooStock2File(ticker='AAPL', foutname=None):
	# production function call (aka, call this function)
	if ticker==None: return None	
	if foutname==None: foutname='%s-history.dat' % ticker
	#
	#stockData=yahooStock2List(ticker)
	stockData=reformatStockList(yahooStock2List(ticker))
	
	#stockData.reverse()	# just a choice...
	fout=open(foutname, 'w')
	for rw in stockData:
		# first row will be the header.
		strRw=''
		if type(rw[0]).__name__=='str':
			# header
			#continue	# skip the header...
			strRw='#'
		for elem in rw: strRw+='%s\t' % str(elem)
		strRw=strRw[:-1]
		while strRw[-1]=='\n': strRw=strRw[:-1]
		strRw += '\n'
		fout.write(strRw)
	fout.close()
		
	return len(stockData)

#############################################
# tools, bits, and components:

def parseYahooStockFile(finName=None, foutName=None):
	# read a yahoo stock history file, reformat, output to foutname.
	# file is not properly commented, but should not contain anything besides a header and data.
	# (this is functional, but it's also a first effort, so it is no optimized to integrate with other parts; it is basically stand-alone at this point)
	fin=open(finName)
	fout=open(foutName, 'w')
	nrows=0
	#
	for rw in fin:
		strOutRw=''
		# dctRw={}	# one approach would be to make a dictionary from the header row, but we still have to know which bits we want.
		rws=rw.split(',')
		if rws[0]=='Date':
			# header row.
			strOutRw='#fdate\t'	# adding the float date col.
			for elem in rws:
				strOutRw+='%s\t' % elem
			strOutRw=strOutRw[:-1]+'\n'
			#
			fout.write(strOutRw)
			continue
		#
		# otherwise, it's data and we need to interpret the date -> float date.
		thisdtStr=rws[0]	# all other vals are floats.
		dts=thisdtStr.split('-')
		print thisdtStr
		thisDtobj=dtm.date(int(dts[0]), int(dts[1]), int(dts[2]))	# day-level resolution. we'll need a new (or imporved) script for higher time-res.
		fdt=thisDtobj.toordinal()	# this is not really a float date, but since we only have day-resolution, this will do. if we have time-resolution:
		# fdt=pylab.date2num(thisDtobj)
		strOutRw='%d\t%s\t' % (fdt, thisdtStr)
		#
		for i in xrange(1, len(rws)):
		#for elem in rws:
			elem=rws[i]
			strOutRw+='%s\t' % float(elem)
		strOutRw=strOutRw[:-1]+'\n'
		#
		fout.write(strOutRw)
		nrows+=1
	fin.close()
	fout.close()
	#
	return nrows

def yahooStock2List(ticker='AAPL'):
	# data come from: 
	# http://ichart.finance.yahoo.com/table.csv?s=AAPL&d=10&e=5&f=2010&g=d&a=8&b=7&c=1984&ignore=.csv
	# http://ichart.finance.yahoo.com/table.csv?s=USL&d=10&e=5&f=2010&g=d&a=11&b=6&c=2007&ignore=.csv
	# but i think this (for example) will suffice (probably from some server-side default settings):
	# http://ichart.finance.yahoo.com/table.csv?s=AAPL
	#
	# load the stock data to a file handler:
	# f = urllib.urlopen('http://www.ncedc.org/cgi-bin/catalog-search2.pl', urllib.urlencode(anssPrams))
	f= urllib.urlopen('http://ichart.finance.yahoo.com/table.csv?s=%s' % ticker)
	lout=[]
	nrows=0
	for rw in f:
		#print rw
		rws=rw.split(',')
		if rws[0]=='Date':
			# header row.
			#strOutRw='#fdate\t'	# adding the float date col.
			thisRw=['fdate']
			for elem in rws:
				#strOutRw+='%s\t' % elem
				thisRw+=[elem]
			thisRw[-1].replace('\n', '')
			lout+=[thisRw]
			#
			continue
		# otherwise, it's data and we need to interpret the date -> float date.
		thisdtStr=rws[0]	# make date-stuff from this; all other vals are floats.
		#print 
		dts=thisdtStr.split('-')
		#print thisdtStr
		thisDtobj=dtm.date(int(dts[0]), int(dts[1]), int(dts[2]))	# day-level resolution. we'll need a new (or imporved) script for higher time-res.
		fdt=thisDtobj.toordinal()	# this is not really a float date, but since we only have day-resolution, this will do. if we have time-resolution:
		# fdt=pylab.date2num(thisDtobj)
		#strOutRw='%d\t%s\t' % (fdt, thisdtStr)
		thisRw=[fdt, thisDtobj]	# maybe use thisdtstr, but for now, we'll let python handle the date object.
		#
		for i in xrange(1, len(rws)):
		#for elem in rws:
			thisRw+=[float(rws[i])]
		lout+=[thisRw]
		nrows+=1
	f.close()
	#
	return lout

def getYearFrac(indt):
	if type(indt).__name__ in ('int', 'float'):
		# a date-number. eventually, we'll code for datetimes, datetime-strings, etc. note that this will auto-work for a date or datetime parameter.
		indt=dtm.date.fromordinal(indt)
	#
	myttuple=indt.timetuple()
	yr=myttuple[0]
	dys=myttuple[7]
	maxdt=dtm.date(yr, 12,31)
	Ndays=maxdt.timetuple()[7]
	#
	outdate=yr+float(dys)/Ndays
	return outdate

def reformatStockList(stockList):
	# oldest to newest, no header, no adj-close, datetime all the way to the right, and float-datetime->year.fracYear (aka, 2010-06-01 ~> 2010.5)
	newList=[]
	for rw in stockList:
		if type(rw[0]).__name__=='str': continue	# skip the header
		#newRw=[]
		#rws=rw.split()
		dtnum=rw[0]
		thisdt=rw[1]
		opn=rw[2]
		hi=rw[3]
		lo=rw[4]
		cls=rw[5]
		vl=rw[6]
		adjcl=rw[7]
		yearFrac=getYearFrac(dtnum)
		#
		#newList+=[[yearFrac, opn, hi, lo, cls, vl, str(thisdt).replace('-', '/')]]
		newList+=[[yearFrac, opn, hi, lo, cls, vl, adjcl, str(thisdt).replace('-', '/')]]
		#newList+=[[yearFrac, opn, hi, lo, cls, vl, thisdt]]
	if newList[0][0]>newList[-1][0]: newList.reverse()
	return newList



