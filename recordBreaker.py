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

class recordbreaker():
	catfile=''
	catalog=[[],[],[],[],[],[]]		# base source of events: [[evnum], [date], [lat], [lon], [mag]]
	activecat=[[],[],[],[],[],[]]	# perform necessary filters, transforms, etc.; use this catalog. number/date formats may vary from [catalog].
									# format: [[evnum], [dt-float], [lat], [lon], [mag], [interval]]
	meanInterval=0
	meanIntSqr=0
	
	def __init__(self, catname='cats/cmt7705.cat', minmag=5.5, minLat=-90, maxLat=90, minLon=-180, maxLon=360):
		self.initialize(catname, minmag, minLat, maxLat, minLon, maxLon)
	
	def initialize(self, catname='cats/cmt7705.cat', minmag=5.5, minLat=-90, maxLat=90, minLon=-180, maxLon=360):
		self.loadCat(catname, minmag, minLat, maxLat, minLon, maxLon)
		self.catfile=catname
	
	def sayHello(self, mystr):
		print 'hello %s.' % mystr
	
	def getSubcat(self, start, end):
		scat=[]
		for col in self.catalog:
			scat+=[col[start:end]]
		return scat
			
	
	def loadCat(self, catname='cats/cmt7705.cat', minmag=5.5, minLat=-90, maxLat=90, minLon=-180, maxLon=360):
		# for the most part, we start with a partially parsed catalog (aka, a regional catalog). we might regularly filter on minmag.
		f=open(catname)
		rownum=0
		self.catalog=[[],[],[],[],[],[]]
		#
		for rw in f:
			if rw[0]=='#': continue
			rws=rw.split('\t')
			if len(rws)==1: print rws
			mag=float(rws[4])
			if mag<minmag: continue
			dt=rws[0]
			tm=rws[1]
			if rownum==0: thisdt=datetimeToFloat(datetimeFromStrings(dt, tm, '/'))
			#
			prevdt=thisdt
			thisdt=datetimeToFloat(datetimeFromStrings(dt, tm, '/'))
			interval=thisdt-prevdt
			#
			self.catalog[0]+=[rownum]
			self.catalog[1]+=[thisdt]
			self.catalog[2]+=[float(rws[2])]	#lat
			self.catalog[3]+=[float(rws[3])] #lon
			self.catalog[4]+=[mag]
			self.catalog[5]+=[interval]
			#
			self.meanInterval+=interval
			self.meanIntSqr+=interval*interval
			
			if rownum%1000==0: print "loadcat rownum: %d" % rownum
			rownum+=1
		f.close()
		#
		self.meanInterval/=rownum
		self.meanIntSqr/=rownum
		self.activecat=self.catalog[:]	# by default...
		#
	#
	# wherever possible, translate all dates to float values to facilitate simple arithmetic, etc.
	def getIntervals(self):
		# get intervals from catalog -> [[ndays], [date-float], [interval]]
		retAry=[self.activecat[0][:], self.activecat[1][:], self.activecat[5][:]]

		minIndex=0
		while retAry[2][minIndex]==0:
			minIndex+=1
			
		return [self.activecat[0][minIndex:], self.activecat[1][minIndex:], self.activecat[5][minIndex:]]
	def getMags(self):
		# get magnitudes from catalog -> [[ndays], [date-float], [interval]]
		#mags=[[],[],[]]]
		return [self.activecat[0], self.activecat[1], self.activecat[4]]
	
	def getPoissonIntervals(self, catLen, startDt=1, type='intervals'):
		# make a poisson-analog catalog.
		# from starting "date", generate a poisson distributed interval; insert row(s), date[i+1]=date[i]+interval.
		# return [[nth row], [date-float], {[interval] or [magnitude]}]
		intervals=[[],[],[]]
		r1=random.Random()
		
		fdt=startDt	# float-date
		for i in xrange(startDt, catLen+startDt, 1):
			interval=r1.expovariate(1.0/self.meanInterval)	# be sure this meanInterval is set correctly to the active catalog.
			#
			intervals[0]+=[i]
			intervals[1]+=[fdt]
			intervals[2]+=[interval]
			#
			fdt+=interval
		
		#print "poisson interval Len: %d" % len(intervals[0])
		return intervals
	
	# now, a NonHomogeneous Poisson Process (HPP) catalog:
	def getNHPPintervalsOmori1(self, catLen, initialRate=1, startDt=0, t0=0):
		# NHHP for omori's law, alpha=1: N(t)=A/(1+t)
		# see: http://fedc.wiwi.hu-berlin.de/xplore/tutorials/stfhtmlnode90.html
		# basically, we generate poisson numbers like:
		# F(x) = 1-exp(-L*t)
		# we draw random numbers u from U (unform distribution, 0:1)
		# u = F(x) = 1-exp(-L*t)
		# 1-u = exp(-L*t)
		# note: statistically, 1-u = u
		# u = exp(-L*t)
		# ln(u) = -L*t
		# t = -ln(u)/L
		# where t are poisson distributed "time" intervals.
		#
		# we take this a step farther: t*L = integral(L(s+v)dv)[0:t]
		# s is the starting time, t is the stopping time (maybe better written as: integral(L(t_(i-1) + t)dt)[0:t_i]
		# picking up in the middle, we have:
		# ln(u)=-integral(L(s+v)dv)[0:t]
		# now, input some function for L(x) and integerate to get f(t), solve for t; we get some:
		# t = g(u) , g(u) of course includes the original ln(u)
		#
		l0=1.0/self.meanInterval
		intervals=[[],[],[]]
		r1=random.Random()
		
		# omori: initialRate/(t0 + t)**alpha
		#
		# alpha=1; we use ln() to itegrate.
		fdt=startDt	# float-date
		for i in xrange(startDt, catLen+startDt, 1):
			interval=(t0 + fdt)*( pow(r1.random(),-1.0/initialRate) - 1.0 )
			
			#print interval
			#
			intervals[0]+=[i+1]
			intervals[1]+=[fdt]
			intervals[2]+=[interval]
			#
			fdt+=interval
		
		return intervals
		
	#def getNHPPintervalsOmori1b(self, Tmax, initialRate=1, startT=0, t0=0):
	def getNHPPintervalsOmori1b(self, Tmax=5, tao=1.0, cm=1.0):
		# NHHP for omori's law, alpha=1: N(t)=A/(1+t)
		# see: http://fedc.wiwi.hu-berlin.de/xplore/tutorials/stfhtmlnode90.html
		# basically, we generate poisson numbers like:
		# F(x) = 1-exp(-L*t)
		# we draw random numbers u from U (unform distribution, 0:1)
		# u = F(x) = 1-exp(-L*t)
		# 1-u = exp(-L*t)
		# note: statistically, 1-u = u
		# u = exp(-L*t)
		# ln(u) = -L*t
		# t = -ln(u)/L
		# where t are poisson distributed "time" intervals.
		#
		# we take this a step farther: t*L = integral(L(s+v)dv)[0:t]
		# s is the starting time, t is the stopping time (maybe better written as: integral(L(t_(i-1) + t)dt)[0:t_i]
		# picking up in the middle, we have:
		# ln(u)=-integral(L(s+v)dv)[0:t]
		# now, input some function for L(x) and integerate to get f(t), solve for t; we get some:
		# t = g(u) , g(u) of course includes the original ln(u)
		# for L(x), use some version of omori's law:
		# r(t) = r0/(t0+t)^p.
		# don likes the "modified" form (see Shcherbakov, "Scaling Properties"):
		# r(t) = (1/tao(m))*(1+t/c(m))**-p
		# one can get this by pulling the c (or t0) out of the denominator, (1/tao) -> r0/c(m)**p. note: c(m) implies some "constant" c as a function of
		# the catalog threshold magnitude m.
		#
		intervals=[[],[],[]]		#[[i], [T], [interval]]
		r1=random.Random()		
		# omori: initialRate/(t0 + t)**alpha
		#
		# alpha=1; we use ln() to itegrate.
		#T=startT	# float-date
		T=0
		#cm=float(cm)
		#tao=float(tao)
		#T=float(T)		
		#for i in xrange(startDt, catLen+startDt, 1):
		i=0
		interval=0
		while T<Tmax:
			i+=1
			#interval=(t0 + T)*( pow(r1.random(),-1.0/initialRate) - 1.0 )
			#interval = pow(r1.random(), -tao/cm)*(cm + T) - cm - T
			interval = (pow(r1.random(), -tao/cm)-1.0)*(cm + T)
			
			#if T+interval>Tmax: break
			
			#print interval
			#
			T+=interval
			if T>Tmax: break
			#if (T<=Tmax):
			intervals[0]+=[i]
			intervals[1]+=[T]
			intervals[2]+=[interval]
			#

		#intervals[0].pop()
		#intervals[1].pop()
		#intervals[2].pop()
		
		return intervals		#[[i], [T], [interval]]; T starts at zero; first T will be 0+interval_1.
		
	def getRecordArraysBack(self, dset=None, NT=True, winLen=256, winStep=1):
		# eventually, this will fetch records backwards in time/n. this should give a more correct picture for time series. plots.
		if dset==None: dset=self.getIntervals()
		# get an array of records (log2 binned; winLen should be powers of 2).
		# data-set, natura/callendar time, window length (aka, records over 1024 days), window step (days to advance after each records scan)
		# dataset: [[event count], [date], [val]] (aka, getIntervals() or getMags()
		#
		# this function doesn't care what kind of data it's looking at; it just looks for bigger/smaller magnitude/count of records.
		#
		# setup:
		nRecordsBig=[]	# number of new records. will be [[dt=1, dt=2, dt=4, dt=8, ... ,dt=1024]..[]], one row for each day in the catalog or average along the way...
		nRecordsSmall=[]
		recMagsBig=[] 
		recMagsSmall=[] # record breaking intervals
		dateCol=1
		if NT==True: dateCol=0
		dateVecs=[[],[],[]]		# time, NT (maybe switched), mag
		#
		nbins=1+long(log2(winLen))
		baseRow=[]
		for i in xrange(nbins-1):
			baseRow+=[0]
		#
		for i in xrange(winLen,len(dset[0]), winStep):
			if i>=len(dset[0]): break		# we won't finish all the bins (note: this truncates our catalog).
			biggest=dset[2][i]
			smallest=dset[2][i]
			nRecordsBig+=[[1]+baseRow[:]]
			nRecordsSmall+=[[1]+baseRow[:]]
			recMagsBig+=[[biggest]+baseRow[:]]
			recMagsSmall+=[[smallest]+baseRow[:]]
			dateVecs[0]+=[dset[0][i]]
			dateVecs[1]+=[dset[1][i]]	# return date/NTdate columns for time series plots, should they be required.
			#print "cat-row: %d" % i
			#
			ii=0
			if i%1000==0: print "date: %s" % datetime.date.fromordinal(long(dset[1][i]))
			#while dset[dateCol][i-ii]-dset[dateCol][i]<winLen and (i-ii)<(len(dset[dateCol])-1):
			while abs(dset[dateCol][i]-dset[dateCol][i-ii])<winLen and (i-ii)>0:
				# falling off the end of the catalog? (this should be redundant from above):
				if abs(dset[dateCol][i-ii]-dset[dateCol][i])>=winLen: continue	# note, this window is in Nat.Time. we might want to calc in "regular" time.				
				#
				# calculate intervalbin (bins along winLen, how many 'days' since starting event).
				if dset[dateCol][i]-dset[dateCol][i-ii]<1:
					intervalbin=0
				else:
					intervalbin=long(log2(dset[dateCol][i]-dset[dateCol][i-ii])+1)	# so intevals<1 give a negative bin (we set it to 0), 1->2: 1, 2-4: 2, 4-8: 3, 2**n-1 - 2**n: n
				#				
				if dset[2][i-ii]>biggest:
					# we've broken a big record:
					#print "new big record: %d, %f, %f" % (i-ii, biggest, nRecordsBig[-1][intervalbin])
					biggest=dset[2][i-ii]
					nRecordsBig[-1][intervalbin]+=1
					recMagsBig[-1][intervalbin]=biggest
				#
				if dset[2][i-ii]<smallest:
					#print "new small record: %d, %f, %f" % (i-ii, smallest, nRecordsSmall[-1][intervalbin])
					smallest=dset[2][i-ii]
					nRecordsSmall[-1][intervalbin]+=1
					recMagsSmall[-1][intervalbin]=smallest
				#
				ii+=1
				#if (i-ii)>len(dset[dateCol]): break
				#print i-ii
			#	
			# now, integrate the new row:
			for irow in xrange(1,len(nRecordsBig[-1])):
				# is this record complete? 
				if i + (2**irow)>len(dset[dateCol]):
					# if not, pop out this element. subsequent bins will also delete. alternatively, set to None
					# note that this will work, but it is a sort of stupid way to shorten these rows.
					nRecordsBig[-1].pop()
					nRecordsSmall[-1].pop()
					recMagsBig[-1].pop()
					recMagsSmall[-1].pop()
					continue
				#
				# if the record is complete:		
				nRecordsBig[-1][irow]+=nRecordsBig[-1][irow-1]
				nRecordsSmall[-1][irow]+=nRecordsSmall[-1][irow-1]
				# record intervals:
				if recMagsBig[-1][irow]==0: recMagsBig[-1][irow]=recMagsBig[-1][irow-1]
				if recMagsSmall[-1][irow]==0: recMagsSmall[-1][irow]=recMagsSmall[-1][irow-1]
				#
		print "records calculated."		
		return [nRecordsBig, nRecordsSmall, recMagsBig, recMagsSmall, dateVecs]
		
		
	def getRecordArrays(self, dset=None, NT=True, winLen=256, winStep=1):
		if dset==None: dset=self.getIntervals()
		# get an array of records (log2 binned; winLen should be powers of 2).
		# data-set, natura/callendar time, window length (aka, records over 1024 days), window step (days to advance after each records scan)
		# dataset: [[event count], [date], [val]] (aka, getIntervals() or getMags()
		#
		# this function doesn't care what kind of data it's looking at; it just looks for bigger/smaller magnitude/count of records.
		#
		# setup:
		nRecordsBig=[]	# number of new records. will be [[dt=1, dt=2, dt=4, dt=8, ... ,dt=1024]..[]], one row for each day in the catalog or average along the way...
		nRecordsSmall=[]
		recMagsBig=[] 
		recMagsSmall=[] # record breaking intervals
		dateCol=1
		if NT==True: dateCol=0
		dateVecs=[[],[], []]		# time, NT (maybe switched), mag
		#
		nbins=1+long(log2(winLen))	# this may need to be 2+... to facilitate non-log2-integer window lengths...
		baseRow=[]
		for i in xrange(nbins-1):
			baseRow+=[0]
		#
		# spin through each row of the data-set:
		for i in xrange(0,len(dset[0]), winStep):
			if i>len(dset[0])-winLen: break		# we won't finish all the bins (note: this truncates our catalog).
			#
			biggest=dset[2][i]
			smallest=dset[2][i]
			nRecordsBig+=[[1]+baseRow[:]]
			nRecordsSmall+=[[1]+baseRow[:]]
			recMagsBig+=[[biggest]+baseRow[:]]
			recMagsSmall+=[[smallest]+baseRow[:]]
			dateVecs[0]+=[dset[0][i+winLen-1]]
			dateVecs[1]+=[dset[1][i+winLen-1]]	# return date/NTdate columns for time series plots, should they be required.
			#print "cat-row: %d" % i
			#
			ii=0
			if i%1000==0: print "date: %s" % datetime.date.fromordinal(long(dset[1][i]))
			while dset[dateCol][i+ii]-dset[dateCol][i]<winLen and (i+ii)<(len(dset[dateCol])-1):
				# falling off the end of the catalog? (this should be redundant from above):
				if dset[dateCol][i+ii]-dset[dateCol][i]>=winLen: continue	# note, this window is in Nat.Time. we might want to calc in "regular" time.				
				#
				# calculate intervalbin (bins along winLen, how many 'days' since starting event).
				if dset[dateCol][i+ii]-dset[dateCol][i]<1:
					intervalbin=0
				else:
					intervalbin=long(log2(dset[dateCol][i+ii]-dset[dateCol][i])+1)	# so intevals<1 give a negative bin (we set it to 0), 1->2: 1, 2-4: 2, 4-8: 3, 2**n-1 - 2**n: n
				#				
				if dset[2][i+ii]>biggest:
					# we've broken a big record:
					#print "new big record: %d, %f, %f" % (i+ii, biggest, nRecordsBig[-1][intervalbin])
					biggest=dset[2][i+ii]
					nRecordsBig[-1][intervalbin]+=1
					recMagsBig[-1][intervalbin]=biggest
				#
				if dset[2][i+ii]<smallest:
					#print "new small record: %d, %f, %f" % (i+ii, smallest, nRecordsSmall[-1][intervalbin])
					smallest=dset[2][i+ii]
					#print "intervalbin: %d/%d" % (intervalbin, len(nRecordsSmall[-1]))
					nRecordsSmall[-1][intervalbin]+=1
					recMagsSmall[-1][intervalbin]=smallest
				#
				ii+=1
				#if (i+ii)>len(dset[dateCol]): break
				#print i+ii
			#
			# now, integrate the new row:
			for irow in xrange(1,len(nRecordsBig[-1])):
				# is this record complete? 
				if i + (2**irow)>len(dset[dateCol]):
					# if not, pop out this element. subsequent bins will also delete. alternatively, set to None
					# note that this will work, but it is a sort of stupid way to shorten these rows.
					nRecordsBig[-1].pop()
					nRecordsSmall[-1].pop()
					recMagsBig[-1].pop()
					recMagsSmall[-1].pop()
					continue
				#
				# if the record is complete:		
				nRecordsBig[-1][irow]+=nRecordsBig[-1][irow-1]
				nRecordsSmall[-1][irow]+=nRecordsSmall[-1][irow-1]
				# record intervals:
				if recMagsBig[-1][irow]==0: recMagsBig[-1][irow]=recMagsBig[-1][irow-1]
				if recMagsSmall[-1][irow]==0: recMagsSmall[-1][irow]=recMagsSmall[-1][irow-1]
				#
		print "records calculated."		
		return [nRecordsBig, nRecordsSmall, recMagsBig, recMagsSmall, dateVecs]
					
	def averageRecordCols(self, records):
		# takes nRecordsBig, recMagsBig, etc.
		# aka, array like [[bin values], [],..,[]]. each set of binned values might not be the same length. for the time being, assume
		# they all start from the same value and that the first row is the longest (only rows toward the end of the set will be shorter).
		#
		# build output arrays:
		aveVals=[]
		sdevs=[]
		denoms=[]
		for x in records[0]:
			aveVals+=[0]
			sdevs+=[0]
			denoms+=[0]
		#
		nrows=0
		for row in records:
			for i in xrange(len(row)):
				aveVals[i]+=row[i]
				sdevs[i]+=row[i]*row[i]
				denoms[i]+=1
			nrows+=1
		#
		for i in xrange(len(aveVals)):
			fdenom=float(denoms[i])
			aveVals[i]/=fdenom
			if fdenom<=1.0:
				sdevs[i]=(((sdevs[i]/fdenom)-aveVals[i]**2.0))**.5
			else:
				sdevs[i]=(fdenom/(fdenom-1.0)*((sdevs[i]/fdenom)-aveVals[i]**2.0))**.5	# note: estimating the population stdev: sigP=(N/(N-1))*sigSample
			#print "sdevs[i]: %f, %f, %f" % (sdevs[i], denoms[i], float(denoms[i])/(denoms[i]-1) )
		#
		return [aveVals, sdevs, denoms]
	
	def calcFullSet(self, catname='cats/cmt7705.cat', minmag=5.5, minLat=-90, maxLat=90, minLon=-180, maxLon=360, isNT=True, winLen=128, winStep=1, outputFolder='paperOutput/cmt/intsNT'):
		if outputFolder[-1]=='/': outputFolder=outputFolder[:-1]
		catstring=catname.split('/')[-1]
		catstring=catstring.split('.')[0]
		strNT='time'
		if isNT: strNT='NT'
		# self.getRecordArrays(self, dset=None, NT=True, winLen=256, winStep=1)
		# self.loadCat(catname, minmag, minLat, maxLat, minLon, maxLon)
		#
		# use class functions to compute record breaking magnitudes and intervals for the given catalog.
		# write a gnuplot ready file and return a pyplot ready array
		# for now, everything is in terms of NaturalTime.
		#
		self.loadCat(catname, minmag, minLat, maxLat, minLon, maxLon)
		print "catlen: %d, %d" % (len(self.catalog[0]), len(self.activecat[0]))
		print "get catalog records"
		recset=self.getRecordArrays(self.getIntervals(), isNT, winLen, winStep)
		print "get poisson records"
		recsPoisson=self.getRecordArrays(self.getPoissonIntervals(len(self.activecat[0])), isNT, winLen, winStep)
		print "poisson cat len: %d" % len(recsPoisson[0])
		meanvals=[[],[],[],[], [],[],[],[]]	# aveNbigs, aveNsmalls, aveMagBig, aveMagSmall, {same for poisson cat}
		#
		meanvals[0]=self.averageRecordCols(recset[0])
		meanvals[1]=self.averageRecordCols(recset[1])
		meanvals[2]=self.averageRecordCols(recset[2])
		meanvals[3]=self.averageRecordCols(recset[3])
		#
		meanvals[4]=self.averageRecordCols(recsPoisson[0])
		meanvals[5]=self.averageRecordCols(recsPoisson[1])
		meanvals[6]=self.averageRecordCols(recsPoisson[2])
		meanvals[7]=self.averageRecordCols(recsPoisson[3])
		#
		# now, make a gnuplottable file:
		#outfname="rbdata/rb-%s-gnu.dat" % catname.split('/')[-1]
		outfname="%s/rb-%s-ints-gnu.dat" % (outputFolder, catstring)
		
		fout=open(outfname, 'w')
		fout.write("#record breaking summary data output\n")
		fout.write("# (catname, minmag, minLat, maxLat, minLon, maxLon, isNT, winLen, winStep)\n")
		fout.write("#calcFullSet(%s, %f, %d, %d, %d, %d, %s, %d, %d)\n" % (catname, minmag, minLat, maxLat, minLon, maxLon, isNT, winLen, winStep))
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
		
		#
		#
		# let's get some time series data as well.
		#outfname="rbdata/rb-%s-tsBig.dat" % catname.split('/')[-1]
		outfname="%s/rb-%s-tsBig.dat" % (outputFolder, catstring)
		fout=open(outfname,'w')
		nbins=len(recset[0][0])
		fout.write("#record breaking time series: number of large records broken before Nmax=2**%d days\n" % nbins)
		fout.write("#NT\tfTime\tNbig\tNsmall\tpoisNbig\tpoisNsmall\n")
		#
		rbtsBig=[]
		rbtsSmall=[]
		rbTime=recset[4][1]		# these time/NT columns might be reversed.
		rbNT=recset[4][0]
		lenBig=len(recset[0][0])
		lenSmall=len(recset[1][0])	# lenBig and lenSmall should almost always have the same value.
		#print "lenPlots: %d, %d, %d" % (len(rbNT), len(recset[0]), len(recsPoisson[0]))
		for rownum in xrange(len(recset[0])):
			if len(recset[0][rownum])!=lenBig: continue	# we could break...
			# if len(recset[0][rownum])==lenBig: 
			rbtsBig+=[recset[0][rownum][-1]]
			#if len(recset[1][rownum])==lenSmall: 
			rbtsSmall+=[recset[1][rownum][-1]]
			# datafile:
			fout.write("%f\t%f" % (rbTime[rownum], rbNT[rownum]))
			fout.write("\t%d\t%d\t%d\t%d\n" % (recset[0][rownum][-1], recset[1][rownum][-1], recsPoisson[0][rownum][-1], recsPoisson[1][rownum][-1]))
			
			#for val in recset[0][rownum]:
			#	fout.write("\t%f" % val)
		fout.close()
		
		# now, pyplot the mean values. we assume/know that the values are log_2 binned; each meanvals[i] is a 3-tuple: [[means], [sdevs], [denoms]]. we only want the first two cols.
		X=[[],[]]
		for i in xrange(len(meanvals[0][0])):
			X[0]+=[i]
			X[1]+=[2**i]
		#
		f0=plt.figure(0)
		#ax0=f0.gci()
		#axo.set_yscale('log')
		f0.set_yscale('log')
		plt.title("Record Breaking Intervals\n(Global CMT, 1976-2005, m>=5.5")
		plt.xlabel("Number of Events (Natural Time)")
		plt.ylabel("Number of Record Breaking Events")
		plt.errorbar(X[0],meanvals[0][0], yerr=meanvals[0][1], fmt='.-', label="CMT Large")
		plt.errorbar(X[0],meanvals[1][0], yerr=meanvals[1][1], fmt='.-', label="CMT Small")
		plt.errorbar(X[0],meanvals[4][0], yerr=meanvals[4][1], fmt='.-', label="Poisson Large")
		plt.errorbar(X[0],meanvals[5][0], yerr=meanvals[5][1], fmt='.-', label="Poisson Small")
		plt.xticks(X[0], X[1])
		plt.legend(loc='lower right')
		
		#plt.show()
		plt.savefig("%s/rb-%s-ints%s-nrbBigSmall.pdf" % (outputFolder, catstring, strNT))
		plt.clf()
		
		plt.figure(1)
		plt.title("Record Breaking Intervals\n(Global CMT, 1976-2005, m>=5.5")
		plt.xlabel("Number of Events (Natural Time)")
		plt.ylabel("Duration of Record Breaking Interval (days)")
		plt.errorbar(X[0],meanvals[2][0], yerr=meanvals[2][1], fmt='.-', label="CMT Large")
		plt.errorbar(X[0],meanvals[3][0], yerr=meanvals[3][1], fmt='.-', label="CMT Small")
		plt.errorbar(X[0],meanvals[6][0], yerr=meanvals[6][1], fmt='.-', label="Poisson Large")
		plt.errorbar(X[0],meanvals[7][0], yerr=meanvals[7][1], fmt='.-', label="Poisson Small")
		plt.xticks(X[0], X[1])
		plt.legend(loc='upper left')
		#plt.show()
		plt.savefig("%s/rb-%s-ints%s-magsBigSmall.pdf" % (outputFolder, catstring, strNT))
		plt.clf()
		
		######
		# and a back-looking time series:
		recset=self.getRecordArraysBack(self.getIntervals(), isNT, winLen, winStep)
		#outfname="rbdata/rb-%s-tsBig.dat" % catname.split('/')[-1]
		outfname="%s/rb-%s-tsBigBack.dat" % (outputFolder, catstring)
		fout=open(outfname,'w')
		nbins=len(recset[0][0])
		fout.write("#back-looking record breaking time series: number of large records broken before Nmax=2**%d days\n" % nbins)
		fout.write("#NT\tfTime\tNbig\tNsmall\tpoisNbig\tpoisNsmall\n")
		#
		rbtsBig=[]
		rbtsSmall=[]
		rbTime=recset[4][1]		# these time/NT columns might be reversed.
		rbNT=recset[4][0]
		lenBig=len(recset[0][0])
		lenSmall=len(recset[1][0])	# lenBig and lenSmall should almost always have the same value.
		#print "lenPlots: %d, %d, %d" % (len(rbNT), len(recset[0]), len(recsPoisson[0]))
		for rownum in xrange(len(recset[0])):
			if len(recset[0][rownum])!=lenBig: continue	# we could break...
			# if len(recset[0][rownum])==lenBig: 
			rbtsBig+=[recset[0][rownum][-1]]
			#if len(recset[1][rownum])==lenSmall: 
			rbtsSmall+=[recset[1][rownum][-1]]
			# datafile:
			fout.write("%f\t%f" % (rbTime[rownum], rbNT[rownum]))
			fout.write("\t%d\t%d\t%d\t%d\n" % (recset[0][rownum][-1], recset[1][rownum][-1], recsPoisson[0][rownum][-1], recsPoisson[1][rownum][-1]))
			
			#for val in recset[0][rownum]:
			#	fout.write("\t%f" % val)
		fout.close()
		
		
		return meanvals


	def calcFullSetMags(self, catname='cats/cmt7705.cat', minmag=5.5, minLat=-90, maxLat=90, minLon=-180, maxLon=360, isNT=True, winLen=128, winStep=1, outputFolder='paperOutput/cmt/magsNT'):
		if outputFolder[-1]=='/': outputFolder=outputFolder[:-1]
		catstring=catname.split('/')[-1]
		catstring=catstring.split('.')[0]
		strNT='time'
		if isNT: strNT='NT'
		# self.getRecordArrays(self, dset=None, NT=True, winLen=256, winStep=1)
		# self.loadCat(catname, minmag, minLat, maxLat, minLon, maxLon)
		#
		# use class functions to compute record breaking magnitudes and intervals for the given catalog.
		# write a gnuplot ready file and return a pyplot ready array
		# for now, everything is in terms of NaturalTime.
		#
		self.loadCat(catname, minmag, minLat, maxLat, minLon, maxLon)
		print "catlen: %d, %d" % (len(self.catalog[0]), len(self.activecat[0]))
		print "get catalog records"
		recset=self.getRecordArrays(self.getMags(), isNT, winLen, winStep)

		meanvals=[[],[],[],[]]	# aveNbigs, aveNsmalls, aveMagBig, aveMagSmall, {same for poisson cat}
		#
		meanvals[0]=self.averageRecordCols(recset[0])
		meanvals[1]=self.averageRecordCols(recset[1])
		meanvals[2]=self.averageRecordCols(recset[2])
		meanvals[3]=self.averageRecordCols(recset[3])
		#
		# now, make a gnuplottable file:
		#outfname="rbdata/rb-%s-mags-gnu.dat" % catname.split('/')[-1]
		outfname="%s/rb-%s-mags-gnu.dat" % (outputFolder, catstring)
		#
		fout=open(outfname, 'w')
		fout.write("#record breaking magnitudes summary data output\n")
		fout.write("# (catname, minmag, minLat, maxLat, minLon, maxLon, isNT, winLen, winStep)\n")
		fout.write("#calcFullSet(%s, %f, %d, %d, %d, %d, %s, %d, %d)" % (catname, minmag, minLat, maxLat, minLon, maxLon, isNT, winLen, winStep))
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
		#
		# now, pyplot the mean values. we assume/know that the values are log_2 binned; each meanvals[i] is a 3-tuple: [[means], [sdevs], [denoms]]. we only want the first two cols.
		X=[[],[]]
		for i in xrange(len(meanvals[0][0])):
			X[0]+=[i]
			X[1]+=[2**i]
		#
		plt.figure(0)
		plt.clf()
		plt.title("Record Breaking Magnitudes\n(Global CMT, 1976-2005, m>=5.5")
		plt.xlabel("Number of Events (Natural Time)")
		plt.ylabel("Number of Record Breaking Events")
		plt.errorbar(X[0],meanvals[0][0], yerr=meanvals[0][1], fmt='.-', label="Large mag")
		plt.errorbar(X[0],meanvals[1][0], yerr=meanvals[1][1], fmt='.-', label="Small mag")
		plt.xticks(X[0], X[1])
		plt.legend(loc='lower right')
		#plt.show()
		plt.savefig("%s/rb-%s-mags%s-nrbBigSmall.pdf" % (outputFolder, catstring, strNT))
		plt.clf()
		
		plt.figure(1)
		plt.clf()
		plt.title("Record Breaking Magnitudes\n(Global CMT, 1976-2005, m>=5.5")
		plt.xlabel("Number of Events (Natural Time)")
		plt.ylabel("Magnitude of Record Breaking Event")
		plt.errorbar(X[0],meanvals[2][0], yerr=meanvals[2][1], fmt='.-', label="Large mag")
		plt.errorbar(X[0],meanvals[3][0], yerr=meanvals[3][1], fmt='.-', label="Small mag")
		plt.xticks(X[0], X[1])
		plt.legend(loc='upper left')
		#plt.show()
		plt.savefig("%s/rb-%s-mags%s-magsBigSmall.pdf" % (outputFolder, catstring, strNT))
		plt.clf()
		
		return meanvals
			
#################
# some other functions...

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

def poiss(x,k):
	return ((x**k)*math.exp(-x))/yfact(k)

def runningmeanNT(dfile, xcol, ycol, width):
	# return a running mean of y, presumably a (natural)time series.
	# for now, assume we're using natural time...
	# open file:
	print "running means..."
	f=open(dfile)
	
	irow=0
	retvals=[[],[], []]	# ready to pyplot... ( (natural)date, <val>, stdev )
	rawcat=[[], []]	# array of raw values. we'll use this to subtract old vals off running mean.
	rmean=0
	rmeansqr=0
	for rw in f:
		if rw=='\n' or rw[0]=='#': continue
		#
		irow+=1
		rws=rw.split('\t')
		rmean+=float(rws[ycol])
		rmeansqr+=float(rws[ycol])**2
		rawcat[0]+=[float(rws[xcol])]
		rawcat[1]+=[float(rws[ycol])]
		#
		# buffer starting end of the set
		if irow<width: continue
		#
		# otherwise, we're adding onto one end of the average and subtracting off the other.
		retvals[0]+=[float(rws[xcol])]
		retvals[1]+=[rmean/float(width)]
		retvals[2]+=[( (width/(width-1.0))*((rmeansqr/float(width))-(rmean/float(width))**2))**.5]		# sdevs[i]=(fdenom/(fdenom-1.0)*((sdevs[i]/fdenom)-aveVals[i]**2.0))**.5
		rmean-=rawcat[1][irow-width]
		rmeansqr-=rawcat[1][irow-width]**2
		#
	#
	#plt.figure(0)
	#plt.plot(retvals[0], retvals[1])
	##plt.errorbar(retvals[0], retvals[1], yerr=retvals[2], fmt='.-')
	#plt.show()
	#
	return retvals

def plotMeanRatios(dfile,width,rb=None, fignum=0, doShow=True, bigeq=5.8, timecol=1, forwardBack=0):
	# note: SoB means "Small over Big" (and conversely BoS -> Big over Small), not it's more contemporary implication.
	#
	# forwardBack = {0: forward, 1: backward}
	rms=meanRatios(dfile, width)	# width is averaging-window, right?
	# find the rb-window length:
	# it's nested stupidly in (some of) the comments like "blah, blah, blah, Nmax=2**7 days" or something.
	f=open(dfile)
	rbwinlen=0
	for line in f:
		if "#" in line and "Nmax=" in line:
			# use "split()" to get the length...
			split1=line.split("Nmax=")[1]
			split1=split1.split("**")
			vbase=split1[0]
			vexpon=split1[1].split(" ")[0]
			
			rbwinlen=float(vbase)**(float(vexpon)-1)
			#if forwardBack==1: rbwinlen=-rbwinlen		# also, for backward RB, plot small/big, not big/small... or maybe we reverse that and flip the colors too...
			
	#
	# make axis ticks:
	minFdt=rms[0][0]
	maxFdt=rms[0][-1]
	minDt=floatToDateTime(minFdt)
	maxDt=floatToDateTime(maxFdt)
	minyr=minDt.year
	maxyr=maxDt.year
	maxBoS=0
	maxSoB=0
	#
	# timecol:
	#timecol=0	# time-time
	#timecol=1 	# natural time
	#
	while minyr%5!=0:
		minyr+=1
	while maxyr%5!=0:
		maxyr+=1

	#
	# our tick labels:
	#print "range: %d, %d" % (minyr, maxyr)
	xyrs=range(minyr,maxyr+1,5)
	print xyrs
	# the ticks must be float values:
	xtks=[]
	for elem in xyrs:
		xtks+=[datetime.date(elem,1,1).toordinal()]
	#
	# tick marks for natural time showing dates:
	# how do we want this to look? messy...
			
	#
	# now, make ratios. we'll plot the above/below mean separately, so construct
	# separate lists for each:
	BoSup=[]			# big/small values > (threshval = {some mean value, 1, etc.} (aka, color1 above, color2 below.
	BoSdown=[]		#
	meanlist=[]
	SoBup=[]			# same, but Small/Big
	SoBdown=[]
	meanlistSoB=[]
	#
	BoSmean=0
	BoSmeanSqr=0
	SoBmean=0
	SoBmeanSqr=0
	
	#
	for i in xrange(len(rms[0])):
		#BoS+=[rms[2][i]/rms[3][i]]
		#SoB+=[1/BoS[-1]]
		#
		BoSmean+=rms[2][i]/rms[3][i]
		BoSmeanSqr+=(rms[2][i]/rms[3][i])**2
		SoBmean+=rms[3][i]/rms[2][i]
		SoBmeanSqr+=(rms[3][i]/rms[2][i])**2
	Nev=float(len(rms[0]))
	BoSmean/=Nev
	SoBmean/=Nev
	SdevFactor=1		# population stdev factor
	if Nev>1: SdevFactor=(Nev/(Nev-1))
	varBoS=SdevFactor*((BoSmeanSqr/Nev)-BoSmean)
	varSoB=SdevFactor*((SoBmeanSqr/Nev)-SoBmean)
	print "mean values: %f, %f, %f, %f" % (BoSmean, varBoS, SoBmean, varSoB)
	#
	threshval=1.0
	BoSmean=threshval
	SoBmean=threshval
	for i in xrange(len(rms[0])):
		# assign values>mean to the "up" plot, values<mean to the "down" plot
		BoSval=float(rms[2][i])/rms[3][i]
		SoBval=1/BoSval
		# establish chart limits, default values.
		if BoSval>maxBoS: maxBoS=BoSval
		if SoBval>maxSoB: maxSoB=SoBval
		#
		# by default, each list gets the mean value.
		valup=BoSmean
		valdown=BoSmean
		valupSoB=SoBmean
		valdownSoB=SoBmean
		#
		# if the value at hand is greter/less than the threshold, put it in the upper/lower list for plotting. note, for values above threshld
		# we plot the threshold value for the "down" plot and vice versa for the "up" plot.
		# then, we just keep a list of the mean (threshold) values to draw a nice little divider along that value.
		if BoSval>BoSmean: valup=BoSval
		if BoSval<BoSmean: valdown=BoSval
		BoSup+=[valup]
		BoSdown+=[valdown]
		meanlist+=[BoSmean]
		#
		if SoBval>SoBmean: valupSoB=SoBval
		if SoBval<SoBmean: valdownSoB=SoBval
		SoBup+=[valupSoB]
		SoBdown+=[valdownSoB]
		meanlistSoB+=[BoSmean]
	#
	# finally, plot:
	plt.figure(fignum)
	plt.clf()
	plt.cla()
	#
	# default:
	lstUp=BoSup
	lstDn=BoSdown
	lstMean=meanlist
	meanval=BoSmean
	maxval=maxBoS
	mainLabel='Nbig/Nsmall'
	#
	# or, if we're doing SoB plots (say we're going backwards):
	if forwardBack==1:
		print "doing backwards protocols..."
		lstUp=SoBup
		lstDn=SoBdown
		#for i in xrange(len(lstUp)): lstUp[i]=1/lstUp[i]
		#for i in xrange(len(lstDn)): lstDn[i]=1/lstDn[i]
		lstMean=meanlistSoB
		meanval=SoBmean
		maxval=maxSoB
		mainLabel='Nsmall/Nbig'
	'''
	SoBup=[]
	SoBdown=[]
	for i in xrange(len(BoSup)):
		SoBup+=[1/BoSup[i]]
		SoBdown+=[1/BoSup[i]]
	lstUp=SoBup
	lstDn=SoBdown
	'''
	#
	'''
	lstUp=SoBup
	lstDn=SoBdown
	lstMean=meanlistSoB
	meanval=SoBmean
	maxval=maxSoB
	'''
	#
	'''
	# make an error envelope:
	erPlus=[[], []]
	erMinus=[[], []]
	for i in xrange(len(rms[0])):
		erPlus[0]+=[rms[timecol][i]]
		erMinus[0]+=[rms[timecol][i]]
		erPlus[1]+=[(rms[2][i]/width)+rms[4][i]]
		erMinus[1]+=[(rms[2][i]/width)-rms[4][i]]
	plt.semilogy(erPlus[0], erPlus[1], '-')
	plt.semilogy(erMinus[0], erMinus[1], '-')
	'''
	#
	#plt.semilogy(rms[timecol], BoSup, label='BoS')
	#plt.semilogy(rms[timecol], BoSdown, label=None)
	plt.semilogy(rms[timecol], lstUp, '-b', label=mainLabel)
	plt.semilogy(rms[timecol], lstDn, '-r', label=None)

	#
	# fill:
	#plt.fill([rms[timecol][0]]+rms[timecol]+[rms[timecol][0]], [BoSmean]+BoSup+[BoSmean], 'b')
	#plt.fill([rms[timecol][0]]+rms[timecol]+[rms[timecol][-1]], [BoSmean]+BoSdown+[BoSmean], 'g')
	plt.fill([rms[timecol][0]]+rms[timecol]+[rms[timecol][-1]], [meanval]+lstUp+[meanval], 'b')
	plt.fill([rms[timecol][0]]+rms[timecol]+[rms[timecol][-1]], [meanval]+lstDn+[meanval], 'r')
	
	# horizontal line(s):
	#plt.semilogy(rms[timecol], getNlist(BoSmean, len(rms[timecol])), label='mean=%f' % BoSmean)
	
	#plt.semilogy([rms[timecol][0], rms[timecol][-1]], [1,1], label='mean=%f' % meanval)
	plt.semilogy([rms[timecol][0], rms[timecol][-1]], [1,1], label=None)
	#
	#plt.semilogy(rms[timecol], SoB, label='SoB')
	# get catalog events?
	if rb!=None:
		mags=[[],[],[], []]
		for i in xrange(len(rb.activecat[4])):
			if rb.activecat[4][i]>=bigeq:
				# activecat is [[nt], [tim], [lat], [lon], [mag], [interval]]
				mags[0]+=[rb.activecat[1][i]]	# time time
				mags[1]+=[rb.activecat[0][i]]	# natural time
				mags[2]+=[(rb.activecat[4][i]/bigeq)*maxval*.5]	# magnitude (scaled for presentation)
				mags[3]+=[rb.activecat[0][i]+rbwinlen]	# natural time
		#plt.plot(mags[0], mags[2], 'o')
		plt.plot(mags[timecol], mags[2], 'o-', label='events m>=%s' % str(bigeq))
		if timecol==1:	# it's natural time; plot the event-shadow (event position + winLen)
			plt.plot(mags[3], mags[2], 'd', label="winlen=%d" % abs(rbwinlen))
		#print "stuff: %s, %s" % (maxBoS, bigeq)
		plt.semilogy([rms[timecol][0], rms[timecol][-1]], [6.0*maxval/(bigeq*2),6.0*maxval/(bigeq*2)], label='mag=%s' % "6.0")	# (6.0*maxBoS/(bigeq*2)
		#plt.vlines(mags[timecol], 1/maxval, maxval)
	#
	# format for time. we could achieve this by checking for a datetime.date or datetime.datetime type, except of course
	# that we're expecting a float-date-time. so let's assume that float means datetime; integers mean natural time.
	if timecol==0:
		#plt.xticks(xtks, xyrs)
		axis=pyplot.gca()
		fig = pyplot.gcf()
		#axis.xaxis.set_major_formatter(dates.DateFormatter('%Y-%b-%d %H:%M'))
		axis.xaxis.set_major_formatter(dates.DateFormatter('%Y-%b-%d'))
		fig.autofmt_xdate() 
	
	plt.legend(loc='lower left')
	if doShow: plt.show()

def getNlist(val, llen):
	# return a list of llen val's:
	rlist=[]
	for i in xrange(llen):
		rlist+=[val]
	return rlist

def meanRatios(dfile, width):
	# for now, assume a standard time-series datafile:
	# #fTime	NT	Nbig	Nsmall	poisNbig	poisNsmall		# (NT and fTime (float-time) might be switched)
	# so reconstruct the "get running mean bit" get:
	# n, t, x1, x2, sig1, sig2.
	f=open(dfile)
	
	irow=0
	retvals=[[],[],[], [],[],[]]		# n,t, Bigs, smalls, SigBig, SigSmall
	rawcat=[[],[],[],[]]
	rmean1=0
	rmean2=0
	rmeansqr1=0
	rmeansqr2=0
	#
	for rw in f:
		if rw=='\n' or rw[0]=='#': continue
		#
		irow+=1
		rws=rw.split('\t')
		rmean1+=float(rws[2])
		rmean2+=float(rws[3])
		#
		rmeansqr1+=float(rws[2])**2
		rmeansqr2+=float(rws[3])**2
		#
		rawcat[0]+=[float(rws[0])]
		rawcat[1]+=[float(rws[1])]
		rawcat[2]+=[float(rws[2])]
		rawcat[3]+=[float(rws[3])]

		#
		# buffer starting end of the set
		if irow<width: continue
		#
		# otherwise, we're adding onto one end of the average and subtracting off the other.
		retvals[0]+=[float(rws[0])]	#float-date
		retvals[1]+=[float(rws[1])]	# nat-time
		#
		retvals[2]+=[rmean1/float(width)]	#NRB_big
		retvals[3]+=[rmean2/float(width)]	#NRB_small
		#
		# note: we use the factor width/(width-1.0) to get population stdev. if width==1: use sample stdev.
		if width>1.0:
			retvals[4]+=[( (width/(width-1.0))*((rmeansqr1/float(width))-(rmean1/float(width))**2))**.5]		# sdevs[i]=(fdenom/(fdenom-1.0)*((sdevs[i]/fdenom)-aveVals[i]**2.0))**.5
			retvals[5]+=[( (width/(width-1.0))*((rmeansqr2/float(width))-(rmean2/float(width))**2))**.5]		# sdevs[i]=(fdenom/(fdenom-1.0)*((sdevs[i]/fdenom)-aveVals[i]**2.0))**.5
		else:
			retvals[4]+=[( ((rmeansqr1/float(width))-(rmean1/float(width))**2))**.5]		# sdevs[i]=(fdenom/(fdenom-1.0)*((sdevs[i]/fdenom)-aveVals[i]**2.0))**.5
			retvals[5]+=[( ((rmeansqr2/float(width))-(rmean2/float(width))**2))**.5]		# sdevs[i]=(fdenom/(fdenom-1.0)*((sdevs[i]/fdenom)-aveVals[i]**2.0))**.5
		#		
		rmean1-=rawcat[2][irow-width]
		rmean2-=rawcat[3][irow-width]
		rmeansqr1-=rawcat[2][irow-width]**2
		rmeansqr2-=rawcat[3][irow-width]**2
		#
	#
	return retvals
	

def rmeanSet(dfile, xcol, ycol, widths, rb=None, outfile="rbdata/meansdata.dat"):
	# get a set of running means averaged over a range of "widths" (aka, over 1, 16, 64, 128 events).
	# plot each smoothed series together.
	plt.figure(0)
	#for width in widths:
	#	thismeans=runningmeanNT(dfile, xcol, ycol, width)
	#	plt.plot(thismeans[0], thismeans[1], label="AvWidth=%d"%width)
	
	thismeans=runningmeanNT(dfile,xcol, ycol, widths[-1])
	plt.plot(thismeans[0], thismeans[1], label="bigs (%d)" % widths[-1])

	thismeans=runningmeanNT(dfile,xcol, ycol+1, widths[-1])
	plt.plot(thismeans[0], thismeans[1], label="smalls")
	# get large earthquakes:
	if rb!=None:
		bigEvents=[[],[]]
		for ivent in xrange(len(rb.activecat[0])):
			if rb.activecat[4][ivent]>=5.5:
				bigEvents[0]+=[rb.activecat[0][ivent]]
				#bigEvents[1]+=[rb.activecat[4][ivent]]
				bigEvents[1]+=[4]
			#
		#
		plt.plot(bigEvents[0], bigEvents[1], 'o')
	plt.legend()
	plt.savefig("images/runningMeans.pdf")
	plt.show()
		#
	# and save a dataset of the time series with big events:
	meanfile=open(outfile, "w")
	
	meanfile.close()
	
	####
	#
	
def calcPaperSets(objRB=None):
	if objRB==None: objRB=recordbreaker()
	#
	winlen=1024
	# socal ints,mags:
	objRB.calcFullSet('cats/socalRB19842010.cat', 3.0, -90, 90, -360,360,True, winlen,1,'paperOutput/socal/intsNT')
	objRB.calcFullSetMags('cats/socalRB19842010.cat', 3.0, -90, 90, -360,360,True, winlen,1,'paperOutput/socal/magsNT')
	#cmt ints,mags:
	objRB.calcFullSet('cats/cmt7705.cat', 5.5, -90, 90, -180, 360, True, winlen, 1, 'paperOutput/cmt/intsNT')
	objRB.calcFullSetMags('cats/cmt7705.cat', 5.5, -90, 90, -180, 360, True, winlen, 1, 'paperOutput/cmt/magsNT')
	#
	# now, the parkfield aftershocks dataset:


###########################3
# workspace: get
# * wlen={32,16}
# * maybe smoothe, but over not more than wlen elements.
# * socal and CMT.

def doEverything():
	print "getting data sets."
	currDataSets()
	#NHPP3set(Nits=10, binsize=1.0
	
	#print "data sets acquired. do plots."
	#getCurrentPlots()
	print "done."
	
	
def currDataSets(objRB=None):
	if objRB==None: objRB=recordbreaker()
	#
	# set up output folders as necessary...
	rootDir='currentOutput'
	socaldir='currentOutput/socal'
	cmtdir='currentOutput/cmt'
	parkfieldshockdir='currentOutput/parkfieldshock'
	parkfieldfulldir='currentOutput/parkfieldfull'
	
	parkfieldshockdir2='currentOutput/parkfieldshockMainEvent'
	parkfieldfulldir2='currentOutput/parkfieldfullMainEvent'

	hminefulldir='currentOutput/hminefull'
	hmineshockdir='currentOutput/hmineshock'
	#
	if os.system('ls %s' % rootDir)!=0: os.system('mkdir %s' % rootDir)
	
	# make a list-n-loop outo of this whole process...
	if os.system('ls %s' % socaldir)!=0: os.system('mkdir %s' % socaldir)
	if os.system('ls %s' % cmtdir)!=0: os.system('mkdir %s' % cmtdir)
	if os.system('ls %s' % parkfieldshockdir)!=0: os.system('mkdir %s' % parkfieldshockdir)	# parkfield before and after mainshock
	if os.system('ls %s' % parkfieldfulldir)!=0: os.system('mkdir %s' % parkfieldfulldir)		# parkfield before and after mainshock
	if os.system('ls %s' % parkfieldshockdir2)!=0: os.system('mkdir %s' % parkfieldshockdir2)	# parkfield after mainshock
	if os.system('ls %s' % parkfieldfulldir2)!=0: os.system('mkdir %s' % parkfieldfulldir2)	# parkfield after mainshock
	if os.system('ls %s' % hminefulldir)!=0: os.system('mkdir %s' % hminefulldir)	# hmine square cat, before and after
	if os.system('ls %s' % hmineshockdir)!=0: os.system('mkdir %s' % hmineshockdir)	# hmine ellipse before, after
	
	#
	'''
	winlen=16
	# socal ints,mags:
	print "calc for winlen %d" % winlen
	objRB.calcFullSet('cats/socalRB19842010.cat', 3.0, -90, 90, -360,360,True, winlen,1,'currentOutput/socal/intsNT16')
	objRB.calcFullSetMags('cats/socalRB19842010.cat', 3.0, -90, 90, -360,360,True, winlen,1, 'currentOutput/socal/magsNT16')
	#cmt ints,mags:
	objRB.calcFullSet('cats/cmt7705.cat', 5.5, -90, 90, -180, 360, True, winlen, 1, 'currentOutput/cmt/intsNT16')
	objRB.calcFullSetMags('cats/cmt7705.cat', 5.5, -90, 90, -180, 360, True, winlen, 1, 'currentOutput/cmt/magsNT16')
	#
	winlen=32
	print "calc for winlen %d" % winlen
	# socal ints,mags:
	objRB.calcFullSet('cats/socalRB19842010.cat', 3.0, -90, 90, -360,360,True, winlen,1,'currentOutput/socal/intsNT32')
	objRB.calcFullSetMags('cats/socalRB19842010.cat', 3.0, -90, 90, -360,360,True, winlen,1, 'currentOutput/socal/magsNT32')
	#cmt ints,mags:
	objRB.calcFullSet('cats/cmt7705.cat', 5.5, -90, 90, -180, 360, True, winlen, 1, 'currentOutput/cmt/intsNT32')
	objRB.calcFullSetMags('cats/cmt7705.cat', 5.5, -90, 90, -180, 360, True, winlen, 1, 'currentOutput/cmt/magsNT32')
	#
	winlen=128
	print "calc for winlen %d" % winlen
	# socal ints,mags:
	objRB.calcFullSet('cats/socalRB19842010.cat', 3.0, -90, 90, -360,360,True, winlen,1,'currentOutput/socal/intsNT%d' % winlen)
	objRB.calcFullSetMags('cats/socalRB19842010.cat', 3.0, -90, 90, -360,360,True, winlen,1, 'currentOutput/socal/magsNT%d' % winlen)
	#cmt ints,mags:
	objRB.calcFullSet('cats/cmt7705.cat', 5.5, -90, 90, -180, 360, True, winlen, 1, 'currentOutput/cmt/intsNT%d' % winlen)
	objRB.calcFullSetMags('cats/cmt7705.cat', 5.5, -90, 90, -180, 360, True, winlen, 1, 'currentOutput/cmt/magsNT%d' % winlen)
	#
	winlen=256
	print "calc for winlen %d" % winlen
	# socal ints,mags:
	objRB.calcFullSet('cats/socalRB19842010.cat', 3.0, -90, 90, -360,360,True, winlen,1,'currentOutput/socal/intsNT%d' % winlen)
	objRB.calcFullSetMags('cats/socalRB19842010.cat', 3.0, -90, 90, -360,360,True, winlen,1, 'currentOutput/socal/magsNT%d' % winlen)
	#cmt ints,mags:
	objRB.calcFullSet('cats/cmt7705.cat', 5.5, -90, 90, -180, 360, True, winlen, 1, 'currentOutput/cmt/intsNT%d' % winlen)
	objRB.calcFullSetMags('cats/cmt7705.cat', 5.5, -90, 90, -180, 360, True, winlen, 1, 'currentOutput/cmt/magsNT%d' % winlen)
	
	#
	'''
	
	
	# now, the parkfield aftershocks dataset:
	# (these sets look at 2000-01-01 to 2009-09-29, roughtly event +/- 5 years)
	# (this should be along fault, +/- 5 years, more or less)
	#for i in xrange(4,8):
	for i in xrange(4,10):
		winlen=2**i
		
		# note: os.listdir() might be a better way to do this.
		if os.system('ls %s/intsNT%d' % (parkfieldshockdir, winlen))!=0: os.system('mkdir %s/intsNT%d' % (parkfieldshockdir, winlen))
		if os.system('ls %s/magsNT%d' % (parkfieldshockdir, winlen))!=0: os.system('mkdir %s/magsNT%d' % (parkfieldshockdir, winlen))
		
		print "calc for Parkfield (aftershocks), winlen %d" % winlen
		objRB.calcFullSet('cats/parkfield10yrs.cat', 1.25, -90, 90, -360,360,True, winlen,1,'%s/intsNT%d' % (parkfieldshockdir, winlen))
		objRB.calcFullSetMags('cats/parkfield10yrs.cat', 1.25, -90, 90, -360,360,True, winlen,1, '%s/magsNT%d' % (parkfieldshockdir, winlen))
	#
	
	
	
	
	# parkfield "full" catalog (not along fault, 4x4 around parkfield
	#for i in xrange(4,10):
	for i in xrange(4,10):
		winlen=2**i
		
		# note: os.listdir() might be a better way to do this.
		if os.system('ls %s/intsNT%d' % (parkfieldfulldir, winlen))!=0: os.system('mkdir %s/intsNT%d' % (parkfieldfulldir, winlen))
		if os.system('ls %s/magsNT%d' % (parkfieldfulldir, winlen))!=0: os.system('mkdir %s/magsNT%d' % (parkfieldfulldir, winlen))
		
		print "calc for Parkfield (full), winlen %d" % winlen
		objRB.calcFullSet('cats/parkfieldfull10yrs.cat', 1.25, -90, 90, -360,360,True, winlen,1,'%s/intsNT%d' % (parkfieldfulldir, winlen))
		objRB.calcFullSetMags('cats/parkfieldfull10yrs.cat', 1.25, -90, 90, -360,360,True, winlen,1, '%s/magsNT%d' % (parkfieldfulldir, winlen))
	
	#
	'''
	# parkfield only after the main event:
	for i in xrange(4,10):
		winlen=2**i
		
		# note: os.listdir() might be a better way to do this.
		if os.system('ls %s/intsNT%d' % (parkfieldshockdir2, winlen))!=0: os.system('mkdir %s/intsNT%d' % (parkfieldshockdir2, winlen))
		if os.system('ls %s/magsNT%d' % (parkfieldshockdir2, winlen))!=0: os.system('mkdir %s/magsNT%d' % (parkfieldshockdir2, winlen))
		
		print "calc for Parkfield (aftershocks), winlen %d" % winlen
		objRB.calcFullSet('cats/parkfieldShocks.cat', 1.25, -90, 90, -360,360,True, winlen,1,'%s/intsNT%d' % (parkfieldshockdir2, winlen))
		objRB.calcFullSetMags('cats/parkfieldShocks.cat', 1.25, -90, 90, -360,360,True, winlen,1, '%s/magsNT%d' % (parkfieldshockdir2, winlen))
		#
	#	
	# parkfield "full" catalog (not along fault, 4x4 around parkfield
	for i in xrange(4,10):
		winlen=2**i
		
		# note: os.listdir() might be a better way to do this.
		if os.system('ls %s/intsNT%d' % (parkfieldfulldir2, winlen))!=0: os.system('mkdir %s/intsNT%d' % (parkfieldfulldir2, winlen))
		if os.system('ls %s/magsNT%d' % (parkfieldfulldir2, winlen))!=0: os.system('mkdir %s/magsNT%d' % (parkfieldfulldir2, winlen))
		
		print "calc for Parkfield (full), winlen %d" % winlen
		objRB.calcFullSet('cats/parkfieldfullShocks.cat', 1.25, -90, 90, -360,360,True, winlen,1,'%s/intsNT%d' % (parkfieldfulldir2, winlen))
		objRB.calcFullSetMags('cats/parkfieldfullShocks.cat', 1.25, -90, 90, -360,360,True, winlen,1, '%s/magsNT%d' % (parkfieldfulldir2, winlen))
	#
	'''
	
	'''
	# hector mine "full" catalog:
	for i in xrange(4,10):
		winlen=2**i
		
		# note: os.listdir() might be a better way to do this.
		if os.system('ls %s/intsNT%d' % (hminefulldir, winlen))!=0: os.system('mkdir %s/intsNT%d' % (hminefulldir, winlen))
		if os.system('ls %s/magsNT%d' % (hminefulldir, winlen))!=0: os.system('mkdir %s/magsNT%d' % (hminefulldir, winlen))
		
		print "calc for Hector-Mine (full), winlen %d" % winlen
		objRB.calcFullSet('cats/hminefull.cat', 2.5, -90, 90, -360,360,True, winlen,1,'%s/intsNT%d' % (hminefulldir, winlen))
		objRB.calcFullSetMags('cats/hminefull.cat', 2.5, -90, 90, -360,360,True, winlen,1, '%s/magsNT%d' % (hminefulldir, winlen))
	#
	# hector mine "shock" catalog:
	for i in xrange(4,10):
		winlen=2**i
		
		# note: os.listdir() might be a better way to do this.
		if os.system('ls %s/intsNT%d' % (hmineshockdir, winlen))!=0: os.system('mkdir %s/intsNT%d' % (hmineshockdir, winlen))
		if os.system('ls %s/magsNT%d' % (hmineshockdir, winlen))!=0: os.system('mkdir %s/magsNT%d' % (hmineshockdir, winlen))
		
		print "calc for Hector-Mine (elliptical), winlen %d" % winlen
		objRB.calcFullSet('cats/hmineshock.cat', 2.5, -90, 90, -360,360,True, winlen,1,'%s/intsNT%d' % (hmineshockdir, winlen))
		objRB.calcFullSetMags('cats/hmineshock.cat', 2.5, -90, 90, -360,360,True, winlen,1, '%s/magsNT%d' % (hmineshockdir, winlen))
	#
	'''
	
	getCurrentPlots(0)

def getCurrentPlots(timecol=0, ntcol=1):
	avlen=16
	#timecol=0		# 0: time-time, 1: natural time
	#
	fignum=0
	plt.clf()
	plt.cla()
	
	'''
	print "socals..."
	objRB=recordbreaker('cats/socalRB19842010.cat', 3.0)
	bigeq=5.0
	plotMeanRatios('currentOutput/socal/intsNT16/rb-socalRB19842010-tsBig.dat', avlen, objRB, fignum, False, bigeq, timecol)	
	plt.title("Socal, winLen=16, avlen=%d" % avlen)
	plt.savefig('currentOutput/socal/intsNT16/socalBoSNT16.pdf')
	fignum+=1
	#
	avlen=16
	plotMeanRatios('currentOutput/socal/intsNT32/rb-socalRB19842010-tsBig.dat', avlen, objRB, fignum, False, bigeq, timecol)
	plt.title("Socal, winlen=32, avlen=%d" % avlen)
	plt.savefig('currentOutput/socal/intsNT32/socalBoSNT32.pdf')
	fignum+=1
	#
	avlen=16
	plotMeanRatios('currentOutput/socal/intsNT128/rb-socalRB19842010-tsBig.dat', avlen, objRB, fignum, False, bigeq, timecol)
	plt.title("Socal, winlen=128, avlen=%d" % avlen)
	plt.savefig('currentOutput/socal/intsNT128/socalBoSNT128.pdf')
	fignum+=1
	#
	avlen=16
	plotMeanRatios('currentOutput/socal/intsNT256/rb-socalRB19842010-tsBig.dat', avlen, objRB, fignum, False, bigeq, timecol)
	plt.title("Socal, winlen=256, avlen=%d" % avlen)
	plt.savefig('currentOutput/socal/intsNT256/socalBoSNT256.pdf')
	fignum+=1
	
	objRB=recordbreaker('cats/cmt7705.cat', 5.5)
	bigeq=7.5
	plotMeanRatios('currentOutput/cmt/intsNT16/rb-cmt7605-tsBig.dat', avlen, objRB, fignum, False, bigeq, timecol)
	plt.title("CMT World, winlen=16, avlen=%d" % avlen)
	plt.savefig('currentOutput/cmt/intsNT16/cmtBoSNT16.pdf')
	fignum+=1
	#
	'''
	
	'''
	# cmt globals:
	print "CMT globals..."
	plotMeanRatios('currentOutput/cmt/intsNT32/rb-cmt7605-tsBig.dat', avlen, objRB, fignum, False, bigeq, timecol)
	plt.title("CMT World, winlen=32, avlen=%d" % avlen)
	plt.savefig('currentOutput/cmt/intsNT16/cmtBoSNT32.pdf')
	fignum+=1
	#
	plotMeanRatios('currentOutput/cmt/intsNT128/rb-cmt7605-tsBig.dat', avlen, objRB, fignum, False, bigeq, timecol)
	plt.title("CMT World, winlen=128, avlen=%d" % avlen)
	plt.savefig('currentOutput/cmt/intsNT16/cmtBoSNT128.pdf')
	fignum+=1
	#
	plotMeanRatios('currentOutput/cmt/intsNT256/rb-cmt7605-tsBig.dat', avlen, objRB, fignum, False, bigeq, timecol)
	plt.title("cmt world, winlen=256, avlen=%d" % avlen)
	plt.savefig('currentOutput/cmt/intsNT16/cmtBoSNT256.pdf')
	fignum+=1
	#
	'''
	#
	# parkfield +/- 5yrs:
	print "plot parkfield shock time series..."
	#parkfieldFdate=datetimeToFloat(datetime.datetime(2004, 9, 28, 17, 15, 24))	#=731852.71902777778
	parkfieldFdate=731852.71902777778
	myfiles=os.listdir('currentOutput/parkfieldshock')
	objRB=recordbreaker('cats/parkfield10yrs.cat', 1.5)
	#bigeq=5.0
	bigeq=3.0
	imageExtensions=['pdf', 'eps', 'png']
	fsize=18
	for myfile in myfiles:
		if myfile[0:6]!='intsNT': continue
		#
		print "folder: %s" % myfile
		if os.system('ls currentOutput/parkfieldshock/%s/rb-parkfield10yrs-tsBig.dat' % myfile)!=0: continue
		#plt.clf()
		# time-time plot:
		plotMeanRatios('currentOutput/parkfieldshock/%s/rb-parkfield10yrs-tsBig.dat' % myfile, avlen, objRB, fignum, False, bigeq, timecol, 0)
		plt.annotate(" Parkfield", (parkfieldFdate, 5), (parkfieldFdate+.001, 8), arrowprops=dict(width=2), fontsize=fsize)
		#plt.title("Parkfield (catalog along fault), winlen=%s, avlen=%s" %(myfile[6:], avlen))
		#plt.title("Record Breaking Interval Ratios\n(Parkfield Aftershock Region, %s Events)" % myfile[6:])
		plt.title('')
		#plt.ylabel("n_big/n_small")
		plt.ylabel("r(t)", fontsize=fsize)
		plt.xlabel("t", fontsize=fsize)
		#plt.savefig('currentOutput/parkfieldshock/%s/parkfieldBoSNT%s.pdf' % (myfile, myfile[6:]))
		#plt.savefig('currentOutput/parkfieldshock/%s/parkfieldBoSNT%s.eps' % (myfile, myfile[6:]))
		#plt.savefig('currentOutput/parkfieldshock/%s/parkfieldBoSNT%s.png' % (myfile, myfile[6:]))
		for ext in imageExtensions:
			plt.savefig('currentOutput/parkfieldshock/%s/parkfieldBoSNT%s.%s' % (myfile, myfile[6:], ext))
			
		fignum+=1
		#
		# natural time plot:
		plotMeanRatios('currentOutput/parkfieldshock/%s/rb-parkfield10yrs-tsBig.dat' % myfile, avlen, objRB, fignum, False, bigeq, ntcol, 0)
		plt.annotate(" Parkfield", (parkfieldFdate, 5), (parkfieldFdate+.001, 8), arrowprops=dict(width=2), fontsize=fsize)
		#plt.title("Parkfield (catalog along fault) NT, winlen=%s, avlen=%s" %(myfile[6:], avlen))
		#plt.title("Record Interval Ratios\n(Parkfield Aftershock Region, %s Events)" % myfile[6:])
		plt.title('')
		#plt.ylabel("N_big/N_small")
		plt.ylabel("r(t)", fontsize=fsize)
		plt.xlabel("Number of Events", fontsize=fsize)
		#plt.savefig('currentOutput/parkfieldshock/%s/parkfieldBoSNTNT%s.pdf' % (myfile, myfile[6:]))
		#plt.savefig('currentOutput/parkfieldshock/%s/parkfieldBoSNTNT%s.eps' % (myfile, myfile[6:]))
		for ext in imageExtensions:
			plt.savefig('currentOutput/parkfieldshock/%s/parkfieldBoSNTNT%s.%s' % (myfile, myfile[6:], ext))
		
		fignum+=1
		
		# time-time plot:
		plotMeanRatios('currentOutput/parkfieldshock/%s/rb-parkfield10yrs-tsBigBack.dat' % myfile, avlen, objRB, fignum, False, bigeq, timecol, 1)
		plt.annotate(" Parkfield", (parkfieldFdate, 5), (parkfieldFdate+.001, 8), arrowprops=dict(width=2), fontsize=fsize)
		#plt.title("Parkfield (catalog along fault), winlen=%s, avlen=%s" %(myfile[6:], avlen))
		#plt.title("Backward Looking Record Interval Ratios\n(Parkfield Aftershock Region, %s Events)" % myfile[6:])
		plt.title('')
		#plt.ylabel("n_big/n_small")
		plt.ylabel("r(t)", fontsize=fsize)
		plt.xlabel("t", fontsize=fsize)
		plt.savefig('currentOutput/parkfieldshock/%s/parkfieldBoSNTBack%s.pdf' % (myfile, myfile[6:]))
		plt.savefig('currentOutput/parkfieldshock/%s/parkfieldBoSNTBack%s.eps' % (myfile, myfile[6:]))
		for ext in imageExtensions:
			plt.savefig('currentOutput/parkfieldshock/%s/parkfieldBoSNTBack%s.%s' % (myfile, myfile[6:], ext))
		fignum+=1
		#
		# natural time plot:
		plotMeanRatios('currentOutput/parkfieldshock/%s/rb-parkfield10yrs-tsBigBack.dat' % myfile, avlen, objRB, fignum, False, bigeq, ntcol, 1)
		plt.annotate(" Parkfield", (parkfieldFdate, 5), (parkfieldFdate+.001, 8), arrowprops=dict(width=2), fontsize=fsize)
		#plt.title("Parkfield (catalog along fault) NT, winlen=%s, avlen=%s" %(myfile[6:], avlen))
		#plt.title("Backward Looking Record Interval Ratios\n(Parkfield Aftershock Region, %s Events)" % myfile[6:])
		plt.title('')
		#plt.ylabel("n_big/n_small")
		plt.ylabel("r(t)", fontsize=fsize)
		plt.xlabel("Number of Events", fontsize=fsize)
		#plt.savefig('currentOutput/parkfieldshock/%s/parkfieldBoSNTNTBack%s.pdf' % (myfile, myfile[6:]))
		#plt.savefig('currentOutput/parkfieldshock/%s/parkfieldBoSNTNTBack%s.eps' % (myfile, myfile[6:]))
		for ext in imageExtensions:
			plt.savefig('currentOutput/parkfieldshock/%s/parkfieldBoSNTNTBack%s.%s' % (myfile, myfile[6:], ext))
		fignum+=1

		

	print "plot parkfield (full) time series..."
	myfiles=os.listdir('currentOutput/parkfieldfull')
	objRB=recordbreaker('cats/parkfieldfull10yrs.cat', 1.5)
	#bigeq=5.0
	bigeq=4.0
	for myfile in myfiles:
		if myfile[0:6]!='intsNT': continue
		#
		print "folder: %s" % myfile
		if os.system('ls currentOutput/parkfieldfull/%s/rb-parkfieldfull10yrs-tsBig.dat' % myfile)!=0:
			skipstr='ls currentOutput/parkfieldfull/%s/rb-parkfieldfull10yrs-tsBig.dat' % myfile
			print "skipping: %s" % skipstr
			skipstr=None
			continue
		#plt.clf()
		plotMeanRatios('currentOutput/parkfieldfull/%s/rb-parkfieldfull10yrs-tsBig.dat' % myfile, avlen, objRB, fignum, False, bigeq, timecol, 0)
		# now, draw an arrow over parkfield:
		plt.annotate(" Parkfield", (parkfieldFdate, 5), (parkfieldFdate+.001, 8), arrowprops=dict(width=2), fontsize=fsize)
		#plt.title("Record Breaking Interval Ratios\n(Parkfield Regional Seismicity, %s Events)" % myfile[6:])
		plt.title('')
		#plt.title("Parkfield (full area catalog), winlen=%s, avlen=%s" %(myfile[6:], avlen))
		#plt.ylabel("n_big/n_small")
		plt.ylabel("r(t)", fontsize=fsize)
		plt.xlabel("t", fontsize=fsize)
		#plt.savefig('currentOutput/parkfieldfull/%s/parkfieldBoSNT%s.pdf' % (myfile, myfile[6:]))
		#plt.savefig('currentOutput/parkfieldfull/%s/parkfieldBoSNT%s.eps' % (myfile, myfile[6:]))
		for ext in imageExtensions:
			plt.savefig('currentOutput/parkfieldfull/%s/parkfieldBoSNT%s.%s' % (myfile, myfile[6:], ext))
		fignum+=1
		#
		if os.system('ls currentOutput/parkfieldfull/%s/rb-parkfieldfull10yrs-tsBig.dat' % myfile)!=0:
			skipstr='ls currentOcleanAllNHPPdata()utput/parkfieldfull/%s/rb-parkfieldfull10yrs-tsBig.dat' % myfile
			print "skipping: %s" % skipstr
			skipstr=None
			continue
		#plt.clf()
		plotMeanRatios('currentOutput/parkfieldfull/%s/rb-parkfieldfull10yrs-tsBig.dat' % myfile, avlen, objRB, fignum, False, bigeq, ntcol, 0)
		# now, draw an arrow over parkfield:
		plt.annotate(" Parkfield", (parkfieldFdate, 5), (parkfieldFdate+.001, 8), arrowprops=dict(width=2), fontsize=fsize)
		#plt.title("Forward Looking Record Interval Ratios\n(Parkfield Regional Seismicity, %s Events)" % myfile[6:])
		plt.title('')
		#plt.title("Parkfield (full area catalog), winlen=%s, avlen=%s" %(myfile[6:], avlen))
		#plt.ylabel("n_big/n_small")
		plt.ylabel("r(t)", fontsize=fsize)
		plt.xlabel("Number of Events", fontsize=fsize)
		#plt.savefig('currentOutput/parkfieldfull/%s/parkfieldBoSNTNT%s.pdf' % (myfile, myfile[6:]))
		#plt.savefig('currentOutput/parkfieldfull/%s/parkfieldBoSNTNT%s.eps' % (myfile, myfile[6:]))
		for ext in imageExtensions:
			plt.savefig('currentOutput/parkfieldfull/%s/parkfieldBoSNTNT%s.%s' % (myfile, myfile[6:], ext))
		fignum+=1
		#
		#
		if os.system('ls currentOutput/parkfieldfull/%s/rb-parkfieldfull10yrs-tsBigBack.dat' % myfile)!=0:
			skipstr='ls currentOutput/parkfieldfull/%s/rb-parkfieldfull10yrs-tsBigBack.dat' % myfile
			print "skipping: %s" % skipstr
			skipstr=None
			continue
		#plt.clf()
		plotMeanRatios('currentOutput/parkfieldfull/%s/rb-parkfieldfull10yrs-tsBigBack.dat' % myfile, avlen, objRB, fignum, False, bigeq, timecol, 1)
		# now, draw an arrow over parkfield:
		plt.annotate(" Parkfield", (parkfieldFdate, 5), (parkfieldFdate+.001, 8), arrowprops=dict(width=2), fontsize=fsize)
		#plt.title("Backward Looking Record Interval Ratios\n(Parkfield Regional Seismicity, %s Events)" % myfile[6:])
		plt.title('')
		#plt.title("Parkfield (full area catalog, back-looking), winlen=%s, avlen=%s" %(myfile[6:], avlen))
		#plt.ylabel("n_big/n_small")
		plt.ylabel("r(t)", fontsize=fsize)
		plt.xlabel("t", fontsize=fsize)
		#plt.savefig('currentOutput/parkfieldfull/%s/parkfieldBoSNTBack%s.pdf' % (myfile, myfile[6:]))
		#plt.savefig('currentOutput/parkfieldfull/%s/parkfieldBoSNTBack%s.eps' % (myfile, myfile[6:]))
		for ext in imageExtensions:
			plt.savefig('currentOutput/parkfieldfull/%s/parkfieldBoSNTBack%s.%s' % (myfile, myfile[6:], ext))
		fignum+=1
		
		if os.system('ls currentOutput/parkfieldfull/%s/rb-parkfieldfull10yrs-tsBigBack.dat' % myfile)!=0:
			skipstr='ls currentOutput/parkfieldfull/%s/rb-parkfieldfull10yrs-tsBigBack.dat' % myfile
			print "skipping: %s" % skipstr
			skipstr=None
			continue
		#plt.clf()
		plotMeanRatios('currentOutput/parkfieldfull/%s/rb-parkfieldfull10yrs-tsBigBack.dat' % myfile, avlen, objRB, fignum, False, bigeq, ntcol, 1)
		# now, draw an arrow over parkfield:
		plt.annotate(" Parkfield", (parkfieldFdate, 5), (parkfieldFdate+.001, 8), arrowprops=dict(width=2))
		#plt.title("Backward Looking Record Interval Ratios\n(Parkfield Regional Seismicity, %s Events)" % myfile[6:])
		plt.title('')
		#plt.ylabel("n_big/n_small")
		plt.ylabel("r(t)")
		plt.xlabel("Number of Events")
		#plt.title("Parkfield (full area catalog, back-looking), winlen=%s, avlen=%s" %(myfile[6:], avlen))
		#plt.savefig('currentOutput/parkfieldfull/%s/parkfieldBoSNTNTBack%s.pdf' % (myfile, myfile[6:]))
		#plt.savefig('currentOutput/parkfieldfull/%s/parkfieldBoSNTNTBack%s.eps' % (myfile, myfile[6:]))
		for ext in imageExtensions:
			plt.savefig('currentOutput/parkfieldfull/%s/parkfieldBoSNTNTBack%s.%s' % (myfile, myfile[6:], ext))		
		fignum+=1
		
	
	#
	#############################
	# parkfield post-main event:
	print "plot parkfield shock time series, post main event..."
	myfiles=os.listdir('currentOutput/parkfieldshockMainEvent')
	objRB=recordbreaker('cats/parkfieldShocks.cat', 1.25)
	#bigeq=5.0
	bigeq=4.0
	for myfile in myfiles:
		if myfile[0:6]!='intsNT': continue
		#
		if os.system('ls currentOutput/parkfieldshockMainEvent/%s/rb-parkfieldShocks-tsBig.dat' % myfile)!=0: continue
		#print "folder: %s" % myfile
		#plt.clf()
		plotMeanRatios('currentOutput/parkfieldshockMainEvent/%s/rb-parkfieldShocks-tsBig.dat' % myfile, avlen, objRB, fignum, False, bigeq, timecol, 0)
		plt.title("Post-Parkfield (catalog along fault), winlen=%s, avlen=%s" %(myfile[6:], avlen))
		plt.savefig('currentOutput/parkfieldshockMainEvent/%s/parkfieldBoSNT%s.pdf' % (myfile, myfile[6:]))
		fignum+=1
	#
	print "plot parkfield (full) time series, post main event..."
	myfiles=os.listdir('currentOutput/parkfieldfullMainEvent')
	objRB=recordbreaker('cats/parkfieldfullShocks.cat', 1.25)
	#bigeq=5.0
	bigeq=4.0
	for myfile in myfiles:
		if myfile[0:6]!='intsNT': continue
		#
		if os.system('ls currentOutput/parkfieldfullMainEvent/%s/rb-parkfieldfullShocks-tsBig.dat' % myfile)!=0: continue
		#print "folder: %s" % myfile
		#plt.clf()
		plotMeanRatios('currentOutput/parkfieldfullMainEvent/%s/rb-parkfieldfullShocks-tsBig.dat' % myfile, avlen, objRB, fignum, False, bigeq, timecol)
		plt.title("Post-Parkfield (full area catalog), winlen=%s, avlen=%s" %(myfile[6:], avlen))
		plt.savefig('currentOutput/parkfieldfullMainEvent/%s/parkfieldBoSNT%s.pdf' % (myfile, myfile[6:]))
		fignum+=1
	
	#
	# hector mine:
		#
	print "plot Hector Mine (full) time series..."
	myfiles=os.listdir('currentOutput/hminefull')
	objRB=recordbreaker('cats/hminefull.cat', 2.5)
	#bigeq=5.0
	bigeq=4.0
	for myfile in myfiles:
		if myfile[0:6]!='intsNT': continue
		#
		if os.system('ls currentOutput/hminefull/%s/rb-hminefull-tsBig.dat' % myfile)!=0: continue
		print "folder: %s" % myfile
		#plt.clf()
		plotMeanRatios('currentOutput/hminefull/%s/rb-hminefull-tsBig.dat' % myfile, avlen, objRB, fignum, False, bigeq, timecol)
		plt.title("Hector Mine (square catalog), winlen=%s, avlen=%s" %(myfile[6:], avlen))
		#plt.savefig('currentOutput/hminefull/%s/hminefullBoSNT%s.pdf' % (myfile, myfile[6:]))
		for ext in imageExtensions:
			plt.savefig('currentOutput/hminefull/%s/hminefullBoSNT%s.%s' % (myfile, myfile[6:], ext))	
		fignum+=1
	print "plot Hector Mine elliptical time series..."
	myfiles=os.listdir('currentOutput/hmineshock')
	objRB=recordbreaker('cats/hmineshock.cat', 2.5)
	#bigeq=5.0
	bigeq=4.0
	for myfile in myfiles:
		if myfile[0:6]!='intsNT': continue
		#
		if os.system('ls currentOutput/hmineshock/%s/rb-hmineshock-tsBig.dat' % myfile)!=0: continue
		print "folder: %s" % myfile
		#plt.clf()
		plotMeanRatios('currentOutput/hmineshock/%s/rb-hmineshock-tsBig.dat' % myfile, avlen, objRB, fignum, False, bigeq, timecol)
		plt.title("Hector Mine (elliptical catalog), winlen=%s, avlen=%s" %(myfile[6:], avlen))
		#plt.savefig('currentOutput/hmineshock/%s/hmineshockBoSNT%s.pdf' % (myfile, myfile[6:]))
		for ext in imageExtensions:
			plt.savefig('currentOutput/hmineshock/%s/hmineshockBoSNT%s.%s' % (myfile, myfile[6:], ext))	
		fignum+=1
	#
	# replot the means:
	replotRBmeans()
	#	
	#plt.show()
	#plt.clf()
#
def plotRBMeanMagsFile(fname, catalog='socal', saveDir='paperOut/socal'):
	# plot mean values (aka, mean NRB(nEvents)) from an output file such as rb-socalRB19842010-ints-gnu.dat
	#if catalog not in ('socal', 'cmt', 'parkfieldsquare', 'parkfieldsquareshock', 'parkfieldfault', 'parkfieldfaultshock'): catalog='socal'
	strNT='NT'	# time-time? natural time?
	#
	if catalog=='socal':
		titleString="Record Breaking Earthquake Magnitudes\n(Southern California, 1984-2010, m>=3.0)"	# obviously, mag, etc. could change...
		catStr='SoCal'
		xlbl="Number of Events (Natural Time)"
		ylblNrb="Number of Record Breaking Events (Magnitudes)"
		ylblMag="Magnitude of Record Breaking Event"
		
	elif catalog=='cmt':
		titleString="Record Breaking Earthquake Magnitudes\n(Global CMT, 1976-2005, m>=5.5)"	# obviously, mag, etc. could change...
		catStr='CMT'
		xlbl="Number of Events (Natural Time)"
		ylblNrb="Number of Record Breaking Events (Magnitudes)"
		ylblMag="Magnitude of Record Breaking Event"
	
	#elif catalog=='parkfieldsquare':
	elif catalog in ['parkfieldsquare', 'parkfieldfull']:
		# square catalog 4x4 degrees centered on parkfield, +/- 5 years (roughly)
		titleString="Record Breaking Earthquake Magnitudes \n(4x4 deg centered on parkfield, 2000-01-01 - 2009-29-09)"
		catStr='Parkfield'
		xlbl="Number of Events (Natural Time)"
		ylblNrb="Number of Record Breaking Events (Intervals)"
		ylblMag="Magnitude of Record Breaking Event"
		
	#if catalog=='parkfieldsquareshock':
	elif catalog in ['parkfieldsquareshock', 'parkfieldfullMainEvent']:
		# square catalog 4x4 degrees centered on parkfield, after main event
		titleString="Record Breaking Earthquake Magnitudes \n(4x4 deg centered on parkfield, after 2004-09-28 mainshock)"
		catStr='Parkfield'
		xlbl="Number of Events (Natural Time)"
		ylblNrb="Number of Record Breaking Events (Intervals)"
		ylblMag="Magnitude of Record Breaking Event"
		
	#if catalog=='parkfieldfault':
	elif catalog in ['parkfieldfault', 'parkfieldshock']:
		# along the fault, +/- 5 years (roughly)
		titleString="Record Breaking Earthquake Magnitudes \n(Ellipse centered on parkfield,  2000-01-01 - 2009-29-09)"
		catStr='Parkfield'
		xlbl="Number of Events (Natural Time)"
		ylblNrb="Number of Record Breaking Events (Intervals)"
		ylblMag="Magnitude of Record Breaking Event"
		
	#if catalog=='parkfieldfaultshock':
	elif catalog in ['parkfieldfaultshock', ' parkfieldshockMainEvent']:
		# along the fault, after main event (28 september 2004
		titleString="Record Breaking Earthquake Magnitudes \n(Ellipse centered on parkfield, after 2004-09-28 mainshock)"
		catStr='Parkfield'
		xlbl="Number of Events (Natural Time)"
		ylblNrb="Number of Record Breaking Events (Intervals)"
		ylblMag="Magnitude of Record Breaking Event"
	elif catalog in ['hminefull']:
		# along the fault, after main event (28 september 2004
		titleString="Record Breaking Earthquake Magnitudes \n(square catalog, Hector Mine)"
		catStr='H.Mine'
		xlbl="Number of Events (Natural Time)"
		ylblNrb="Number of Record Breaking Events (Intervals)"
		ylblMag="Magnitude of Record Breaking Event"
	elif catalog in ['hmineshock']:
		# along the fault, after main event (28 september 2004
		titleString="Record Breaking Earthquake Magnitudes \n(Elliptical catalog, Hector Mine)"
		catStr='H.Mine'
		xlbl="Number of Events (Natural Time)"
		ylblNrb="Number of Record Breaking Events (Intervals)"
		ylblMag="Magnitude of Record Breaking Event"
	elif catalog in ['hminefull']:
		# along the fault, after main event (28 september 2004
		titleString="Record Breaking Earthquake Magnitudes \n(square catalog, Hector Mine)"
		catStr='H.Mine'
		xlbl="Number of Events (Natural Time)"
		ylblNrb="Number of Record Breaking Events (Intervals)"
		ylblMag="Magnitude of Record Breaking Event"
	else:
		titleString="Record Breaking Earthquake Magnitudes\n(secret catalog)"	# obviously, mag, etc. could change...
		catStr='SoCal'
		xlbl="Number of Events (Natural Time)"
		ylblNrb="Number of Record Breaking Events (Magnitudes)"
		ylblMag="Magnitude of Record Breaking Event"
	#
	#
	# get arrays from the data file:
	dataWidth=9	# eventually this might be more dynamic
	f=open(fname)
	dset=[[]]
	for i in xrange(dataWidth):
		dset+=[[]]
	for rw in f:
		if rw=='\n': continue
		if rw[0]=='#' or rw[0]=='\t' or rw[0]=='\n': continue	#comments...
		rws=rw.split('\t')
		#
		dset[0]+=[long(rws[0])]
		dset[1]+=[long(2.0**long(rws[0]))]
		for j in range(1, dataWidth, 1):
			dset[j+1]+=[float(rws[j])]
		#
	#
	#	
	plt.figure(0)
	plt.title(titleString)
	plt.xlabel(xlbl)
	plt.ylabel(ylblNrb)
	plt.errorbar(dset[0],dset[2], yerr=dset[3], fmt='.-', label="%s Large" % catStr)
	plt.errorbar(dset[0],dset[4], yerr=dset[5], fmt='.-', label="%s Small" % catStr)
	plt.xticks(dset[0], dset[1])
	plt.legend(loc='upper left')
	#plt.show()
	plt.savefig("%s/rb-%s-ints%s-nrbBigSmall.pdf" % (saveDir, catalog, strNT))
	plt.clf()
	
	plt.figure(1)
	plt.title(titleString)
	plt.xlabel(xlbl)
	plt.ylabel(ylblMag)
	plt.errorbar(dset[0],dset[6], yerr=dset[7], fmt='.-', label="%s Large" % catStr)
	plt.errorbar(dset[0],dset[8], yerr=dset[9], fmt='.-', label="%s Small" % catStr)
	plt.xticks(dset[0], dset[1])
	plt.legend(loc='upper left')
	#plt.show()
	plt.savefig("%s/rb-%s-ints%s-magsBigSmall.pdf" % (saveDir, catalog, strNT))
	plt.clf()

def replotRBmeans():
	# main folders (this could be automated with better organization)
	rootdir='currentOutput'
	folders1=['cmt', 'parkfieldfull', 'parkfieldshock', 'parkfieldfullMainEvent', 'parkfieldshockMainEvent', 'socal', 'hminefull', 'hmineshock']
	for thisdir in folders1:
		fls=os.listdir('%s/%s' % (rootdir,thisdir))
		for subdir in fls:
			currdir='%s/%s/%s' % (rootdir, thisdir, subdir)
			#print currdir
			# look for the file that's like 'rb-*-gnu.dat'
			thesefiles=os.listdir(currdir)
			for thisfile in thesefiles:
				if thisfile[0:3]=='rb-' and (thisfile[-13:]=='-ints-gnu.dat' or thisfile[-13:]=='-mags-gnu.dat'):
					workingfile="%s/%s" % (currdir, thisfile)
					break
			#
			strCat=thisdir
			#if 'parkfield' in strCat: strCat='parkfield'
			#plotRBMeanFile(workingfile, strCat, currdir)
			if subdir[0:4]=='ints':
				print "plotIntsPrams: %s, %s, %s" % (workingfile, strCat, currdir)
				plotRBMeanFile(workingfile, strCat, currdir)
			elif subdir[0:4]=='mags':
				print "plotMagsPrams: %s, %s, %s" % (workingfile, strCat, currdir)
				#plotRBMeanMagsFile(fname, catalog='socal', saveDir='paperOut/socal'):
				plotRBMeanMagsFile(workingfile, strCat, currdir)
			else:
				print "unknonwPrams:  %s, %s, %s" % (workingfile, strCat, currdir)

def plotRBMeanFile(fname, catalog='socal', saveDir='paperOutput/socal'):
	# plot mean values (aka, mean NRB(nEvents)) from an output file such as rb-socalRB19842010-ints-gnu.dat
	#if catalog not in ('socal', 'cmt', 'parkfieldsquare', 'parkfieldsquareshock', 'parkfieldfault', 'parkfieldfaultshock'): catalog='socal'
	strNT='NT'	# time-time? natural time?
	#
	xlbl="Number of Events (Natural Time)"
	ylblNrb="Number of Record Breaking Events"
	ylblDuration="Duration of Record Breaking Interval (days)"
	if catalog=='socal':
		titleString="Record Breaking Earthquake Intervals\n(Southern California, 1984-2010, m>=3.0)"	# obviously, mag, etc. could change...
		catStr='SoCal'
		xlbl="Number of Events (Natural Time)"
		ylblNrb="Number of Record Breaking Events (Intervals)"
		ylblMag="Duration of Record Breaking Interval (days)"
	
	elif catalog=='cmt':
		titleString="Record Breaking Earthquake Intervals\n(Global CMT, 1977-2005, m>=5.5)"	# obviously, mag, etc. could change...
		catStr='CMT'
		xlbl="Number of Events (Natural Time)"
		ylblNrb="Number of Record Breaking Events"
		ylblMag="Duration of Record Breaking Interval (days)"
	
	# the parkfield nomenclature is a little confusing. we look at parkfield: along-fault, square around event, +/- 5 years, just after main event.
	# sometimes we call "shock" events after the main 2004 event, sometimes we refer to the spatial distribution (ellipse along fault)
	elif catalog in ['parkfieldsquare', 'parkfieldfull']:
		# square catalog 4x4 degrees centered on parkfield, +/- 5 years (roughly)
		titleString="Record Breaking Earthquake Intervals \n(4x4 deg centered on parkfield, 2000-01-01 - 2009-29-09)"
		catStr='Parkfield'
		xlbl="Number of Events (Natural Time)"
		ylblNrb="Number of Record Breaking Events (Intervals)"
		ylblMag="Duration of Record Breaking Interval (days)"
		
	elif catalog in ['parkfieldsquareshock', 'parkfieldfullMainEvent']:
		# square catalog 4x4 degrees centered on parkfield, after main event
		titleString="Record Breaking Earthquake Intervals \n(4x4 deg centered on parkfield, after 2004-09-28 mainshock)"
		catStr='Parkfield'
		xlbl="Number of Events (Natural Time)"
		ylblNrb="Number of Record Breaking Events (Intervals)"
		ylblMag="Duration of Record Breaking Interval (days)"
		
	elif catalog in ['parkfieldfault', 'parkfieldshock']:
		# along the fault, +/- 5 years (roughly)
		titleString="Record Breaking Earthquake Intervals \n(Ellipse centered on parkfield,  2000-01-01 - 2009-29-09)"
		catStr='Parkfield'
		xlbl="Number of Events (Natural Time)"
		ylblNrb="Number of Record Breaking Events (Intervals)"
		ylblMag="Duration of Record Breaking Interval (days)"
		
	elif catalog in ['parkfieldfaultshock', 'parkfieldshockMainEvent']:
		# along the fault, after main event (28 september 2004
		titleString="Record Breaking Earthquake Intervals \n(Ellipse centered on parkfield, after 2004-09-28 mainshock)"
		catStr='Parkfield'
		xlbl="Number of Events (Natural Time)"
		ylblNrb="Number of Record Breaking Events (Intervals)"
		ylblMag="Duration of Record Breaking Interval (days)"
	elif catalog in ['hminefull']:
		# along the fault, after main event (28 september 2004
		titleString="Record Breaking Earthquake Intervals \n(Hector Mine, square catalog)"
		catStr='Hmine'
		xlbl="Number of Events (Natural Time)"
		ylblNrb="Number of Record Breaking Events (Intervals)"
		ylblMag="Duration of Record Breaking Interval (days)"
	elif catalog in ['hmineshock']:
		# along the fault, after main event (28 september 2004
		titleString="Record Breaking Earthquake Intervals \n(Hector Mine, Elliptical catalog)"
		catStr='Hmine'
		xlbl="Number of Events (Natural Time)"
		ylblNrb="Number of Record Breaking Events (Intervals)"
		ylblMag="Duration of Record Breaking Interval (days)"	
	else:
		titleString="Record Breaking Earthquake Intervals\n(secret catalog)"	# obviously, mag, etc. could change...
		catStr='SoCal'
		xlbl="Number of Events (Natural Time)"
		ylblNrb="Number of Record Breaking Events (Intervals)"
		ylblMag="Duration of Record Breaking Interval (days)"
		
	#
	# get arrays from the data file:
	dataWidth=17	# eventually this might be more dynamic
	f=open(fname)
	dset=[[]]
	for i in xrange(dataWidth):
		dset+=[[]]
	for rw in f:
		if rw[0]=='#' or rw[0]=='\t' or rw[0]=='\n': continue	#comments...
		rws=rw.split('\t')
		#
		dset[0]+=[long(rws[0])]
		dset[1]+=[long(2.0**long(rws[0]))]
		for j in range(1, dataWidth, 1):
			dset[j+1]+=[float(rws[j])]
		#
	#
	#	
	plt.figure(0)
	plt.clf()
	plt.title(titleString)
	plt.xlabel(xlbl)
	plt.ylabel(ylblNrb)
	plt.errorbar(dset[0],dset[2], yerr=dset[3], fmt='.-', label="%s Large" % catStr)
	plt.errorbar(dset[0],dset[4], yerr=dset[5], fmt='.-', label="%s Small" % catStr)
	plt.errorbar(dset[0],dset[10], yerr=dset[11], fmt='.-', label="Poisson Large")
	plt.errorbar(dset[0],dset[12], yerr=dset[13], fmt='.-', label="Poisson Small")
	plt.xticks(dset[0], dset[1])
	plt.legend(loc='upper left')
	#plt.show()
	plt.savefig("%s/rb-%s-ints%s-nrbBigSmall.pdf" % (saveDir, catalog, strNT))
	plt.clf()
	
	plt.figure(1)
	plt.clf()
	plt.title(titleString)
	plt.xlabel(xlbl)
	plt.ylabel(ylblMag)
	plt.errorbar(dset[0],dset[6], yerr=dset[7], fmt='.-', label="%s Large" % catStr)
	plt.errorbar(dset[0],dset[8], yerr=dset[9], fmt='.-', label="%s Small" % catStr)
	plt.errorbar(dset[0],dset[14], yerr=dset[15], fmt='.-', label="Poisson Large")
	plt.errorbar(dset[0],dset[16], yerr=dset[17], fmt='.-', label="Poisson Small")
	plt.xticks(dset[0], dset[1])
	plt.legend(loc='upper left')
	#plt.show()
	plt.savefig("%s/rb-%s-ints%s-magsBigSmall.pdf" % (saveDir, catalog, strNT))
	plt.clf()

def cleanAllNHPPdata():
	'''
	cleanUpNHPPdata('images/NHPPcat3/mc10/binnedNRBnt.dat')
	cleanUpNHPPdata('images/NHPPcat3/mc15/binnedNRBnt.dat')
	cleanUpNHPPdata('images/NHPPcat3/mc20/binnedNRBnt.dat')
	cleanUpNHPPdata('images/NHPPcat3/mc25/binnedNRBnt.dat')
	cleanUpNHPPdata('images/NHPPcat3/mc30/binnedNRBnt.dat')
	
	cleanUpNHPPdata('images/NHPPcat3/mc10/binnedintNT.dat')
	cleanUpNHPPdata('images/NHPPcat3/mc15/binnedintNT.dat')
	cleanUpNHPPdata('images/NHPPcat3/mc20/binnedintNT.dat')
	cleanUpNHPPdata('images/NHPPcat3/mc25/binnedintNT.dat')
	cleanUpNHPPdata('images/NHPPcat3/mc30/binnedintNT.dat')

	cleanUpNHPPdata('images/NHPPcat3/mc10/binnedNRB.dat')
	cleanUpNHPPdata('images/NHPPcat3/mc15/binnedNRB.dat')
	cleanUpNHPPdata('images/NHPPcat3/mc20/binnedNRB.dat')
	cleanUpNHPPdata('images/NHPPcat3/mc25/binnedNRB.dat')
	cleanUpNHPPdata('images/NHPPcat3/mc30/binnedNRB.dat')
	
	cleanUpNHPPdata('images/NHPPcat3/mc10/binnedint.dat')
	cleanUpNHPPdata('images/NHPPcat3/mc15/binnedint.dat')
	cleanUpNHPPdata('images/NHPPcat3/mc20/binnedint.dat')
	cleanUpNHPPdata('images/NHPPcat3/mc25/binnedint.dat')
	cleanUpNHPPdata('images/NHPPcat3/mc30/binnedint.dat')
	'''
	
	flist=['images/NHPPcat3/mc10/binnedNRBnt.dat', 'images/NHPPcat3/mc15/binnedNRBnt.dat', 'images/NHPPcat3/mc20/binnedNRBnt.dat', 'images/NHPPcat3/mc25/binnedNRBnt.dat', 'images/NHPPcat3/mc30/binnedNRBnt.dat', 'images/NHPPcat3/mc10/binnedintNT.dat', 'images/NHPPcat3/mc15/binnedintNT.dat', 'images/NHPPcat3/mc20/binnedintNT.dat', 'images/NHPPcat3/mc25/binnedintNT.dat', 'images/NHPPcat3/mc30/binnedintNT.dat', 'images/NHPPcat3/mc10/binnedNRB.dat', 'images/NHPPcat3/mc15/binnedNRB.dat', 'images/NHPPcat3/mc20/binnedNRB.dat', 'images/NHPPcat3/mc25/binnedNRB.dat', 'images/NHPPcat3/mc30/binnedNRB.dat', 'images/NHPPcat3/mc10/binnedint.dat', 'images/NHPPcat3/mc15/binnedint.dat', 'images/NHPPcat3/mc20/binnedint.dat', 'images/NHPPcat3/mc25/binnedint.dat', 'images/NHPPcat3/mc30/binnedint.dat']
	for fl in flist:
		if os.system("ls %s" % fl)!=0: continue
		cleanUpNHPPdata(fl)

	
def cleanUpNHPPdata(dfile, dfilelog2=None):
	# when we make NHPP data, we end up artifically filling the end of the series. this looks like crap when we plot it.
	# start at the end of the file, throw out the last row until tha value changes (end of the filler).
	# also, export the final (real) value and log2 values to the log2-data file. this is used to produce uncrowded errorbars
	# in log-scaled plots.
	# nominally, we should save the original file. maybe we just comment away the ends...
	#
	if dfilelog2==None:
		fnames=dfile.split('.', 1)
		dfilelog2="%s-log2.%s" % (fnames[0], fnames[1])
	
	f1=open(dfile, 'r')	# load into an array, then rewrite.
	aryData=[]
	aryDataLog=[]
	newData=[]
	comments=[]
	#
	# first, read data into arrays:
	for rw in f1:
		# data are [time-ish, nrb or rbInt, stdev]
		#if rw[0]=="#" or rw[0]==" " or rw[0]=="\t":
		if rw[0] in ["#", " ", "\t", "\n"]:
			# a comment:
			comments+=[rw]
			continue
		rws=rw.split('\t')
		# below is code to convert to floats, but we don't actually care. string comparison shoudl suffice.
		#for i in xrange (len(rew)):
		#	rws[i]=float(rws[i]))
		#
		aryData+=[rws]
		if log2(float(rws[0]))%1==0: aryDataLog+=[rws]
	f1.close()
	#
	# so the log-binned array should be ready to go, unless we want to append the last good value from the raw data.
	# now, start at the end of the raw-data and walk up until we find a new value.
	# output commented rows until we find a new value.
	#fnames=dfile.split('.', 1)
	#dfileClean="%s-plot.%s" % (fnames[0], fnames[1])
	fnames=None
	#
	aryLen=len(aryData)
	thisVals=[aryData[aryLen-1][1], aryData[aryLen-1][2]]	# val and stdev.
	isFirstNewVal=0
	for i in xrange(1, aryLen):
		prevVals=thisVals
		thisVals=[aryData[aryLen-1-i][1], aryData[aryLen-1-i][2]]
		#
		if (thisVals==prevVals and isFirstNewVal==0):
			newData+=[["#%s" % aryData[aryLen-1-i][0], aryData[aryLen-1-i][1], aryData[aryLen-1-i][2]]]
		else:
			# it's one of many new-vals (allowing for repeats during the realy sequence):
			if isFirstNewVal==0:
				# we've found the first new value:
				isFirstNewVal=1
				# record the last of the filler-vals in log2 (previous row):
				#aryDataLog+=[[aryData[aryLen-i][0], aryData[aryLen-i][1], aryData[aryLen-i][0]]]	# this error bar usually looks like crap...
			newData+=[[aryData[aryLen-1-i][0], aryData[aryLen-1-i][1], aryData[aryLen-1-i][2]]]
	#
	#f1.open(dfileClean, 'w')
	fnames=dfile.split('.', 1)
	dfileRaw="%s-raw.%s" % (fnames[0], fnames[1])
	if os.system('ls %s' % dfileRaw)!=0: os.system('mv %s %s' % (dfile, dfileRaw))
	f1=open(dfile, 'w')
	f2=open(dfilelog2, 'w')
	#
	f1.write("#cleaned up NHPP data\n")
	f2.write("#log2-binned NHPP data\n")
	for rw in comments:
		if rw in ["\t", "\n"]: continue
		f1.write(rw)
		f2.write(rw)
	
	for rw in aryDataLog:
		strRow=''
		for elem in rw:
			strRow+="%s\t" % elem
		f2.write("%s" % strRow[:-1])
	f2.close()
		
	#
	# the main data array is inverted. start from the bottom:
	for i in xrange(len(newData)):
		rw=newData[len(newData)-1-i]
		strRow=''
		for elem in rw:
			strRow+="%s\t" % elem
		f1.write("%s" % strRow[:-1])
	f1.close()

def plotMagDists():
	fnames=["cats/cmt7705.cat", "cats/parkfield10yrs.cat", "cats/parkfieldfull10yrs.cat", "cats/socalRB19822009.cat", "cats/hminefull.cat", "cats/hmineshock.cat"]
	for fname in fnames:
		plotMagDist(fname, False)
	print "mag dists finished."

def plotMagDist(catfile="cats/cmt7705.cat", doshow=True):
	# assume format like: date, time, lat, lon, mag, interval(?)	# note: columns are tab separated; commas are just for show here.
	colDelim='\t'
	#
	# no binning. just get the magnitudes, sort wee-tiny to big, then integrate N...
	#fname="writeup/cmtMagDistribution"
	# get fname from catfile:
	fnameRoot=catfile.split("/")[-1].split(".")[0] + "MagDist"
	outdir="writeup/"
	#
	mags=[]
	f=open(catfile, 'r')
	for row in f:
		if row[0]=="#": continue	# comment
		mags+=[float(row.split('\t')[4])]
	#	#
	#
	mags.sort()	# this should sort wee-tiny to big, but just in case the standard changes one day:
	if mags[0]<mags[-1]: mags.reverse()
	numEvents=range(1,len(mags)+1)
	print "lens: %d, %d" % (len(mags), len(numEvents))
	#
	plt.figure(0)
	plt.clf()
	plt.cla()
	plt.semilogy(mags, numEvents, '.-')
	plt.title("%s Catalog Magnitude Distribution\n(1977-2009)" % fnameRoot)
	plt.xlabel("Magnitude")
	plt.ylabel("Number of Events")
	plt.savefig("%s%s.jpg" % (outdir, fnameRoot))
	plt.savefig("%s%s.eps" % (outdir, fnameRoot))

	if doshow: plt.show()


def getNHPPmeans(Nits=1000, Tmax=90.0, tao=9.65*10**(-7), cm=5.33*10**(-3)):
	binWidth=.1
	NHPPs=[[], []]	# t, <y>
	t=0
	print "initialize..."
	while t<Tmax:
		NHPPs[0]+=[t]
		NHPPs[1]+=[0]	# counts
		t+=binWidth
	print "initialized... %d" % len(NHPPs[0])
	#
	rb1=recordbreaker()
	for i in xrange(Nits):
		intervals=rb1.getNHPPintervalsOmori1b(Tmax, tao, cm)
		print "intervals len: %d" % len(intervals[0])
		for ii in xrange(len(intervals[0])):
			tbin=int(intervals[1][ii]/(binWidth))
			NHPPs[1][tbin]+=intervals[2][ii]/float(Nits)
	
	# now, plot distribution...
	Nvals=range(1, len(NHPPs[1])+1)
	Nvals.reverse()
	Xvals=NHPPs[1]
	Xvals.sort()
	plt.figure(1)
	plt.clf()
	plt.loglog(Xvals, Nvals)
	plt.show()
	
	return NHPPs
	
	
