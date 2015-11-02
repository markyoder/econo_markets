import math
from math import *	
#from scipy import *
import scipy
from pylab import *
from matplotlib import *
#import numpy.fft as nft
import scipy.optimize as spo
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse

# maping bits:
import matplotlib	# note that we've tome from ... import *. we should probably eventually get rid of that and use the matplotlib namespace.
matplotlib.use('Agg')
#from matplotlib.toolkits.basemap import Basemap
from mpl_toolkits.basemap import Basemap as Basemap
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
####
#
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

import operator

import yodapy as yp

#import omoribits as ob

###############################################
# rbIntervals.py                              #
#
# mark yoder
# UC Davis
# april 2009
#
# compute record breaking intervals between aftershocks
# namely, we will use events following the Parkfield Eq.
# catalog, etc. will be described inline.
###############################################

class intervalRecordBreaker:
	fullCat=[]
	shockCat=[]		# [[evDateTime, lat, lon, mag, a, b], [row1], ...]. where a,b are the eliptical x,y coordinates
	#catname='parkcat.cat'
	catname='cats/parkfieldfull10yrs.cat'
	eventDtTime=None
	
	# coordinate transformation variables:
	tLat=None # 35.9
	tLon=None # -120.5
	tTheta=None # 40.0		#47?	note: tTheta is the angle CCW of the x' (transformed) axis from the x axis.
	tA=None # .4		# ellipse axes
	tB=None # .15
	
	
	def __init__(self, catFname='cats/parkfieldfull10yrs.cat', theta=40.0, clat=35.9, clon=-120.5, ra=.4, rb=.15, eventDate=datetime.datetime(2004,9,28, 17,15,24), maxDate=datetime.datetime(2009,9,28, 17,15,24), skipSeconds=0):
		# don't know yet. we'll need a catalog and bits like that...
		self.initialize(catFname, theta, clat, clon, ra, rb)
	
	def initialize(self, catFname='cats/parkfieldfull10yrs.cat', theta=40.0, clat=35.9, clon=-120.5, ra=.4, rb=.15, eventDate=datetime.datetime(2004,9,28, 17,15,24), maxDate=datetime.datetime(2009,9,28, 17,15,24), skipSeconds=0):
		self.catname=catFname
		self.tLat=clat
		self.tLon=clon
		self.tTheta=theta
		self.tA=ra
		self.tB=rb
		
		if catFname==None: return None
		
		if theta!=None:
			print "set aftershock catalog"
			self.setAftershockCatalog(catFname, theta, clat, clon, ra, rb, eventDate, maxDate, skipSeconds )
		else:
			self.setNormalCat(catFname)	# default to socal...
	
	# 29 apr 2009 yoder: we want to do RB for a 2 year or so period in socal. create a normal active catalog.
	# (but then we abandoned this for a renewed approach; see recordBreaker.py)
	def setNormalCat(self, catFname=None, minlat=32, maxlat=36.5, minlon=-125, maxlon=-115, minDt=datetime.datetime(1984,01,01), maxDt=datetime.datetime(1985,12,31)):
		# set a normal (rectangular) catalog, no ellipse.
		if catFname!=None:
			self.setFullCat(catFname)
			
		tempcat=self.fullCat
		if minDt==None: minDt=tempcat[0][0]
		if maxDt==None: maxDt=tempcat[0][-1]
		
		for row in tempcat:
			# no rotations; just transfer the events in the space,time,mag space.
			if row[0]<=minDt: continue
			if row[0]>maxDt: break
		self.shockCat+=[[row[0], row[1], row[2], row[3], row[2], row[1]]]
		#
	#
	##################
	
	#def setAftershockCatalog(self, catFname=None, theta=tTheta, clat=35.9, clon=-120.5, ra=tA, rb=tB, eventDate=datetime.datetime(2004,9,28, 17,15,24), maxDate=datetime.datetime(2009,9,28, 17,15,24), skipNevents=0):
	def setAftershockCatalog(self, catFname=None, theta=tTheta, clat=35.9, clon=-120.5, ra=tA, rb=tB, eventDate=datetime.datetime(2004,9,28, 17,15,24), maxDate=datetime.datetime(2009,9,28, 17,15,24), skipSeconds=0):
		# keep all events in catFname defined by the following ellipse:
		# skipNevents: number of events after main-shock to skip (avoid messiness during mainshock)
		self.eventDtTime=eventDate
		self.catname=catFname
		self.shockCat=[]
		#print "event (start) date, catname: %s, %s, %s" % (eventDate, catFname, self.catname)
		#
		if catFname!=None:
			print "setting catalog from [%s]" % catFname
			self.setFullCat(catFname)
		#print "first cat row: %s" % self.fullCat[0]
		#
		tempCat=self.fullCat
		if eventDate==None:
			# we want the whole catalog (no min date):
			#print tempCat[-1][0]
			eventDate=tempCat[0][0]
		
		#nEventsSinceMS=0
		for row in tempCat:
			# rotate each element into our aftershock axis, is it in the ellipse?
			#print row[0]
			# 2009-06-04 yoder: skipping nEvents since mainshock... (note, at this time we're already excluding up to and including MS)
			#nEventsSinceMS+=1
			#if nEventsSinceMS<=skipNevents: continue
			# or some period of time-time:
			#if row[0]<=eventDate: continue
			if maxDate>=eventDate and row[0]<=eventDate+datetime.timedelta(seconds=skipSeconds): continue		# 1_day = 86400 seconds
			if maxDate<=eventDate and row[0]>=eventDate-datetime.timedelta(seconds=skipSeconds): continue		# 22 july 2009 yoder: facilitat reverse RB
			
			# remove the .01 day immediately after the event from the data-set:
			#if row[0]<=eventDate+datetime.timedelta(seconds=864): continue
			#if row[0]<=eventDate+datetime.timedelta(seconds=1600): continue
			if maxDate>=eventDate and row[0]>maxDate: break
			if maxDate<=eventDate and row[0]<maxDate: break		# 22 july 2009 yoder: facilitat reverse RB
			
			newVec=self.faultTransform(row[2], row[1], clat, clon, theta, ra, rb)
			#
			# is the rotated vector in our ellipse?
			if abs(newVec[0])>ra: continue
			Y=ellipseY(newVec[0], ra, rb)
			if abs(newVec[1])>Y: continue
			# dtm, lat, lon, mag, tX, tY 		(note this is like y,x, x`, y` for the space coordinates).
			self.shockCat+=[[row[0], row[1], row[2], row[3], newVec[0], newVec[1]]]	
		#
			
	def setFullCat(self, catFname=None, minmag=1.25):
		if catFname==None: catFname=self.catname
		self.fullCat=[]
		
		#print "fullcat catname: %s" % catFname
		f=open(catFname)
		nrows=0
		for row in f:
			if row=='\n' or row[0]=='#': continue
			#if nrows==0: print row
			nrows+=1
			delim=' '
			if '\t' in row: delim='\t'
			#
			rowl=row.split(delim)
			# 2008 12 2 -117.318 35.97 3.06
			mag=float(rowl[4])
			if mag<minmag: continue
			delim=rowl[0][4]
			thisDtm=datetimeFromStrings(rowl[0], rowl[1], delim)
			#
			# skip (probable) duplicate entries.
			if nrows>1:
				if thisDtm==self.fullCat[-1][0]: continue
				
			# dtm, lat, lon, mag
			self.fullCat+=[[thisDtm, float(rowl[2]), float(rowl[3]), mag]]
			#print 'catrow: %s' % self.fullCat[-1]
		print "len of fullcat: %d, %s" % (len(self.fullCat), self.fullCat[0])
	
	def getMainEvent(self, cat=None):
		# return catalog row of max magnitude (epicenter location (more or less)). note, by default we use ths shock-cat because it will
		# be faster than using the fullCat AND, the fullCat is likely to have other large earthquakes.
		if cat==None: cat=self.shockCat
		if len(cat)==0: cat=self.fullCat
		maxMag=cat[0][3]
		maxIndex=0
		for i in xrange(len(cat)):
			#print i, maxMag, maxIndex, cat[i][3]
			if cat[i][3]>maxMag:
				maxIndex=i
				maxMag=cat[i][3]
		#
		return cat[maxIndex] + [maxIndex]
			
	def fitAftershocksToOmori(self, minmag=1.2, p=None):
		# in any case, this current form seems to work pretty well so long as we guess the initial prams within an order of magnitude or so
		# for example: (50, .001, 1.5)
		#
		# p=scipy.array([0,1])
		# plsq=spo.leastsq(linRes, p, args=(scipy.array(y), scipy.array(x)), full_output=0, maxfev=20000)
		# omori's law: n(t) = k/(c+t)**p --> p[0]/(p[1]+t)**p[2]
		#p=[1.0/4000, 1/1000, 1]		# on the first day, we get 12 events; figure it takes a few hours to set up; this should give us reasonalbe starting prams.
		if p==None or len(p)<3:
			p=[50, .02, 1.5]	# in days...
			
		x=[]
		y=[]
		nthEvent=0
		for rw in self.shockCat:
			if rw[3]<minmag: continue
			#print "use %s" % x
			x+=[datetimeToFloat(rw[0])-datetimeToFloat(self.eventDtTime)]
			y+=[nthEvent]
			nthEvent+=1
		
		plsq=spo.leastsq(omoriCumRes, p, args=(scipy.array(y), scipy.array(x)), full_output=0, maxfev=200000)
		p=plsq[0]
		for i in xrange(1000):
			# loop pseudo-recursively to converge:
			plsq=spo.leastsq(omoriCumRes, p, args=(scipy.array(y), scipy.array(x)), full_output=0, maxfev=200000)
			p=plsq[0]
		#print plsq
		#return plsq[0]
		yfit=[]
		chiSqr=0
		ndof=-3
		
		for X in x:
			yfit+=[fomoriCum(X,p)]
			#
			chiSqr+=(yfit[ndof+3]-y[ndof+3])
			ndof+=1
		chiSqr/=ndof
		
		print p[0], p[1], p[2], chiSqr
		return [p.tolist()+[chiSqr], [x,y,yfit]]

	def plotfitAftershocksToOmori(self, minmag=1.2, p=None):
		prams=self.fitAftershocksToOmori(minmag, p)
		p=prams[0]
		x=prams[1][0]
		y=prams[1][1]
		yfit=prams[1][2]
		
		plt.figure(1)
		plt.plot(x,y,'.-', label='data')
		plt.plot(x,yfit, '.-', label='fit')
		plt.legend(loc='lower right')
		plt.show()
				
		return [p, [x,y,yfit]]
	
	def xyPlotShocks(self):
		#lats=map(operator.itemgetter(1), rbi.shockCat)
		plt.figure(0)
		plt.clf()
		plt.plot(map(operator.itemgetter(2), self.shockCat), map(operator.itemgetter(1), self.shockCat), '.')
		plt.show()
	
	def xyPlotFull(self):
		#lats=map(operator.itemgetter(1), rbi.shockCat)
		plt.figure(0)
		plt.clf()
		plt.plot(map(operator.itemgetter(2), self.fullCat), map(operator.itemgetter(1), self.fullCat), '.')
		plt.show()
	
	def getLatLonRange(self, cat=None, latloncols=[1,2]):
		if cat==None: cat=self.fullCat
		if latloncols==None: latloncols=[1,2]	# latitude, lon cols of catalog (order is lat, lon).
		#
		minLat=cat[0][latloncols[0]]
		maxLat=cat[0][latloncols[0]]
		minLon=cat[0][latloncols[1]]
		maxLon=cat[0][latloncols[1]]
		#
		for rw in cat:
			thisLat=rw[latloncols[0]]
			thisLon=rw[latloncols[1]]
			#
			if thisLat>maxLat: maxLat=thisLat
			if thisLat<minLat: minLat=thisLat
			if thisLon>maxLon: maxLon=thisLon
			if thisLon<minLon: minLon=thisLon
		#
		return [[minLat, minLon], [maxLat, maxLon]]
	
	def xyPlotCatsMap(self, doShow=True, doSave=False, saveName='catalogPlot.png', epicenter=None, legendLoc='upper left'):
		
		if epicenter==None: epicenter=[self.tLon, self.tLat]
		fcat=[]
		scat=[]
		for rw in self.shockCat:
			scat+=[rw[0:4]]
		for rw in self.fullCat:
			if rw not in scat: fcat+=[rw]
		#return [scat, fcat]
	
		f0=plt.figure(0)	
		plt.clf()
		#		
		#set up map:
		llr=self.getLatLonRange()	# latLonRange
		cntr=[float(llr[0][0])+(llr[1][0]-float(llr[0][0]))/2.0, float(llr[0][1])+(llr[1][1]-float(llr[0][1]))/2.0]
		catmap=Basemap(llcrnrlon=llr[0][1], llcrnrlat=llr[0][0], urcrnrlon=llr[1][1], urcrnrlat=llr[1][0], resolution ='l', projection='tmerc', lon_0=cntr[1], lat_0=cntr[0])
		canvas=FigureCanvas(f0)
		catmap.ax=f0.add_axes([0,0,1,1])
		f0.set_figsize_inches((8/catmap.aspect,8.))
		#
		catmap.drawcoastlines(color='gray')
		catmap.drawcountries(color='gray')
		catmap.fillcontinents(color='beige')
		xfull, yfull=catmap(map(operator.itemgetter(2), fcat), map(operator.itemgetter(1), fcat))
		xshock, yshock=catmap(map(operator.itemgetter(2), scat), map(operator.itemgetter(1), scat))
		epx, epy=catmap(epicenter[0], epicenter[1])
		catmap.plot(xfull, yfull, 'g+', label='Full Catalog')
		catmap.plot(xshock, yshock, 'b.', label='Aftershock zone')
		catmap.plot(epx, epy, 'ro')
		
		canvas.print_figure(saveName)
		
		#
		#ax=plt.gca()
		el = Ellipse((self.tLon, self.tLat), 2.0*self.tA, 2.0*self.tB, -self.tTheta, facecolor='b', alpha=0.4)
		#catmap.ax.add_artist(el)
		#ax.add_artist(el)
		#
		#plt.plot(map(operator.itemgetter(2), self.fullCat), map(operator.itemgetter(1), self.fullCat), '+')
		#plt.plot(map(operator.itemgetter(2), self.shockCat), map(operator.itemgetter(1), self.shockCat), '.')
		#plt.plot(map(operator.itemgetter(2), fcat), map(operator.itemgetter(1), fcat), '+', label='Full Catalog')
		#plt.plot(map(operator.itemgetter(2), scat), map(operator.itemgetter(1), scat), '.', label='Aftershock zone')
		#plt.plot([epicenter[0]], [epicenter[1]], 'ro', label='epicenter')
		plt.legend(loc=legendLoc, numpoints=1)
		if doSave: plt.savefig(saveName)
		if doShow: plt.show()

	def xyPlotCatalogs(self, doShow=True, doSave=False, saveName='catalogPlot.png', epicenter=None, legendLoc='upper left'):
		
		if epicenter==None: epicenter=[self.tLon, self.tLat]
		fcat=[]
		scat=[]
		for rw in self.shockCat:
			scat+=[rw[0:4]]
		for rw in self.fullCat:
			if rw not in scat: fcat+=[rw]
		#return [scat, fcat]
	
		plt.figure(0)	
		plt.clf()
		#
		ax=plt.gca()
		#el = Ellipse((-72.533,18.457), 1.25, .25, 15, facecolor='r', alpha=0.5)
		el = Ellipse((self.tLon, self.tLat), 2.0*self.tA, 2.0*self.tB, -self.tTheta, facecolor='b', alpha=0.4)
		ax.add_artist(el)
		#
		#plt.plot(map(operator.itemgetter(2), self.fullCat), map(operator.itemgetter(1), self.fullCat), '+')
		#plt.plot(map(operator.itemgetter(2), self.shockCat), map(operator.itemgetter(1), self.shockCat), '.')
		plt.plot(map(operator.itemgetter(2), fcat), map(operator.itemgetter(1), fcat), '+', label='Full Catalog')
		plt.plot(map(operator.itemgetter(2), scat), map(operator.itemgetter(1), scat), '.', label='Aftershock zone')
		plt.plot([epicenter[0]], [epicenter[1]], 'ro', label='epicenter')
		plt.legend(loc=legendLoc, numpoints=1)
		if doSave: plt.savefig(saveName)
		if doShow: plt.show()	
	
	def plotMagsIntervals(self, doShow=True, doSave=False, saveName='catalogMagsInts', pltTitle='Catalog Seismicity', nTime=False, minmag=1.5):
		fnameFull="%sFull.png"
		fnameShock="%sShock.png"
		
		cats=[self.fullCat, self.shockCat]
		fignums=xrange(len(cats))
		#
		fsize=20
		#minmag=1.1	# we should get this from the data
		for fnum in fignums:
			currcat=cats[fnum]
			
			plt.figure(fnum)
			plt.clf()
			ax0=plt.axes([.1,.1,.85, .20])
			
			#
			#plt.title("Magnitudes")
			if nTime:
				plt.xlabel("Number of Events (n)", fontsize=fsize)
			else:
				plt.xlabel("date", fontsize=fsize)
			
			plt.ylabel("Mags", fontsize=fsize)
			#
			#ax1=plt.axes([.1, .4, .85, .50], sharex=ax0)
			ax1=plt.axes([.1, .55, .85, .35], sharex=ax0)
			#plt.title(pltTitle, fontsize=fsize)
			if currcat==self.fullCat: plt.title("Parkfield 4x4 catalog", fontsize=fsize)
			if currcat==self.shockCat: plt.title("Parkfield aftershock catalog", fontsize=fsize)
			plt.ylabel("intervals (days)", fontsize=fsize)
			plt.xlabel("")
			#
			intervals=[]
			intbars=[]
			magbars=[]
			for i in xrange(len(currcat)):
				thisx=currcat[i][0]
				magbars+=[[thisx,minmag], [thisx, currcat[i][3]], [thisx,minmag]]
				if i==0: continue
				intbars+=[[thisx,0], [thisx, date2num(currcat[i][0])-date2num(currcat[i-1][0])], [thisx,0]]
			
			#zeros=[0]
			#for i in xrange(1, len(currcat)):
			#	intervals+=[date2num(currcat[i][0])-date2num(currcat[i-1][0])]
			#	#zeros+=[0]
			#
			#ax0([datetime.datetime(2000,01,01), datetime.datetime(2010,12,31)])
			
			#ax0.plot_date(map(operator.itemgetter(0), self.fullCat), map(operator.itemgetter(3), self.fullCat), '-')
			#ax1.plot_date(map(operator.itemgetter(0), currcat[1:]), intervals, '-')
			
			if nTime==False:
				ax0.plot_date(map(operator.itemgetter(0), magbars), map(operator.itemgetter(1), magbars), '-')
				ax1.plot_date(map(operator.itemgetter(0), intbars), map(operator.itemgetter(1), intbars), '-')
				ax1.set_yscale('log')
				a = plt.gca()
				#a.set_xlim([currcat[0][0]-datetime.timedelta(days=20), currcat[-1][0]+datetime.timedelta(days=20)])
				
			if nTime==True:
				y1=map(operator.itemgetter(1), magbars)
				y2=map(operator.itemgetter(1), intbars)
				ax0.plot(xrange(len(y1)), y1, '-')
				ax1.plot(xrange(len(y2)), y2, '-')
				ax1.set_yscale('log')
				a = plt.gca()
				a.set_xlim([0, int(1.1*len(y2))])		
			#
			if doSave:
				if currcat==self.fullCat: catTag='-full'
				if currcat==self.shockCat: catTag='-shock'
				plt.savefig('%s-%s.png' % (saveName, catTag))
		if doShow: plt.show()
		
		
	def GRshock(self, doShow=True, fname='GRdist.png', plotTitle="Aftershock Region Magnitude Distribution"):
		# [[evDateTime, lat, lon, mag, a, b], [row1], ...]
		mags=map(operator.itemgetter(3), self.shockCat)
		return self.GRdist(mags, doShow, fname, plotTitle)
		
	def GRfullcat(self, doShow=True, fname='GRdist.png', plotTitle="Full Catalog Magnitude Distribution"):
		mags=map(operator.itemgetter(3), self.fullCat)
		return self.GRdist(mags, doShow, fname, plotTitle)
	
	def GRdist(self, mags, doShow=True, fname='GRdist.png', plotTitle="Magnitude Distribution"):
		# cat: a 1D array of magnitudes
		mags.sort()
		# get rid of biggest event (probably a large off-GR earthquake):
		mags.pop()
		#mags.reverse()
		#print mags
		#print len(mags)
		
		if doShow==True or fname!=None:
			# make a plot and show and/or save
			#Y=range(1, len(mags)+1)
			Y=yp.frange(1, len(mags), -1)
			#print Y
			#print len(Y)
			plt.figure(0)
			plt.clf()
			plt.semilogy(mags, Y, '.-')
			plt.xlabel("Magnitude, m")
			plt.ylabel("Number of Events, n")
			plt.title(plotTitle)
			if fname!=None: plt.savefig(fname)
			if doShow: plt.show()
		
		return mags
		#
		
	def fitOmoriRange(self, minmag=1.25, maxmag=4.5, dmag=.1, p=None):
		pramses=[[], [], [], [], []]

		fname="parkfieldOmoriFit.dat"
		f=open(fname,'w')
		f.write("#k\tc\tp\trChiSqr\n")
		f.close()
		
		while minmag<=maxmag:
			print "fitting for mag %f." % minmag
			prams=self.fitAftershocksToOmori(minmag,p)[0]
			print prams
			f=open(fname,'a')
			f.write("%f\t%f\t%f\t%f\t%f\n" % (minmag, prams[0], prams[1], prams[2], prams[3]))
			f.close()
			pramses[0]+=[minmag]
			pramses[1]+=[prams[0]]
			pramses[2]+=[prams[1]]
			pramses[3]+=[prams[2]]
			pramses[4]+=[prams[3]]
			
			minmag+=dmag
		return pramses
		
	def faultTransform(self, x, y, Lat=tLat, Lon=tLon, theta=tTheta, A=tA, B=tB):
		# x,y to transform via blah, blah.
		#
		theta=deg2rad(float(theta))
		xprime = (x-Lon)*cos(theta) - (y-Lat)*sin(theta)
		yprime = (x-Lon)*sin(theta) + (y-Lat)*cos(theta)
		
		return [xprime, yprime]
	#
	def getBigShocks(self, minmag=1.5, bigmag=5.0, cat=None):
		if cat==None: cat=self.shockCat
		rnum=0
		bigShocks=[]
		for rw in cat:
			if rw[3]<minmag: continue	# don't count rnum...
			if rw[3]>=bigmag: bigShocks+=[[rnum] + rw]
			rnum+=1
		return bigShocks
		
	def getIntervalRatios(self, minmag=3.0, windowLen=10, cat0=None, deltaipos=1):
		# deltaipos: event resolution; how many events to advance each step.
		if cat0==None: cat=self.shockCat
		#
		# problem: if you give this function a catalog with magnitudes below the min-mag, the getRBintervals() function can break (it has a bit
		# that strips out sub-m0 events). SO, do that here.
		cat=[]
		#icount=0
		for rw in cat0:
			#print icount, rw[3], minmag, rw[3]>=minmag
			#icount+=1
			if rw[3]>=minmag: cat+=[rw]
		#
		ipos=0
		rbRatios=[]
		#print ipos, len(cat)-windowLen, minmag, len(cat), len(cat0)
		cat0=None
		while ipos<(len(cat)-windowLen):
			thisRatios=self.getRBintervals(minmag, cat[ipos:(ipos+windowLen)])	# note: if we want to do backwards RB, using cat[i:i+l].reverse() is a nice trick.
			rbRatios+=[[ipos+windowLen, cat[(ipos+windowLen)][0], float(thisRatios[0][2][-1])/float(thisRatios[1][2][-1])]]
			#ipos+=deltaipos
			ipos+=1
		# return [[n, dt, r]]
		thisRatios=None
		return rbRatios

	def plotIntervalRatiosAx(self, minmag=3.0, windowLen=10, cat0=None, hitThreshold=1.0, bigmag=5.0, thisAx=None, ratios=None, deltaipos=1, avlen=1, mainEV=None):
		# avlen=10
		# eventually, this will probably be the sole version of this function....
		# this creates a plot of interval ratios in/on a specified axis. use this function by itself or to make complex figures.
		#
		if thisAx==None:
			f0=plt.figure()
			thisAx=f0.gca()
		#	
		legLoc='upper left'
		eventName="Event RB ratios"
		if cat0==None: cat0=self.shockCat
		# getIntervalRatios(self, minmag=3.0, windowLen=10, cat0=None, deltaipos=1):
		if ratios==None: ratios=self.getIntervalRatios(minmag, windowLen, cat0, deltaipos)
		fdts=[]
		for rw in ratios:
			 fdts+=[rw[1].toordinal() + float(rw[1].hour)/24 + rw[1].minute/(24*60) + rw[1].second/(24*3600) + rw[1].microsecond/(24*3600000000)]
		plaindts=map(operator.itemgetter(1), ratios)
		
		if mainEV==None: mainEV=self.getMainEvent(cat0)
		eventFloatDt=yp.datetimeToFloat(mainEV[0])
		#
		#f=plt.figure(fignum)
		#theseFigs+=[f]
		thisAx.cla()	
		thisAx.axvline(x=eventFloatDt, color='c', lw=3, label='mainshock' )
		#plt.fill_between(fdts, scipy.ones(len(fdts),int), map(operator.itemgetter(2), hratios), color='b', where=scipy.array([val>=1 for val in map(operator.itemgetter(2), hratios)]))
		#plt.fill_between(fdts, scipy.ones(len(fdts),int), map(operator.itemgetter(2), hratios), color='r', where=scipy.array([val<=1 for val in map(operator.itemgetter(2), hratios)]))
		reload(yp)
		#
		# ultimately, averaging the actual ratios is not meaningful; we need to average the logs.
		# so, <log(x)> = (1/N)(x1+x2+...xn) = log(Prod(x_i)**(1/N))
		#ploty=yp.averageOver(map(operator.itemgetter(2), ratios), avlen)
		ploty=yp.logaverageOver(map(operator.itemgetter(2), ratios), avlen)
		#
		#thisAx.fill_between(fdts, hitThreshold*scipy.ones(len(fdts),int), ploty[0], color='b', where=scipy.array([val>=hitThreshold for val in ploty[0]]))
		#thisAx.fill_between(fdts, hitThreshold*scipy.ones(len(fdts),int), ploty[0], color='r', where=scipy.array([val<=hitThreshold for val in ploty[0]]))
		thisAx.fill_between(plaindts, hitThreshold*scipy.ones(len(fdts),int), ploty[0], color='b', where=scipy.array([val>=hitThreshold for val in ploty[0]]))
		thisAx.fill_between(plaindts, hitThreshold*scipy.ones(len(fdts),int), ploty[0], color='r', where=scipy.array([val<=hitThreshold for val in ploty[0]]))
		#
		maxy=math.log10(max(ploty[0]))
		#miny=min(ploty[0])
		#print "maxy, miny: %d, %d" % (maxy, miny)
		
		# note: we don't set the y-log scale in these "fill()" commands. we can do that with axis.set_yscale('log') i think.
		# we achieve this by doing semilogy() plots below.
		fg=plt.gcf()
		#plt.title("%s rupture area, time-time, wLen=%d" % (eventName, windowLen))
		#plt.xlabel('time')
		#plt.ylabel('$r=N_{rb-long} / N_{rb-short}$')
		thisAx.axvline(x=eventFloatDt)
		nbigshocks=0
		bigShocks=self.getBigShocks(minmag, bigmag, cat0)
		for rw in bigShocks:
			if nbigshocks==0:
				#thisAx.axvline(x=yp.datetimeToFloat(rw[1]), color='g', label='m > %f' % bigmag)
				thisAx.axvline(x=rw[1], color='g', label='m > %f' % bigmag)
				nbigshocks+=1
			else:
				#thisAx.axvline(x=yp.datetimeToFloat(rw[1]), color='g')
				#thisAx.plot([yp.datetimeToFloat(rw[1])], rw[4], '*')
				thisAx.axvline(x=rw[1], color='g')
				thisAx.plot([rw[1]], rw[4]*maxy, '*')
			
		#thisAx.semilogy([eventFloatDt], [1], 'r^', ms=10)
		thisAx.semilogy([mainEV[0]], [1], 'r^', ms=10)
		
		thisAx.axhline(y=1, color='k')
		thisAx.legend(loc=legLoc, numpoints=2)
		#ax=plt.gca()
		#fg=plt.gcf()
	#	thisAx.xaxis.set_major_formatter(dates.DateFormatter('%Y-%b-%d'))
		thisAx.set_ylim([.1,10])
	#	fg.autofmt_xdate()
		#plt.savefig('images/%sRuptureTimeTime-Wlen%d-mc%d.png' % (eventName, wLen, int(10*minMag)))
		#plt.show()		
	
	def plotIntervalRatios(self, minmag=3.0, windowLen=10, cat0=None, hitThreshold=1.0, bigmag=5.0, fignum=0, ratios=None, deltaipos=1, avlen=1):
		# avlen=10
		#			
		legLoc='upper left'
		eventName="Event RB ratios"
		if cat0==None: cat0=self.shockCat
		# getIntervalRatios(self, minmag=3.0, windowLen=10, cat0=None, deltaipos=1):
		if ratios==None: ratios=self.getIntervalRatios(minmag, windowLen, cat0, deltaipos)
		fdts=[]
		for rw in ratios:
			 fdts+=[rw[1].toordinal() + float(rw[1].hour)/24 + rw[1].minute/(24*60) + rw[1].second/(24*3600) + rw[1].microsecond/(24*3600000000)]
		mainEv=self.getMainEvent(cat0)
		eventFloatDt=yp.datetimeToFloat(mainEv[0])
		#
		f=plt.figure(fignum)
		#theseFigs+=[f]
		plt.clf()	
		plt.axvline(x=eventFloatDt, color='c', lw=3, label='mainshock' )
		#plt.fill_between(fdts, scipy.ones(len(fdts),int), map(operator.itemgetter(2), hratios), color='b', where=scipy.array([val>=1 for val in map(operator.itemgetter(2), hratios)]))
		#plt.fill_between(fdts, scipy.ones(len(fdts),int), map(operator.itemgetter(2), hratios), color='r', where=scipy.array([val<=1 for val in map(operator.itemgetter(2), hratios)]))
		reload(yp)
		
		ploty=yp.averageOver(map(operator.itemgetter(2), ratios), avlen)
		plt.fill_between(fdts, hitThreshold*scipy.ones(len(fdts),int), ploty[0], color='b', where=scipy.array([val>=hitThreshold for val in ploty[0]]))
		plt.fill_between(fdts, hitThreshold*scipy.ones(len(fdts),int), ploty[0], color='r', where=scipy.array([val<=hitThreshold for val in ploty[0]]))
	
		# note: we don't set the y-log scale in these "fill()" commands. we can do that with axis.set_yscale('log') i think.
		# we achieve this by doing semilogy() plots below.
		plt.title("%s rupture area, time-time, wLen=%d" % (eventName, windowLen))
		plt.xlabel('time')
		plt.ylabel('$r=N_{rb-long} / N_{rb-short}$')
		plt.axvline(x=eventFloatDt)
		nbigshocks=0
		bigShocks=self.getBigShocks(minmag, bigmag, cat0)
		for rw in bigShocks:
			if nbigshocks==0:
				plt.axvline(x=yp.datetimeToFloat(rw[1]), color='g', label='m > %f' % bigmag)
				nbigshocks+=1
			else:
				plt.axvline(x=yp.datetimeToFloat(rw[1]), color='g')
				plt.plot([yp.datetimeToFloat(rw[1])], rw[4], '*')
			
		plt.semilogy([eventFloatDt], [1], 'r^', ms=10)
		
		plt.axhline(y=1, color='k')
		plt.legend(loc=legLoc, numpoints=2)
		ax=plt.gca()
		fg=plt.gcf()
		ax.xaxis.set_major_formatter(dates.DateFormatter('%Y-%b-%d'))
		ax.set_ylim([.1,10])
		fg.autofmt_xdate()
		#plt.savefig('images/%sRuptureTimeTime-Wlen%d-mc%d.png' % (eventName, wLen, int(10*minMag)))
		plt.show()
	
	def getMagSubset(self, mag=5.0, cat=None):
		# get a subset of a catalog with m>mmag
		if cat==None: cat=self.shockCat
		outCat=[]
		for rw in cat:
			if rw[3]>=mag: outCat+=[rw]
		return outCat
	
	def getLargeAftershocks(self, mag=5.0, cat=None):
		if cat==None: cat=self.shockCat
		mainShock=self.getMainEvent()
		outcat=[]
		for rw in cat:
			if rw[0]>mainShock[0] and rw[3]>=mag: outcat+=[rw]
		mainShock=None
		
		return outcat
		#
	def getEarthquakeRatioScore(self, ratios=None, earthquakes=None):
		# the idea here is to get the "current" value of r(t) when an earthquake occurs. was the event predicted?
		# so, for each earthquake in earthquakes[], what was the most recent value of r, in ratios[].
		# if ratios and earthquakes are not provided, use getIntervalRatios(self), and m5.0 events from self.shockCat, respectively.
		#
		if ratios==None: ratios=self.getIntervalRatios()	#[n, date, r]
		if earthquakes==None: earthquakes=self.getMagSubset(5.0, self.shockCat)		#[date, lat, lon, mag, a, b]
		returnQuakes=[]
		#
		# now, assign a ratio value to each earthquake:
		#
		#for eqrow in earthquakes:
		for i in xrange(len(earthquakes)):
			#if eqrow[0]<ratios[0][1]: continue	# there is a winLen lag; we dont' have a forecast yet.
			# add the r value of the closest ratios entry to each earthquake:
			for irat in xrange(1, len(ratios)):
				if earthquakes[i][0]==ratios[irat][1]:
				#if ratios[irat][1]>=earthquakes[i][0]:
					returnQuakes+=[[earthquakes[i]+[ratios[irat-1][2]]]]
					continue
		#		
		return returnQuakes
		
	# science:
	def getRBintervals(self, minmag=1.0, useCat=None):
		# default record-breaking. walk forward in the shockCat; look for the largest/smallest intervals between events:
		cat=[]
		# print "shockcat[0]: %s, %s, %s, %s" % (self.shockCat[0][0], self.shockCat[1][0], self.shockCat[2][0], self.shockCat[3][0])
		if useCat==None: useCat=self.shockCat
		#for row in self.shockCat:
		for row in useCat:
			#print "rw3, minmag: %s, %s" % (row[3], minmag)
			if row[3]>=minmag: cat+=[row]
		#
		#print "catlen: %d" % len(cat)
		biggest=abs(datetimeToFloat(cat[1][0])-datetimeToFloat(cat[0][0]))	# interval between second and first events.
		#print "biggest interval: %f" % biggest
		smallest=biggest
		nbigger=1
		nsmaller=1
		
		# these arrays will be returned plot-ready:
		biggers=[[cat[1][0]], [biggest], [nbigger], [1]]		# [[date], [bigInts], [NRB_big], [i]]
		smallers=[[cat[1][0]],[smallest], [nsmaller], [1]]
		
		#prevDtmBig=cat[0][0]
		#prevDtmSmall=cat[0][0]
		for i in xrange(1,len(cat)):
			#dT=(cat[i][0]-cat[i-1][0]).seconds
			dT=abs(datetimeToFloat(cat[i][0])-datetimeToFloat(cat[i-1][0]))		# aka, the interval between the current and previous event...
			# note: i guess this WILL work backwards. intrinsically, it will work in reverse - intervals will be negative. we could
			# use this as is, but it probably makes sense to use the absolute value. after all, we want the magnitude of the interval:
			#
			#if cat[i][0]<cat[i-1][0]:
			#	print cat[i]
			#	a=input("type something")
			
			#print cat[i][0], cat[i-1][0], cat[i][3], dT
			if dT>biggest:
				nbigger+=1
				biggers[0]+=[cat[i][0]]		# date record was broken (datetime object)
				biggers[1]+=[dT]				# interval since last event (RBID)
				biggers[2]+=[nbigger]		# number of records broken (NRB)
				biggers[3]+=[i]				# record broken on i'th earthquake since mainshock (natural time).	[i+1] ??
				biggest=dT
				#prevDtmBig=cat[i][0]
				#
			if dT<smallest:
				nsmaller+=1
				smallers[0]+=[cat[i][0]]
				smallers[1]+=[dT]
				smallers[2]+=[nsmaller]
				smallers[3]+=[i]				# record broken on i'th earthquake since mainshock.
				smallest=dT
				#prevDtmSmall=cat[i][0]
				#
			#
		#
		#print biggers[0]
		#print biggers[1]
		
		cat=None
		return[biggers, smallers]
	
	def plotRBintervalSet(self, minmag=1.0, maxmag=5.0, dmag=.1, outdir='images/parkfield/'):
		plt.clf()
		plt.cla()
		
		curmag=minmag
		intervalSet=[]
		mags=[]
		if outdir[-1]!='/': outdir="%s/" % outdir
		
		print "begin plotRBintervalSet() {%f}" % curmag
		
		while curmag<=maxmag:
			curRecords=self.getRBintervals(curmag)	# curr records -> [[bigs], [smalls]]
			#bigsX=curRecords[0][0]
			#bigsY=curRecorss[0][1]
			# big records from currRecords
			intervalSet+=[[curRecords[0][0], curRecords[0][1], curRecords[0][2], curRecords[0][3], curmag]]		# "biggest" records, [date occured, interval, Nth record, mag-bin]
			mags+=[curmag]
			curmag+=dmag
		
		print "intervalSet(s) assigned..."
			
		#print len(intervalSet)
		#print len(intervalSet[0][0])
		
		startTime=datetimeToFloat(self.eventDtTime)
		
		# for giggles, get all the intervals:
		Xev=[0]
		Yev=[0]
		XevLog=[1]
		NevLog=[1]
		YevLog=[1]
		for iev in xrange(1,len(self.shockCat)):
			if (datetimeToFloat(self.shockCat[iev][0]) - datetimeToFloat(self.shockCat[iev-1][0])==0) : continue
			#
			Xev+=[datetimeToFloat(self.shockCat[iev][0]) - datetimeToFloat(self.shockCat[0][0])]
			Yev+=[datetimeToFloat(self.shockCat[iev][0]) - datetimeToFloat(self.shockCat[iev-1][0])]
			#
			#if abs(log10(datetimeToFloat(self.shockCat[iev][0]) - datetimeToFloat(self.shockCat[iev-1][0])))>3: continue
			XevLog+=[log10(datetimeToFloat(self.shockCat[iev][0]) - datetimeToFloat(self.shockCat[0][0]))]
			NevLog+=[log10(iev)+1]
			YevLog+=[log10(datetimeToFloat(self.shockCat[iev][0]) - datetimeToFloat(self.shockCat[iev-1][0]))]
		
		fignum=0
		plt.figure(fignum)
		plt.clf()
		fignum+=1
		plt.loglog(Xev, Yev, '.')
		plt.title("All Intervals")
		plt.xlabel("days since earthquake")
		plt.ylabel("interval")
		
		plt.figure(fignum)
		plt.clf()
		fignum+=1
		plt.loglog(Xev, Yev, '.')
		for pset in intervalSet:
			#plt.plot_date(pset[0], pset[1], '.-', label=str(minmag))
			x=[]
			#x2=[]
			for elem in pset[0]:
				x+=[datetimeToFloat(elem)-startTime]
				#x2+=[datetimeToFloat(elem)]
			plt.loglog(x, pset[1], '-', label=str(pset[3]))
			#plt.loglog(x2, pset[1], '-')
			
		#
		plt.title("RB interval")
		plt.xlabel("days since earthquake")
		plt.ylabel("interval (days)")
		plt.savefig("%sRBintervals.pdf" % outdir)
		#plt.legend(loc='upper left')
		
		#print "saveFig: %sRBintervals.pdf" % outdir
		
		plt.figure(fignum)
		plt.clf()
		fignum+=1
		for pset in intervalSet:
			#print len(pset)
			#plt.plot_date(pset[0], pset[2], '.-', label=str(minmag))
			x=[]
			for elem in pset[0]:
				x+=[datetimeToFloat(elem)-startTime]
			plt.loglog(x, pset[2], '.-', label=str(pset[4]))
		
		plt.title("Number of New RB Intervals")
		plt.xlabel("days since earthquake")
		plt.ylabel("Number of Broken Records")
		plt.savefig("%snewRBintervals.pdf" % outdir)
		#plt.legend(loc='upper left')
		
		# and let's do some data fitting:
		# for now, fit t>=10**-1.5 in the Nrecords plot, maybe interval>10**-2 for RBinterval??
		# alternatively, always trim the first (maybe first 2) data points.
		#
		# we want the log/log fit, so we can eithe rconstruct a log-based error or we can take the log-log
		# of the data and do a linear fit.
		#
		fitPlotsInt=[]
		fitPlotsN=[]
		slopes=[]
		#print YevLog
		#print "lens: %d, %d" % (len(XevLog), len(YevLog))
		fbooga=open("%sXevLog.dat" % outdir, 'w')
		for booga in xrange(len(XevLog)):
			fbooga.write("%f\t%f\t%f\t%f\n" %(XevLog[booga], YevLog[booga], Xev[booga], Yev[booga]))
			#fbooga.write("%f\t%f\n" %(XevLog[booga], YevLog[booga]))
		fbooga.close()
		
		plt.figure(fignum)
		plt.clf()
		fignum+=1
		plt.plot(XevLog, YevLog, '.')
		for pset in intervalSet:
			x=[]
			y=[]
			yfit=[]
			for i in xrange(len(pset[0])):
				#if i<1: continue		# (skip first element)
				#print "dtm arg: %s" % (datetimeToFloat(pset[0][i])-startTime)
				if log10(datetimeToFloat(pset[0][i])-startTime)<-1.0: continue
				x+=[log10(datetimeToFloat(pset[0][i])-startTime)]
				y+=[log10(pset[1][i])]
				
			#
			# now, we have a log-log set of one dataset. fit it to a line...
			#print "x: %s" % str(x)
			#print "y: %s" % str(y)
			
			p=scipy.array([0,1])
			plsq=spo.leastsq(linRes, p, args=(scipy.array(y), scipy.array(x)), full_output=0, maxfev=20000)	# (function to minimize, initial prams, argument-arrays, max-iterations)
			slopes+=[plsq[0][1]]
			for X in x:
				fitval=plsq[0][0] + X*plsq[0][1]
				yfit+=[fitval]
			plt.plot(x,yfit, '.-', label=str(pset[3]))
		plt.title("Fit to RB Intervals")
		plt.xlabel("log10(days since mainshock)")
		plt.ylabel("log10(Days Since Prev eq m>=m0)")
		plt.savefig("%sRBintervalFits.pdf" % outdir)
		#plt.legend(loc='upper left')
		#
		plt.figure(fignum)
		plt.clf()
		fignum+=1
		plt.title("Interval Slopes")
		plt.xlabel("mag")
		plt.ylabel("slope")
		plt.plot(mags, slopes, '.-')
		plt.savefig("%sRBintervalSlopes.pdf" % outdir)
		
		#############
		fitPlotsInt=[]
		fitPlotsN=[]
		slopes=[]
		plt.figure(fignum)
		plt.clf()
		fignum+=1
		for pset in intervalSet:
			x=[]
			y=[]
			yfit=[]
			for i in xrange(len(pset[0])):
				#if i<1: continue		# (skip first element)
				if log10(datetimeToFloat(pset[0][i])-startTime)<-1.0: continue
				x+=[log10(datetimeToFloat(pset[0][i])-startTime)]
				y+=[log10(pset[2][i])]
				
			#
			# now, we have a log-log set of one dataset. fit it to a line...
			#print "x: %s" % str(x)
			#print "y: %s" % str(y)
			p=scipy.array([0,1])
			plsq=spo.leastsq(linRes, p, args=(scipy.array(y), scipy.array(x)), full_output=0, maxfev=50000)	# (function to minimize, initial prams, argument-arrays, max-iterations)
			slopes+=[plsq[0][1]]
			for X in x:
				#print X
				fitval=plsq[0][0] + X*plsq[0][1]
				yfit+=[fitval]
			plt.plot(x,yfit, '.-', label=str(pset[4]))
		plt.title("fit to Number of New Records")
		plt.xlabel("log10(days since event)")
		plt.ylabel("log10(Nrecord breaking Intervals)")
		plt.savefig("%sRBfitNewRecords.pdf" % outdir)
		#plt.legend(loc='upper left')
		
		plt.figure(fignum)
		plt.clf()
		fignum+=1
		plt.title("Slopes Nrecords")
		plt.xlabel("mag")
		plt.ylabel("slope")
		plt.plot(mags, slopes, '.-')
		plt.savefig("%sRBslopesNewRecords.pdf" % outdir)
		
		###########################################################
		###########################################################
		print "natural time plots..."
		# natrual time plots:
		plt.figure(fignum)
		plt.clf()
		fignum+=1
		plt.loglog(xrange(1,len(Yev)+1), Yev, '.')	# all the events...
		#plt.semilogy(range(1,len(Yev)+1), Yev, '.')
		for pset in intervalSet:
			#print len(pset), pset[3]
			#print len(pset[2])
			#plt.plot_date(pset[0], pset[1], '.-', label=str(minmag))
			#
			plt.loglog(pset[3], pset[1], '+-', label=str(pset[4]))
			
			#plt.semilogy(pset[3], pset[1], '-', label=str(pset[4]))
		#
		plt.title("RB interval (natural time)")
		plt.xlabel("n events since mainshock")
		plt.ylabel("interval (days)")
		plt.savefig("%sRBintervalsNT.pdf" % outdir)
		#plt.legend(loc='upper left')
		#
		# the RB line does not appear to line up with the data-data. what are the values of the last elements?
		#print "elements: %d, %f, %f, %f" % (len(Yev)+1, Yev[-1], pset[3][-1], pset[1][-1])
		# output to text-file:
		fout=open('%sNRBNTdata.dat' % outdir, 'w')
		fout.write("#nthEvent\tinterval\n")
		for ii in xrange(len(Yev)):
			fout.write("%d\t%f\n" % (ii, Yev[ii]))
		fout.close()
		fout=open('%sNRBNTrecs.dat' % outdir, 'w')
		fout.write("#nEvents\tRBinterval\n")
		for ii in xrange(len(pset[3])):
			fout.write("%d\t%f\n" % (pset[3][ii], pset[1][ii]))
		fout.close()
		
		
		plt.figure(fignum)
		plt.clf()
		fignum+=1
		for pset in intervalSet:
			#print len(pset)
			#plt.plot_date(pset[0], pset[2], '.-', label=str(minmag))
			#
			plt.loglog(pset[3], pset[2], '.-', label=str(pset[4]))
			#plt.semilogy(pset[3], pset[2], '.-', label=str(pset[4]))
		
		plt.title("Number of RB Events (NT)")
		plt.xlabel("n events since mainshock")
		plt.ylabel("Number of Record Breaking Events")
		plt.legend(loc='lower right')
		plt.savefig("%snewRBintervalsNT.pdf" % outdir)
		
		# and let's do some data fitting:
		# for now, fit t>=10**-1.5 in the Nrecords plot, maybe interval>10**-2 for RBinterval??
		# alternatively, always trim the first (maybe first 2) data points.
		#
		# we want the log/log fit, so we can either construct a log-based error or we can take the log-log
		# of the data and do a linear fit.
		#
		fitPlotsInt=[]
		fitPlotsN=[]
		slopes=[]
		plt.figure(fignum)
		plt.clf()
		fignum+=1
		#plt.plot(NevLog, YevLog, '.')
		for pset in intervalSet:
			x=[]
			y=[]
			yfit=[]
			#print len(pset[0])
			for i in xrange(len(pset[0])):
				#if i<1: continue		# (skip first element)
				if log10(datetimeToFloat(pset[0][i])-startTime)<-1.0: continue
				x+=[log10(pset[3][i])]
				y+=[log10(pset[1][i])]
				
			#
			# now, we have a log-log set of one dataset. fit it to a line...
			#print "x: %s" % str(x)
			#print "y: %s" % str(y)
			p=scipy.array([0,1])
			plsq=spo.leastsq(linRes, p, args=(scipy.array(y), scipy.array(x)), full_output=0, maxfev=20000)	# (function to minimize, initial prams, argument-arrays, max-iterations)
			slopes+=[plsq[0][1]]
			for X in x:
				fitval=plsq[0][0] + X*plsq[0][1]
				yfit+=[fitval]
			#plt.plot(x,yfit, '.-', label=str(pset[4]))
		plt.title("Fit to RB Intervals (NT)")
		plt.xlabel("log10(nevents since mainshock)")
		plt.ylabel("log10(Days Since Prev eq m>=m0)")
		plt.savefig("%sRBintervalFitsNT.pdf" % outdir)
		#plt.legend(loc='upper left')
		#
		plt.figure(fignum)
		plt.clf()
		fignum+=1
		plt.title("Interval Slopes (NT)")
		plt.xlabel("mag")
		plt.ylabel("slope")
		plt.plot(mags, slopes, '.-')
		plt.savefig("%sRBintervalSlopesNT.pdf" % outdir)
		
		#############
		fitPlotsInt=[]
		fitPlotsN=[]
		slopes=[]
		plt.figure(fignum)
		plt.clf()
		fignum+=1
		for pset in intervalSet:
			x=[]
			y=[]
			yfit=[]
			for i in xrange(len(pset[0])):
				#if i<1: continue		# (skip first element)
				if log10(datetimeToFloat(pset[0][i])-startTime)<-1.0: continue
				x+=[log10(pset[3][i])]
				y+=[log10(pset[2][i])]
			#
			# now, we have a log-log set of one dataset. fit it to a line...
			#print "x: %s" % str(x)
			#print "y: %s" % str(y)
			p=scipy.array([0,1])
			plsq=spo.leastsq(linRes, p, args=(scipy.array(y), scipy.array(x)), full_output=0, maxfev=50000)	# (function to minimize, initial prams, argument-arrays, max-iterations)
			slopes+=[plsq[0][1]]
			for X in x:
				#print X
				fitval=plsq[0][0] + X*plsq[0][1]
				yfit+=[fitval]
			plt.plot(x,yfit, '.-', label=str(pset[4]))
		plt.title("fit to Number of New Records (NT)")
		plt.xlabel("log10(nevents since event)")
		plt.ylabel("log10(Nrecord breaking Intervals)")
		plt.savefig("%sRBfitNewRecordsNT.pdf" % outdir)
		#plt.legend(loc='upper left')
		
		plt.figure(fignum)
		plt.clf()
		fignum+=1
		plt.title("Slopes Nrecords (NT)")
		plt.xlabel("mag")
		plt.ylabel("slope")
		plt.plot(mags, slopes, '.-')
		plt.savefig("%sRBslopesNewRecordsNT.pdf" % outdir)
		#
		#plt.show()
	#####################	
	
	# utilities:
	def scatterShockCat(self):
		# create a scatter plot of the shock-catalog.
		vecs=[[],[],[],[]]
		for row in self.shockCat:
			vecs[0]+=[row[1]]
			vecs[1]+=[row[2]]
			vecs[2]+=[row[4]]
			vecs[3]+=[row[5]]
		#
		#print vecs[1]
		plt.figure(0)
		plt.plot(vecs[1], vecs[0], '.')
		
		plt.figure(1)
		plt.plot(vecs[2], vecs[3], '.')
		
		plt.show()

	#
def deg2rad(theta):
	return 2*pi*theta/360.0
	
def ellipseY(x, a, b):
	#print b, x, a
	return b*(1-x*x/(a*a))**.5

def datetimeFromStrings(strDt, strTm, dtdelim='/'):
	# looks like date[time].strptime(date_string, format) does this...
	if strTm=='': strTm='00:00:00.0'
	#
	ldt=strDt.split(dtdelim)
	ltm=strTm.split(':')
	
	lsecs=ltm[2].split('.')
	secs = long(lsecs[0])
	msecs=0
	if len(lsecs)>1:
		msecs = long(float("."+lsecs[1])*10**6)
	#
	#return datetime.datetime(long(ldt[0]), long(ldt[1]), long(ldt[2]), long(ltm[0]), long(ltm[1]), long(ltm[2]))
	return datetime.datetime(long(ldt[0]), long(ldt[1]), long(ldt[2]), long(ltm[0]), long(ltm[1]), secs, msecs)

def datetimeToFloat(dtm):
	# return number of days, including fractional bit.
	return float(dtm.toordinal()) + float(dtm.hour)/24.0 + float(dtm.minute)/(24.0*60.0) + float(dtm.second)/(24.0*60.0*60.0) + float(dtm.microsecond)/(24.0*60*60*10**6)

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


def doRBintervals1():
	rbi=intervalRecordBreaker()
	rbi.scatterShockCat()

#################
# fit functions:
#################
def linf(x,p):
	return (p[0] + x*p[1])
	
def fomori(t,p):
	return p[0]/(p[1]+t)**p[2]

def fomoriCum(t,p):
	# nominally, we like to convert all the parameters into the same functional form, aka a*x**b --> (10**a`)*(x**b).
	# however, since p[2] is both an exponential and an (anti)linear term (as per integration), we get a big fat mess if we do this.
	# visually, the convergance appears to be better with the first form, though paremeters must be chosen more carefully.
	
	# in any case, this current form seems to work pretty well so long as we guess the initial prams within an order of magnitude or so
	# for example: (50, .001, 1.5)
	p0=p[0]
	#p0=10**p[0]
	p1=p[1]
	p2=p[2]
	
	#foc= p[0]*((p[1]+t)**(1-p[2]) - p[1]**(1-p[2]))/(1-p[2])
	foc= p0*((p1+t)**(1-p2) - p1**(1-p2))/(1-p2)

	#foc= (10**p[0])*((10**p[1] + t)**(1-p[2]) - 10**(p[1]*(1-p[2])))/(1-10**p[2])
	#print p
	return foc

def linRes(p, y, x):
	# "linear residual"
	#print p, y, x
	err=y-linf(x,p)
	errs=scipy.array(err).tolist()	# once upon a time, "errs" was a class attribute; we set it here so we could sum for chi-sqr.

	return err

def omoriRes(p, y, t):
	err=y-fomori(t,p)
	return err

def omoriCumRes(p, y, t):
	err=y-fomoriCum(t,p)
	return err

def standardPFrecords():
	rbp=intervalRecordBreaker()
	rbp.plotRBintervalSet(1.5,3.5, .5)	

def standardSocalRecs():
	# def setNormalCat(self, catFname=None, minlat=32, maxlat=36.5, minlon=-125, maxlon=-115, minDt=datetime.datetime(1984,01,01), maxDt=datetime.datetime(1985,12,31)):
	rpb=intervalRecordBreaker()
	rpb.setNormalCat('cats/socalRB19841985.cat')
	rpb.plotRBintervalSet(2,4.5,.25)

def doHminePlots(incat='cats/hminefull.cat', minmag=2.5, maxmag=5.0):
	# ultimately, this will be a singel function call to produce the pertinent hector mine (data?) and plots. (right now, we start
	# from a catalog. we could produce the catalog from sql here as well for posterity).
	# 
	# catFname='parkcat.cat', theta=tTheta, clat=tLat, clon=tLon, ra=tA, rb=tB)
	# makeHmineCats(dtMinus=365.24*8, dtPlus=365.24*5, rx=.5, ry=.5, ra=.65, rb=.2, theta=67.4, fullcatout='cats/hminefull.cat', shockcat='cats/hmineshock.cat'):
	# t0=datetime.datetime(1999, 10, 16, 02, 46, 44)
	
	# first, make a catalog from sql. start, for now, immediately after HM.
	catname='cats/hminesquareshocks.cat'
	shcatname='cats/hmineelshocks.cat'	#"shock" cat. we'll make a new one inside the intervalRecordBreaker() object.
	#
	#makeHmineCats(dtMinus=0, dtPlus=365.24*5, rx=.75, ry=.75, ra=.65, rb=.2, theta=67.4, fullcatout='cats/hminefull.cat', shockcat='cats/hmineshock.cat')
	#makeHmineCats(0, 365.24*5, .75, .75, .65, .2, 67.4, catname, shcatname)
	
	rbi=intervalRecordBreaker(incat, 67.4, 34.594, -116.271, .65, .2)
	# reset the catalog. note we're getting only events after the mainshock.
	rbi.setAftershockCatalog(incat, 67.4, 34.594, -116.271, .65, .2, datetime.datetime(1999, 10, 16, 02, 46, 45), datetime.datetime(2007,10,17) )
		
	rbi.plotRBintervalSet(minmag, maxmag, .25, 'images/hmineshock/')
	
	return rbi

#def doHminePlots2(incat='cats/hminefull.cat', minmag=2.75, maxmag=2.75, dmag=.25):
def doHminePlots2(incat='cats/hminefull.cat', minmag=2.75, outdir='images/hmineshock'):
	# create second set of plot-types. here, we'll look at one or two mag-bins at different starting times, aka, 0, 10, 100, 200, x[0]
	# events after the mainshock.
	# note, we'll create a new "doAllMyPlots" type function.
	#
	# ultimately, this will be a singel function call to produce the pertinent hector mine (data?) and plots. (right now, we start
	# from a catalog. we could produce the catalog from sql here as well for posterity).
	# 
	# catFname='parkcat.cat', theta=tTheta, clat=tLat, clon=tLon, ra=tA, rb=tB)
	# makeHmineCats(dtMinus=365.24*8, dtPlus=365.24*5, rx=.5, ry=.5, ra=.65, rb=.2, theta=67.4, fullcatout='cats/hminefull.cat', shockcat='cats/hmineshock.cat'):
	# t0=datetime.datetime(1999, 10, 16, 02, 46, 44)
	
	# first, make a catalog from sql. start, for now, immediately after HM.
	catname='cats/hminesquareshocks.cat'
	shcatname='cats/hmineelshocks.cat'	#"shock" cat. we'll make a new one inside the intervalRecordBreaker() object.
	#outdir='images/hmineshock/'
	#
	rbPlots=[]	# a list to hold our plot data for fitting. this will be like: [ [[x1],[y1]], [[x2],[y2]], ...].
					# a data set will be rbPlots[i]; replotting a series goes like: plot(rbPlots[i][0], rpPlots[i][1])
	#
	#makeHmineCats(dtMinus=0, dtPlus=365.24*5, rx=.75, ry=.75, ra=.65, rb=.2, theta=67.4, fullcatout='cats/hminefull.cat', shockcat='cats/hmineshock.cat')
	#makeHmineCats(0, 365.24*5, .75, .75, .65, .2, 67.4, catname, shcatname)
	#
	# shock-limited hmine catalog...
	rbi=intervalRecordBreaker(incat, 67.4, 34.594, -116.271, .65, .2)
	# reset the catalog. note we're getting only events after the mainshock.
	# setAftershockCatalog(self, catFname=None, theta=tTheta, clat=35.9, clon=-120.5, ra=tA, rb=tB, eventDate=datetime.datetime(2004,9,28, 17,15,24), maxDate=datetime.datetime(2009,9,28, 17,15,24), skipSeconds=0)
	rbi.setAftershockCatalog(incat, 67.4, 34.594, -116.271, .65, .2, datetime.datetime(1999, 10, 16, 02, 46, 45), datetime.datetime(2009,10,17), 0)
	# this is our base catalog and RB object. now, run a RB series. when do we get our first RB event?
	recs1=rbi.getRBintervals(minmag) # returns [biggers, smallers] -> [ [ [dates],[interval],[nRecords],[nth earthquake (natural time)] ], [ [],[],[],[] ] ]
	dNmax=recs1[0][3][1]	#NT of the first (non-trivial) RB interval.
	
	#dNlist=[0, 4, 16, 64, 256, 1024]
	secdays=86400.0
	currMag=4.0
	intCurrMag=long(10*currMag)
	strmag=str(currMag).replace('.', '')
	
	#maglags = [[2.5, [2*secdays, 3*secdays], '.-'], [3.0, [2*secdays, secdays], '+-'], [3.5, [2*secdays, 1*secdays, .1*secdays], '-^'], [4.0, [2*secdays, 1*secdays,.1*secdays], '-h'] ]
	#maglags = [[currMag, [20*secdays, 10*secdays, 5*secdays, 3*secdays, 2*secdays, 1*secdays,.1*secdays], '-h'] ]
	maglags = [[currMag, [20*secdays, 10*secdays, 5*secdays, 3*secdays, 2*secdays, 1*secdays,.5*secdays,.1*secdays], '-h'] ]
	# clear figures:
	for i in xrange(4):
		plt.figure(i)
		plt.clf()
		plt.cla()
		
	minx=1.0
	for rw in maglags:
		thismag=rw[0]
		thismarker=rw[2]
		for dT in rw[1]:
			# it's messy but simple; make a whole bunch of data-files.
			dfileInts=open('%s/hmineInts%s-dT%s.dat' % (outdir, strmag, int(dT*10/secdays)), 'w')
			dfileIntsNT=open('%s/hmineIntsNT%s-dT%s.dat' % (outdir, strmag, int(dT*10/secdays)), 'w')
			dfileNRB=open('%s/hmineNRB%s-dT%s.dat' % (outdir, strmag, int(dT*10/secdays)), 'w')
			dfileNRBnt=open('%s/hmineNRBnt%s-dT%s.dat' % (outdir, strmag, int(dT*10/secdays)), 'w')
			#
			dfileInts.write("#hmine event RB data, RB Interval Duration (time-time)\n#t\tinterval\n")
			dfileIntsNT.write("#hmine event RB data, RB Interval Duration (NT)\n#n\tinterval\n")
			dfileNRB.write("#hmine event RB data, NRB (time-time)\n#t\tnrb\n")
			dfileNRBnt.write("#hmine event RB data, NRB (NT)\n#n\tnrb\n")
			#
			thismarker=rw[2]
			if dT==10*secdays: thismarker='-+'
			rbi.setAftershockCatalog(incat, 67.4, 34.594, -116.271, .65, .2, datetime.datetime(1999, 10, 16, 02, 46, 45), datetime.datetime(2009,10,17), dT)
			recs2=rbi.getRBintervals(thismag)
			#
			# NT-NRB:
			plt.figure(0)
			x=[]
			y=[]
			for i in xrange(len(recs2[0][3])):
				if log10(recs2[0][3][i])<minx: continue
				x+=[recs2[0][3][i]]
				y+=[recs2[0][2][i]]
				dfileNRBnt.write("%d\t%d\n" % (recs2[0][3][i], recs2[0][2][i]))
			#
			rbPlots+=[[recs2[0][3], recs2[0][2], [thismag, dT, thismarker]]]
			plt.loglog(x, y, thismarker, label='m%s,dT%s' % (str(thismag)[0:4], float(dT)/secdays))
			#
			# NT-IntervalDuration:
			plt.figure(1)
			x=[]
			y=[]
			for i in xrange(len(recs2[0][3])):
				if log10(recs2[0][3][i])<minx: continue
				x+=[recs2[0][3][i]]
				y+=[recs2[0][1][i]]
				dfileIntsNT.write("%d\t%f\n" % (recs2[0][3][i], recs2[0][1][i]))
			plt.loglog(x, y, thismarker, label='m%s,dT%s' % (str(thismag)[0:4], float(dT)/secdays))
			#
			# NRB time-time:
			plt.figure(2)
			x=[]
			y=[]
			firstdt=datetimeToFloat(recs2[0][0][0])
			print "firstdt: %f" % firstdt
			for i in xrange(len(recs2[0][3])):
				#if log10(recs2[0][3][i])<minx: continue
				x+=[datetimeToFloat(recs2[0][0][i])-firstdt]
				y+=[recs2[0][2][i]]
				dfileNRB.write("%f\t%d\n" % (datetimeToFloat(recs2[0][0][i])-firstdt, recs2[0][2][i]))
			plt.loglog(x, y, thismarker, label='m%s,dT%s' % (str(thismag)[0:4], float(dT)/secdays))
			#
			# interval duration, time-time
			plt.figure(3)
			x=[]
			y=[]
			firstdt=datetimeToFloat(recs2[0][0][0])
			for i in xrange(len(recs2[0][3])):
				#if log10(recs2[0][3][i])<minx: continue
				x+=[datetimeToFloat(recs2[0][0][i])-firstdt]
				y+=[recs2[0][1][i]]
				dfileInts.write("%f\t%d\n" % (datetimeToFloat(recs2[0][0][i])-firstdt, recs2[0][1][i]))
			plt.loglog(x, y, thismarker, label='m%s,dT%s' % (str(thismag)[0:4], float(dT)/secdays))
			#
			dfileInts.close()
			dfileIntsNT.close()
			dfileNRB.close()
			dfileNRBnt.close()
			#
			
			print "SeriesLen(%d): %d" % (dT, len(recs2[0][2]))
	
	strFmt='eps'
	plt.figure(0)
	plt.title("Number of RB Events (NT)\n(Hector Mine)")
	#plt.xlabel("n events since mainshock")
	plt.xlabel("n events since RB start")
	plt.ylabel("Number of Record Breaking Events")
	plt.legend(loc='upper left')
	#plt.legend(loc='lower right')
	plt.savefig("%s/nrbNTHminedTdM-%d-%d.pdf" % (outdir, dT, intCurrMag), format='%s' % strFmt)
	#
	#
	plt.figure(1)
	plt.title("RB Interval Duration (NT)\n(Hector Mine)")
	#plt.xlabel("n events since mainshock")
	plt.xlabel("n events since RB start")
	plt.ylabel("Interval Duration (days)")
	plt.legend(loc='upper left')
	#plt.legend(loc='lower right')
	plt.savefig("%s/RBintervalDurationsNTHminedTdM-%d-%d.pdf" % (outdir, dT, intCurrMag), format='%s' % strFmt)
	#
	plt.figure(2)
	plt.title("Number of RB Events (time)\n(Hector Mine)")
	#plt.xlabel("n events since mainshock")
	plt.xlabel("days since main event")
	plt.ylabel("Number of Record Breaking Events")
	plt.legend(loc='upper left')
	#plt.legend(loc='lower right')
	plt.savefig("%s/nrbHminedTdM-%d-%d.%s" % (outdir, dT, intCurrMag, strFmt), format='%s' % strFmt)
	#
	plt.figure(3)
	plt.title("RB Interval Duration (time)\n(Hector Mine)")
	#plt.xlabel("n events since mainshock")
	plt.xlabel("days since main event")
	plt.ylabel("Interval Duration (days)")
	plt.legend(loc='upper left')
	#plt.legend(loc='lower right')
	plt.savefig("%s/RBintervalDurationsHminedTdM-%d-%d.%s" % (outdir, dT, intCurrMag, strFmt), format='%s' % strFmt)
	#
	
	# for now, let's take a break on data-fitting:
	'''
	# and let's do some data fitting:
	# for now, fit t>=10**-1.5 in the Nrecords plot, maybe interval>10**-2 for RBinterval??
	# alternatively, always trim the first (maybe first 2) data points.
	#
	# we want the log/log fit, so we can either construct a log-based error or we can take the log-log
	# of the data and do a linear fit.
	#
	fitPlotsInt=[]	# fit RB interval (mags)
	fitPlotsN=[]	# fit to Nrb plots.
	slopes=[]
	plt.figure(1)
	plt.clf()
	plt.cla()
	#plt.plot(NevLog, YevLog, '.')
	for pset in rbPlots:
		x=[]
		y=[]
		yfit=[]
		#print len(pset[0])
		for i in xrange(len(pset[0])):
			# let's hard-code some fit-range conditions. maybe we can be more dynamic later...
			#if log10(pset[1][i])<.6: continue
			if log10(pset[0][i])<minx: continue
			x+=[log10(pset[0][i])]
			y+=[log10(pset[1][i])]
		if len(x)<=1: continue	
		#
		# now, we have a log-log set of one dataset. fit it to a line...
		#print "x: %s" % str(x)
		#print "y: %s" % str(y)
		p=scipy.array([0,1])
		#print "x: %s" % str(x)
		
		plsq=spo.leastsq(linRes, p, args=(scipy.array(y), scipy.array(x)), full_output=0, maxfev=20000)	# (function to minimize, initial prams, argument-arrays, max-iterations)
		slopes+=[plsq[0][1]]
		# make an array of the fit function. note we can plot this or get a chi-sqr bc we have data-points.
		for X in x:
			fitval=plsq[0][0] + X*plsq[0][1]
			yfit+=[fitval]
		plt.plot(x,yfit, pset[2][2], label='M%s,dT=%s,m=%s' % (pset[2][0], float(pset[2][1])/secdays, str(slopes[-1])[0:4] ))
	plt.title("Fit to Nrb (NT)\n(Hector Mine)")
	plt.xlabel("log10(nevents since mainshock)")
	plt.ylabel("log10(Nrb)")
	#plt.legend(loc='lower right')
	plt.legend(loc='upper left')
	plt.savefig("%sNrbFitsNTHmineNTdTdM-%d-%d.pdf" % (outdir, dT, intCurrMag), format='pdf')
	#plt.legend(loc='upper left')
	#
	#plt.figure(2)
	#plt.title("Hector Mine Interval Slopes (NT)")
	#plt.xlabel("mag")
	#plt.ylabel("slope")
	#plt.plot(mags, slopes, '.-')
	#plt.savefig("%sRBintervalSlopesNT.pdf" % outdir)
	print "slopes: %s" % slopes
	#
	#return recs1
	#
	#rbi.plotRBintervalSet(minmag, maxmag, dmag, 'images/hmineshock/')
	#
	'''
	
def doNHPPrecords():
	# basically, copy doHminePlots2() but for a NHPP catalog.
	# getNHPPintervalsOmori1(self, catLen, initialRate=1, startDt=0, t0=0):
	# secdays=86400.0
	import recordBreaker as rb2	# maybe it's time to consolidate these two record breaking codes???
	
	rbp=rb2.recordbreaker()	# get any recordBreaker() object. this is the (typically) poisson RB tool from which we get our NHPP catalog.
	rbi = intervalRecordBreaker()	# this is the "one event" record breaker object (this module, original coded to study parkfield)
	# see: getNHPPintervalsOmori1(self, catLen, initialRate=1, startDt=0, t0=0) --> f_omori = initialRate/(t0+t). note, p=1. then, the startDate pram just starts the process
	# at a later time. see getNHPPintervalsOmori1(); we integrate L(t) in the numerator, so we get T (aka, startDt), current total time, in the random number generator.
	NHPPints = rbp.getNHPPintervalsOmori1(1000000, 100000, 0, 1)
	#NHPPints = rbp.getNHPPintervalsOmori1(10000, 1000, 0, 1)	# the initial rate here is like a full m0 interval with data spanning 3 magnitudes.
	#																			# aka, a M6 with m_min=3. here, startDt is a way of skipping to the middle of the sequence;
	#																			# t0 is the constant in the omori denominator: N=A/(t0 + t)**p. using t0=1 gives a simple
	#																			# GR distribution initial rate (aka, N0=N_GR (M_mainshock))
	#																			# NHHPints -> [[event number], [ft], [interval]]
	#
	# record the catalog:
	fout=open('images/NHPPcat/NHPPseries.dat', 'w')
	fout.write("#n\tt\tdt\n")
	for i in xrange(len(NHPPints[0])):
		fout.write("%d\t%f\t%f\n" % (NHPPints[0][i], NHPPints[1][i], NHPPints[2][i]))
	fout.close()
		
	#getRecordArrays(self, dset=None, NT=True, winLen=256, winStep=1)
	#newrecs=[[],[],[],[]]
	newrecs=[]
	for i in xrange(len(NHPPints[0])):
		#newrecs[0]+=[floatToDateTime(1+NHPPints[1][i])]
		#newrecs[1]+=[42]
		#newrecs[2]+=[42]
		#newrecs[3]+=[5.5]	# these 3 columns have arbitrary values. we just have to fill the list.
		#
		newrecs+=[[floatToDateTime(1+NHPPints[1][i]), 42, 42, 5.5]]
		#
	#
	rbi.fullCat=newrecs[:]
	rbi.shockCat=newrecs[:]
	rbi.eventDtTime=rbi.fullCat[0][0]
	# self.eventDtTime
		
	print "len(shockCat): %d, %s" % (len(rbi.shockCat[3]), rbi.shockCat[3][0])
	rbi.plotRBintervalSet(2.5, 2.5, .1, 'images/NHPPcat/')
	
	return rbi



def doNHPPrecords3(Nits=100, Tmax=365, tao=.01, cm=1.0, binsize=1.0, binsizeNT=1.0, maxLen=1000, doshow=True):
	maxLenSmall=0
	# this one is a little different. we have to average over many realizations. getRBintervals() returns a set of big and small RB data.
	# so,
	# - get an RB set
	# - bin.
	# - add to cumulative bin-array.
	# - repeat.
	#
	#  getRBintervals(self, minmag=1.0)
	# returns [biggers, smallers]
	# returns [biggers, smallers] -> [ [ [dates],[interval],[nRecords],[nth earthquake (natural time)] ], [ [],[],[],[] ] ]
	#
	#binsize=1.0
	#binsizeNT=1.0
	#
	import recordBreaker as rb2
	reload(rb2)
	
	maxXbin=float(Tmax)
	if os.system('ls images/NHPPcat3')!=0: os.system('mkdir images/NHPPcat3')
	
	rbp = rb2.recordbreaker()	# get any recordBreaker() object. this is the (typically) poisson RB tool from which we get our NHPP catalog.
	rbi = intervalRecordBreaker()	# this is the "one event" record breaker object (this module, original coded to study parkfield)
	#
	bigRecInts=[[binsize],[0],[0]]		# [[bin], [<val>], [<val**2>]
	bigRecIntsNT=[[binsizeNT],[0],[0]]
	bigNRB=[[binsize],[0],[0]]
	bigNRBnt=[[binsizeNT],[0],[0]]
	#
	# and the smallies:
	smallRecInts=[[binsize],[0],[0]]		# [[bin], [<val>], [<val**2>]
	smallRecIntsNT=[[binsizeNT],[0],[0]]
	smallNRB=[[binsize],[0],[0]]
	smallNRBnt=[[binsizeNT],[0],[0]]
	
	for x in xrange(4):
		plt.figure(x)
		plt.clf()
		plt.cla()
	#
	# write data-file headers:
	fs="images/NHPPcat3/scatter.dat"
	fb1="images/NHPPcat3/binnedint.dat"
	fb2="images/NHPPcat3/binnedintNT.dat"
	fb3="images/NHPPcat3/binnedNRB.dat"
	fb4="images/NHPPcat3/binnedNRBnt.dat"
	fb5="images/NHPPcat3/binnedNRBnt-log2.dat"
	fb6="images/NHPPcat3/binnedintNT-log2.dat"
	fscat=open(fs, "w")
	fbin1=open(fb1, "w")
	fbin2=open(fb2, "w")
	fbin3=open(fb3, "w")
	fbin4=open(fb4, "w")
	fbin5=open(fb5, 'w')
	fbin6=open(fb6, 'w')
	#
	fscat.write("#scatter-data from NHPP3\n#Nits=%d;Tmax=%f;tao=%f;cm=%f;binsize=%f;binsizeNT=%f\n#n\tfdt\tinterval\tnrb\n" % (Nits, Tmax, tao, cm, binsize, binsizeNT))
	fbin1.write("#binned-data from NHPP3, intervals\n#Nits=%d;Tmax=%f;tao=%f;cm=%f;binsize=%f;binsizeNT=%f\n#fdt\tinterval\n" % (Nits, Tmax, tao, cm, binsize, binsizeNT))
	fbin1.close()
	fbin2.write("#binned-data from NHPP3, intervals NT\n#Nits=%d;Tmax=%f;tao=%f;cm=%f;binsize=%f;binsizeNT=%f\n#n\tinterval\n" % (Nits, Tmax, tao, cm, binsize, binsizeNT))
	fbin2.close()
	fbin3.write("#binned-data from NHPP3\n#Nits=%d;Tmax=%f;tao=%f;cm=%f;binsize=%f;binsizeNT=%f\n#fdt\tnrb\n" % (Nits, Tmax, tao, cm, binsize, binsizeNT))
	fbin3.close()
	fbin4.write("#binned-data from NHPP3\n#Nits=%d;Tmax=%f;tao=%f;cm=%f;binsize=%f;binsizeNT=%f\n#n\tnrb\n"% (Nits, Tmax, tao, cm, binsize, binsizeNT))
	fbin4.close()
	fbin5.write("#binned-data from NHPP3\n#Nits=%d;Tmax=%f;tao=%f;cm=%f;binsize=%f;binsizeNT=%f\n#fdt\tnrb\n" % (Nits, Tmax, tao, cm, binsize, binsizeNT))
	fbin5.close()
	fbin6.write("#binned-data from NHPP3\n#Nits=%d;Tmax=%f;tao=%f;cm=%f;binsize=%f;binsizeNT=%f\n#n\tnrb\n"% (Nits, Tmax, tao, cm, binsize, binsizeNT))
	fbin6.close()
	############
	#
	for i in xrange(Nits):
		intervals=rbp.getNHPPintervalsOmori1b(Tmax, tao, cm)	#-> [[n], [T], [interval]]
		#if len(intervals[0])>maxLen: maxLen=len(intervals[0])
		#
		# put intervals "catalog" into rbi object:
		newrecs=[]
		for ii in xrange(len(intervals[0])):
			newrecs+=[[floatToDateTime(1+intervals[1][ii]), 42, 42, 5.5]]
		#
		rbi.fullCat=newrecs[:]
		rbi.shockCat=newrecs[:]
		rbi.eventDtTime=rbi.fullCat[0][0]
		#
		# get record-breaking intervals (returns [[biggers], [smallers]], aka, big RBevents, small RBevents)
		if len(rbi.shockCat)<=1: continue		# sometimes, when we have a high rate, we get a single event (or two) that jumps out of the interval, and we get an error.
		allrecs=rbi.getRBintervals()
		if len(allrecs[0][0])>maxLen: maxLen=len(allrecs[0][0])
		if len(allrecs[1][0])>maxLenSmall: maxLenSmall=len(allrecs[1][0])
		# "lots of points" plots:
		if Nits<100 or i%(Nits/100)==0:
			# note on datetimes->float: times start at 1, so our running time since event is shifted by 1. next time, don't use datetimes; use
			# float-time, convert to datetimes for pretty plots, etc. as necessary. for now, subtract 1 from running time.
			for ifile in xrange(len(allrecs[0][3])): fscat.write("%d\t%f\t%f\t%d\n" % (allrecs[0][3][ifile], datetimeToFloat(allrecs[0][0][ifile])-1, allrecs[0][1][ifile], allrecs[0][2][ifile]))
			print "NHPP3 Nit: %d" % i
			plt.figure(0)
			currX=[]
			for dt in allrecs[0][0]: currX+=[datetimeToFloat(dt)-1]
			#
			plt.loglog(currX, allrecs[0][1], '.')
			#plt.title("NHPP RB Interval Durations")
			#
			plt.figure(1)
			#for dt in allrecs[0][0]: currX+=[datetimeToFloat(dt)]
			plt.semilogy(allrecs[0][3], allrecs[0][1], '.')
			#plt.title("NHPP RB Interval Durations, NT")
			#
			plt.figure(2)
			#for dt in allrecs[0][0]: currX+=[datetimeToFloat(dt)]
			plt.loglog(currX, allrecs[0][2], '.')
			#plt.title("NHPP NRB")
			#
			plt.figure(3)
			#for dt in allrecs[0][0]: currX+=[datetimeToFloat(dt)]
			plt.semilogy(allrecs[0][3], allrecs[0][2], '.')
			#plt.title("NHPP NRB, NT")
			#
		
		# bin cumulative NRB. for interval magnitudes, use maximum value up to the end of that bin. then, add bins and average.
		# note: recbinVals() converts the datetime column (allrecs[0][0] to a float.
		thisBigRecInts=recbinVals(allrecs[0][0][:], allrecs[0][1][:], binsize, maxXbin)		# rbIntervals (note: allrecs[0][0] is a column of datetime objects.)
		thisBigNRB=recbinVals(allrecs[0][0][:], allrecs[0][2][:], binsize, maxXbin)			# NRB
		thisBigRecIntsNT=recbinVals(allrecs[0][3][:], allrecs[0][1][:], binsizeNT, maxLen)	# RBintervals-NT
		thisBigNRBnt=recbinVals(allrecs[0][3][:], allrecs[0][2][:], binsizeNT, maxLen)		#NRBnt
		#
		thisSmallRecInts=recbinVals(allrecs[1][0][:], allrecs[1][1][:], binsize, maxXbin)		# rbIntervals (note: allrecs[0][0] is a column of datetime objects.)
		thisSmallNRB=recbinVals(allrecs[1][0][:], allrecs[1][2][:], binsize, maxXbin)			# NRB
		thisSmallRecIntsNT=recbinVals(allrecs[1][3][:], allrecs[1][1][:], binsizeNT, maxLenSmall)	# RBintervals-NT
		thisSmallNRBnt=recbinVals(allrecs[1][3][:], allrecs[1][2][:], binsizeNT, maxLenSmall)		#NRBnt
		#
		# now, add "this" histogram to the total histogram. all X arrays will be continuous (by binsize) up to their maximum bin-size. add Y values; we'll normalize later.
		# normalization factor:
		normfact=float(Nits)
		normfactNT=float(Nits)
		ii=0
		while ii<len(thisBigRecInts[0]):
		#for ii in xrange(len(thisBigRecInts[0])):
			if ii>(len(bigRecInts[0])-1):
				bigRecInts[0]+=[bigRecInts[0][-1]+binsize]
				bigRecInts[1]+=[0]
				bigRecInts[2]+=[0]
				#bigRecInts[1]+=[bigRecInts[1][-1]]
			bigRecInts[1][ii]+=thisBigRecInts[1][ii]/normfact
			bigRecInts[2][ii]+=(thisBigRecInts[1][ii]**2)/(normfact)
			ii+=1
		while ii<len(bigRecInts[0]):
			bigRecInts[1][ii]+=thisBigRecInts[1][-1]/normfact
			bigRecInts[2][ii]+=(thisBigRecInts[1][-1]**2)/(normfact)
			ii+=1
			#
		#
		#
		# adjust arrays to be the same size:
		ii=0
		while ii<len(thisBigRecIntsNT[0]):
		#for ii in xrange(len(thisBigRecIntsNT[0])):
			if ii>(len(bigRecIntsNT[0])-1):
				bigRecIntsNT[0]+=[bigRecIntsNT[0][-1]+binsizeNT]
				bigRecIntsNT[1]+=[0]
				bigRecIntsNT[2]+=[0]
				#bigRecIntsNT[1]+=[bigRecIntsNT[1][-1]]
			bigRecIntsNT[1][ii]+=thisBigRecIntsNT[1][ii]/normfactNT
			bigRecIntsNT[2][ii]+=(thisBigRecIntsNT[1][ii]**2)/(normfactNT)
			ii+=1
		while ii<len(bigRecIntsNT[0]):
			bigRecIntsNT[1][ii]+=thisBigRecIntsNT[1][-1]/normfactNT
			bigRecIntsNT[2][ii]+=(thisBigRecIntsNT[1][-1]**2)/(normfactNT)
			ii+=1
			#
		#
		#
		ii=0
		while ii<len(thisBigNRB[0]):
		#for ii in xrange(len(thisBigNRB[0])):
			if ii>(len(bigNRB[0])-1):
				bigNRB[0]+=[bigNRB[0][-1]+binsize]
				bigNRB[1]+=[0]
				bigNRB[2]+=[0]
				#bigNRB[1]+=[bigNRB[1][-1]]
			bigNRB[1][ii]+=thisBigNRB[1][ii]/float(Nits)
			bigNRB[2][ii]+=(thisBigNRB[1][ii]**2)/(float(Nits))
			ii+=1
		while ii<len(bigNRB[0]):
			bigNRB[1][ii]+=thisBigNRB[1][-1]/float(Nits)
			bigNRB[2][ii]+=(thisBigNRB[1][-1]**2)/(float(Nits))
			ii+=1
			#
		#	
		#
		ii=0
		while ii<len(thisBigNRBnt[0]):
		#for ii in xrange(len(thisBigNRBnt[0])):
			if ii>(len(bigNRBnt[0])-1):
				bigNRBnt[0]+=[bigNRBnt[0][-1]+binsizeNT]
				bigNRBnt[1]+=[0]
				bigNRBnt[2]+=[0]
				#bigNRBnt[1]+=[bigNRBnt[1][-1]]
			bigNRBnt[1][ii]+=thisBigNRBnt[1][ii]/float(Nits)
			bigNRBnt[2][ii]+=(thisBigNRBnt[1][ii]**2)/(float(Nits))
			ii+=1
		while ii<len(bigNRBnt[0]):
			bigNRBnt[1][ii]+=thisBigNRBnt[1][-1]/float(Nits)
			bigNRBnt[2][ii]+=(thisBigNRBnt[1][-1]**2)/(float(Nits))
			ii+=1
			#
		#
		#
	fbin1=open(fb1, "a")
	for ii in xrange(len(bigRecInts[0])): fbin1.write("%f\t%f\t%f\n" % (bigRecInts[0][ii], bigRecInts[1][ii], bigRecInts[2][ii]-bigRecInts[1][ii]**2))
	fbin1.close()
	#
	fbin2=open(fb2, "a")
	fbin6=open(fb5, 'a')
	for ii in xrange(len(bigRecIntsNT[0])):
		fbin2.write("%d\t%f\t%f\n" % (bigRecIntsNT[0][ii], bigRecIntsNT[1][ii], bigRecIntsNT[2][ii]-bigRecIntsNT[1][ii]**2))
		if math.log(bigRecIntsNT[0][ii],2)%1==0:
			fbin6.write("%d\t%f\t%f\n" % (bigRecIntsNT[0][ii], bigRecIntsNT[1][ii], bigRecIntsNT[2][ii]-bigRecIntsNT[1][ii]**2))
	fbin6.close()
	fbin2.close()
	#
	fbin3=open(fb3, "a")
	for ii in xrange(len(bigNRB[0])): fbin3.write("%f\t%f\t%f\n" % (bigNRB[0][ii], bigNRB[1][ii], bigNRB[2][ii]-bigNRB[1][ii]**2))
	fbin3.close()
	#
	fbin4=open(fb4, "a")
	fbin5=open(fb6, "a")
	for ii in xrange(len(bigNRBnt[0])):
		fbin4.write("%d\t%f\t%f\n" % (bigNRBnt[0][ii], bigNRBnt[1][ii], bigNRBnt[2][ii]-bigNRBnt[1][ii]**2))
		if math.log(bigNRBnt[0][ii],2)%1==0:
			fbin5.write("%d\t%f\t%f\n" % (bigNRBnt[0][ii], bigNRBnt[1][ii], bigNRBnt[2][ii]-bigNRBnt[1][ii]**2))
	fbin5.close()
	fbin6.close()
	
	#
	# binned, cumulative type plots:
	# , ls='-', marker='^', lw=2
	plt.figure(0)
	plt.loglog(bigRecInts[0], bigRecInts[1], '-')
	plt.title("NHPP Large Record Intervals")
	plt.xlabel("days")
	plt.ylabel("Mean RB Interval Duration (days)")
	plt.savefig("images/NHPPcat3/NHPPints.svg")
	
	plt.figure(1)
	plt.semilogy(bigRecIntsNT[0], bigRecIntsNT[1], '-')
	plt.title("NHPP Large Record Intervals (NT)")
	plt.xlabel("nEvents")
	plt.ylabel("Mean RB Interval Duration (days)")
	plt.savefig("images/NHPPcat3/NHPPintsNT.svg")
	
	plt.figure(2)
	plt.loglog(bigNRB[0], bigNRB[1], '-')
	plt.title("NHPP Large Records: NRB")
	plt.xlabel("days")
	plt.ylabel("Mean NRB (big)")
	plt.savefig("images/NHPPcat3/NHPP-NRB.svg")
	
	plt.figure(3)
	plt.semilogy(bigNRBnt[0], bigNRBnt[1], '-')
	plt.title("NHPP Large Records: NRB (NT)")
	plt.xlabel("nEvents")
	plt.ylabel("Mean NRB (big)")
	plt.savefig("images/NHPPcat3/NHPP-NRBnt.svg")
	
	fscat.close()
	if doshow==True: plt.show()
		
def recbinVals(X, Y, binsize=1.0, maxBin=None):
	#maxBin=None
	# note: this assumes the list is in order...
	# X: x-vals
	# Y: y-vals
	# binsize: bin size.
	# for the time being, return floats for X.
	# note: this process will satisfy magnitude and NRB type 'binning'.
	# returns an array like: [[X], [Y]]. where X are the upper, non-inclusive limit. aka, max Y up-to,not including X[i]
	#
	# 6 july 2009 yoder:
	# add maxBins thing. if we want our series to always extend to some max time/N/whatever.
	# naturally, series come in with varying lengths. pseudo-bins at the end of the observation period
	# can be unoccupied, so the bin-value, averaged over many realizations, declines as an artifact of 
	# the observation window and binning method.
	#
	#print "xy lens0: %d, %d" % (len(X), len(Y))
	if type(X[0]).__name__ in ['datetime', 'date']:
		# convert X to floats:
		for i in xrange(len(X)):
			X[i]=datetimeToFloat(X[i])-1		# note: typically, we want time since an event. datetimeToFloat() starts at minDate=1.0 (as per dtime limitations).
														# never again use datetime objecs to hold low level variables. float->datetime is too easy...
	#print "xy lens0: %d, %d" % (len(X), len(Y))
	#
	# now, is the dataset long enough? if we have a maxBins requirement, extend the DS as necessary:
	if maxBin!=None:
		while X[-1]<=maxBin+binsize:
			#if X[-1]+binsize>maxBin: break
			X+=[X[-1]+binsize]
			Y+=[0]
	#print "xy lens: %d, %d" % (len(X), len(Y))	
	####
	# cum-bin the data:
	retAry=[[binsize],[0]]
	for i in xrange(len(X)):
		thisbin=int(X[i]/binsize)
		# what about exactly (inclusive)? NT binning gets screwed up; (X[i]=1)/1 -> thisbin=1, but clearly we want thisbin=0. this is not so problematic in time-time
		# since time-time is continuous. let's make our bins inclusive of their upper value:
		if X[i]%binsize==0: thisbin-=1	# so, for NT, n=1 -> bin_0; n=2->bin_1. for time-time, binsize=.1, .1-> bin_0, .2->bin_1, etc.
		#
		while len(retAry[0])<=(thisbin+1):
			retAry[0]+=[retAry[0][-1]+binsize]
			retAry[1]+=[retAry[1][-1]]
			#
		#print "i, thisbin: %f, %f, %f, %d, %d, %f, %f" % (i, X[i], thisbin, len(retAry[1]), len(Y), X[-1], Y[-1])
		if Y[i]>retAry[1][thisbin]: retAry[1][thisbin]=Y[i]
	#
	# now, spin through the list. any empty sites (0) are equal to the previous bin.
	for i in xrange(1, len(retAry[0])):
		if retAry[1][i]==0: retAry[1][i]=retAry[1][i-1]	# arguably, a better condition is: if retAry[1][i]<retAry[1][i-1]: retAry[1][i]=retAry[1][i-1]. BUT, this should never happen...
	#
	ftemp=open("temp.dat", 'w')
	for i in xrange(len(retAry[0])):
		ftemp.write("%f\t%f\n" % (retAry[0][i], retAry[1][i]))
	ftemp.close()
	
	return retAry

def reclogbinVals(X, Y, logbase=10, maxBin=None):
	#maxBin=None
	# note: this assumes the list is in order...
	# X: x-vals
	# Y: y-vals
	# binsize: bin size.
	# for the time being, return floats for X.
	# note: this process will satisfy magnitude and NRB type 'binning'.
	# returns an array like: [[X], [Y]]. where X are the upper, non-inclusive limit. aka, max Y up-to,not including X[i]
	#
	# 6 july 2009 yoder:
	# add maxBins thing. if we want our series to always extend to some max time/N/whatever.
	# naturally, series come in with varying lengths. pseudo-bins at the end of the observation period
	# can be unoccupied, so the bin-value, averaged over many realizations, declines as an artifact of 
	# the observation window and binning method.
	#
	#print "xy lens0: %d, %d" % (len(X), len(Y))
	if type(X[0]).__name__ in ['datetime', 'date']:
		# convert X to floats:
		for i in xrange(len(X)):
			X[i]=datetimeToFloat(X[i])-1		# note: typically, we want time since an event. datetimeToFloat() starts at minDate=1.0 (as per dtime limitations).
														# never again use datetime objecs to hold low level variables. float->datetime is too easy...
	#print "xy lens0: %d, %d" % (len(X), len(Y))
	#
	# now, is the dataset long enough? if we have a maxBins requirement, extend the DS as necessary:
	if maxBin!=None:
		#while X[-1]<=maxBin+binsize:
		#	#if X[-1]+binsize>maxBin: break
		#	X+=[X[-1]+binsize]
		#	Y+=[0]
		X+=[maxBin]
		Y+=[0]
	#print "xy lens: %d, %d" % (len(X), len(Y))	
	####
	# cum-bin the data:
	print "starters: %f" % (log10(27))
	minx=int(math.log(X[0],logbase))
	retAry=[[minx],[0]]
	for i in xrange(len(X)):
		#thisbin=int(X[i]/binsize)
		thisbin=int(math.log(X[i],logbase))-minx
		while len(retAry[0])<=(thisbin+1):
			retAry[0]+=[retAry[0][-1]*logbase]
			retAry[1]+=[retAry[1][-1]]
			#
		#print "i, thisbin: %f, %f, %f, %d, %d, %f, %f" % (i, X[i], thisbin, len(retAry[1]), len(Y), X[-1], Y[-1])
		if Y[i]>retAry[1][thisbin]: retAry[1][thisbin]=Y[i]
	#
	# now, spin through the list. any empty sites (0) are equal to the previous bin.
	for i in xrange(1, len(retAry[0])):
		if retAry[1][i]==0: retAry[1][i]=retAry[1][i-1]	# arguably, a better condition is: if retAry[1][i]<retAry[1][i-1]: retAry[1][i]=retAry[1][i-1]. BUT, this should never happen...
	#
#	ftemp=open("templog.dat", 'w')
#	for i in xrange(len(retAry[0])):
#		ftemp.write("%f\t%f\n" % (retAry[0][i], retAry[1][i]))
#	ftemp.close()
	
	return retAry
		
def doNHPPrecords2(Tmax=365, tao=1.0, cm=1.0, returnSeries=False):
	# secdays=86400.0
	# basically, copy doHminePlots2() but for a NHPP catalog.
	# getNHPPintervalsOmori1(self, catLen, initialRate=1, startDt=0, t0=0):
	# maybe it's time to consolidate these two record breaking codes???
	import recordBreaker as rb2
	reload(rb2)
	
	rbp=rb2.recordbreaker()	# get any recordBreaker() object. this is the (typically) poisson RB tool from which we get our NHPP catalog.
	rbi = intervalRecordBreaker()	# this is the "one event" record breaker object (this module, original coded to study parkfield)
	# see: getNHPPintervalsOmori1(self, catLen, initialRate=1, startDt=0, t0=0) --> f_omori = initialRate/(t0+t). note, p=1. then, the startDate pram just starts the process
	# at a later time. see getNHPPintervalsOmori1(); we integrate L(t) in the numerator, so we get T (aka, startDt), current total time, in the random number generator.
	#NHPPints = rbp.getNHPPintervalsOmori1(1000000, 100000, 0, 1)
	#NHPPints = rbp.getNHPPintervalsOmori1(10000, 1000, 0, 1)	# the initial rate here is like a full m0 interval with data spanning 3 magnitudes.
	#																			# aka, a M6 with m_min=3. here, startDt is a way of skipping to the middle of the sequence;
	#																			# t0 is the constant in the omori denominator: N=A/(t0 + t)**p. using t0=1 gives a simple
	#																			# GR distribution initial rate (aka, N0=N_GR (M_mainshock))
	#																			# NHHPints -> [[event number], [ft], [interval]]
	# getNHPPintervalsOmori1b(self, Tmax=10, tao=1, cm=1)	# time units must be consistent. let's use days...
	NHPPints = rbp.getNHPPintervalsOmori1b(Tmax, tao, cm)	# same as above (getNHPPintervalsOmori1(), except we limit T instead of Nevents.
	print "NHPPints len: %d" % len(NHPPints[0])
	#
	# record the catalog:
	fout=open('images/NHPPcat2/NHPPseries.dat', 'w')
	fout.write("#n\tt\tdt\n")
	for i in xrange(len(NHPPints[0])):
		fout.write("%d\t%f\t%f\n" % (NHPPints[0][i], NHPPints[1][i], NHPPints[2][i]))
	fout.close()
		
	#getRecordArrays(self, dset=None, NT=True, winLen=256, winStep=1)
	#newrecs=[[],[],[],[]]
	newrecs=[]
	for i in xrange(len(NHPPints[0])):
		#newrecs[0]+=[floatToDateTime(1+NHPPints[1][i])]
		#newrecs[1]+=[42]
		#newrecs[2]+=[42]
		#newrecs[3]+=[5.5]	# these 3 columns have arbitrary values. we just have to fill the list.
		#
		newrecs+=[[floatToDateTime(1+NHPPints[1][i]), 42, 42, 5.5]]
		#
	#
	rbi.fullCat=newrecs[:]
	rbi.shockCat=newrecs[:]
	rbi.eventDtTime=rbi.fullCat[0][0]
	# self.eventDtTime
		
	print "len(shockCat): %d, %s" % (len(rbi.shockCat[3]), rbi.shockCat[3][0])
	rbi.plotRBintervalSet(2.5, 2.5, .1, 'images/NHPPcat2/')
	if returnSeries:
		return NHPPints
	else:
		return rbi
		

def checkNHPP(Nits=10, Tmax=365, tao=1.0, cm=1.0, doShow=True, binres=1.0):
	# secdays=86400.0
	# 2009 june 29: this thing explodes in memory as Nits gets big. basically, we produce a lot of individual NHPP runs
	# and then try to pdf them all at once. we need to produce one pdf at a time.
	import recordBreaker as rb2
	reload(rb2)
	#doNHPPrecords2(Tmax=365, tao=1, cm=1, returnSeries=False)
	#
	# fetch a bunch of NHPP series and assemble a pdf.
	# start with Parkfield mc=1: cm=.1313, tao=38.066 sec.
	ints=[]
	rbp=rb2.recordbreaker()	# get any recordBreaker() object. this is the (typically) poisson RB tool from which we get our NHPP catalog.	
	# we'll call rbp.getNHPPintervalsOmori1b to return [[n], [T], [dt]]:
	
	pdf=[[],[]]
	thispdf=[[],[]]
	#
	# get the first series/pdf for base x-value set:
	thisInts=rbp.getNHPPintervalsOmori1b(Tmax, tao, cm)
	#print "get pdf: %d, %f" % (len(thisInts), binres)
	pdf=getpdf(thisInts[2], binres)
	
	for i in xrange(Nits-1):
		#print "Nit: %d" % (i+1)
		if i%(Nits/10)==0: print "doing checkNHPP(), Nit: %d of %d" % (i, Nits)
		#if i%1000==0: print "doing checkNHPP(), Nit: %d" % i
		#
		thisInts=rbp.getNHPPintervalsOmori1b(Tmax, tao, cm)
		#ints+=thisInts[2]
		thispdf=getpdf(thisInts[2], binres)
		
		#rint "lens: %d, %d, %f, %f, %f, %f" % (len(thispdf[0]), len(pdf[0]), thispdf[0][0], thispdf[0][-1], pdf[0][0], pdf[0][-1])
		#if len(thispdf[0])!=len(pdf[0]):
		#	print "length error."
		#	break
		#
		#while thispdf[0][-1]>pdf[0][-1]:
		while len(thispdf[0])>len(pdf[0]):
			#note: the above 2 statements should be functionally equivalent.
			#print "[[%f, %f]], %f" % (pdf[0][-1], thispdf[0][-1], pdf[0][-1]+binres)
			pdf[0]+=[pdf[0][-1]+binres]
			pdf[1]+=[0]
			
		#
		#print "%%%%%%%%"
		#print "pdf: %s :: %d" % (pdf[0], len(pdf[0]))
		#print "thispdf: %s :: %d" % (thispdf[0], len(thispdf[0]))
		#print "lens: %d, %d" % (len(pdf[0]), len(thispdf[0]))
		#
		for ii in xrange(len(thispdf[1])):
			#if ii>=len(pdf[0]):
			#	# our new pdf is longer than the previous pdf.
			#	print "(%d) extending pdf (%d, %d; %f, %f)" % (i, len(pdf[0]), ii, pdf[0][-1], thispdf[0][-1])
			#	pdf[0]+=[pdf[0][-1]+binres]
			#	pdf[1]+=[0]
			#
			pdf[1][ii]+=thispdf[1][ii]
			#
		#
		
	#myPdf=getlogpdf(ints, .1)
	#myPdf=getpdf(ints, .1)
	#
	# normalize the pdf:
	totalY=0.0
	for y in pdf[1]:
		totalY+=y
	for i in xrange(len(pdf[1])):
		pdf[1][i]/=(float(totalY)*binres)
	myPdf=pdf
	
	if doShow:
		plt.figure(0)
		plt.semilogy(myPdf[0], myPdf[1])
		plt.show()
	#
	return myPdf

def glebspdf(T=2.0,c=1.0):
	X=[]
	Y=[]
	N=1000
	x=0
	for i in xrange(N):
		X+=[x]
		Y+=[( ((1.0-x/(1.0+T))**c) - ((1.0+x)**(-1.0-c)))/(x*log(1.0+T))]
		
		x+=T/float(N)
	
	plt.figure(0)
	plt.semilogy(X,Y)
	plt.show()
	
	return [X,Y]

def checkAnalyticalNHPP(Nits=100, Tmax=1.25, tao=1.0, cm=1.0, binres=.1):
	# note: there are a few screwy disaster mitigators in this script. the point here is just to check the NHPP sequence generator against
	# gleb's analytical results. depending on the input parameters, NHPP can produce a HUGE data-set and blows up from time to time.
	# the sequence appears to check out.
	#
	# note: the analytical result comes from gleb's paper: Yakovlev et al,
	# feb 2 2008, cond-mat.stat-mech (??), "Inter-arrival time distributios for the non-homogeneous Poisson Process"
	# the paper studies lambda(t) = c/(1+t). in our variables: tao -> 1/c, cm=1.0
	# obviously, the analytical 'result' will not line up for cm!=1.0.
	#
	T=Tmax
	c=1.0/float(tao)
	print "c,tao: %f, %f" % (c, tao)
	# now, check against analytic solution (as per gleb's paper):
	# pdf(x) = [(1-x/(1+T))**c - (1+x)**(-(1+c)] / (x*ln(1+T))
	# where L(t) = c/(1+t).
	# so, c - 1/tao	(rate)
	# cm = 1 (exponent)
	#
	myPdf=checkNHPP(Nits, Tmax, tao, cm, False, binres)
	# now, trim it up. let's say no more than 2*Tmax
	nmax=0
	for x in myPdf[0]:
		if x>Tmax*2.0: break
		nmax+=1
		
	anaPDF=[]
	#anaPDF2=[]
	for x in myPdf[0]:
		#if x==0:
		#	anaPDF+=[0]
		#	continue
		#val=(((1.0-x/(1.0+T) )**c) - ((1.0+x)**(-1.0-c)) )/(x*log(1.0+T))
		#anaPDF+=[val]
		if (1.0-x/(1.0+T))<0:
			anaPDF+=[0]
			continue
		anaPDF+=[( ((1.0-x/(1.0+T))**c) - ((1.0+x)**(-1.0-c)))/(x*log(1.0+T))]
		#anaPDF2+=[( ((1.0-x/(1.0+T))**tao) - ((1.0+x)**(-1.0-tao)))/(x*log(1.0+T))]
		

	#print "ingegralval: %f" % pdfIntegral
	#for i in xrange(len(anaPDF)):
	#	anaPDF[i]/=pdfIntegral
	
	fout=open('NHPP/pdf1.dat', 'w')
	for i in xrange(len(myPdf[0])):
		fout.write('%f\t%f\t%f\n' % (myPdf[0][i], myPdf[1][i], anaPDF[i]))
	fout.close()
	
	#
	plt.figure(0)
	plt.semilogy(myPdf[0][0:nmax], myPdf[1][0:nmax], '.')
	plt.semilogy(myPdf[0][0:nmax], anaPDF[0:nmax])
	#plt.semilogy(myPdf[0][0:nmax], anaPDF2[0:nmax], '-.')
	plt.show()
	return myPdf+[anaPDF]
	#return [myPdf[0], myPdf[1], anaPDF]

def getpdfFlat(data, binres=.1, ycol=1, maxbin=None):
	# wrapper for getpdf when data are in row format ([[x,y]...]
	X=[]
	for rw in data:
		X+=[rw[ycol]]

	return getpdf(X, binres, maxbin)
	

def getpdf(data, binres=.1, maxbin=None):
	#binres=.1
	#first, spin through data to get max value (and bin-size)
	# data are in [y] type list. if original data are in an array like: [[x,y,z], ...], then use getpdfFlat()
	#print "datalen: %d" % len(data)
	if len(data)==0: return [[0], [0]]
	
	maxy=data[0]
	miny=maxy
	for y in data:
		if y>maxy: maxy=float(y)
		if y<miny: miny=float(y)
	#print "minY, maxY: %f, %f" % (miny, maxy)
	# but for now, keep this simple:
	miny=0
	nBins=1+long((maxy-miny)/binres)
	#print "nBins: %d" % nBins
	#
	# yThreshold: the minimum possible y-value, aka where n=1; y_min=n/Nevents. exclude these values,
	# maybe exclude all vals < a/Nevents, where a is some constant of our choosing.
	yThreshold=2.0*1.0/float(len(data))
	#print "yThreshold: %f" % yThreshold
	#print "bins: %f, %f, %d" % (miny, maxy, nBins)
	#
	# initialize an array:
	pdf=[[],[]]
	#pdf[1]=[0]*nBins
	#print "initializing pdf with (%d) elements." % nBins
	#print "nBins: %d" % nBins
	if nBins>10**7: nBins=10**7
	for i in xrange(nBins):
		pdf[0]+=[miny + i*binres]
		pdf[1]+=[0]
		# if pdf[0][-1]>10000: break
	#
	totalEvents=0
	#print "make raw bins."
	for x in data:
		#if x>10000: continue
		thisbin=long(x/binres)
		if thisbin>=len(pdf[1]): continue
		pdf[1][thisbin]+=1
		totalEvents+=1
	totalEvents=float(totalEvents)*binres	# now, total events is an integral (multiply by bin width)
	#
	## now, normalize:
	#for i in xrange(len(pdf[0])):
	#	pdf[1][i]/=totalEvents
	##
	'''
	# now, trim the return set. pop off all the 0 value elements off the end. (???)
	# or maybe just remove all 0 value elements.
	i=0
	thisLen=len(pdf[1])
	newLen=0
	for i in xrange(len(pdf[0])):
		if pdf[1][i]>yThreshold: newLen+=1
	#
	print "reduced array length: %d" % newLen
	newPDF=[[0]*newLen,[0]*newLen]	# short pdf list...
	#
	newIndex=0
	#while i<thisLen:
	print "assign values to short PDF"
	for i in xrange(len(pdf[0])):
		if pdf[1][i]<=yThreshold:
			#pdf[0].pop(i)
			#pdf[1].pop(i)
			#thisLen=len(pdf[1])
			continue
		newPDF[0][newIndex]=pdf[0][i]
		newPDF[1][newIndex]=pdf[1][i]
		newIndex+=1
	pdf=None
	
	print "return."
	'''
	
	return pdf
	#return newPDF

def getlogpdf(data, binres=.1, maxbin=None):
	#first, spin through data to get max value (and bin-size)
	# data are in [y] type list. if original data are in an array like: [[x,y,z], ...], then use getpdfFlat()
	maxy=data[0]
	miny=maxy
	for y in data:
		if y>maxy: maxy=float(y)
		if y<miny: miny=float(y)
	# but for now, keep this simple:
	miny=0
	nBins=1+long(log(maxy-miny)/binres)
	print "bins: %f, %f, %d" % (miny, maxy, nBins)
	#
	# initialize an array:
	pdf=[[],[]]
	for i in xrange(long(nBins)+1):
		pdf[0]+=[exp(miny + i*binres)]
		pdf[1]+=[0]
	#
	totalEvents=0
	for x in data:
		thisbin=long(log(x)/binres)
		if thisbin<0: thisbin=0
		#print "len(pdf), thisbin: %d, %d" % (len(pdf[1]), thisbin)
		pdf[1][thisbin]+=1
		#if x>1.1: continue
		totalEvents+=1
	#
	# now, normalize:
	for i in xrange(1,len(pdf[0])):
		pdf[1][i]/=(totalEvents*(pdf[0][i]-pdf[0][i-1]))
	
	return pdf

def xyplotEvents(catarray, xcol=0, ycol=1):
	# assumes flat array structure: [[a,b,c,d], [a,b,c,d]]
	X=[]
	Y=[]
	for rw in catarray:
		X+=[rw[xcol]]
		Y+=[rw[ycol]]
	plt.figure(0)
	plt.clf()
	plt.plot(X,Y, '.')
	plt.show()

def testDoParkfield2(incat='cats/parkfieldfullShocks.cat', minmag=1.5, outdir='images/parkfieldDev', doshow=False):
	# use this script to diagnose doParkfield2(). more or less, copy doParkfield2() code into this function; exit, return, plot, etc. as
	# desired.
	if os.system('ls %s' % outdir)!=0: os.system('mkdir %s' % outdir)
	
	rbi=intervalRecordBreaker(incat)
	doReverse=False
	if rbi.fullCat[0][0]>rbi.fullCat[-1][0]: doReverse=True
	evdt=datetime.datetime(2004,9,28, 17,15,24)
	maxFwdDate=datetime.datetime(2010,10,17)
	maxRevDate=datetime.datetime(1999,10,17)
	#rbi.setAftershockCatalog(incat, 40, 35.9, -120.5, .4, .15, datetime.datetime(2004,9,28, 17,15,24), datetime.datetime(2010,10,17), 0)
	#
	if doReverse==False:
		rbi.setAftershockCatalog(incat, 40, 35.9, -120.5, .4, .15, evdt, maxFwdDate, 0)
		print "Parkfield T=%d" % ((maxFwdDate-evdt).days)
	if doReverse==True:
		rbi.setAftershockCatalog(incat, 40, 35.9, -120.5, .4, .15, evdt, maxRevDate, 0)
		print "Parkfield T=%d" % ((maxRevDate-evdt).days)
	
	xyplotEvents(rbi.fullCat, 1,2)
	xyplotEvents(rbi.shockCat,1,2)
	#
	recs1=rbi.getRBintervals(minmag)

	return recs1
	

def doParkfield2(incat='cats/parkfieldfullShocks.cat', minmag=1.5, outdir='images/parkfield', doshow=False):
	# remember, we can give this a full (square) catalog. the rb-object, given parameters, can select elements along the fault after the event.
	#
	# create a second set of parkfield plots akin the the Hmine plots
	# clearly, this is a copy of doHminePlots2() with parkfield catalogs, coordinates, etc substituded where appropriate.
	#
	# here, we'll look at one or two mag-bins at different starting times, aka, 0, 10, 100, 200, x[0]
	# events after the mainshock.
	#
	# ultimately, this will be a singel function call to produce the pertinent Parkfield (data?) and plots. (right now, we start
	# from a catalog. we could produce the catalog from sql here as well for posterity).
	# 
	# catFname='parkcat.cat', theta=tTheta, clat=tLat, clon=tLon, ra=tA, rb=tB)
	# makeHmineCats(dtMinus=365.24*8, dtPlus=365.24*5, rx=.5, ry=.5, ra=.65, rb=.2, theta=67.4, fullcatout='cats/hminefull.cat', shockcat='cats/hmineshock.cat'):
	# t0=datetime.datetime(1999, 10, 16, 02, 46, 44)
	#outdir='images/parkfield/mag%d' % int(10*minmag)
	#
	rbPlots=[]	# a list to hold our plot data for fitting. this will be like: [ [[x1],[y1]], [[x2],[y2]], ...].
					# a data set will be rbPlots[i]; replotting a series goes like: plot(rbPlots[i][0], rpPlots[i][1])
	#
	#
	# note: parkfield is the default setting of the intervalRecordBreaker() class:
	rbi=intervalRecordBreaker(incat)
	# reset the catalog. note we're getting only events after the mainshock.
	# setAftershockCatalog(self, catFname=None, theta=tTheta, clat=35.9, clon=-120.5, ra=tA, rb=tB, eventDate=datetime.datetime(2004,9,28, 17,15,24), maxDate=datetime.datetime(2009,9,28, 17,15,24), skipSeconds=0)
	#
	'''
	tLat=35.9
	tLon=-120.5
	tTheta=40.0		#47?
	tA=.4		# ellipse axes
	tB=.15
	datetime.datetime(2004,9,28, 17,15,24)
	'''
	# set the shock catalog for immediately after parkfield, along fault.
	# note: we insert some max/min datetimes here. we could get rid of these and alway read to the end of the catalog.
	# catalog completeness might be a problem.
	doReverse=False
	if rbi.fullCat[0][0]>rbi.fullCat[-1][0]: doReverse=True
	evdt=datetime.datetime(2004,9,28, 17,15,24)
	maxFwdDate=datetime.datetime(2010,10,17)
	maxRevDate=datetime.datetime(1999,10,17)
	#rbi.setAftershockCatalog(incat, 40, 35.9, -120.5, .4, .15, datetime.datetime(2004,9,28, 17,15,24), datetime.datetime(2010,10,17), 0)
	#
	if doReverse==False:
		rbi.setAftershockCatalog(incat, 40, 35.9, -120.5, .4, .15, evdt, maxFwdDate, 0)
		#rbi.setAftershockCatalog(incat, 40, 35.9, -120.5, .15, .15, evdt, maxFwdDate, 0)
		print "Parkfield T=%d" % ((maxFwdDate-evdt).days)
	if doReverse==True:
		rbi.setAftershockCatalog(incat, 40, 35.9, -120.5, .4, .15, evdt, maxRevDate, 0)
		#rbi.setAftershockCatalog(incat, 40, 35.9, -120.5, .15, .15, evdt, maxRevDate, 0)
		print "Parkfield T=%d" % ((maxRevDate-evdt).days)
	
	#xyplotEvents(rbi.fullCat, 1,2)
	#xyplotEvents(rbi.shockCat,1,2)
	
	
	# this is our base catalog and RB object. now, run a RB series. when do we get our first RB event?
	recs1=rbi.getRBintervals(minmag) # returns [biggers, smallers] -> [ [ [dates],[interval],[nRecords],[nth earthquake (natural time)] ], [ [],[],[],[] ] ]
	dNmax=recs1[0][3][1]	#NT of the first (non-trivial) RB interval.
	
	#dNlist=[0, 4, 16, 64, 256, 1024]
	secdays=86400.0
	#currMag=1.25
	currMag=minmag
	intCurrMag=long(10*currMag)
	strmag=str(currMag).replace('.', '')
	
	#maglags = [[2.5, [2*secdays, 3*secdays], '.-'], [3.0, [2*secdays, secdays], '+-'], [3.5, [2*secdays, 1*secdays, .1*secdays], '-^'], [4.0, [2*secdays, 1*secdays,.1*secdays], '-h'] ]
	# 20*secdays -> 20*(number of seconds in a day) -> 20 days.
	maglags = [[currMag, [20*secdays, 10*secdays, 5*secdays, 3*secdays, 2*secdays, 1*secdays,.1*secdays, .01*secdays, 0.0], '-h'] ]
	maglags = [[currMag, [1*secdays,.1*secdays, .01*secdays, 0.0], '-h'] ]
	
	# clear figures:
	for i in xrange(8):
		plt.figure(i)
		plt.clf()
		plt.cla()
	#
	minx=0.0
	for rw in maglags:
		#
		thismag=rw[0]
		thismarker=rw[2]
		for dT in rw[1]:
			# it's messy but simple; make a whole bunch of data-files.
			dfileInts=open('%s/parkfieldInts%s-dT%s.dat' % (outdir, strmag, int(dT*100/secdays)), 'w')
			dfileIntsNT=open('%s/parkfieldIntsNT%s-dT%s.dat' % (outdir, strmag, int(dT*100/secdays)), 'w')
			dfileNRB=open('%s/parkfieldNRB%s-dT%s.dat' % (outdir, strmag, int(dT*100/secdays)), 'w')
			dfileNRBnt=open('%s/parkfieldNRBnt%s-dT%s.dat' % (outdir, strmag, int(dT*100/secdays)), 'w')
			#
			dfileInts.write("#parkfield event RB data, RB Interval Duration (time-time)\n#incat, minmag, dT\n#%s, %f, %f\n#t\tinterval\n" % (incat, minmag, dT/secdays))
			dfileIntsNT.write("#parkfield event RB data, RB Interval Duration (NT)\n#incat, minmag, dT\n#%s, %f, %f\n#n\tinterval\n" % (incat, minmag, dT/secdays))
			dfileNRB.write("#parkfield event RB data, NRB (time-time)\n#incat, minmag, dT\n#%s, %f, %f\n#t\tnrb\n" % (incat, minmag, dT/secdays))
			dfileNRBnt.write("#parkfield event RB data, NRB (NT)\n#incat, minmag, dT\n#%s, %f, %f\n#n\tnrb\n" % (incat, minmag, dT/secdays))
			
			dfileIntsSmall=open('%s/parkfieldIntsSmall%s-dT%s.dat' % (outdir, strmag, int(dT*100/secdays)), 'w')
			dfileIntsNTSmall=open('%s/parkfieldIntsNTSmall%s-dT%s.dat' % (outdir, strmag, int(dT*100/secdays)), 'w')
			dfileNRBSmall=open('%s/parkfieldNRBSmall%s-dT%s.dat' % (outdir, strmag, int(dT*100/secdays)), 'w')
			dfileNRBntSmall=open('%s/parkfieldNRBntSmall%s-dT%s.dat' % (outdir, strmag, int(dT*100/secdays)), 'w')
			#
			dfileIntsSmall.write("#parkfield event RBsmall data, RB Interval Duration (time-time)\n#incat, minmag, dT\n#%s, %f, %f\n#t\tinterval\n" % (incat, minmag, dT/secdays))
			dfileIntsNTSmall.write("#parkfield event RBsmall data, RB Interval Duration (NT)\n#incat, minmag, dT\n#%s, %f, %f\n#n\tinterval\n" % (incat, minmag, dT/secdays))
			dfileNRBSmall.write("#parkfield event RBsmall data, NRB (time-time)\n#incat, minmag, dT\n#%s, %f, %f\n#t\tnrb\n" % (incat, minmag, dT/secdays))
			dfileNRBntSmall.write("#parkfield event RBsmall data, NRB (NT)\n#incat, minmag, dT\n#%s, %f, %f\n#n\tnrb\n" % (incat, minmag, dT/secdays))
			#
			#
			thismarker=rw[2]
			#if dT==10*secdays: thismarker='-+'
			rbi.setAftershockCatalog(incat, 40, 35.9, -120.5, .4, .15, datetime.datetime(2004,9,28, 17,15,24), datetime.datetime(2010,10,17), dT)
			recs2=rbi.getRBintervals(thismag)
			# NT-NRB:
			plt.figure(0)
			x=[]
			y=[]
			for i in xrange(len(recs2[0][3])):
				if log10(recs2[0][3][i])<minx: continue
				x+=[recs2[0][3][i]]
				y+=[recs2[0][2][i]]
				dfileNRBnt.write("%d\t%d\n" % (x[-1], y[-1]))
			#
			rbPlots+=[[recs2[0][3], recs2[0][2], [thismag, dT, thismarker]]]
			plt.loglog(x, y, thismarker, label='m%s,dT%s' % (str(thismag)[0:4], float(dT)/secdays))
			#
			# now, get NT-RBintervals (interval duration):
			plt.figure(1)
			x=[]
			y=[]
			for i in xrange(len(recs2[0][3])):
				if log10(recs2[0][3][i])<minx: continue
				x+=[recs2[0][3][i]]
				y+=[recs2[0][1][i]]
				dfileIntsNT.write("%d\t%f\n" % (x[-1], y[-1]))
			plt.loglog(x, y, thismarker, label='m%s,dT%s' % (str(thismag)[0:4], float(dT)/secdays))
			#
			# NRB time-time:
			plt.figure(2)
			x=[]
			y=[]
			firstdt=datetimeToFloat(recs2[0][0][0])
			print "firstdt: %f" % firstdt
			for i in xrange(len(recs2[0][3])):
				#if log10(recs2[0][3][i])<minx: continue
				x+=[abs(datetimeToFloat(recs2[0][0][i])-firstdt)]
				y+=[abs(recs2[0][2][i])]
				dfileNRB.write("%f\t%d\n" % (x[-1], y[-1]))
			plt.loglog(x, y, thismarker, label='m%s,dT%s' % (str(thismag)[0:4], float(dT)/secdays))
			#
			# interval duration, time-time
			plt.figure(3)
			x=[]
			y=[]
			firstdt=datetimeToFloat(recs2[0][0][0])
			for i in xrange(len(recs2[0][3])):
				#if log10(recs2[0][3][i])<minx: continue
				x+=[abs(datetimeToFloat(recs2[0][0][i])-firstdt)]
				y+=[recs2[0][1][i]]
				dfileInts.write("%f\t%f\n" % (x[-1], y[-1]))
			plt.loglog(x, y, thismarker, label='m%s,dT%s' % (str(thismag)[0:4], float(dT)/secdays))
			#
			# small records:
			# NT-NRB:
			plt.figure(4)
			x=[]
			y=[]
			for i in xrange(len(recs2[1][3])):
				if log10(recs2[1][3][i])<minx: continue
				x+=[recs2[1][3][i]]
				y+=[recs2[1][2][i]]
				dfileNRBntSmall.write("%d\t%d\n" % (x[-1], y[-1]))
			#
			rbPlots+=[[recs2[1][3], recs2[1][2], [thismag, dT, thismarker]]]
			plt.loglog(x, y, thismarker, label='m%s,dT%s' % (str(thismag)[0:4], float(dT)/secdays))
			#
			# now, get NT-RBintervals (interval duration):
			plt.figure(5)
			x=[]
			y=[]
			for i in xrange(len(recs2[1][3])):
				if log10(recs2[1][3][i])<minx: continue
				x+=[recs2[1][3][i]]
				y+=[recs2[1][1][i]]
				dfileIntsNTSmall.write("%d\t%f\n" % (x[-1], y[-1]))
			plt.loglog(x, y, thismarker, label='m%s,dT%s' % (str(thismag)[0:4], float(dT)/secdays))
			#
			# NRB time-time:
			plt.figure(6)
			x=[]
			y=[]
			firstdt=datetimeToFloat(recs2[1][0][0])
			print "firstdt: %f" % firstdt
			for i in xrange(len(recs2[1][3])):
				#if log10(recs2[1][3][i])<minx: continue
				x+=[abs(datetimeToFloat(recs2[1][0][i])-firstdt)]
				y+=[recs2[1][2][i]]
				dfileNRBSmall.write("%f\t%d\n" % (x[-1], y[-1]))
			plt.loglog(x, y, thismarker, label='m%s,dT%s' % (str(thismag)[0:4], float(dT)/secdays))
			#
			# interval duration, time-time
			plt.figure(7)
			x=[]
			y=[]
			firstdt=datetimeToFloat(recs2[1][0][0])
			for i in xrange(len(recs2[1][3])):
				#if log10(recs2[1][3][i])<minx: continue
				x+=[abs(datetimeToFloat(recs2[1][0][i])-firstdt)]
				y+=[recs2[1][1][i]]
				dfileIntsSmall.write("%f\t%f\n" % (x[-1], y[-1]))
			plt.loglog(x, y, thismarker, label='m%s,dT%s' % (str(thismag)[0:4], float(dT)/secdays))
			#
			dfileInts.close()
			dfileIntsNT.close()
			dfileNRB.close()
			dfileNRBnt.close()
			#
			dfileIntsSmall.close()
			dfileIntsNTSmall.close()
			dfileNRBSmall.close()
			dfileNRBntSmall.close()
			#
			print "SeriesLen(%d): %d" % (dT, len(recs2[0][2]))
	
	strFmt='eps'
	plt.figure(0)	
	plt.title("Number of Large RB Events (NT)\n(Parkfield)")
	#plt.xlabel("n events since mainshock")
	plt.xlabel("n events since RB start")
	plt.ylabel("Number of Record Breaking Events")
	plt.legend(loc='upper left')
	#plt.legend(loc='lower right')
	plt.savefig("%s/newRBintervalsNTParkfielddTdM-%d-%d.%s" % (outdir, dT, intCurrMag, strFmt), format='%s' % strFmt)
	#
	plt.figure(1)
	plt.title("Large RB Interval Duration (NT)\n(Parkfield)")
	#plt.xlabel("n events since mainshock")
	plt.xlabel("n events since RB start")
	plt.ylabel("Interval Duration (days)")
	plt.legend(loc='upper left')
	#plt.legend(loc='lower right')
	plt.savefig("%s/RBintervalDurationsNTParkfielddTdM-%d-%d.%s" % (outdir, dT, intCurrMag, strFmt), format='%s' % strFmt)
	#
	plt.figure(2)
	plt.title("Number of Large RB Events (time)\n(Parkfield)")
	#plt.xlabel("n events since mainshock")
	plt.xlabel("days since main event")
	plt.ylabel("Number of Record Breaking Events")
	plt.legend(loc='upper left')
	#plt.legend(loc='lower right')
	plt.savefig("%s/newRBintervalsParkfielddTdM-%d-%d.%s" % (outdir, dT, intCurrMag, strFmt), format='%s' % strFmt)
	#
	plt.figure(3)
	plt.title("Large RB Interval Duration (time)\n(Parkfield)")
	#plt.xlabel("n events since mainshock")
	plt.xlabel("days since main event")
	plt.ylabel("Interval Duration (days)")
	plt.legend(loc='upper left')
	#plt.legend(loc='lower right')
	#strFmt='svg'
	plt.savefig("%s/RBintervalDurationsParkfielddTdM-%d-%d.%s" % (outdir, dT, intCurrMag, strFmt), format='%s' % strFmt)
	#
	plt.figure(4)	
	plt.title("Number of Small RB Events (NT)\n(Parkfield)")
	#plt.xlabel("n events since mainshock")
	plt.xlabel("n events since RB start")
	plt.ylabel("Number of Record Breaking Events")
	plt.legend(loc='upper left')
	#plt.legend(loc='lower right')
	plt.savefig("%s/newRBintervalsNTSmallParkfielddTdM-%d-%d.%s" % (outdir, dT, intCurrMag, strFmt), format='%s' % strFmt)
	#
	plt.figure(5)
	plt.title("Small RB Interval Duration (NT)\n(Parkfield)")
	#plt.xlabel("n events since mainshock")
	plt.xlabel("n events since RB start")
	plt.ylabel("Interval Duration (days)")
	plt.legend(loc='lower left')
	#plt.legend(loc='lower right')
	plt.savefig("%s/RBintervalDurationsNTSmallParkfielddTdM-%d-%d.%s" % (outdir, dT, intCurrMag, strFmt), format='%s' % strFmt)
	#
	plt.figure(6)
	plt.title("Number of Small RB Events (time)\n(Parkfield)")
	#plt.xlabel("n events since mainshock")
	plt.xlabel("days since main event")
	plt.ylabel("Number of Record Breaking Events")
	plt.legend(loc='upper left')
	#plt.legend(loc='lower right')
	plt.savefig("%s/newRBintervalsSmallParkfielddTdM-%d-%d.%s" % (outdir, dT, intCurrMag, strFmt), format='%s' % strFmt)
	#
	plt.figure(7)
	plt.title("Small RB Interval Duration (time)\n(Parkfield)")
	#plt.xlabel("n events since mainshock")
	plt.xlabel("days since main event")
	plt.ylabel("Interval Duration (days)")
	plt.legend(loc='lower left')
	#plt.legend(loc='lower right')
	plt.savefig("%s/RBintervalDurationsSmallParkfielddTdM-%d-%d.%s" % (outdir, dT, intCurrMag, strFmt), format='%s' % strFmt)
	#
	#
	
	# and let's do some data fitting:
	# for now, fit t>=10**-1.5 in the Nrecords plot, maybe interval>10**-2 for RBinterval??
	# alternatively, always trim the first (maybe first 2) data points.
	#
	# we want the log/log fit, so we can either construct a log-based error or we can take the log-log
	# of the data and do a linear fit.
	#
	# for now, let's take a break on the fits:
	'''
	fitPlotsInt=[]	# fit RB interval (mags)
	fitPlotsN=[]	# fit to Nrb plots.
	slopes=[]
	plt.figure(0)
	plt.clf()
	plt.cla()
	#plt.plot(NevLog, YevLog, '.')
	for pset in rbPlots:
		x=[]
		y=[]
		yfit=[]
		#print len(pset[0])
		for i in xrange(len(pset[0])):
			# let's hard-code some fit-range conditions. maybe we can be more dynamic later...
			#if log10(pset[1][i])<.6: continue
			if log10(pset[0][i])<minx: continue
			x+=[log10(pset[0][i])]
			y+=[log10(pset[1][i])]
		if len(x)<=1: continue	
		#
		# now, we have a log-log set of one dataset. fit it to a line...
		#print "x: %s" % str(x)
		#print "y: %s" % str(y)
		p=scipy.array([0,1])
		#print "x: %s" % str(x)
		
		plsq=spo.leastsq(linRes, p, args=(scipy.array(y), scipy.array(x)), full_output=0, maxfev=20000)	# (function to minimize, initial prams, argument-arrays, max-iterations)
		slopes+=[plsq[0][1]]
		# make an array of the fit function. note we can plot this or get a chi-sqr bc we have data-points.
		for X in x:
			fitval=plsq[0][0] + X*plsq[0][1]
			yfit+=[fitval]
		plt.plot(x,yfit, pset[2][2], label='M%s,dT=%s,m=%s' % (pset[2][0], float(pset[2][1])/secdays, str(slopes[-1])[0:4] ))
	plt.title("Fit to Nrb (NT)\n(Parkfield)")
	plt.xlabel("log10(nevents since mainshock)")
	plt.ylabel("log10(Nrb)")
	#plt.legend(loc='lower right')
	plt.legend(loc='upper left')
	plt.savefig("%sNrbFitsNTParkfieldNTdTdM-%d-%d.pdf" % (outdir, dT, intCurrMag), format='pdf')
	#plt.legend(loc='upper left')
	#
	#plt.figure(2)
	#plt.title("Hector Mine Interval Slopes (NT)")
	#plt.xlabel("mag")
	#plt.ylabel("slope")
	#plt.plot(mags, slopes, '.-')
	#plt.savefig("%sRBintervalSlopesNT.pdf" % outdir)
	print "slopes: %s" % slopes
	'''
	#
	#return recs1
	#
	
		
	#rbi.plotRBintervalSet(minmag, maxmag, dmag, 'images/hmineshock/')
	
	#

	
def makeParkfieldFaultCatalogSql(startDate=datetime.datetime(1999,9,28, 17,15,24), endDate=datetime.datetime(2009,9,28, 17,15,24), fullcatout='cats/parkfieldfull10yrs.cat', shockcatout='cats/parkfield10yrs.cat', doRev=False):
	# def setAftershockCatalog(self, catFname=None, theta=tTheta, clat=35.9, clon=-120.5, ra=tA, rb=tB, eventDate=datetime.datetime(2004,9,28, 17,15,24), maxDate=datetime.datetime(2009,9,28, 17,15,24)):
	# get a catalog, make a catalog file:
	import _mysql
	import MySQLdb
	#
	sqlHost = 'localhost'
	sqlUser = 'myoder'
	sqlPassword = 'yoda'
	sqlPort = 3306
	sqlDB = 'QuakeData'
	con=MySQLdb.connect(host=sqlHost, user=sqlUser, passwd=sqlPassword, port=sqlPort, db=sqlDB)
	c1=con.cursor()
	sqlstr='select date(eventDateTime) as dt, time(eventDateTime) as tm, lat, lon, mag from parkfieldquakes where eventDateTime between \'1999-09-27\' and \'2009-09-29\' order by eventDateTime asc'
	if doRev==True: sqlstr='select date(eventDateTime) as dt, time(eventDateTime) as tm, lat, lon, mag from parkfieldquakes where eventDateTime between \'1999-09-27\' and \'2009-09-29\' order by eventDateTime desc'
	c1.execute(sqlstr)
	fout=open(fullcatout, 'w')
	fout.write("#parkfield catalog generated from:\n#%s\n" % sqlstr)
	rw1=c1.fetchone()
	while rw1!=None:
		# spin through the cursor; write a catalog. note formatting choices...
		thisdt=str(rw1[0]).replace('-', '/')
		#thistm=str(rw1[1])
		fout.write("%s\t%s\t%f\t%f\t%f\n" % (thisdt, str(rw1[1]), rw1[2], rw1[3], rw1[4]))
		rw1=c1.fetchone()
	fout.close()
	c1.close()
	con.close()
	# now we have a catalog of the parkfield area (note, it is partially defined by our "parkfieldquakes" MySQL view.
	#
	#makeShockCat(incat, outcat)
	makeShockCat(fullcatout, shockcatout)

def makeHmineCats(dtMinus=365.24*8, dtPlus=365.24*5, rx=.75, ry=.75, ra=.65, rb=.2, theta=67.4, fullcatout='cats/hminefull.cat', shockcat='cats/hmineshock.cat'):
	#
	# hmine event occured at: 34.36, -116.116, 1999-10-16 02:46:44
	t0=datetime.datetime(1999, 10, 16, 02, 46, 44)
	tstart=floatToDateTime(datetimeToFloat(t0)-dtMinus)
	tend=floatToDateTime(datetimeToFloat(t0)+dtPlus)
	lat0=34.594
	lon0=-116.271
	#
	# get a catalog centered around this coordinate -dtMinus, +dtPlus, and a rectangle with sides 2*rx x 2*ry
	# then get a "shock cat" from the ellipse centered on the epicenter and dfined by a=ra, b=rb, theta.
	#
	# get the catalog from sql.
	import _mysql
	import MySQLdb
	#
	sqlHost = 'localhost'
	sqlUser = 'myoder'
	sqlPassword = 'yoda'
	sqlPort = 3306
	sqlDB = 'QuakeData'
	con=MySQLdb.connect(host=sqlHost, user=sqlUser, passwd=sqlPassword, port=sqlPort, db=sqlDB)
	c1=con.cursor()
	dtmStartDate=datetime.datetime
	sqlstr='select date(eventDateTime) as dt, time(eventDateTime) as tm, lat, lon, mag from Earthquakes where catalogID=517 and eventDateTime between \'%s\' and \'%s\' and lat between %f and %f and lon between %f and %f and mag>=2.5 order by eventDateTime asc' % (str(tstart), str(tend), lat0-ry, lat0+ry, lon0-rx, lon0+rx)
	#
	#print sqlstr
	# make a square-catalog for starters.
	fullout=open(fullcatout,'w')
	c1.execute(sqlstr)
	fullout.write("#Hector Mine Catalog created from:\n#%s\n" % sqlstr)
	rw1=c1.fetchone()
	while rw1!=None:
		thisdt=str(rw1[0]).replace('-', '/')
		fullout.write("%s\t%s\t%f\t%f\t%f\n" % (thisdt, str(rw1[1]), rw1[2], rw1[3], rw1[4]))
		rw1=c1.fetchone()	
	fullout.close()
	c1.close()
	con.close()
	#
	# make an elliptical shock-cat:
	# eyeballing the catalog, hmine aftershocks live along a line (-116.5,35.2); (-116,34) which gives us an angle of 67.4 deg cw from west,
	# ra=.65 deg (or so); rb=.3 (just eyeballing for now)
	# copy "makeShockCat()" to get an ellipitical dataset...
	#theta=67.4
	#ra=.65
	#rb=.3
	#
	rbp=intervalRecordBreaker()
	rbp.setAftershockCatalog(fullcatout, theta, lat0, lon0, ra, rb, tstart, tend)
	fout=open(shockcat, 'w')
	fout.write("#hector mine shock-cat (events along fault defined by x0=(%f, %f), theta=%f, ra=%f, rb=%f\n" % (lon0, lat0, theta, ra, rb))
	for rw in rbp.shockCat:
		fout.write("%d/%d/%d\t%d:%d:%d.%d\t%f\t%f\t%f\n" %  (rw[0].year, rw[0].month, rw[0].day, rw[0].hour, rw[0].minute,rw[0].second, rw[0].microsecond, float(rw[1]), float(rw[2]), float(rw[3])))
		
	fout.close()
	
	return rbp
	
	#

def makeShockCat(incat, outcat):
	# now, load this catalog into an IntervalRecordBreaker() object, set the shockCat, and output a shock-catalog to shockcatout
	#catFname='parkcat.cat', theta=tTheta, clat=tLat, clon=tLon, ra=tA, rb=tB
	rbp=intervalRecordBreaker()
	#rbp=intervalRecordBreaker(fullcatout, 40.0, 35.9, -120.5, .4, .15)
	#print "aftershock prams: %s, %f, %f, %f, %f, %f, %s, %s" % (fullcatout, rbp.tTheta, rbp.tLat, rbp.tLon, rbp.tA, rbp.tB, datetime.datetime(1999,9,28), datetime.datetime(2009,9,29))
	#rbp.setAftershockCatalog(fullcatout, rbp.tTheta, rbp.tLat, rbp.tLon, rbp.tA, rbp.tB, datetime.datetime(1999,9,28), datetime.datetime(2009,9,29))
	rbp.setAftershockCatalog(incat, rbp.tTheta, rbp.tLat, rbp.tLon, rbp.tA, rbp.tB, datetime.datetime(1999,9,28), datetime.datetime(2009,9,29))
	
	# note, we have now set the full catalog, then transformed the catalog along the fault and selected only those in the ellipse defined by center lat,lon and ra, rb.
	# the latter, transformed catalog, is rbp.shockCat; the format is [[datetime(eventdate), lat, lon, mag, x`, y`]..]
	# where x` and y` are the x,y coordinates along the transformed (fault) axis, and may be listed here in reverse order.
	# 
	# now, write a new catalog file:
	fout=open(outcat, 'w')
	fout.write("#parkfield \'shock-cat\', events along fault defined by: x0(lat,lon)=%f, %f, theta=%f, ra=%f, rb=%f\n" % (rbp.tLat, rbp.tLon, rbp.tTheta, rbp.tA, rbp.tB))
	for rw in rbp.shockCat:
		fout.write("%d/%d/%d\t%d:%d:%d.%d\t%f\t%f\t%f\n" %  (rw[0].year, rw[0].month, rw[0].day, rw[0].hour, rw[0].minute,rw[0].second, rw[0].microsecond, float(rw[1]), float(rw[2]), float(rw[3])))
		
	fout.close()
	
	return rbp
	

#
def plotPFNHPP():
	import recordBreaker as nhrb2
	sdays=86400.0
	#
	rbp=nhrb2.recordbreaker()
	ints1=rbp.getNHPPintervalsOmori1b(365, 38.066/sdays, .1313)
	ints2=rbp.getNHPPintervalsOmori1b(365, 43.6874/sdays, .0547)
	ints3=rbp.getNHPPintervalsOmori1b(365, 50.139/sdays, .0228)
	ints4=rbp.getNHPPintervalsOmori1b(365, 57.5433/sdays, .0095)
	ints5=rbp.getNHPPintervalsOmori1b(365, 66.0411/sdays, .00395)
	#
	import matplotlib.pyplot as plt
	# ints# -> [[n], [T], [interval]]
	#
	# time-time plots:
	thisPlt0=plt.semilogy
	thisPlt1=plt.plot
	thisPlt2=plt.plot
	
	plt.figure(0)
	plt.clf()
	thisPlt0(ints1[1], ints1[2], label=".8.066, .1313")
	thisPlt0(ints2[1], ints2[2], label="43.6874, .0547")
	thisPlt0(ints3[1], ints3[2], label="50.139, .0228")
	#thisPlt0(ints4[1], ints4[2], label="57.5433, .0095")
	#thisPlt0(ints5[1], ints5[2], label="66.0411, .00395")
	plt.title("NHPP Interoccurence Intervals (time-time)")
	plt.xlabel("elapsed time (days)")
	plt.ylabel("Interval Duration (days)")
	plt.legend(loc='upper left')
	#
	plt.figure(1)
	plt.clf()
	plt.semilogy(ints1[0], ints1[2], label=".8.066, .1313")
	plt.semilogy(ints2[0], ints2[2], label="43.6874, .0547")
	plt.semilogy(ints3[0], ints3[2], label="50.139, .0228")
	plt.semilogy(ints4[0], ints4[2], label="57.5433, .0095")
	plt.semilogy(ints5[0], ints5[2], label="66.0411, .00395")
	plt.title("NHPP Interoccurence Intervals (Natural Time)")
	plt.xlabel("n-events (natural time)")
	plt.ylabel("Interval Duration (days)")
	plt.legend(loc='lower right')
	#
	plt.figure(2)
	plt.clf()
	thisPlt2(ints1[0], ints1[1], label=".8.066, .1313")
	thisPlt2(ints2[0], ints2[1], label="43.6874, .0547")
	thisPlt2(ints3[0], ints3[1], label="50.139, .0228")
	thisPlt2(ints4[0], ints4[1], label="57.5433, .0095")
	thisPlt2(ints5[0], ints5[1], label="66.0411, .00395")
	plt.title("NHPP Interoccurence Intervals")
	plt.xlabel("n-events (natural time)")
	plt.ylabel("Total Elapsed Time (days)")
	plt.legend(loc='lower right')
	#
	plt.show()

############################################
############################################
# paper-production scripts:
def doPF2set():
	# the whole doParkfield2() set. this should be everything we need for publication.
	mags=[1.0, 1.5, 2.0, 2.5, 3.0]
	thiscat='cats/parkfieldfullShocks.cat'
	#
	for mg in mags:
		print "getting parkfiled RB events for mag=%f" % mg
		magstr=str(mg).replace('.', '')
		outdir='images/parkfield/mag%s' % magstr
		if os.system('ls %s' % outdir)!=0: os.system('mkdir %s' % outdir)
		doParkfield2(thiscat, mg, outdir, False)
	#
	# now, time-reverse the catalog and do it again.
	# we need a full or reverse-shock cat. use the full-cat; setshock() will cull events as necessary.
	print "*************************\ndo parkfield-reverse\n*************************"
	thiscat='cats/parkcatRev.cat'
	#if os.system('ls images/parkfield/reverse')!=0: os.system('mkdir images/parkfield/reverse')
	for mg in mags:
		print "getting parkfiled reverse RB events for mag=%f" % mg
		magstr=str(mg).replace('.', '')
		outdir='images/parkfield/revmag%s' % magstr
		if os.system('ls %s' % outdir)!=0: os.system('mkdir %s' % outdir)
		doParkfield2(thiscat, mg, outdir, False)
	

def doHmine2set():
	# the whole doParkfield2() set. this should be everything we need for publication.
	mags=[2.75, 3.0, 3.25, 3.5]
	thiscat='cats/hminefull.cat'
	#
	for mg in mags:
		print "getting Hmine RB events for mag=%f" % mg
		magstr=str(mg).replace('.', '')
		outdir='images/hmineshock/mag%s' % magstr
		if os.system('ls %s' % outdir)!=0: os.system('mkdir %s' % outdir)
		doHminePlots2(thiscat, mg, outdir)



def NHPP3set(Nits=10, binsize=1.0):
	# get NHPP RB series plots for the full range of Robert's parkfield m_c/tao values.
	#Nits=10
	if os.system('ls images/NHPPcat3')!=0: os.system('mkdir images/NHPPcat3')
	#
	#Tmax=1000
	#Tmax=365.0
	Tmax=2209.28	#T when we monitor parkfield up to (more or less) today...
	#binsize=1.0
	#binsize = .5
	binsizeNT=1.0
	secdays=86400.0
	#
	
	# from robert's parkfield:
	tao=38.066/secdays			#86400 sec/day
	cm=.1313
	print "tao, cm: %f, %f" % (tao, cm)
	newdir='mc10'
	doNHPPrecords3(Nits, Tmax, tao, cm, binsize, binsizeNT, 2700.0, False)
	if os.system('ls images/NHPPcat3/%s' % newdir)!=0: os.system('mkdir images/NHPPcat3/%s' % newdir)
	os.system("mv images/NHPPcat3/NHPPints.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPPintsNT.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPP-NRB.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPP-NRBnt.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/scatter.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedint.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedintNT.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRB.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRBnt.dat images/NHPPcat3/%s" % newdir)
	#
	os.system("mv images/NHPPcat3/binnedintNT-log2.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRBnt-log2.dat images/NHPPcat3/%s" % newdir)
	#	
	tao=43.6874/secdays
	cm=.0547
	print "tao, cm: %f, %f" % (tao, cm)
	doNHPPrecords3(Nits, Tmax, tao, cm, binsize, binsizeNT, 1700, False)
	newdir='mc15'
	if os.system('ls images/NHPPcat3/%s' % newdir)!=0: os.system('mkdir images/NHPPcat3/%s' % newdir)
	os.system("mv images/NHPPcat3/NHPPints.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPPintsNT.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPP-NRB.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPP-NRBnt.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/scatter.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedint.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedintNT.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRB.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRBnt.dat images/NHPPcat3/%s" % newdir)
	#
	os.system("mv images/NHPPcat3/binnedintNT-log2.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRBnt-log2.dat images/NHPPcat3/%s" % newdir)
	#
	
	tao=50.139/secdays
	cm=.0228
	print "tao, cm: %f, %f" % (tao, cm)
	doNHPPrecords3(Nits, Tmax, tao, cm, binsize, binsizeNT, 700, False)
	newdir='mc20'
	if os.system('ls images/NHPPcat3/%s' % newdir)!=0: os.system('mkdir images/NHPPcat3/%s' % newdir)
	os.system("mv images/NHPPcat3/NHPPints.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPPintsNT.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPP-NRB.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPP-NRBnt.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/scatter.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedint.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedintNT.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRB.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRBnt.dat images/NHPPcat3/%s" % newdir)
	#
	os.system("mv images/NHPPcat3/binnedintNT-log2.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRBnt-log2.dat images/NHPPcat3/%s" % newdir)
	#
	
	tao=57.5433/secdays
	cm=.0095
	print "tao, cm: %f, %f" % (tao, cm)
	doNHPPrecords3(Nits, Tmax, tao, cm, binsize, binsizeNT, 250, False)
	newdir='mc25'
	if os.system('ls images/NHPPcat3/%s' % newdir)!=0: os.system('mkdir images/NHPPcat3/%s' % newdir)
	os.system("mv images/NHPPcat3/NHPPints.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPPintsNT.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPP-NRB.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPP-NRBnt.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/scatter.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedint.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedintNT.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRB.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRBnt.dat images/NHPPcat3/%s" % newdir)
	#
	os.system("mv images/NHPPcat3/binnedintNT-log2.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRBnt-log2.dat images/NHPPcat3/%s" % newdir)
	
	tao=66.0411/secdays
	cm=.00395
	print "tao, cm: %f, %f" % (tao, cm)
	doNHPPrecords3(Nits, Tmax, tao, cm, binsize, binsizeNT, 75, False)
	newdir='mc30'
	if os.system('ls images/NHPPcat3/%s' % newdir)!=0: os.system('mkdir images/NHPPcat3/%s' % newdir)
	os.system("mv images/NHPPcat3/NHPPints.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPPintsNT.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPP-NRB.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPP-NRBnt.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/scatter.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedint.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedintNT.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRB.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRBnt.dat images/NHPPcat3/%s" % newdir)
	#
	os.system("mv images/NHPPcat3/binnedintNT-log2.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRBnt-log2.dat images/NHPPcat3/%s" % newdir)
	#	
	# and the "dimensionless" version:
	# we want a "dimensionless", universal version, so:
	# r = (1/tao) * 1/(1+t/cm))**p
	# for now p=1
	# rearrange:
	# r*tao = 1/(1+t/c)**p
	# R==r*tao; T=t/c:
	# R = 1/(1+T)**p (but for now, p=1):
	# R = 1/(1+T).
	# but this is not really universal; it just sets the characteristic scale to 1, doesn't it?
	# in any case, we can accomplish this by using parameter values:
	'''
	#
	# some general versions:
	tao=1.0
	cm=1.0
	print "tao, cm: %f, %f" % (tao, cm)
	doNHPPrecords3(Nits, Tmax, tao, cm, binsize, binsizeNT, 50, False)
	newdir='ND1-1'
	if os.system('ls images/NHPPcat3/%s' % newdir)!=0: os.system('mkdir images/NHPPcat3/%s' % newdir)
	os.system("mv images/NHPPcat3/NHPPints.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPPintsNT.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPP-NRB.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPP-NRBnt.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/scatter.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedint.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedintNT.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRB.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRBnt.dat images/NHPPcat3/%s" % newdir)
	'''
	#
	tao=1.0
	cm=10.0
	print "tao, cm: %f, %f" % (tao, cm)
	doNHPPrecords3(Nits, Tmax, tao, cm, binsize, binsizeNT, 75, False)
	newdir='ND1-10'
	if os.system('ls images/NHPPcat3/%s' % newdir)!=0: os.system('mkdir images/NHPPcat3/%s' % newdir)
	os.system("mv images/NHPPcat3/NHPPints.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPPintsNT.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPP-NRB.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPP-NRBnt.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/scatter.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedint.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedintNT.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRB.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRBnt.dat images/NHPPcat3/%s" % newdir)
	#
	os.system("mv images/NHPPcat3/binnedintNT-log2.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRBnt-log2.dat images/NHPPcat3/%s" % newdir)
	#
	#
	tao=1.0
	cm=100.0
	print "tao, cm: %f, %f" % (tao, cm)
	doNHPPrecords3(Nits, Tmax, tao, cm, binsize, binsizeNT, 100, False)
	newdir='ND1-100'
	if os.system('ls images/NHPPcat3/%s' % newdir)!=0: os.system('mkdir images/NHPPcat3/%s' % newdir)
	os.system("mv images/NHPPcat3/NHPPints.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPPintsNT.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPP-NRB.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPP-NRBnt.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/scatter.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedint.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedintNT.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRB.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRBnt.dat images/NHPPcat3/%s" % newdir)
	#
	os.system("mv images/NHPPcat3/binnedintNT-log2.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRBnt-log2.dat images/NHPPcat3/%s" % newdir)
	#
	#
	tao=1.0
	cm=1000.0
	print "tao, cm: %f, %f" % (tao, cm)
	doNHPPrecords3(Nits, Tmax, tao, cm, binsize, binsizeNT, 500, False)
	newdir='ND1-1000'
	if os.system('ls images/NHPPcat3/%s' % newdir)!=0: os.system('mkdir images/NHPPcat3/%s' % newdir)
	os.system("mv images/NHPPcat3/NHPPints.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPPintsNT.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPP-NRB.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPP-NRBnt.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/scatter.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedint.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedintNT.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRB.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRBnt.dat images/NHPPcat3/%s" % newdir)
	#
	os.system("mv images/NHPPcat3/binnedintNT-log2.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRBnt-log2.dat images/NHPPcat3/%s" % newdir)
	#
	tao=.01
	cm=.01
	print "tao, cm: %f, %f" % (tao, cm)
	doNHPPrecords3(Nits, Tmax, tao, cm, binsize, binsizeNT, 30, False)
	newdir='ND001-001'
	if os.system('ls images/NHPPcat3/%s' % newdir)!=0: os.system('mkdir images/NHPPcat3/%s' % newdir)
	os.system("mv images/NHPPcat3/NHPPints.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPPintsNT.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPP-NRB.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPP-NRBnt.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/scatter.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedint.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedintNT.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRB.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRBnt.dat images/NHPPcat3/%s" % newdir)
	#
	os.system("mv images/NHPPcat3/binnedintNT-log2.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRBnt-log2.dat images/NHPPcat3/%s" % newdir)
	#
	tao=.01
	cm=.1
	print "tao, cm: %f, %f" % (tao, cm)
	doNHPPrecords3(Nits, Tmax, tao, cm, binsize, binsizeNT, 160, False)
	newdir='ND001-01'
	if os.system('ls images/NHPPcat3/%s' % newdir)!=0: os.system('mkdir images/NHPPcat3/%s' % newdir)
	os.system("mv images/NHPPcat3/NHPPints.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPPintsNT.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPP-NRB.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPP-NRBnt.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/scatter.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedint.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedintNT.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRB.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRBnt.dat images/NHPPcat3/%s" % newdir)
	#
	os.system("mv images/NHPPcat3/binnedintNT-log2.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRBnt-log2.dat images/NHPPcat3/%s" % newdir)
	#
	tao=.01
	cm=1.0
	print "tao, cm: %f, %f" % (tao, cm)
	doNHPPrecords3(Nits, Tmax, tao, cm, binsize, binsizeNT, 500, False)
	newdir='ND001-1'
	if os.system('ls images/NHPPcat3/%s' % newdir)!=0: os.system('mkdir images/NHPPcat3/%s' % newdir)
	os.system("mv images/NHPPcat3/NHPPints.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPPintsNT.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPP-NRB.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPP-NRBnt.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/scatter.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedint.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedintNT.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRB.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRBnt.dat images/NHPPcat3/%s" % newdir)
	#
	os.system("mv images/NHPPcat3/binnedintNT-log2.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRBnt-log2.dat images/NHPPcat3/%s" % newdir)
	#
	tao=.01
	cm=10.0
	print "tao, cm: %f, %f" % (tao, cm)
	doNHPPrecords3(Nits, Tmax, tao, cm, binsize, binsizeNT, 500, False)
	newdir='ND001-10'
	if os.system('ls images/NHPPcat3/%s' % newdir)!=0: os.system('mkdir images/NHPPcat3/%s' % newdir)
	os.system("mv images/NHPPcat3/NHPPints.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPPintsNT.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPP-NRB.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/NHPP-NRBnt.svg images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/scatter.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedint.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedintNT.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRB.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRBnt.dat images/NHPPcat3/%s" % newdir)
	#
	os.system("mv images/NHPPcat3/binnedintNT-log2.dat images/NHPPcat3/%s" % newdir)
	os.system("mv images/NHPPcat3/binnedNRBnt-log2.dat images/NHPPcat3/%s" % newdir)
	print "finished."
#
#
# an experimental PF run:
def doPF2setExp():
	# the whole doParkfield2() set. this should be everything we need for publication.
	mags=[1.0, 1.5, 2.0, 2.5, 3.0]
	thiscat='cats/parkfieldfullShocks.cat'
	#
	for mg in mags:
		print "getting parkfiled RB events for mag=%f" % mg
		magstr=str(mg).replace('.', '')
		outdir='images/parkfieldExp/mag%s' % magstr
		if os.system('ls %s' % outdir)!=0: os.system('mkdir %s' % outdir)
		doParkfield2(thiscat, mg, outdir, False)
	#
	# now, time-reverse the catalog and do it again.
	# we need a full or reverse-shock cat. use the full-cat; setshock() will cull events as necessary.
	print "*************************\ndo parkfield-reverse\n*************************"
	thiscat='cats/parkcatRev.cat'
	#if os.system('ls images/parkfield/reverse')!=0: os.system('mkdir images/parkfield/reverse')
	for mg in mags:
		print "getting parkfiled reverse RB events for mag=%f" % mg
		magstr=str(mg).replace('.', '')
		outdir='images/parkfieldExp/revmag%s' % magstr
		if os.system('ls %s' % outdir)!=0: os.system('mkdir %s' % outdir)
		doParkfield2(thiscat, mg, outdir, False)


#############################################################
# these 3 functions may have been strictly part of the Parkfield characteristic eq. project, now moved to parkfieldChar.py
def fPL(x,p):
	return (10**p[0])*(x**p[1])

def fitRes(p,y,x):
	# "residuals" function for scipy.optimize fitting.
	# assume unform weight.
	#err = (y-fPL(x,p))*wt	
	#err = (y-(10**p[0])*(x**p[1]))
	err=y-(p[0]+x*p[1])
	return err
	
def plotParkCat():
	# def plotMagsIntervals(self, doShow=True, doSave=False, saveName='catalogMagsInts', pltTitle='Parkfield Seismicity', nTime=False):
	# make a time series of intervals and magnitudes...
	rb1=intervalRecordBreaker()
	rb1.setAftershockCatalog('cats/parkfieldfull10yrs.cat', 40.0, 35.9, -120.5, .4, .15, datetime.datetime(2000,01,01))
	rb1.plotMagsIntervals(False, True, 'images/parkfield/parkfieldSeismicity-timetime', 'Parkfield Seismicity (time)', False)
	rb1.plotMagsIntervals(False, True, 'images/parkfield/parkfieldSeismicity-NT', 'Parkfield Seismicity (natural time)', True)

#############################################################
	
