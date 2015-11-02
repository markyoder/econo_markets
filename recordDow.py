from Scientific.Geometry import *
from math import *
import string
#from matplotlib import *
from pylab import *
import os
import random
import time

#
import _mysql
import MySQLdb

# make a sql connection:
sqlHost='localhost'
sqlPort=3306
#sqlUser='defaultInter'
#sqlPW = 'Int3r@ct1Ve'
sqlUser='myoder'
sqlPW='yoda'
#mycon = _mysql.connect(host=sqlHost, user=sqlUser, passwd=sqlPW, port=sqlPort, db='QuakeData')
con1 = MySQLdb.connect(host=sqlHost, user=sqlUser, passwd=sqlPW, port=sqlPort, db='markets')
con2  = MySQLdb.connect(host=sqlHost, user=sqlUser, passwd=sqlPW, port=sqlPort, db='markets')
c1=con1.cursor()
c2=con2.cursor()

# get the djia history:
djiaSelect='select intDt, seq, open, close, (open-close)/open from stocks where ticker=\'djia\' order by seq asc'
c1.execute(djiaSelect)
#
recordCrashes=[]	# this will be rows of tuples [date, deltaT, seq, mag]
recordMag = 0
lastSeq = 0

print 'spin djia cursor...'
while (1==1):
	row = c1.fetchone()		# or use .fetchall()
	if row==None:
		break		# end of cursor
	if abs((row[4]))>recordMag:
		print 'seq, dP, dT: %d, %f, %d' % (row[1], row[4], row[1]-lastSeq)
		recordCrashes=recordCrashes + [[row[0], row[1], row[1]-lastSeq, row[4]]]
		recordMag=row[4]
		lastSeq=row[1]
		
print recordCrashes

c1.close()
c2.close()
con1.close()
con2.close()
		

