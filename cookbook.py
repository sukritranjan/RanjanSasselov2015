# -*- coding: iso-8859-1 -*-

"""
This file contains functions used for manipulation of binned data. Contained here are both prev
"""

def rebin_uneven(lefts_old, rights_old, vals_old, lefts_new, rights_new):
	"""
	Purpose of this script is to rebin a histogram onto a new x-axis. It's goal is to accomodate irregular bins. 
	
	This method AVERAGES bins together. So if you have two bins containing the numbers 2 and 3 and want to rebin them into one bin, the value of that bin would be 2.5. This is well suited to dealing with things that are rates. For example, suppose you have those two bins. Suppose each contains the flux in ergs/cm2/s/nm: 10 and 20 for 200-210 and 210-220 nm, respectively. Then if you want the flux average from 200-220 nm, you need to compute (10*10+20*10)/(10+10)=15 erg/cm2/s/nm. 
	
	lefts_old=left edges of old bins
	rights_old=right edges of old bins
	vals_old=values of old bins
	lefts_new=left edges of new bins
	rights_new=right edges of new bins
	
	returns: values of new bins
	Uses linear interpolation to split up bins.
	ASSUMES: abscissa monotonically increasing ###
	The bins defined by lefts_new and rights_new must be entirely enclosed in the bins defined by left_old and right_old
	"""
	import numpy as np
	import pdb
	vals_new=np.zeros(np.shape(lefts_new))
	binwidths_old=rights_old-lefts_old#widths of old bins. May be nonuniform.
	
	for ind in range(0, len(vals_new)):
		#Loop over each of the new bins
		val_new=0. #initialize the value of the new bin to zero.
		left_new=lefts_new[ind]
		right_new=rights_new[ind]

		#is the new bin contained entirely within an existing bin?
		fullycontained_ind=np.where((left_new>=lefts_old) & (right_new<=rights_old)) #if there is a bin for which both these conditions hold, then this bin is contained entirely within it.

		if np.size(fullycontained_ind)>0: #if such a bin exists...
			val_new=vals_old[fullycontained_ind]# then the value of the new bin is just the enclosing old bin.
		else: #otherwise, we know that the new bin must span at least 1 junction of bins, meaning there will be a left partial part and a right partial part
			val_new_numerator=0 #adding up multiple contributions so need to keep track of the numerators and denominators of the weighted sum.
			val_new_denominator=0
			
			#Given that we are not wholely contained within a single bin, there may exist a partial bin enclosed by the left edge, i.e. left_edge must be in the middle of some bin. Only case in which this is *not* true is if the left edge _perfectly_ overlaps with one of the old bin edges.
			leftedgemiddle_ind=np.where((left_new>lefts_old) & (left_new<rights_old))
			if np.size(leftedgemiddle_ind)>0: #This condition is met as long as the left edge does not perfectly overlap with another edge
				val_new_numerator=val_new_numerator+vals_old[leftedgemiddle_ind]*(rights_old[leftedgemiddle_ind]-left_new)
				val_new_denominator=val_new_denominator+rights_old[leftedgemiddle_ind]-left_new
			
			#As there must be a left edge in the middle of a bin somwhere, there must be a right edge in the middle of a bin somewhere, unless it perfectly overlaps with an old edge.
			rightedgemiddle_ind=np.where((right_new<rights_old) & (right_new>lefts_old))
			if np.size(rightedgemiddle_ind)>0: #condition met as long as there is no perfect overlap in the edges
				val_new_numerator=val_new_numerator+vals_old[rightedgemiddle_ind]*(right_new-lefts_old[rightedgemiddle_ind])
				val_new_denominator=val_new_denominator+right_new-lefts_old[rightedgemiddle_ind]
			
			#Lastly, there may exist whole old bins encapsulated within the new bin. Add those in as well.
			oldencapsulatedinds=np.where((left_new<=lefts_old) & (right_new>=rights_old)) #indices of old bins which are completely encapsulated by the new bins. Note this captures the cases where a new edge aligns with an old edge
			if np.size(oldencapsulatedinds) > 0: #if such encapsulated bins exist then...
			      val_new_numerator=val_new_numerator+np.sum(vals_old[oldencapsulatedinds]*binwidths_old[oldencapsulatedinds])
			      val_new_denominator=val_new_denominator+np.sum(binwidths_old[oldencapsulatedinds])
			
			val_new=val_new_numerator/val_new_denominator
		vals_new[ind]=val_new
		#pdb.set_trace()
	return vals_new

#"""Test the uneven rebinning function we have written. """
#import numpy as np

#########First family of tests: evenly-binned input (old data)
#lefts_old=np.array([0., 1., 2., 3., 4., 5.])
#rights_old=np.array([1., 2., 3., 4., 5., 6.])
#vals_old=np.array([100., 110., 120., 130., 140., 150.]) 

##First test: ensure we recover the same bins if we return the same edges
#lefts_new=lefts_old
#rights_new=rights_old
#vals_new_1=rebin_uneven(lefts_old, rights_old, vals_old, lefts_new, rights_new)
#print vals_new_1-vals_old #should be 0 if we did it right

##Second test: bin down the data by a factor of two.
#lefts_new=np.array([0., 2., 4.])
#rights_new=np.array([2., 4., 6.])
#vals_new_2=rebin_uneven(lefts_old, rights_old, vals_old, lefts_new, rights_new)
#print vals_new_2-np.array([105., 125., 145.]) #should be 0 if we did it right

##Third test: evenly spaced values, but displaced from edges
#lefts_new=np.array([0.5, 1.5, 2.5, 3.5, 4.5])
#rights_new=np.array([1.5, 2.5, 3.5, 4.5, 5.5])
#vals_new_3=rebin_uneven(lefts_old, rights_old, vals_old, lefts_new, rights_new)
#print vals_new_3-np.array([105., 115., 125., 135., 145.])

##Fourth test: evenly spaced values, displaced from edges, reduce resolution by factor of 2
#lefts_new=np.array([0.7, 2.7])
#rights_new=np.array([2.7, 4.7])
#vals_new_4=rebin_uneven(lefts_old, rights_old, vals_old, lefts_new, rights_new)
#print vals_new_4-np.array([(0.3*100+1.0*110.+0.7*120.)/(0.3+1.0+0.7), (0.3*120+1.0*130+0.7*140.)/(0.3+1+0.7)]) #Test reveals numpy floating-point errors at the <1e-14 level.

##Fifth test: unevenly spaced values
#lefts_new=np.array([0.0, 1.0, 1.3, 1.5, 2.0, 3.5]) #same, left overlap/right in same bin, same bin, left in same bin/right overlap, left overlap/right in next bin, left in one bin/right overlap in next bin.
#rights_new=np.array([1.0,1.3, 1.5, 2.0, 3.5, 5.0])
#vals_new_5=rebin_uneven(lefts_old, rights_old, vals_old, lefts_new, rights_new)
#print vals_new_5-np.array([100.,110, 110, 110., (120+0.5*130)/1.5, (0.5*130.+140)/1.5])

#lefts_new=np.array([0.3, 1.4]) #both in different bins, both in different bins with multiple bins in the middle
#rights_new=np.array([1.4, 4.1])
#vals_new_5=rebin_uneven(lefts_old, rights_old, vals_old, lefts_new, rights_new)
#print vals_new_5-np.array([(0.7*100+0.4*110)/(1.1), (0.6*110+120+130+0.1*140)/(0.6+2+0.1)]) #Test reveals numpy floating-point errors at the <1e-14 level.




###########Second family of tests: unevenly-binned input
#lefts_old=np.array([0., 1., 3., 3.5,4.5])
#rights_old=np.array([1.,3., 3.5,4.5, 6.])
#vals_old=np.array([100., 115., 130., 135., 146.666666667]) 

##First test: even output
#lefts_new=np.array([0., 1., 2., 3., 4., 5.])
#rights_new=np.array([1., 2., 3., 4., 5., 6.])
#vals_new_6=rebin_uneven(lefts_old, rights_old, vals_old, lefts_new, rights_new)
#shouldbe_6=np.array([100., 115., 115., (0.5*130+0.5*135), (0.5*135.+0.5*146.666666667), 146.666666667])
#print vals_new_6-shouldbe_6
##print np.sum(vals_old*(rights_old-lefts_old))-np.sum(vals_new_6) #should equal zero if everything conserved.

##Last test: uneven input, uneven output
#lefts_old=np.array
#lefts_old=np.array([ 0., 1., 3.,  4.5, 4.9, 5.5, 6.0])
#rights_old=np.array([1., 3., 4.5, 4.9, 5.5, 6.0, 10.])
#vals_old=np.array([  0, 10., 20., 30., 40,  50,  60.]) 

#lefts_new=np.array([ 0.5, 1.5, 3.1, 3.9, 5.0, 7.5])
#rights_new=np.array([1.5, 3.1, 3.9, 5.0, 7.5, 10.])
#vals_new_7=rebin_uneven(lefts_old, rights_old, vals_old, lefts_new, rights_new)
#shouldbe_7=np.array([(0.5*0+0.5*10), (1.5*10+0.1*20)/1.6, 20., (0.6*20.+0.4*30+0.1*40.)/(1.1), (0.5*40+0.5*50+1.5*60)/2.5, 60.])
#print vals_new_7-shouldbe_7

def get_bin_edges(x_centers):
	"""
	This function is used to generate bin edges for histogram-type data for which only the bin centers are available (e.g., data extracted from plots). It does so using the following algorithm:
	-The boundary between points x_i and x_i+1 is set at the midpoint between them.
	-The left edge of bin x_0 is equidistant from the center of the bin as the right edge
	-The right edge of bin x_n is equidistant from the center of the bin as the left edge
	
	Input: x_centers, bin centers
	
	Output: x_lefts, x_rights
	"""
	import numpy as np
	x_lefts=np.zeros(np.shape(x_centers))
	x_rights=np.zeros(np.shape(x_centers))
	
	for ind in range(0, len(x_centers)-1):
		boundaryval=0.5*(x_centers[ind]+x_centers[ind+1])
		
		x_lefts[ind+1]=boundaryval
		x_rights[ind]=boundaryval
	
	x_lefts[0]=x_centers[0]-(x_rights[0]-x_centers[0])
	x_rights[-1]=x_centers[-1]+(x_centers[-1]-x_lefts[-1])
	
	return x_lefts, x_rights
