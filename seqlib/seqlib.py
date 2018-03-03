#!/usr/bin/env python

"""
Seqlib class project for 6-scientific-python
"""

import numpy as np
import pandas as pd
    
class Seqlib:

    def __init__(self, ninds, nsites, maxfreq, minfreq):
        self.ninds = ninds
        self.nsites = nsites
        self.seqs = self.simulate()
        self.filter = self.filter_missing(), self.filter_maf()
    
    def simulate(self):  #ninds is a number of rows and nsites is a number of columns
        pass 
        """
        The function simulate generate variable sequence data by creating mutations and sites with missing data. 
        """
        
        #choose a random letter from the list "ACGT" for a number of times indicated in nsites
        oseq = np.random.choice(list("ACGT"), size=self.nsites) 
    
        #construct an array by creating rows of oseq
        arr = np.array([oseq for i in range(self.ninds)]) 
    
        #Binomial sampling through the array where the probability of one outcome is 0.1 (p=0.1). 
        #This will return an array of binary integers  
        muts = np.random.binomial(1, 0.1, (self.ninds, self.nsites)) 
    
        for col in range(self.nsites):      
            newbase = mutate(arr[0, col])    # creating a random mutation in the coulmns through the interation. 
            mask = muts[:, col].astype(bool) # muts flips a coin to assign outcome in binary integers(e.g. 0 or 1) which will 
                                                # then be converted to a boolean type using the astype() call.
            arr[:, col][mask] = newbase      # the arr[:,col] part pulls out a full column from the array. Then, the [mask] index pulls 
                                                # out some indices (defined above sell) from that column and stores in newbase
    
        # generate random missing data by using binomial sampling with probability of 0.1
        missing = np.random.binomial(1, 0.1, (self.ninds, self.nsites)) 
        #indicate missing data as "N"
        arr[missing.astype(bool)] = "N"  
        return arr
        
    def filter_missing(self):  #arr is the variable seq generated using the function, simulate.
        
        """
        The filter_missing function removes sequences that have more than centain frequency of missing data.
        """  
        # freqmissing is obtained by deviding the numner of N within the row by length of columns (shape[0]).
        freqmissing = np.sum(self.arr == "N", axis=0) / self.arr.shape[0]    
    
        #return rows with missing frequency less than maxfreq 
        return self.arr[:, freqmissing <= maxfreq]  
    
    def filter_maf(self):
    
        """
        The filter_maf function removes sequences that has too little minor allele frequencies. It uses copy becasue we do not want 
        to make changes in the original sequences.
        """  
    
        # The sum of variables in each column is divided by the length of columns (shape[0]) to obtein frequency and this view 
        # is stored in freqs. A number of variables is found by comparing the first row against the rest. 
        freqs = np.sum(self.arr != self.arr[0], axis=0) / self.arr.shape[0] 
    
        # store a copy to avoid modifying the original array 'arr'
        maf = freqs.copy()
    
        # subselect sites (columns) with major freq (>0.5) and modify to be 1-value (e.g 0.875 to 0.125)
        maf[maf > 0.5] = 1 - maf[maf > 0.5]
    
        # print only columns of minor allele frequencies greater than the minimum frequency (sequeces with too little mutations are eliminated here)
        return self.arr[:, maf > minfreq]  #print all row 
  
    
    def calculcate_statistics(self):
        """
        A. The calculcate_statistics function computes mean nucleotide diversity, mean minor allele frequency, variant sites, 
        variable sites
        """     
        # find mean nucleotide diversity: arr == arr[0] compares the first row against the rest of the rows. The similarity of 
        # each column is calculated and the similarity values from each column are used to compute the mean. 
        nd = np.var(self.arr == self.arr[0], axis=0).mean()
    
        # find mean minor allele frequency:  arr != arr[0] finds a number of differences between the first row against the rest
        # of the rows. The number of differeces is calculated for each column and divided by the total length of each column. 
        # The sum of the frequencies of all columns is used to compute the mean.
        mf = np.mean(np.sum(self.arr != self.arr[0], axis=0) / self.arr.shape[0])
   
        # find variant sites: create boolean mask for whether any sites in a column is different from the first row. Sum the 
        # number of columns with some variability 
        var = np.any(self.arr != self.arr[0], axis=0).sum()
    
        # find variable sites: substract the number of variant sites from the total number of columns
        inv = self.arr.shape[1] - var
    
        #Use Pandas to return all the values 
        return pd.Series(
            {"mean nucleotide diversity": nd,
             "mean minor allele frequency": mf,
             "invariant sites": inv,
             "variable sites": var,
            })
    