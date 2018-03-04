#!/usr/bin/env python

"""
Seqlib library for class project: 6-scientific-python
"""
import copy
import numpy as np
import pandas as pd
    
class Seqlib:

    """
    A seqlib object for class.
    """
    def __init__(self, ninds, nsites):

        ## generate the full sequence array
        self.ninds = ninds
        self.nsites = nsites
        self.seqs = self._simulate()
        
        ## store maf of the full seq array
        self.maf = self._get_maf()   

   
    ## private functions used only during init -----
    def _mutate(self, base):
        "converts a base to another base"
        diff = set("ACTG") - set(base)
        return np.random.choice(list(diff)) 
   

    def _simulate(self):  #ninds is a number of rows and nsites is a number of columns
       
        """
        The function simulate generate variable sequence data by creating mutations and sites with missing data. 
        """
        
        "returns a random array of DNA bases with mutations"
        oseq = np.random.choice(list("ACGT"), size=self.nsites) 
    
        #construct an array by creating rows of oseq
        arr = np.array([oseq for i in range(self.ninds)]) 
    
        #Binomial sampling through the array where the probability of one outcome is 0.1 (p=0.1). 
        #This will return an array of binary integers  
        muts = np.random.binomial(1, 0.1, (self.ninds, self.nsites)) 
    
        for col in range(self.nsites):      
            newbase = self._mutate(arr[0, col])    # creating a random mutation in the coulmns through the interation. 
            mask = muts[:, col].astype(bool) # muts flips a coin to assign outcome in binary integers(e.g. 0 or 1) which will 
                                                # then be converted to a boolean type using the astype() call.
            arr[:, col][mask] = newbase      # the arr[:,col] part pulls out a full column from the array. Then, the [mask] index pulls 
                                                # out some indices (defined above sell) from that column and stores in newbase
    
        # generate random missing data by using binomial sampling with probability of 0.1
        missing = np.random.binomial(1, 0.1, (self.ninds, self.nsites)) 
        #indicate missing data as "N"
        arr[missing.astype(bool)] = "N"  
        return arr
        
    def _get_maf(self):
        "returns the maf of the full seqarray while not counting Ns"
        ## init an array to fill and iterate over columns
        maf = np.zeros(self.nsites)
        for col in range(self.nsites):
            ## select this column of bases
            thiscol = self.seqs[:, col]

            ## mask "N" bases and get new length
            nmask = thiscol != "N"
            no_n_len = np.sum(nmask)

            ## mask "N" bases and get the first base
            first_non_n_base = thiscol[nmask][0]

            ## calculate maf of "N" masked bases
            freq = np.sum(thiscol[nmask] != first_non_n_base) / no_n_len
            if freq > 0.5:
                maf[col] = 1 - freq
            else:
                maf[col] = freq
        return maf


    def _filter_missing(self, maxmissing): 
          
        "returns a boolean filter True for columns with Ns > maxmissing"
        freqmissing = np.sum(self.seq == "N", axis=0) / self.seq.shape[0]    
        return freqmissing > maxmissing
    
    def _filter_maf(self, minmaf):
    
        "returns a boolean filter True for columns with maf < minmaf"
        return self.maf < minmaf

    ## public functions -----
    def filter(self, minmaf, maxmissing):
        """
        Applies maf and missing filters to the array 
        Parameters
        ----------
        minmaf: float
            The minimum minor allele frequency. Filter columns below this.
        maxmissing: float
            The maximum prop. missing data. Filter columns with prop Ns > this.
        """
        filter1 = self._filter_maf(minmaf)
        filter2 = self._filter_missing(maxmissing)
        fullfilter = filter1 + filter2
        return self.seqs[:, np.invert(fullfilter)]
    
    def calculate_statistics(self):
        """ 
        Returns a dataframe of statistics on the seqs array. The earlier 
        example from the notebook had a bug where var and inv were switched.
        """
        if self.seqs.size:
            nd = np.var(self.seqs == self.seqs[0], axis=0).mean()
            mf = np.mean(
                np.sum(self.seqs != self.seqs[0], axis=0) / self.seqs.shape[0])
            inv = np.all(self.seqs == self.seqs[0], axis=0).sum()
            var = self.seqs.shape[1] - inv
            return pd.Series(
                {"mean nucleotide diversity": nd,
                 "mean minor allele frequency": mf,
                 "invariant sites": inv,
                 "variable sites": var,
                })
        else:
            print("seqs array is empty")