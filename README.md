# seqlib

The seqlib library is a module developed as a practice exercise for the PDSB class.

### Installation 
Clone the repo and use pip to install it. 
```bash
git clone https://github.com/kuratanp/seqlib.git
cd seqlib/
pip install .
```

### API Usage
To use the seqlib Python API:
```python

## import the seqlib library
import seqlib

## generate a seqlib class object with a sequence array of shape (n, m) 
s = seqlib.Seqlib(10, 100)

## return sequence array
print(s.seqs)

## return maf of sequence array
print(s.maf)

## return a filtered view of the seqarray filter based on maf
## and missing (N) sites
print(s.filter(minmaf=0.1, maxmissing=0.0)

## return a new copy of seqlib object with modified seqarray 
n = s.filter_seqlib(minmaf=0.1, maxmissing=0.0)

## view stats on the full seqarray
s.calculate_stats()

## view stats on the modified seqarray
n.calculate_stats()

## or do the same in one shot
s.filter_seqlib(minmaf=0.1, maxmissing=0.0).calculate_stats()




```

