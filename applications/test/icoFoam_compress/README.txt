20180706: compression

- each time step p                      : 320000
- so over 10 timesteps                  : 3200000
- compressed 10 (lossy, 1e-3 presicion) : 926040

Binary writing (losless)
- 0.055/p               : 321153
- 0.055/p.gz            : 289072
- so lossy improves factor 3

Ascii writing (6 digits, lossy)
- with gzip single time : 121302 Jul  6 15:37 0.006/p.gz
- so lossy only improves 1.3 compared to gzip
