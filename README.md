# OFDM-Sim
A quick MATLAB simulation of an OFDM communications system. The script was written in MATLAB to facilitate it getting running quickly as an OFDM refresher exercise. As such, I didn't take any particular efforts to keep this code as clean as I normally do. However, if you're here to learn about OFDM systems, this code should work nicely for your purposes.

The system works by generating a pulse vector and a random bit stream. I've included several pulse generation scripts, but I use the root raised cosine pulse in the simulation (https://en.wikipedia.org/wiki/Root-raised-cosine_filter). The bit stream is either subjected to either an IFFT or FFT depending on whether we're dealing with the transmitter or receiver. The root raised cosine pulse is applied to the result to form the baseband signal. The main advantage is to doing this is that once you take into account the sending and receiving ends, the spectrum is that of a raised cosine, which is one of several options that result in no intersymbol interference. There are two carriers that can convey +1 and -1, resulting in QPSK.

The script outputs several plots to show experimental and theoretical results. These include the spectrum of the signal and bit error rates as a function of SNR. You can also see the signal at various points in the system, which might be educationally useful to someone down the road.



