"""Class doc"""

import math
import wave
import numpy
from warnings import warn

# Based on Kalker, Haitsma - "A Highly Robust Audio Fingerprinting System"
# Journal of New Music Research, Volume 32, Issue 2, 2003
# We assume that the amplitudes in wave files are linearly encoded


#Author - Srijan Magapu - Symphonium Project

class Fingerprinter(object):
    


    # actually calculate all sub-fingerprints
    def init_fingerprints(self):
        # read all bytes from wavefile
        data = numpy.frombuffer(self.wavefile.readframes(self.nsamples), dtype="i" + str(self.sample_width))
        
        # loop through the file, calculating the sub-fingerprint for each frame
        position = 0 + self.nsamples_between_frames
        there_are_frames_left = (position < (self.nsamples - self.nsamples_per_frame) )
        
        last_frame = data[0:self.nsamples_per_frame]
        fft_last_frame = numpy.fft.rfft(numpy.hamming(self.nsamples_per_frame) * last_frame)
        #fft_last_frame = numpy.fft.rfft( last_frame)
        index_fingerprint = 0
        
        powers = 2 ** numpy.arange(0, self.FINGERPRINT_NBITS)
        while there_are_frames_left:
            current_frame = data[position:position + self.nsamples_per_frame]
            # fourier transform current frame, weighted in time by hamming window
            fft_current_frame = numpy.fft.rfft(numpy.hamming(self.nsamples_per_frame) * current_frame)
            #fft_current_frame = numpy.fft.rfft(current_frame)
            
            # crunch crunch crunch
            self.fingerprints_binary[index_fingerprint] = self.generate_binary_fingerprint_from_frames(fft_last_frame, fft_current_frame)
            self.fingerprints[index_fingerprint] = numpy.sum(powers * self.fingerprints_binary[index_fingerprint])
            
            position = position + self.nsamples_between_frames
            there_are_frames_left = (position < (self.nsamples - self.nsamples_per_frame) )
            last_frame = current_frame
            fft_last_frame = fft_current_frame
            index_fingerprint = index_fingerprint + 1
            
            #print position
            #print current_frame.shape
            #print fft_current_frame.shape
            # this should print the index with the maximum pressure in this frame
            # index(frequency) =   length_of_interval_in_seconds * frequency
            #print numpy.argmax(fft_current_frame)
            #print current_frame
            #print "***"
            
            # Keep on going as long as there are enough samples left for a whole frame
        self.wavefile.close();
    
    
    
    ########################
    def print_info(self):
        """
        Prints some info about the file this Wrapper wraps and the wrapper itself.
        """
        print " * * * * * * * "
        print "Width of frames in seconds: " + str(self.framewidth)
        print "Overlap between frames: " + str(self.overlap)
        print "# of samples per frame: " + str(self.nsamples_per_frame)
        print "# of samples between beginning of adjacent frames: " + str(self.nsamples_between_frames)
        print " -> Resolution approx. " + str(self.resolution)
        print "# of fourier components between 300 and 2000 Hz: " + str(self.index_width_lower_to_upper)
        print " * * * * * * * "
        
    
        
    ##############################
    def calculate_frequency_bands(self, nsamples_per_frame, frame_length):
        # Since the paper cited above does not give the explicit spacing, we are forced to chose our own spacing.
        # Since apparently humans perceive loudness on a log scale, we use frequency bands that grow exponentially
        # in frequency to have an even "loudness" distribution in the bands
        
        # index of a frequency component in the FFT array =   length_of_frame_in_seconds * frequency
        # length of fft array : even n/2 + 1       odd (n+1)/2
        
        # length of fft array.
        if( nsamples_per_frame % 2 == 0): 
            nsamples_fft = 1 + nsamples_per_frame/2
        else:
            nsamples_fft = (nsamples_per_frame + 1)/2
        
        index_lower_limit = int(math.floor( frame_length * self.FREQUENCY_BAND_LOWER_LIMIT )) # index of 300Hz freq component
        index_upper_limit = int(math.ceil( frame_length * self.FREQUENCY_BAND_UPPER_LIMIT )) # index of 2000Hz freq component
        index_width = index_upper_limit - index_lower_limit + 1 # how many entries lie within 300Hz..2000Hz
        self.index_width_lower_to_upper = index_width
        
        # we must partition index_width into 33 subsets of width proportional to the exp() of the frequency they lie next to
        # we use f(n) = a + exp(b * n) where n is element [0,33] and f(0) == index_lower_limit, f(33) = index_upper_limit
        # a and b given here solve these equations
        a = float(index_lower_limit)
        b = numpy.log( float(index_upper_limit)/a ) / float(self.FINGERPRINT_NBITS + 1)
        # return indices for the frequency bands
        return numpy.round( a * numpy.exp(b * numpy.arange(0,self.FINGERPRINT_NBITS + 2)) )
        
        
        
    def generate_binary_fingerprint_from_frames(self, frame1, frame2):
        """
        Generate a sub-fingerprint for a frame
        """
        fingerprint = numpy.zeros(self.FINGERPRINT_NBITS, dtype=numpy.bool)
        
        energies1 = self.calculate_band_energies_from_frame(frame1)
        energies2 = self.calculate_band_energies_from_frame(frame2)
        for band in range(self.FINGERPRINT_NBITS):
            fingerprint[band] = energies2[band] - energies2[band+1] > energies1[band] - energies1[band+1]
            
        return fingerprint
        #self.frequency_band_boundary_indices
      
      
      
    def calculate_band_energies_from_frame(self, frame):
        """
        Calculate the sum of fourier components of a frame, by bands
        """
        energies = numpy.zeros(self.FINGERPRINT_NBITS + 1)
        for band in range(self.FINGERPRINT_NBITS + 1):
            energies[band] = numpy.sum ( numpy.abs( frame[ self.frequency_band_boundary_indices[band] : self.frequency_band_boundary_indices[band+1]  ] ) ** 2 )
        return energies
    
    
    
    def binary_distance(self, print1, print2):
        return numpy.sum(numpy.logical_xor(print1, print2))
    
    
    
    def block_distance(self, block1, block2):
        if block1.shape != block2.shape:
            raise "You are trying to calculate the distance between two different objects."
        block_dist = 0
        for i in range(len(block1)):
            block_dist += self.binary_distance(block1[i], block2[i])
        return block_dist
        
        
        
    # Find the best match between this fingerprinter and a given fingerprinter
    def find_position(self, fingerprinted):
        compare = fingerprinted.fingerprints_binary
        block_length = len(compare)
        if block_length > 256:
            block_length = 256
        
        min = 10000000000
        index = 0
        for i in range(len(self.fingerprints) - block_length):
            diff = self.block_distance(compare[:block_length], self.fingerprints_binary[i:i+block_length])
            if diff < min:
                print "found new minimum at " + str(i) + " (" + str(float(i) * self.nsamples_between_frames / self.sample_rate) + ") - " + str(diff)
                index = i
                min = diff
                
        print index # index of frame with lowest block distance
        print float(index) * self.nsamples_between_frames / self.sample_rate
