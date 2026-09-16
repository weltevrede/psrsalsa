# Copyright (c) 2015, Patrick Weltevrede
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# 
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from psrsalsa import io
import psrchive
import numpy as np

def convert_from_psrchive(psrchive_archive):

    nrbins = psrchive_archive.get_nbin()
    nsub = psrchive_archive.get_nsubint()
    nchan = psrchive_archive.get_nchan()
    npol = psrchive_archive.get_npol()

    psrsalsa_archive = io.create_empty(nsub, nrbins, nrfreqchans=nchan, nrpols=npol)

    psrname = psrchive_archive.get_source()
    io.set_psrname(psrsalsa_archive, psrname)


    # PSRCHIVE is broken. It returns galactic longitude for both l and b.
    #print(archive.get_Integration(0).get_galactic_longitude())
    #print(archive.get_Integration(0).get_galactic_latitude())
    #gal_l = archive.get_Integration(0).get_galactic_longitude()
    #gal_b = archive.get_Integration(0).get_galactic_latitude()
    #gc = SkyCoord(l=gal_l*u.degree, b=gal_b*u.degree, frame='galactic')
    #print(gc)
    #print("      Pulsar is at l=%f, b=%f"%(gal_l, gal_b))
    #print(archive.execute("edit int:point:gl"))
    #print(archive.execute("edit int:point:gb"))
    # Getting the RA/DEC is not supported in the python interface.
    # Can run psrsh commands via python interface, which allows the source coordinates to be extracted.
    coord = psrchive_archive.execute("edit coord").strip()
    ra = coord.split("=")[1].replace("+", "-").split("-")[0]
    ra_hour = int(ra.split(":")[0])
    ra_min = int(ra.split(":")[1])
    ra_sec = float(ra.split(":")[2])
    if "+" in coord:
         dec_sign = +1
    else:
         dec_sign = -1
    dec = coord.split("=")[1].replace("+", "-").split("-")[1]
    dec_deg = int(dec.split(":")[0])
    dec_min = int(dec.split(":")[1])
    dec_sec = float(dec.split(":")[2])
    #print("      %s ra=%s=%d:%d:%f dec=%s=%d:%d:%f"%(coord,ra,ra_hour,ra_min,ra_sec,dec,dec_sign*dec_deg,dec_min,dec_sec))
    ra = ra_hour + ra_min/60.0 + ra_sec/3600.0
    io.set_ra(psrsalsa_archive, ra*2*np.pi/24.0)
    dec = dec_deg + dec_min/60.0 + dec_sec/3600.0
    dec *= dec_sign
    io.set_dec(psrsalsa_archive, dec*np.pi/180.0)
    
  
    dm  = psrchive_archive.get_dispersion_measure()
    io.set_dm(psrsalsa_archive, dm)
    
    rm  = psrchive_archive.get_rotation_measure()
    io.set_rm(psrsalsa_archive, rm)

    isdedisp = psrchive_archive.get_dedispersed()
    io.set_dedisp_state(psrsalsa_archive, isdedisp)

    isdefarad = psrchive_archive.get_faraday_corrected()
    io.set_defarad_state(psrsalsa_archive, isdefarad)

    mjd = psrchive_archive.get_Integration(0).get_start_time().intday()
    mjd_frac = psrchive_archive.get_Integration(0).get_start_time().fracday()
    io.set_mjd_start(psrsalsa_archive, mjd+mjd_frac)
    
    #parang = archive.get_Integration(0).get_parallactic_angle()
    
    #obslength = psrchive_archive.integration_length()

    cfreq = psrchive_archive.get_centre_frequency()
    io.set_centre_freq(psrsalsa_archive, cfreq)
    bw = psrchive_archive.get_bandwidth()
    io.set_bw(psrsalsa_archive, bw)
    
    freq_table = []
    for subint in range(0, nsub):
        freq_table.append([])
        for chan in range(0, nchan):
            freq_table[subint].append(psrchive_archive.get_Integration(subint).weighted_frequency(chan,chan+1))
    #print(freq_table)
    io.set_freq_table(psrsalsa_archive, freq_table)

    tsub_table = []
    for subint in range(0, nsub):
        tsub_table.append(psrchive_archive.get_Integration(subint).get_duration())
    #print(tsub_table)
    io.set_tsub_table(psrsalsa_archive, tsub_table)

    polstate = psrchive_archive.get_state()
    if polstate == "Coherence":
        io.set_poltype(psrsalsa_archive, "Coherency")
    if polstate == "Stokes":
        io.set_poltype(psrsalsa_archive, "Stokes")

    basis = psrchive_archive.execute("edit rcvr:basis").strip()
    #print(basis)
    basis_str = basis.split("=")[1]
    if basis_str == "lin":
        io.set_feedtype(psrsalsa_archive, "Linear")
    if basis_str == "cir":
        io.set_feedtype(psrsalsa_archive, "Circular")

    period = psrchive_archive.get_Integration(0).get_folding_period()
    io.set_fold_period(psrsalsa_archive, period)
    io.set_tsamp(psrsalsa_archive, period/nrbins)

    ##cfreq = archive.get_Integration(0).get_centre_frequency()
    #cfreq = archive.get_Integration(0).weighted_frequency(0,archive.get_nchan())

    # bw = archive.get_Integration(0).get_bandwidth()
    #freq_low   = archive.get_Integration(0).get_centre_frequency()-0.5*bw
    #freq_high = archive.get_Integration(0).get_centre_frequency()+0.5*bw

    # Hmmm, should be a nicer way of doing this, but this is MEMORY_format, the data is loaded in memory.
    # But should be set to MEMORY_format, as otherwise can't use io.set_pulse_data()
    io.set_format(psrsalsa_archive, 99) 

    for subint in range(0, nsub):
        for chan in range(0, nchan):
            for pol in range(0, npol):
                profile = psrchive_archive.get_Profile(subint, pol, chan).get_amps()
                io.set_pulse_data(psrsalsa_archive, profile.tolist(), subint, pol, chan)
                
    return psrsalsa_archive
    
def load(filename):
    #print("Opening %s"%(filename))
    if "Archive_load" in dir(psrchive):
        psrchive_archive = psrchive.Archive_load(filename)
    else:
        psrchive_archive = psrchive.Archive.load(filename)
    return convert_from_psrchive(psrchive_archive)
