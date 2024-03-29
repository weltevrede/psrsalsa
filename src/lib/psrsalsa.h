/*
Copyright (c) 2015, Patrick Weltevrede
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <fitsio.h>
#include <psrsalsa_defines.h>
#include <psrsalsa_typedefs.h>
int isValidPSRDATA_format(int format);
void printPSRDataFormats(FILE *printdevice, int nrspaces);
int parsePSRDataFormats(char *cmd);
void cleanPSRData(datafile_definition *datafile, verbose_definition verbose);
int copy_params_PSRData(datafile_definition datafile_source, datafile_definition *datafile_dest, verbose_definition verbose);
int set_filename_PSRData(datafile_definition *datafile_dest, char *filename, verbose_definition verbose);
int set_psrname_PSRData(datafile_definition *datafile_dest, char *psrname, verbose_definition verbose);
int set_observatory_PSRData(datafile_definition *datafile_dest, char *observatory, verbose_definition verbose);
int set_institute_PSRData(datafile_definition *datafile_dest, char *institute, verbose_definition verbose);
int set_instrument_PSRData(datafile_definition *datafile_dest, char *instrument, verbose_definition verbose);
int set_scanID_PSRData(datafile_definition *datafile_dest, char *scanID, verbose_definition verbose);
int set_observer_PSRData(datafile_definition *datafile_dest, char *observer, verbose_definition verbose);
int set_projectID_PSRData(datafile_definition *datafile_dest, char *projectID, verbose_definition verbose);
int guessPSRData_format(char *filename, int noerror, verbose_definition verbose);
int openPSRData(datafile_definition *datafile, char *filename, int format, int enable_write, int read_in_memory, int nowarnings, verbose_definition verbose);
int closePSRData(datafile_definition *datafile, int perserve_header, int perserve_data, verbose_definition verbose);
void printHeaderPSRData(datafile_definition datafile, int update, verbose_definition verbose);
int readHeaderPSRData(datafile_definition *datafile, int readnoscales, int nowarnings, verbose_definition verbose);
int writeHeaderPSRData(datafile_definition *datafile, int argc, char **argv, int cmdOnly, char *notes, verbose_definition verbose);
int get_pointer_PulsePSRData(datafile_definition *datafile, long pulsenr, int polarization, int freq, int binnr, float **pulse_ptr, verbose_definition verbose);
int readPulsePSRData(datafile_definition *datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse, verbose_definition verbose);
int writePulsePSRData(datafile_definition *datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse, verbose_definition verbose);
int readPSRData(datafile_definition *datafile, float *data, verbose_definition verbose);
int writePSRData(datafile_definition *datafile, float *data, verbose_definition verbose);
int read_profilePSRData(datafile_definition datafile, float *profileI, int *zapMask, int polchan, verbose_definition verbose);
int read_partprofilePSRData(datafile_definition datafile, float *profileI, int *zapMask, int polchan, long nskip, long nread, verbose_definition verbose);
int read_rmsPSRData(datafile_definition datafile, float *rms, float *avrg, int *zapMask, pulselongitude_regions_definition *regions, int invert, int polchan, int freqchan, verbose_definition verbose);
int convert_if_uniform_frequency_spacing(datafile_definition *datafile, int nowarnings, verbose_definition verbose);
int force_uniform_frequency_spacing(datafile_definition *datafile, verbose_definition verbose);
void cleanVerboseState(verbose_definition *verbose_state);
void copyVerboseState(verbose_definition verbose_state_src, verbose_definition *verbose_state_dst);
int PSRDataHeader_parse_commandline(datafile_definition *psrdata, int argc, char **argv, verbose_definition verbose);
void printHeaderCommandlineOptions(FILE *printdevice);
void printHeaderGentypeOptions(FILE *printdevice);
int make_clone(datafile_definition original, datafile_definition *clone, verbose_definition verbose);
void swap_orig_clone(datafile_definition *original, datafile_definition *clone, verbose_definition verbose);
void printGenType(int gentype, FILE *destination);
char *returnGenType_str(int gentype);
char *returnFileFormat_str(int format);
int showHistory(datafile_definition datafile, verbose_definition verbose);
int skipallhashedlines(datafile_definition *datafile);
int setITRFlocation_by_name(datafile_definition *datafile, char *observatory, verbose_definition verbose);
int get_period(datafile_definition datafile, long subint, double *period, verbose_definition verbose);
double get_tsamp(datafile_definition datafile, long subint, verbose_definition verbose);
int convert_to_fixed_tsamp(datafile_definition *datafile, verbose_definition verbose);
double get_pulse_longitude(datafile_definition datafile, long subint, long binnr, verbose_definition verbose);
double get_tsub(datafile_definition datafile, long subint, verbose_definition verbose);
double get_tobs(datafile_definition datafile, verbose_definition verbose);
long double get_mjd_subint(datafile_definition datafile, long subint, verbose_definition verbose);
int get_channelbandwidth(datafile_definition datafile, double *channelbw, verbose_definition verbose);
void set_channelbandwidth(datafile_definition *datafile, double channelbw, verbose_definition verbose);
double get_bandwidth(datafile_definition datafile, verbose_definition verbose);
int set_bandwidth(datafile_definition *datafile, double bw, verbose_definition verbose);
double get_centre_frequency(datafile_definition datafile, verbose_definition verbose);
void set_centre_frequency(datafile_definition *datafile, double freq, verbose_definition verbose);
double get_nonweighted_channel_freq(datafile_definition psrdata, long channel, verbose_definition verbose);
double get_weighted_channel_freq(datafile_definition psrdata, long subint, long channel, verbose_definition verbose);
int set_weighted_channel_freq(datafile_definition *psrdata, long subint, long channel, double freq, verbose_definition verbose);
char * get_history_notes_last(datafile_definition *data);
int rebinPulse(float *Ipulse, long NrBins, float *Ipulse2, long NrBins2, int noDependencyWarning, verbose_definition verbose);
int continuous_shift(datafile_definition fin, datafile_definition *fout, int shift, int circularShift, char *output_name, int oformat, int argc, char **argv, verbose_definition verbose, int verbose2);
int data_parang(datafile_definition data, long subintnr, double *parang, verbose_definition verbose);
int check_baseline_subtracted(datafile_definition data, verbose_definition verbose);
char *str_replace_header_params(datafile_definition data, char *text, verbose_definition verbose);
void str_list_replace_keys(int nrspaces);
int rmSynthesis(datafile_definition data, float rm_low, float rm_high, float **rmsynth_array, int nrrmsteps, pulselongitude_regions_definition *onpulse, verbose_definition verbose);
void collapseRMSynthesisArray(float *rmsynth_array, int nrrmsteps, int nrBins, pulselongitude_regions_definition onpulse, float *singlespectrum, verbose_definition verbose);
int rmSynthesis_instrument_responds(int nrFreqChan, double chanbw, double cfreq0, double rm_low, double rm_high, float **rmsynth_responds, double **rmsynth_responds_double, int usedouble, int nrrmsteps, double rmshift, int callmulti, verbose_definition verbose);
int rmSynthesis_fitInstrumentalResponds(float *singlespectrum, float rmmin, float rmmax, int nrrmsteps, float *rm, float *offset, float *scale, int nrFreqChan, float chanbw, float cfreq0, float ftol, verbose_definition verbose);
void print_fitsio_version_used(FILE *stream);
void psrfits_set_noweights(int val);
void psrfits_set_absweights(int val);
void psrfits_set_use_weighted_freq(int val);
int filterPApoints(datafile_definition *datafile, verbose_definition verbose);
int readPPOLfile(datafile_definition *datafile, float *data, int extended, float add_longitude_shift, verbose_definition verbose);
int writePPOLfile(datafile_definition datafile, float *data, int extended, int onlysignificantPA, int twoprofiles, float PAoffset, verbose_definition verbose);
int make_paswing_fromIQUV(datafile_definition *datafile, int extended, int spstat, float sigma_limit, int sigmaI, pulselongitude_regions_definition onpulse, int normalize, int correctLbias, int correctPbias, float correctQV, float correctV, int nolongitudes, float loffset, float paoffset, datafile_definition *rms_file, float rebin_factor, int onpulseonly, verbose_definition verbose);
int paswing_remove_observed_PA_swing(datafile_definition *datafile, datafile_definition datafile_reference, int add, verbose_definition verbose);
int make_polarization_projection_map(datafile_definition datafile, float *map, int nrx, int nry, float background, int binnr, pulselongitude_regions_definition onpulse, int weighting, float threshold, int projection, float rot_long, float rot_lat, float conalselection, datafile_definition *subtract_pa_data, verbose_definition verbose);
int make_pa_distribution(datafile_definition datain, datafile_definition *dataout, int nrbins, int normalise, int weighttype, datafile_definition *pamask, float pamask_value, int ellipticity, verbose_definition verbose);
int preprocess_make_profile(datafile_definition original, datafile_definition *profile, int stokesI, verbose_definition verbose);
int preprocess_addNoise(datafile_definition original, datafile_definition *clone, float rms, verbose_definition verbose);
int preprocess_shuffle(datafile_definition original, datafile_definition *clone, int fixseed, verbose_definition verbose);
int preprocess_stokes(datafile_definition *original, verbose_definition verbose);
int preprocess_coherency(datafile_definition *original, verbose_definition verbose);
int preprocess_rotateStokes(datafile_definition *original, datafile_definition *clone, int inplace, int subint, float angle, float *angle_array, int stokes1, int stokes2, verbose_definition verbose);
int preprocess_addsuccessivepulses(datafile_definition original, datafile_definition *clone, long nrpulses, int complete, verbose_definition verbose);
int preprocess_dedisperse(datafile_definition *original, int undo, int update, double freq_ref, verbose_definition verbose);
int preprocess_deFaraday(datafile_definition *original, int undo, int update, double freq_ref, double *rm_table, verbose_definition verbose);
int preprocess_changeRefFreq(datafile_definition *original, double freq_ref_new, verbose_definition verbose);
int preprocess_addsuccessiveFreqChans(datafile_definition original, datafile_definition *clone, long nrfreq, int *fzapMask, verbose_definition verbose);
int preprocess_rebin(datafile_definition original, datafile_definition *clone, long NrBins, verbose_definition verbose);
int preprocess_debase(datafile_definition *original, pulselongitude_regions_definition *onpulse, float **baseline, int remove_shape, verbose_definition verbose);
int preprocess_channelselect(datafile_definition original, datafile_definition *clone, long chanelnr, verbose_definition verbose);
int preprocess_pulsesselect(datafile_definition original, datafile_definition *clone, long nskip, long nread, verbose_definition verbose);
int preprocess_blocksize(datafile_definition original, datafile_definition *clone, int blocksize, verbose_definition verbose);
int preprocess_fftshift(datafile_definition original, long singlesubint, float shiftPhase, int addslope, float slope, verbose_definition verbose);
int preprocess_polselect(datafile_definition original, datafile_definition *clone, long polnr, verbose_definition verbose);
int preprocess_transposeRawFBdata(datafile_definition original, datafile_definition *clone, verbose_definition verbose);
int preprocess_norm(datafile_definition original, float normvalue, pulselongitude_regions_definition *onpulse, int global, verbose_definition verbose);
int preprocess_clip(datafile_definition original, float clipvalue, verbose_definition verbose);
int preprocess_scale(datafile_definition original, float factor, float offset, verbose_definition verbose);
int preprocess_checknan(datafile_definition original, int generate_warning, verbose_definition verbose);
int preprocess_removenan(datafile_definition original, verbose_definition verbose);
int preprocess_checkinf(datafile_definition original, int generate_warning, verbose_definition verbose);
int preprocess_corrParAng(datafile_definition *original, datafile_definition *clone, int undo, verbose_definition verbose);
void print_fftw_version_used(FILE *stream);
int rotateSinglepulse(float *data, int npts, float epsilon, verbose_definition verbose);
int crosscorrelation_fft(float *data1, float *data2, int ndata, float *cc, int remove_baseline, verbose_definition verbose);
int crosscorrelation_fft_padding_cclength(int ndata, int extrazeropad);
int crosscorrelation_fft_padding(float *data1, float *data2, int ndata, int extrazeropad, float **cc, int *cclength, int remove_baseline, verbose_definition verbose);
int calcLRFS(float *data, long nry, long nrx, unsigned long fft_size, float *lrfs, int subtractDC, float *avrg_offpulse_lrfs_power, float *phase_track, float *phase_track_phases, int calcPhaseTrack, float freq_min, float freq_max, int track_only_first_region, float *subpulseAmplitude, int calcsubpulseAmplitude, int mask_freqs, int inverseFFT, pulselongitude_regions_definition *regions, float *var_rms, int argc, char **argv, verbose_definition verbose);
void calcModindex(float *lrfs, float *profile, long nrx, unsigned long fft_size, unsigned long nrpulses, float *sigma, float *rms_sigma, float *modind, float *rms_modind, pulselongitude_regions_definition *regions, float var_rms, float *avrg_offpulse_lrfs_power, float *avrg_mod, int avrg_mod_squared, verbose_definition verbose);
int calc2DFS(float *data, long nry, long nrx, unsigned long fft_size, float *twodfs, pulselongitude_regions_definition *onpulse, int region, int allow_cat_offpulse, int noise_subtract_mode, float *var_rms, verbose_definition verbose);
int foldP3(float *data, long nry, long nrx, float *map, int nr_p3_bins, float foldp3, int refine, int cyclesperblock, int noSmooth, float smoothWidth, float slope, float subpulse_offset, pulselongitude_regions_definition *onpulse
    , verbose_definition verbose);
long double calcDMDelay(long double freq, long double freq_ref, int inffrq, long double dm);
float calcRMAngle(float freq, float freq_ref, int inffrq, float rm);
void tempo2_ITRF_to_GRS80(double obs_X, double obs_Y, double obs_Z, double *longitude, double *latitude, double *height);
void
tempo2_GRS80_to_ITRF(double longitude, double latitude, double height, double *obs_X, double *obs_Y, double *obs_Z);
double observatory_long_geocentric(datafile_definition datafile);
double observatory_lat_geocentric(datafile_definition datafile);
double observatory_long_geodetic(datafile_definition datafile);
double observatory_lat_geodetic(datafile_definition datafile);
double observatory_height_geodetic(datafile_definition datafile);
void mjd2dateString(long double mjd, char *string, int precision, int type, char *separator);
void converthms(char *hms, double *h);
void converthms_string(char *hms, long double number, int precision, int type);
double calc_parang(double longitude, double latitude, double ra, double dec, double mjd, int precess);
int calc_precess_nut_ab(char system, double mjd, double *ra, double *dec, int nutation, int aberration, int verbose);
void projectionHammerAitoff_xy(float longitude, float latitude, float dlongitude, float dlatitude, float *x, float *y);
int projection_sphere_xy(float longitude, float latitude, float dlongitude, float dlatitude, float *x, float *y, float *weight);
int projection_sphere_longlat(float x, float y, float dlongitude, float dlatitude, float *longitude, float *latitude);
void projection_longlat_xy(float longitude, float latitude, float dlongitude, float dlatitude, float *x, float *y);
float derotate_deg(float a);
float derotate_180(float a);
double derotate_180_double(double a);
double derotate_180_rad_double(double a);
float derotate_90(float a);
double derotate_90_double(double a);
double derotate_180_small_double(double a);
double derotate_90_small_double(double a);
float polar_angle_rad(float x, float y);
float paswing(float alpha, float beta, float l, float pa0, float l0, int nrJumps, float *jump_longitude, float *jump_offset, float add_height_longitude, float add_height_shift, int add_height_shift_bcw_only);
double paswing_double(double alpha, double beta, double l, double pa0, double l0, int nrJumps, double *jump_longitude, double *jump_offset, double add_height_longitude, double add_height_shift, int add_height_shift_bcw_only);
void print_pgplot_version_used(FILE *stream);
void pgplot_clear_viewport_def(pgplot_viewport_definition *viewport);
void clear_pgplot_box(pgplot_box_definition *box);
void pgplot_clear_options(pgplot_options_definition *pgplot);
int initPulselongitudeRegion(pulselongitude_regions_definition *region, verbose_definition verbose);
void clearPulselongitudeRegion(pulselongitude_regions_definition *region);
void freePulselongitudeRegion(pulselongitude_regions_definition *region);
void copyPulselongitudeRegion(pulselongitude_regions_definition source, pulselongitude_regions_definition *destination);
void pgplot_setWindowsize(int windowwidth, int windowheight, float aspectratio);
void clear_pgplot_frame(pgplot_frame_def_internal *frame);
void pgplot_drawbox(pgplot_box_definition *box);
void clearRegion(pulselongitude_regions_definition *region);
void region_frac_to_int(pulselongitude_regions_definition *region, float scale, float offset);
void region_int_to_frac(pulselongitude_regions_definition *region, float scale, float offset);
int region_make_even(int regionnr, pulselongitude_regions_definition *region, int nrx);
void regionShowNextTimeUse(pulselongitude_regions_definition region, char *option, char *optionFrac, FILE *where);
int pgplot_opendevice(pgplot_viewport_definition *viewport, int *deviceID, verbose_definition verbose);
int pgplotGraph1(pgplot_options_definition *pgplot, float *data, float *datax, float *sigma, int nrx, float xmin, float xmax, int dontsetranges, float xmin_show, float xmax_show, float ymin_show, float ymax_show, int forceMinZero, int hist, int noline, int linewidth, int pointtype, int color, int boxcolor, pulselongitude_regions_definition *regions, int onpulsecolor, verbose_definition verbose);
int pgplotMap(pgplot_options_definition *pgplot, float *cmap, int nrx, int nry, float xmin, float xmax, float xminshow, float xmaxshow, float ymin, float ymax, float yminshow, float ymaxshow, int maptype, int itf, int nogray, int nrcontours, float *contours, int contourlw, int forceMinZero, float saturize, int levelset, float levelmin, float levelmax, int levelInversion, int onlyData, int sideright, int forceMinZeroRight, int sidetop, int forceMinZeroTop, int sidelw, int showwedge, int showwedge_actualmax, int plotSubset, int showTwice, verbose_definition verbose);
int pgplotMapCoordinate(float x, float y, int *nx, int *ny);
int pgplotMapCoordinate_dbl(double x, double y, int *nx, int *ny);
void pgplotMapCoordinateInverse(float *x, float *y, int nx, int ny);
void pgplotMapCoordinateInverse_dbl(double *x, double *y, int nx, int ny);
int pgplot_process_text_option(char *textoption, int outline_txt, int outline_lw, int outline_color, int argc, char **argv, int xunit_type, datafile_definition data, verbose_definition verbose);
void pgplotMapCoordinateBinSize(float *dx, float *dy);
int selectRegions(float *profileI, int nrBins, pgplot_options_definition *pgplot, int onlyOne, int powerTwo, int evenNumber, pulselongitude_regions_definition *regions, verbose_definition verbose);
int checkRegions(int bin, pulselongitude_regions_definition *regions, int whichregion, verbose_definition verbose);
int pgplot_device_type(char *devicename, verbose_definition verbose);
int pgplotPAplot(datafile_definition data, int showtotpol, int nopaswing, int showEll, pgplot_options_definition *pgplot, char *xlabel, char *ylabel, char *ylabel_pa, char *ylabel_ell, char *ylabel_spstatfrac, float longitude_left, float longitude_right, int xunit_type, float loffset, float Imin, float Imax, float pa_bottom, float pa_top, float PAoffset, float sigma_limit, float datalinewidth, float ysize2, int dashed, int noynumbers, char *textoption, char *textoption_pa, float ytick_pa, int nysub_pa, char *textoption_ell, float ytick_ell, int nysub_ell, char *textoption_padist, char *textoption_elldist, char *textoption_spstat, char *herrorbaroption, char *herrorbaroptionpa, char *herrorbaroptionpa2, char *herrorbaroptionell, char *herrorbaroptionpadist, char *herrorbaroptionelldist, char *verrorbaroption, char *verrorbaroptionpa, char *verrorbaroptionpa2, char *verrorbaroptionell, int argc, char **argv, int outline_txt, int outline_lw, int outline_color, int overlayPA, float overlayalpha, float overlaybeta, float overlaypa0, float overlayl0, int overlayPAfine, int nrJumps, float *jump_longitudes, float *jump_offsets, datafile_definition *padist, float padist_pamin, float padist_pamax, float padist_saturize, int padist_overlayavpa, int padist_paswing, datafile_definition *elldist, float elldist_saturize, int nowedge, datafile_definition *spstatfrac, verbose_definition verbose);
void drawSphericalGrid(float dlat, float dlong, float rot_long, float rot_lat, int lw, int projection);
void printCMAPCommandlineOptions(FILE *printdevice);
int cmap_parse_commandline(int argc, char **argv, int debug);
int ppgopen(const char *device);
void ppgclos(void);
void ppgend(void);
void ppgsvp(float xleft, float xright, float ybot, float ytop);
void ppgswin(float x1, float x2, float y1, float y2);
void ppgsci(int ci);
void ppgbbuf(void);
void ppgebuf(void);
void ppgmove(float x, float y);
void ppgdraw(float x, float y);
void ppgpt1(float xpt, float ypt, int symbol);
void ppgerr1(int dir, float x, float y, float e, float t);
void ppgask(int flag);
void ppgpage(void);
void ppgslw(int lw);
void ppgsclp(int state);
void ppgsls(int ls);
void ppgsfs(int fs);
void ppgsch(float size);
void ppgqcs(int units, float *xch, float *ych);
void ppglab(const char *xlbl, const char *ylbl, const char *toplbl);
void ppgbox(const char *xopt, float xtick, int nxsub, const char *yopt, float ytick, int nysub);
void ppgaxis(const char *opt, float x1, float y1, float x2, float y2, float v1, float v2, float step, int nsub, float dmajl, float dmajr, float fmin, float disp, float orient);
void ppgmtxt(const char *side, float disp, float coord, float fjust, const char *text);
void ppgtick(float x1, float y1, float x2, float y2, float v, float tikl, float tikr, float disp, float orient, const char *str);
void ppgscr(int ci, float cr, float cg, float cb);
void ppgsitf(int itf);
void ppgarro(float x1, float y1, float x2, float y2);
void ppgcirc(float xcent, float ycent, float radius);
void ppgtext(float x, float y, const char *text);
void ppgptxt(float x, float y, float angle, float fjust, const char *text);
void ppgscf(int font);
void ppgshs(float angle, float sepn, float phase);
void ppgrect(float x1, float x2, float y1, float y2);
void ppgerry(int n, const float *x, const float *y1, const float *y2, float t);
void ppggray(const float *a, int idim, int jdim, int i1, int i2, int j1, int j2, float fg, float bg, const float *tr);
void ppgimag(const float *a, int idim, int jdim, int i1, int i2, int j1, int j2, float a1, float a2, const float *tr);
void ppgcont(const float *a, int idim, int jdim, int i1, int i2, int j1, int j2, const float *c, int nc, const float *tr);
void ppgpoly(int n, const float *xpts, const float *ypts);
void ppgconl(const float *a, int idim, int jdim, int i1, int i2, int j1, int j2, float c, const float *tr, const char *label, int intval, int minint);
int ppgband(int mode, int posn, float xref, float yref, float *x, float *y, char *ch_scalar);
int ppgcurs(float *x, float *y, char *ch);
void ppgscir(int icilo, int icihi);
void ppgctab(float *l, float *r, float *g, float *b, int nc, float contra, float bright);
void ppgwedg(const char *side, float disp, float width, float fg, float bg, const char *label);
int ppgqid(int *id);
int ppgslct(int id);
void ppgpap(float width, float aspect);
int ppgqvp(int units, float *xleft, float *xright, float *ybot, float *ytop);
int ppgqwin(float *xleft, float *xright, float *ybot, float *ytop);
void ppgqinf(const char *item, char *value, int *value_length);
int readVonMisesModel(char *filename, vonMises_collection_definition *components, verbose_definition verbose);
int writeVonMisesModel(char *filename, vonMises_collection_definition *components, verbose_definition verbose);
double calcVonMisesFunction(vonMises_collection_definition *components, double phase, double shift);
double calcVonMisesFunction2(double centre, double concentration, double height, double phase, double shift);
double integrateVonMisesFunction(vonMises_collection_definition *components);
void calcVonMisesProfile(vonMises_collection_definition *components, int nrbins, float *profile, double shift, int normalize);
void calcVonMisesProfile_resid_rms(vonMises_collection_definition *components, int nrbins, float *profile, double shift, double *rms);
void calcVonMisesProfile_shape_parameter(vonMises_collection_definition *components, double shift, double phaseprecision, int shapepar, double *shapepar_aux, double *measurement, verbose_definition verbose);
void print_shape_par(FILE *fout, int showdescr, int shapepar, double measurement, double error);
float correlateVonMisesFunction(vonMises_collection_definition *components, int nrbins, float *profile, verbose_definition verbose);
void find_boundaries(float *profile, int nrbins, float y, pulselongitude_regions_definition *regions);
int fitvonmises_refine_model(datafile_definition psrdata, vonMises_collection_definition *components, int fitbaseline, int avoid_neg_components, float *baseline, int fixamp, int fixwidth, int fixphase, int fixrelamp, int fixrelphase, verbose_definition verbose);
int find_peak_correlation(float *data1, float *data2, int ndata, int zeropad, int circularpad, int duplicate, int remove_baseline, int *lag, float *correl_max, verbose_definition verbose);
void randomize_idnum(long *idnum);
long calculate_bin_number(double x, double dx, double min_x, int centered_at_zero, double extra_phase);
double calculate_bin_location(long binnr, double dx, double min_x, int centered_at_zero, double extra_phase);
double calculate_required_bin_width(double x, long binnr, double min_x, int centered_at_zero, double extra_phase, verbose_definition verbose);
int set_binning_histogram(double min_x_data, double max_x_data, int rangex_set, double rangex_min, double rangex_max, int nrbins_specified, long nrbins, int centered_at_zero, double extra_phase, double *min_x, double *max_x, double *dx, verbose_definition verbose);
void kstest(double *data1, long n1, double *data2, long n2, int cdf_type, double input_value1, double input_value2, double (*cdf)(double), double *max_diff, double *prob, verbose_definition verbose);
void print_gsl_version_used(FILE *stream);
int print_fitfunctions(fitfunc_collection_type *function, int novalue, int showerror, int index, verbose_definition verbose);
int fit_levmar(int algorithm, fitfunc_collection_type *function, double *data_x, double *data_y, double *data_sigma, long ndata, int oneatatime, int force_chi2_1, double epsabs, double epsrel, int maxiter, int *status, int showresults, int showcovariance, verbose_definition verbose);
double evaluate_fitfunc_collection(fitfunc_collection_type *function, double x, verbose_definition verbose);
int minimize_1D_double(int findroot, double (*funk)(double, void *), void *params, double x_lower, double x_upper, int gridsearch, int investigateLocalMinima, int nested, double *x_minimum, int max_iter, double epsabs, double epsrel, int verbose, int debug_verbose);
int find_1D_error(double (*funk)(double *, void *), double *xminimum, int paramnr, int nrparameters, double dx, double dxmax, void *params, double sigma, double chi2min, int max_itr, double epsabs, double epsrel, double *errorbar, int verbose);
int doAmoeba(int algorithm, float *xstart, float *dx, int *fixed, float *xfit, float *yfit, int nrparams, float (*funk)(float []), float ftol, int *nfunk, int verbose, int finderrors, float sigma, float *dplus, float *dmin);
int doAmoeba_d(int algorithm, double *xstart, double *dx, int *fixed, double *xfit, double *yfit, int nrparams, double (*funk)(double []), double ftol, int *nfunk, int verbose, int finderrors, double sigma, double *dplus, double *dmin);
int find_errors_amoeba(int algorithm, float *dx, int *fixed, float *xfit, float yfit, int nrparams, float (*funk)(float []), float ftol, int paramnr, float *dplus, float *dmin, float sigma);
int find_errors_amoeba_d(int algorithm, double *dx, int *fixed, double *xfit, double yfit, int nrparams, double (*funk)(double []), double ftol, int paramnr, double *dplus, double *dmin, double sigma);
int boxcarFindpeak(float *pulse, int nrBins, pulselongitude_regions_definition *onpulse, int *bin, int *pulsewidth, float *snrbest, float *E_best, int squared, int posOrNeg, int allwidths, int refine, int maxwidth, int only_onpulse, int nodebase, verbose_definition verbose);
void initApplication(psrsalsaApplication *application, char *name, char *genusage);
void terminateApplication(psrsalsaApplication *application);
void printApplicationHelp(psrsalsaApplication *application);
int processCommandLine(psrsalsaApplication *application, int argc, char **argv, int *index);
void printCitationInfo();
int parse_command_string(verbose_definition verbose, int argc, char **argv, int argv_index, int check_only, int minrequestedparameters, char *format, ...);
int preprocessApplication(psrsalsaApplication *application, datafile_definition *psrdata);
int applicationAddFilename(int argi, verbose_definition verbose);
int applicationFilenameList_checkConsecutive(char **argv, verbose_definition verbose);
int numberInApplicationFilenameList(psrsalsaApplication *application, char **argv, verbose_definition verbose);
char *getNextFilenameFromList(psrsalsaApplication *application, char **argv, verbose_definition verbose);
void rewindFilenameList(psrsalsaApplication *application);
int getOutputName(psrsalsaApplication *application, char *filename, char *outputname, verbose_definition verbose);
void showlibraryversioninformation(FILE *stream);
int pgetch(void);
void fprintf_color(FILE *destination, int color, const char *format, ...);
int pgetch_macro(psrsalsaApplication *application, verbose_definition verbose);
int getUsername(char **username, verbose_definition verbose);
void constructCommandLineString(char *txt, int length, int argc, char **argv, verbose_definition verbose);
int getMachinename(char *hostname, int size, verbose_definition verbose);
char * pickWordFromString(char *string, int n, int *nrwords, int replacetabs, char separator, verbose_definition verbose);
char *str_replace(char *orig, char *rep, char *with, verbose_definition verbose);
int ascii_file_get_next_line(FILE *fin, char *txt, int maxlinelength, int skipChar, verbose_definition verbose);
int ascii_file_stats(FILE *fin, char skipChar, long *nrlines, int maxlinelength, int autoNrColumns, int *nrColumns, verbose_definition verbose);
int change_filename_extension(char *inputname, char *outputname, char *extension, int outputnamelength, verbose_definition verbose);
int read_ascii_column(char *fname, int skiplines, char skipChar, int nrColumns, int autoNrColumns, long *nrdatapoints, long nrdatapoints_max, int colnum, double scale, double norm, int read_log, float **data, float *mindata, float *maxdata, float *avdata, verbose_definition verbose, int verbose_stderr);
int read_ascii_column_double(char *fname, int skiplines, char skipChar, int nrColumns, int autoNrColumns, long *nrdatapoints, long nrdatapoints_max, int colnum, double scale, double norm, int read_log, double **data, double *mindata, double *maxdata, double *avdata, verbose_definition verbose, int verbose_stderr);
int read_ascii_column_int(char *fname, int skiplines, char skipChar, int nrColumns, int autoNrColumns, long *nrdatapoints, long nrdatapoints_max, int colnum, int **data, int *mindata, int *maxdata, double *avdata, verbose_definition verbose, int verbose_stderr);
int read_ascii_column_str(char *fname, int skiplines, char skipChar, int nrColumns, int autoNrColumns, long *nrdatapoints, long nrdatapoints_max, int colnum, char ***data, verbose_definition verbose, int verbose_stderr);
int linalg_solve_matrix_eq_gauss_jordan(double *matrixa, double *matrixb, int n, int m, verbose_definition verbose);
