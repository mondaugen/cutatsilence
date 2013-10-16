/* Cut a soundfile at its silences. This implimentation is not optimized for
 * memory requirements and will have at some point in the program twice your
 * soundfile's length's worth of doubles in memory. Keep this in mind when using
 * the program on low memory systems. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h> 
#include <vector>
#include <utility> 
#include <stk/Fir.h>

#define START 1
#define ABOVE 2
#define BELOW 3
#define UP    4
#define DOWN  5 

/* Find the mid-point between two integers */
#define int_mid_point(a,b) (((a) + (b)) / (2))

/* Half-wave rectify an array of doubles (make them all positive) */

int
cas_rectify(double *samples, size_t length)
{
  while (length--)
    if (samples[length] < 0.)
      samples[length] *= -1.;
}

/* Create a windowed-sinc FIR filter with a given cutoff frequency and window
 * length. The cutoff frequency and sampling rate are in Hertz. The resulting
 * filter is entirely causal, which means x[n]'s coefficient is the first filter
 * coefficient, x[n-1] is the second filter coefficient and so on (you might
 * think that x[n]'s coefficient is the middle filter coefficient... but here
 * it's not). The window is a Hann window (for now). */
stk::Fir
cas_build_sinc_filter(size_t size, double cutoff, double sr)
{
  std::vector coefs;
  double rad_freq = cutoff / sr * 2 * M_PI;
  double x_max = ((double)(size - 1)) / 2.0;
  for (double x = -x_max; x <= x_max; ++x) {
    coefs.push_back((cos(x / (double)size * 2.0 * M_PI) + 1.0) / 2.
        * sin(rad_freq * x) / (M_PI * x));
  }
  stk::Fir result = stk::Fir(coefs);
  return result;
}

/* Create a windowed-unit-step FIR filter with a given cutoff frequency and
 * window length. The resulting filter is entirely causal, which means x[n]'s
 * coefficient is the first filter coefficient, x[n-1] is the second filter
 * coefficient and so on (you might think that x[n]'s coefficient is the middle
 * filter coefficient... but here it's not). The window is a Hann window (for
 * now). */
stk::Fir
cas_build_unitstep_filter(size_t size)
{
  std::vector coefs;
  double x_max = ((double)(size - 1)) / 2.0;
  for (double x = -x_max; x <= x_max; ++x) {
    coefs.push_back((cos(x / (double)size * 2.0 * M_PI) + 1.0) / 2.);
  }
  stk::Fir result = stk::Fir(coefs);
  return result;
}

/* Returns a vector of <size_t, int> pairs representing the sample-points at
 * which transitions occur in the signal sig. thresh is the transition
 * threshold. */
std::vector< std::pair<size_t, int> >
cas_find_trans_points(std::vector<double> &sig, size_t length, double thresh)
{
  std::vector< std::pair<size_t, int> > trans;
  int state = START;
  for (size_t i = 0; i < length; i++) {
    switch (state) {
      case START:
        if (sig[i] >= thresh)
          state = ABOVE;
        else
          state = BELOW;
        continue;
      case ABOVE:
        if (v[i] < thresh) {
          state = BELOW;
          trans.push_back( std::pair<size_t, int>(i, DOWN) );
        }
        continue;
      case BELOW:
        if (v[i] >= THRESH) {
          state = ABOVE;
          trans.push_back( std::pair<size_t, int>(i, UP) );
        }
        continue;
    }
  }
  return trans;
}

/* Returns a vector of <size_t>s representing the points at which to cut the
 * sound file. Needs a vector of pairs from cas_find_trans_points(). Does not
 * add the first or last point as a cut point. */
std::vector< size_t >
cas_find_cut_points( std::vector< std::pair<size_t, int> > &trans )
{
  std::vector< size_t > cutpoints;
  /* The last transition point, so we can calculate the mid-point which is gonna
   * be the cut-point */
  size_t lastsamp = 0;
  int state = START;
  for (size_t i = 0; i < trans.size(); i++) {
    switch (state) {
      case START:
        state = trans[i].second;
        lastsamp = trans[i].first;
        continue;
      case UP:
        state = trans[i].second;
        lastsamp = trans[i].first;
        continue;
      case DOWN:
        state = trans[i].second;
        if (state == UP)
          cutpoints.push_back( int_mid_point(lastsamp, trans[i].first) );
        lastsamp = trans[i].first;
        continue;
    }
  }
  return cutpoints;
}

/* Using a vector of sorted cut-points and a file name format, write all the
 * cut-out sections of a signal to files with unique names. The cut-points
 * vector must be sorted and we assume that all the cut-points are valid
 * addresses in sig. We also assume that there are at least the addresses of the
 * first and last sample in the soundfile. The filename format is in that class
 * c-style and may only contain 1 %d field. It should also have the proper file
 * ending corresponding to the format you chose. See sndfile.h in libsndfile
 * sources for different valid formats. sr is the sampling rate of the
 * soundfiles in hz*/
int
cas_write_signal_sections(std::vector<double> &sig, const char *path_form,
    std::vector< size_t > &cutpoints, int sf_format, int sr)
{
  size_t lastcut = cutpoints[0];
  char *path;
  for (size_t i = 1; i < cutpoints.size(); i++) {
    if (asprintf(&path, path_form, (int)i) < 0)
      return -2; /* couldn't allocate string for the path name */
    SndfileHandle sfHandle(path, SFM_WRITE, sf_format, 1, sr);
    if (sfHandle.write(&(sig.front()) + lastcut,
          cutpoints[i] - lastcut) < (cutpoints[i] - lastcut))
      return -1 /* couldn't write all the samples! */
    lastcut = cutpoints[i];
    free(path); /* deallocate path name */
  }
  return 0;
}

int
main (int argc, char **argv)
{
  /* Argument order:
   *  path to sound file
   *  output path string format
   *  threshold (in dB)
   *  cut-off frequency (in Hz)
   *  window size
   */
  char *inpath          = argv[1];
  char *outpathstrform  = argv[2];
  double thresh         = pow(2.0, atof(argv[3]) / 6.0);
  double cofreq         = atof(argv[4]);
  size_t winsize        = atoi(argv[5]); /* assume positive */

  if (!(winsize % 2)) {
    fprintf("Window size must be odd.\n");
    exit(1);
  }

  /* open sound file */
  SndfileHandle infile(path);

  if (infile.channels() > 1) {
    fprintf(stderr, "Mono files only.\n");
    exit(1);
  }

  std::vector<double> origbuff(infile.frames() + winsize, 0);
  std::vector<double> procbuff(infile.frames() + winsize, 0);

  /* load file into origbuff. We delay by half a window to align with the
   * filtered signal */
  if (infile.frames() > infile.read((&(origbuff.front()))
        + (winsize / 2) ,
      infile.frames())) {
    fprintf("Failed to load soundfile %s into buffer.\n", inpath);
    exit(1);
  }

  /* copy file to buffer that will be filtered */
  memcpy(&(procbuff.front()), ((void*)(&(origbuff.front())))
      + (winsize / 2 * sizeof(double)),
      infile.frames() * sizeof(double));

  /* build a filter */
  stk::Fir filter = cas_build_unitstep_filter(winsize);

  /* do filtering to procbuf */
  for ( std::vector<double>::iterator it = procbuff.begin();
        it != procbuff.end();
        ++it ) {
    *it = filter.tick( StkFrames( *it ) );
  }

  /* find the transition points in soundfile */
  std::vector< std::pair<size_t, int> > tpoints = 
    cas_find_trans_points(procbuff, procbuff.size(), thresh);

  /* find the cut-points in the soundfile */
  std::vector<size_t> cpoints = cas_find_cut_points(tpoints);

  /* write the signal sections. Default format is 16 bit wave for now. */
  cas_write_signal_sections(origbuff, outpathstrform, cpoints,
      (SF_FORMAT_WAV, SF_FORMAT_PCM_16), infile.samplerate());

  exit(0);
}

