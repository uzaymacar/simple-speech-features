"""
Helper functions for extracting features from sound waves based on Praat and parselmouth.
Some methods are inspired from: http://www.fon.hum.uva.nl/rob/NKI_TEVA/TEVA/HTML/Analysis.html
"""

import math
import parselmouth

from parselmouth.praat import call


def get_intensity_attributes(sound,
                             time_step=0., min_time=0., max_time=0.,
                             pitch_floor=75.,
                             interpolation_method='Parabolic',
                             return_values=False, replacement_for_nan=0.):
    """
    Function to get intensity attributes such as minimum intensity, maximum intensity, mean
    intensity, and standard deviation of intensity.

    NOTE: Notice that we don't need a unit parameter for intensity as intensity is consistently
    reported as dB SPL throughout Praat. dB SPL is simply dB relative to the normative auditory
    threshold for a 1000-Hz sine wave, 2 x 10^(-5) Pascal.

    NOTE: The standard interpolation method is 'Parabolic' because of the usual non-linearity
    (logarithm) in the computation of intensity; sinc interpolation would be too stiff and may
    give unexpected results.

    :param (parselmouth.Sound) sound: sound waveform
    :param (float) time_step: the measurement interval (frame duration), in seconds (default: 0.)
           NOTE: The default 0. value corresponds to a time step of 0.75 / pitch floor
    :param (float) min_time: minimum time value considered for time range (t1, t2) (default: 0.)
    :param (float) max_time: maximum time value considered for time range (t1, t2) (default: 0.)
           NOTE: If max_time <= min_time, the entire time domain is considered
    :param pitch_floor: minimum pitch (default: 75.)
    :param (str) interpolation_method: method of sampling new data points with a discrete set of
           known data points, 'None', 'Parabolic', 'Cubic', or 'Sinc' (default: 'Parabolic')
    :param (bool) return_values: whether to return a continuous list of intensity values
           from all frames or not
    :param (float) replacement_for_nan: a float number that will represent frames with NaN values
    :return: (a dictionary of mentioned attributes, a list of intensity values OR None)
    """
    # Get total duration of the sound
    duration = call(sound, 'Get end time')

    # Create Intensity object
    intensity = call(sound, 'To Intensity', pitch_floor, time_step)

    attributes = dict()

    attributes['min_intensity'] = call(intensity, 'Get minimum',
                                       min_time, max_time,
                                       interpolation_method)
    attributes['relative_min_intensity_time'] = call(intensity, 'Get time of minimum',
                                                     min_time, max_time,
                                                     interpolation_method) / duration
    attributes['max_intensity'] = call(intensity, 'Get maximum',
                                       min_time, max_time,
                                       interpolation_method)
    attributes['relative_max_intensity_time'] = call(intensity, 'Get time of maximum',
                                                     min_time, max_time,
                                                     interpolation_method) / duration

    attributes['mean_intensity'] = call(intensity, 'Get mean',
                                        min_time, max_time)
    attributes['stddev_intensity'] = call(intensity, 'Get standard deviation',
                                          min_time, max_time)

    attributes['q1_intensity'] = call(intensity, 'Get quantile',
                                      min_time, max_time,
                                      0.25)
    attributes['median_intensity'] = call(intensity, 'Get quantile',
                                          min_time, max_time,
                                          0.50)
    attributes['q3_intensity'] = call(intensity, 'Get quantile',
                                      min_time, max_time,
                                      0.75)

    intensity_values = None

    if return_values:
        intensity_values = [call(intensity, 'Get value in frame', frame_no)
                            for frame_no in range(len(intensity))]
        # Convert NaN values to floats (default: 0)
        intensity_values = [value if not math.isnan(value) else replacement_for_nan
                            for value in intensity_values]

    return attributes,  intensity_values


def get_pitch_attributes(sound, pitch_type='preferred',
                         time_step=0., min_time=0., max_time=0.,
                         pitch_floor=75., pitch_ceiling=600.,
                         unit='Hertz', interpolation_method='Parabolic',
                         return_values=False, replacement_for_nan=0.):
    """
    Function to get pitch attributes such as minimum pitch, maximum pitch, mean pitch, and
    standard deviation of pitch.

    :param (parselmouth.Sound) sound: sound waveform
    :param (str) pitch_type: the type of pitch analysis to be performed; values include 'preferred'
           optimized for speech based on auto-correlation method, and 'cc' for performing acoustic
           periodicity detection based on cross-correlation method
           NOTE: Praat also includes an option for type 'ac', a variation of 'preferred' that
           requires several more parameters. We are not including this for simplification.
    :param (float) time_step: the measurement interval (frame duration), in seconds (default: 0.)
           NOTE: The default 0. value corresponds to a time step of 0.75 / pitch floor
    :param (float) min_time: minimum time value considered for time range (t1, t2) (default: 0.)
    :param (float) max_time: maximum time value considered for time range (t1, t2) (default: 0.)
           NOTE: If max_time <= min_time, the entire time domain is considered
    :param (float) pitch_floor: minimum pitch (default: 75.)
    :param (float) pitch_ceiling: maximum pitch (default: 600.)
    :param (str) unit: units of the result, 'Hertz' or 'Bark' (default: 'Hertz)
    :param (str) interpolation_method: method of sampling new data points with a discrete set of
           known data points, 'None' or 'Parabolic' (default: 'Parabolic')
    :param (bool) return_values: whether to return a continuous list of pitch values from all frames
           or not
    :param (float) replacement_for_nan: a float number that will represent frames with NaN values
    :return: (a dictionary of mentioned attributes, a list of pitch values OR None)
    """
    # Get total duration of the sound
    duration = call(sound, 'Get end time')

    # Create pitch object
    if pitch_type == 'preferred':
        pitch = call(sound, 'To Pitch', time_step, pitch_floor, pitch_ceiling)
    elif pitch_type == 'cc':
        pitch = call(sound, 'To Pitch (cc)', time_step, pitch_floor, pitch_ceiling)
    else:
        raise ValueError('Argument for @pitch_type not recognized!')

    attributes = dict()

    attributes['voiced_fraction'] = call(pitch, 'Count voiced frames') / len(pitch)

    attributes['min_pitch'] = call(pitch, 'Get minimum',
                                   min_time, max_time,
                                   unit,
                                   interpolation_method)
    attributes['relative_min_pitch_time'] = call(pitch, 'Get time of minimum',
                                                 min_time, max_time,
                                                 unit,
                                                 interpolation_method) / duration
    attributes['max_pitch'] = call(pitch, 'Get maximum',
                                   min_time, max_time,
                                   unit,
                                   interpolation_method)
    attributes['relative_max_pitch_time'] = call(pitch, 'Get time of maximum',
                                                 min_time, max_time,
                                                 unit,
                                                 interpolation_method) / duration

    attributes['mean_pitch'] = call(pitch, 'Get mean',
                                    min_time, max_time,
                                    unit)
    attributes['stddev_pitch'] = call(pitch, 'Get standard deviation',
                                      min_time, max_time,
                                      unit)

    attributes['q1_pitch'] = call(pitch, 'Get quantile',
                                  min_time, max_time,
                                  0.25)
    attributes['median_intensity'] = call(pitch, 'Get quantile',
                                          min_time, max_time,
                                          0.50)
    attributes['q3_pitch'] = call(pitch, 'Get quantile',
                                  min_time, max_time,
                                  0.75)

    attributes['mean_absolute_pitch_slope'] = call(pitch, 'Get mean absolute slope',
                                                   unit)
    attributes['pitch_slope_without_octave_jumps'] = call(pitch, 'Get slope without octave jumps')

    pitch_values = None

    if return_values:
        pitch_values = [call(pitch, 'Get value in frame', frame_no, unit)
                        for frame_no in range(len(pitch))]
        # Convert NaN values to floats (default: 0)
        pitch_values = [value if not math.isnan(value) else replacement_for_nan
                        for value in pitch_values]

    return attributes, pitch_values


def get_harmonics_to_noise_ratio_attributes(sound,
                                            harmonics_type='preferred',
                                            time_step=0.01, min_time=0., max_time=0.,
                                            minimum_pitch=75.,
                                            silence_threshold=0.1, num_periods_per_window=1.0,
                                            interpolation_method='Parabolic',
                                            return_values=False,
                                            replacement_for_nan=0.):
    """
    Function to get Harmonics-to-Noise Ratio (HNR) attributes such as minimum HNR, maximum HNR,
    mean HNR, and standard deviation of HNR. HNR is defined as a measure that quantifies the amount
    of additive noise in a voice signal.

    NOTE: Harmonicity object represents the degree of acoustic periodicity, also called
    Harmonics-to-Noise Ratio (HNR). Harmonicity is expressed in dB: if 99% of the energy of the
    signal is in the periodic part, and 1% is noise, the HNR is 10*log10(99/1) = 20 dB. A HNR of
    0 dB means that there is equal energy in the harmonics and in the noise.

    :param (parselmouth.Sound) sound: sound waveform
    :param (str) harmonics_type: the type of harmonicity analysis to be performed; values include
           'preferred' for short-term analysis on cross-correlation method, and 'ac' for performing
           acoustic periodicity detection based on an accurate auto-correlation method
    :param (float) time_step: the measurement interval (frame duration), in seconds (default: 0.01)
    :param (float) min_time: minimum time value considered for time range (t1, t2) (default: 0.)
    :param (float) max_time: maximum time value considered for time range (t1, t2) (default: 0.)
           NOTE: If max_time <= min_time, the entire time domain is considered
    :param (float) minimum_pitch: determines the length of the analysis window (default: 75.)
    :param (float) silence_threshold: frames that do not contain amplitudes above this threshold
           (relative to the global maximum amplitude), are considered silent (default: 0.1)
    :param (float) num_periods_per_window: 4.5 is usually best for speech; HNR values up to 37 dB
           are guaranteed to be detected reliably; 6 periods per window raises this figure to more
           than 60 dB, but the algorithm becomes more sensitive to dynamic changes in the signal
           (default: 1.0)
    :param (str) interpolation_method: method of sampling new data points with a discrete set of
           known data points, 'None', 'Parabolic', 'Cubic', 'Sinc70', or 'Sinc700'
           (default: 'Parabolic')
    :param (bool) return_values: whether to return a continuous list of harmonicity values
           from all frames or not
    :param (float) replacement_for_nan: a float number that will represent frames with NaN values
    :return: (a dictionary of mentioned attributes, a list of harmonicity values OR None)
    """
    # Get total duration of the sound
    duration = call(sound, 'Get end time')

    # Create a Harmonicity object
    if harmonics_type == 'preferred':
        harmonicity = call(sound, "To Harmonicity (cc)",
                           time_step,
                           minimum_pitch,
                           silence_threshold, num_periods_per_window)
    elif harmonics_type == 'ac':
        harmonicity = call(sound, "To Harmonicity (ac)",
                           time_step,
                           minimum_pitch,
                           silence_threshold, num_periods_per_window)
    else:
        raise ValueError('Argument for @harmonics_type is not recognized!')

    attributes = dict()

    attributes['min_hnr'] = call(harmonicity, 'Get minimum',
                                 min_time, max_time,
                                 interpolation_method)
    attributes['relative_min_hnr_time'] = call(harmonicity, 'Get time of minimum',
                                               min_time, max_time,
                                               interpolation_method) / duration
    attributes['max_hnr'] = call(harmonicity, 'Get maximum',
                                 min_time, max_time,
                                 interpolation_method)
    attributes['relative_max_hnr_time'] = call(harmonicity, 'Get time of maximum',
                                               min_time, max_time,
                                               interpolation_method) / duration

    attributes['mean_hnr'] = call(harmonicity, 'Get mean',
                                  min_time, max_time)
    attributes['stddev_hnr'] = call(harmonicity, 'Get standard deviation',
                                    min_time, max_time)

    harmonicity_values = None

    if return_values:
        harmonicity_values = [call(harmonicity, 'Get value in frame', frame_no)
                              for frame_no in range(len(harmonicity))]
        # Convert NaN values to floats (default: 0)
        harmonicity_values = [value if not math.isnan(value) else replacement_for_nan
                              for value in harmonicity_values]

    return attributes, harmonicity_values


def get_glottal_to_noise_ratio_attributes(sound,
                                          minimum_frequency=500., maximum_frequency=4500.,
                                          bandwidth=1000., step=80.):
    """
    Function to get Glottal-to-Noise Ratio (GNE) attributes such as minimum GNE, maximum GNE,
    mean GNE, standard deviation of GNE, and sum of GNE. GNE is a measure that indicates whether a
    given voice signal originates from vibrations of the vocal folds or from turbulent noise
    generated in the vocal tract and is thus related to (but not a direct measure of) breathiness.
    (D.Michaelis et al. 1997)

    NOTE: The default units for the operations performed in this function are all 'Hertz'.

    :param (parselmouth.Sound) sound: sound waveform
    :param (float) minimum_frequency: minimum frequency for analysis (default: 500.)
    :param (float) maximum_frequency: maximum frequency for analysis (default: 4500.)
    :param (float) bandwidth: frequency difference between upper and lower signals (default: 1000.)
    :param (float) step: frequency steps for intervals (default: 80.)
    :return: a dictionary of mentioned attributes
    """
    # Create a Matrix object that represents GNE
    matrix = call(sound, "To Harmonicity (gne)",
                  minimum_frequency, maximum_frequency,
                  bandwidth, step)

    attributes = dict()

    attributes['min_gne'] = call(matrix, 'Get minimum')
    attributes['max_gne'] = call(matrix, 'Get maximum')

    attributes['mean_gne'] = call(matrix, 'Get mean')
    attributes['stddev_gne'] = call(matrix, 'Get standard deviation')

    attributes['sum_gne'] = call(matrix, 'Get sum')

    return attributes


def get_local_jitter(sound,
                     min_time=0., max_time=0.,
                     pitch_floor=75., pitch_ceiling=600., period_floor=0.0001, period_ceiling=0.02,
                     max_period_factor=1.3):
    """
    Function to calculate (local) jitter from a periodic PointProcess.

    :param (parselmouth.Sound) sound: sound waveform
    :param (float) min_time: minimum time value considered for time range (t1, t2) (default: 0.)
    :param (float) max_time: maximum time value considered for time range (t1, t2) (default: 0.)
           NOTE: If max_time <= min_time, the entire time domain is considered
    :param (float) pitch_floor: minimum pitch (default: 75.)
    :param (float) pitch_ceiling: maximum pitch (default: 600.)
    :param (float) period_floor: the shortest possible interval that will be used in the computation
           of jitter, in seconds (default: 0.0001)
    :param (float) period_ceiling: the longest possible interval that will be used in the
           computation of jitter, in seconds (default: 0.02)
    :param (float) max_period_factor: the largest possible difference between consecutive intervals
           that will be used in the computation of jitter (default: 1.3)
    :return: value of (local) jitter
    """
    # Create a PointProcess object
    point_process = call(sound, 'To PointProcess (periodic, cc)',
                         pitch_floor, pitch_ceiling)

    local_jitter = call(point_process, 'Get jitter (local)',
                        min_time, max_time,
                        period_floor, period_ceiling, max_period_factor)

    return local_jitter


def get_local_shimmer(sound,
                      min_time=0., max_time=0.,
                      pitch_floor=75., pitch_ceiling=600., period_floor=0.0001, period_ceiling=0.02,
                      max_period_factor=1.3, max_amplitude_factor=1.6):
    """
    Function to calculate (local) shimmer from a periodic PointProcess.

    :param (parselmouth.Sound) sound: sound waveform
    :param (float) min_time: minimum time value considered for time range (t1, t2) (default: 0.)
    :param (float) max_time: maximum time value considered for time range (t1, t2) (default: 0.)
           NOTE: If max_time <= min_time, the entire time domain is considered
    :param (float) pitch_floor: minimum pitch (default: 75.)
    :param (float) pitch_ceiling: maximum pitch (default: 600.)
    :param (float) period_floor: the shortest possible interval that will be used in the computation
           of shimmer, in seconds (default: 0.0001)
    :param (float) period_ceiling: the longest possible interval that will be used in the
           computation of shimmer, in seconds (default: 0.02)
    :param (float) max_period_factor: the largest possible difference between consecutive intervals
           that will be used in the computation of shimmer (default: 1.3)
    :param (float) max_amplitude_factor: maximum amplitude factor for shimmer (default: 1.6)
    :return: value of (local) shimmer
    """
    # Create a PointProcess object
    point_process = call(sound, 'To PointProcess (periodic, cc)',
                         pitch_floor, pitch_ceiling)

    local_shimmer = call([sound, point_process], 'Get shimmer (local)',
                         min_time, max_time,
                         period_floor, period_ceiling, max_period_factor, max_amplitude_factor)

    return local_shimmer


def get_spectrum_attributes(sound,
                            band_floor=200., band_ceiling=1000.,
                            low_band_floor=0., low_band_ceiling=500.,
                            high_band_floor=500., high_band_ceiling=4000.,
                            power=2., moment=3.,
                            return_values=False, replacement_for_nan=0.):
    """
    Function to get pitch attributes such as minimum pitch, maximum pitch, mean pitch, and
    standard deviation of pitch.

    :param (parselmouth.Sound) sound: sound waveform
    :param (str) pitch_type: the type of pitch analysis to be performed; values include 'preferred'
           optimized for speech based on auto-correlation method, and 'cc' for performing acoustic
           periodicity detection based on cross-correlation method
           NOTE: Praat also includes an option for type 'ac', a variation of 'preferred' that
           requires several more parameters. We are not including this for simplification.
    :param (float) time_step: the measurement interval (frame duration), in seconds (default: 0.)
           NOTE: The default 0. value corresponds to a time step of 0.75 / pitch floor
    :param (float) min_time: minimum time value considered for time range (t1, t2) (default: 0.)
    :param (float) max_time: maximum time value considered for time range (t1, t2) (default: 0.)
           NOTE: If max_time <= min_time, the entire time domain is considered
    :param (float) pitch_floor: minimum pitch (default: 75.)
    :param (float) pitch_ceiling: maximum pitch (default: 600.)
    :param (str) unit: units of the result, 'Hertz' or 'Bark' (default: 'Hertz)
    :param (str) interpolation_method: method of sampling new data points with a discrete set of
           known data points, 'None' or 'Parabolic' (default: 'Parabolic')
    :param (bool) return_values: whether to return a continuous list of pitch values from all frames
           or not
    :param (float) replacement_for_nan: a float number that will represent frames with NaN values
    :return: (a dictionary of mentioned attributes, a list of pitch values OR None)
    """
    # Create a Spectrum object
    spectrum = call(sound, 'To spectrum',
                    'yes')

    attributes = dict()

    attributes['band_energy'] = call(spectrum, 'Get band energy',
                                     band_floor, band_ceiling)
    attributes['band_density'] = call(spectrum, 'Get band density',
                                      band_floor, band_ceiling)
    attributes['band_energy_difference'] = call(spectrum, 'Get band energy difference',
                                                low_band_floor, low_band_ceiling,
                                                high_band_floor, high_band_ceiling)
    attributes['band_density_difference'] = call(spectrum, 'Get band density difference',
                                                 low_band_floor, low_band_ceiling,
                                                 high_band_floor, high_band_ceiling)

    attributes['center_of_gravity_spectrum'] = call(spectrum, 'Get centre of gravity',
                                                    power)
    attributes['stddev_spectrum'] = call(spectrum, 'Get standard deviation',
                                         power)
    attributes['skewness_spectrum'] = call(spectrum, 'Get skewness',
                                           power)
    attributes['kurtosis_spectrum'] = call(spectrum, 'Get kurtosis',
                                           power)

    # TODO: Dig deeper into what the moments represent here and what their values can be!
    attributes['central_moment_spectrum'] = call(spectrum, 'Get central moment',
                                                 moment)

    spectrum_values = None

    if return_values:
        spectrum_values = [call(spectrum, 'Get real value in bin', bin_no)
                           for bin_no in range(len(spectrum))]
        # Convert NaN values to floats (default: 0)
        spectrum_values = [val if not math.isnan(val) else replacement_for_nan
                           for val in spectrum_values]

    return attributes, spectrum_values



# TODO: Have a formant method based on the link on top of this script! (Formant quality factors
# TODO: and much more through the Formant object)

def get_speaking_rate(sound, text_filepath):
    """
    Function to get speaking rate, approximated as number of words divided by total duration.

    :param (parselmouth.Sound) sound: sound waveform
    :param (str) text_filepath: path to text file that annotates the given speech, @sound
    :return: speaking rate
    """
    # Get total duration of the sound
    duration = call(sound, 'Get end time')

    with open(text_filepath, mode='r') as f:
        text = f.read()

    # Approximate speaking rate as #words / duration
    speaking_rate = len(text.split()) / duration

    return speaking_rate
