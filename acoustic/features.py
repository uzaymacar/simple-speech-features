"""Helper functions for extracting features from sound waves based on Praat and parselmouth."""

from parselmouth.praat import call


def get_pitch_attributes(sound, pitch_type='preferred',
                         time_step=0., min_time=0., max_time=0.,
                         pitch_floor=75., pitch_ceiling=600.,
                         unit='Hertz', interpolation_method='Parabolic'):
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
    :return: a dictionary of mentioned attributes
    """

    if pitch_type == 'preferred':
        pitch = call(sound, 'To Pitch', time_step, pitch_floor, pitch_ceiling)
    elif pitch_type == 'cc':
        pitch = call(sound, 'To Pitch (cc)', time_step, pitch_floor, pitch_ceiling)
    else:
        raise ValueError('Argument for @pitch_type not recognized!')

    attributes = dict()

    attributes['min_pitch'] = call(pitch, 'Get minimum',
                                   min_time, max_time,
                                   unit, interpolation_method)
    attributes['max_pitch'] = call(pitch, 'Get maximum',
                                   min_time, max_time,
                                   unit, interpolation_method)
    attributes['mean_pitch'] = call(pitch, 'Get mean',
                                    min_time, max_time,
                                    unit)
    attributes['stddev_pitch'] = call(pitch, 'Get standard deviation',
                                      min_time, max_time, unit)

    return attributes


def get_intensity_attributes(sound,
                             time_step=0., min_time=0., max_time=0.,
                             pitch_floor=75.,
                             interpolation_method='Parabolic'):
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
    :return: a dictionary of mentioned attributes
    """
    intensity = call(sound, 'To Intensity', pitch_floor, time_step)

    attributes = dict()

    attributes['min_intensity'] = call(intensity, 'Get minimum',
                                       min_time, max_time,
                                       interpolation_method)
    attributes['max_intensity'] = call(intensity, 'Get maximum',
                                       min_time, max_time,
                                       interpolation_method)
    attributes['mean_intensity'] = call(intensity, 'Get mean',
                                        min_time, max_time)
    attributes['stddev_intensity'] = call(intensity, 'Get standard deviation',
                                          min_time, max_time)

    return attributes


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


def get_harmonics_to_noise_ratio_attributes(sound,
                                            harmonics_type='preferred',
                                            time_step=0.01, min_time=0., max_time=0.,
                                            minimum_pitch=75.,
                                            silence_threshold=0.1, num_periods_per_window=1.0,
                                            interpolation_method='Parabolic'):
    """
    Function to get Harmonics-to-Noise Ratio (HNR) attributes such as minimum HNR, maximum HNR,
    mean HNR, and standard deviation of HNR.

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
    :return: a dictionary of mentioned attributes
    """
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
    attributes['max_hnr'] = call(harmonicity, 'Get maximum',
                                 min_time, max_time,
                                 interpolation_method)
    attributes['mean_hnr'] = call(harmonicity, 'Get mean',
                                  min_time, max_time)
    attributes['stddev_hnr'] = call(harmonicity, 'Get standard deviation',
                                    min_time, max_time)

    return attributes
