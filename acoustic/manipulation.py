"""Helper functions to manipulate sound waves based on Praat and parselmouth."""

from parselmouth.praat import call


def replace_pitch_tier(original_sound, replicated_sound,
                       time_step=0.01,
                       pitch_floor=75., pitch_ceiling=600.,
                       output_filepath=None):
    """
    Function to replace pitch tier of a sound with another given sound.

    :param (parselmouth.Sound) original_sound: original sound waveform
    :param (parselmouth.Sound) replicated_sound: to-be-replicated sound waveform
    :param (float) time_step: the measurement interval (frame duration), in seconds (default: 0.01)
    :param (float) pitch_floor: minimum pitch (default: 75.)
    :param (float) pitch_ceiling:: maximum pitch (default: 600.)
    :param output_filepath: path to the output file as the new sound waveform (default: None)
    :return: new Sound (parselmouth.Sound) object with replaced pitch tier
    """
    # Create Manipulation object for original sound
    manipulation = call(original_sound, 'To Manipulation',
                        time_step,
                        pitch_floor, pitch_ceiling)

    # Create Manipulation object for replicated sound
    replicated_manipulation = call(replicated_sound, 'To Manipulation',
                                   time_step,
                                   pitch_floor, pitch_ceiling)

    # Create a PitchTier object for replicated Manipulation object
    replicated_pitch_tier = call(replicated_manipulation, 'Extract pitch tier')

    call([replicated_pitch_tier, manipulation], 'Replace pitch tier')

    # Resynthesize the new sound
    new_sound = call(manipulation, 'Get resynthesis (overlap-add)')

    if output_filepath is not None:
        call(new_sound, 'Save as WAV file', output_filepath)

    return new_sound


def replace_intensity_tier(original_sound, replicated_sound,
                           time_step=0.01,
                           pitch_floor=75.,
                           output_filepath=None):
    """
    Function to replace pitch tier of a sound with another given sound.

    :param (parselmouth.Sound) original_sound: original sound waveform
    :param (parselmouth.Sound) replicated_sound: to-be-replicated sound waveform
    :param (float) time_step: the measurement interval (frame duration), in seconds (default: 0.01)
    :param (float) pitch_floor: minimum pitch (default: 75.)
    :param output_filepath: path to the output file as the new sound waveform (default: None)
    :return: new Sound (parselmouth.Sound) object with replaced intensity tier
    """
    # Create Intensity object for replicated sound
    replicated_intensity = call(replicated_sound, 'To Intensity',
                                pitch_floor, time_step)

    # Create IntensityTier object for replicated Intensity object
    replicated_intensity_tier = call(replicated_intensity, 'To IntensityTier (peaks)')

    # Multiply intensities of the original sound to effectively replace intensity tiers.
    # NOTE: The boolean argument 'yes' below scales the multiplication by a 0.9 factor
    new_sound = call([replicated_intensity_tier, original_sound], 'Multiply...',
                     'yes')

    if output_filepath is not None:
        call(new_sound, 'Save as WAV file', output_filepath)

    return new_sound


def replace_pulses(original_sound, replicated_sound,
                   time_step=0.01,
                   pitch_floor=75., pitch_ceiling=600.,
                   output_filepath=None):
    """
    Function to replace pulses of a sound with another given sound.
    NOTE: This often leads to unnatural speech synthesis, unlike the previous 2 functions!

    :param (parselmouth.Sound) original_sound: original sound waveform
    :param (parselmouth.Sound) replicated_sound: to-be-replicated sound waveform
    :param (float) time_step: the measurement interval (frame duration), in seconds (default: 0.01)
    :param (float) pitch_floor: minimum pitch (default: 75.)
    :param (float) pitch_ceiling: maximum pitch (default: 600.)
    :param output_filepath: path to the output file as the new sound waveform (default: None)
    :return: new Sound (parselmouth.Sound) object with replaced pulses
    """
    # Create Manipulation object from original sound
    manipulation = call(original_sound, 'To Manipulation',
                        time_step,
                        pitch_floor, pitch_ceiling)

    # Create Manipulation object from replicated sound
    replicated_manipulation = call(replicated_sound, 'To Manipulation',
                                   time_step,
                                   pitch_floor, pitch_ceiling)

    # Create PointProcess object of pulses from replicated Manipulation object
    replicated_pulses = call(replicated_manipulation, 'Extract pulses')

    call([replicated_pulses, manipulation], 'Replace pulses')

    # Resynthesize the new sound
    new_sound = call(manipulation, 'Get resynthesis (overlap-add)')

    if output_filepath is not None:
        call(new_sound, 'Save as WAV file', output_filepath)

    return new_sound
